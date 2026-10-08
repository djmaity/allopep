#!/usr/bin/env python3
"""AlloPep: allosteric peptide prediction pipeline."""

import argparse
import json
import os
import shutil
import subprocess
from pathlib import Path

import pandas as pd

import count_clashes
import run_gaps
from pdb_utils import split_protein_peptide


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog='python allopep.py',
        description='Pipeline for allosteric peptide design',
    )
    parser.add_argument('pdb_file')
    parser.add_argument('--max_rank', default=None, type=int,
                        help='Maximum ranked GAPS patch to use for peptide design')
    parser.add_argument('--length_min', default=5, type=int,
                        help='Minimum length of the predicted peptide')
    parser.add_argument('--length_max', default=10, type=int,
                        help='Maximum length of the predicted peptide')
    parser.add_argument('--num_samples', default=10, type=int,
                        help='Number of predicted peptides per patch')
    parser.add_argument('--peptide-chain', default=None,
                        help='Known peptide chain ID; automatically detected if omitted')
    return parser.parse_args(argv)


def run_pepglad(args, pdb_file, out_dir, prefix, patches):
    sites_dir = out_dir / f'{prefix}_allosteric_sites'
    sites_dir.mkdir(parents=True, exist_ok=True)

    env = os.environ.copy()
    env['RAY_ACCEL_ENV_VAR_OVERRIDE_ON_ZERO'] = '0'

    max_rank = len(patches)
    if args.max_rank is not None and 0 < args.max_rank <= max_rank:
        max_rank = args.max_rank

    generated = []
    for rank, patch in enumerate(patches[:max_rank], start=1):
        pocket_file = sites_dir / f'{prefix}_allosteric_site_{rank}.json'
        with pocket_file.open('w', encoding='utf-8') as handle:
            json.dump(patch.as_pepglad_pocket(), handle)

        codesign_dir = (
            out_dir / f'{prefix}_PepGLAD_outputs'
            / f'{prefix}_allosteric_site_{rank}_codesign'
        )
        subprocess.run([
            'python', '-m', 'api.run', '--mode', 'codesign',
            '--pdb', pdb_file,
            '--pocket', f'../{pocket_file}',
            '--out_dir', f'../{codesign_dir}',
            '--length_min', str(args.length_min),
            '--length_max', str(args.length_max),
            '--n_samples', str(args.num_samples),
        ], cwd='PepGLAD', env=env, check=True)

        for sample_index in range(args.num_samples):
            output_pdb = codesign_dir / f'{prefix}_{sample_index}.pdb'
            if not output_pdb.is_file():
                raise FileNotFoundError(
                    f'PepGLAD did not produce expected output: {output_pdb}'
                )
            generated.append((output_pdb, rank, sample_index))
    return generated


def prepare_rosetta_inputs(generated, score_dir, prefix, peptide_chain):
    score_dir.mkdir(parents=True, exist_ok=True)
    list_file = score_dir / f'{prefix}_pdb_file_list.txt'

    with list_file.open('w', encoding='utf-8') as handle:
        for output_pdb, rank, sample_index in generated:
            base = score_dir / f'{prefix}_allosteric_site_{rank}_sample_{sample_index}'
            complex_pdb = Path(f'{base}.pdb')
            shutil.copy(output_pdb, complex_pdb)

            if count_clashes.count_clashes(complex_pdb) > 0:
                print(f'{prefix}_{sample_index}', 'Steric Clashes')
                continue

            split_protein_peptide(complex_pdb, peptide_chain)
            complex_pdb.rename(f'{base}_complex.pdb')
            for structure in ('complex', 'protein', 'peptide'):
                handle.write(f'{base}_{structure}.pdb\n')
    return list_file


def read_interface_scores(score_file):
    scores = pd.read_table(
        score_file, skiprows=1, sep=r'\s+',
        usecols=['description', 'total_score'],
    )
    scores['description'] = scores['description'].str.slice(stop=-5)
    scores[['prefix', 'site', 'peptide', 'structure']] = (
        scores['description'].str.extract(
            r'(\S+)_allosteric_site_(\d+)_sample_(\d+)_([a-zA-Z]+)'
        )
    )
    scores = scores[scores['total_score'] > 0]

    by_structure = {
        structure: scores[scores['structure'] == structure]
        .set_index(['prefix', 'site', 'peptide'])[['total_score']]
        for structure in ('complex', 'protein', 'peptide')
    }
    interface = (
        by_structure['protein']
        + by_structure['peptide']
        - by_structure['complex']
    )
    return interface.dropna()


def main(argv=None):
    args = parse_args(argv)
    prefix = Path(args.pdb_file).stem
    pdb_file = os.path.abspath(args.pdb_file)
    out_dir = Path('output') / f'{prefix}_output'

    run_gaps.run_gaps(pdb_file, out_dir=str(out_dir))
    patches = run_gaps.find_high_bfactor_spatial_patches(
        out_dir / f'{prefix}_GAPS_output.pdb'
    )

    pocket_scores = pd.DataFrame([
        {
            'pocket': patch.as_pepglad_pocket(),
            'mean_gaps_score': patch.mean_bfactor,
            'max_gaps_score': patch.max_bfactor,
        }
        for patch in patches
    ])
    print(pocket_scores)
    pocket_scores.to_csv(out_dir / 'pocket_scores.csv', index=False)

    generated = run_pepglad(args, pdb_file, out_dir, prefix, patches)
    score_dir = out_dir / f'{prefix}_rosetta_score'
    list_file = prepare_rosetta_inputs(
        generated, score_dir, prefix, args.peptide_chain
    )

    score_file = score_dir / f'{prefix}_rosetta_scores.tsv'
    score_file.unlink(missing_ok=True)
    subprocess.run([
        'score_jd2', '-in:file:l', str(list_file),
        '-score:weights', 'rosetta_pdbbind_interface_regression.wts',
        '-out:file:scorefile', str(score_file),
    ], check=True)

    interface_scores = read_interface_scores(score_file)
    interface_scores.to_csv(
        score_dir / f'{prefix}_rosetta_interface_scores.csv'
    )
    print(interface_scores)


if __name__ == '__main__':
    main()
