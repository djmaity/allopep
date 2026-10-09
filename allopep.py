#!/usr/bin/env python3
"""AlloPep: allosteric peptide prediction pipeline."""

import argparse
import csv
import json
import os
import re
import shutil
import subprocess
from pathlib import Path

import pandas as pd

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
    parser.add_argument('--validate-existing', action='store_true',
                        help='Check existing PepGLAD PDBs without rerunning design or scoring')
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


def run_openstructure_validation(generated, score_dir, prefix, peptide_chain):
    """Check all generated complexes and return paths that pass validation."""
    ost = os.environ.get('ALLOPEP_OST', 'ost')
    if shutil.which(ost) is None:
        raise FileNotFoundError(
            'OpenStructure ost executable not found. Set ALLOPEP_OST to the '
            'executable in a separate OpenStructure environment.'
        )
    manifest = score_dir / f'{prefix}_pepglad_validation_manifest.json'
    csv_file = score_dir / f'{prefix}_pepglad_validation.csv'
    jsonl_file = score_dir / f'{prefix}_pepglad_validation.jsonl'
    manifest.write_text(json.dumps([
        {'path': str(path.resolve()), 'site': rank, 'sample': sample}
        for path, rank, sample in generated
    ], indent=2) + '\n', encoding='utf-8')
    command = [ost, str(Path(__file__).with_name('pepglad_validation.py')),
               '--manifest', str(manifest.resolve()), '--csv', str(csv_file.resolve()),
               '--jsonl', str(jsonl_file.resolve())]
    if peptide_chain:
        command.extend(['--peptide-chain', peptide_chain])
    subprocess.run(command, check=True)
    with csv_file.open(newline='', encoding='utf-8') as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != len(generated):
        raise RuntimeError('OpenStructure result count differs from PepGLAD output count')
    return {row['path'] for row in rows if row['status'] == 'valid'}


def prepare_rosetta_inputs(generated, score_dir, prefix, peptide_chain):
    score_dir.mkdir(parents=True, exist_ok=True)
    list_file = score_dir / f'{prefix}_pdb_file_list.txt'
    valid = run_openstructure_validation(generated, score_dir, prefix, peptide_chain)

    with list_file.open('w', encoding='utf-8') as handle:
        for output_pdb, rank, sample_index in generated:
            if str(output_pdb.resolve()) not in valid:
                continue
            base = score_dir / f'{prefix}_allosteric_site_{rank}_sample_{sample_index}'
            complex_pdb = Path(f'{base}.pdb')
            shutil.copy(output_pdb, complex_pdb)
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


def find_existing_pepglad_outputs(out_dir, prefix):
    generated = []
    pattern = re.compile(rf'^{re.escape(prefix)}_allosteric_site_(\d+)_codesign$')
    for site_dir in (out_dir / f'{prefix}_PepGLAD_outputs').glob('*_codesign'):
        site = pattern.fullmatch(site_dir.name)
        if not site:
            continue
        paths = {pdb.stem: pdb for pdb in site_dir.glob('*.pdb')}
        summary = site_dir / 'summary.jsonl'
        if summary.is_file():
            expected = {json.loads(line)['id'] for line in summary.read_text(
                encoding='utf-8').splitlines() if line.strip()}
            if expected != set(paths):
                raise ValueError(f'PepGLAD PDBs do not match {summary}: '
                                 f'missing={sorted(expected - set(paths))}, '
                                 f'extra={sorted(set(paths) - expected)}')
        for stem, pdb in paths.items():
            sample = re.fullmatch(rf'{re.escape(prefix)}_(\d+)', stem)
            if sample:
                generated.append((pdb, int(site.group(1)), int(sample.group(1))))
    return sorted(generated, key=lambda item: (item[1], item[2]))


def main(argv=None):
    args = parse_args(argv)
    prefix = Path(args.pdb_file).stem
    pdb_file = os.path.abspath(args.pdb_file)
    out_dir = Path('output') / f'{prefix}_output'
    if args.validate_existing:
        generated = find_existing_pepglad_outputs(out_dir, prefix)
        if not generated:
            raise FileNotFoundError(f'No PepGLAD PDBs found under {out_dir}')
        score_dir = out_dir / f'{prefix}_rosetta_score'
        score_dir.mkdir(parents=True, exist_ok=True)
        run_openstructure_validation(generated, score_dir, prefix, args.peptide_chain)
        return

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
    if not list_file.read_text(encoding='utf-8').strip():
        print('No PepGLAD structures passed OpenStructure validation; '
              f'see {score_dir / (prefix + "_pepglad_validation.csv")}')
        return

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
