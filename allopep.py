#!/usr/bin/env python3
"""AlloPep: allosteric peptide prediction pipeline."""

import argparse
import csv
import filecmp
import json
import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

import run_gaps


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
                        help='Known peptide chain ID; uses PepGLAD summary if omitted')
    existing = parser.add_mutually_exclusive_group()
    existing.add_argument('--validate-existing', action='store_true',
                          help='Check existing PepGLAD PDBs without rerunning design or scoring')
    existing.add_argument('--rank-existing', action='store_true',
                          help='Validate and rank existing PepGLAD PDBs with Vina')
    parser.add_argument('--vina-box-padding', default=5.0, type=float,
                        help='Vina box padding on each side of the peptide in angstroms (default: 5)')
    parser.add_argument('--vina-hydrogen-policy', choices=('keep', 'rebuild'),
                        default='keep', help='Meeko hydrogen policy (default: keep)')
    return parser.parse_args(argv)


def prepare_input_pdb(input_pdb, out_dir, validate_existing=False):
    """Keep the input PDB with its output and reuse it for later validation."""
    input_pdb = Path(input_pdb)
    saved_pdb = out_dir / input_pdb.name
    if validate_existing and saved_pdb.is_file():
        return str(saved_pdb.resolve())

    source = input_pdb.resolve()
    if not source.is_file():
        raise FileNotFoundError(source)
    if not validate_existing:
        out_dir.mkdir(parents=True, exist_ok=True)
        if source != saved_pdb.resolve():
            if saved_pdb.is_file():
                if not filecmp.cmp(source, saved_pdb, shallow=False):
                    raise ValueError(
                        f'{out_dir} already contains a different input PDB: '
                        f'{saved_pdb}'
                    )
            else:
                shutil.copy2(source, saved_pdb)
        return str(saved_pdb.resolve())
    return str(source)


def run_pepglad(args, pdb_file, out_dir, prefix, patches):
    sites_dir = out_dir / f'{prefix}_allosteric_sites'
    sites_dir.mkdir(parents=True, exist_ok=True)

    env = os.environ.copy()
    env['RAY_ACCEL_ENV_VAR_OVERRIDE_ON_ZERO'] = '0'
    # PepGLAD pins OpenMM 8.0 and uses PDBFixer 1.9, which still imports
    # pkg_resources. Limit this upstream warning filter to PepGLAD workers.
    pdbfixer_warning = (
        'ignore:pkg_resources is deprecated as an API:'
        'UserWarning:pdbfixer.pdbfixer'
    )
    env['PYTHONWARNINGS'] = ','.join(filter(None, (
        env.get('PYTHONWARNINGS'), pdbfixer_warning,
    )))

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
            'conda', 'run', '--no-capture-output', '-n', 'PepGLAD',
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


def run_openstructure_validation(generated, score_dir, prefix, peptide_chain, input_pdb):
    """Check all generated complexes and return paths that pass validation."""
    ost = shutil.which('ost')
    if ost is None:
        raise FileNotFoundError(
            'OpenStructure ost executable not found in the AlloPep environment.'
        )
    manifest = score_dir / f'{prefix}_pepglad_validation_manifest.json'
    csv_file = score_dir / f'{prefix}_pepglad_validation.csv'
    jsonl_file = score_dir / f'{prefix}_pepglad_validation.jsonl'
    summaries = {}
    entries = []
    for path, rank, sample in generated:
        summary_path = path.parent / 'summary.jsonl'
        if summary_path not in summaries:
            records = [json.loads(line) for line in summary_path.read_text(
                encoding='utf-8').splitlines() if line.strip()]
            summaries[summary_path] = {record['id']: record for record in records}
            if len(summaries[summary_path]) != len(records):
                raise ValueError(f'Duplicate PepGLAD IDs in {summary_path}')
        if path.stem not in summaries[summary_path]:
            raise ValueError(f'{path.stem} is missing from {summary_path}')
        record = summaries[summary_path][path.stem]
        entries.append({
            'path': str(path.resolve()), 'site': rank, 'sample': sample,
            'input_pdb': str(Path(input_pdb).resolve()),
            'pep_seq': record['pep_seq'], 'pep_chain': record['pep_chain'],
            'rec_chains': record['rec_chains'],
        })
    manifest.write_text(json.dumps(entries, indent=2) + '\n', encoding='utf-8')
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
    errors = [row for row in rows if row['status'] == 'error']
    if errors:
        raise RuntimeError(f'{len(errors)} OpenStructure checks failed; see {csv_file}')
    return {row['path']: row['peptide_chain'] for row in rows
            if row['status'] == 'valid'}


def rank_pepglad_outputs(generated, valid_chains, score_dir, prefix,
                         box_padding=5.0, hydrogen_policy='keep'):
    """Run the bundled rigid-pose scorer and rank successful total scores."""
    score_file = score_dir / f'{prefix}_vina_ranked.csv'
    score_file.unlink(missing_ok=True)
    manifest_file = score_dir / f'{prefix}_vina_manifest.csv'
    result_file = score_dir / f'{prefix}_vina_results.csv'
    result_file.unlink(missing_ok=True)
    temporary_results = score_dir / f'{prefix}_vina_results.csv.tmp'
    temporary_results.unlink(missing_ok=True)
    samples = {}
    for output_pdb, site, sample in generated:
        path = str(output_pdb.resolve())
        if path not in valid_chains:
            continue
        samples[f'site_{site}_sample_{sample}'] = {
            'site': site, 'sample': sample, 'path': path,
            'peptide_chain': valid_chains[path],
        }
    if not samples:
        raise ValueError('No validated PepGLAD structures to score')

    with manifest_file.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=('id', 'complex_pdb',
                                                     'peptide_chain'))
        writer.writeheader()
        for job_id, sample in samples.items():
            writer.writerow({'id': job_id, 'complex_pdb': sample['path'],
                             'peptide_chain': sample['peptide_chain']})

    scorer = Path(__file__).with_name('vina_scoring') / 'vina_score.py'
    subprocess.run([
        sys.executable, str(scorer), '--manifest', str(manifest_file),
        '--workers', '1', '--box-padding', str(box_padding),
        '--hydrogen-policy', hydrogen_policy,
        '--output', str(temporary_results),
    ], check=True)
    temporary_results.replace(result_file)

    with result_file.open(newline='', encoding='utf-8') as handle:
        results = list(csv.DictReader(handle))
    if len(results) != len(samples) or {row['id'] for row in results} != set(samples):
        raise RuntimeError(f'Vina result count or IDs differ from {manifest_file}')
    rows = []
    for result in results:
        if result['status'] == 'error':
            continue
        if result['status'] != 'ok':
            raise RuntimeError(f'Unexpected Vina status for {result["id"]}: '
                               f'{result["status"]}')
        score = float(result['score_total'])
        if not math.isfinite(score):
            raise RuntimeError(f'Non-finite Vina score for {result["id"]}')
        rows.append({
            **samples[result['id']],
            'vina_score_kcal_mol': score,
            'vina_ligand_inter_kcal_mol': result['score_ligand_inter'],
            'vina_box_padding_a': box_padding,
            'vina_hydrogen_policy': hydrogen_policy,
            'vina_version': result['vina_version'],
            'meeko_version': result['meeko_version'],
        })
    if not rows:
        raise RuntimeError(f'All {len(samples)} Vina scoring jobs failed; '
                           f'see {result_file}')
    rows.sort(key=lambda row: (row['vina_score_kcal_mol'],
                               row['site'], row['sample']))
    for rank, row in enumerate(rows, start=1):
        row['rank'] = rank

    with score_file.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=(
            'rank', 'site', 'sample', 'path', 'peptide_chain',
            'vina_score_kcal_mol', 'vina_ligand_inter_kcal_mol',
            'vina_box_padding_a', 'vina_hydrogen_policy',
            'vina_version', 'meeko_version'))
        writer.writeheader()
        writer.writerows(rows)
    return rows


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


def write_pymol_script(generated, input_pdb, out_dir, prefix):
    """Load every PepGLAD complex with its receptor aligned to the input PDB."""
    summaries = {}
    models = []
    for path, site, sample in generated:
        summary_path = path.parent / 'summary.jsonl'
        if summary_path not in summaries:
            records = [json.loads(line) for line in summary_path.read_text(
                encoding='utf-8').splitlines() if line.strip()]
            summaries[summary_path] = {record['id']: record for record in records}
        record = summaries[summary_path][path.stem]
        chains = record['rec_chains']
        peptide_chain = record['pep_chain']
        if (not chains or any(not re.fullmatch(r'[A-Za-z0-9]', chain)
                              for chain in [*chains, peptide_chain])):
            raise ValueError(f'Unsupported PepGLAD chain ID in {summary_path}')
        models.append((str(path.resolve()), site, sample, chains, peptide_chain))

    object_prefix = 'allopep_' + re.sub(r'[^A-Za-z0-9_]', '_', prefix)
    reference = f'{object_prefix}_input'
    script_path = out_dir / f'{prefix}_pepglad_aligned.pml'
    lines = [
        '# Generated by AlloPep for PyMOL',
        f'load {json.dumps(str(Path(input_pdb).resolve()))}, {reference}',
    ]
    sites = {}
    for path, site, sample, chains, peptide_chain in models:
        name = f'{object_prefix}_site_{site}_sample_{sample}'
        receptor = ' or '.join(f'chain {chain}' for chain in chains)
        mobile = f'({name} and name CA and ({receptor}))'
        target = f'({reference} and name CA and ({receptor}))'
        lines.extend([
            f'load {json.dumps(path)}, {name}',
            f'align {mobile}, {target}, cycles=0',
        ])
        sites.setdefault(site, []).append((name, peptide_chain))
    lines.extend([
        f'hide everything, {object_prefix}_*',
        f'show cartoon, {reference}',
        f'color gray, {reference}',
    ])
    colors = ('yellow', 'cyan', 'magenta', 'orange', 'green', 'salmon')
    for site, site_models in sorted(sites.items()):
        pattern = f'{object_prefix}_site_{site}_*'
        lines.append(f'group {object_prefix}_site_{site}, {pattern}')
        peptide_chains = {chain for _, chain in site_models}
        if len(peptide_chains) == 1:
            selections = [f'{pattern} and chain {peptide_chains.pop()}']
        else:
            selections = [f'{name} and chain {chain}'
                          for name, chain in site_models]
        for selection in selections:
            lines.extend([
                f'show sticks, {selection}',
                f'color {colors[(site - 1) % len(colors)]}, {selection}',
            ])
    lines.extend([f'orient {reference}', ''])
    script_path.write_text('\n'.join(lines), encoding='utf-8')
    return script_path


def main(argv=None):
    args = parse_args(argv)
    if not math.isfinite(args.vina_box_padding) or args.vina_box_padding <= 0:
        raise ValueError('--vina-box-padding must be finite and positive')
    prefix = Path(args.pdb_file).stem
    out_dir = Path('output') / f'{prefix}_output'
    reuse_existing = args.validate_existing or args.rank_existing
    pdb_file = prepare_input_pdb(args.pdb_file, out_dir, reuse_existing)
    score_dir = out_dir / f'{prefix}_vina_score'
    if reuse_existing:
        generated = find_existing_pepglad_outputs(out_dir, prefix)
        if not generated:
            raise FileNotFoundError(f'No PepGLAD PDBs found under {out_dir}')
        score_dir.mkdir(parents=True, exist_ok=True)
    else:
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
        score_dir.mkdir(parents=True, exist_ok=True)

    pymol_script = write_pymol_script(generated, pdb_file, out_dir, prefix)
    print(f'PyMOL session script: {pymol_script}')
    valid_chains = run_openstructure_validation(
        generated, score_dir, prefix, args.peptide_chain, pdb_file)
    if args.validate_existing:
        return
    if not valid_chains:
        (score_dir / f'{prefix}_vina_ranked.csv').unlink(missing_ok=True)
        (score_dir / f'{prefix}_vina_results.csv').unlink(missing_ok=True)
        (score_dir / f'{prefix}_vina_manifest.csv').unlink(missing_ok=True)
        print('No PepGLAD structures passed OpenStructure validation; '
              f'see {score_dir / (prefix + "_pepglad_validation.csv")}')
        return

    ranked = rank_pepglad_outputs(
        generated, valid_chains, score_dir, prefix,
        box_padding=args.vina_box_padding,
        hydrogen_policy=args.vina_hydrogen_policy)
    failures = len(valid_chains) - len(ranked)
    print(f'Ranked {len(ranked)} PepGLAD structures by Vina score '
          f'({failures} scoring failures); '
          f'see {score_dir / (prefix + "_vina_ranked.csv")} and '
          f'{score_dir / (prefix + "_vina_results.csv")}')


if __name__ == '__main__':
    main()
