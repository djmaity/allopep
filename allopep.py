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
import tempfile
from pathlib import Path

from Bio.PDB import PDBParser

from pepglad_outputs import load_design_outputs
from rosetta_scoring import rank_pepglad_outputs, resolve_rosetta_binary
import run_gaps


def _max_rank(value):
    if value == 'all':
        return value
    try:
        rank = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError('expected a positive integer or "all"') from exc
    if rank < 1:
        raise argparse.ArgumentTypeError('expected a positive integer or "all"')
    return rank


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog='python allopep.py',
        description='Pipeline for allosteric peptide design',
    )
    parser.add_argument('pdb_file')
    parser.add_argument('-r', '--max_rank', default='all', type=_max_rank,
                        help='Number of ranked GAPS pockets to design (default: "all"; or a positive integer)')
    parser.add_argument('-l', '--length_min', default=5, type=int,
                        help='Minimum length of the predicted peptide')
    parser.add_argument('-L', '--length_max', default=10, type=int,
                        help='Maximum length of the predicted peptide')
    parser.add_argument('-n', '--num_samples', default=10, type=int,
                        help='Number of predicted peptides per site')
    parser.add_argument('-c', '--peptide-chain', default=None,
                        help='Known peptide chain ID; uses PepGLAD summary if omitted')
    existing = parser.add_mutually_exclusive_group()
    existing.add_argument('-v', '--validate-existing', action='store_true',
                          help='Check existing PepGLAD PDBs without rerunning design or scoring')
    existing.add_argument('-R', '--rank-existing', action='store_true',
                          help='Validate and rank existing PepGLAD PDBs with Rosetta')
    existing.add_argument('--resume', action='store_true',
                          help='Reuse completed PepGLAD sites and design missing sites')
    parser.add_argument('--rosetta-bin', default=None,
                        help='InterfaceAnalyzer executable (default: find on PATH)')
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


def run_pepglad(args, pdb_file, out_dir, prefix, contexts):
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

    generated = []
    for rank, context in enumerate(contexts, start=1):
        pocket_file = sites_dir / f'{prefix}_allosteric_site_{rank}.json'
        pocket = context.as_pepglad_pocket()
        codesign_dir = (
            out_dir / f'{prefix}_PepGLAD_outputs'
            / f'{prefix}_allosteric_site_{rank}_codesign'
        )
        if getattr(args, 'resume', False) and codesign_dir.exists():
            if not pocket_file.is_file() or json.loads(pocket_file.read_text()) != pocket:
                raise ValueError(f'Saved PepGLAD pocket differs for site {rank}: {pocket_file}')
            summary_path = codesign_dir / 'summary.jsonl'
            records = [json.loads(line) for line in summary_path.read_text().splitlines()
                       if line.strip()]
            expected = {f'{prefix}_{i}' for i in range(args.num_samples)}
            if len(records) != args.num_samples or {r['id'] for r in records} != expected:
                raise ValueError(f'Saved PepGLAD sample IDs differ: {summary_path}')
            site_generated = [(codesign_dir / f'{prefix}_{i}.pdb', rank, i)
                              for i in range(args.num_samples)
                              if (codesign_dir / f'{prefix}_{i}.pdb').is_file()]
            designs = load_design_outputs(
                site_generated, exact=True, summary_paths=[summary_path])
            for design in designs:
                with design.path.open(encoding='utf-8') as handle:
                    if not handle.readline().startswith('REMARK   1 CREATED WITH OPENMM'):
                        raise RuntimeError(
                            f'PepGLAD site {rank} has unrelaxed samples; '
                            f'run recover_pepglad_relaxation.py on {codesign_dir}')
            generated.extend(site_generated)
            print(f'Reusing {len(designs)} relaxed designs for site {rank}')
            continue
        with pocket_file.open('w', encoding='utf-8') as handle:
            json.dump(pocket, handle)
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
            if output_pdb.is_file():
                generated.append((output_pdb, rank, sample_index))
        load_design_outputs(
            [entry for entry in generated if entry[1] == rank],
            exact=True, summary_paths=[codesign_dir / 'summary.jsonl'])
    return load_design_outputs(generated)


def run_openstructure_validation(generated, score_dir, prefix, peptide_chain, input_pdb):
    """Check generated complexes and return the designs that pass."""
    ost = shutil.which('ost')
    if ost is None:
        raise FileNotFoundError(
            'OpenStructure ost executable not found in the AlloPep environment.'
        )
    manifest = score_dir / f'{prefix}_pepglad_validation_manifest.json'
    csv_file = score_dir / f'{prefix}_pepglad_validation.csv'
    jsonl_file = score_dir / f'{prefix}_pepglad_validation.jsonl'
    entries = []
    for design in generated:
        entries.append({
            'path': design.resolved_path, 'site': design.site,
            'sample': design.sample,
            'input_pdb': str(Path(input_pdb).resolve()),
            'pep_seq': design.sequence, 'pep_chain': design.peptide_chain,
            'rec_chains': list(design.receptor_chains),
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
    by_path = {design.resolved_path: design for design in generated}
    if (len(rows) != len(generated) or len(by_path) != len(generated)
            or {row['path'] for row in rows} != set(by_path)):
        raise RuntimeError('OpenStructure result paths differ from PepGLAD outputs')
    errors = [row for row in rows if row['status'] == 'error']
    if errors:
        raise RuntimeError(f'{len(errors)} OpenStructure checks failed; see {csv_file}')
    if any(row['status'] not in ('valid', 'invalid') for row in rows):
        raise RuntimeError(f'Unexpected OpenStructure status; see {csv_file}')
    valid = []
    for row in rows:
        design = by_path[row['path']]
        if row['peptide_chain'] != design.peptide_chain:
            raise RuntimeError(f'OpenStructure peptide chain differs for {row["path"]}')
        if row['status'] == 'valid':
            valid.append(design)
    return valid


def find_existing_pepglad_outputs(out_dir, prefix):
    generated = []
    summary_paths = []
    pattern = re.compile(rf'^{re.escape(prefix)}_allosteric_site_(\d+)_codesign$')
    for site_dir in (out_dir / f'{prefix}_PepGLAD_outputs').glob('*_codesign'):
        site = pattern.fullmatch(site_dir.name)
        if not site:
            continue
        summary_paths.append(site_dir / 'summary.jsonl')
        for pdb in site_dir.glob('*.pdb'):
            stem = pdb.stem
            sample = re.fullmatch(rf'{re.escape(prefix)}_(\d+)', stem)
            if sample:
                generated.append((pdb, int(site.group(1)), int(sample.group(1))))
    generated.sort(key=lambda item: (item[1], item[2]))
    return load_design_outputs(generated, exact=True, summary_paths=summary_paths)


def _pymol_score_commands(out_dir, prefix, reference):
    """Apply saved GAPS scores to the original receptor's residue identities."""
    scores_file = out_dir / f'{prefix}_GAPS_scores.json'
    records = json.loads(scores_file.read_text(encoding='utf-8'))['residues']
    structure = PDBParser(QUIET=True).get_structure(
        prefix, str(out_dir / f'{prefix}.pdb'))
    model = next(structure.get_models())
    lines = [f'alter {reference}, b=-1', f'color gray70, {reference}']
    seen = set()
    if not records:
        raise ValueError(f'No GAPS scores found in {scores_file}')
    for record in records:
        key = record['chain'], record['number'], record['icode']
        chain, number, icode = key
        try:
            residue = model[chain][(' ', number, icode)]
        except KeyError as exc:
            raise ValueError(f'GAPS score residue missing from original PDB: {key}') from exc
        score = record['score']
        if (key in seen or residue.resname != record['resname']
                or type(score) not in (float, int)
                or not math.isfinite(score) or not 0 <= score <= 1):
            raise ValueError(f'Invalid or mismatched GAPS score for {key}')
        if record['ca'] is not None and (
                'CA' not in residue
                or math.dist(residue['CA'].coord, record['ca']) >= 0.01):
            raise ValueError(f'GAPS score coordinates differ from original PDB: {key}')
        seen.add(key)
        lines.append(
            f'alter {reference} and polymer.protein and chain {json.dumps(chain.strip())} '
            f'and resi {number}{icode.strip()}, b={score!r}')
    return lines


def _pymol_commands(generated, out_dir, prefix, *, include_all_sites=False):
    """Show the original receptor colored by GAPS score and aligned peptides."""
    reference = re.sub(r'[^A-Za-z0-9_]', '_', prefix)
    input_pdb = out_dir / f'{prefix}.pdb'
    if not input_pdb.is_file():
        raise FileNotFoundError(f'Original receptor structure not found: {input_pdb}')
    lines = [
        '# Generated by AlloPep for PyMOL',
        f'load {json.dumps(str(input_pdb.resolve()))}, {reference}',
        *_pymol_score_commands(out_dir, prefix, reference),
    ]
    sites = {}
    for design in generated:
        name = f'site_{design.site}_peptide_{design.sample + 1}'
        receptor = ' or '.join(f'chain {chain}' for chain in design.receptor_chains)
        mobile = f'({name} and name CA and ({receptor}))'
        target = f'({reference} and name CA and ({receptor}))'
        lines.extend([
            f'load {json.dumps(design.resolved_path)}, {name}',
            f'align {mobile}, {target}, cycles=0',
            f'remove {name} and not chain {design.peptide_chain}',
        ])
        sites.setdefault(design.site, []).append((name, design.peptide_chain))
    if include_all_sites:
        sites_dir = out_dir / f'{prefix}_allosteric_sites'
        pattern = re.compile(rf'^{re.escape(prefix)}_allosteric_site_([1-9]\d*)\.json$')
        for pocket_file in sites_dir.glob('*.json'):
            match = pattern.fullmatch(pocket_file.name)
            if match:
                sites.setdefault(int(match.group(1)), [])
        if not sites:
            raise FileNotFoundError(f'No allosteric site definitions found in {sites_dir}')
    lines.extend([
        'hide everything, all',
        f'show cartoon, {reference}',
        f'spectrum b, blue_white_red, {reference} and polymer.protein and b > -0.5, '
        'minimum=0, maximum=1',
    ])
    # PyMOL's default object color cycle, with site 1 set to green.
    colors = (
        'green', 'cyan', 'lightmagenta', 'yellow', 'salmon', 'hydrogen',
        'slate', 'orange', 'lime', 'deepteal', 'hotpink', 'yelloworange',
        'violetpurple', 'grey70', 'marine', 'olive', 'smudge', 'teal',
        'dirtyviolet', 'wheat', 'deepsalmon', 'lightpink', 'aquamarine',
        'paleyellow', 'limegreen', 'skyblue', 'warmpink', 'limon',
        'violet', 'bluewhite', 'greencyan', 'sand', 'forest', 'lightteal',
        'darksalmon', 'splitpea', 'raspberry', 'grey50', 'deepblue', 'brown',
    )
    for site, site_models in sorted(sites.items()):
        peptide_pattern = f'site_{site}_peptide_*'
        color = colors[(site - 1) % len(colors)]
        pocket_file = (out_dir / f'{prefix}_allosteric_sites'
                       / f'{prefix}_allosteric_site_{site}.json')
        if pocket_file.is_file():
            residues_by_chain = {}
            for entry in json.loads(pocket_file.read_text(encoding='utf-8')):
                chain, (number, icode) = entry[:2]
                residues_by_chain.setdefault(chain.strip() or '""', []).append(
                    f'{number}{icode.strip()}')
            surface = f'site_{site}'
            site_selection = ' or '.join(
                f'(chain {chain} and resi {"+".join(residues)})'
                for chain, residues in sorted(residues_by_chain.items()))
            lines.extend([
                f'create {surface}, {reference} and polymer.protein',
                f'hide everything, {surface}',
                f'show surface, {surface} and ({site_selection})',
                f'color {color}, {surface}',
                f'set transparency, 0.4, {surface}',
                f'group sites, {surface}',
            ])
        if not site_models:
            continue
        lines.append(f'group site_{site}_peptides, {peptide_pattern}')
        peptide_chains = {chain for _, chain in site_models}
        if len(peptide_chains) == 1:
            selections = [f'{peptide_pattern} and chain {peptide_chains.pop()}']
        else:
            selections = [f'{name} and chain {chain}'
                          for name, chain in site_models]
        for selection in selections:
            lines.extend([
                f'show sticks, {selection}',
                f'color {color}, {selection}',
                f'color atomic, {selection} and not elem C',
            ])
        lines.append(f'group peptides, site_{site}_peptides')
    lines.extend([f'orient {reference}', ''])
    return lines


def _save_pymol_session(commands, out_dir, stem, *, keep_script=False):
    """Build a session in isolation and publish it only after PyMOL succeeds."""
    session_path = out_dir / f'{stem}.pse'
    saved_script = out_dir / f'{stem}.pml'
    with tempfile.TemporaryDirectory(prefix='.allopep-pymol-', dir=out_dir) as temporary:
        temporary = Path(temporary)
        script_path = temporary / 'build.pml'
        temporary_session = temporary / 'session.pse'
        script_path.write_text(
            '\n'.join([*commands, 'save session.pse', '']),
            encoding='utf-8')
        try:
            result = subprocess.run(
                ['pymol', '-cqQ', script_path.name], cwd=temporary,
                capture_output=True, text=True)
        except FileNotFoundError as exc:
            raise RuntimeError('PyMOL is required to create the .pse session') from exc
        output = '\n'.join((result.stdout, result.stderr))
        errors = [line.strip() for line in output.splitlines()
                  if line.lstrip().startswith(('Error:', 'PyMOL-Error:', 'Traceback'))]
        if (result.returncode or errors or not temporary_session.is_file()
                or not temporary_session.stat().st_size):
            details = '\n'.join([*errors[-5:], *output.splitlines()[-8:]])
            raise RuntimeError(f'PyMOL session creation failed: {details}')
        if keep_script:
            script_path.write_text('\n'.join(commands), encoding='utf-8')
        os.replace(temporary_session, session_path)
        if keep_script:
            os.replace(script_path, saved_script)
    if not keep_script:
        saved_script.unlink(missing_ok=True)
    return session_path


def write_pymol_session(generated, out_dir, prefix):
    """Save a ready-to-open PyMOL session without leaving a PML script."""
    return _save_pymol_session(
        _pymol_commands(generated, out_dir, prefix), out_dir,
        f'{prefix}_pepglad_aligned')


def write_allosteric_site_visualization(out_dir, prefix):
    """Save a native PML script and session for every selected GAPS pocket."""
    out_dir = Path(out_dir)
    return _save_pymol_session(
        _pymol_commands([], out_dir, prefix, include_all_sites=True), out_dir,
        f'{prefix}_allosteric_site', keep_script=True)


def main(argv=None):
    args = parse_args(argv)
    if args.length_min < 1 or args.length_max < args.length_min:
        raise ValueError('--length_min must be positive and no greater than --length_max')
    if args.num_samples < 1:
        raise ValueError('--num_samples must be positive')
    if args.peptide_chain is not None and not re.fullmatch(r'[A-Za-z0-9]', args.peptide_chain):
        raise ValueError('--peptide-chain must be one alphanumeric chain ID')
    rosetta_bin = None if args.validate_existing else resolve_rosetta_binary(args.rosetta_bin)
    prefix = Path(args.pdb_file).stem
    out_dir = Path('output') / f'{prefix}_output'
    reuse_existing = args.validate_existing or args.rank_existing
    pdb_file = prepare_input_pdb(args.pdb_file, out_dir, reuse_existing)
    score_dir = out_dir / f'{prefix}_rosetta_score'
    if reuse_existing:
        generated = find_existing_pepglad_outputs(out_dir, prefix)
        if not generated:
            raise FileNotFoundError(f'No PepGLAD PDBs found under {out_dir}')
        score_dir.mkdir(parents=True, exist_ok=True)
    else:
        run_gaps.run_gaps(pdb_file, out_dir=str(out_dir))
        contexts = run_gaps.select_guidance_contexts(
            pdb_file, out_dir / f'{prefix}_GAPS_scores.json',
            max_rank=args.max_rank,
        )
        print(f'Designing against {len(contexts)} GAPS-guided sites')
        generated = run_pepglad(args, pdb_file, out_dir, prefix, contexts)
        score_dir.mkdir(parents=True, exist_ok=True)

    pymol_session = write_pymol_session(generated, out_dir, prefix)
    print(f'PyMOL session: {pymol_session}')
    valid_designs = run_openstructure_validation(
        generated, score_dir, prefix, args.peptide_chain, pdb_file)
    if args.validate_existing:
        return
    if not valid_designs:
        (score_dir / f'{prefix}_rosetta_ranked.csv').unlink(missing_ok=True)
        (score_dir / f'{prefix}_rosetta_results.csv').unlink(missing_ok=True)
        (score_dir / f'{prefix}_rosetta_manifest.csv').unlink(missing_ok=True)
        print('No PepGLAD structures passed OpenStructure validation; '
              f'see {score_dir / (prefix + "_pepglad_validation.csv")}')
        return

    ranked = rank_pepglad_outputs(valid_designs, score_dir, prefix, rosetta_bin)
    failures = len(valid_designs) - len(ranked)
    print(f'Ranked {len(ranked)} PepGLAD structures by Rosetta dG_separated '
          f'({failures} scoring failures); '
          f'see {score_dir / (prefix + "_rosetta_ranked.csv")} and '
          f'{score_dir / (prefix + "_rosetta_results.csv")}')


if __name__ == '__main__':
    main()
