"""Command line, stage orchestration, and validation invocation."""

import argparse
import filecmp
import os
import shutil
from pathlib import Path

from . import gaps
from .layout import REPO_ROOT, OutputLayout
from .outputs import load_round_outputs, site_rounds
from .process_control import run_logged
from .pymol_session import (write_allosteric_site_visualization,
                            write_pymol_session)
from .ranking import check_vina_dependencies, rank_pepglad_outputs
from .sampling import (SamplingConfig, SamplingController, require_resume,
                       sampling_lock)
from .validation_reports import read_rows, run_ost, write_manifest


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
        prog='python -m allopep',
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
                        help='Number of physically valid structures per site (default: 10)')
    parser.add_argument('--attempt-multiplier', default=1000, type=int,
                        help='Attempt budget per site is this positive integer times -n (default: 1000)')
    parser.add_argument('--max_num_peptides', default=None, type=int,
                        help='Absolute attempt cap per site; overrides --attempt-multiplier times -n')
    parser.add_argument('--round-timeout-seconds', default=3600, type=int,
                        help='Deadline for generation, relaxation and validation of one round (default: 3600)')
    parser.add_argument('-c', '--peptide-chain', default=None,
                        help='Known peptide chain ID; uses PepGLAD summary if omitted')
    existing = parser.add_mutually_exclusive_group()
    existing.add_argument('-v', '--validate-existing', action='store_true',
                          help='Check existing PepGLAD PDBs without rerunning design or scoring')
    existing.add_argument('-R', '--rank-existing', action='store_true',
                          help='Validate and rank existing PepGLAD PDBs with Vina')
    existing.add_argument('--resume', action='store_true',
                          help='Resume sampling using saved rounds and cumulative attempt budgets')
    return parser.parse_args(argv)


def prepare_input_pdb(input_pdb, layout, validate_existing=False):
    """Keep the input PDB with its output and reuse it for later validation."""
    input_pdb = Path(input_pdb)
    out_dir = layout.root
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


def run_pepglad_round(config, pdb_file, pocket_file, round_dir, count, *, deadline):
    """Generate one isolated round; the controller owns IDs and attempt budgets."""
    env = os.environ.copy()
    env['RAY_ACCEL_ENV_VAR_OVERRIDE_ON_ZERO'] = '0'
    pdbfixer_warning = (
        'ignore:pkg_resources is deprecated as an API:'
        'UserWarning:pdbfixer.pdbfixer'
    )
    env['PYTHONWARNINGS'] = ','.join(filter(None, (
        env.get('PYTHONWARNINGS'), pdbfixer_warning,
    )))
    run_logged([
        'conda', 'run', '--no-capture-output', '-n', 'PepGLAD',
        'python', '-m', 'api.run', '--mode', 'codesign',
        '--pdb', str(Path(pdb_file).resolve()),
        '--pocket', str(Path(pocket_file).resolve()),
        '--out_dir', str(Path(round_dir).resolve()),
        '--length_min', str(config.length_min),
        '--length_max', str(config.length_max),
        '--n_samples', str(count),
    ], cwd=REPO_ROOT / 'PepGLAD', env=env,
        log_file=round_dir / 'pepglad.log', deadline=deadline)


def run_pepglad(config, pdb_file, layout, contexts):
    with sampling_lock(layout):
        controller = SamplingController(
            config, pdb_file, layout, contexts,
            run_pepglad_round, run_openstructure_validation)
        session = write_allosteric_site_visualization(layout)
        print(f'Allosteric site visualization: {session.with_suffix(".pml")} and {session}')
        return controller.run()


def run_openstructure_validation(generated, directory, layout, peptide_chain, input_pdb,
                                *, deadline=None, allow_errors=False):
    """Check generated complexes and return the designs that pass."""
    csv_file, jsonl_file, manifest = layout.validation_reports(directory)
    write_manifest(generated, manifest, input_pdb)
    run_ost(manifest, csv_file, jsonl_file, peptide_chain, deadline=deadline,
            log_file=directory / 'validation.log')
    _, rows, by_path = read_rows(csv_file, generated)
    errors = [row for row in rows if row['status'] == 'error']
    if errors and not allow_errors:
        raise RuntimeError(f'{len(errors)} OpenStructure checks failed; see {csv_file}')
    valid = []
    for row in rows:
        design = by_path[row['path']]
        if row['status'] == 'error':
            continue
        if row['peptide_chain'] != design.peptide_chain:
            raise RuntimeError(f'OpenStructure peptide chain differs for {row["path"]}')
        if row['status'] == 'valid':
            valid.append(design)
    return valid


def find_existing_pepglad_outputs(layout):
    designs = []
    for site_dir in layout.designs_dir.glob('*_codesign'):
        site = layout.site_number(site_dir)
        if site is None:
            continue
        for entry in site_rounds(site_dir, layout.prefix):
            outputs, _ = load_round_outputs(
                entry['directory'], site, layout.prefix,
                entry['start'], entry['count'], partial=True)
            designs.extend(outputs)
    designs.sort(key=lambda design: (design.site, design.sample))
    return designs


def sampling_config(args):
    """Resolve the per-site attempt budget and validate the settings."""
    if args.attempt_multiplier < 1:
        raise ValueError('--attempt-multiplier must be positive')
    return SamplingConfig(
        length_min=args.length_min, length_max=args.length_max,
        num_samples=args.num_samples,
        max_attempts=(args.max_num_peptides if args.max_num_peptides is not None
                      else args.attempt_multiplier * args.num_samples),
        round_timeout_seconds=args.round_timeout_seconds,
        peptide_chain=args.peptide_chain, resume=args.resume,
    )


def main(argv=None):
    args = parse_args(argv)
    config = sampling_config(args)
    if not args.validate_existing:
        check_vina_dependencies()
    layout = OutputLayout.for_input(args.pdb_file)
    reuse_existing = args.validate_existing or args.rank_existing
    if not reuse_existing:
        require_resume(layout, config.resume)
    pdb_file = prepare_input_pdb(args.pdb_file, layout, reuse_existing)
    layout.score_dir.mkdir(parents=True, exist_ok=True)
    validation_csv = layout.validation_reports(layout.score_dir)[0]
    shortfall = False
    if reuse_existing:
        generated = find_existing_pepglad_outputs(layout)
        if not generated:
            raise FileNotFoundError(f'No PepGLAD PDBs found under {layout.root}')
        session = write_allosteric_site_visualization(layout)
        print(f'Allosteric site visualization: {session.with_suffix(".pml")} and {session}')
        valid_designs = run_openstructure_validation(
            generated, layout.score_dir, layout, config.peptide_chain, pdb_file)
    else:
        gaps.run_gaps(pdb_file, out_dir=str(layout.root))
        contexts = gaps.select_guidance_contexts(
            pdb_file, layout.gaps_scores, max_rank=args.max_rank,
        )
        print(f'Designing against {len(contexts)} GAPS-guided sites')
        result = run_pepglad(config, pdb_file, layout, contexts)
        valid_designs = result.selected
        shortfall = result.shortfall
        print(f'Sampling summary: {result.summary_path}')

    pymol_session = write_pymol_session(valid_designs, layout)
    print(f'PyMOL session: {pymol_session}')
    if args.validate_existing:
        return 0
    if not valid_designs:
        for report in (layout.vina_ranked, layout.vina_results,
                       layout.vina_manifest):
            report.unlink(missing_ok=True)
        print('No PepGLAD structures passed OpenStructure validation; '
              f'see {validation_csv}')
        if shortfall:
            print('Sampling stopped at --max_num_peptides; exiting with status 2.')
        return 2 if shortfall else 0

    ranked = rank_pepglad_outputs(valid_designs, layout)
    failures = len(valid_designs) - len(ranked)
    print(f'Ranked {len(ranked)} PepGLAD structures by Vina score_total '
          f'({failures} scoring failures); '
          f'see {layout.vina_ranked} and {layout.vina_results}')
    if shortfall:
        print('One or more sites exhausted --max_num_peptides before reaching '
              'the valid-structure target; partial results were retained.')
    return 2 if shortfall else 0
