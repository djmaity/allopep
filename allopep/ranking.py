"""Score physically valid PepGLAD complexes with AlloPep's Vina scorer."""

import importlib
import math
from pathlib import Path

from .fileio import atomic_csv, sha256

IDENTITY_FIELDS = ('site', 'sample', 'path')


def check_vina_dependencies():
    """Check scoring dependencies before spending time on generation."""
    try:
        for module in ('gemmi', 'meeko', 'psutil', 'rdkit', 'scipy', 'vina'):
            importlib.import_module(module)
        from .scoring.preparation import find_meeko_receptor_preparer
        find_meeko_receptor_preparer()
    except (ImportError, OSError, RuntimeError) as exc:
        raise RuntimeError(
            'Vina scoring dependencies are unavailable; update the allopep '
            'environment with conda env update -n allopep -f environment.yml: '
            f'{exc}') from exc


def _check_chains(design):
    """Require the coordinate chains recorded in the PepGLAD summary."""
    with design.path.open(encoding='utf-8', errors='replace') as handle:
        actual = {line[21] for line in handle
                  if line.startswith(('ATOM  ', 'HETATM')) and len(line) >= 22}
    expected = set(design.receptor_chains) | {design.peptide_chain}
    if actual != expected:
        raise ValueError(f'PDB chains {sorted(actual)} differ from '
                         f'PepGLAD chains {sorted(expected)}')


def rank_pepglad_outputs(valid_designs, layout):
    """Score selected valid designs and rank finite Vina totals ascending.

    Each pose runs in a disposable worker through the shared scorer. A
    scoring failure never changes the physical-validity or sampling policy.
    """
    from .scoring.service import BatchRunner
    from .scoring.models import CombinedStructure, ScoreJob, WorkflowConfig
    from .scoring.records import (
        PROTOCOL_ID, RESULT_FIELDS, failure, result_record,
    )

    layout.score_dir.mkdir(parents=True, exist_ok=True)
    manifest_file = layout.vina_manifest
    result_file = layout.vina_results
    ranked_file = layout.vina_ranked
    for report in (manifest_file, result_file, ranked_file):
        report.unlink(missing_ok=True)
    if not valid_designs:
        raise ValueError('No validated PepGLAD structures to score')
    identifiers = [design.job_id for design in valid_designs]
    if len(set(identifiers)) != len(identifiers):
        raise ValueError('Duplicate PepGLAD job ID')

    config = WorkflowConfig()
    jobs, preflight_failures, manifest = [], [], []
    for index, design in enumerate(valid_designs):
        job = ScoreJob(index, design.job_id, CombinedStructure(
            Path(design.resolved_path), design.peptide_chain,
            tuple(design.receptor_chains)))
        digest = ''
        try:
            _check_chains(design)
            digest = sha256(design.path)
        except (OSError, ValueError) as exc:
            preflight_failures.append(failure(job, config, exc))
        else:
            jobs.append(job)
        manifest.append({
            'id': job.job_id, 'site': design.site, 'sample': design.sample,
            'path': design.resolved_path, 'complex_pdb': design.resolved_path,
            'peptide_chain': design.peptide_chain,
            'receptor_chains': ';'.join(design.receptor_chains),
            'complex_sha256': digest, 'protocol_id': PROTOCOL_ID,
            'hydrogen_policy': config.hydrogen_policy, 'ph': config.ph,
            'box_padding': config.box_padding, 'spacing': config.spacing,
        })
    atomic_csv(manifest_file, tuple(manifest[0]), manifest)

    outcomes = preflight_failures
    if jobs:
        outcomes += BatchRunner(config).run(jobs)
    designs = {design.job_id: design for design in valid_designs}
    records = []
    for outcome in sorted(outcomes, key=lambda item: item.job.index):
        row = result_record(outcome)
        design = designs[row['id']]
        row.update(site=design.site, sample=design.sample,
                   path=design.resolved_path)
        if row['status'] == 'ok':
            try:
                score = float(row['score_total'])
                if not math.isfinite(score):
                    raise ValueError('non-finite score_total')
            except (TypeError, ValueError) as exc:
                row.update(status='error', error=f'Invalid Vina total: {exc}')
            else:
                row['score_total'] = score
        records.append(row)
    fields = (*IDENTITY_FIELDS, *RESULT_FIELDS)
    atomic_csv(result_file, fields, records)
    ranked = [row.copy() for row in records if row['status'] == 'ok']
    if not ranked:
        raise RuntimeError(f'All {len(records)} Vina scoring jobs failed; '
                           f'see {result_file}')
    ranked.sort(key=lambda row: (
        row['score_total'], row['site'], row['sample']))
    for rank, row in enumerate(ranked, 1):
        row['rank'] = rank
    atomic_csv(ranked_file, ('rank', *fields), ranked)
    return ranked
