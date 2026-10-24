"""Copy selected valid complexes and collect their user-facing results."""

import csv
import json
import math
import shutil

from . import gaps
from .fileio import atomic_csv, atomic_path


FIELDS = (
    'site', 'sample', 'site_residues', 'site_residue_count',
    'gaps_pocket_score', 'pocket_sasa_a2', 'peptide_sequence',
    'peptide_length', 'peptide_chain', 'receptor_chains',
    'source_path', 'complex_path', 'validation_status',
    'vina_score', 'vina_units', 'vina_rank', 'vina_status', 'vina_error',
    'vina_protocol_id', 'buried_sasa_a2', 'receptor_contact_residues_5a',
    'minimum_interface_distance_a',
)


def _csv_rows(path):
    with path.open(newline='', encoding='utf-8') as handle:
        return list(csv.DictReader(handle))


def export_peptides(designs, layout, input_pdb, *, scored=False,
                    scoring_error=''):
    """Publish selected complexes, pocket properties, and score outcomes.

    Pocket SASA sums residue areas measured in the original receptor.
    Only this helper's copied PDBs are refreshed; source PDBs stay put.
    """
    layout = layout.resolved()
    layout.valid_complexes_dir.mkdir(parents=True, exist_ok=True)
    rows, current, pockets = [], set(), {}
    validation, scores, ranks = {}, {}, {}
    if designs:
        validation = {row['path']: row for row in _csv_rows(
            layout.validation_reports(layout.score_dir)[0])}
        residues = {residue.key: residue for residue in
                    gaps.load_scored_residues(input_pdb, layout.gaps_scores)}
        names = {(r['chain'], r['number'], r['icode']): r['resname']
                 for r in gaps.load_score_records(layout.gaps_scores)}
    if scored:
        if layout.vina_results.exists():
            scores = {row['id']: row for row in _csv_rows(layout.vina_results)}
        if layout.vina_ranked.exists():
            ranks = {row['id']: row['rank']
                     for row in _csv_rows(layout.vina_ranked)}
    for design in sorted(designs, key=lambda d: (d.site, d.sample)):
        checked = validation[design.resolved_path]
        if checked['status'] != 'valid':
            raise ValueError(
                f'Cannot export unvalidated complex: {design.path}')
        if design.site not in pockets:
            entries = json.loads(layout.pocket_file(design.site).read_text())
            keys = sorted((chain, number, icode)
                          for chain, (number, icode) in entries)
            pockets[design.site] = {
                'site_residues': ';'.join(
                    f'{chain}-{names[chain, number, icode]}-'
                    f'{number}{icode.strip()}'
                    for chain, number, icode in keys),
                'site_residue_count': len(keys),
                'gaps_pocket_score': math.fsum(
                    residues[key].score for key in keys) / len(keys),
                'pocket_sasa_a2': math.fsum(
                    residues[key].sasa for key in keys),
            }
        destination = layout.valid_complex_file(design.site, design.sample)
        with atomic_path(destination) as temporary:
            shutil.copy2(design.path, temporary)
        current.add(destination)
        score = scores.get(design.job_id, {})
        status = score.get('status', 'error' if scored else 'not_run')
        rows.append({
            'site': design.site, 'sample': design.sample,
            **pockets[design.site],
            'peptide_sequence': design.sequence,
            'peptide_length': len(design.sequence),
            'peptide_chain': design.peptide_chain,
            'receptor_chains': ';'.join(design.receptor_chains),
            'source_path': str(design.path.resolve().relative_to(layout.root)),
            'complex_path': str(destination.relative_to(layout.root)),
            'validation_status': checked['status'],
            'vina_score': (score.get('score_total', '')
                           if status == 'ok' else ''),
            'vina_units': score.get('score_units', ''),
            'vina_rank': (ranks.get(design.job_id, '')
                          if status == 'ok' else ''),
            'vina_status': status,
            'vina_error': score.get('error', scoring_error or
                                    'Vina scoring produced no result')
            if status == 'error' else '',
            'vina_protocol_id': score.get('protocol_id', ''),
            **{name: checked.get(name, '') for name in (
                'buried_sasa_a2', 'receptor_contact_residues_5a',
                'minimum_interface_distance_a')},
        })
    # Remove stale copies when a smaller selection is requested.
    for path in layout.valid_complexes_dir.glob('site_*/*.pdb'):
        if path not in current:
            path.unlink()
    atomic_csv(layout.peptides_csv, FIELDS, rows)
    return layout.peptides_csv
