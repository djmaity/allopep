"""Read PepGLAD design metadata once for the downstream pipeline stages."""

import csv
import json
import re
from dataclasses import dataclass
from pathlib import Path

# A design PepGLAD itself reported as unusable; it keeps no PDB.
REJECTED = ('invalid', 'error')
# OpenMM stamps every design it relaxes, so a missing stamp means the
# relaxation step never finished for that file.
RELAXED_REMARK = 'REMARK   1 CREATED WITH OPENMM'


def is_relaxed(path):
    """Whether OpenMM finished relaxing this design."""
    with Path(path).open(encoding='utf-8') as handle:
        return handle.readline().startswith(RELAXED_REMARK)


@dataclass(frozen=True)
class DesignOutput:
    path: Path
    site: int
    sample: int
    sequence: str
    peptide_chain: str
    receptor_chains: tuple[str, ...]

    @property
    def job_id(self):
        return f'site_{self.site}_sample_{self.sample}'

    @property
    def resolved_path(self):
        return str(self.path.resolve())


def summary_records(summary_path):
    """Index a PepGLAD summary by design ID, rejecting repeats."""
    records = {}
    for line_number, line in enumerate(summary_path.read_text(encoding='utf-8').splitlines(), 1):
        if not line.strip():
            continue
        record = json.loads(line)
        if not isinstance(record, dict) or not isinstance(record.get('id'), str):
            raise ValueError(f'Invalid PepGLAD entry in {summary_path}:{line_number}')
        identifier = record['id']
        if identifier in records:
            raise ValueError(f'Duplicate PepGLAD IDs in {summary_path}: {identifier}')
        records[identifier] = record
    return records


def removed_structures(directory, records=None):
    """Read the deletion journal, allowing only structures with a recorded rejection."""
    directory = Path(directory)
    journal = directory / 'removed_structures.json'
    if not journal.exists():
        return {}
    removed = json.loads(journal.read_text(encoding='utf-8'))
    if not isinstance(removed, dict):
        raise ValueError(f'Invalid removed-structure journal: {journal}')
    records = summary_records(directory / 'summary.jsonl') if records is None else records
    validation = None
    for relative, item in removed.items():
        path = Path(relative)
        if (path.is_absolute() or path.suffix != '.pdb'
                or path.parts not in ((path.name,), ('invalid_designs', path.name))
                or not (directory / path).resolve().is_relative_to(directory.resolve())
                or path.stem not in records or not isinstance(item, dict)
                or item.get('status') not in REJECTED
                or not isinstance(item.get('sha256'), str)
                or not re.fullmatch('[0-9a-f]{64}', item['sha256'])):
            raise ValueError(f'Invalid removed-structure entry in {journal}: {relative}')
        record = records[path.stem]
        if record.get('status') in REJECTED:
            status = record['status']
        else:
            if validation is None:
                reports = list(directory.glob('*_pepglad_validation.csv'))
                if len(reports) != 1:
                    raise ValueError(f'Missing rejection report for {journal}')
                with reports[0].open(newline='', encoding='utf-8') as handle:
                    validation = {row['path']: row for row in csv.DictReader(handle)}
            row = validation.get(str((directory / path).resolve()), {})
            status = row.get('status')
            if row.get('sha256') and row['sha256'] != item['sha256']:
                raise ValueError(f'Removed-structure hash differs from validation: {relative}')
        if status != item['status']:
            raise ValueError(f'Structure has no matching rejection in {journal}: {relative}')
    return removed


def removed_top_level_stems(directory, records=None):
    """Stems of rejected PDBs journalled as removed from the directory itself."""
    return {Path(relative).stem
            for relative in removed_structures(directory, records)
            if len(Path(relative).parts) == 1}


def check_summary_matches_pdbs(directory, records):
    """Require the PDBs present to be exactly the designs not rejected."""
    directory = Path(directory)
    removed = removed_top_level_stems(directory, records)
    actual = {path.stem for path in directory.glob('*.pdb')} - removed
    expected = {identifier for identifier, record in records.items()
                if record.get('status') not in REJECTED} - removed
    if actual != expected:
        raise ValueError(f'PepGLAD PDBs do not match {directory}: '
                         f'missing={sorted(expected - actual)}, '
                         f'extra={sorted(actual - expected)}')


def load_design_outputs(generated, *, exact=False, summary_paths=()):
    """Pair output PDBs with their summary records, caching each summary."""
    generated = [(Path(path), site, sample) for path, site, sample in generated]
    summaries = {}
    for summary_path in [*summary_paths, *(path.parent / 'summary.jsonl' for path, _, _ in generated)]:
        summary_path = Path(summary_path)
        if summary_path not in summaries:
            summaries[summary_path] = summary_records(summary_path)
    if exact:
        for summary_path, records in summaries.items():
            check_summary_matches_pdbs(summary_path.parent, records)

    designs = []
    seen = set()
    for path, site, sample in generated:
        summary_path = path.parent / 'summary.jsonl'
        records = summaries[summary_path]
        if path.stem not in records:
            raise ValueError(f'{path.stem} is missing from {summary_path}')
        if (site, sample) in seen:
            raise ValueError(f'Duplicate PepGLAD site/sample: {site}/{sample}')
        seen.add((site, sample))
        record = records[path.stem]
        if record.get('status') in REJECTED:
            raise ValueError(f'Invalid PepGLAD design still has a PDB: {path}')
        try:
            sequence = record['pep_seq']
            peptide_chain = record['pep_chain']
            receptor_chains = record['rec_chains']
        except KeyError as exc:
            raise ValueError(f'Incomplete PepGLAD entry {path.stem} in {summary_path}') from exc
        chains = [*receptor_chains, peptide_chain] if isinstance(receptor_chains, list) else []
        if (not isinstance(sequence, str) or not sequence
                or not isinstance(receptor_chains, list) or not receptor_chains
                or not isinstance(peptide_chain, str)
                or any(not isinstance(chain, str)
                       or not re.fullmatch(r'[A-Za-z0-9]', chain) for chain in chains)
                or len(set(chains)) != len(chains)):
            raise ValueError(f'Invalid PepGLAD chains or sequence in {summary_path}')
        designs.append(DesignOutput(
            path, site, sample, sequence, peptide_chain, tuple(receptor_chains)))
    return designs


def load_round_outputs(directory, site, prefix, start, count, *, partial=False,
                       include_removed=False):
    """Map PepGLAD's local IDs to the reserved per-site sample range."""
    directory = Path(directory)
    summary = directory / 'summary.jsonl'
    if not summary.exists() and partial and not list(directory.glob('*.pdb')):
        return [], {}
    records = summary_records(summary)
    expected = {f'{prefix}_{index}' for index in range(count)}
    if not set(records) <= expected or (not partial and set(records) != expected):
        raise ValueError(f'Saved PepGLAD sample IDs differ: {summary}')
    removed = removed_structures(directory, records)
    generated = [(directory / f'{prefix}_{index}.pdb', site, start + index)
                 for index in range(count)
                 if (records.get(f'{prefix}_{index}', {}).get('status')
                     not in REJECTED)
                 and ((directory / f'{prefix}_{index}.pdb').is_file()
                      or (include_removed and f'{prefix}_{index}.pdb' in removed))
                 and (include_removed or f'{prefix}_{index}.pdb' not in removed)]
    designs = load_design_outputs(generated, exact=True, summary_paths=[summary])
    return designs, records


def site_rounds(site_dir, prefix):
    """Discover reserved round directories without guessing sample IDs.

    Designs live only in ``round_NNNNNN`` subdirectories. PDBs or a
    summary directly in the site directory are the pre-round layout,
    which is no longer read.
    """
    site_dir = Path(site_dir)
    rounds = []
    if (site_dir / 'summary.jsonl').is_file() or list(site_dir.glob('*.pdb')):
        raise ValueError(
            f'Unsupported PepGLAD outputs directly in {site_dir}; designs '
            'belong in round_NNNNNN subdirectories')
    for directory in sorted(site_dir.glob('round_*')):
        if not directory.is_dir():
            continue
        metadata = json.loads((directory / 'round.json').read_text())
        start, count = metadata['start'], metadata['count']
        if type(start) is not int or type(count) is not int or start < 0 or count < 1:
            raise ValueError(f'Invalid reserved sample range in {directory}')
        rounds.append({'directory': directory, 'start': start, 'count': count})
    end = 0
    for entry in sorted(rounds, key=lambda entry: entry['start']):
        if entry['start'] < end:
            raise ValueError(f'Overlapping PepGLAD sample ranges in {site_dir}')
        end = entry['start'] + entry['count']
    return rounds
