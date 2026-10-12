"""Read PepGLAD design metadata once for the downstream pipeline stages."""

import json
import re
from dataclasses import dataclass
from pathlib import Path


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


def _summary_records(summary_path):
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


def load_design_outputs(generated, *, exact=False, summary_paths=()):
    """Pair output PDBs with their summary records, caching each summary."""
    generated = [(Path(path), site, sample) for path, site, sample in generated]
    summaries = {}
    for summary_path in [*summary_paths, *(path.parent / 'summary.jsonl' for path, _, _ in generated)]:
        summary_path = Path(summary_path)
        if summary_path not in summaries:
            summaries[summary_path] = _summary_records(summary_path)
    if exact:
        for summary_path, records in summaries.items():
            actual = {path.stem for path in summary_path.parent.glob('*.pdb')}
            expected = {identifier for identifier, record in records.items()
                        if record.get('status') != 'invalid'}
            if actual != expected:
                raise ValueError(f'PepGLAD PDBs do not match {summary_path}: '
                                 f'missing={sorted(expected - actual)}, '
                                 f'extra={sorted(actual - expected)}')

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
        if record.get('status') == 'invalid':
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
