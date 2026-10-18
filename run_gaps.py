"""Run GAPS and group high-scoring surface residues into connected pockets."""

import json
import math
import subprocess
from dataclasses import dataclass
from pathlib import Path

import freesasa
import numpy as np
from Bio.PDB import PDBParser


_MINIMUM_GAPS_SCORE = 0.5
_MINIMUM_POCKET_SASA = 0.0
_POCKET_CONNECTION_RADIUS = 8.0


@dataclass(frozen=True)
class ScoredResidue:
    chain: str
    number: int
    icode: str
    score: float
    ca: tuple[float, float, float]
    sasa: float  # Absolute solvent-accessible area in Å².

    @property
    def key(self):
        return self.chain, self.number, self.icode

    def as_pepglad_entry(self):
        return [self.chain, [self.number, self.icode]]


@dataclass(frozen=True)
class GuidanceContext:
    """A connected pocket and its highest-scoring member as a reference."""

    seed: ScoredResidue
    residues: tuple[ScoredResidue, ...]

    def as_pepglad_pocket(self):
        return [residue.as_pepglad_entry() for residue in self.residues]


def _residue_sasa(residue, areas):
    chain = residue.get_parent().id
    _, number, icode = residue.id
    try:
        value = areas[chain][f'{number}{icode.strip()}'].total
    except KeyError:
        return 0.0
    return value if math.isfinite(value) else 0.0


def _read_score_records(scores_json):
    payload = json.loads(Path(scores_json).read_text(encoding='utf-8'))
    if (not isinstance(payload, dict) or type(payload.get('schema_version')) is not int
            or payload['schema_version'] != 1
            or not isinstance(payload.get('residues'), list)):
        raise ValueError(f'Invalid GAPS score sidecar: {scores_json}')
    return payload['residues']


def _scored_residue(record, model, areas, seen):
    """Verify one sidecar record against the unmodified input structure."""
    try:
        chain = record['chain']
        number = record['number']
        icode = record['icode']
        score = record['score']
        resname = record['resname']
        ca = record['ca']
    except (KeyError, TypeError) as exc:
        raise ValueError(f'Malformed GAPS score record: {record}') from exc
    if (not isinstance(chain, str) or len(chain) != 1
            or type(number) is not int or not isinstance(icode, str)
            or len(icode) != 1 or not isinstance(resname, str)
            or type(score) not in (float, int)
            or not math.isfinite(score) or not 0 <= score <= 1):
        raise ValueError(f'Invalid GAPS residue identity or score: {record}')
    key = (chain, number, icode)
    if key in seen:
        raise ValueError(f'Duplicate GAPS residue: {key}')
    seen.add(key)
    try:
        residue = model[chain][(' ', number, icode)]
    except KeyError as exc:
        raise ValueError(f'GAPS residue missing from input PDB: {key}') from exc
    if residue.resname != resname:
        raise ValueError(f'GAPS residue type differs from input PDB: {key}')
    if ca is None:
        return None
    if 'CA' not in residue or not isinstance(ca, list) or len(ca) != 3:
        raise ValueError(f'Invalid C-alpha mapping for GAPS residue: {key}')
    if not all(type(x) in (int, float) and math.isfinite(x) for x in ca):
        raise ValueError(f'Invalid C-alpha coordinates for GAPS residue: {key}')
    if np.linalg.norm(residue['CA'].coord - np.asarray(ca)) > 0.05:
        raise ValueError(f'GAPS C-alpha differs from input PDB: {key}')
    return ScoredResidue(
        chain, number, icode, float(score), tuple(map(float, ca)),
        _residue_sasa(residue, areas),
    )


def load_scored_residues(input_pdb, scores_json):
    """Check full-precision GAPS records against the original receptor PDB."""
    input_pdb = Path(input_pdb)
    records = _read_score_records(scores_json)

    structure = PDBParser(QUIET=True).get_structure(input_pdb.stem, str(input_pdb))
    model = next(structure.get_models(), None)
    if model is None:
        raise ValueError(f'No model found in {input_pdb}')
    areas = freesasa.calcBioPDB(structure)[0].residueAreas()

    scored = []
    seen = set()
    for record in records:
        residue = _scored_residue(record, model, areas, seen)
        if residue is not None:
            scored.append(residue)
    if not scored:
        raise ValueError(f'No C-alpha-bearing GAPS residues in {scores_json}')
    return scored


def select_contexts(residues, max_rank='all'):
    """Rank groups with GAPS score >0.5, absolute SASA >0 and C-alpha links <=8 Å."""
    if max_rank != 'all' and (type(max_rank) is not int or max_rank < 1):
        raise ValueError('max_rank must be a positive integer or "all"')
    eligible = sorted(
        (r for r in residues if r.score > _MINIMUM_GAPS_SCORE
         and r.sasa > _MINIMUM_POCKET_SASA),
        key=lambda r: r.key,
    )
    coords = np.asarray([r.ca for r in eligible], dtype=float)
    visited = np.zeros(len(eligible), dtype=bool)
    contexts = []
    for start in range(len(eligible)):
        if visited[start]:
            continue
        visited[start] = True
        frontier = [start]
        members = []
        while frontier:
            index = frontier.pop()
            members.append(eligible[index])
            distances = np.linalg.norm(coords - coords[index], axis=1)
            neighbors = np.flatnonzero(
                ~visited & (distances <= _POCKET_CONNECTION_RADIUS))
            visited[neighbors] = True
            frontier.extend(neighbors.tolist())
        if len(members) < 2:
            continue
        members.sort(key=lambda r: r.key)
        seed = min(members, key=lambda r: (-r.score, r.key))
        contexts.append(GuidanceContext(seed, tuple(members)))
    if not contexts:
        raise ValueError('No usable GAPS-guided receptor contexts found')
    contexts.sort(key=lambda ctx: (
        -math.fsum(r.score for r in ctx.residues) / len(ctx.residues),
        ctx.residues[0].key,
    ))
    return contexts if max_rank == 'all' else contexts[:max_rank]


def select_guidance_contexts(input_pdb, scores_json, max_rank='all'):
    """Load full-precision scores and select connected surface pockets."""
    return select_contexts(load_scored_residues(input_pdb, scores_json), max_rank)


def run_gaps(pdb_filepath, out_dir='.', pytorch_device='cuda'):
    """Run GAPS prediction in its dedicated Conda environment."""
    subprocess.run([
        'conda', 'run', '--no-capture-output', '-n', 'gaps', 'python',
        str(Path(__file__).with_name('gaps_predict.py')),
        '--pdb', str(Path(pdb_filepath).resolve()),
        '--out-dir', str(Path(out_dir).resolve()),
        '--device', pytorch_device,
    ], check=True)
