"""Run GAPS and choose local receptor contexts from its residue scores."""

import json
import math
import subprocess
from dataclasses import dataclass
from pathlib import Path

import freesasa
import numpy as np
from Bio.PDB import PDBParser


@dataclass(frozen=True)
class ScoredResidue:
    chain: str
    number: int
    icode: str
    score: float
    ca: tuple[float, float, float]
    relative_sasa: float

    @property
    def key(self):
        return self.chain, self.number, self.icode

    def as_pepglad_entry(self):
        return [self.chain, [self.number, self.icode], self.score]


@dataclass(frozen=True)
class GuidanceContext:
    seed: ScoredResidue
    residues: tuple[ScoredResidue, ...]

    def as_pepglad_pocket(self):
        return [residue.as_pepglad_entry() for residue in self.residues]


@dataclass(frozen=True)
class GuidanceSelection:
    max_rank: int | str = 3
    surface_rsa_threshold: float = 0.2
    seed_separation: float = 15.0
    context_radius: float = 10.0
    max_context_size: int = 48

    def __post_init__(self):
        if self.max_rank != 'all' and (type(self.max_rank) is not int or self.max_rank < 1):
            raise ValueError('max_rank must be a positive integer or "all"')
        if (not math.isfinite(self.surface_rsa_threshold)
                or not 0 <= self.surface_rsa_threshold <= 1):
            raise ValueError('surface_rsa_threshold must be between 0 and 1')
        if (not math.isfinite(self.seed_separation) or self.seed_separation <= 0
                or not math.isfinite(self.context_radius) or self.context_radius <= 0
                or type(self.max_context_size) is not int or self.max_context_size < 3):
            raise ValueError('Invalid spatial context parameters')


def _relative_sasa(residue, areas):
    chain = residue.get_parent().id
    _, number, icode = residue.id
    try:
        value = areas[chain][f'{number}{icode.strip()}'].relativeTotal
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
        _relative_sasa(residue, areas),
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


def select_contexts(residues, selection):
    """Choose scored local contexts without reading files or running GAPS."""
    candidates = sorted(
        (r for r in residues if r.relative_sasa >= selection.surface_rsa_threshold),
        key=lambda r: (-r.score, r.key),
    )
    contexts = []
    for seed in candidates:
        if any(np.linalg.norm(np.subtract(seed.ca, ctx.seed.ca)) < selection.seed_separation
               for ctx in contexts):
            continue
        nearby = sorted(
            (r for r in residues
             if np.linalg.norm(np.subtract(seed.ca, r.ca)) <= selection.context_radius),
            key=lambda r: (np.linalg.norm(np.subtract(seed.ca, r.ca)),
                           -r.score, r.key),
        )[:selection.max_context_size]
        if len(nearby) < 3 or not any(r.score > 0 for r in nearby):
            continue
        contexts.append(GuidanceContext(seed, tuple(sorted(nearby, key=lambda r: r.key))))
        if selection.max_rank != 'all' and len(contexts) >= selection.max_rank:
            break
    if not contexts:
        raise ValueError('No usable GAPS-guided receptor contexts found')
    return contexts


def select_guidance_contexts(
    input_pdb, scores_json, max_rank=3, surface_rsa_threshold=0.2,
    seed_separation=15.0, context_radius=10.0, max_context_size=48,
):
    """Load GAPS scores, then choose separated surface contexts."""
    selection = GuidanceSelection(
        max_rank, surface_rsa_threshold, seed_separation,
        context_radius, max_context_size)
    return select_contexts(load_scored_residues(input_pdb, scores_json), selection)


def run_gaps(pdb_filepath, out_dir='.', pytorch_device='cuda'):
    """Run GAPS prediction in its dedicated Conda environment."""
    subprocess.run([
        'conda', 'run', '--no-capture-output', '-n', 'gaps', 'python',
        str(Path(__file__).with_name('gaps_predict.py')),
        '--pdb', str(Path(pdb_filepath).resolve()),
        '--out-dir', str(Path(out_dir).resolve()),
        '--device', pytorch_device,
    ], check=True)
