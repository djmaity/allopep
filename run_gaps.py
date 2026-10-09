""" Run GAPS to predict peptide binding sites """

import math
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple, Union

import freesasa
from Bio.PDB import NeighborSearch, PDBParser
from Bio.PDB.Residue import Residue


ResidueId = Tuple[str, int, str]


@dataclass(frozen=True)
class SpatialPatch:
    """A connected group of high-scoring residues."""

    residues: Tuple[ResidueId, ...]
    mean_bfactor: float
    max_bfactor: float

    def as_pepglad_pocket(self) -> List[list]:
        """Return residues in PepGLAD's ``[chain, [number, icode]]`` format."""

        return [[chain, [number, icode]]
                for chain, number, icode in self.residues]


def _mean_bfactor(residue: Residue) -> float:
    b_factors = [atom.get_bfactor() for atom in residue.get_atoms()]
    if not b_factors:
        raise ValueError(f"Residue {residue.get_full_id()} contains no atoms")
    return sum(b_factors) / len(b_factors)


def _residue_id(residue: Residue) -> ResidueId:
    chain = residue.get_parent()
    _, number, insertion_code = residue.id
    return chain.id, number, insertion_code


def _relative_sasa(residue: Residue, residue_areas: dict) -> float:
    """Return the FreeSASA relative accessibility for a residue."""

    chain, number, insertion_code = _residue_id(residue)
    residue_number = f"{number}{insertion_code.strip()}"
    try:
        relative_sasa = residue_areas[chain][residue_number].relativeTotal
    except KeyError:
        return 0.0
    return relative_sasa if math.isfinite(relative_sasa) else 0.0


def find_high_bfactor_spatial_patches(
    pdb_file: Union[str, Path],
    bfactor_threshold: float = 0.5,
    distance_cutoff: float = 5.0,
    min_patch_size: int = 3,
    surface_rsa_threshold: float = 0.2,
) -> List[SpatialPatch]:
    """Find connected surface patches of residues with high GAPS scores.

    A residue is retained when its mean atomic B-factor is greater than or
    equal to ``bfactor_threshold`` and its relative solvent-accessible surface
    area (RSA) is at least ``surface_rsa_threshold``. RSA is calculated over
    the complete first model so that atoms in neighboring chains can occlude
    one another. A threshold of 0.2 therefore requires at least 20 percent of
    the residue's reference surface area to be exposed. Two retained surface
    residues are connected when any pair of their non-hydrogen atoms is within
    ``distance_cutoff`` Angstrom. Connections may cross chain boundaries.

    Patches are ranked by mean B-factor, then by patch size. Hetero residues
    and water are ignored.
    """

    if not 0.0 <= bfactor_threshold <= 1.0:
        raise ValueError("bfactor_threshold must be between 0 and 1")
    if distance_cutoff <= 0.0:
        raise ValueError("distance_cutoff must be greater than 0")
    if min_patch_size < 1:
        raise ValueError("min_patch_size must be at least 1")
    if not 0.0 <= surface_rsa_threshold <= 1.0:
        raise ValueError("surface_rsa_threshold must be between 0 and 1")

    pdb_path = Path(pdb_file)
    if not pdb_path.is_file():
        raise FileNotFoundError(pdb_path)

    structure = PDBParser(QUIET=True).get_structure(pdb_path.stem, str(pdb_path))
    try:
        model = next(structure.get_models())
    except StopIteration as exc:
        raise ValueError(f"No model found in {pdb_path}") from exc

    residue_scores = {
        residue: _mean_bfactor(residue)
        for residue in model.get_residues()
        if residue.id[0] == " "
    }
    score_candidates = {
        residue for residue, score in residue_scores.items()
        if score >= bfactor_threshold
    }
    if not score_candidates:
        return []

    sasa_result = freesasa.calcBioPDB(structure)[0]
    residue_areas = sasa_result.residueAreas()
    selected = {
        residue for residue in score_candidates
        if _relative_sasa(residue, residue_areas) >= surface_rsa_threshold
    }
    if not selected:
        return []

    adjacency = {residue: set() for residue in selected}
    heavy_atoms = [
        atom
        for residue in selected
        for atom in residue.get_atoms()
        if atom.element != "H"
    ]
    for left, right in NeighborSearch(heavy_atoms).search_all(
            distance_cutoff, level="R"):
        if left is right:
            continue
        adjacency[left].add(right)
        adjacency[right].add(left)

    patches = []
    unseen = set(selected)
    while unseen:
        seed = unseen.pop()
        component = {seed}
        stack = [seed]
        while stack:
            neighbors = adjacency[stack.pop()] & unseen
            unseen.difference_update(neighbors)
            component.update(neighbors)
            stack.extend(neighbors)

        if len(component) < min_patch_size:
            continue

        ordered_residues = tuple(sorted(_residue_id(residue)
                                        for residue in component))
        scores = [residue_scores[residue] for residue in component]
        patches.append(SpatialPatch(
            residues=ordered_residues,
            mean_bfactor=sum(scores) / len(scores),
            max_bfactor=max(scores),
        ))

    patches.sort(
        key=lambda patch: (patch.mean_bfactor, len(patch.residues)),
        reverse=True,
    )
    return patches


def select_contiguous_high_bfactor_residues(
    pdb_file: Union[str, Path],
    bfactor_threshold: float = 0.5,
    distance_cutoff: float = 5.0,
    min_patch_size: int = 3,
    surface_rsa_threshold: float = 0.2,
) -> List[list]:
    """Return the best high-GAPS-score surface patch in PepGLAD format."""

    patches = find_high_bfactor_spatial_patches(
        pdb_file=pdb_file,
        bfactor_threshold=bfactor_threshold,
        distance_cutoff=distance_cutoff,
        min_patch_size=min_patch_size,
        surface_rsa_threshold=surface_rsa_threshold,
    )
    return patches[0].as_pepglad_pocket() if patches else []


def run_gaps(pdb_filepath, out_dir='.', pytorch_device='cuda'):
    """Run GAPS prediction in its dedicated Conda environment."""

    subprocess.run([
        'conda', 'run', '--no-capture-output', '-n', 'gaps', 'python',
        str(Path(__file__).with_name('gaps_predict.py')),
        '--pdb', str(Path(pdb_filepath).resolve()),
        '--out-dir', str(Path(out_dir).resolve()),
        '--device', pytorch_device,
    ], check=True)
