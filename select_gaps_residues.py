"""Select spatial patches with high GAPS scores from a PDB file.

GAPS writes its prediction score to the B-factor column.  This module first
filters residues by that score and then builds connected components using the
minimum heavy-atom distance between residues.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple, Union

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

        return [[chain, [number, icode]] for chain, number, icode in self.residues]


def _mean_bfactor(residue: Residue) -> float:
    b_factors = [atom.get_bfactor() for atom in residue.get_atoms()]
    if not b_factors:
        raise ValueError(f"Residue {residue.get_full_id()} contains no atoms")
    return sum(b_factors) / len(b_factors)


def _residue_id(residue: Residue) -> ResidueId:
    chain = residue.get_parent()
    _, number, insertion_code = residue.id
    return chain.id, number, insertion_code


def find_high_bfactor_spatial_patches(
    pdb_file: Union[str, Path],
    bfactor_threshold: float = 0.5,
    distance_cutoff: float = 5.0,
    min_patch_size: int = 3,
) -> List[SpatialPatch]:
    """Find connected patches of residues with high B-factor values.

    A residue is retained when its mean atomic B-factor is greater than or
    equal to ``bfactor_threshold``. Two retained residues are connected when
    any pair of their non-hydrogen atoms is within ``distance_cutoff`` Angstrom.
    Connections may cross chain boundaries. Only the first model in the PDB is
    considered.

    Patches are ranked by mean B-factor, then by patch size. Hetero residues
    and water are ignored.

    Args:
        pdb_file: A ``*_GAPS_output.pdb`` file.
        bfactor_threshold: Minimum mean residue B-factor (GAPS probability).
        distance_cutoff: Maximum heavy-atom distance defining an edge.
        min_patch_size: Minimum number of residues in a returned patch.

    Returns:
        Spatial patches ordered from highest to lowest mean B-factor.
    """

    if not 0.0 <= bfactor_threshold <= 1.0:
        raise ValueError("bfactor_threshold must be between 0 and 1")
    if distance_cutoff <= 0.0:
        raise ValueError("distance_cutoff must be greater than 0")
    if min_patch_size < 1:
        raise ValueError("min_patch_size must be at least 1")

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
    selected = {
        residue for residue, score in residue_scores.items()
        if score >= bfactor_threshold
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
) -> List[list]:
    """Return the best high-B-factor spatial patch in PepGLAD format.

    An empty list is returned when no patch meets the threshold and minimum
    size. Use :func:`find_high_bfactor_spatial_patches` when all patches or
    their scores are needed.
    """

    patches = find_high_bfactor_spatial_patches(
        pdb_file=pdb_file,
        bfactor_threshold=bfactor_threshold,
        distance_cutoff=distance_cutoff,
        min_patch_size=min_patch_size,
    )
    return patches[0].as_pepglad_pocket() if patches else []
