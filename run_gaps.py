""" Run GAPS to predict peptide binding sites """

import os
from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple, Union

import torch as pt
from Bio.PDB import NeighborSearch, PDBParser
from Bio.PDB.Residue import Residue
from tqdm import tqdm

from GAPS.model.model import Model
from GAPS.src.config import config_model
from GAPS.src.data_encoding import en_features, en_structure, ext_topology
from GAPS.src.dataset import StructuresDataset, col_batch
from GAPS.src.structure import (concatenate_chains, encode_bfactor,
                                split_by_chain)
from GAPS.src.structure_io import save_pdb


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


def find_high_bfactor_spatial_patches(
    pdb_file: Union[str, Path],
    bfactor_threshold: float = 0.5,
    distance_cutoff: float = 5.0,
    min_patch_size: int = 3,
) -> List[SpatialPatch]:
    """Find connected patches of residues with high GAPS scores.

    A residue is retained when its mean atomic B-factor is greater than or
    equal to ``bfactor_threshold``. Two retained residues are connected when
    any pair of their non-hydrogen atoms is within ``distance_cutoff`` Angstrom.
    Connections may cross chain boundaries. Only the first model in the PDB is
    considered.

    Patches are ranked by mean B-factor, then by patch size. Hetero residues
    and water are ignored.
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
    """Return the best high-GAPS-score patch in PepGLAD format."""

    patches = find_high_bfactor_spatial_patches(
        pdb_file=pdb_file,
        bfactor_threshold=bfactor_threshold,
        distance_cutoff=distance_cutoff,
        min_patch_size=min_patch_size,
    )
    return patches[0].as_pepglad_pocket() if patches else []


def run_gaps(pdb_filepaths, out_dir='.', pytorch_device='cuda'):
    """ Get peptide binging site prediction using GAPS """

    model_weights = 'GAPS/checkpoints/model.pt'
    device = pt.device(pytorch_device)

    model = Model(config_model)
    model.load_state_dict(pt.load(model_weights,
                                  map_location=pt.device(pytorch_device),
                                  weights_only=True)
    )
    model = model.eval().to(device)

    if not isinstance(pdb_filepaths, list):
        pdb_filepaths = [pdb_filepaths]
    dataset = StructuresDataset(pdb_filepaths, with_preprocessing=True)

    os.makedirs(out_dir, exist_ok=True)
    with pt.no_grad():
        for subunits, filepath in tqdm(dataset):
            structure = concatenate_chains(subunits)
            X, M = en_structure(structure)
            q = en_features(structure)[0]
            ids_topk, _, _, _, _ = ext_topology(X, 64)
            X, ids_topk, q, M = col_batch([[X, ids_topk, q, M]])
            z = model(X.to(device), ids_topk.to(device), q.to(device), M.float().to(device))
            p = pt.sigmoid(z)
            structure = encode_bfactor(structure, p.cpu().numpy())

            filename = os.path.splitext(os.path.basename(filepath))[0]
            out_file = os.path.join(out_dir, filename + '_GAPS_output.pdb')
            save_pdb(split_by_chain(structure), out_file)


if __name__ == '__main__':
    # TODO: add argument parser
    run_gaps(['input/2VH7.pdb', 'input/5lvp_holo.pdb'], out_dir='output')
