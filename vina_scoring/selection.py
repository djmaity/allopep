"""PDB/mmCIF atom selection and connectivity preservation."""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Any

MMCIF_SUFFIXES = {".cif", ".mmcif"}
WATER_RESIDUES = {"DOD", "H2O", "HOH", "WAT"}


def _chosen_altloc_atom_ids(model: Any) -> set[int]:
    """Choose the highest-occupancy conformer for each disordered atom.
    """
    selected: set[int] = set()
    for chain in model:
        for residue in chain:
            by_name: dict[str, list[Any]] = {}
            for atom in residue.get_unpacked_list():
                by_name.setdefault(atom.get_name(), []).append(atom)
            for atoms in by_name.values():
                chosen = max(atoms, key=lambda atom: (atom.get_occupancy()
                        if atom.get_occupancy() is not None else 0.0,
                        atom.get_altloc() == "A", atom.get_altloc(),),)
                selected.add(id(chosen))
    return selected


def _is_hydrogen(atom: Any) -> bool:
    """Return whether a Biopython atom represents hydrogen."""
    element = str(getattr(atom, "element", "") or "").strip().upper()
    return element == "H" or (
        not element and atom.get_name().upper().startswith("H"))


class StructureSelector:
    """Select PDB/mmCIF atoms and retain applicable connectivity."""

    def __init__(self,
        *, model_index: int = 0, keep_hetero_residues: Sequence[str] = (),
        include_all_hetero: bool = False,
        hydrogen_policy: str = "keep",) -> None:
        if model_index < 0:
            raise ValueError("Model index cannot be negative")
        if hydrogen_policy not in {"keep", "rebuild"}:
            raise ValueError("Hydrogen policy must be 'keep' or 'rebuild'")
        self.model_index = model_index
        self.keep_hetero_residues = {
            name.strip().upper()
            for name in keep_hetero_residues if name.strip()}
        self.include_all_hetero = include_all_hetero
        self.hydrogen_policy = hydrogen_policy

    def _load_structure(self, source: Path, *, model_index: int | None = None
    ) -> tuple[Any, Any]:
        try:
            from Bio.PDB import MMCIFParser, PDBParser
        except ImportError as exc:
            raise RuntimeError(
                "Biopython is required for structure preparation") from exc

        if not source.is_file():
            raise ValueError(f"Input structure does not exist: {source}")
        if source.suffix.lower() in MMCIF_SUFFIXES:
            parser = MMCIFParser(QUIET=True)
        else:
            parser = PDBParser(QUIET=True)
        structure = parser.get_structure(source.stem, str(source))
        models = list(structure.get_models())
        selected_model_index = (
            self.model_index if model_index is None else model_index)
        if selected_model_index >= len(models):
            raise ValueError(
                f"Model index {selected_model_index} is invalid for {source}; "
                f"the structure contains {len(models)} model(s)")
        return structure, models[selected_model_index]

    def available_chains(self, source: Path) -> tuple[str, ...]:
        """Return chain IDs in the selected PDB/mmCIF model."""
        _, model = self._load_structure(source)
        return tuple(chain.get_id() for chain in model)

    def _validate_chains(
        self, model: Any, source: Path, chain_ids: Sequence[str] | None
    ) -> tuple[str, ...]:
        available = tuple(chain.get_id() for chain in model)
        selected = tuple(dict.fromkeys(chain_ids or available))
        unknown = sorted(set(selected) - set(available))
        if unknown:
            raise ValueError(
                f"Unknown chain(s) in {source}: {', '.join(unknown)}; "
                f"available chains: {', '.join(available)}")
        if not selected:
            raise ValueError(f"No chains were selected from {source}")
        return selected

    def _accept_residue(self, residue: Any, *, as_ligand: bool) -> bool:
        residue_name = residue.get_resname().strip().upper()
        if as_ligand:
            # The selected peptide chain is itself the scored ligand.
            # Modified residues and caps are therefore part of the pose,
            # not incidental receptor hetero compounds. Waters are never
            # part of that pose.
            return residue_name not in WATER_RESIDUES
        if residue.get_id()[0] == " ":
            return True
        return self.include_all_hetero or (
            residue_name in self.keep_hetero_residues)

    def selected_coordinates(
        self, source: Path, chain_ids: Sequence[str] | None = None,
        *, as_ligand: bool = False,) -> list[tuple[float, float, float]]:
        """Read selected PDB/mmCIF coordinates for box and memory
        estimation."""
        _, model = self._load_structure(source)
        selected_chains = set(self._validate_chains(model, source, chain_ids))
        selected_atoms = _chosen_altloc_atom_ids(model)
        coordinates: list[tuple[float, float, float]] = []
        for chain in model:
            if chain.get_id() not in selected_chains:
                continue
            for residue in chain:
                if not self._accept_residue(residue, as_ligand=as_ligand):
                    continue
                for atom in residue.get_unpacked_list():
                    if id(atom) in selected_atoms:
                        coordinates.append(
                            tuple(float(value) for value in atom.coord))
        if not coordinates:
            raise ValueError(f"Selection from {source} contains no atoms")
        return coordinates

    def _write_selected_structure(
        self, source: Path, destination: Path, chain_ids: Sequence[str] | None,
        *, as_ligand: bool,) -> tuple[tuple[str, ...], list[str], list[str]]:
        try:
            from Bio.PDB import MMCIFIO, PDBIO, Select
        except ImportError as exc:
            raise RuntimeError(
                "Biopython is required for structure preparation") from exc

        structure, model = self._load_structure(source)
        selected_chains = self._validate_chains(model, source, chain_ids)
        selected_chain_set = set(selected_chains)
        selected_atoms = _chosen_altloc_atom_ids(model)
        hetero_present = sorted({residue.get_resname().strip().upper()
                for chain in model if chain.get_id() in selected_chain_set
                for residue in chain if residue.get_id()[0] != " "})
        excluded_hetero = sorted({residue.get_resname().strip().upper()
                for chain in model if chain.get_id() in selected_chain_set
                for residue in chain if residue.get_id()[0] != " "
                and not self._accept_residue(residue, as_ligand=as_ligand)})

        preparer = self

        class ComponentSelect(Select):
            def accept_model(self, candidate: Any) -> bool:
                return candidate is model

            def accept_chain(self, chain: Any) -> bool:
                return chain.get_id() in selected_chain_set

            def accept_residue(self, residue: Any) -> bool:
                return preparer._accept_residue(residue, as_ligand=as_ligand)

            def accept_atom(self, atom: Any) -> bool:
                return id(atom) in selected_atoms and not (
                    preparer.hydrogen_policy == "rebuild"
                    and _is_hydrogen(atom))

        destination.parent.mkdir(parents=True, exist_ok=True)
        writer = (
            MMCIFIO()
            if destination.suffix.lower() in MMCIF_SUFFIXES else PDBIO())
        writer.set_structure(structure)
        writer.save(str(destination), ComponentSelect())
        self._restore_connections(source, destination)
        selected_atom_count = sum(1 for chain in model
            if chain.get_id() in selected_chain_set for residue in chain
            if self._accept_residue(residue, as_ligand=as_ligand)
            for atom in residue.get_unpacked_list()
            if id(atom) in selected_atoms
            and not (self.hydrogen_policy == "rebuild" and _is_hydrogen(atom)))
        if selected_atom_count == 0:
            raise ValueError(f"Selection from {source} contains no atoms")
        return selected_chains, hetero_present, excluded_hetero

    @staticmethod
    def _restore_connections(source: Path, destination: Path) -> None:
        """Transfer source LINK/struct_conn records retained by the
        selection."""
        try:
            import gemmi
        except ImportError as exc:
            raise RuntimeError(
                "Gemmi is required to preserve connectivity") from exc

        source_structure = gemmi.read_structure(str(source))
        selected_structure = gemmi.read_structure(str(destination))
        if not selected_structure:
            raise ValueError(
                f"Selection from {source} contains no coordinate model")
        selected_atoms = {
            (chain.name, residue.name, str(residue.seqid), atom.name)
            for model in selected_structure
            for chain in model for residue in chain for atom in residue}
        selected_structure.connections.clear()
        for connection in source_structure.connections:
            partners = (connection.partner1, connection.partner2)
            if all((partner.chain_name, partner.res_id.name,
                    str(partner.res_id.seqid), partner.atom_name,
                ) in selected_atoms for partner in partners):
                selected_structure.connections.append(connection)

        if destination.suffix.lower() in MMCIF_SUFFIXES:
            selected_structure.make_mmcif_document().write_file(
                str(destination))
        else:
            selected_structure.write_pdb(str(destination))
