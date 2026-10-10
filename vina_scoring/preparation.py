"""PDB/mmCIF selection, hydrogen handling, and rigid PDBQT preparation.
"""

from __future__ import annotations

import hashlib
import importlib.metadata
import math
import shutil
import subprocess
import sys
from collections import Counter
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

MMCIF_SUFFIXES = {".cif", ".mmcif"}
WATER_RESIDUES = {"DOD", "H2O", "HOH", "WAT"}


@dataclass(frozen=True)
class PreparedComponent:
    """A temporary rigid PDBQT and the provenance used to create it."""

    pdbqt: Path
    metadata: dict[str, Any]


def sha256_file(filename: str | Path) -> str:
    """Return a streaming SHA-256 digest for an input structure."""
    digest = hashlib.sha256()
    with open(filename, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def find_meeko_receptor_preparer() -> Path:
    """Locate Meeko's receptor preparer, including beside virtualenv
    Python."""
    executable = shutil.which("mk_prepare_receptor.py")
    if executable:
        return Path(executable)
    virtualenv_executable = Path(sys.executable).with_name(
        "mk_prepare_receptor.py")
    if virtualenv_executable.is_file():
        return virtualenv_executable
    raise RuntimeError(
        "Meeko executable 'mk_prepare_receptor.py' was not found")


def meeko_version() -> str | None:
    """Return the installed Meeko distribution version."""
    try:
        return importlib.metadata.version("meeko")
    except importlib.metadata.PackageNotFoundError:
        return None


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


class PdbPreparer:
    """Select structure atoms and convert receptor or peptide to rigid
    PDBQT."""

    def __init__(self,
        *, model_index: int = 0, keep_hetero_residues: Sequence[str] = (),
        include_all_hetero: bool = False,
        ph: float = 7.4, hydrogen_policy: str = "keep",
        meeko_templates: Path | None = None,) -> None:
        if model_index < 0:
            raise ValueError("Model index cannot be negative")
        if not math.isfinite(ph):
            raise ValueError("pH must be finite")
        if hydrogen_policy not in {"keep", "rebuild"}:
            raise ValueError("Hydrogen policy must be 'keep' or 'rebuild'")
        self.model_index = model_index
        self.keep_hetero_residues = {
            name.strip().upper()
            for name in keep_hetero_residues if name.strip()}
        self.include_all_hetero = include_all_hetero
        self.ph = float(ph)
        self.hydrogen_policy = hydrogen_policy
        self.meeko_templates = meeko_templates
        try:
            from meeko import ResidueChemTemplates
        except ImportError as exc:
            raise RuntimeError(
                "Meeko is required for PDBQT preparation") from exc
        templates = ResidueChemTemplates.create_from_defaults()
        if self.meeko_templates is not None:
            if not self.meeko_templates.is_file():
                raise ValueError("Meeko template file does not exist: "
                    f"{self.meeko_templates}")
            templates.add_json_file(str(self.meeko_templates))
        self._supported_residue_names = set(templates.residue_templates) | set(
            templates.ambiguous)

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

    def _selected_atom_metrics(self, source: Path) -> tuple[
        int, list[tuple[float, float, float]], Counter[tuple[str, int]],]:
        """Return hydrogen, coordinate, and residue metrics for a
        component."""
        # A selected temporary structure contains exactly one model,
        # even when a later model was selected from the original input.
        _, model = self._load_structure(source, model_index=0)
        hydrogen_count = 0
        heavy_coordinates: list[tuple[float, float, float]] = []
        residue_signatures: Counter[tuple[str, int]] = Counter()
        unsupported: list[str] = []
        for chain in model:
            for residue in chain:
                residue_name = residue.get_resname().strip().upper()
                residue_id = residue.get_id()
                if residue_name not in self._supported_residue_names:
                    insertion = str(residue_id[2]).strip()
                    unsupported.append(f"{residue_name} {chain.get_id()}:"
                        f"{residue_id[1]}{insertion}")
                residue_heavy_atoms = 0
                for atom in residue.get_atoms():
                    if _is_hydrogen(atom):
                        hydrogen_count += 1
                    else:
                        residue_heavy_atoms += 1
                        heavy_coordinates.append(
                            tuple(float(value) for value in atom.coord))
                residue_signatures[(residue_name, residue_heavy_atoms)] += 1
        if unsupported:
            raise ValueError("Unsupported Meeko residue template(s): "
                + ", ".join(unsupported))
        return hydrogen_count, heavy_coordinates, residue_signatures

    @staticmethod
    def _validate_heavy_atoms(source_coordinates: Sequence[Sequence[float]],
        source_residues: Counter[tuple[str, int]], destination: Path,) -> None:
        """Require Meeko conversion to preserve heavy atoms and
        residues."""
        destination_coordinates: list[tuple[float, float, float]] = []
        residue_heavy_counts: dict[tuple[str, str, str], int] = {}
        with destination.open(encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, start=1):
                if not line.startswith(("ATOM  ", "HETATM")):
                    continue
                fields = line.split()
                atom_type = fields[-1].upper() if fields else ""
                if atom_type in {"H", "HD", "HS"}:
                    continue
                residue_key = (
                    line[17:20].strip().upper(),
                    line[21:22].strip(), line[22:27].strip(),)
                residue_heavy_counts[residue_key] = (
                    residue_heavy_counts.get(residue_key, 0) + 1)
                try:
                    coordinate = (
                        float(line[30:38]),
                        float(line[38:46]), float(line[46:54]),)
                except ValueError as exc:
                    raise RuntimeError(
                        f"Invalid PDBQT coordinates in {destination} at line "
                        f"{line_number}") from exc
                destination_coordinates.append(coordinate)

        if len(source_coordinates) != len(destination_coordinates):
            raise RuntimeError(
                "Meeko changed the heavy-atom count during PDBQT conversion "
                f"({len(source_coordinates)} input, "
                f"{len(destination_coordinates)} output): {destination}")
        destination_residues = Counter(
            (residue_name, heavy_count) for (residue_name,
                _chain, _number,), heavy_count in residue_heavy_counts.items())
        if source_residues != destination_residues:
            raise RuntimeError(
                "Meeko changed the residue composition during PDBQT "
                "conversion: " f"{destination}")

        # Meeko writes PDBQT coordinates to three decimals. Match each
        # output atom to an unused source position rather than relying
        # on atom order or names, which the converter may legitimately
        # change.
        unmatched = [
            tuple(map(float, coordinate)) for coordinate in source_coordinates]
        tolerance = 0.002
        for coordinate in destination_coordinates:
            nearest_index = min(range(len(unmatched)),
                key=lambda index: math.dist(coordinate, unmatched[index]),)
            displacement = math.dist(coordinate, unmatched[nearest_index])
            if displacement > tolerance:
                raise RuntimeError("Meeko changed a heavy-atom coordinate by "
                    f"{displacement:.4f} A during PDBQT conversion: "
                    f"{destination}")
            unmatched.pop(nearest_index)

    @staticmethod
    def _format_rigid_ligand(destination: Path) -> None:
        """Wrap rigid PDBQT atoms in the ligand records required by
        Vina."""
        lines = destination.read_text(encoding="utf-8").splitlines()
        remarks = [line for line in lines if line.startswith("REMARK")]
        atoms = [
            line for line in lines if line.startswith(("ATOM  ", "HETATM"))]
        if not atoms:
            raise RuntimeError(f"Meeko wrote no atoms to {destination}")

        # Meeko's receptor-style rigid output can contain repeated atom
        # serials for disconnected fragments. A rigid Vina ligand is one
        # ROOT, so give every atom a unique serial while preserving its
        # atom typing.
        renumbered = [
            f"{line[:6]}{serial:5d}{line[11:]}"
            for serial, line in enumerate(atoms, start=1)]
        rigid_ligand = [
            *remarks[:3], "ROOT", *renumbered, "ENDROOT", "TORSDOF 0",]
        destination.write_text(
            "\n".join(rigid_ligand) + "\n", encoding="utf-8")

    def _convert_to_rigid_pdbqt(
        self, source: Path, destination: Path, *, as_ligand: bool) -> str:
        command = [
            str(find_meeko_receptor_preparer()),
            "--read_pdb", str(source), "--write_pdbqt", str(destination),]
        if self.meeko_templates is not None:
            command.extend(["--add_templates", str(self.meeko_templates)])
        try:
            result = subprocess.run(
                command, check=True, capture_output=True, text=True,)
        except subprocess.CalledProcessError as exc:
            detail = (
                exc.stderr or exc.stdout or "no diagnostic output").strip()
            raise RuntimeError(f"Meeko preparation failed for {source}: "
                f"{' '.join(detail.split())}") from exc
        diagnostic = "\n".join(
            text for text in (result.stdout, result.stderr) if text).strip()
        if "warning" in diagnostic.lower() or "error" in diagnostic.lower():
            raise RuntimeError(
                f"Meeko reported unsafe chemistry for {source}: "
                f"{' '.join(diagnostic.split())}")
        if not destination.is_file() or destination.stat().st_size == 0:
            raise RuntimeError(f"Meeko did not create {destination}")
        if as_ligand:
            self._format_rigid_ligand(destination)
        return diagnostic

    def prepare_component(self, source: Path,
        output_prefix: Path, chain_ids: Sequence[str] | None = None,
        *, as_ligand: bool = False,) -> PreparedComponent:
        """Select one PDB/mmCIF component and create a temporary rigid
        PDBQT."""
        # Meeko's native PDB reader is more reliable than its optional
        # ProDy mmCIF path. mmCIF remains the preferred input because
        # Gemmi converts its struct_conn records into LINK records in
        # this selected PDB.
        selected_structure = output_prefix.with_suffix(".pdb")
        pdbqt = output_prefix.with_suffix(".pdbqt")
        selected_chains, hetero_present, excluded_hetero = (
            self._write_selected_structure(
                source, selected_structure, chain_ids, as_ligand=as_ligand,))
        (
            input_hydrogens, input_heavy_coordinates,
            input_residues,) = self._selected_atom_metrics(selected_structure)
        diagnostic = self._convert_to_rigid_pdbqt(
            selected_structure, pdbqt, as_ligand=as_ligand)
        self._validate_heavy_atoms(
            input_heavy_coordinates, input_residues, pdbqt,)
        diagnostic = " ".join(
            diagnostic.replace(str(selected_structure), str(source)).split())

        warnings: list[str] = []
        if excluded_hetero:
            warnings.append(
                "Excluded HETATM residue names: " + ", ".join(excluded_hetero))
        return PreparedComponent(
            pdbqt=pdbqt, metadata={"source_pdb": str(source),
                "source_sha256": sha256_file(source),
                "model_index": self.model_index,
                "chains": list(selected_chains),
                "hetero_residues_present": hetero_present,
                "included_hetero_residues": (
                    "all_non_water" if as_ligand else (
                        "all" if self.include_all_hetero
                        else sorted(self.keep_hetero_residues)
                    )), "input_hydrogen_count": input_hydrogens,
                "hydrogen_policy": self.hydrogen_policy, "hydrogen_action": (
                    "input_removed_then_meeko_templates"
                    if self.hydrogen_policy == "rebuild"
                    else "input_retained_then_meeko_templates"
                ), "ph": self.ph, "pdbqt_role": "rigid_ligand"
                if as_ligand else "rigid_receptor",
                "preparation_tool": "meeko", "meeko_templates": (
                    str(self.meeko_templates) if self.meeko_templates
                    else "built_in"), "warnings": warnings,},)
