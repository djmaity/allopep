"""Meeko rigid PDBQT preparation and heavy-atom validation."""

from __future__ import annotations

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

from ..fileio import sha256
from .selection import StructureSelector, _is_hydrogen


@dataclass(frozen=True)
class PreparedComponent:
    """A temporary rigid PDBQT and the provenance used to create it."""

    pdbqt: Path
    metadata: dict[str, Any]


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


def installed_version(distribution: str) -> str:
    """Return an installed distribution version for result provenance."""
    try:
        return importlib.metadata.version(distribution)
    except importlib.metadata.PackageNotFoundError:
        return ""


class MeekoPreparer:
    """Convert selected structures to rigid PDBQT and verify the result."""

    def __init__(self,
        *, model_index: int = 0, keep_hetero_residues: Sequence[str] = (),
        include_all_hetero: bool = False,
        ph: float = 7.4, hydrogen_policy: str = "keep",
        meeko_templates: Path | None = None,) -> None:
        self.selector = StructureSelector(
            model_index=model_index,
            keep_hetero_residues=keep_hetero_residues,
            include_all_hetero=include_all_hetero,
            hydrogen_policy=hydrogen_policy)
        self.model_index = model_index
        self.keep_hetero_residues = self.selector.keep_hetero_residues
        self.include_all_hetero = include_all_hetero
        self.ph = float(ph)
        self.hydrogen_policy = hydrogen_policy
        self.meeko_templates = meeko_templates
        try:
            from meeko import ResidueChemTemplates
        except ImportError as exc:
            raise RuntimeError("Meeko is required for PDBQT preparation") from exc
        templates = ResidueChemTemplates.create_from_defaults()
        if meeko_templates is not None:
            templates.add_json_file(str(meeko_templates))
        self._supported_residue_names = (
            set(templates.residue_templates) | set(templates.ambiguous))

    def _selected_atom_metrics(self, source: Path) -> tuple[
        int, list[tuple[float, float, float]], Counter[tuple[str, int]],]:
        """Return hydrogen, coordinate, and residue metrics for a
        component."""
        # A selected temporary structure contains exactly one model,
        # even when a later model was selected from the original input.
        _, model = self.selector._load_structure(source, model_index=0)
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
            self.selector._write_selected_structure(
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
                "source_sha256": sha256(source),
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
