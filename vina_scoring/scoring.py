"""Rigid-pose AutoDock Vina scoring and score diagnostics."""

from __future__ import annotations

import gc
import math
import os
from collections import OrderedDict
from collections.abc import Iterable, Sequence
from pathlib import Path
from typing import Any

import numpy as np

from models import BoxSpec
from preparation import sha256_file

SCORE_TERMS = (
    "total", "ligand_inter",
    "flex_receptor_inter", "other_inter", "flex_receptor_intra",
    "ligand_intra", "torsional", "ligand_intra_best_pose",)


def named_energies(values: Iterable[float]) -> dict[str, float]:
    """Attach names to components returned by the supported Vina API."""
    energies = [float(value) for value in values]
    if len(energies) != len(SCORE_TERMS):
        raise ValueError(f"Vina score() returned {len(energies)} components; "
            f"expected {len(SCORE_TERMS)}")
    return dict(zip(SCORE_TERMS, energies))


def size_normalized_diagnostics(
    energies: dict[str, float], heavy_atoms: int, residues: int
) -> dict[str, float | None]:
    """Calculate the requested size diagnostics without changing the raw
    score."""
    interaction = energies.get("ligand_inter")
    total = energies.get("total")

    def divide(value: float | None, denominator: int) -> float | None:
        return (
            value / denominator
            if value is not None and denominator > 0 else None)

    return {
        "total_per_heavy_atom": divide(total, heavy_atoms),
        "interaction_per_heavy_atom": divide(interaction, heavy_atoms),
        "total_per_residue": divide(total, residues),
        "interaction_per_residue": divide(interaction, residues),}


def parse_pdbqt(filename: str | os.PathLike[str]) -> dict[str, Any]:
    """Extract coordinates and peptide size metrics from a PDBQT file.
    """
    coordinates: list[tuple[float, float, float]] = []
    residues: set[tuple[str, str, str, str]] = set()
    heavy_atoms = 0
    branch_count = 0
    torsdof: int | None = None

    with open(filename, encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            record = line[:6].strip()
            if record in {"ATOM", "HETATM"}:
                try:
                    coordinate = (
                        float(line[30:38]),
                        float(line[38:46]), float(line[46:54]),)
                except ValueError as exc:
                    raise ValueError(
                        f"Invalid coordinates in {filename} at line "
                        f"{line_number}") from exc
                coordinates.append(coordinate)
                residues.add((line[21:22].strip(), line[22:26].strip(),
                        line[26:27].strip(), line[17:20].strip(),))
                fields = line.split()
                atom_type = fields[-1].upper() if fields else ""
                if atom_type not in {"H", "HD", "HS"}:
                    heavy_atoms += 1
            elif record == "BRANCH":
                branch_count += 1
            elif record == "TORSDO":  # first six columns of TORSDOF
                fields = line.split()
                if len(fields) == 2:
                    try:
                        torsdof = int(fields[1])
                    except ValueError:
                        pass

    if not coordinates:
        raise ValueError(f"No ATOM/HETATM coordinates found in {filename}")
    return {
        "coordinates": coordinates, "atom_count": len(coordinates),
        "heavy_atom_count": heavy_atoms, "residue_count": len(residues),
        "active_torsions": torsdof if torsdof is not None else branch_count,}


def rigid_residue_ligands(
    filename: str | os.PathLike[str],
) -> list[str]:
    """Split a rigid peptide PDBQT into residue-sized ligands.

    Vina 1.2.7 allocates a quadratic per-atom lookup table even when a
    ligand has no active torsions or intramolecular score terms. Residue
    fragments preserve the local atom typing needed by protein residues
    while keeping that redundant allocation small.
    """
    residues: OrderedDict[tuple[str, str, str, str], list[str]] = (
        OrderedDict())
    with open(filename, encoding="utf-8") as handle:
        for line in handle:
            if line[:6].strip() not in {"ATOM", "HETATM"}:
                continue
            key = (
                line[21:22], line[22:26], line[26:27], line[17:20],)
            residues.setdefault(key, []).append(line.rstrip("\n"))
    if not residues:
        raise ValueError(f"No ATOM/HETATM records found in {filename}")
    return [
        "ROOT\n" + "\n".join(atom_lines)
        + "\nENDROOT\nTORSDOF 0\n"
        for atom_lines in residues.values()]


def box_from_coordinates(
    coordinate_sets: Iterable[Iterable[Sequence[float]]], padding: float
) -> tuple[list[float], list[float]]:
    """Calculate a box around the union of peptide coordinate sets."""
    if not math.isfinite(padding) or padding < 0:
        raise ValueError("Box padding must be a finite, non-negative number")
    coordinates = [
        tuple(float(value) for value in coordinate)
        for coordinate_set in coordinate_sets for coordinate in coordinate_set]
    if not coordinates:
        raise ValueError("Cannot calculate a box without peptide coordinates")
    if any(len(coordinate) != 3 for coordinate in coordinates):
        raise ValueError("Every coordinate must contain exactly three values")
    if not all(math.isfinite(value)
        for coordinate in coordinates for value in coordinate):
        raise ValueError("Coordinates must be finite")

    minima = [
        min(coordinate[axis] for coordinate in coordinates)
        for axis in range(3)]
    maxima = [
        max(coordinate[axis] for coordinate in coordinates)
        for axis in range(3)]
    center = [(low + high) / 2.0 for low, high in zip(minima, maxima)]
    size = [(high - low) + 2.0 * padding for low, high in zip(minima, maxima)]
    if any(dimension <= 0 for dimension in size):
        raise ValueError(
            "Box dimensions must be positive; increase box padding")
    return center, size


def resolve_box(metrics: Sequence[dict[str, Any]], box: BoxSpec, padding: float
) -> tuple[list[float], list[float], str]:
    """Resolve and validate an explicit box or the union of peptide
    bounds."""
    coordinate_sets = [item["coordinates"] for item in metrics]
    if box.is_explicit:
        assert box.center is not None and box.size is not None
        center = [float(value) for value in box.center]
        size = [float(value) for value in box.size]
        source = "explicit"
    else:
        center, size = box_from_coordinates(coordinate_sets, padding)
        source = "union_of_peptide_bounds"
    validate_box(center, size, coordinate_sets)
    return center, size, source


def validate_box(center: Sequence[float], size: Sequence[float],
    coordinate_sets: Iterable[Iterable[Sequence[float]]],
    tolerance: float = 1e-6,) -> None:
    """Reject malformed boxes and peptide atoms outside their bounds."""
    if len(center) != 3 or len(size) != 3:
        raise ValueError("Box center and size must each contain three values")
    if not all(math.isfinite(float(value)) for value in (*center, *size)):
        raise ValueError("Box center and size must be finite")
    if any(float(dimension) <= 0 for dimension in size):
        raise ValueError("Box dimensions must be positive")

    lower = [
        float(c) - float(s) / 2.0 - tolerance for c, s in zip(center, size)]
    upper = [
        float(c) + float(s) / 2.0 + tolerance for c, s in zip(center, size)]
    for coordinate_set in coordinate_sets:
        for coordinate in coordinate_set:
            if any(float(coordinate[axis]) < lower[axis]
                or float(coordinate[axis]) > upper[axis] for axis in range(3)):
                raise ValueError(
                    "At least one peptide atom is outside the map box")


class VinaPoseScorer:
    """Score one rigid peptide pose as additive residue interactions."""

    def __init__(self, receptor_pdbqt: str | os.PathLike[str],
        center: Sequence[float], box_size: Sequence[float],
        *, spacing: float = 0.375, cpu: int = 1,) -> None:
        try:
            from vina import Vina
        except ImportError as exc:
            raise RuntimeError(
                "The Python package 'vina' is required") from exc
        if not math.isfinite(spacing) or spacing <= 0:
            raise ValueError("Grid spacing must be a finite, positive number")
        if cpu < 1:
            raise ValueError("CPU count must be at least 1")

        self.receptor_pdbqt = Path(receptor_pdbqt)
        if not self.receptor_pdbqt.is_file():
            raise ValueError(
                f"Receptor PDBQT does not exist: {self.receptor_pdbqt}")
        self.center = [float(value) for value in center]
        self.box_size = [float(value) for value in box_size]
        self.spacing = float(spacing)

        self._vina: Any | None = Vina(
            sf_name="vina", cpu=cpu, verbosity=0)
        self._vina.set_receptor(str(self.receptor_pdbqt))
        # Build all Vina atom-type maps once so every residue fragment
        # can reuse them. Fragmenting avoids Vina 1.2.7's quadratic
        # whole-ligand lookup table without approximating the score.
        self._vina.compute_vina_maps(
            center=self.center, box_size=self.box_size,
            spacing=self.spacing, force_even_voxels=True,)
        self._used = False

    def score_pose(
        self, peptide_pdbqt: str | os.PathLike[str], metrics: dict[str, Any]
    ) -> dict[str, Any]:
        """Score one rigid, supplied peptide pose without optimization.
        """
        if self._used:
            raise RuntimeError(
                "A VinaPoseScorer can score only one peptide pose")
        if self._vina is None:
            raise RuntimeError("The Vina scorer has already been closed")
        self._used = True
        peptide_path = Path(peptide_pdbqt)
        if metrics["active_torsions"] != 0:
            raise ValueError("Rigid peptide PDBQT has "
                f"{metrics['active_torsions']} active torsion(s): "
                f"{peptide_path}")
        accumulated = [0.0] * len(SCORE_TERMS)
        for residue_ligand in rigid_residue_ligands(peptide_path):
            self._vina.set_ligand_from_string(residue_ligand)
            # The public Python wrapper rounds each score to three
            # decimals. Sum the native values first so fragment rounding
            # cannot accumulate across a peptide.
            native_vina = getattr(self._vina, "_vina", None)
            if native_vina is None:
                raise RuntimeError(
                    "The installed Vina package does not expose raw scores")
            values = [float(value) for value in native_vina.score()]
            if len(values) != len(SCORE_TERMS):
                raise ValueError(
                    f"Vina score() returned {len(values)} components; "
                    f"expected {len(SCORE_TERMS)}")
            for index, value in enumerate(values):
                accumulated[index] += value
        raw = named_energies(np.around(accumulated, decimals=3))
        return {
            "peptide_pdbqt_sha256": sha256_file(peptide_path),
            "pose_metrics": {
                key: value
                for key, value in metrics.items() if key != "coordinates"
            }, "scores": raw, "size_normalized": size_normalized_diagnostics(
                raw, heavy_atoms=metrics["heavy_atom_count"],
                residues=metrics["residue_count"],),}

    def close(self) -> None:
        """Release the native Vina object and its affinity maps."""
        self._vina = None
        gc.collect()


def estimate_vina_memory_bytes(
    box_size: Sequence[float], spacing: float) -> int:
    """Return a conservative memory reservation for one Vina map worker.
    """
    if spacing <= 0:
        raise ValueError("Grid spacing must be positive")
    voxel_counts: list[int] = []
    for dimension in box_size:
        voxels = math.ceil(float(dimension) / spacing)
        if voxels % 2:
            voxels += 1
        voxel_counts.append(voxels + 1)
    grid_points = math.prod(voxel_counts)
    # Vina creates maps for all supported atom types. Eight-byte values
    # and a 3x overhead factor conservatively cover map copies and
    # native work arrays.
    map_bytes = grid_points * 22 * 8 * 3
    process_baseline = 1024**3
    return process_baseline + map_bytes
