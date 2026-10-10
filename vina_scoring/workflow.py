"""High-level preparation and scoring workflow for one complex."""

from __future__ import annotations

import gc
import importlib.metadata
import tempfile
from pathlib import Path
from typing import Any

from models import BoxSpec, ScoreGroup, ScoreJob, WorkflowConfig
from preparation import PdbPreparer, meeko_version
from scoring import (
    SCORE_TERMS,
    VinaPoseScorer,
    box_from_coordinates,
    estimate_vina_memory_bytes,
    parse_pdbqt,
    resolve_box,
)

NORMALIZED_TERMS = (
    "total_per_heavy_atom", "interaction_per_heavy_atom",
    "total_per_residue", "interaction_per_residue",)


def hetero_policy_description(config: WorkflowConfig) -> str:
    """Describe the role-aware hetero-residue policy for output
    provenance."""
    if config.include_all_hetero:
        receptor_policy = "all"
    elif config.keep_hetero_residues:
        receptor_policy = ",".join(sorted(config.keep_hetero_residues))
    else:
        receptor_policy = "none"
    return f"peptide:all_non_water;receptor:{receptor_policy}"


RESULT_FIELDS = (
    "id", "status", "error", "input_mode",
    "complex_pdb", "receptor_pdb", "peptide_pdb", "complex_sha256",
    "receptor_sha256", "peptide_sha256", "model_index", "peptide_chain",
    "receptor_chains", "hydrogen_policy", "ph", "included_hetero_residues",
    "warnings", "box_source", "box_center_x", "box_center_y",
    "box_center_z", "box_size_x", "box_size_y", "box_size_z", "spacing",
    "poses_in_map_group", "peptide_atom_count", "peptide_heavy_atom_count",
    "peptide_residue_count", "peptide_active_torsions",
    "score_units", *(f"score_{term}" for term in SCORE_TERMS),
    *(f"normalized_{term}" for term in NORMALIZED_TERMS),
    "vina_version", "meeko_version",)


def _package_version(distribution: str) -> str | None:
    try:
        return importlib.metadata.version(distribution)
    except importlib.metadata.PackageNotFoundError:
        return None


def _preparer(config: WorkflowConfig) -> PdbPreparer:
    return PdbPreparer(model_index=config.model_index,
        keep_hetero_residues=config.keep_hetero_residues,
        include_all_hetero=config.include_all_hetero,
        ph=config.ph, hydrogen_policy=config.hydrogen_policy,
        meeko_templates=config.meeko_templates,)


def _actual_receptor_chains(
    job: ScoreJob, preparer: PdbPreparer) -> tuple[str, ...]:
    if job.input_mode == "separate":
        assert job.receptor_pdb is not None
        return preparer.available_chains(job.receptor_pdb)
    assert job.complex_pdb is not None and job.peptide_chain is not None
    if job.receptor_chains:
        return job.receptor_chains
    return tuple(chain for chain in preparer.available_chains(job.complex_pdb)
        if chain != job.peptide_chain)


def failure_result(
    job: ScoreJob, error: str, config: WorkflowConfig) -> dict[str, Any]:
    """Create a complete CSV-compatible failure record."""
    record: dict[str, Any] = {field: "" for field in RESULT_FIELDS}
    record.update({"id": job.job_id,
            "status": "error", "error": error, "input_mode": job.input_mode,
            "complex_pdb": str(job.complex_pdb or ""),
            "receptor_pdb": str(job.receptor_pdb or ""),
            "peptide_pdb": str(job.peptide_pdb or ""),
            "model_index": config.model_index,
            "peptide_chain": job.peptide_chain or "",
            "receptor_chains": ";".join(job.receptor_chains),
            "hydrogen_policy": config.hydrogen_policy, "ph": config.ph,
            "included_hetero_residues": hetero_policy_description(config),
            "spacing": config.spacing,
            "vina_version": _package_version("vina") or "",
            "meeko_version": meeko_version() or "",})
    return record


class VinaScoringWorkflow:
    """Prepare and score isolated protein-peptide complexes."""

    def __init__(self, config: WorkflowConfig) -> None:
        self.config = config
        self.preparer = _preparer(config)

    def estimate_group_memory(self, group: ScoreGroup) -> int:
        """Estimate map memory from the group's source peptide
        coordinates."""
        coordinate_sets: list[list[tuple[float, float, float]]] = []
        for job in group.jobs:
            if job.input_mode == "combined":
                assert (
                    job.complex_pdb is not None
                    and job.peptide_chain is not None)
                coordinates = self.preparer.selected_coordinates(
                    job.complex_pdb, (job.peptide_chain,), as_ligand=True,)
            else:
                assert job.peptide_pdb is not None
                coordinates = self.preparer.selected_coordinates(
                    job.peptide_pdb, as_ligand=True)
            coordinate_sets.append(coordinates)

        box = group.jobs[0].box
        if box.is_explicit:
            assert box.size is not None
            size = list(box.size)
        else:
            _, size = box_from_coordinates(
                coordinate_sets, self.config.box_padding)
            # Rebuilt hydrogen coordinates can extend modestly beyond
            # the input heavy-atom bounds. Reserve for a box two
            # Angstrom larger per axis.
            size = [dimension + 2.0 for dimension in size]
        return estimate_vina_memory_bytes(size, self.config.spacing)

    def _prepare_receptor(self, job: ScoreJob, directory: Path) -> Any:
        receptor_chains = _actual_receptor_chains(job, self.preparer)
        if not receptor_chains:
            raise ValueError(f"Job {job.job_id!r} selected no receptor chains")
        source = job.complex_pdb or job.receptor_pdb
        assert source is not None
        return self.preparer.prepare_component(
            source, directory / "receptor", receptor_chains)

    def _prepare_peptide(self, job: ScoreJob, output_prefix: Path) -> Any:
        if job.input_mode == "combined":
            assert (
                job.complex_pdb is not None and job.peptide_chain is not None)
            return self.preparer.prepare_component(job.complex_pdb,
                output_prefix, (job.peptide_chain,), as_ligand=True,)
        assert job.peptide_pdb is not None
        return self.preparer.prepare_component(
            job.peptide_pdb, output_prefix, as_ligand=True)

    def _success_result(self, job: ScoreJob, receptor: Any, peptide: Any,
        scored: dict[str, Any], center: list[float], size: list[float],
        box_source: str, group_size: int, vina_version: str | None,
        installed_meeko_version: str | None,) -> dict[str, Any]:
        receptor_metadata = receptor.metadata
        peptide_metadata = peptide.metadata
        warnings = list(dict.fromkeys(
                receptor_metadata["warnings"] + peptide_metadata["warnings"]))
        if job.input_mode == "combined" and not job.receptor_chains:
            warnings.insert(0,
                "Receptor chains were inferred as every non-peptide chain.",)
        record: dict[str, Any] = {field: "" for field in RESULT_FIELDS}
        record.update(
            {"id": job.job_id, "status": "ok", "input_mode": job.input_mode,
                "complex_pdb": str(job.complex_pdb or ""),
                "receptor_pdb": str(job.receptor_pdb or ""),
                "peptide_pdb": str(job.peptide_pdb or ""),
                "model_index": self.config.model_index,
                "peptide_chain": job.peptide_chain or "",
                "receptor_chains": ";".join(receptor_metadata["chains"]),
                "hydrogen_policy": self.config.hydrogen_policy,
                "ph": self.config.ph,
                "included_hetero_residues": hetero_policy_description(
                    self.config), "warnings": " | ".join(warnings),
                "box_source": box_source, "box_center_x": center[0],
                "box_center_y": center[1], "box_center_z": center[2],
                "box_size_x": size[0], "box_size_y": size[1],
                "box_size_z": size[2], "spacing": self.config.spacing,
                "poses_in_map_group": group_size,
                "score_units": "nominal_kcal_per_mol",
                "vina_version": vina_version or "",
                "meeko_version": installed_meeko_version or "",})
        if job.input_mode == "combined":
            record["complex_sha256"] = receptor_metadata["source_sha256"]
        else:
            record["receptor_sha256"] = receptor_metadata["source_sha256"]
            record["peptide_sha256"] = peptide_metadata["source_sha256"]

        metrics = scored["pose_metrics"]
        record.update({"peptide_atom_count": metrics["atom_count"],
                "peptide_heavy_atom_count": metrics["heavy_atom_count"],
                "peptide_residue_count": metrics["residue_count"],
                "peptide_active_torsions": metrics["active_torsions"],})
        record.update({f"score_{name}": value
                for name, value in scored["scores"].items()})
        record.update({f"normalized_{name}": value
                for name, value in scored["size_normalized"].items()})
        return record

    def score_group(self, group: ScoreGroup) -> list[dict[str, Any]]:
        """Prepare and score one pose in an isolated worker group."""
        if not group.jobs:
            return []
        if len(group.jobs) != 1:
            raise ValueError(
                "Each Vina worker group must contain exactly one job")
        job = group.jobs[0]
        with tempfile.TemporaryDirectory(
            prefix="vina-rigid-score-") as temporary:
            work_directory = Path(temporary)
            receptor = self._prepare_receptor(job, work_directory)
            peptide = self._prepare_peptide(
                job, work_directory / "peptide")
            metrics = parse_pdbqt(peptide.pdbqt)
            center, size, box_source = resolve_box(
                [metrics], job.box, self.config.box_padding,)
            scorer: VinaPoseScorer | None = None
            try:
                scorer = VinaPoseScorer(receptor.pdbqt,
                    center, size, spacing=self.config.spacing, cpu=1,)
                scored = scorer.score_pose(peptide.pdbqt, metrics)
                return [self._success_result(job, receptor, peptide, scored,
                    center, size, box_source, 1, _package_version("vina"),
                    meeko_version(),)]
            finally:
                if scorer is not None:
                    scorer.close()
                gc.collect()


def score_group_safely(
    group: ScoreGroup, config: WorkflowConfig) -> list[dict[str, Any]]:
    """Process-pool entry point that turns per-group failures into
    result rows."""
    try:
        return VinaScoringWorkflow(config).score_group(group)
    except (
        Exception
    ) as exc:  # Batch jobs must not stop unrelated receptor groups.
        message = f"{type(exc).__name__}: {exc}"
        return [failure_result(job, message, config) for job in group.jobs]
