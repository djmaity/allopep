"""Typed score outcomes and the single supported result schema."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from .models import ScoreJob, WorkflowConfig
from .preparation import installed_version
from .vina import SCORE_TERMS, PoseMetrics, size_normalized_diagnostics

PROTOCOL_ID = "rigid_meeko_vina_residue_sum_v1"
NORMALIZED_TERMS = (
    "total_per_heavy_atom", "interaction_per_heavy_atom",
    "total_per_residue", "interaction_per_residue",
)
RESULT_FIELDS = (
    "protocol_id", "id", "status", "error",
    "complex_pdb", "complex_sha256", "model_index", "peptide_chain",
    "receptor_chains", "hydrogen_policy", "ph", "included_hetero_residues",
    "warnings", "box_source", "box_center_x", "box_center_y",
    "box_center_z", "box_size_x", "box_size_y", "box_size_z", "spacing",
    "peptide_atom_count", "peptide_heavy_atom_count",
    "peptide_residue_count", "peptide_active_torsions", "score_units",
    *(f"score_{term}" for term in SCORE_TERMS),
    *(f"normalized_{term}" for term in NORMALIZED_TERMS),
    "vina_version", "meeko_version",
)


def hetero_policy_description(config: WorkflowConfig) -> str:
    """Describe the role-aware hetero-residue policy."""
    if config.include_all_hetero:
        receptor_policy = "all"
    elif config.keep_hetero_residues:
        receptor_policy = ",".join(sorted(config.keep_hetero_residues))
    else:
        receptor_policy = "none"
    return f"peptide:all_non_water;receptor:{receptor_policy}"


@dataclass(frozen=True)
class ScoreSuccess:
    """A fully prepared and scored rigid pose."""

    job: ScoreJob
    config: WorkflowConfig
    receptor_metadata: dict[str, Any]
    peptide_metadata: dict[str, Any]
    metrics: PoseMetrics
    scores: dict[str, float]
    center: tuple[float, float, float]
    size: tuple[float, float, float]
    box_source: str
    vina_version: str
    meeko_version: str


@dataclass(frozen=True)
class ScoreFailure:
    """An isolated job failure that can be serialized beside successes."""

    job: ScoreJob
    config: WorkflowConfig
    error: str


ScoreOutcome = ScoreSuccess | ScoreFailure


def failure(job: ScoreJob, config: WorkflowConfig, exc: Exception) -> ScoreFailure:
    """Keep exception type in a per-job diagnostic."""
    return ScoreFailure(job, config, f"{type(exc).__name__}: {exc}")


def result_record(outcome: ScoreOutcome) -> dict[str, Any]:
    """Serialize either outcome with identical fields and ordering."""
    job, config = outcome.job, outcome.config
    source = job.source
    record: dict[str, Any] = {name: "" for name in RESULT_FIELDS}
    record.update({
        "protocol_id": PROTOCOL_ID,
        "id": job.job_id,
        "status": "error" if isinstance(outcome, ScoreFailure) else "ok",
        "complex_pdb": str(source.path),
        "model_index": config.model_index,
        "peptide_chain": source.peptide_chain,
        "receptor_chains": ";".join(source.receptor_chains),
        "hydrogen_policy": config.hydrogen_policy,
        "ph": config.ph,
        "included_hetero_residues": hetero_policy_description(config),
        "spacing": config.spacing,
        "vina_version": installed_version("vina"),
        "meeko_version": installed_version("meeko"),
    })
    if isinstance(outcome, ScoreFailure):
        record["error"] = outcome.error
        return record

    receptor_metadata = outcome.receptor_metadata
    peptide_metadata = outcome.peptide_metadata
    warnings = list(dict.fromkeys(
        receptor_metadata["warnings"] + peptide_metadata["warnings"]))
    if not source.receptor_chains:
        warnings.insert(0, "Receptor chains were inferred as every non-peptide chain.")
    record.update({
        "receptor_chains": ";".join(receptor_metadata["chains"]),
        "warnings": " | ".join(warnings),
        "box_source": outcome.box_source,
        "box_center_x": outcome.center[0],
        "box_center_y": outcome.center[1],
        "box_center_z": outcome.center[2],
        "box_size_x": outcome.size[0],
        "box_size_y": outcome.size[1],
        "box_size_z": outcome.size[2],
        "score_units": "nominal_kcal_per_mol",
        "vina_version": outcome.vina_version,
        "meeko_version": outcome.meeko_version,
    })
    record["complex_sha256"] = receptor_metadata["source_sha256"]
    metrics = outcome.metrics
    record.update({
        "peptide_atom_count": metrics.atom_count,
        "peptide_heavy_atom_count": metrics.heavy_atom_count,
        "peptide_residue_count": metrics.residue_count,
        "peptide_active_torsions": metrics.active_torsions,
    })
    record.update({f"score_{name}": value
                   for name, value in outcome.scores.items()})
    record.update({f"normalized_{name}": value
                   for name, value in size_normalized_diagnostics(
                       outcome.scores, metrics.heavy_atom_count,
                       metrics.residue_count).items()})
    return record
