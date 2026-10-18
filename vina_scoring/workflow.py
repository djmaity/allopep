"""Coordinate one rigid preparation and scoring job."""

from __future__ import annotations

import tempfile
from pathlib import Path

from .models import CombinedStructure, ScoreJob, WorkflowConfig
from .preparation import PdbPreparer, PreparedComponent, meeko_version
from .results import ScoreOutcome, ScoreSuccess, failure, installed_version
from .scoring import (
    VinaPoseScorer,
    box_from_coordinates,
    estimate_vina_memory_bytes,
    parse_pdbqt,
    resolve_box,
)


class ScoringService:
    """Select, prepare, validate, and score one supplied complex."""

    def __init__(self, config: WorkflowConfig) -> None:
        self.config = config
        self.preparer = PdbPreparer(
            model_index=config.model_index,
            keep_hetero_residues=config.keep_hetero_residues,
            include_all_hetero=config.include_all_hetero,
            ph=config.ph,
            hydrogen_policy=config.hydrogen_policy,
            meeko_templates=config.meeko_templates,
        )

    def estimate_job_memory(self, job: ScoreJob) -> int:
        """Reserve for a source-coordinate box before launching Vina."""
        job.validate_files()
        source = job.source
        if isinstance(source, CombinedStructure):
            coordinates = self.preparer.selected_coordinates(
                source.path, (source.peptide_chain,), as_ligand=True)
        else:
            coordinates = self.preparer.selected_coordinates(
                source.peptide, as_ligand=True)
        if job.box.is_explicit:
            assert job.box.size is not None
            size = list(job.box.size)
        else:
            _, size = box_from_coordinates([coordinates], self.config.box_padding)
            # Rebuilt hydrogens can extend beyond source heavy-atom bounds.
            size = [dimension + 2.0 for dimension in size]
        return estimate_vina_memory_bytes(size, self.config.spacing)

    def _receptor_chains(self, source: CombinedStructure) -> tuple[str, ...]:
        if source.receptor_chains:
            return source.receptor_chains
        return tuple(chain for chain in self.preparer.available_chains(source.path)
                     if chain != source.peptide_chain)

    def _prepare_receptor(self, job: ScoreJob, directory: Path) -> PreparedComponent:
        source = job.source
        if isinstance(source, CombinedStructure):
            chains = self._receptor_chains(source)
            if not chains:
                raise ValueError(f"Job {job.job_id!r} selected no receptor chains")
            return self.preparer.prepare_component(source.path, directory, chains)
        return self.preparer.prepare_component(source.receptor, directory)

    def _prepare_peptide(self, job: ScoreJob, directory: Path) -> PreparedComponent:
        source = job.source
        if isinstance(source, CombinedStructure):
            return self.preparer.prepare_component(
                source.path, directory, (source.peptide_chain,), as_ligand=True)
        return self.preparer.prepare_component(
            source.peptide, directory, as_ligand=True)

    def score_job(self, job: ScoreJob) -> ScoreSuccess:
        """Score a pose without docking or changing supplied coordinates."""
        job.validate_files()
        with tempfile.TemporaryDirectory(prefix="vina-rigid-score-") as temporary:
            directory = Path(temporary)
            receptor = self._prepare_receptor(job, directory / "receptor")
            peptide = self._prepare_peptide(job, directory / "peptide")
            metrics = parse_pdbqt(peptide.pdbqt)
            center, size, box_source = resolve_box(
                [metrics], job.box, self.config.box_padding)
            scorer: VinaPoseScorer | None = None
            try:
                scorer = VinaPoseScorer(receptor.pdbqt, center, size,
                                        spacing=self.config.spacing, cpu=1)
                scores = scorer.score_pose(peptide.pdbqt, metrics)
                return ScoreSuccess(
                    job=job, config=self.config,
                    receptor_metadata=receptor.metadata,
                    peptide_metadata=peptide.metadata,
                    metrics=metrics, scores=scores,
                    center=tuple(center), size=tuple(size),
                    box_source=box_source,
                    vina_version=installed_version("vina"),
                    meeko_version=meeko_version() or "",
                )
            finally:
                if scorer is not None:
                    scorer.close()


def score_job_safely(job: ScoreJob, config: WorkflowConfig) -> ScoreOutcome:
    """Turn a worker exception into a result for that one job."""
    try:
        return ScoringService(config).score_job(job)
    except Exception as exc:  # noqa: BLE001 - isolate independent jobs
        return failure(job, config, exc)
