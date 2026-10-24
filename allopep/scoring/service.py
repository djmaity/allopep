"""Prepare, score and schedule rigid-pose scoring jobs.

Every job runs in a disposable child process. Vina allocates native
affinity maps that it does not fully release, so a fresh interpreter per
job is what keeps a long batch's memory flat.
"""

from __future__ import annotations

import multiprocessing
import tempfile
from collections.abc import Callable
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import psutil

from .models import ScoreJob, WorkflowConfig
from .preparation import MeekoPreparer, PreparedComponent, installed_version
from .records import ScoreOutcome, ScoreSuccess, failure
from .vina import (
    VinaPoseScorer,
    box_from_coordinates,
    estimate_vina_memory_bytes,
    parse_pdbqt,
    resolve_box,
)


def default_memory_budget_bytes() -> int:
    """Leave 25 percent of currently available RAM for other work."""
    return int(psutil.virtual_memory().available * 0.75)


class ScoringService:
    """Select, prepare, validate, and score one supplied complex."""

    def __init__(self, config: WorkflowConfig) -> None:
        self.config = config
        self.preparer = MeekoPreparer(
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
        coordinates = self.preparer.selector.selected_coordinates(
            source.path, (source.peptide_chain,), as_ligand=True)
        _, size = box_from_coordinates([coordinates], self.config.box_padding)
        # Rebuilt hydrogens can extend beyond source heavy-atom bounds.
        size = [dimension + 2.0 for dimension in size]
        return estimate_vina_memory_bytes(size, self.config.spacing)

    def _receptor_chains(self, job: ScoreJob) -> tuple[str, ...]:
        source = job.source
        if source.receptor_chains:
            return source.receptor_chains
        chains = tuple(
            chain for chain in self.preparer.selector.available_chains(source.path)
            if chain != source.peptide_chain)
        if not chains:
            raise ValueError(f"Job {job.job_id!r} selected no receptor chains")
        return chains

    def _prepare_receptor(self, job: ScoreJob, directory: Path) -> PreparedComponent:
        return self.preparer.prepare_component(
            job.source.path, directory, self._receptor_chains(job))

    def _prepare_peptide(self, job: ScoreJob, directory: Path) -> PreparedComponent:
        return self.preparer.prepare_component(
            job.source.path, directory, (job.source.peptide_chain,),
            as_ligand=True)

    def score_job(self, job: ScoreJob) -> ScoreSuccess:
        """Score a pose without docking or changing supplied coordinates."""
        job.validate_files()
        with tempfile.TemporaryDirectory(prefix="vina-rigid-score-") as temporary:
            directory = Path(temporary)
            receptor = self._prepare_receptor(job, directory / "receptor")
            peptide = self._prepare_peptide(job, directory / "peptide")
            metrics = parse_pdbqt(peptide.pdbqt)
            center, size, box_source = resolve_box(
                [metrics], self.config.box_padding)
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
                    meeko_version=installed_version("meeko"),
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


class BatchRunner:
    """Score each job in a disposable child process, one at a time."""

    def __init__(self, config: WorkflowConfig, *,
                 memory_budget_bytes: int | None = None) -> None:
        self.config = config
        self.memory_budget_bytes = (
            default_memory_budget_bytes() if memory_budget_bytes is None
            else memory_budget_bytes)
        if self.memory_budget_bytes <= 0:
            raise ValueError("Memory budget must be positive")

    def _preflight(self, job: ScoreJob) -> ScoreOutcome | None:
        """Reject a job whose maps cannot fit, instead of being killed."""
        try:
            estimate = ScoringService(self.config).estimate_job_memory(job)
            if estimate > self.memory_budget_bytes:
                raise ValueError(
                    f"Estimated memory {estimate / 1024**3:.2f} GiB exceeds the "
                    f"batch budget of {self.memory_budget_bytes / 1024**3:.2f} GiB")
        except Exception as exc:  # noqa: BLE001 - isolate each job
            return failure(job, self.config, exc)
        return None

    def run(self, jobs: list[ScoreJob],
            on_result: Callable[[ScoreOutcome], None] | None = None,
            ) -> list[ScoreOutcome]:
        """Return outcomes in job-index order."""
        results: list[ScoreOutcome] = []
        runnable: list[ScoreJob] = []
        for job in jobs:
            rejected = self._preflight(job)
            if rejected is None:
                runnable.append(job)
            else:
                results.append(rejected)
                if on_result is not None:
                    on_result(rejected)
        if runnable:
            context = multiprocessing.get_context("spawn")
            with ProcessPoolExecutor(max_workers=1, mp_context=context,
                                     max_tasks_per_child=1) as executor:
                futures = [(job, executor.submit(score_job_safely, job, self.config))
                           for job in runnable]
                for job, future in futures:
                    try:
                        outcome = future.result()
                    except Exception as exc:  # noqa: BLE001 - child failure
                        outcome = failure(job, self.config,
                                          RuntimeError(f"Worker failure: {exc}"))
                    results.append(outcome)
                    if on_result is not None:
                        on_result(outcome)
        return sorted(results, key=lambda outcome: outcome.job.index)
