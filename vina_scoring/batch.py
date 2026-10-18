"""Manifest parsing and memory-aware disposable-worker scheduling."""

from __future__ import annotations

import csv
import math
import multiprocessing
import os
from collections.abc import Callable
from concurrent.futures import FIRST_COMPLETED, Future, ProcessPoolExecutor, wait
from pathlib import Path

import psutil

from .models import (
    BoxSpec,
    CombinedStructure,
    ScoreJob,
    SeparateStructures,
    WorkflowConfig,
)
from .results import RESULT_FIELDS, ScoreOutcome, failure, result_record
from .workflow import ScoringService, score_job_safely

MANIFEST_FIELDS = {
    "id", "complex_pdb", "receptor_pdb", "peptide_pdb", "peptide_chain",
    "receptor_chains", "box_center_x", "box_center_y", "box_center_z",
    "box_size_x", "box_size_y", "box_size_z",
}
BOX_FIELDS = (
    "box_center_x", "box_center_y", "box_center_z",
    "box_size_x", "box_size_y", "box_size_z",
)


def _manifest_path(value: str, manifest_directory: Path) -> Path | None:
    text = value.strip()
    if not text:
        return None
    path = Path(text)
    return (path if path.is_absolute()
            else Path(os.path.normpath(manifest_directory / path)))


def _manifest_box(row: dict[str, str], row_number: int) -> BoxSpec:
    values = [(row.get(name) or "").strip() for name in BOX_FIELDS]
    if not any(values):
        return BoxSpec()
    if not all(values):
        raise ValueError(
            f"Manifest row {row_number} must provide all six box values or none")
    try:
        numbers = tuple(float(value) for value in values)
    except ValueError as exc:
        raise ValueError(
            f"Manifest row {row_number} contains a non-numeric box value") from exc
    if not all(math.isfinite(number) for number in numbers):
        raise ValueError(
            f"Manifest row {row_number} contains a non-finite box value")
    return BoxSpec(numbers[:3], numbers[3:])


def read_manifest(filename: Path) -> list[ScoreJob]:
    """Parse job shapes; missing files are isolated during batch preflight."""
    if not filename.is_file():
        raise ValueError(f"Manifest does not exist: {filename}")
    jobs: list[ScoreJob] = []
    seen_ids: set[str] = set()
    with filename.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"Manifest has no header: {filename}")
        headers = {header.strip() for header in reader.fieldnames}
        unknown = sorted(headers - MANIFEST_FIELDS)
        if unknown:
            raise ValueError(f"Unknown manifest column(s): {', '.join(unknown)}")
        if "id" not in headers:
            raise ValueError("Manifest requires an 'id' column")
        for row_number, raw_row in enumerate(reader, start=2):
            row = {(key or "").strip(): value or ""
                   for key, value in raw_row.items()}
            job_id = row.get("id", "").strip()
            if not job_id:
                raise ValueError(f"Manifest row {row_number} has an empty id")
            if job_id in seen_ids:
                raise ValueError(f"Manifest contains duplicate id {job_id!r}")
            seen_ids.add(job_id)
            base = filename.parent
            complex_path = _manifest_path(row.get("complex_pdb", ""), base)
            receptor_path = _manifest_path(row.get("receptor_pdb", ""), base)
            peptide_path = _manifest_path(row.get("peptide_pdb", ""), base)
            peptide_chain = row.get("peptide_chain", "").strip()
            receptor_chains = tuple(value.strip()
                                    for value in row.get("receptor_chains", "").split(";")
                                    if value.strip())
            if complex_path is not None:
                if receptor_path is not None or peptide_path is not None:
                    raise ValueError(
                        f"Manifest row {row_number} mixes input modes")
                if not peptide_chain:
                    raise ValueError(
                        f"Manifest row {row_number} requires peptide_chain")
                source = CombinedStructure(complex_path, peptide_chain,
                                           receptor_chains)
            else:
                if receptor_path is None or peptide_path is None:
                    raise ValueError(
                        f"Manifest row {row_number} requires receptor_pdb "
                        "and peptide_pdb")
                if peptide_chain or receptor_chains:
                    raise ValueError(
                        f"Manifest row {row_number} cannot select chains "
                        "for separate files")
                source = SeparateStructures(receptor_path, peptide_path)
            jobs.append(ScoreJob(len(jobs), job_id, source,
                                 _manifest_box(row, row_number)))
    if not jobs:
        raise ValueError("Manifest contains no jobs")
    return jobs


def default_memory_budget_bytes() -> int:
    """Leave 25 percent of currently available RAM for other work."""
    return int(psutil.virtual_memory().available * 0.75)


class BatchRunner:
    """Reserve memory and score each job in a disposable child process."""

    def __init__(self, config: WorkflowConfig, *, max_workers: int,
                 memory_budget_bytes: int | None = None) -> None:
        if max_workers < 1:
            raise ValueError("Worker count must be at least 1")
        self.config = config
        self.max_workers = max_workers
        self.memory_budget_bytes = (
            default_memory_budget_bytes() if memory_budget_bytes is None
            else memory_budget_bytes)
        if self.memory_budget_bytes <= 0:
            raise ValueError("Memory budget must be positive")

    def _estimated_jobs(self, jobs: list[ScoreJob]
                        ) -> tuple[list[tuple[ScoreJob, int]], list[ScoreOutcome]]:
        service = ScoringService(self.config)
        runnable: list[tuple[ScoreJob, int]] = []
        failures: list[ScoreOutcome] = []
        for job in jobs:
            try:
                estimate = service.estimate_job_memory(job)
                if estimate > self.memory_budget_bytes:
                    gib = estimate / 1024**3
                    budget_gib = self.memory_budget_bytes / 1024**3
                    raise ValueError(
                        f"Estimated memory {gib:.2f} GiB exceeds the batch budget "
                        f"of {budget_gib:.2f} GiB")
            except Exception as exc:  # noqa: BLE001 - isolate each job
                failures.append(failure(job, self.config, exc))
            else:
                runnable.append((job, estimate))
        return runnable, failures

    def run(self, jobs: list[ScoreJob],
            on_result: Callable[[ScoreOutcome], None] | None = None,
            ) -> list[ScoreOutcome]:
        """Return manifest-order outcomes; callbacks follow completion order."""
        pending, results = self._estimated_jobs(jobs)
        if on_result is not None:
            for outcome in results:
                on_result(outcome)
        if not pending:
            return sorted(results, key=lambda outcome: outcome.job.index)

        running: dict[Future[ScoreOutcome], tuple[ScoreJob, int]] = {}
        reserved = 0
        worker_count = min(self.max_workers, len(pending))
        context = multiprocessing.get_context("spawn")
        with ProcessPoolExecutor(max_workers=worker_count, mp_context=context,
                                 max_tasks_per_child=1) as executor:
            while pending or running:
                for job, estimate in list(pending):
                    if len(running) >= worker_count:
                        break
                    if estimate > self.memory_budget_bytes - reserved:
                        continue
                    future = executor.submit(score_job_safely, job, self.config)
                    running[future] = (job, estimate)
                    pending.remove((job, estimate))
                    reserved += estimate
                if not running:
                    # Preflight removed all oversized jobs; this is an invariant.
                    raise RuntimeError(
                        "No memory reservation available for a pending job")
                done, _ = wait(running, return_when=FIRST_COMPLETED)
                for future in done:
                    job, estimate = running.pop(future)
                    reserved -= estimate
                    try:
                        outcome = future.result()
                    except Exception as exc:  # noqa: BLE001 - child process failure
                        outcome = failure(job, self.config,
                                          RuntimeError(f"Worker failure: {exc}"))
                    results.append(outcome)
                    if on_result is not None:
                        on_result(outcome)
        return sorted(results, key=lambda outcome: outcome.job.index)


def write_results_csv(filename: Path, results: list[ScoreOutcome]) -> None:
    """Write result rows without replacing an existing output."""
    if filename.exists():
        raise ValueError(f"Output already exists: {filename}")
    filename.parent.mkdir(parents=True, exist_ok=True)
    with filename.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=RESULT_FIELDS)
        writer.writeheader()
        writer.writerows(result_record(outcome) for outcome in results)
