"""CSV manifest parsing and memory-aware parallel Vina batch execution.
"""

from __future__ import annotations

import csv
import math
import multiprocessing
import os
from collections.abc import Callable
from concurrent.futures import (
    FIRST_COMPLETED,
    Future,
    ProcessPoolExecutor,
    wait,
)
from pathlib import Path
from typing import Any

import psutil
from models import BoxSpec, ScoreGroup, ScoreJob, WorkflowConfig, group_jobs
from workflow import (
    RESULT_FIELDS,
    VinaScoringWorkflow,
    failure_result,
    score_group_safely,
)

MANIFEST_FIELDS = {
    "id", "complex_pdb", "receptor_pdb",
    "peptide_pdb", "peptide_chain", "receptor_chains", "box_center_x",
    "box_center_y", "box_center_z", "box_size_x", "box_size_y", "box_size_z",}


def _manifest_path(value: str, manifest_directory: Path) -> Path | None:
    text = value.strip()
    if not text:
        return None
    path = Path(text)
    if path.is_absolute():
        return path
    return Path(os.path.normpath(manifest_directory / path))


def _manifest_box(row: dict[str, str], row_number: int) -> BoxSpec:
    names = (
        "box_center_x", "box_center_y",
        "box_center_z", "box_size_x", "box_size_y", "box_size_z",)
    values = [(row.get(name) or "").strip() for name in names]
    if not any(values):
        return BoxSpec()
    if not all(values):
        raise ValueError(f"Manifest row {row_number} must provide all six box "
            "values or none")
    try:
        numbers = [float(value) for value in values]
    except ValueError as exc:
        raise ValueError(
            f"Manifest row {row_number} contains a non-numeric box value"
        ) from exc
    if not all(math.isfinite(number) for number in numbers):
        raise ValueError(
            f"Manifest row {row_number} contains a non-finite box value")
    return BoxSpec(tuple(numbers[:3]), tuple(numbers[3:]))


def read_manifest(filename: Path) -> list[ScoreJob]:
    """Read and validate a CSV manifest with paths relative to the
    manifest."""
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
            raise ValueError(
                f"Unknown manifest column(s): {', '.join(unknown)}")
        if "id" not in headers:
            raise ValueError("Manifest requires an 'id' column")

        for row_number, raw_row in enumerate(reader, start=2):
            row = {
                (key or "").strip(): value or ""
                for key, value in raw_row.items()}
            job_id = row.get("id", "").strip()
            if not job_id:
                raise ValueError(f"Manifest row {row_number} has an empty id")
            if job_id in seen_ids:
                raise ValueError(f"Manifest contains duplicate id {job_id!r}")
            seen_ids.add(job_id)
            base = filename.parent
            complex_pdb = _manifest_path(row.get("complex_pdb", ""), base)
            receptor_pdb = _manifest_path(row.get("receptor_pdb", ""), base)
            peptide_pdb = _manifest_path(row.get("peptide_pdb", ""), base)
            receptor_chains = tuple(value.strip()
                for value in row.get("receptor_chains", "").split(";")
                if value.strip())
            job = ScoreJob(
                index=len(jobs), job_id=job_id, complex_pdb=complex_pdb,
                receptor_pdb=receptor_pdb, peptide_pdb=peptide_pdb,
                peptide_chain=row.get("peptide_chain", "").strip() or None,
                receptor_chains=receptor_chains,
                box=_manifest_box(row, row_number),)
            job.validate()
            jobs.append(job)
    if not jobs:
        raise ValueError("Manifest contains no jobs")
    return jobs


def default_memory_budget_bytes() -> int:
    """Reserve 25% of currently available RAM for the OS and other work.
    """
    return int(psutil.virtual_memory().available * 0.75)


class MemoryAwareBatchRunner:
    """Run complexes in disposable parallel worker processes.

    Each worker handles one complex before exiting so the
    operating system reclaims native Vina map allocations.
    """

    def __init__(self, config: WorkflowConfig,
        *, max_workers: int, memory_budget_bytes: int | None = None,) -> None:
        if max_workers < 1:
            raise ValueError("Worker count must be at least 1")
        self.config = config
        self.max_workers = max_workers
        self.memory_budget_bytes = (
            memory_budget_bytes if memory_budget_bytes is not None
            else default_memory_budget_bytes())
        if self.memory_budget_bytes <= 0:
            raise ValueError("Memory budget must be positive")

    def _estimated_groups(self, groups: list[ScoreGroup]
    ) -> tuple[list[ScoreGroup], list[dict[str, Any]]]:
        workflow = VinaScoringWorkflow(self.config)
        runnable: list[ScoreGroup] = []
        failures: list[dict[str, Any]] = []
        for group in groups:
            try:
                estimate = workflow.estimate_group_memory(group)
            except Exception as exc:
                message = f"{type(exc).__name__}: {exc}"
                failures.extend(failure_result(job, message, self.config)
                    for job in group.jobs)
                continue
            estimated_group = group.with_memory_estimate(estimate)
            if estimate > self.memory_budget_bytes:
                gib = estimate / 1024**3
                budget_gib = self.memory_budget_bytes / 1024**3
                message = (
                    f"Estimated memory {gib:.2f} GiB exceeds the batch budget "
                    f"of {budget_gib:.2f} GiB")
                failures.extend(failure_result(job, message, self.config)
                    for job in group.jobs)
            else:
                runnable.append(estimated_group)
        return runnable, failures

    def run(self, jobs: list[ScoreJob],
        on_results: Callable[[list[dict[str, Any]]], None] | None = None,
    ) -> list[dict[str, Any]]:
        """Score all jobs, optionally reporting each completed result
        group."""
        groups, results = self._estimated_groups(group_jobs(jobs))
        if on_results is not None and results:
            on_results(list(results))
        if not groups:
            return sorted(
                results, key=lambda record: jobs_by_id(jobs)[record["id"]])

        pending = list(groups)
        running: dict[Future[list[dict[str, Any]]], ScoreGroup] = {}
        reserved = 0
        context = multiprocessing.get_context("spawn")
        worker_count = min(self.max_workers, len(groups))
        with ProcessPoolExecutor(
            max_workers=worker_count, mp_context=context,
            max_tasks_per_child=1,
        ) as executor:
            while pending or running:
                submitted = False
                for group in list(pending):
                    if len(running) >= worker_count:
                        break
                    available_reservation = self.memory_budget_bytes - reserved
                    if group.estimated_memory_bytes > available_reservation:
                        continue
                    future = executor.submit(
                        score_group_safely, group, self.config)
                    running[future] = group
                    pending.remove(group)
                    reserved += group.estimated_memory_bytes
                    submitted = True

                if running:
                    # Wait after filling the available worker or memory
                    # slots.
                    done, _ = wait(running, return_when=FIRST_COMPLETED)
                    for future in done:
                        group = running.pop(future)
                        reserved -= group.estimated_memory_bytes
                        try:
                            completed_results = future.result()
                        except Exception as exc:
                            message = (
                                f"Worker failure: {type(exc).__name__}: {exc}")
                            completed_results = [
                                failure_result(job, message, self.config)
                                for job in group.jobs]
                        results.extend(completed_results)
                        if on_results is not None:
                            on_results(completed_results)
                elif pending and not submitted:
                    # Oversized groups are removed during preflight, so
                    # reaching this branch indicates an internal
                    # scheduler inconsistency.
                    group = pending.pop(0)
                    completed_results = [
                        failure_result(job,
                            "Internal error: no memory reservation available",
                            self.config,) for job in group.jobs]
                    results.extend(completed_results)
                    if on_results is not None:
                        on_results(completed_results)

        order = jobs_by_id(jobs)
        return sorted(results, key=lambda record: order[record["id"]])


def jobs_by_id(jobs: list[ScoreJob]) -> dict[str, int]:
    """Return manifest ordering for deterministic parallel output."""
    return {job.job_id: job.index for job in jobs}


def write_results_csv(filename: Path, results: list[dict[str, Any]]) -> None:
    """Write the batch summary without silently replacing an existing
    file."""
    if filename.exists():
        raise ValueError(f"Output already exists: {filename}")
    filename.parent.mkdir(parents=True, exist_ok=True)
    with filename.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=RESULT_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(results)
