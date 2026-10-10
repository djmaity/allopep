"""Data models shared by the Vina preparation and batch workflows."""

from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path


@dataclass(frozen=True)
class BoxSpec:
    """An optional explicit Vina affinity-map box."""

    center: tuple[float, float, float] | None = None
    size: tuple[float, float, float] | None = None

    def __post_init__(self) -> None:
        if (self.center is None) != (self.size is None):
            raise ValueError("Box center and size must be supplied together")

    @property
    def is_explicit(self) -> bool:
        """Return whether fixed box coordinates were supplied."""
        return self.center is not None


@dataclass(frozen=True)
class ScoreJob:
    """One protein-peptide PDB/mmCIF pose requested by the CLI or
    manifest."""

    index: int
    job_id: str
    complex_pdb: Path | None = None
    receptor_pdb: Path | None = None
    peptide_pdb: Path | None = None
    peptide_chain: str | None = None
    receptor_chains: tuple[str, ...] = ()
    box: BoxSpec = BoxSpec()

    @property
    def input_mode(self) -> str:
        """Return ``combined`` or ``separate`` for output metadata."""
        return "combined" if self.complex_pdb is not None else "separate"

    def validate(self) -> None:
        """Validate mutually exclusive structure inputs and required
        files."""
        if self.complex_pdb is not None:
            if self.receptor_pdb is not None or self.peptide_pdb is not None:
                raise ValueError(
                    f"Job {self.job_id!r} mixes combined and separate inputs")
            if not self.peptide_chain:
                raise ValueError(
                    f"Job {self.job_id!r} requires a peptide chain "
                    "for its complex PDB")
            if self.peptide_chain in self.receptor_chains:
                raise ValueError(
                    f"Job {self.job_id!r} selects its peptide chain "
                    "as receptor")
            paths = (self.complex_pdb,)
        else:
            if self.receptor_pdb is None or self.peptide_pdb is None:
                raise ValueError(f"Job {self.job_id!r} requires receptor and "
                    "peptide structures")
            if self.peptide_chain or self.receptor_chains:
                raise ValueError(
                    f"Job {self.job_id!r} cannot use chain selection "
                    "with separate PDBs")
            paths = (self.receptor_pdb, self.peptide_pdb)

        for path in paths:
            if not path.is_file():
                raise ValueError(
                    f"Job {self.job_id!r} input does not exist: {path}")

    def group_key(self) -> tuple[object, ...]:
        """Return a key for jobs that may share one receptor map
        calculation."""
        box_key = (self.box.center, self.box.size)
        if self.input_mode == "separate":
            assert self.receptor_pdb is not None
            return ("separate", str(self.receptor_pdb.resolve()), box_key)
        # Combined PDBs are deliberately independent: their receptor
        # extraction depends on the peptide and receptor chain
        # selections in that same file.
        return ("combined", self.index)


@dataclass(frozen=True)
class ScoreGroup:
    """One job assigned to a disposable Vina worker."""

    jobs: tuple[ScoreJob, ...]
    estimated_memory_bytes: int = 0

    def with_memory_estimate(self, memory_bytes: int) -> "ScoreGroup":
        """Return a copy carrying its scheduler memory reservation."""
        return replace(self, estimated_memory_bytes=memory_bytes)


@dataclass(frozen=True)
class WorkflowConfig:
    """Preparation and scoring settings common to all jobs in one
    invocation."""

    model_index: int = 0
    keep_hetero_residues: tuple[str, ...] = ()
    include_all_hetero: bool = False
    ph: float = 7.4
    hydrogen_policy: str = "keep"
    box_padding: float = 5.0
    spacing: float = 0.375
    meeko_templates: Path | None = None


def group_jobs(jobs: list[ScoreJob]) -> list[ScoreGroup]:
    """Return one isolated score group per input job.

    Vina owns native affinity-map allocations that are most reliably
    returned to the operating system when the worker exits. Keeping every
    job in its own group ensures a disposable worker handles exactly one
    complex.
    """
    return [ScoreGroup((job,)) for job in jobs]
