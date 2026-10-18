"""Validated inputs and settings for rigid protein-peptide scoring."""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path


@dataclass(frozen=True)
class BoxSpec:
    """An optional explicit Vina affinity-map box."""

    center: tuple[float, float, float] | None = None
    size: tuple[float, float, float] | None = None

    def __post_init__(self) -> None:
        if (self.center is None) != (self.size is None):
            raise ValueError("Box center and size must be supplied together")
        if self.center is not None and self.size is not None:
            if len(self.center) != 3 or len(self.size) != 3:
                raise ValueError("Box center and size must each contain three values")
            if not all(math.isfinite(value) for value in (*self.center, *self.size)):
                raise ValueError("Box center and size must be finite")
            if any(value <= 0 for value in self.size):
                raise ValueError("Box dimensions must be positive")

    @property
    def is_explicit(self) -> bool:
        return self.center is not None


@dataclass(frozen=True)
class CombinedStructure:
    """A complex with a selected peptide chain and receptor chains."""

    path: Path
    peptide_chain: str
    receptor_chains: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if not self.peptide_chain:
            raise ValueError("A combined structure requires a peptide chain")
        if self.peptide_chain in self.receptor_chains:
            raise ValueError("The peptide chain cannot also be a receptor chain")


@dataclass(frozen=True)
class SeparateStructures:
    """Receptor and peptide structures in the same coordinate frame."""

    receptor: Path
    peptide: Path


@dataclass(frozen=True)
class ScoreJob:
    """One independent rigid pose to score."""

    index: int
    job_id: str
    source: CombinedStructure | SeparateStructures
    box: BoxSpec = field(default_factory=BoxSpec)

    def __post_init__(self) -> None:
        if self.index < 0:
            raise ValueError("Job index cannot be negative")
        if not self.job_id:
            raise ValueError("Job ID cannot be empty")
        if not isinstance(self.source, (CombinedStructure, SeparateStructures)):
            raise TypeError("Job source must be a combined or separate structure")

    @property
    def input_mode(self) -> str:
        return "combined" if isinstance(self.source, CombinedStructure) else "separate"

    def validate_files(self) -> None:
        """Check existence in preflight so one missing file fails one job."""
        paths = ((self.source.path,) if isinstance(self.source, CombinedStructure)
                 else (self.source.receptor, self.source.peptide))
        for path in paths:
            if not path.is_file():
                raise ValueError(f"Job {self.job_id!r} input does not exist: {path}")


@dataclass(frozen=True)
class WorkflowConfig:
    """Preparation and scoring settings shared by one invocation."""

    model_index: int = 0
    keep_hetero_residues: tuple[str, ...] = ()
    include_all_hetero: bool = False
    ph: float = 7.4
    hydrogen_policy: str = "keep"
    box_padding: float = 5.0
    spacing: float = 0.375
    meeko_templates: Path | None = None

    def __post_init__(self) -> None:
        if self.model_index < 0:
            raise ValueError("Model index cannot be negative")
        if not math.isfinite(self.ph):
            raise ValueError("pH must be finite")
        if self.hydrogen_policy not in {"keep", "rebuild"}:
            raise ValueError("Hydrogen policy must be 'keep' or 'rebuild'")
        if self.include_all_hetero and self.keep_hetero_residues:
            raise ValueError("Choose all hetero residues or selected names")
        if not math.isfinite(self.box_padding) or self.box_padding < 0:
            raise ValueError("Box padding must be finite and non-negative")
        if not math.isfinite(self.spacing) or self.spacing <= 0:
            raise ValueError("Grid spacing must be finite and positive")
        if self.meeko_templates is not None and not self.meeko_templates.is_file():
            raise ValueError(
                f"Meeko template file does not exist: {self.meeko_templates}")
