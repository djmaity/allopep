"""Prepare and score rigid protein-peptide poses with AutoDock Vina."""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from collections.abc import Sequence
from pathlib import Path

# Worker count must be controlled by this CLI rather than nested numerical
# libraries that may otherwise create large thread-local allocations.
os.environ.update({"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1", "NUMEXPR_NUM_THREADS": "1",})

from batch import MemoryAwareBatchRunner, read_manifest, write_results_csv
from models import BoxSpec, ScoreGroup, ScoreJob, WorkflowConfig
from workflow import score_group_safely


def argument_parser(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse combined, separate, or CSV-manifest PDB/mmCIF inputs."""
    parser = argparse.ArgumentParser(prog="vina_score.py", description=(
            "Score supplied protein-peptide coordinates with Vina. Receptor "
            "and peptide coordinates remain rigid; docking is not performed."
        ),)
    parser.add_argument(
        "complex_pdb", nargs="?", type=Path, metavar="COMPLEX_STRUCTURE",
        help="PDB or mmCIF containing both receptor and peptide",)
    parser.add_argument("-c", "--chain",
        dest="peptide_chain", help="Peptide chain in a combined structure",)
    parser.add_argument(
        "-r", "--receptor-chain", dest="receptor_chains", action="append",
        help=("Receptor chain in a combined structure; repeat for "
            "multiple chains"),)

    separate = parser.add_argument_group("separate structure input")
    separate.add_argument("--receptor-pdb", "--receptor-structure",
        dest="receptor_pdb", type=Path, help="Receptor-only PDB or mmCIF",)
    separate.add_argument("--peptide-pdb",
        "--peptide-structure", dest="peptide_pdb", type=Path, nargs="+", help=(
            "One or more peptide-only PDB/mmCIF poses. They are scored "
            "independently "
            "and must use the same coordinate frame as the receptor."),)

    batch = parser.add_argument_group("batch input")
    batch.add_argument("--manifest", type=Path,
        help="CSV manifest containing combined or separate structure jobs",)
    batch.add_argument(
        "--workers", type=int, default=max(1, os.cpu_count() or 1),
        help="Maximum parallel receptor groups (default: available CPUs)",)
    batch.add_argument("--memory-budget-gb", type=float,
        help="RAM available to workers; default is 75%% of available memory",)

    preparation = parser.add_argument_group("structure preparation")
    preparation.add_argument("--model-index", type=int,
        default=0, help="Zero-based coordinate model index (default: 0)",)
    hetero = preparation.add_mutually_exclusive_group()
    hetero.add_argument("--keep-hetero-residue",
        action="append", default=[], metavar="RESNAME",
        help="Retain a receptor HETATM residue name; repeat as needed",)
    hetero.add_argument("--include-all-hetero", action="store_true",
        help="Retain every receptor HETATM residue, including waters",)
    preparation.add_argument("--hydrogen-policy",
        choices=("keep", "rebuild"), default="keep", help=(
            "Keep input hydrogens or rebuild them from Meeko residue "
            "templates " "(default: keep)"),)
    preparation.add_argument("--ph", type=float, default=7.4, help=(
            "Intended protonation pH recorded as provenance; Meeko "
            "protonation "
            "is controlled by residue templates (default: 7.4)"),)
    preparation.add_argument("--meeko-templates", type=Path,
        help="Additional cached/custom Meeko residue-template JSON",)

    box = parser.add_argument_group("affinity-map box")
    box.add_argument("--box-padding", type=float, default=5.0,
        help=("Padding on every side of automatic peptide bounds "
            "(default: 5 A)"),)
    box.add_argument(
        "--box-center", nargs=3, type=float, metavar=("X", "Y", "Z"),
        help="Fixed box center for a single-input invocation",)
    box.add_argument(
        "--box-size", nargs=3, type=float, metavar=("X", "Y", "Z"),
        help="Fixed box dimensions; requires --box-center",)
    box.add_argument("--spacing",
        type=float, default=0.375, help="Vina map spacing (default: 0.375 A)",)

    parser.add_argument("-o", "--output", type=Path,
        help=("JSON output for single input or required CSV output for "
            "a manifest"),)
    parser.add_argument("--compact", action="store_true",
        help="Use compact JSON for single-input standard output",)
    return parser.parse_args(argv)


def _box_spec(args: argparse.Namespace) -> BoxSpec:
    if (args.box_center is None) != (args.box_size is None):
        raise ValueError(
            "--box-center and --box-size must be supplied together")
    if args.box_center is None:
        return BoxSpec()
    return BoxSpec(tuple(args.box_center), tuple(args.box_size))


def _config(args: argparse.Namespace) -> WorkflowConfig:
    if args.model_index < 0:
        raise ValueError("--model-index cannot be negative")
    if not math.isfinite(args.ph):
        raise ValueError("--ph must be finite")
    if not math.isfinite(args.box_padding) or args.box_padding < 0:
        raise ValueError("--box-padding must be finite and non-negative")
    if not math.isfinite(args.spacing) or args.spacing <= 0:
        raise ValueError("--spacing must be finite and positive")
    if args.workers < 1:
        raise ValueError("--workers must be at least 1")
    if args.memory_budget_gb is not None and (
        not math.isfinite(args.memory_budget_gb) or args.memory_budget_gb <= 0
    ):
        raise ValueError("--memory-budget-gb must be finite and positive")
    if args.meeko_templates is not None and not args.meeko_templates.is_file():
        raise ValueError(
            f"--meeko-templates does not exist: {args.meeko_templates}")
    return WorkflowConfig(model_index=args.model_index,
        keep_hetero_residues=tuple(args.keep_hetero_residue),
        include_all_hetero=args.include_all_hetero, ph=args.ph,
        hydrogen_policy=args.hydrogen_policy, box_padding=args.box_padding,
        spacing=args.spacing, meeko_templates=args.meeko_templates,)


def _single_jobs(args: argparse.Namespace) -> list[ScoreJob]:
    separate_requested = (
        args.receptor_pdb is not None or args.peptide_pdb is not None)
    if args.complex_pdb is not None and separate_requested:
        raise ValueError(
            "Do not combine a complex structure with separate inputs")
    if args.complex_pdb is not None:
        if not args.peptide_chain:
            raise ValueError("--chain is required for a combined complex PDB")
        return [
            ScoreJob(index=0, job_id=args.complex_pdb.stem,
                complex_pdb=args.complex_pdb, peptide_chain=args.peptide_chain,
                receptor_chains=tuple(args.receptor_chains or ()),
                box=_box_spec(args),)]
    if separate_requested:
        if args.receptor_pdb is None or args.peptide_pdb is None:
            raise ValueError(
                "--receptor-pdb and --peptide-pdb must be supplied together")
        if args.peptide_chain or args.receptor_chains:
            raise ValueError("Chain selection applies only to a combined PDB")
        return [
            ScoreJob(
                index=index, job_id=(peptide.stem if len(args.peptide_pdb) > 1
                    else f"{args.receptor_pdb.stem}__{peptide.stem}"
                ), receptor_pdb=args.receptor_pdb,
                peptide_pdb=peptide, box=_box_spec(args),
            ) for index, peptide in enumerate(args.peptide_pdb)]
    raise ValueError(
        "Supply a complex structure, separate receptor/peptide structures, "
        "or --manifest")


def _validate_and_select_mode(
    args: argparse.Namespace,) -> tuple[str, list[ScoreJob]]:
    if args.manifest is not None:
        if any(value is not None for value in (
                args.complex_pdb, args.receptor_pdb, args.peptide_pdb,)):
            raise ValueError(
                "Do not combine --manifest with single-input PDBs")
        if args.peptide_chain or args.receptor_chains:
            raise ValueError(
                "Manifest rows provide combined-PDB chain selection")
        if args.box_center is not None or args.box_size is not None:
            raise ValueError(
                "Manifest rows provide any explicit box coordinates")
        if args.output is None:
            raise ValueError("--output is required with --manifest")
        return "batch", read_manifest(args.manifest)
    return "single", _single_jobs(args)


def _write_single_json(
    results: list[dict[str, object]], output: Path | None, compact: bool
) -> None:
    document: object = (
        results[0] if len(results) == 1 else {"results": results})
    serialized = json.dumps(
        document, indent=None if compact else 2, sort_keys=True,)
    if output is None:
        print(serialized)
        return
    if output.exists():
        raise ValueError(f"Output already exists: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(serialized + "\n", encoding="utf-8")


def main(argv: Sequence[str] | None = None) -> int:
    """Run one rigid-pose score or a memory-aware parallel batch."""
    args = argument_parser(argv)
    try:
        config = _config(args)
        mode, jobs = _validate_and_select_mode(args)
        for job in jobs:
            job.validate()
        if args.output is not None and args.output.exists():
            raise ValueError(f"Output already exists: {args.output}")

        if mode == "batch":
            budget = (
                int(args.memory_budget_gb * 1024**3)
                if args.memory_budget_gb is not None else None)
            results = MemoryAwareBatchRunner(config, max_workers=args.workers,
                memory_budget_bytes=budget,).run(jobs)
            assert args.output is not None
            write_results_csv(args.output, results)
        else:
            results = []
            for job in jobs:
                results.extend(
                    score_group_safely(ScoreGroup((job,)), config))
            _write_single_json(results, args.output, args.compact)
            if any(result["status"] != "ok" for result in results):
                return 1
    except (OSError, RuntimeError, ValueError) as exc:
        print(f"vina_score.py: error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
