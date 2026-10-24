"""Score a supplied rigid protein-peptide complex with AutoDock Vina."""

from __future__ import annotations

import argparse
import json
import os
import sys
from collections.abc import Sequence
from pathlib import Path

from .models import CombinedStructure, ScoreJob, WorkflowConfig


def argument_parser(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse one combined PDB/mmCIF complex and its preparation settings."""
    parser = argparse.ArgumentParser(
        prog="python -m allopep.scoring", description=(
            "Score supplied protein-peptide coordinates with Vina. Receptor "
            "and peptide coordinates remain rigid; docking is not performed."
        ),)
    parser.add_argument(
        "complex_pdb", type=Path, metavar="COMPLEX_STRUCTURE",
        help="PDB or mmCIF containing both receptor and peptide",)
    parser.add_argument("-c", "--chain", dest="peptide_chain", required=True,
        help="Peptide chain in the combined structure",)
    parser.add_argument(
        "-r", "--receptor-chain", dest="receptor_chains", action="append",
        help=("Receptor chain; repeat for multiple chains. The default is "
              "every non-peptide chain."),)

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
            "templates (default: keep)"),)
    preparation.add_argument("--ph", type=float, default=7.4, help=(
            "Intended protonation pH recorded as provenance; Meeko "
            "protonation is controlled by residue templates (default: 7.4)"),)
    preparation.add_argument("--meeko-templates", type=Path,
        help="Additional cached/custom Meeko residue-template JSON",)

    box = parser.add_argument_group("affinity-map box")
    box.add_argument("--box-padding", type=float, default=5.0,
        help=("Padding on every side of the peptide's own bounds "
              "(default: 5 A)"),)
    box.add_argument("--spacing",
        type=float, default=0.375, help="Vina map spacing (default: 0.375 A)",)

    parser.add_argument("-o", "--output", type=Path,
        help="JSON output path; the default is standard output",)
    parser.add_argument("--compact", action="store_true",
        help="Use compact JSON for standard output",)
    return parser.parse_args(argv)


def _config(args: argparse.Namespace) -> WorkflowConfig:
    return WorkflowConfig(model_index=args.model_index,
        keep_hetero_residues=tuple(args.keep_hetero_residue),
        include_all_hetero=args.include_all_hetero, ph=args.ph,
        hydrogen_policy=args.hydrogen_policy, box_padding=args.box_padding,
        spacing=args.spacing, meeko_templates=args.meeko_templates,)


def _write_json(record: dict[str, object], output: Path | None,
                compact: bool) -> None:
    serialized = json.dumps(
        record, indent=None if compact else 2, sort_keys=True,)
    if output is None:
        print(serialized)
        return
    if output.exists():
        raise ValueError(f"Output already exists: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(serialized + "\n", encoding="utf-8")


def main(argv: Sequence[str] | None = None) -> int:
    """Score one rigid pose and report it as JSON."""
    args = argument_parser(argv)
    try:
        # Set thread limits before loading NumPy or native scoring libraries.
        os.environ.update({"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
                           "MKL_NUM_THREADS": "1", "NUMEXPR_NUM_THREADS": "1"})
        from .records import ScoreFailure, result_record
        from .service import score_job_safely

        if args.output is not None and args.output.exists():
            raise ValueError(f"Output already exists: {args.output}")
        job = ScoreJob(
            index=0, job_id=args.complex_pdb.stem,
            source=CombinedStructure(args.complex_pdb, args.peptide_chain,
                                     tuple(args.receptor_chains or ())))
        outcome = score_job_safely(job, _config(args))
        _write_json(result_record(outcome), args.output, args.compact)
        if isinstance(outcome, ScoreFailure):
            return 1
    except (OSError, RuntimeError, ValueError) as exc:
        print(f"allopep.scoring: error: {exc}", file=sys.stderr)
        return 1
    return 0
