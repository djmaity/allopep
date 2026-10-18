# Pocket selection

AlloPep groups high-scoring, solvent-accessible GAPS residues into connected
pockets. The selected residue identities guide PepGLAD's generation; their
numerical GAPS scores are not supplied to the generative model.

## Current selection rules

| Property | Rule |
| --- | --- |
| GAPS score | Strictly **>0.5**. A score exactly equal to 0.5 is excluded. |
| Accessibility | Absolute per-residue **SASA >0 Å²**, calculated with FreeSASA on the original receptor. |
| Missing or non-finite accessibility | Treated as zero and excluded. |
| Coordinates | A mapped Cα atom is required. |
| Connection | Eligible residues connect when their Cα distance is **≤8 Å**, including across chains. |
| Pocket | A connected component containing **at least two residues**. |
| Upper size limit | None; a retained component is passed in full. |
| Ranking | Mean member GAPS score, highest first. |
| Ranking ties | Smallest member identity: chain, residue number, insertion code. |
| Number of pockets used | All by default; `--max_rank N` selects the first N. |

SASA here is an **absolute area**, not relative solvent accessibility. There
is no reference-area normalization or percentage-exposure cutoff.

## Connectivity is transitive

If residues A and B are 7 Å apart, and B and C are 7 Å apart, all three can be
in one pocket even when A and C are more than 8 Å apart. A pocket therefore
need not fit inside an 8 Å sphere.

Only eligible residues participate in connectivity. A residue excluded by the
score or accessibility threshold cannot bridge two components. Every eligible
residue belongs to at most one component, and distinct retained pockets have
no cross-pocket residue pair within 8 Å.

There is no center-radius restriction, neighborhood-score blend, 48-residue cap,
or separate seed-separation rule in this selector. `GuidanceContext.seed` is
the highest-scoring member, with identity used to break ties; it is a reference
member and does not determine the pocket's extent or rank.

## Full-precision scores and residue mapping

GAPS writes `structure_GAPS_scores.json` with schema version `1` and a
`residues` list. Records contain the chain, residue number, insertion code,
residue name, Cα coordinates, and full-precision sigmoid score.

Before selection, AlloPep checks the records against the original receptor:

- Identities must be unique and residue names must match.
- Scores must be finite and between zero and one.
- A supplied Cα coordinate must agree with the receptor within 0.05 Å.
- Residues without a usable Cα mapping cannot enter a pocket.

Selection uses the JSON scores, avoiding the reduced precision of PDB B factors.
The PDB and PyMOL coloring provide a visual representation of those scores.

## PepGLAD pocket format

Each pocket is a JSON list of two-field residue entries:

```json
[
  ["A", [42, " "]],
  ["A", [43, "A"]]
]
```

The second field contains the original residue number and insertion code. A
blank insertion code is represented by a space. Members are sorted by chain,
number, and insertion code. GAPS scores stay in the sidecar and AlloPep's
selection logic; PepGLAD uses its original geometry guidance.

Two residues are the minimum supported by the patched PepGLAD input preparation.
Producing a finite pocket covariance at that size is a preprocessing check,
not evidence that a two-residue pocket yields high-quality designs.

## Historical benchmarks and resume

`GAPS_DEFAULT_SELECTION_REPORT.md` in the repository records historical
optimization evidence. Its radius-limited, score-blended selector and the
earlier connected-pocket threshold experiments differ from the current rules.
Those results do not measure recovery for the current algorithm.

The benchmark workspace is `~/pocket_detection`; consult its own README for
cohorts and saved artifacts. Historical scripts importing the retired
`GuidanceSelection` API require migration before use with this checkout.

Resume requires the same selected site count and exact saved pocket JSON.
Changed thresholds, changed pocket membership, or old three-field scored
entries can cause a mismatch. See [Sampling and resume](sampling-and-resume.md).

Implementation: `gaps_predict.py`, `run_gaps.py::load_scored_residues()`, and
`run_gaps.py::select_contexts()`.
