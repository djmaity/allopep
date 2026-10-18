# Scientific interpretation

AlloPep provides candidate sites, generated complexes, physical-validation
records, and empirical rankings. Keep the evidence for those stages distinct
when choosing candidates or reporting results.

## What each result supports

| Result | Supported interpretation | Further evidence needed |
| --- | --- | --- |
| High GAPS residue or pocket score | The model prioritizes this receptor region under the current selector. | Experimental allosteric relevance and independence from model training data. |
| Successful generation or OpenMM relaxation | A candidate was produced and processed. | Passing physical validation and any claim of stability. |
| AlloPep `valid` status | The candidate passed the implemented physical checks. | Native pose accuracy, binding, function, comprehensive stereochemical assessment, and dynamics. |
| Favorable fitted Rosetta interface score | The fixed-pose scorer ranks the candidate favorably among the scored cohort. | Calibrated affinity prediction and held-out performance on generated peptides. |
| Native-pocket recovery | A selected receptor set overlaps a specified native contact set under an explicit metric. | Allosteric activity, peptide design success, and generalization beyond the evaluation cohort. |

## Physical validity and repair

The gate checks specific peptide geometry, conformation, completeness, and
contact requirements. It is stricter than a basic evaluator flag that checks
only peptide length and Cα completeness. Name the validator and policy when
reporting a valid fraction.

Local minimization is a preparation step, not evidence of dynamic stability.
Favorable contacts or scores cannot compensate for failed chirality, peptide
bond, completeness, or clash criteria. Modified structures require full
revalidation. Do not describe relaxed thresholds as repaired coordinates.

## Pocket-recovery evidence

The repository's `GAPS_DEFAULT_SELECTION_REPORT.md` preserves saved historical
selector optimization evidence. Its radius-limited pockets, score cutoff,
relative-accessibility settings, and neighborhood blending differ from the
current connected-pocket selector.

That historical search tested the first pocket only. Its full-cohort best
setting was an **in-sample selection result**; the grouped selection analysis
was **internal evaluation**, not an external validation of one fixed setting.
Neither establishes the current selector's recovery or the recovery of later
ranked pockets. Earlier connected-pocket experiments also used different
thresholds.

Benchmark code and saved artifacts live in `~/pocket_detection`. Record the
exact selector, cohort, native-label definition, rank convention, and metric
threshold when using those results. A 5 Å heavy-atom native contact definition
and an 8 Å Cα definition are different labels; strict `>` and inclusive `>=`
cutoffs also produce different success counts.

## Fair comparisons and reproducible reports

When comparing models or settings:

1. Pair receptor/target IDs and sample indices, using the same number of
   structures per target for each method.
2. State whether results concern all attempts, available complexes, valid
   complexes, selected structures, or successfully scored structures.
3. Report invalidity separately from generation, relaxation, validation, and
   scoring errors. Rejection categories can overlap.
4. Separate full-cohort fitting or tuning from grouped held-out evaluation.
   Record sequence-similarity grouping and leakage checks when relevant.
5. Preserve the input, code/checkpoint/weight hashes, selected pocket identities,
   raw reports, and selection rules. Keep the original cohort reports when
   rechecking or reranking available structures.

The artifacts described in [Outputs and visualization](outputs.md) support
these audits. Experimental affinity, allosteric activity, and dynamic stability
of the generated peptides have not been established by this pipeline alone.
