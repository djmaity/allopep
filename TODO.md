# AlloPep improvement checklist

## Goals, in priority order

- [ ] **Benchmark known allosteric site recovery.** Curate the 44-record AlloBench peptide set. Run GAPS on peptide-stripped bound receptors and on matched apo receptors where available. Report top-1 and top-3 known-site recovery per target.
- [ ] **Verify benchmark interfaces.** Check receptor and peptide chain mappings, derive receptor contact residues from the complexes, and document exclusions. Treat this as a positive-only dataset; do not claim allosteric specificity from it.
- [ ] **Diagnose structural validation failures.** For recovered sites, count generated and valid PepGLAD designs and summarize each rejection reason. Compare validation results with native complexes before changing thresholds.
- [ ] **Evaluate Vina ranking.** Use native poses and controlled decoys from the 1,048-complex structural set. Report ranking performance by target and peptide length.
- [ ] **Connect candidates to function.** Keep predicted binding separate from demonstrated changes in target activity.

## First deliverable

- [ ] Produce a reproducible known-site recovery report with dataset coverage, exclusions, per-target top-1 and top-3 results, apo comparisons where available, and an exposed-surface baseline.
