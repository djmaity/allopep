# Vina ranking

AlloPep ranks selected physically valid complexes with its bundled rigid-pose
Meeko/Vina scorer. **Lower `score_total` ranks first.** Validation establishes
eligibility before scoring; site and per-site sample IDs break score ties.

## Scoring protocol

The scorer selects explicit receptor and peptide chains from PepGLAD's
metadata. Coordinate chain IDs must match the summary. Meeko prepares both
components as rigid PDBQT structures, preserving peptide residues, heavy
atoms, and their coordinates. Missing templates, preparation warnings, or
changed/lost heavy atoms fail the complete pose.

The calculation scores the saved coordinates without docking or optimization.
It computes one complete set of receptor affinity maps, scores each peptide
residue against those maps, and sums raw interaction energies before Vina's
three-decimal rounding. The peptide has zero active torsions. This avoids
Vina's redundant quadratic allocation for a complete rigid peptide.

| Setting | Value |
| --- | --- |
| Protocol ID | `rigid_meeko_vina_residue_sum_v1` |
| Meeko / Vina | `0.8.0` / `1.2.7` in the main environment |
| Hydrogen policy | Keep input hydrogens |
| Automatic box | Peptide bounds plus 5 Å on each side |
| Grid spacing | 0.375 Å |
| Intended pH provenance | 7.4; Meeko templates control protonation |
| Receptor hetero residues | Excluded by default |
| Peptide residues | All non-water residues, including supported caps/modifications |
| Pipeline workers | One disposable scoring process at a time |

The scheduler reserves memory using the affinity-map box, with a default
budget of 75% of available memory. This reservation is a scheduling estimate,
not an operating-system memory limit.

## Results and failures

Reports are written under `output/structure_output/structure_vina_score/`:

| File | Meaning |
| --- | --- |
| `structure_vina_manifest.csv` | Selected valid designs, chain IDs, input hashes, and preparation settings. |
| `structure_vina_results.csv` | Every submitted design, including failures, all score components, normalized diagnostics, hashes, software versions, and errors. |
| `structure_vina_ranked.csv` | Successful finite totals, ordered ascending, with rank starting at one. |

Each scoring job runs in a disposable worker. Temporary prepared files are
removed after scoring; diagnostics and provenance remain in the CSV reports.
Failed jobs are excluded from ranking. All jobs failing causes an operational
failure. Scoring errors do not trigger replacement generation. Reranking
clears the previous Vina scoring reports before writing the current cohort.

## Rerank available complexes

```bash
conda activate allopep
python allopep.py input/structure.pdb --rank-existing
```

This validates available saved candidates and scores those passing the gate.
It can include valid extras retained from the final sampling round. It
refreshes cumulative validation and scoring reports without generating
peptides. Archive reports first when preserving a particular analysis.
Earlier Rosetta reports remain historical artifacts; compatible saved sampling
rounds can be resumed with the new scorer.

## Source and interpretation

`vina_scoring/` is the canonical source, imported directly by the comparison
repository through its `run_scorers.vina_backend` bridge. It was copied from
`peptide-vina-scoring` revision
`4a97b5ce7c9df8600d63fc2b135efda8be9cc1db`, retaining the scoring protocol and
result schema. No separate scorer installation is required. The standalone
CLI is available from AlloPep's root:

```bash
python -m vina_scoring.vina_score --help
```

Vina totals have nominal kcal/mol units. They are empirical ranking signals;
they do not establish experimental affinity, allosteric activity, or dynamic
stability. Compare scores prepared with the same protocol. See
[Scientific interpretation](scientific-scope.md).
