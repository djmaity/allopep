# Rigid-peptide Vina scoring

AlloPep owns this scorer. Its source modules were copied unchanged from
`~/repos/peptide-vina-scoring`, revision
`4a97b5ce7c9df8600d63fc2b135efda8be9cc1db` (package version 0.2.0).
The protocol remains `rigid_meeko_vina_residue_sum_v1`.
The comparison repository imports this directory directly through its
`run_scorers.vina_backend` bridge. No separately installed scorer distribution
is required.

`vina_score.py` evaluates supplied protein-peptide coordinates with the
AutoDock Vina scoring function. It prepares temporary PDBQT files but does not
dock, move, or optimize the peptide. Both receptor and peptide are rigid, so
the peptide has zero active torsions and no torsional flexibility penalty.

Both PDB and mmCIF coordinate inputs are supported. Prefer mmCIF when explicit
connectivity records are available. The selected component is converted to a
temporary PDB with applicable connections represented as `LINK` records before
Meeko's native macromolecule parser is invoked.

Install the dependencies through AlloPep's `environment.yml` and run the
commands below from the AlloPep repository root with `allopep` activated.

## Single complex structure

Supply the peptide chain and optionally enumerate the receptor chains:

```bash
python -m vina_scoring.vina_score complex.pdb \
  --chain B \
  --receptor-chain A \
  --output score.json
```

When `--receptor-chain` is omitted, every non-peptide chain is used as the
receptor.

## Separate receptor and peptide structures

Separate files use all chains in their respective PDB or mmCIF files:

```bash
python -m vina_scoring.vina_score \
  --receptor-pdb receptor.pdb \
  --peptide-pdb peptide.pdb \
  --output score.json
```

Multiple peptide poses can be supplied together. Every pose is scored in an
independent worker and must already use the receptor's coordinate frame:

```bash
python -m vina_scoring.vina_score \
  --receptor-pdb receptor.pdb \
  --peptide-pdb pose_1.pdb pose_2.pdb \
  --output scores.json
```

## Batch manifest

Batch input is a CSV manifest. Paths are resolved relative to the manifest.
Each row uses either `complex_pdb` plus `peptide_chain`, or both `receptor_pdb`
and `peptide_pdb`. These legacy column names accept PDB or mmCIF paths.
Receptor chains in a combined structure are separated with
semicolons.

```csv
id,complex_pdb,receptor_pdb,peptide_pdb,peptide_chain,receptor_chains
complex_1,structures/complex_1.pdb,,,B,A
pose_2,,structures/receptor_2.pdb,structures/peptide_2.pdb,,
pose_3,,structures/receptor_2.pdb,structures/peptide_3.pdb,,
```

Every manifest row is scored in an isolated disposable worker. Distinct jobs
can be processed in parallel:

```bash
python -m vina_scoring.vina_score \
  --manifest jobs.csv \
  --workers 8 \
  --memory-budget-gb 24 \
  --output vina_scores.csv
```

The scheduler estimates memory from each affinity-map box and starts a job
only when its reservation fits the configured budget. Without
`--memory-budget-gb`, 75% of currently available memory is used. This is a
conservative scheduling policy, not an operating-system hard memory limit.
Each complex runs in a disposable worker process so native Vina map memory is
returned to the operating system before that worker is replaced. Vina 1.2.7
creates a quadratic atom-pair table even for a fully rigid peptide. To avoid
that redundant allocation, the workflow builds one complete set of receptor
maps, scores each rigid peptide residue against those maps, and sums the raw
interaction energies before Vina's normal three-decimal rounding. The
requested number of workers can still process different complexes in
parallel.
Failures, including missing input files, are written to their result rows
without stopping other jobs. Output rows retain manifest order.

An optional explicit box can be supplied in a manifest with all six columns:

```text
box_center_x,box_center_y,box_center_z,box_size_x,box_size_y,box_size_z
```

When these columns are empty, the box is calculated from the peptide
coordinates plus `--box-padding` on every side.

## Preparation policy

Meeko's macromolecule residue templates assign the Vina atom types and prepare
both rigid components. The default `--hydrogen-policy keep` retains the common
all-atom input hydrogens; `--hydrogen-policy rebuild` removes them before
template preparation.
`--ph` records the intended protonation pH as provenance, but protonation is
controlled by the selected Meeko residue templates. Additional cached or custom
templates can be supplied with `--meeko-templates`.

Unknown residue names are rejected before Meeko runs, so scoring never depends
on a runtime CCD download. Any additional template JSON should be kept under
version control and supplied explicitly for reproducible reruns.

Every non-water peptide residue is retained, including modified residues and
caps. A missing Meeko template, preparation warning, lost residue or heavy
atom, or changed heavy-atom coordinate makes the complete job fail explicitly.
Receptor hetero residues remain opt-in: repeat `--keep-hetero-residue RESNAME`
for required cofactors or use `--include-all-hetero` after reviewing the input.

## Results

Single-input execution writes JSON. Batch execution writes one CSV row per
manifest row. Every result includes
`protocol_id=rigid_meeko_vina_residue_sum_v1`. Results include:

- all eight values returned by Vina's `score()` call;
- total and interaction scores normalized by peptide heavy atoms and residues;
- peptide atom, heavy-atom, residue, and active-torsion counts;
- box coordinates, preparation settings, checksums, and software versions;
- preparation warnings and per-row errors.

The Vina total is a rigid-pose score in nominal kcal/mol, not an experimental
binding free energy. Compare only structures prepared with the same protocol.

Library code can construct a `ScoreJob` from `CombinedStructure` or
`SeparateStructures`, run it through `ScoringService` or `BatchRunner`, and
serialize its `ScoreSuccess` or `ScoreFailure` with `result_record`. The
`PdbPreparer` facade remains available for independent rigid-PDBQT preparation.

Run the focused tests from the repository root:

```bash
python -m unittest discover -s tests -p 'test_vina_score.py' -v
```
