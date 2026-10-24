# Structure validation

Every usable PepGLAD complex must pass AlloPep's OpenStructure and CCTBX
physical-validity policy before it can be selected for Vina ranking.
**Any one of the 16 rejection reasons below makes a checked structure invalid.**
There is no allowance for a percentage of outliers.

## Status meanings

| Status | Meaning | Eligible for selection and ranking? |
| --- | --- | --- |
| `valid` | Measurements completed and no rejection reason was found. | Yes. |
| `invalid` | Measurements completed and at least one rejection criterion failed. | No. |
| `error` | The complex could not be checked successfully. | No. |

Generation and relaxation failures without a usable PDB are accounted for in
the [sampling summary](sampling-and-resume.md#attempts-are-charged-before-launch),
rather than being forced into a physical-validation row.

## All 16 rejection criteria

The first 15 metrics must be zero. The final criterion requires at least one
peptide residue contacting the receptor within 5 Å.

| Rejection reason | Failure condition |
| --- | --- |
| `peptide_clashes` | At least one detected nonbonded clash entirely within the peptide. |
| `interface_clashes` | At least one detected clash involving both peptide and receptor. |
| `peptide_bad_bonds` | A peptide-involving bond length differs from its reference mean by **>12 standard deviations**. |
| `peptide_bad_angles` | A peptide-involving bond angle differs from its reference mean by **>12 standard deviations**. |
| `inverted_ca` | A non-glycine peptide residue has signed Cα volume **<−0.5 Å³**, indicating inverted chirality in the implemented convention. |
| `flat_ca` | A non-glycine peptide residue has signed Cα volume **−0.5 to +0.5 Å³**, inclusive, or lacks an atom needed to assess that center. |
| `missing_backbone_atoms` | A peptide residue lacks **N, CA, C, or O**. |
| `backbone_breaks` | Consecutive peptide residues have a C–N distance outside **1.1–1.7 Å**, inclusive. Missing atoms are handled by the completeness checks. |
| `sequence_mismatches` | The observed peptide sequence differs from `pep_seq` in PepGLAD's summary. |
| `missing_heavy_atoms` | A peptide residue lacks a required non-hydrogen atom from OpenStructure's compound definition. Terminal **OXT is optional**. |
| `rama_outliers` | CCTBX identifies at least one peptide Ramachandran outlier. |
| `rotamer_outliers` | CCTBX identifies at least one peptide side-chain rotamer outlier. |
| `omega_unassessable` | The ω torsion of a consecutive peptide bond cannot be evaluated. |
| `twisted_peptide_bonds` | A peptide-bond ω torsion deviates **>30°** from the nearest planar cis or trans configuration. |
| `nonpro_cis_peptide_bonds` | A bond has **absolute ω <30°** and its following residue is not proline. |
| `no_receptor_contact_5a` | Zero peptide residues have a heavy atom within **≤5 Å** of a receptor heavy atom. |

Rejection reasons are nonexclusive. A structure with inverted chirality,
Ramachandran outliers, and no receptor contact receives all three reasons;
criterion counts therefore need not sum to the number of invalid structures.

## How the measurements are defined

### Heavy-atom stereochemistry

Hydrogen and deuterium are excluded from OpenStructure's clash, bond, and angle
checks. This avoids comparing OpenMM hydrogen geometry against different
reference tables. Bond and angle violations entirely within the receptor are
not included in the peptide rejection metrics.

OpenStructure's ordinary clash rule is:

```text
nonbonded distance < van der Waals radius 1 + radius 2 − 1.5 Å
```

The installed default radii cover C, N, O, P, and S; bonded atom pairs are
excluded, and potential disulfides have a separate threshold. This is the
specific OpenStructure clash check used here, not an all-atom MolProbity
clashscore. The bond and angle tolerances are **12 standard deviations**,
not 12 Å or 12 degrees.

### Cα chirality

For each non-glycine peptide residue, the code uses CA, N, C, and CB coordinates:

```text
V = ((N − CA) × (C − CA)) · (CB − CA)
```

The L-amino-acid convention in this atom order requires `V > 0.5 Å³` to avoid
the inverted/flat flags. Glycine is excluded. This check concerns peptide Cα
centers; it does not constitute a comprehensive side-chain stereochemistry
audit.

### Ramachandran, rotamers, and peptide bonds

CCTBX `ramalyze` and `rotalyze` supply their outlier classifications. AlloPep
filters their results to the peptide chain and records the number evaluated.
It uses those classifications directly and adds no outlier-percentage
threshold.

For ω normalized to `[-180°, 180°)`, planarity deviation is the smaller distance
to 0° or ±180°. A deviation exactly equal to 30° does not trigger the twisted
flag; cis classification uses a strict absolute ω <30°. Cis bonds immediately
preceding proline are permitted if all other checks pass.

### Contacts and buried area

Distances use peptide and receptor heavy atoms. Reports include atom-pair,
peptide-residue, and receptor-residue contact counts at both 4 Å and 5 Å.
Only the existence of a **5 Å peptide-residue contact** is an acceptance rule.

Buried SASA is reported in Å² as:

```text
(SASA receptor + SASA peptide − SASA complex) / 2
```

There is no positive minimum buried-SASA threshold. A calculated value below
−1 Å² raises an error; small negative values are clamped to zero.

## Reported diagnostics without rejection cutoffs

The following measurements provide context but are not additional policy gates:

- Receptor-only clash count.
- Total cis-bond count, when cis bonds precede proline.
- Four-angstrom contact counts and the number of contacting residues beyond
  the required nonzero 5 Å contact count.
- Minimum interface distance and buried-SASA magnitude.
- Receptor Cα RMSD and maximum displacement after rigid-body alignment.

The receptor comparison still requires matching receptor chains and residue
names, a Cα in each compared residue, and at least three receptor Cα atoms.
A structural inconsistency can therefore produce an error even though RMSD
itself has no acceptance cutoff.

## What produces an error

Examples include an unreadable or malformed PDB, unexpected chains, disagreement
with a peptide-chain override, missing receptor Cα atoms, a receptor sequence
different from the saved input, an unknown peptide compound definition, or a
failed measurement. Errors remain separate from measured invalidity.

The sampling controller can account for individual validation errors and keep
sampling within the budget. The top-level existing-output modes treat a
validation error as an operational failure. A missing validation dependency
can prevent the validation subprocess from running at all.

## Inspect the reports

For an input named `structure.pdb`, cumulative reports are under
`output/structure_output/structure_vina_score/`:

| File | Contents |
| --- | --- |
| `structure_pepglad_validation.csv` | Status, semicolon-separated reasons, scalar metrics, paths, hashes, site, and sample. |
| `structure_pepglad_validation.jsonl` | One JSON object per complex containing `result` and per-atom/per-residue `details`. |
| `structure_pepglad_validation_manifest.json` | Complex paths, expected sequences/chains, site/sample IDs, and reference input paths. |

Each round also retains its own reports. The CSV records candidate and input
SHA-256 values and the OpenStructure version. Sampling runs use these records
and hashes to check resume compatibility.

For a read-only summary of an existing cumulative CSV:

```python
import csv
from collections import Counter
from pathlib import Path

path = Path("output/structure_output/structure_vina_score/"
            "structure_pepglad_validation.csv")
with path.open(newline="") as handle:
    rows = list(csv.DictReader(handle))

print("Status:", Counter(row["status"] for row in rows))
print("Rejections:", Counter(
    reason
    for row in rows if row["status"] == "invalid"
    for reason in row["reasons"].split(";") if reason
))
```

See [report retention](outputs.md#report-and-structure-retention) before
rechecking existing outputs: rejected coordinates may already have been deleted.

## Validation, repair, and scientific limits

Relaxing a threshold changes which structures pass; it does not repair their
coordinates. Idealization and minimization may change local geometry, but a
modified complex must undergo the complete validation policy again. Missing
atoms, sequence differences, chirality inversions, chain breaks, and absent
interface contact can require rebuilding or repositioning rather than another
local minimization.

Passing these checks is an operational physical-validity result. It does not
establish native pose recovery, affinity, allosteric function, or stability.
See [Scientific interpretation](scientific-scope.md).

Implementation: `allopep/validation.py::validation_reasons()` defines the gate;
the measurement helpers in that module define each metric. Threshold meanings
also depend on the installed `ost.mol.alg.stereochemistry` implementation.
