# AlloPep

AlloPep finds candidate allosteric sites with GAPS, designs peptides with
PepGLAD, checks their structures, and ranks valid designs with Vina.

## Setup

Install Conda, Git, `unzip`, `curl` or `wget`, and `sha256sum` or `shasum`.
Then clone AlloPep and run its installer:

```bash
git clone https://github.com/djmaity/allopep.git
cd allopep
bash install.sh
```

The installer fetches pinned GAPS and PepGLAD code, applies the GAPS patch,
downloads PepGLAD checkpoints, and creates or updates the `gaps`, `PepGLAD`,
and `allopep` environments. It also installs Ray's metrics dependencies and
builds the rotarama cache. Rerun `bash install.sh` to complete an interrupted
installation or update existing environments. The Vina scorer is bundled here.

## Run

From the repository directory, activate the main environment and run the
pipeline on a PDB file:

```bash
conda activate allopep
python allopep.py input/structure.pdb
```

AlloPep starts GAPS and PepGLAD in their own environments. To reuse PepGLAD
structures already generated for this input:

```bash
python allopep.py input/structure.pdb --validate-existing
python allopep.py input/structure.pdb --rank-existing
```

`--validate-existing` checks structures without scoring them.
`--rank-existing` validates and scores them without rerunning GAPS or PepGLAD.
For Vina scoring, `--vina-box-padding` sets the box padding in ångströms and
`--vina-hydrogen-policy {keep,rebuild}` selects how Meeko handles hydrogens.

## Results

For `input/structure.pdb`, results go in `output/structure_output/`:

| File | Contents |
| --- | --- |
| `structure.pdb` | Saved input used when checking existing outputs. A new full run rejects a different PDB with the same filename. |
| `structure_PepGLAD_outputs/` | Generated complex PDBs and PepGLAD summaries. |
| `structure_pepglad_aligned.pml` | PyMOL script showing generated peptides aligned to the input receptor. |
| `structure_vina_score/structure_pepglad_validation.csv` | Status, counts, and rejection reasons for every PepGLAD structure. |
| `structure_vina_score/structure_pepglad_validation.jsonl` | Per-atom validation details. |
| `structure_vina_score/structure_vina_manifest.csv` | Valid structures submitted to Vina. |
| `structure_vina_score/structure_vina_results.csv` | Full Vina scores and errors for valid structures submitted to the scorer. |
| `structure_vina_score/structure_vina_ranked.csv` | Successfully scored structures, ordered from lowest to highest Vina score. |

Open the PyMOL script with
`pymol output/structure_output/structure_pepglad_aligned.pml`. It uses absolute
PDB paths, so regenerate it after moving the output directory.

Validation rejects inverted chirality, severe peptide geometry or continuity
problems, missing peptide atoms, Ramachandran or rotamer outliers, peptide or
interface clashes, and absent receptor contact. Hydrogen atoms are excluded
from heavy-atom geometry checks. If no structures pass, AlloPep stops before
Vina scoring. Inspect the validation CSV and JSONL for the reasons.

Vina scores the generated pose without docking and ranks lower scores first.
These scores and geometry checks do not establish experimental affinity or
dynamic stability. See [the bundled scorer documentation](vina_scoring/README.md)
for scoring details.
