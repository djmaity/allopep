# AlloPep

AlloPep finds candidate allosteric sites with GAPS, designs peptides with
PepGLAD, checks their structures, and ranks valid designs with Rosetta
InterfaceAnalyzer.

## Setup

Install Conda, Git, `unzip`, `curl` or `wget`, and `sha256sum` or `shasum`.
Run the installer from `~/repos/allopep`:

```bash
cd ~/repos
git clone https://github.com/djmaity/allopep.git
cd allopep
bash install.sh
```

The installer fetches pinned GAPS and PepGLAD code, applies their tracked patches,
downloads PepGLAD checkpoints, and creates or updates the `gaps`, `PepGLAD`,
and `allopep` environments. It also installs Ray's metrics dependencies and
builds the rotarama cache. Rerun `bash install.sh` to complete an interrupted
installation or update existing environments.

Install Rosetta separately and make its `InterfaceAnalyzer` executable
available on `PATH` when running the pipeline. Rosetta is not installed by
`install.sh` or `environment.yml`. You can instead pass
`--rosetta-bin /path/to/InterfaceAnalyzer` to use a specific executable.
See the [Rosetta installation instructions](https://docs.rosettacommons.org/docs/latest/getting_started/Getting-Started).

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
For new designs, GAPS saves full-precision residue scores in a JSON sidecar.
AlloPep chooses local receptor contexts around separated exposed high-scoring
residues and passes each residue's score to PepGLAD for spatial guidance.
The three highest-scoring sites are used by default; `--max_rank all` designs
against every separated site, or pass a positive number. Set
`--score-guidance-strength 0` for a matched run without the score attraction
term. The receptor context is still chosen from GAPS scores.

## Results

For `input/structure.pdb`, results go in `output/structure_output/`:

| File | Contents |
| --- | --- |
| `structure.pdb` | Saved input used when checking existing outputs. A new full run rejects a different PDB with the same filename. |
| `structure_GAPS_scores.json` | Full-precision GAPS scores keyed by the input PDB's chain, residue number, and insertion code. |
| `structure_allosteric_sites/` | Local receptor contexts with a score for every included residue. |
| `structure_PepGLAD_outputs/` | Generated complex PDBs and PepGLAD summaries. |
| `structure_pepglad_aligned.pml` | PyMOL script showing generated peptides aligned to the input receptor. |
| `structure_rosetta_score/structure_pepglad_validation.csv` | Status, counts, and rejection reasons for every PepGLAD structure. |
| `structure_rosetta_score/structure_pepglad_validation.jsonl` | Per-atom validation details. |
| `structure_rosetta_score/structure_rosetta_manifest.csv` | Valid structures submitted to InterfaceAnalyzer. |
| `structure_rosetta_score/structure_rosetta_results.csv` | Interface metrics and per-design errors for every submitted structure. |
| `structure_rosetta_score/structure_rosetta_ranked.csv` | Successfully scored structures, ordered from lowest to highest `rosetta_dG_separated`. |
| `structure_rosetta_score/rosetta_runs/` | Per-design raw scorefiles and logs. |

Open the PyMOL script with
`pymol output/structure_output/structure_pepglad_aligned.pml`. It uses absolute
PDB paths, so regenerate it after moving the output directory.

Validation rejects inverted chirality, severe peptide geometry or continuity
problems, missing peptide atoms, Ramachandran or rotamer outliers, peptide or
interface clashes, and absent receptor contact. Hydrogen atoms are excluded
from heavy-atom geometry checks. If no structures pass, AlloPep stops before
Rosetta scoring. Inspect the validation CSV and JSONL for the reasons.

InterfaceAnalyzer scores each validated, unchanged complex with explicit
receptor and peptide chains. It neither relaxes nor repacks the bound or
separated structures. AlloPep ranks lower `dG_separated` values first using
the bundled 26-term ridge weights. The weights came from
`compare_docking_scoring_function/output/rosetta_affinity/optimized_rosetta_26_terms.wts`;
their SHA-256 is
`e193a7acc0cfe0918a15662cac9d021f232797da81d57aa89e6a1c6005983044`.
The results CSV records the same checksum, executable path, and each design's
raw scorefile and log. Failed Rosetta jobs stay in the results CSV but are
omitted from ranking.
Existing Vina output directories remain historical and are not converted.

These fitted scores are empirical ranking signals, not measured binding free
energies. GAPS-guided sampling and the fitted weights have not established
experimental affinity, allosteric activity, or dynamic stability for the
generated peptides.
