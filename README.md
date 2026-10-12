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
builds the rotarama cache. The `allopep` environment includes PyMOL to build
the saved visualization session. Rerun `bash install.sh` to complete an
interrupted installation or update existing environments.

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

The same commands accept `-v` and `-R`, respectively. Run each script with
`-h` to see its command-line options.

`--validate-existing` checks structures without scoring them.
`--rank-existing` validates and scores them without rerunning GAPS or PepGLAD.
If a design run stops after some PepGLAD sites finish, rerun it with `--resume`
to reuse sites whose saved pocket, sample IDs, and relaxed PDBs match. For a
site interrupted during OpenMM relaxation, inspect its codesign directory with
`conda run -n PepGLAD python recover_pepglad_relaxation.py SITE_DIR`; add
`--apply` to quarantine malformed PDBs and relax unfinished designs, then
rerun AlloPep with `--resume`.
Validation checks independent complexes in parallel, automatically choosing a
worker count from the available CPUs, memory, and number of complexes.
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
| `structure_GAPS_output.pdb` | GAPS receptor structure with residue scores in its B factors. |
| `structure_allosteric_sites/` | Local receptor contexts with a score for every included residue. |
| `structure_PepGLAD_outputs/` | Generated complex PDBs and PepGLAD summaries. |
| `structure_pepglad_aligned.pse` | PyMOL session showing the GAPS receptor colored by score, selected site surfaces, and aligned PepGLAD peptides. |
| `structure_rosetta_score/structure_pepglad_validation.csv` | Status, counts, and rejection reasons for every PepGLAD structure. |
| `structure_rosetta_score/structure_pepglad_validation.jsonl` | Per-atom validation details. |
| `structure_rosetta_score/structure_rosetta_manifest.csv` | Valid structures submitted to InterfaceAnalyzer. |
| `structure_rosetta_score/structure_rosetta_results.csv` | Interface metrics and per-design errors for every submitted structure. |
| `structure_rosetta_score/structure_rosetta_ranked.csv` | Successfully scored structures, ordered from lowest to highest `rosetta_dG_separated`. |
| `structure_rosetta_score/rosetta_runs/` | Per-design raw scorefiles and logs. |

Open the saved session with
`pymol output/structure_output/structure_pepglad_aligned.pse`.
The `sites` group contains pocket surfaces named `site_1`, `site_2`, and so on.
The `peptides` group contains separate `site_1_peptides`, `site_2_peptides`,
and other per-site peptide groups. The GAPS
receptor is named `<structure>`; peptide objects use names such as
`site_1_peptide_1`. Only peptide chains from PepGLAD complexes are shown.
Peptides and pocket surfaces share each site's color, following PyMOL's
default object color sequence starting with green for site 1. Peptide stick
models use the site color for carbon and element colors for all other atoms.
The receptor cartoon uses blue for low B-factor scores
and red for high scores (0–1).

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
