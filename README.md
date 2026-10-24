# AlloPep

AlloPep finds candidate allosteric sites with GAPS, designs peptides with
PepGLAD, checks their structures, and ranks valid designs with rigid-pose
Vina scoring.

## Documentation

The [documentation](docs/index.md) covers installation, command-line usage,
pocket selection, resumable sampling, all 16 [structure-validation
criteria](docs/validation.md), Vina ranking, outputs, and troubleshooting.
To preview the searchable MkDocs site locally:

```bash
python -m venv .venv-docs
.venv-docs/bin/python -m pip install -r requirements-docs.txt
.venv-docs/bin/python -m mkdocs serve
```

Open `http://127.0.0.1:8000/`. See the [development guide](docs/development.md)
for strict builds and documentation maintenance.

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

The Vina scorer is included in `allopep/scoring/`. The `allopep` environment
installs Meeko 0.8.0 and Vina 1.2.7 directly; no separate scorer package or
checkout is needed to run AlloPep.

## Run

From the repository directory, activate the main environment and run the
pipeline on a PDB file:

```bash
conda activate allopep
python -m allopep input/structure.pdb
```

`-n / --num_samples` requests **physically valid structures per selected site**
(default: 10). AlloPep keeps sampling until each site reaches that target or
its own `--max_num_peptides` attempt cap. The cap defaults to **1,000 times
`-n` per site**; for example, `-n 10` permits up to 10,000 candidates at each
site. Change the multiplier with `--attempt-multiplier`, a positive integer:

```bash
python -m allopep input/structure.pdb -n 10 --attempt-multiplier 200
```

This permits up to **200 × 10 = 2,000 attempts per site**. An explicit
`--max_num_peptides` takes precedence over the multiplier and sets an absolute
per-site total. It may be higher or lower than the default, but must be at
least `-n`. For example, set a smaller budget:

```bash
python -m allopep input/structure.pdb -n 10 --max_num_peptides 1000
```

Every candidate consumes an attempt, including failed coordinates, failed
relaxation, and rejected structures. Different poses with the same sequence
count separately. Validity uses the existing OpenStructure/CCTBX gate;
Vina scoring failures do not trigger replacement sampling.

Unfinished sites receive one round each in GAPS-rank order. Each round runs a
separate PepGLAD process with GPU minibatches of eight. Validation rounds start
at eight candidates, double when no structures pass, and otherwise estimate
the candidates needed from that site's cumulative acceptance rate. They grow
by at most a factor of two, stay between eight and 64, and shrink to the
remaining budget when necessary. Once a site reaches its target, no further
rounds are scheduled there. The first `n` valid structures in sample order are
selected; any additional valid structures in the final round remain on disk
and in the validation reports.

After each completed round at each site, AlloPep deletes rejected PDBs,
including validation failures, validation errors, and quarantined generation
or relaxation failures. It keeps every valid PDB, the rejection reports,
PepGLAD summaries, error files, and logs. A `removed_structures.json` journal
records intentional deletions before cleanup, allowing interrupted cleanup
to finish on resume without resetting the attempt budget. Unfinished rounds
are kept for recovery. Resuming older validated rounds also performs this
cleanup using their saved rejection results.

`--round-timeout-seconds` (default: 3,600) covers one round's generation,
relaxation, and validation. A timeout terminates the subprocess and its
descendants, saves progress and logs, and reports an operational failure.
Individual OpenMM sample failures are reported and temporarily quarantined
until round cleanup, without aborting other samples. A site's sampling
shortfall does not stop the other sites: AlloPep retains and ranks the partial
valid results. Exit statuses are
`0` when all sites reach their targets, `2` for a sampling shortfall, and `1`
for an operational failure. Existing Vina scoring-error reports remain
available; all scoring jobs failing is an operational failure.

AlloPep starts GAPS and PepGLAD in their own environments.
PepGLAD idealizes each generated peptide's backbone and side chains before
writing and relaxing the complex. This is the default for new AlloPep designs.
To reuse PepGLAD structures already generated for this input:

```bash
python -m allopep input/structure.pdb --validate-existing
python -m allopep input/structure.pdb --rank-existing
```

The same commands accept `-v` and `-R`, respectively. Run each script with
`-h` to see its command-line options.

`--validate-existing` checks structures without scoring them.
`--rank-existing` validates and scores them without rerunning GAPS or PepGLAD.
Fresh generation refuses existing candidate outputs. Rerun with `--resume`
to reuse checked rounds and continue unfinished sites. Resume verifies the
input, selected pockets, peptide-length settings, peptide-chain override,
generator/checkpoint, validation policy, and saved artifact hashes. Already
validated rounds are not validated again. Legacy flat outputs are checked and
charged against the budget before new rounds are generated. You may adjust
`-n`, the attempt multiplier, and the absolute attempt cap on resume; the
resulting cap cannot be below attempts already charged, and the cumulative
budget never resets.

Each site's `round_000001/`, `round_000002/`, and later directories keep their
own PDBs, summary, logs, and validation reports. Round-local filenames repeat,
but reserved sample ranges give every candidate a unique per-site sample ID
in validation and ranking reports. The complete round budget is reserved
before launch. An interrupted round remains charged even if some candidates
were never written; resume can validate its saved relaxed candidates without
resampling that reserved range.

For a round interrupted during OpenMM relaxation, inspect its directory with
`conda run -n PepGLAD python -m allopep.recover_relaxation ROUND_DIR`; add
`--apply` to quarantine malformed PDBs and relax unfinished designs, then
rerun AlloPep with `--resume`. Use the site directory for legacy flat outputs.
Unrelaxed saved structures require this recovery before resume can continue.
After a selector default change, `--resume` checks that the selected site count
matches any saved sampling checkpoint and compares each site's residue
identities with its saved pocket JSON before reusing that site's PepGLAD
outputs. A changed site count or pocket raises a mismatch error; matching
pockets can still be reused. Existing scored three-field pockets also fail
this exact match against the current two-field format.
Validation checks independent complexes in parallel, automatically choosing a
worker count from the available CPUs, memory, and number of complexes.
For new designs, GAPS saves full-precision residue scores in a JSON sidecar.
AlloPep keeps solvent-accessible residues with **GAPS score >0.5** and
**absolute solvent-accessible surface area (SASA) >0 Å²**. Any positive SASA
qualifies; residues with SASA zero or missing accessibility are excluded.
It groups these residues into connected pockets
using **Cα distances ≤8 Å**, including connections across chains. Connections
are transitive: A and C share a pocket when both connect to B, even if A and C
are more than 8 Å apart. Filtered-out residues cannot act as bridges.
Each pocket contains at least **two residues**, matching the minimum supported
by the current PepGLAD input preparation, and has no upper size limit or
center-radius restriction. Every residue belongs to at most one pocket;
every residue pair from separate pockets is more than 8 Å apart.
Pockets are ranked by **mean GAPS score**, highest first. Ties use the
smallest residue identity (chain, residue number, insertion code); pocket
members are sorted by the same identity. All pockets are used by default;
`--max_rank N` limits design to the first N pockets. AlloPep passes only
residue identities to PepGLAD, which uses its original geometry guidance
without GAPS scores.
The `GAPS_DEFAULT_SELECTION_REPORT.md` report records the saved optimization
evidence and its limits. That search used a ≥0.4 score cutoff for centers
and pocket members and measured only the first pocket.
It does not establish recovery for the current connected pockets or for later
ranked pockets. The old selector tuning API has been retired; historical
`pocket_detection` scripts that import it require a separate migration.

### Check pocket recovery on bound protein-peptide structures

The benchmark code, instructions, and saved full and pilot results are in
`~/pocket_detection`. It compares the first GAPS pocket with receptor residues
contacting each bound peptide chain; see `~/pocket_detection/README.md` for the
cohort, evaluation method, and commands.

## Results

For `input/structure.pdb`, results go in `output/structure_output/`:

| File | Contents |
| --- | --- |
| `structure.pdb` | Saved input used when checking existing outputs. A new full run rejects a different PDB with the same filename. |
| `structure_GAPS_scores.json` | Full-precision GAPS scores keyed by the input PDB's chain, residue number, and insertion code. |
| `structure_GAPS_output.pdb` | GAPS receptor structure with residue scores in its B factors. |
| `structure_allosteric_sites/` | Local receptor pocket definitions containing chain, residue number, and insertion code. |
| `structure_allosteric_site.pml` / `structure_allosteric_site.pse` | PyMOL script and session showing the original input receptor colored by GAPS score (0–1) and every selected allosteric site as a surface. Created before peptide sampling. |
| `structure_PepGLAD_outputs/` | Generated complex PDBs and PepGLAD summaries. |
| `structure_sampling_summary.json` | Per-site target, charged attempts, valid count, selected paths, failure counts, reserved rounds, and completion reason; also the resume checkpoint. |
| `structure_PepGLAD_outputs/*_codesign/round_*/` | Isolated round outputs, `round.json` sample range, `pepglad.log`, `validation.log`, and per-round validation reports. |
| `structure_PepGLAD_outputs/*_codesign/round_*/removed_structures.json` | Rejected PDB paths, statuses, and original hashes; permits safe resume after deleting those structures. |
| `structure_pepglad_aligned.pse` | PyMOL session showing the original input receptor colored by GAPS score, selected site surfaces, and selected valid PepGLAD peptides. |
| `structure_vina_score/structure_pepglad_validation.csv` | Status, counts, and rejection reasons for every PepGLAD structure. |
| `structure_vina_score/structure_pepglad_validation.jsonl` | Per-atom validation details. |
| `structure_vina_score/structure_vina_manifest.csv` | Selected valid structures submitted to Vina. |
| `structure_vina_score/structure_vina_results.csv` | Vina components, provenance, and per-design errors for every submitted structure. |
| `structure_vina_score/structure_vina_ranked.csv` | Successfully scored structures, ordered from lowest to highest `score_total`. |

Open the allosteric site view with
`pymol output/structure_output/structure_allosteric_site.pse`, or run the
matching `.pml` script in PyMOL. Each pocket surface can be toggled in the
`sites` group. The view loads `structure.pdb` and applies full-precision scores
from `structure_GAPS_scores.json` by chain, residue number, and insertion code.
These files are also refreshed by `--validate-existing` and
`--rank-existing` using the saved site definitions.

Open the session with selected valid peptides using
`pymol output/structure_output/structure_pepglad_aligned.pse`.
The `sites` group contains pocket surfaces named `site_1`, `site_2`, and so on.
The `peptides` group contains separate `site_1_peptides`, `site_2_peptides`,
and other per-site peptide groups. The original input
receptor is named `<structure>`; peptide objects use names such as
`site_1_peptide_1`. Only peptide chains from PepGLAD complexes are shown.
Peptides and pocket surfaces share each site's color, following PyMOL's
default object color sequence starting with green for site 1. Peptide stick
models use the site color for carbon and element colors for all other atoms.
The receptor cartoon uses blue for low GAPS scores and red for high scores
(0–1). Scores are assigned to B factors inside the visualization; residues
without predictions remain gray.

For new runs, cumulative validation CSV/JSONL reports retain the results for
every checked structure, including deleted rejected structures and final-round
extras. Report paths for deleted structures remain as historical identifiers;
their PDB files are removed. Generation and relaxation failures without usable
PDBs remain in the round summaries,
error files, and sampling-summary counts. `--validate-existing` and
`--rank-existing` discover both legacy flat outputs and round directories;
they check all available structures and perform no generation.

Validation rejects inverted chirality, severe peptide geometry or continuity
problems, missing peptide atoms, Ramachandran or rotamer outliers, peptide or
interface clashes, and absent receptor contact. Hydrogen atoms are excluded
from heavy-atom geometry checks. If no structures pass, AlloPep stops before
Vina scoring. Inspect the validation CSV and JSONL for the reasons.

The bundled scorer prepares each selected complex with Meeko, using explicit
receptor and peptide chains from PepGLAD. It scores the saved rigid pose
without docking or coordinate optimization. One set of receptor maps is used
to score peptide residue fragments; their raw interaction energies are summed
before rounding. Lower `score_total` ranks first, with site and sample IDs
breaking ties.

The protocol is `rigid_meeko_vina_residue_sum_v1`. Defaults retain input
hydrogens, pad the peptide bounds by 5 Å per side, and use 0.375 Å grid
spacing. Results record score components, normalized diagnostics, input
hashes, preparation settings, and Vina/Meeko versions. Failed jobs remain in
the results CSV and are excluded from ranking. See [Vina ranking](docs/scoring.md).

Earlier Rosetta reports remain historical artifacts. Saved generation rounds
can be resumed or reranked with Vina; their original reports are preserved.
Vina scores are empirical ranking signals in nominal kcal/mol. Experimental
affinity, allosteric activity, and dynamic stability of generated peptides
have not been established by this pipeline.
