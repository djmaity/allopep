# AlloPep

AlloPep is a pipeline for de novo design of allosteric peptides.

## Installation

AlloPep requires Conda and the following command-line tools:

- Git
- `unzip`
- either `curl` or `wget`
- either `sha256sum` or `shasum`

Clone this repository and enter its directory:

```bash
git clone https://github.com/djmaity/allopep.git
cd allopep
```

Run the dependency installer to fetch the pinned GAPS and PepGLAD revisions:

```bash
bash install.sh
```

The installer also applies AlloPep's GAPS compatibility patch and downloads
and verifies the PepGLAD v1.0 model checkpoints.

The script may be run again after an interrupted or completed installation.
Already-installed dependencies and checkpoints are reused. If an existing
dependency checkout is at an incompatible revision, the script stops with
instructions instead of overwriting it.

Create the three Conda environments. GAPS uses Python 3.10, PyTorch 1.13.1,
and CUDA 11.7 as specified in its README. NumPy 1.26.4 and MKL 2024.0.0
keep that PyTorch build importable with current Conda packages:

```bash
conda create -n gaps python=3.10
conda install -n gaps pytorch=1.13.1 torchvision=0.14.1 \
  torchaudio=0.13.1 pytorch-cuda=11.7 -c pytorch -c nvidia
conda install -n gaps pandas scikit-learn tqdm h5py gemmi \
  numpy=1.26.4 mkl=2024.0.0 -c conda-forge
conda env create -f PepGLAD/env.yaml
conda env create -f environment.yml -n allopep
conda run -n allopep mmtbx.rebuild_rotarama_cache
```

The PepGLAD environment follows its upstream `env.yaml` and is named
`PepGLAD`. The `allopep` environment contains the remaining pipeline tools,
including Rosetta and structural validation. OpenStructure provides geometry,
contacts, and solvent accessibility; CCTBX and `chem_data` provide reference
data for Ramachandran and rotamer checks. The cache command above prepares
those reference tables. Rosetta 2025.12 is pinned because the earlier 2025.06
Conda build cannot be solved with OpenStructure and CCTBX in one environment.

Because the script resolves paths relative to its own location, it can also be
launched from another directory:

```bash
bash /path/to/allopep/install.sh
```

## Running

Activate the main environment and run AlloPep from the repository directory.
The script starts GAPS and PepGLAD with their dedicated Conda environments:

```bash
conda activate allopep
python allopep.py input/structure.pdb
```

Replace `input/structure.pdb` with the path to the input PDB file.
To check already generated PepGLAD structures without rerunning design or
Rosetta scoring, run:

```bash
python allopep.py input/structure.pdb --validate-existing
```

## Output

Results for `input/structure.pdb` are written beneath
`output/structure_output` in the AlloPep repository.
The original PDB is copied to `output/structure_output/structure.pdb` before
the pipeline runs. `--validate-existing` uses this saved copy when present, so
its receptor comparison uses the input associated with those outputs. Runs
with the same input filename use the same output directory. The pipeline stops
if that directory already contains an input PDB with different contents.

The workflow checks every PepGLAD complex before Rosetta scoring. Invalidating
checks cover peptide heavy-atom clashes, severe bond or angle violations,
inverted or flat C-alpha chirality, missing or broken backbones, sequence
mismatches against PepGLAD's `summary.jsonl`, missing standard heavy atoms,
Ramachandran or rotamer outliers, peptide bonds twisted more than 30° from
planar cis/trans geometry, non-proline cis peptide bonds, and no receptor
heavy-atom contact within 5 Å. Cis-proline bonds are reported but not rejected.
OpenStructure uses a 1.5 Å clash tolerance and flags bonds or angles beyond 12
reference standard deviations. Hydrogen atoms and optional terminal OXT are
excluded from completeness checks; hydrogens are also omitted from geometry
checks because OpenMM-relaxed hydrogen bond lengths can differ from the
OpenStructure reference values.

The report also contains 4 and 5 Å heavy-atom contact-pair counts, the number
of peptide and receptor residues contacted at each cutoff, and buried solvent-
accessible surface area calculated as `(SASA_receptor + SASA_peptide -
SASA_complex) / 2`. Receptor displacement from the input PDB is reported as
RMSD and maximum C-alpha displacement after aligning the receptor C-alpha
atoms; the report records the input PDB hash, and no universal displacement
cutoff is applied. Receptor-only clashes are
reported separately. The per-sample counts and reasons are saved as
`*_pepglad_validation.csv`, with atom-level details in
`*_pepglad_validation.jsonl` (including each measured omega angle), under
`output/structure_output/structure_rosetta_score/`.
Validation errors stop the workflow; if no designs pass, it stops before
Rosetta scoring.
These geometric checks do not establish binding affinity or dynamic stability.
