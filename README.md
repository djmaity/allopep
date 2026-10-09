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

Create and activate the AlloPep Conda environment:

```bash
conda env create -f environment.yml -n allopep
conda activate allopep
```

Run the dependency installer:

```bash
bash install.sh
```

The installer:

- installs the compatible GAPS and PepGLAD revisions under the
  repository directory;
- applies AlloPep's GAPS compatibility patch; and
- downloads and verifies the PepGLAD v1.0 model checkpoints.

The script may be run again after an interrupted or completed installation.
Already-installed dependencies and checkpoints are reused. If an existing
dependency checkout is at an incompatible revision, the script stops with
instructions instead of overwriting it.

OpenStructure validation runs in a separate Conda environment so installing it
does not replace the CUDA packages used by PepGLAD:

```bash
conda create -n allopep-ost -c conda-forge -c bioconda openstructure=2.11.1
export ALLOPEP_OST="$(conda info --base)/envs/allopep-ost/bin/ost"
```

Set `ALLOPEP_OST` in each shell used to run AlloPep, or put the OpenStructure
`ost` executable on `PATH`.

Because the script resolves paths relative to its own location, it can also be
launched from another directory:

```bash
bash /path/to/allopep/install.sh
```

## Running

Activate the environment and run AlloPep from the repository directory:

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

The workflow checks every PepGLAD complex with OpenStructure before Rosetta
scoring. It excludes structures with peptide heavy-atom clashes, severe bond or
angle violations, inverted or flat C-alpha chirality, missing or broken peptide
backbones, or no receptor heavy-atom contact within 5 Å. The check also reports
receptor-only clashes separately. OpenStructure uses a 1.5 Å clash tolerance
and flags bonds or angles beyond 12 reference standard deviations. Hydrogen
atoms are omitted from geometry checks because the OpenMM-relaxed hydrogen
bond lengths can differ from the OpenStructure reference values. The per-sample
counts and reasons are saved as
`*_pepglad_validation.csv`, with atom-level details in
`*_pepglad_validation.jsonl`, under `output/structure_output/structure_rosetta_score/`.
If no designs pass, the workflow stops before Rosetta scoring.
These geometric checks do not establish binding affinity or dynamic stability.
