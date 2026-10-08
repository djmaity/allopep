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

## Output

Results for `input/structure.pdb` are written beneath
`output/structure_output` in the AlloPep repository.
