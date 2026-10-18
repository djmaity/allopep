# Installation

AlloPep's installer prepares three separate Conda environments and patched,
pinned GAPS and PepGLAD checkouts. The Vina scorer is bundled in this repository.

## Prerequisites

Use a Linux environment with Conda, Git, `unzip`, either `curl` or `wget`, and
either `sha256sum` or `shasum`. The current sampling lock and validation worker
implementation use Linux/Unix facilities; native Windows execution is not
covered by this workflow.

The standard design pipeline uses CUDA for GAPS, PepGLAD generation, and
PepGLAD's current OpenMM minimizer. Ensure the GPU and driver are available to
the processes running the Conda environments. The top-level CLI has no device
or CPU-only switch.

## Install AlloPep and its dependencies

```bash
git clone https://github.com/djmaity/allopep.git
cd allopep
bash install.sh
```

The installer:

1. Fetches GAPS and PepGLAD at the revisions below.
2. Applies `patches/GAPS.patch` and `patches/PepGLAD.patch`.
3. Downloads and verifies the PepGLAD checkpoint archive when checkpoints are
   missing.
4. Creates or updates the `gaps`, `PepGLAD`, and `allopep` environments.
5. Installs `ray[default]==2.51.2` in `PepGLAD` and rebuilds the CCTBX rotarama
   cache in `allopep`.

| Component | Pinned revision or artifact |
| --- | --- |
| GAPS | `c67cba6ca3ef45470bee8bbb2f0c832680c18044` |
| PepGLAD | `bad015ca50c312a89482adb5220c3d907f13df5c` |
| PepGLAD checkpoint release | `v1.0/checkpoints.zip` |
| Checkpoint archive SHA-256 | `a610e079492a1b1ceab213aea3a0ab875846415e16470897bb54850557f462a0` |

The checksum applies to the downloaded archive. Existing nonempty checkpoint
files are reused by the installer; the sampling controller separately records
the hash of the actual `codesign.ckpt` used for a run.

## Environment responsibilities

| Environment | Responsibility | Selected dependencies |
| --- | --- | --- |
| `allopep` | Pipeline control, accessibility, validation, visualization, Vina scoring | Python 3.12, OpenStructure 2.11.1, CCTBX 2025.9, chem_data 2026.9, FreeSASA, PyMOL 3.1.0, Meeko 0.8.0, Vina 1.2.7 |
| `gaps` | GAPS inference | Python 3.10, PyTorch 1.13.1, PyTorch CUDA 11.7 |
| `PepGLAD` | Peptide generation and OpenMM relaxation | The patched upstream `PepGLAD/env.yaml`, plus Ray's default extras |

Activate `allopep` for the main CLI. AlloPep launches the other two environments
through `conda run`; their names and capitalization matter.

## Bundled Vina scoring

`vina_scoring/` contains the rigid-pose scorer. `environment.yml` installs
Meeko 0.8.0 and Vina 1.2.7 directly. A design or reranking run checks these
dependencies and Meeko's receptor-preparation executable before processing
structures. `--validate-existing` skips scoring dependency checks.

An existing installation can update just the main environment:

```bash
conda env update -n allopep -f environment.yml
```

## Check the installation

From the repository root:

```bash
conda run --no-capture-output -n allopep python allopep.py --help
conda run --no-capture-output -n allopep python -c \
  "import ost, mmtbx, freesasa, pymol; print('Validation and visualization imports OK')"
conda run --no-capture-output -n allopep python -c \
  "from vina_ranking import check_vina_dependencies; check_vina_dependencies()"
conda run --no-capture-output -n gaps python -c \
  "import torch; print('GAPS CUDA available:', torch.cuda.is_available())"
conda run --no-capture-output -n PepGLAD python -c \
  "import torch, openmm; print('PepGLAD CUDA available:', torch.cuda.is_available()); print('OpenMM platforms:', [openmm.Platform.getPlatform(i).getName() for i in range(openmm.Platform.getNumPlatforms())])"
```

These are environment checks, not a complete design or scientific validation
run. Vina scoring uses the bundled source and dependencies in `allopep`.

## Update or repair an installation

Rerun `bash install.sh` to update declared environments or finish an interrupted
installation. The installer accepts an already-applied patch and refuses an
unexpected dependency revision or a patch that cannot be applied cleanly.
Inspect local dependency changes before replacing a checkout.

After environment changes, rebuild the conformation-data cache if needed:

```bash
conda run --no-capture-output -n allopep mmtbx.rebuild_rotarama_cache
```

See [Troubleshooting](troubleshooting.md) for missing data, PyMOL, and dependency
mismatch errors. To work on the documentation alone, use the separate
[documentation environment](development.md#build-and-preview-the-documentation).

Implementation: `install.sh` and `environment.yml`.
