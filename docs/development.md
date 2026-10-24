# Development

The repository coordinates patched external models and local validation,
sampling, visualization, and scoring code. Use the environment appropriate to
each component, and keep documentation aligned with the actual CLI and saved
artifact formats.

## Source map

The pipeline is the `allopep` package, run as `python -m allopep`.

| File or directory | Responsibility |
| --- | --- |
| `allopep/__main__.py` | Entry point for `python -m allopep`. |
| `allopep/cli.py` | Command line, stage orchestration, and validation invocation. |
| `allopep/layout.py` | Repository anchors and every output path, derived in one place. |
| `allopep/fileio.py` | Atomic file replacement and content hashing. |
| `allopep/gaps_predict.py` | GAPS inference and full-precision residue score export. |
| `allopep/gaps.py` | Residue mapping, accessibility, connected pockets, and ranking. |
| `allopep/sampling.py` | Shared or per-site budgets, adaptive rounds, selection, and integrity checks. |
| `allopep/peptide_exports.py` | Selected-complex copies and the peptide/pocket/scoring CSV. |
| `allopep/process_control.py` | Logged subprocess execution and deadline/process-tree handling. |
| `allopep/outputs.py` | Summary parsing, chain/sequence metadata, round IDs, and deletion journals. |
| `allopep/validation.py` | OpenStructure/CCTBX measurements, policy, parallel workers, and reports. |
| `allopep/validation_reports.py` | Validation manifest, `ost` invocation, and report-to-design matching. |
| `allopep/pymol_session.py` | PyMOL scripts and saved sessions. |
| `allopep/ranking.py` | Valid-design submission, diagnostics, and Vina ranked CSVs. |
| `allopep/recover_relaxation.py` | Inspection and optional recovery of interrupted relaxation. |
| `allopep/scoring/` | Rigid-pose Meeko/Vina scorer, also usable on its own as `python -m allopep.scoring`. |
| `install.sh`, `environment.yml` | Pinned dependency setup and the main environment declaration. |
| `patches/` | Changes applied to the pinned external checkouts. |
| `tests/` | Focused regression tests and fixtures. |
| `docs/`, `mkdocs.yml` | This documentation site. |
| `normal_mode_analysis.ipynb` | Reference material kept for planned normal-mode support; not part of the pipeline. |

Two constraints are easy to break by accident:

- **`allopep/__init__.py` imports nothing.** It is imported from the `gaps`
  environment for `allopep.gaps_predict` and from the `PepGLAD`
  environment for `allopep.recover_relaxation`, and neither has the
  `allopep` environment's dependencies.
- **`allopep/validation.py` imports nothing from its own package.** The
  `ost` interpreter runs it by absolute path rather than importing it.

## Run the focused tests

From the repository root:

```bash
conda run --no-capture-output -n allopep \
  python -m unittest discover -s tests -v
```

PepGLAD-specific cases use its environment and their explicit test filename:

```bash
conda run --no-capture-output -n PepGLAD \
  python tests/pepglad_guidance_case.py
```

For installer syntax:

```bash
bash -n install.sh
```

Because the package is imported from three environments, check that the two
foreign ones can still load it after changing imports:

```bash
conda run -n gaps python -c "import allopep, allopep.gaps_predict"
conda run -n PepGLAD python -c "import allopep, allopep.recover_relaxation"
```

Focused tests verify their particular code paths. They do not constitute a
complete GPU design run or establish physical-validity improvements, pocket
recovery, binding, or allosteric function.

## Maintain dependency patches

The installer checks pinned upstream revisions before applying patches. When
changing a patched dependency, keep the installed checkout and tracked patch
synchronized. A reverse check can confirm that the patch is represented in an
already-patched checkout:

```bash
git -C PepGLAD apply --reverse --check --ignore-space-change ../patches/PepGLAD.patch
git -C GAPS apply --reverse --check --ignore-space-change ../patches/GAPS.patch
```

These commands check rather than modify the checkout. Also verify forward
application against a clean checkout at the pinned revision when updating a
patch. Preserve unrelated local changes and expensive experiment artifacts.

## Build and preview the documentation

Build and preview the documentation in the `allopep` Conda environment.
MkDocs and its Material theme are declared in `environment.yml`. Update
an existing environment before building:

```bash
conda env update -n allopep -f environment.yml
conda run --no-capture-output -n allopep \
  python -m mkdocs build --strict
```

The static HTML site is written to `site/`, with direct page links such as
`validation.html`. After building, you can open `site/index.html` to browse
the pages directly without directory listings. The `site/` directory is
ignored by Git. For the full preview, including search, use the local server:

```bash
conda run --no-capture-output -n allopep python -m mkdocs serve
```

Open `http://127.0.0.1:8000/`. Material for MkDocs provides navigation and
search with its default settings. Publishing is a separate step; the
configuration does not presume a hosted site URL.

The environment declaration constrains compatible MkDocs and Material
major versions; it is not a complete dependency lock. Record
`conda env export -n allopep` when toolchain reproduction is needed.

## Documentation maintenance

Use relative Markdown links between pages under `docs/`. MkDocs discovers
pages automatically and uses their first heading as the navigation title;
`mkdocs.yml` has no explicit navigation entries. Add new Markdown pages
directly under `docs/` and update the guide links on the [home page](index.md).
Use standard Markdown formatting and run `python -m mkdocs build --strict`
before committing documentation changes.

Update related documentation together when behavior changes:

| Code change | Pages to review |
| --- | --- |
| CLI names/defaults or length handling | [Usage](usage.md), [Home](index.md), README |
| Pocket thresholds, grouping, or limits | [Pocket selection](pocket-selection.md), [Scientific interpretation](scientific-scope.md) |
| Budgets, round layout, resume, or cleanup | [Sampling](sampling-and-resume.md), [Outputs](outputs.md), [Troubleshooting](troubleshooting.md) |
| Validation thresholds or rejection policy | [Validation](validation.md), report interpretation, relevant tests |
| Scoring preparation or protocol | [Vina ranking](scoring.md), source provenance, relevant tests |
| Environment or external dependency pins | [Installation](installation.md), installer, environment declaration |

Use `git diff --check --` with the changed paths for whitespace validation.
Generated patches can contain unrelated pre-existing whitespace; assess those
separately from documentation edits.

## Upstream references

- [GAPS source](https://github.com/hongliangduan/GAPS)
- [PepGLAD source](https://github.com/THUNLP-MT/PepGLAD)
- [OpenStructure documentation](https://openstructure.org/docs/)
- [CCTBX source](https://github.com/cctbx/cctbx_project)
- [AutoDock Vina documentation](https://autodock-vina.readthedocs.io/)
- [MkDocs configuration](https://www.mkdocs.org/user-guide/configuration/)
- [Material for MkDocs](https://squidfunk.github.io/mkdocs-material/)
