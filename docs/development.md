# Development

The repository coordinates patched external models and local validation,
sampling, visualization, and scoring code. Use the environment appropriate to
each component, and keep documentation aligned with the actual CLI and saved
artifact formats.

## Source map

| File or directory | Responsibility |
| --- | --- |
| `allopep.py` | CLI, stage orchestration, validation invocation, and PyMOL artifacts. |
| `gaps_predict.py` | GAPS inference and full-precision residue score export. |
| `run_gaps.py` | Residue mapping, accessibility, connected pockets, and ranking. |
| `peptide_sampling.py` | Per-site budgets, adaptive rounds, selection, integrity checks, and rejection cleanup. |
| `process_control.py` | Logged subprocess execution and deadline/process-tree handling. |
| `pepglad_outputs.py` | Summary parsing, chain/sequence metadata, round IDs, and deletion journals. |
| `pepglad_validation.py` | OpenStructure/CCTBX measurements, policy, parallel workers, and reports. |
| `vina_ranking.py` | Valid-design submission, diagnostics, and Vina ranked CSVs. |
| `recover_pepglad_relaxation.py` | Inspection and optional recovery of interrupted relaxation. |
| `install.sh`, `environment.yml` | Pinned dependency setup and the main environment declaration. |
| `patches/` | Changes applied to the pinned external checkouts. |
| `vina_scoring/` | Canonical rigid-pose Meeko/Vina scorer, also imported by the comparison repository. |
| `tests/` | Focused regression tests and fixtures. |
| `docs/`, `mkdocs.yml` | This documentation site. |

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

The documentation needs only a separate Python environment; it does not need
GAPS, PepGLAD, Vina, or the scientific Conda dependencies:

```bash
python -m venv .venv-docs
.venv-docs/bin/python -m pip install -r requirements-docs.txt
.venv-docs/bin/python -m mkdocs build --strict
```

The static HTML site is written to `site/`. Both `site/` and `.venv-docs/` are
ignored by Git. To preview it:

```bash
.venv-docs/bin/python -m mkdocs serve
```

Open `http://127.0.0.1:8000/`. Search, navigation, code-copy buttons, and a
light/dark toggle are provided by Material for MkDocs. Publishing is a separate
step; the configuration does not presume a hosted site URL.

The dependency file constrains compatible MkDocs and Material major versions;
it is not a complete dependency lock. Record `pip freeze` when exact toolchain
reproduction is needed.

## Documentation maintenance

Use relative Markdown links between pages under `docs/`. Add new pages to
`mkdocs.yml`, and run the strict build to catch missing pages and anchors.

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
