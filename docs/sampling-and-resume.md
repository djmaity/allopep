# Sampling and resume

AlloPep samples toward a target number of **physically valid structures at each
site**. The target is distinct from the number of candidates attempted or the
number of successful Vina scores.

## Adaptive rounds

Unfinished sites receive one generation round each in GAPS-rank order. Each
round launches a separate PepGLAD process with GPU minibatches of eight.

| Condition | Next round |
| --- | --- |
| No previous attempts | Start with eight candidates. |
| No passing structures yet | Double the previous round size. |
| Some passing structures | Estimate remaining demand from the cumulative acceptance rate, rounded to minibatches. |
| Growth and size limits | Grow by at most a factor of two; normally stay between eight and 64 candidates. |
| Budget nearly exhausted | Shrink to the remaining number of attempts, even below eight. |
| Target reached or cap exhausted | Schedule no further rounds for that site. |

The controller selects the first `-n` valid structures in per-site sample
order. A final round can overshoot the target; additional valid structures
remain on disk and in validation reports, but are not selected for that run's
ranking or peptide visualization.

## Generation and relaxation

The patched PepGLAD API runs in co-design mode using `codesign.ckpt`. It
idealizes generated backbone and side-chain geometry before exporting PDBs,
and screens generated coordinates for finite, representable values.

Usable complexes then undergo PDBFixer preparation and restrained CHARMM36
OpenMM minimization. Individual minimization errors are recorded and separated
from successful candidates. An infrastructure failure, such as a broken Ray
runtime, can still fail the round.

Idealization and minimization do not bypass validation. The controller checks
the OpenMM creation remark before accepting a saved complex for validation.

## Attempts are charged before launch

The full sample range for a round is reserved and saved before its subprocess
starts. Every reserved candidate consumes an attempt, including generation
failures, relaxation failures, rejected structures, and candidates never
written during an interrupted round.

Sample IDs in validation and ranking reports are unique within each site.
Round-local filenames can repeat: `structure_0.pdb` in the second round is a
different sample from `structure_0.pdb` in the first round. Use the recorded
`site`, `sample`, and full path to identify a candidate.

The sampling summary records these failure categories:

| Category | Meaning |
| --- | --- |
| `generation_invalid` | PepGLAD marked a candidate invalid before producing a usable complex. |
| `relaxation_error` | PepGLAD recorded an individual relaxation error. |
| `validation_invalid` | A usable complex failed the physical-validity policy. |
| `validation_error` | A usable complex could not be checked successfully. |
| `unreported` | Reserved attempts have no generation summary record. |

Vina errors occur after selection and do not trigger replacement sampling.

## Resume compatibility

```bash
python -m allopep input/structure.pdb --resume
```

Use the original generation settings, including any nondefault length or pocket
limit. The controller checks:

- Input PDB SHA-256.
- Peptide length bounds and peptide-chain override.
- Validator, PepGLAD patch, API, model source, and co-design checkpoint hashes.
- Selected site count and exact two-field pocket definitions.
- Reserved sample ranges, round metadata, saved reports, and candidate hashes.

Already-validated rounds are reused after integrity checks. An interrupted
round keeps its full reserved budget; resume processes available relaxed
candidates without generating replacements within that old range.

Designs live only in `round_NNNNNN` subdirectories of a site directory. The
earlier flat layout, which wrote PDBs and `summary.jsonl` directly into the
site directory, is no longer adopted; such a directory is reported as an
error instead.

Output directories created before the pipeline became the `allopep` package
cannot be resumed. The saved signature includes a hash of the validation
module, which moved and changed, so resume reports differing validation
policy. Use a fresh output prefix for those inputs. `--validate-existing`
and `--rank-existing` never read the sampling summary and still work on
older output directories.

You may adjust `-n`, `--attempt-multiplier`, and `--max_num_peptides`. The
effective cap must remain at least `-n` and cannot be below attempts already
charged at any site. Budgets accumulate across resumes.

Deleted rejected structures are recognized through their rejection records and
hash journal. Valid structures must still exist and match their saved hashes.
See [report and structure retention](outputs.md#report-and-structure-retention).

## Interrupted OpenMM relaxation

Inspect the affected round from the repository root:

```bash
conda run --no-capture-output -n PepGLAD \
  python -m allopep.recover_relaxation \
  output/structure_output/structure_PepGLAD_outputs/structure_allosteric_site_1_codesign/round_000001
```

Inspection is read-only by default. To quarantine malformed generated PDBs and
relax pending usable designs, rerun with `--apply`:

```bash
conda run --no-capture-output -n PepGLAD \
  python -m allopep.recover_relaxation \
  output/structure_output/structure_PepGLAD_outputs/structure_allosteric_site_1_codesign/round_000001 \
  --apply --timeout-seconds 180
```

The recovery timeout is **per design**, with a default of 180 seconds. Use the
site directory instead for legacy flat outputs. After recovery, rerun AlloPep
with `--resume` and the original settings.

## Timeouts and concurrent runs

`--round-timeout-seconds` covers generation, relaxation, and validation together
for one new round. On timeout, AlloPep terminates the subprocess and its
descendants, saves progress, and exits with an operational error. This is
separate from Vina scoring, which runs after selection in disposable workers.

An advisory `.sampling.lock` prevents two sampling controllers from using the
same output directory simultaneously. A leftover lock file is not evidence of
a live lock; the operating system releases the advisory lock when its holder
exits.

Implementation: `allopep/sampling.py`, `allopep/process_control.py`,
`allopep/outputs.py`, and `allopep/recover_relaxation.py`.
