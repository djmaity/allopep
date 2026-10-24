# Usage and CLI

Run commands from the AlloPep repository root with the `allopep` environment
active. AlloPep accepts a receptor structure in PDB format; provide the receptor
chains intended for design. This CLI does not expose input cleaning or native
peptide removal.

## Start a bounded design run

```bash
conda activate allopep
python -m allopep input/structure.pdb \
  --max_rank 1 \
  --length_min 8 --length_max 13 \
  --num_samples 5 --max_num_peptides 500
```

This designs against the highest-ranked pocket, samples lengths **8–12**, and
seeks five physically valid structures with at most 500 candidate attempts.
Results use the input filename stem: `structure.pdb` maps to
`output/structure_output/`.

!!! note "Peptide length uses an exclusive upper bound"

    PepGLAD samples `length_min <= length < length_max`. The defaults `5` and
    `10` therefore generate lengths **5–9**. To request only length eight, use
    `--length_min 8 --length_max 9`. Keep the upper bound strictly greater than
    the lower bound.

## Command-line reference

```text
python -m allopep PDB_FILE [OPTIONS]
```

| Argument | Default | Meaning |
| --- | --- | --- |
| `pdb_file` | Required | Receptor PDB path; its stem determines the output directory. |
| `-r`, `--max_rank` | `all` | Use every ranked GAPS pocket, or the first positive integer number of pockets. |
| `-l`, `--length_min` | `5` | Minimum generated peptide length, inclusive. |
| `-L`, `--length_max` | `10` | Maximum length, exclusive in PepGLAD. |
| `-n`, `--num_samples` | `10` | Target number of physically valid structures **per site**. |
| `--attempt-multiplier` | `1000` | Positive integer multiplied by `-n` to calculate each site's attempt cap. |
| `--max_num_peptides` | Calculated | Absolute per-site attempt cap; overrides the multiplier and must be at least `-n`. |
| `--round-timeout-seconds` | `3600` | Shared deadline for one round's generation, relaxation, and validation. |
| `-c`, `--peptide-chain` | From summary | Known peptide chain ID; must agree with PepGLAD's summary. One alphanumeric character. |
| `-v`, `--validate-existing` | Off | Validate saved available complexes without generation or Vina scoring. |
| `-R`, `--rank-existing` | Off | Validate and rank saved available complexes without generation. |
| `--resume` | Off | Reuse compatible saved rounds and continue incomplete sites. |
| `-h`, `--help` | — | Display help and exit. |

`--validate-existing`, `--rank-existing`, and `--resume` are mutually exclusive.
The target, multiplier, cap, and round timeout must be positive integers.

## Choose the attempt budget

Without an explicit cap:

```text
per-site attempt cap = num_samples × attempt_multiplier
```

For example, `-n 10 --attempt-multiplier 200` permits 2,000 candidates at each
selected site. `-n 10 --max_num_peptides 1000` instead sets a cap of 1,000.
The budget applies separately to each site, so designing more pockets can
increase total work substantially.

Every candidate consumes an attempt, including invalid coordinates, relaxation
failures, and rejected structures. Different structures with the same peptide
sequence count separately. See [Sampling and resume](sampling-and-resume.md)
for adaptive rounds, reservation accounting, and selected structures.

## Continue sampling

Use the same input and generation settings:

```bash
python -m allopep input/structure.pdb --max_rank 1 \
  --length_min 8 --length_max 13 -n 5 --max_num_peptides 1000 --resume
```

You may change the valid-structure target or attempt budget, but the cap cannot
be less than attempts already charged. Changed pockets, incompatible generation
settings, or altered saved artifacts cause resume to fail. See the
[compatibility rules](sampling-and-resume.md#resume-compatibility).

## Recheck or rerank saved structures

```bash
python -m allopep input/structure.pdb --validate-existing
python -m allopep input/structure.pdb --rank-existing
```

These modes discover available structures in both legacy flat outputs and round
directories. They check every available candidate, including retained valid
extras, rather than applying `-n` as a new selection limit. They do not generate
replacement samples. Rejected structures already deleted by the sampler remain
documented in round reports; their coordinates cannot be rechecked. See
[report retention](outputs.md#report-and-structure-retention).

## Exit status

| Status | Meaning |
| --- | --- |
| `0` | Normal completion. For sampling, every site reached its valid-structure target. Validation-only can finish normally with zero passing structures; inspect the reports. |
| `2` | At least one site exhausted its attempt cap before reaching the target. Partial valid results are retained and ranked when available. |
| `1` | Operational failure, such as a timeout, incompatible resume state, failed dependency, or all Vina jobs failing. |

A partial Vina scoring failure is recorded in the results CSV and does not
cause replacement sampling. Exit status alone does not give the number of
valid structures or successful scores.

Implementation: `allopep/cli.py::parse_args()`,
`allopep/cli.py::sampling_config()`, and `allopep/cli.py::main()`.
