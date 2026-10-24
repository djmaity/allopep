# Outputs and visualization

An input named `input/structure.pdb` uses the output prefix `structure` and
writes under `output/structure_output/`. Output naming depends on the filename
stem rather than the full input path. A different input with the same basename
is rejected if it conflicts with the saved receptor.

## Directory map

The following layout is schematic; optional files depend on how far the run
progressed:

```text
output/structure_output/
├── structure.pdb
├── structure_GAPS_scores.json
├── structure_GAPS_output.pdb
├── structure_allosteric_sites/
│   └── structure_allosteric_site_1.json
├── structure_allosteric_site.pml
├── structure_allosteric_site.pse
├── structure_sampling_summary.json
├── structure_peptides.csv
├── structure_valid_complexes/
│   └── site_1/
│       └── structure_site_1_peptide_0.pdb
├── structure_PepGLAD_outputs/
│   └── structure_allosteric_site_1_codesign/
│       └── round_000001/
│           ├── round.json
│           ├── summary.jsonl
│           ├── structure_0.pdb
│           ├── pepglad.log
│           ├── validation.log
│           ├── structure_pepglad_validation_manifest.json
│           ├── structure_pepglad_validation.csv
│           └── structure_pepglad_validation.jsonl
├── structure_pepglad_aligned.pse
└── structure_vina_score/
    ├── structure_pepglad_validation_manifest.json
    ├── structure_pepglad_validation.csv
    ├── structure_pepglad_validation.jsonl
    ├── structure_vina_manifest.csv
    ├── structure_vina_results.csv
    └── structure_vina_ranked.csv
```

## Which artifact answers which question?

| Question | Inspect |
| --- | --- |
| What receptor was used? | Saved `structure.pdb` and input hashes in reports/checkpoint. |
| Which residues were scored or selected? | `structure_GAPS_scores.json` and the pocket JSON files. |
| Was the overall or per-site target reached? | `structure_sampling_summary.json`. |
| Why was a candidate rejected? | Its round validation CSV and JSONL. |
| Did generation or relaxation fail? | Round `summary.jsonl`, `pepglad.log`, and any per-design error records. |
| Which structures were selected? | `structure_peptides.csv`, the valid-complex copies, and each site's `selected` checkpoint entries. |
| Which structures were scored? | `structure_vina_manifest.csv`. |
| Which score should determine rank? | `score_total` in `structure_vina_ranked.csv`. |
| Why is a valid structure absent from ranking? | Selection limit, then `structure_vina_results.csv` and its per-design error fields. |

## Sampling summary

The summary is both a report and a resume checkpoint. It contains the run
signature, overall status, `peptides_per_site` setting, target and cap,
reserved rounds, total and per-site charged attempts and valid counts,
selected sample/path pairs, and failure counts.

Site completion reasons are:

| Reason | Meaning |
| --- | --- |
| `target_reached` | At least the requested number of valid structures are available. |
| `max_num_peptides` | The attempt cap was reached before the valid target. |
| `pending` | The site has neither reached its target nor exhausted its cap. |

`valid_count` can exceed the requested target because complete rounds are
validated. `selected` contains at most the first `-n` passing samples for each
site when `--peptides-per-site` is set, or at most `-n` across all sites
otherwise. Candidate attempts, accepted structures, selected structures,
and successful scores are different counts. In total mode, site completion
reasons describe the shared target and budget.

## Selected-complex CSV

`structure_peptides.csv` has one row per selected complex in site and sample
order. The copies retain sample IDs, including gaps between IDs, and contain
the complete receptor and peptide. Source and copy paths are relative to the
output directory. Repeating a command refreshes the copies to match its
selection, preserving original predictions.

Site residues use semicolon-separated identifiers such as
`A-ARG-42;A-GLY-43;B-TYR-105A`. Each includes chain, three-letter residue name,
and residue number with any insertion code. `gaps_pocket_score` is their mean
full-precision GAPS score. `pocket_sasa_a2` sums their absolute FreeSASA residue
areas measured in the original receptor before adding the generated peptide.

Other columns give the sequence, length, chains, validation status, Vina score,
units, rank, status and error, buried SASA, receptor contacts at 5 Å, and minimum
interface distance. Scoring failures keep a row with blank score and rank;
validation-only runs use `vina_status=not_run`. No passing complexes produces
a header-only CSV and no copied PDBs.

## Report and structure retention

During sampling, AlloPep retains round-level and cumulative validation records
for every checked complex, including invalid structures and valid extras.
Generation failures without a usable complex remain in generation summaries
and sampling failure counts.

The sampler keeps rejected PDBs and files in `invalid_designs/`, along with
their validation details, summaries, errors, and logs. Selected valid complexes
are copied separately; final-round valid extras also remain in their source
directories.

> **Existing-output modes refresh cumulative reports**
>
> `--validate-existing` and `--rank-existing` check available PDBs. Their
> cumulative reports describe that available cohort and can replace the
> cumulative reports from sampling. Round reports
> preserve the earlier checked cohort. Archive reports before reranking if
> the previous cumulative view is needed for an analysis.

Keep round directories, summaries, reports, and the sampling checkpoint
together. Renaming candidates or editing cached reports can invalidate resume.

## Visualize selected sites

```bash
conda run --no-capture-output -n allopep \
  pymol output/structure_output/structure_allosteric_site.pse
```

The allosteric-site session and matching `.pml` script are created before peptide
sampling. They load the saved original receptor and apply full-precision scores
from `structure_GAPS_scores.json`, colored from blue (0) to red (1), with
unscored residues gray. Selected pocket surfaces appear in the `sites` group.
Residue identities, including insertion codes, control score mapping.
`--validate-existing` and `--rank-existing` also refresh these files using
the saved pocket definitions.

## Visualize selected valid peptides

```bash
conda run --no-capture-output -n allopep \
  pymol output/structure_output/structure_pepglad_aligned.pse
```

The original input receptor is named `structure`. Pocket surfaces use names
such as `site_1`. The `peptides` group contains `site_1_peptides` and the other
per-site groups, with objects such as `site_1_peptide_1`. Complexes are aligned
using receptor Cα atoms; only their peptide chains are shown alongside the
original receptor.

Peptide carbons and their site's surface share a color, following PyMOL's
default object color sequence starting with green for site 1. Other peptide
atoms use element colors. Fresh/resumed sampling shows its selected valid
structures; existing-output modes show the requested selection after
rechecking candidates.

PyMOL coloring changes the visualization's B factors; it does not replace the
full-precision JSON score source used for selection.
