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
│           ├── structure_pepglad_validation.jsonl
│           └── removed_structures.json
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
| Did each site reach its target? | `structure_sampling_summary.json`. |
| Why was a candidate rejected? | Its round validation CSV and JSONL. |
| Did generation or relaxation fail? | Round `summary.jsonl`, `pepglad.log`, and any per-design error records. |
| Why is a previously rejected PDB absent? | The round's `removed_structures.json` and rejection reports. |
| Which structures were selected? | Each site's `selected` entries in the sampling summary. |
| Which structures were scored? | `structure_vina_manifest.csv`. |
| Which score should determine rank? | `score_total` in `structure_vina_ranked.csv`. |
| Why is a valid structure absent from ranking? | Selection limit, then `structure_vina_results.csv` and its per-design error fields. |

## Sampling summary

The summary is both a report and a resume checkpoint. It contains the run
signature, overall status, per-site target and cap, reserved rounds, charged
attempts, valid counts, selected sample/path pairs, and failure counts.

Site completion reasons are:

| Reason | Meaning |
| --- | --- |
| `target_reached` | At least the requested number of valid structures are available. |
| `max_num_peptides` | The attempt cap was reached before the valid target. |
| `pending` | The site has neither reached its target nor exhausted its cap. |

`valid_count` can exceed the requested target because complete rounds are
validated. `selected` contains at most the first `-n` passing samples for each
site. Candidate attempts, accepted structures, selected structures, and
successful scores are different counts.

## Report and structure retention

During sampling, AlloPep retains round-level and cumulative validation records
for every checked complex, including invalid structures and valid extras.
Generation failures without a usable complex remain in generation summaries
and sampling failure counts.

After recording rejection details and hashes, the sampler removes rejected
PDBs, including rejected files in `invalid_designs/`. The round's
`removed_structures.json` journals each deleted relative path, status, and
SHA-256. Valid structures remain on disk. This journal allows resume to verify
an intentional rejection without treating the missing file as lost data.

The journal is a provenance record, not a coordinate backup. Use the retained
JSONL to inspect the original reported violations; a deleted PDB cannot be
reloaded for a new coordinate analysis.

!!! note "Existing-output modes refresh cumulative reports"

    `--validate-existing` and `--rank-existing` check available PDBs. Their
    cumulative reports describe that available cohort and can replace the
    cumulative reports from sampling. Round reports and deletion journals
    preserve the earlier checked cohort. Archive reports before reranking if
    the previous cumulative view is needed for an analysis.

Keep round directories, summaries, reports, journals, and the sampling checkpoint
together. Renaming candidates or editing cached reports can invalidate resume.

## Visualize selected sites

```bash
conda run --no-capture-output -n allopep \
  pymol output/structure_output/structure_allosteric_site.pse
```

The allosteric-site session and matching `.pml` script are created before peptide
sampling. They show the original receptor with GAPS scores from blue (low) to
red (high), with unscored residues gray, and surfaces for selected pockets in
the `sites` group. Residue identities, including insertion codes, control score
mapping.

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

Peptide carbons and their site's surface share a color. Other peptide atoms
use element colors. Fresh/resumed sampling shows its selected valid structures;
existing-output modes show the available structures that pass their recheck.

PyMOL coloring changes the visualization's B factors; it does not replace the
full-precision JSON score source used for selection.
