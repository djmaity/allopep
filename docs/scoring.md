# Rosetta ranking

AlloPep ranks selected physically valid complexes using an external Rosetta
InterfaceAnalyzer executable. **Lower `rosetta_dG_separated` ranks first.**
Validation status is established before scoring.

## Scoring protocol

The scorer constructs an explicit receptor–peptide interface from PepGLAD's
metadata. Receptor chains A and B with peptide chain C, for example, become
`AB_C`. Coordinate chain IDs must match those expected from the summary.

For each complex, the command enables interface shape complementarity and
packing statistics, uses a constant random seed of 42, and specifies:

```text
-pack_input false
-pack_separated false
```

The calculation does not relax the complex or repack its bound or separated
states. It scores the saved coordinates after PepGLAD's earlier preparation
and relaxation. Each design has an isolated run directory and a scoring
timeout of 1,800 seconds in the current implementation; the top-level CLI does
not expose that timeout.

## Weight provenance

| Property | Value |
| --- | --- |
| Bundled weights | `weights/optimized_rosetta_26_terms.wts` |
| Selected fit | 26-term ridge weights |
| Origin | `compare_docking_scoring_function/output/rosetta_affinity/optimized_rosetta_26_terms.wts` |
| SHA-256 | `e193a7acc0cfe0918a15662cac9d021f232797da81d57aa89e6a1c6005983044` |
| Ranking field | `rosetta_dG_separated` |
| Ties | Site, then per-site sample ID. |

The scorer verifies the bundled file's checksum. Changing the weights is a
change to the scoring protocol and provenance contract; a silent replacement
is rejected.

The fitted score is an empirical ranking signal. Its name does not make it a
measured binding free energy, an affinity estimate in kcal/mol, or a validated
predictor for generated peptides. See [Scientific interpretation](scientific-scope.md).

## Results and failures

Files are written under `output/structure_output/structure_rosetta_score/`:

| File | Meaning |
| --- | --- |
| `structure_rosetta_manifest.csv` | Valid selected designs submitted for scoring, with chains, binary path, and weight checksum. |
| `structure_rosetta_results.csv` | Every submitted design, with `status`, error text, available metrics, and raw artifact paths. |
| `structure_rosetta_ranked.csv` | Only successfully scored designs, with rank starting at one. |
| `rosetta_runs/site_N_sample_M/interface.sc` | Raw scorefile for a design. |
| `rosetta_runs/site_N_sample_M/interface.log` | Process output and diagnostics. |

A successful result requires a complete score row and a finite numeric
`dG_separated`. Exit code zero without a scorefile is a scoring error. Failed
jobs remain in the results CSV and are excluded from ranking; all scoring jobs
failing causes an operational failure. Individual scoring failures do not
trigger additional generation.

## Rerank available complexes

```bash
conda activate allopep
python allopep.py input/structure.pdb --rank-existing \
  --rosetta-bin /path/to/InterfaceAnalyzer
```

This validates all available saved candidates and scores those passing the
current gate. It can include valid extras retained from a final sampling round.
It refreshes cumulative validation and scoring reports, so archive reports
first when preserving a previous analysis. No new peptides are generated.

Old Vina directories are historical artifacts; the current pipeline does not
convert them into Rosetta results.

Implementation: `rosetta_scoring.py` and the ranking branch of `allopep.py`.
