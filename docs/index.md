# AlloPep

AlloPep turns a receptor PDB into candidate peptide–receptor complexes. It uses
GAPS to identify candidate sites, PepGLAD to generate peptides, OpenStructure
and CCTBX to check their physical validity, and the bundled Vina scorer to
rank the accepted structures.

## From receptor to ranked designs

| Stage | What happens | Main result |
| --- | --- | --- |
| 1. Predict | GAPS assigns scores to receptor residues. | Full-precision scores and a score-colored PDB. |
| 2. Select | Accessible, high-scoring residues form connected pockets. | Ranked pocket definitions. |
| 3. Generate | PepGLAD samples peptides, idealizes geometry, and applies OpenMM relaxation. | Candidate complexes and generation records. |
| 4. Validate | Each usable complex must pass all 16 physical-validity criteria. | Accepted structures and rejection details. |
| 5. Continue | Each site receives more samples until its valid-structure target or attempt cap is reached. | A resumable sampling summary. |
| 6. Rank | Vina scores the selected valid structures. | A CSV ordered by increasing `score_total`. |

!!! note "What a result establishes"

    A passing structure satisfies the implemented physical checks. A favorable
    Vina score is a ranking signal. Experimental binding, allosteric activity,
    and dynamic stability require separate evidence. See
    [Scientific interpretation](scientific-scope.md).

## Start a run

After [installation](installation.md), run these commands from the repository
root, replacing `input/structure.pdb` with your receptor PDB:

```bash
conda activate allopep
python allopep.py input/structure.pdb
```

The defaults use every selected pocket, seek **10 physically valid structures
per site**, and permit **10,000 candidate attempts per site**. For a smaller
initial run:

```bash
python allopep.py input/structure.pdb --max_rank 1 -n 1 --max_num_peptides 32
```

This bounds the candidate count; it does not guarantee a passing structure.
Results are written under `output/structure_output/`.

## Find the right guide

| Task | Guide |
| --- | --- |
| Install dependencies and check the environments | [Installation](installation.md) |
| Choose lengths, sites, and attempt budgets | [Usage and CLI](usage.md) |
| Understand how a pocket is selected | [Pocket selection](pocket-selection.md) |
| Continue an interrupted or incomplete run | [Sampling and resume](sampling-and-resume.md) |
| Interpret every rejection reason and threshold | [Structure validation](validation.md) |
| Understand the ranking score and its provenance | [Vina ranking](scoring.md) |
| Find reports, structures, and PyMOL sessions | [Outputs and visualization](outputs.md) |
| Diagnose a failed run | [Troubleshooting](troubleshooting.md) |
| Maintain the code or build this documentation | [Development](development.md) |

## Scope of this documentation

These pages describe the current checkout's connected-pocket selector,
bounded sampling controller, and physical-validity policy. Historical pocket
optimization reports describe different selector settings; they do not
validate the current algorithm. The [scientific interpretation
guide](scientific-scope.md#pocket-recovery-evidence) explains that boundary.
