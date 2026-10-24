# Troubleshooting

Start with the console error, `structure_sampling_summary.json`, and the affected
round's `pepglad.log` and `validation.log`. For scoring errors, inspect the
results CSV and the design-specific `interface.log`.

## Missing dependencies or environments

| Symptom | Action |
| --- | --- |
| `Vina scoring dependencies are unavailable` | Activate `allopep` and run `conda env update -n allopep -f environment.yml`; check that Meeko's `mk_prepare_receptor.py` is available. Validation-only skips scorer checks. |
| `OpenStructure ost executable not found` | Run the CLI from the `allopep` environment and check its installation. |
| Missing `ost`, `mmtbx`, `freesasa`, or `pymol` | Use `conda run -n allopep`; rerun the installer to apply the declared dependencies if needed. |
| Missing Ramachandran/rotamer data | Check `chem_data` in `allopep`, then rebuild the rotarama cache. |
| `ModuleNotFoundError: No module named 'api'` for a direct PepGLAD command | Run `python -m api.run` from the `PepGLAD/` checkout in the `PepGLAD` environment. The top-level CLI sets this working directory itself. |
| CUDA unavailable | Check GPU visibility and CUDA availability in both `gaps` and `PepGLAD`; also check OpenMM platforms. See [installation checks](installation.md#check-the-installation). |

Rebuild conformation tables with:

```bash
conda run --no-capture-output -n allopep mmtbx.rebuild_rotarama_cache
```

## No usable GAPS pockets

`No usable GAPS-guided receptor contexts found` means selection produced no
connected component of at least two eligible residues. Inspect the score
sidecar and receptor mapping. Eligible residues need GAPS score **>0.5**,
absolute SASA **>0 Å²**, and a Cα atom; links must be **≤8 Å**.

The CLI exposes a pocket-count limit, not threshold overrides. Older commands
using radius, relative-accessibility, blend, or `GuidanceSelection` settings do
not describe the current [selector](pocket-selection.md).

## Existing outputs or resume mismatch

| Message or condition | Meaning and next step |
| --- | --- |
| `Existing PepGLAD outputs require --resume` | Use `--resume` for the same run or a separate input filename stem for a fresh run. |
| Output contains a different input PDB | Two different inputs share a filename/output prefix. Preserve the existing output and use a distinct prefix. |
| `Saved sampling ... differs` | Input, generation code/checkpoint, settings, or validator no longer match. Reproduce the original settings or start a separate run. |
| `Saved sampling site count differs` | The selected pocket count changed, possibly through `--max_rank` or selector changes. |
| `Saved PepGLAD pocket differs for site N` | The selected residue identities or JSON format changed. Legacy three-field scored pockets do not equal the current two-field format. |
| `Saved sampling artifact changed` | A cached candidate, report, summary, or round record no longer matches its saved hash. |
| Already-charged attempts exceed the cap | Increase the cap to at least the charged count; resume does not reset attempts. |

Do not edit a checkpoint or deletion journal to force compatibility. Use
[resume checks](sampling-and-resume.md#resume-compatibility) to identify the
actual difference.

## Unrelaxed or malformed saved candidates

If resume reports unrelaxed samples, inspect the affected directory with
`python -m allopep.recover_relaxation` in the `PepGLAD` environment. Its default is
read-only; `--apply` quarantines malformed generated PDBs and relaxes pending
designs. The complete [recovery commands](sampling-and-resume.md#interrupted-openmm-relaxation)
include the per-design timeout.

Do not reinterpret overflowed or non-finite PDB coordinates as a normal
geometric outlier. They can prevent measurements and produce an error instead.

## No structures pass validation

Inspect the semicolon-separated reasons in the validation CSV and the detailed
JSONL. The [validation reference](validation.md#all-16-rejection-criteria)
explains every threshold. Reasons overlap, so summing their counts overstates
the number of failed structures.

Check `valid_count`, `attempted`, and `completion_reason` at each site. A
sampling exit status of `2` means the budget was exhausted; increasing the cap
allows more sampling but does not guarantee valid output. A positive validation
result is required before scoring, regardless of interface contacts or an
otherwise favorable energy.

## Missing rejected PDBs

The sampler intentionally removes rejected coordinates after saving their
validation records and hashes. Check the round's `removed_structures.json` and
validation reports. This is distinct from a missing valid candidate, which
breaks resume integrity. See [report retention](outputs.md#report-and-structure-retention).

## Round timeout or apparent stall

Use log timestamps and the round's reserved sample range to determine whether
generation, OpenMM, or validation was active. The round timeout covers all three
stages and terminates the process tree when exceeded. Resume retains the
charged range; recover unrelaxed saved candidates before continuing.

Only one sampling controller can use an output directory at a time. For
`Another sampling process is using ...`, check for the actual live process;
deleting `.sampling.lock` does not safely resolve a held advisory lock.

## Vina errors or missing ranked rows

Inspect the `error` and `warnings` fields in `structure_vina_results.csv`.
Common causes include chain mismatch, missing or unsupported residue templates,
Meeko preparation warnings, changed or lost heavy atoms, insufficient worker
memory, or non-finite `score_total`. Preparation errors fail the complete pose.

The ranked CSV includes only successful jobs. The sampling target counts valid
structures, so a scoring error does not schedule a replacement peptide. If a
valid extra was never selected during sampling, it also need not appear in the
manifest. See [Vina ranking](scoring.md).

## PyMOL launch or session failure

Use the launcher and Python module from the same environment:

```bash
conda run --no-capture-output -n allopep python -c "import pymol"
conda run --no-capture-output -n allopep \
  pymol -cq output/structure_output/structure_allosteric_site.pml
```

The second command checks the generated script without opening a GUI. A GUI
session additionally needs an available display. Updating `environment.yml`
alone does not update an existing environment; rerun the installer or apply
the environment update explicitly.
