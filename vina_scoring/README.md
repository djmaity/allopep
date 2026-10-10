# Bundled rigid-pose Vina scorer

These Python modules are copied from
`compare_docking_scoring_function/programs/calculate_vina_score/` at Git commit
`169903500fd1fe78b9a1cf7b5687188d54b6a163`. The scoring and preparation
modules are unchanged. They run entirely from this directory; the source
checkout is not needed at runtime.

The bundled test file omits one source-repository test that depended on an
external cyclic-peptide fixture absent from AlloPep.

The scorer prepares receptor and peptide PDBQT files with Meeko, keeps the
input pose rigid, scores residue-sized peptide fragments against one receptor
map, and sums their native Vina energies before rounding. Each batch job runs
in a disposable process. See `vina_score.py --help` for standalone usage.

AlloPep supplies a manifest of PepGLAD complexes that passed structural
validation. It ranks successful `score_total` values from lowest to highest
and retains the full scorer result CSV, including per-sample failures.
