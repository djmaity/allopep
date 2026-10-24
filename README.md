# AlloPep

AlloPep finds candidate allosteric sites with GAPS, designs peptides with
PepGLAD, checks their physical validity, and ranks valid protein–peptide
complexes with rigid-pose Vina scoring.

## Quick start

Requires Linux, Conda, and a CUDA-capable GPU. See
[installation](docs/installation.md) for prerequisites and environment setup.

```bash
git clone https://github.com/djmaity/allopep.git
cd allopep
bash install.sh
conda activate allopep
python -m allopep input/structure.pdb -n 10
```

`-n` requests valid complexes across all selected sites. Add
`--peptides-per-site` to request that many at each site.

Results are saved in `output/<input-stem>_output/`, including the peptide CSV,
selected valid complexes, validation reports, and PyMOL sessions.

## Documentation

See the [full documentation](docs/index.md) for detailed instructions.

- [Installation](docs/installation.md)
- [Usage and CLI](docs/usage.md) · [Sampling and resume](docs/sampling-and-resume.md)
- [Pocket selection](docs/pocket-selection.md) · [Validation](docs/validation.md) · [Vina ranking](docs/scoring.md)
- [Outputs and visualization](docs/outputs.md) · [Troubleshooting](docs/troubleshooting.md)
- [Scientific interpretation](docs/scientific-scope.md) · [Development](docs/development.md)

## License

[GNU GPL v3](LICENSE).
