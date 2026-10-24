# 1a08 score regression fixture

These selected receptor and peptide mmCIF files were copied from
`compare_docking_scoring_function/data/processed/fixed_structures/1a08/`.
The source SHA-256 digests are:

- `1a08_protein.cif`: `b3751e7ed10973067d8ded83e164771043b1198a1e75b59c382c758f9a11c354`
- `1a08_peptide.cif`: `fd7070c1b6ea031c6b178d2e5f33d00cb69655f588e98921c61eb8abaeb8bdeb`

`1a08_complex.cif` holds both of them in one file, receptor chain `A` and
peptide chain `C`, in their original coordinate frame. It was written with
Gemmi by appending the peptide chain to the receptor structure, and the
scorer takes combined structures only:

- `1a08_complex.cif`: `0bbc2db21b522f1468d6d312cbb9f180be1df98d32d722dc7380d35bc703eb01`

With Meeko 0.8.0 and Vina 1.2.7, the rigid-pose total and ligand interaction
scores are both -5.640; the other six components are zero. The peptide has
27 atoms, 22 heavy atoms, two residues, and zero active torsions.

Only `1a08_complex.cif` is used by the tests, because the scorer takes
combined structures only. The two separate chain files are kept as the
provenance of the combined fixture.
