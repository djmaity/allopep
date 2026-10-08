""" https://www.blopig.com/blog/2023/05/checking-your-pdb-file-for-clashing-atoms/
Author: Brennan Abanades Kenyon
Modified by Dibyajyoti Maity
"""

import argparse

import numpy as np
from Bio import PDB


# Atomic radii for supported elements.
atom_radii = {
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "S": 1.80,
    "F": 1.47,
    "P": 1.80,
    "CL": 1.75,
    "MG": 1.73,
}


def argument_parser():
    parser = argparse.ArgumentParser(
                        prog='count_clashes.py',
                        description='Count the number of clashes in a protein structure.')
    parser.add_argument('pdb_file', help='Input PDB file')
    return parser.parse_args()


def count_clashes(pdb_file, clash_cutoff=0.63):
    pdb_parser = PDB.PDBParser(QUIET=True)
    structure = pdb_parser.get_structure('structure', pdb_file)

    # Only atoms with known radii can be checked.
    atoms = [x for x in structure.get_atoms() if x.element in atom_radii]
    coords = np.array([a.coord for a in atoms], dtype="d")
    kdt = PDB.kdtrees.KDTree(coords)
    search_radius = 2 * clash_cutoff * max(atom_radii.values())

    clashes = 0
    for index, atom_1 in enumerate(atoms):
        for neighbor in kdt.search(coords[index], search_radius):
            if neighbor.index <= index:
                continue
            atom_2 = atoms[neighbor.index]

            # Exclude clashes from atoms in the same residue
            if atom_1.parent.id == atom_2.parent.id:
                continue

            # Exclude clashes from peptide bonds
            if {atom_1.name, atom_2.name} == {"C", "N"}:
                continue

            # Exclude clashes from disulphide bridges
            if atom_1.name == atom_2.name == "SG" and neighbor.radius > 1.88:
                continue

            cutoff = clash_cutoff * (
                atom_radii[atom_1.element] + atom_radii[atom_2.element]
            )
            if neighbor.radius < cutoff:
                clashes += 1

    return clashes


def main():
    args = argument_parser()
    print(count_clashes(pdb_file=args.pdb_file))


if __name__ == '__main__':
    main()
