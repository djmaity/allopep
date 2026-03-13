#! /usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
import statistics
import subprocess

import numpy as np
from Bio.PDB import PDBIO, PDBParser, Select
from vina import Vina


def argument_parser():
    parser = argparse.ArgumentParser(
        prog='vina_score.py',
        description='Calculate the AutoDock Vina score between one chain and others in a PDB file',
        epilog='The mean of the AutoDock Vina scores is given in kcal/mol.'
               'Multiple runs of openbabel results in minor charge differences'
               'in the third decimal places, which results in differences'
               ' in the Vina Score calculation.')
    parser.add_argument('pdb_file', help='Input PDB file')
    parser.add_argument('-c', '--chain', dest='peptide_chain', type=str,
                        help='Peptide chain name', required=True)
    parser.add_argument('-b', '--box_padding', default=5, type=float,
                        help='Addtional space in angstroms in each dimension '
                             'of the box around the peptide during '
                             'computation of AutoDock Vina map')
    parser.add_argument('-n', '--num_iter', dest='num_iterations', default=1, type=int,
                        help='Number of times to generate the PDBQT files and '
                            'calculate scores')
    return parser.parse_args()


def get_vina_score(protein_filename, peptide_filename, index, center, box_size, pH=7.4):
    subprocess.run(['obabel', f'{protein_filename}.pdb',
                    '-O', f'{protein_filename}_H_{index}.pdb', '-xr', '-p', str(pH)],
                    capture_output=True)
    subprocess.run(['obabel', f'{protein_filename}_H_{index}.pdb',
                    '-O', f'{protein_filename}_{index}.pdbqt',
                    '--partialcharge', 'eem', '-xr'],
                    capture_output=True)

    subprocess.run(['obabel', f'{peptide_filename}.pdb',
                    '-O', f'{peptide_filename}_H_{index}.pdb', '-p', str(pH), '-xr'],
                    capture_output=True)
    subprocess.run(['obabel', f'{peptide_filename}_H_{index}.pdb', '-O',
                    f'{peptide_filename}_{index}.pdbqt', '--partialcharge', 'eem'],
                    capture_output=True)

    energy = None
    try:
        v = Vina(sf_name='vina', verbosity=0)
        v.set_receptor(f'{protein_filename}_{index}.pdbqt')
        v.set_ligand_from_file(f'{peptide_filename}_{index}.pdbqt')
        v.compute_vina_maps(center=center, box_size=box_size)
        # Score the current pose
        energy = v.score()[0]
    except TypeError:
        print('AutoDock Vina score calculation failed')

    os.remove(f'{protein_filename}_H_{index}.pdb')
    os.remove(f'{protein_filename}_{index}.pdbqt')
    os.remove(f'{peptide_filename}_H_{index}.pdb')
    os.remove(f'{peptide_filename}_{index}.pdbqt')
    return energy


def vina_score(pdb_file, peptide_chain, box_padding=5, num_iterations=1, pH=7.4):
    basename, _ = os.path.splitext(os.path.basename(pdb_file))
    protein_filename = f'{basename}_protein'
    peptide_filename = f'{basename}_peptide'

    pdb_parser = PDBParser(QUIET=True)
    structure = pdb_parser.get_structure(basename, pdb_file)


    class SelectPeptideChain(Select):
        def accept_chain(self, chain):
            return chain.get_id() == peptide_chain

    class SelectProteinChain(Select):
        def accept_chain(self, chain):
            return chain.get_id() != peptide_chain

    io = PDBIO()
    io.set_structure(structure)
    io.save(f'{peptide_filename}.pdb', SelectPeptideChain())
    io.save(f'{protein_filename}.pdb', SelectProteinChain())

    atoms = structure[0][peptide_chain].get_atoms()
    coordinates = np.array([atom.get_coord() for atom in atoms])
    center = [float(coordinate) for coordinate in
              (coordinates.max(axis=0) + coordinates.min(axis=0))/2]
    box_size = [float(coordinate) for coordinate in
                coordinates.max(axis=0) - coordinates.min(axis=0) + box_padding]

    scores = []
    for index in range(num_iterations):
        """ Multiple runs of openbabel results in minor charge differences
            in the third decimal places, which results in differences in
            the Vina Score calculation.
        """
        score = get_vina_score(protein_filename=protein_filename,
                               peptide_filename=peptide_filename,
                               index=index, center=center,
                               box_size=box_size, pH=pH)
        if score is not None:
            scores.append(score)

    os.remove(f'{protein_filename}.pdb')
    os.remove(f'{peptide_filename}.pdb')

    # Average energies in kcal/mol
    if len(scores) > 1:
        return basename, statistics.fmean(scores), statistics.stdev(scores)
    elif len(scores) == 1:
        return basename, float(round(scores[0], 3))
    else:
        return


def main():
    args = argument_parser()
    score = vina_score(pdb_file=args.pdb_file,
                       peptide_chain=args.peptide_chain,
                       box_padding=args.box_padding,
                       num_iterations=args.num_iterations)
    print(score)


if __name__ == '__main__':
    main()
