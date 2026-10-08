""" Calculate the mean GAPS score of APOP pockets """

import io
import os
import re
import statistics
import zipfile

import Bio.SeqIO
from Bio.PDB.PDBParser import PDBParser

pocket_pattern = re.compile(r"Rank=(\d+)\nPocket name: (\S+)\nAPOP score: (\S+)\nResidues: (.+)")
residue_pattern = re.compile(r"(\d+)([a-zA-Z]+)")


def first_residue_by_chain(pdb_file):
    return {
        chain.id.split(':')[1]: chain.annotations['start']
        for chain in Bio.SeqIO.parse(pdb_file, 'pdb-atom')
    }


def pocket_mean_gaps_score(input_file, apop_file, gaps_file):
    """ Calculate the mean GAPS score of APOP pockets """
    prefix = os.path.splitext(os.path.basename(input_file))[0]

    input_first_residue = first_residue_by_chain(input_file)
    gaps_first_residue = first_residue_by_chain(gaps_file)

    parser = PDBParser(QUIET=True)
    gaps_structure = parser.get_structure(id=prefix, file=gaps_file)

    pockets = []
    with zipfile.ZipFile(apop_file) as apop_output_zip:
        with apop_output_zip.open('apop_output.txt', 'r') as apop_output:
            apop_data = apop_output.read().decode("utf-8")

        for _, pocket_filename, apop_score, residues in re.findall(pocket_pattern, apop_data):
            apop_pocket_name = os.path.splitext(pocket_filename)[0]
            with apop_output_zip.open(pocket_filename, 'r') as pocket_file:
                pdb_data = pocket_file.read().decode("utf-8")
                apop_pocket = parser.get_structure(
                    id=apop_pocket_name, file=io.StringIO(pdb_data))

                b_factors = []
                for atom in apop_pocket.get_atoms():
                    residue = atom.get_parent()
                    chain = residue.get_parent()
                    model = chain.get_parent()

                    offset = input_first_residue[chain.id] - gaps_first_residue[chain.id]
                    gaps_atom = gaps_structure[model.id][chain.id][residue.id[1] - offset][atom.id]
                    b_factors.append(gaps_atom.get_bfactor())
                gaps_score = statistics.mean(b_factors)

            pocket = [
                [chain, [int(resid), " "]]
                for resid, chain in re.findall(residue_pattern, residues)
            ]
            pockets.append([pocket_filename, pocket, apop_score, gaps_score])

    pockets.sort(key=lambda x: x[3], reverse=True)
    return pockets


if __name__ == '__main__':
    pocket_scores = pocket_mean_gaps_score(
        input_file='./input/2VH7.pdb',
        apop_file='output/2VH7_output/2VH7_apop_output.zip',
        gaps_file='output/2VH7_output/2VH7_GAPS_output.pdb'
        )

    from pprint import pprint
    pprint(pocket_scores)
