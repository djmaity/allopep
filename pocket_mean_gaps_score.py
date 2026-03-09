""" Calculate the mean GAPS score of APOP pockets """

import io
import json
import os
import re
import statistics
import zipfile

import Bio.Align
import Bio.SeqIO
from Bio.PDB.PDBParser import PDBParser

pocket_pattern = re.compile(r"Rank=(\d+)\nPocket name: (\S+)\nAPOP score: (\S+)\nResidues: (.+)")
residue_pattern = re.compile(r"(\d+)([a-zA-Z]+)")


def pocket_mean_gaps_score(input_file, apop_file, gaps_file):
    """ Calculate the mean GAPS score of APOP pockets """
    prefix = os.path.splitext(os.path.basename(input_file))[0]

    input_first_residue = {}
    for input_chain in Bio.SeqIO.parse(input_file, "pdb-atom"):
        input_chain_id = input_chain.id.split(':')[1]
        input_first_residue[input_chain_id] = input_chain.annotations['start']

    gaps_first_residue = {}
    for gaps_chain in Bio.SeqIO.parse(gaps_file, "pdb-atom"):
        gaps_chain_id = gaps_chain.id.split(':')[1]
        gaps_first_residue[gaps_chain_id] = gaps_chain.annotations['start']

    parser = PDBParser(QUIET=True)
    gaps_structure = parser.get_structure(id=prefix, file=gaps_file)

    pockets = []
    with zipfile.ZipFile(apop_file) as apop_output_zip:
        with apop_output_zip.open('apop_output.txt', 'r') as apop_output:
            apop_data = apop_output.read().decode("utf-8")

        for rank, pocket_filename, apop_score, residues in re.findall(pocket_pattern, apop_data):
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

            pocket = []
            for resid, chain in re.findall(residue_pattern, residues):
                pocket.append([chain, [int(resid), " "]])
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
