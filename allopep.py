#!/usr/bin/env python3

import os
import re
import subprocess
import zipfile
import json
import argparse

from glob import glob

parser = argparse.ArgumentParser(
    prog='python mlsim.py', description='Pipeline for allosteric peptide design')
parser.add_argument('pdb_file')
parser.add_argument('--max_rank', default=3, type=int,
    help='Maximum ranked allosteric pocket by APOP to use for peptide design using PepGLAD'
)
args = parser.parse_args()

prefix = os.path.splitext(os.path.basename(args.pdb_file))[0]
full_pdb_path = os.path.abspath(args.pdb_file)

print(prefix)

# Run APOP to get allosteric pockets
subprocess.run(['python', 'apop.py', full_pdb_path,
                '--chain', 'A', '--output',
                os.path.abspath(f'output/{prefix}_apop_output.zip')], cwd='APOP')
# TODO: convert to subpackage and call module directly

with zipfile.ZipFile(f'output/{prefix}_apop_output.zip') as apop_output_zip:
    with apop_output_zip.open('apop_output.txt', 'r') as apop_output:
        apop_data = apop_output.read().decode("utf-8")

pocket_pattern = re.compile(r"Rank=(\d+)\nPocket name: (\S+)\nAPOP score: (\S+)\nResidues: (.+)")
residue_pattern = re.compile(r"(\d+)([a-zA-Z]+)")

pocket_json = []
for rank, name, score, residues in re.findall(pocket_pattern, apop_data):
    if int(rank) <= args.max_rank:
        pocket = []
        for resid, chain in re.findall(residue_pattern, residues):
            pocket.append([chain, [int(resid), " "]])
        json.dump(pocket, open(f'output/{prefix}_allosteric_site_{rank}.json', 'w'))

        # Run PepGLAD
        subprocess.run([
            'python', '-m', 'api.run', '--mode', 'codesign',
            '--pdb', full_pdb_path,
            '--pocket', f'../output/{prefix}_allosteric_site_{rank}.json',
            '--out_dir', f'../output/{prefix}_allosteric_site_{rank}_codesign',
            '--length_min', '8',
            '--length_max', '15',
            '--n_samples', '10'
        ], cwd='PepGLAD')
        # TODO: convert to subpackage and call module directly

# # Run PepGLAD
# # TODO: convert to subpackage and call module directly
# for pocket_file in glob(f'output/{prefix}_allosteric_site_*.json'):
#     pocket_name = os.path.splitext(os.path.basename(pocket_file))[0]
#     subprocess.run([
#         'python', '-m', 'api.run', '--mode', 'codesign',
#         '--pdb', full_pdb_path,
#         '--pocket', os.path.relpath(pocket_file, start='./PepGLAD'),
#         '--out_dir', f'../output/{pocket_name}_codesign',
#         '--length_min', '8',
#         '--length_max', '15',
#         '--n_samples', '10'
#     ], cwd='PepGLAD')


# Calculate Rosetta Score of the Predicted Peptides
with open(f'output/{prefix}_output_pdb_files.txt', 'w') as rosetta_input:
    for output_pdb_file in sorted(glob(
            f'output/{prefix}_allosteric_site_*_codesign/{prefix}_*.pdb')):
        print(output_pdb_file, file=rosetta_input)

subprocess.run(['score_jd2', '-in:file:l', f'output/{prefix}_output_pdb_files.txt',
'-out:file:scorefile', f'output/{prefix}_rosetta_scores.tsv'
])
