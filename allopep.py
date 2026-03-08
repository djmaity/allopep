#!/usr/bin/env python3

import os
import re
import subprocess
import zipfile
import json
import argparse

parser = argparse.ArgumentParser(
    prog='python mlsim.py', description='Pipeline for allosteric peptide design')
parser.add_argument('pdb_file')
# parser.add_argument('--max_rank', default=None, type=int,
#     help='Maximum ranked allosteric pocket by APOP to use for peptide design using PepGLAD')
parser.add_argument('--length_min', default=5, type=int,
    help='Minimum length of the predicted peptide')
parser.add_argument('--length_max', default=10, type=int,
    help='Maximum length of the predicted peptide')
args = parser.parse_args()

prefix = os.path.splitext(os.path.basename(args.pdb_file))[0]
full_pdb_path = os.path.abspath(args.pdb_file)

apop_output_filename = os.path.abspath(f'output/{prefix}_apop_output.zip')

# Run APOP to get allosteric pockets
subprocess.run(['python', '../APOP/apop.py', full_pdb_path,
                '--output', apop_output_filename],
                cwd='output', check=True)

with zipfile.ZipFile(apop_output_filename) as apop_output_zip:
    with apop_output_zip.open('apop_output.txt', 'r') as apop_output:
        apop_data = apop_output.read().decode("utf-8")

pocket_pattern = re.compile(r"Rank=(\d+)\nPocket name: (\S+)\nAPOP score: (\S+)\nResidues: (.+)")
residue_pattern = re.compile(r"(\d+)([a-zA-Z]+)")

os.makedirs(f'output/{prefix}_allosteric_sites/', exist_ok=True)

pocket_json = []
for rank, name, score, residues in re.findall(pocket_pattern, apop_data):
    pocket = []
    for resid, chain in re.findall(residue_pattern, residues):
        pocket.append([chain, [int(resid), " "]])

    # Write the JSON input file for PepGLAD
    json.dump(pocket, open(f'output/{prefix}_allosteric_sites/{prefix}_allosteric_site_{rank}.json', 'w', encoding="utf-8"))

    # Run PepGLAD
    subprocess.run(['RAY_ACCEL_ENV_VAR_OVERRIDE_ON_ZERO=0',
        'python', '-m', 'api.run', '--mode', 'codesign',
        '--pdb', full_pdb_path,
        '--pocket', f'../output/{prefix}_allosteric_site_{rank}.json',
        '--out_dir', f'../output/{prefix}_allosteric_site_{rank}_codesign',
        '--length_min', str(args.length_min),
        '--length_max', str(args.length_max),
        '--n_samples', '10'
    ], cwd='PepGLAD', check=True)

# # Calculate Rosetta Score of the Predicted Peptides
# with open(f'output/{prefix}_output_pdb_files.txt', 'w') as rosetta_input:
#     for output_pdb_file in sorted(glob(
#             f'output/{prefix}_allosteric_site_*_codesign/{prefix}_*.pdb')):
#         print(output_pdb_file, file=rosetta_input)

# subprocess.run(['score_jd2', '-in:file:l', f'output/{prefix}_output_pdb_files.txt',
# '-out:file:scorefile', f'output/{prefix}_rosetta_scores.tsv'
# ])
