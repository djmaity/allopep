"""AlloPep: allosteric peptide prediction pipeline"""

#!/usr/bin/env python3

import argparse
import json
import os
import subprocess

# Local imports
import count_clashes
import pocket_mean_gaps_score
import run_gaps

import vina_score

parser = argparse.ArgumentParser(
    prog='python mlsim.py', description='Pipeline for allosteric peptide design')
parser.add_argument('pdb_file')
parser.add_argument('--max_rank', default=None, type=int,
    help='Maximum ranked allosteric pocket by GAPS to use for peptide design using PepGLAD')
parser.add_argument('--length_min', default=5, type=int,
    help='Minimum length of the predicted peptide')
parser.add_argument('--length_max', default=10, type=int,
    help='Maximum length of the predicted peptide')
parser.add_argument('--num_samples', default=10, type=int,
    help='Maximum length of the predicted peptide')
args = parser.parse_args()

prefix = os.path.splitext(os.path.basename(args.pdb_file))[0]
full_pdb_path = os.path.abspath(args.pdb_file)
OUT_DIR = f'output/{prefix}_output'
apop_out_file = os.path.abspath(f'{OUT_DIR}/{prefix}_apop_output.zip')

# Run GAPS to get peptide binding sites
run_gaps.run_gaps(full_pdb_path, out_dir=OUT_DIR)

# Run APOP to get allosteric pockets
subprocess.run(['python', '../APOP/apop.py', full_pdb_path,
                '--output', apop_out_file],
                cwd='output', check=True)

# Calculate mean GAPS score of the APOP pockets
pocket_scores = pocket_mean_gaps_score.pocket_mean_gaps_score(
    input_file=full_pdb_path,
    apop_file=apop_out_file,
    gaps_file=f'{OUT_DIR}/{prefix}_GAPS_output.pdb'
    )

# TODO: Check if allosteric pocket overlaps active site pocket using Jaccard Index

ALLOSTERIC_SITES_DIR = f'{OUT_DIR}/{prefix}_allosteric_sites/'
os.makedirs(ALLOSTERIC_SITES_DIR, exist_ok=True)

# Add required environment variables to supress warining from Ray
pepgald_env = os.environ.copy()
pepgald_env['RAY_ACCEL_ENV_VAR_OVERRIDE_ON_ZERO'] = '0'

MAX_RANK = len(pocket_scores)
if args.max_rank is not None and 0 < args.max_rank <= len(pocket_scores):
    MAX_RANK = args.max_rank

for rank, pocket in enumerate(pocket_scores[:MAX_RANK], start=1):
    # Write the JSON input file for PepGLAD
    JSON_FILENAME = f'{OUT_DIR}/{prefix}_allosteric_sites/{prefix}_allosteric_site_{rank}.json'
    with open(JSON_FILENAME, 'w', encoding="utf-8") as json_file:
        json.dump(pocket[1], json_file)

    # Run PepGLAD
    PEPGLAD_OUT_DIR = f'../{OUT_DIR}/{prefix}_PepGLAD_outputs/{prefix}_allosteric_site_{rank}_codesign'
    subprocess.run([
        'python', '-m', 'api.run', '--mode', 'codesign',
        '--pdb', full_pdb_path,
        '--pocket', f'../{ALLOSTERIC_SITES_DIR}/{prefix}_allosteric_site_{rank}.json',
        '--out_dir', PEPGLAD_OUT_DIR,
        '--length_min', str(args.length_min),
        '--length_max', str(args.length_max),
        '--n_samples', str(args.num_samples)
    ], cwd='PepGLAD', env=pepgald_env, check=True)

    # Calculate AutoDock Vina Score
    for i in range(args.num_samples):
        pepglad_prediction = f'{OUT_DIR}/{prefix}_PepGLAD_outputs/{prefix}_allosteric_site_{rank}_codesign/{prefix}_{i}.pdb'
        clashes = count_clashes.count_clashes(pdb_file=pepglad_prediction)
        if clashes > 0:
            print(f'{prefix}_{i}', 'Steric Clashes')
        else:
            peptide_chain = 'B'  # TODO: write code to automatically detect peptide chain
            score = vina_score.vina_score(
                pdb_file=pepglad_prediction, peptide_chain=peptide_chain)
            if score is not None:
                print(*score, 'kcal/mol')
            else:
                print([score])

# # Calculate Rosetta Score of the Predicted Peptides
# with open(f'output/{prefix}_output_pdb_files.txt', 'w') as rosetta_input:
#     for output_pdb_file in sorted(glob(
#             f'output/{prefix}_allosteric_site_*_codesign/{prefix}_*.pdb')):
#         print(output_pdb_file, file=rosetta_input)

# subprocess.run(['score_jd2', '-in:file:l', f'output/{prefix}_output_pdb_files.txt',
# '-out:file:scorefile', f'output/{prefix}_rosetta_scores.tsv'
# ])
