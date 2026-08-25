"""AlloPep: allosteric peptide prediction pipeline"""

#!/usr/bin/env python3

import argparse
import json
import os
import shutil
import subprocess
from pathlib import Path
import pandas as pd

from Bio.PDB import PDBIO, PDBParser, Select

# Local imports
import count_clashes
import run_gaps
import vina_score



def detect_peptide_chain(structure):
    shortest_chain_id = None
    shortest_chain_length = None

    for chain in structure[0]:
        chain_length = 0
        for residue in chain:
            if residue.id[0] == ' ':
                chain_length += 1

        if chain_length == 0:
            continue

        is_shorter = (
            shortest_chain_length is None
            or chain_length < shortest_chain_length
        )
        if is_shorter:
            shortest_chain_id = chain.id
            shortest_chain_length = chain_length

    if shortest_chain_id is None:
        raise ValueError('No peptide or protein chains found in the PDB file')

    return shortest_chain_id


def split_protein_peptide(pdb_file, peptide_chain=None):
    file_path, extension = os.path.splitext(pdb_file)
    basename = os.path.basename(file_path)
    protein_filename = f'{file_path}_protein'
    peptide_filename = f'{file_path}_peptide'

    pdb_parser = PDBParser(QUIET=True)
    structure = pdb_parser.get_structure(basename, pdb_file)

    if peptide_chain is None:
        peptide_chain = detect_peptide_chain(structure)

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


parser = argparse.ArgumentParser(
    prog='python allopep.py', description='Pipeline for allosteric peptide design')
parser.add_argument('pdb_file')
parser.add_argument('--max_rank', default=None, type=int,
    help='Maximum ranked allosteric pocket by GAPS to use for peptide design using PepGLAD')
parser.add_argument('--length_min', default=5, type=int,
    help='Minimum length of the predicted peptide')
parser.add_argument('--length_max', default=10, type=int,
    help='Maximum length of the predicted peptide')
parser.add_argument('--num_samples', default=10, type=int,
    help='Maximum length of the predicted peptide')
parser.add_argument('--peptide-chain', default=None,
    help='Known peptide chain ID; automatically detected if omitted')
args = parser.parse_args()

prefix = os.path.splitext(os.path.basename(args.pdb_file))[0]
full_pdb_path = os.path.abspath(args.pdb_file)
OUT_DIR = f'output/{prefix}_output'

# Run GAPS to get peptide binding sites
run_gaps.run_gaps(full_pdb_path, out_dir=OUT_DIR)

# Find spatial patches with high GAPS scores
spatial_patches = run_gaps.find_high_bfactor_spatial_patches(
    pdb_file=f'{OUT_DIR}/{prefix}_GAPS_output.pdb'
)

# TODO: Check if allosteric pocket overlaps active site pocket using Jaccard Index

ALLOSTERIC_SITES_DIR = f'{OUT_DIR}/{prefix}_allosteric_sites/'
os.makedirs(ALLOSTERIC_SITES_DIR, exist_ok=True)

# Add required environment variables to supress warining from Ray
pepgald_env = os.environ.copy()
pepgald_env['RAY_ACCEL_ENV_VAR_OVERRIDE_ON_ZERO'] = '0'

MAX_RANK = len(spatial_patches)
if args.max_rank is not None and 0 < args.max_rank <= len(spatial_patches):
    MAX_RANK = args.max_rank


df = pd.DataFrame([
    {
        'pocket': patch.as_pepglad_pocket(),
        'mean_gaps_score': patch.mean_bfactor,
        'max_gaps_score': patch.max_bfactor,
    }
    for patch in spatial_patches
])
print(df)
df.to_csv(f'{OUT_DIR}/pocket_scores.csv', index=False)


generated_peptides = []
for rank, patch in enumerate(spatial_patches[:MAX_RANK], start=1):
    # Write the JSON input file for PepGLAD
    JSON_FILENAME = f'{OUT_DIR}/{prefix}_allosteric_sites/{prefix}_allosteric_site_{rank}.json'
    with open(JSON_FILENAME, 'w', encoding="utf-8") as json_file:
        json.dump(patch.as_pepglad_pocket(), json_file)

    # Run PepGLAD
    codesign_output_dir = (
        Path(OUT_DIR)
        / f'{prefix}_PepGLAD_outputs'
        / f'{prefix}_allosteric_site_{rank}_codesign'
    )
    PEPGLAD_OUT_DIR = f'../{codesign_output_dir}'
    subprocess.run([
        'python', '-m', 'api.run', '--mode', 'codesign',
        '--pdb', full_pdb_path,
        '--pocket', f'../{ALLOSTERIC_SITES_DIR}/{prefix}_allosteric_site_{rank}.json',
        '--out_dir', PEPGLAD_OUT_DIR,
        '--length_min', str(args.length_min),
        '--length_max', str(args.length_max),
        '--n_samples', str(args.num_samples)
    ], cwd='PepGLAD', env=pepgald_env, check=True)

    for sample_index in range(args.num_samples):
        output_pdb_path = codesign_output_dir / f'{prefix}_{sample_index}.pdb'
        if not output_pdb_path.is_file():
            raise FileNotFoundError(
                f'PepGLAD did not produce expected output: {output_pdb_path}'
            )
        generated_peptides.append(
            (output_pdb_path, prefix, rank, sample_index)
        )

    # # Calculate AutoDock Vina Score
    # for i in range(args.num_samples):
    #     pepglad_prediction = f'{OUT_DIR}/{prefix}_PepGLAD_outputs/{prefix}_allosteric_site_{rank}_codesign/{prefix}_{i}.pdb'
    #     clashes = count_clashes.count_clashes(pdb_file=pepglad_prediction)
    #     if clashes > 0:
    #         print(f'{prefix}_{i}', 'Steric Clashes')
    #     else:
    #         peptide_chain = 'B'  # TODO: write code to automatically detect peptide chain
    #         score = vina_score.vina_score(
    #             pdb_file=pepglad_prediction, peptide_chain=peptide_chain)
    #         if score is not None:
    #             print(*score, 'kcal/mol')
    #         else:
    #             print([score])

ROSETTA_SCORE_DIR = Path(f'{OUT_DIR}/{prefix}_rosetta_score')
ROSETTA_SCORE_DIR.mkdir(parents=True, exist_ok=True)
tmp_array = []
for output_pdb_path, run_prefix, rank, sample_index in generated_peptides:
    rosetta_pdb_filename = (
        ROSETTA_SCORE_DIR
        / f'{run_prefix}_allosteric_site_{rank}_sample_{sample_index}.pdb'
    )
    shutil.copy(output_pdb_path, rosetta_pdb_filename)
    tmp_array.append((run_prefix, rank, sample_index))

# Sort the array by prefix, rank, and index
tmp_array.sort(key=lambda x: (x[0], x[1], x[2]))

# Calculate Rosetta Score of the Predicted Peptides
rosetta_input = f'{ROSETTA_SCORE_DIR}/{prefix}_pdb_file_list.txt'
rosetta_input_file = open(rosetta_input, 'w')

for run_prefix, rank, i in tmp_array:
    rosetta_pdb_filename = f'{ROSETTA_SCORE_DIR}/{run_prefix}_allosteric_site_{rank}_sample_{i}.pdb'
    clashes = count_clashes.count_clashes(pdb_file=rosetta_pdb_filename)
    if clashes > 0:
        print(f'{run_prefix}_{i}', 'Steric Clashes')
    else:
        split_protein_peptide(
            pdb_file=rosetta_pdb_filename,
            peptide_chain=args.peptide_chain)

        rosetta_pdb_path, _ = os.path.splitext(rosetta_pdb_filename)
        os.rename(rosetta_pdb_filename, rosetta_pdb_path + '_complex.pdb')

        rosetta_input_file.write(rosetta_pdb_path + '_complex.pdb\n')
        rosetta_input_file.write(rosetta_pdb_path + '_protein.pdb\n')
        rosetta_input_file.write(rosetta_pdb_path + '_peptide.pdb\n')


rosetta_input_file.close()
rosetta_output = f'{ROSETTA_SCORE_DIR}/{prefix}_rosetta_scores.tsv'

# Remove the old Rosetta Score output file to prevent appending to it
if os.path.isfile(rosetta_output):
    os.remove(rosetta_output)

subprocess.run(['score_jd2', '-in:file:l', rosetta_input, '-score:weights', 'rosetta_pdbbind_interface_regression.wts', '-out:file:scorefile', rosetta_output])

# Calculate Interface Rosetta Score
df_rosetta_score = pd.read_table(rosetta_output, skiprows=1, sep=r'\s+', usecols=['description', 'total_score'])
df_rosetta_score['description'] = df_rosetta_score['description'].str.slice(stop=-5)
df_rosetta_score[['prefix', 'site', 'peptide', 'structure']] = df_rosetta_score['description'].str.extract(r'(\S+)_allosteric_site_(\d+)_sample_(\d+)_([a-zA-Z]+)')
df_rosetta_score.drop(columns='description', inplace=True)
df_rosetta_score = df_rosetta_score[df_rosetta_score['total_score'] > 0]

df_complex = df_rosetta_score[df_rosetta_score['structure'] == 'complex'].copy()
df_complex = df_complex.set_index(['prefix', 'site', 'peptide'], drop=True)
df_complex.drop(columns='structure', inplace=True)

df_protein = df_rosetta_score[df_rosetta_score['structure'] == 'protein'].copy()
df_protein = df_protein.set_index(['prefix', 'site', 'peptide'], drop=True)
df_protein.drop(columns='structure', inplace=True)

df_peptide = df_rosetta_score[df_rosetta_score['structure'] == 'peptide'].copy()
df_peptide = df_peptide.set_index(['prefix', 'site', 'peptide'], drop=True)
df_peptide.drop(columns='structure', inplace=True)

df_interface = df_protein + df_peptide - df_complex
df_interface.dropna(inplace=True)

df_interface.to_csv(f'{ROSETTA_SCORE_DIR}/{prefix}_rosetta_interface_scores.csv')
print(df_interface)
