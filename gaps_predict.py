#!/usr/bin/env python3
"""Run the GAPS model inside the dedicated ``gaps`` environment."""

import argparse
import json
from pathlib import Path

import numpy as np
import torch as pt
from tqdm import tqdm

from GAPS.model.model import Model
from GAPS.src.config import config_model
from GAPS.src.data_encoding import en_features, en_structure, ext_topology
from GAPS.src.dataset import StructuresDataset, col_batch
from GAPS.src.structure import (concatenate_chains, encode_bfactor,
                                split_by_chain)
from GAPS.src.structure_io import save_pdb


def score_records(structure, probabilities):
    """Pair raw residue probabilities with PDB identities before B-factor export."""
    scores = np.asarray(probabilities).reshape(-1)
    internal_ids = np.unique(structure['resid'])
    if len(scores) != len(internal_ids):
        raise ValueError('GAPS score count differs from preprocessed residue count')

    records = []
    seen = set()
    for internal_id, score in zip(internal_ids, scores):
        indices = np.flatnonzero(structure['resid'] == internal_id)
        identities = {
            (str(structure['chain_name'][i]).split(':')[0],
             int(structure['original_resid'][i]),
             str(structure['original_icode'][i]) or ' ',
             str(structure['resname'][i]))
            for i in indices
        }
        if len(identities) != 1:
            raise ValueError(f'GAPS residue {internal_id} has ambiguous PDB identity')
        chain, number, icode, resname = identities.pop()
        key = (chain, number, icode)
        if key in seen:
            raise ValueError(f'Duplicate GAPS PDB residue identity: {key}')
        seen.add(key)
        ca_indices = [i for i in indices
                      if structure['name'][i] == 'CA'
                      and structure['element'][i] == 'C']
        if len(ca_indices) > 1:
            raise ValueError(f'Multiple C-alpha atoms for GAPS residue: {key}')
        ca = (structure['xyz'][ca_indices[0]].tolist()
              if ca_indices else None)
        records.append({
            'chain': chain, 'number': number, 'icode': icode,
            'resname': resname, 'score': float(score), 'ca': ca,
        })
    return records


def run_prediction(pdb_filepath, out_dir, pytorch_device='cuda'):
    device = pt.device(pytorch_device)
    model = Model(config_model)
    model.load_state_dict(pt.load(
        Path(__file__).with_name('GAPS') / 'checkpoints' / 'model.pt',
        map_location=device, weights_only=True,
    ))
    model = model.eval().to(device)

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    dataset = StructuresDataset([str(pdb_filepath)], with_preprocessing=True)
    with pt.no_grad():
        for subunits, filepath in tqdm(dataset):
            structure = concatenate_chains(subunits)
            X, M = en_structure(structure)
            q = en_features(structure)[0]
            ids_topk, _, _, _, _ = ext_topology(X, 64)
            X, ids_topk, q, M = col_batch([[X, ids_topk, q, M]])
            z = model(X.to(device), ids_topk.to(device), q.to(device), M.float().to(device))
            p = pt.sigmoid(z)
            raw_scores = p.cpu().numpy()
            scores_file = out_dir / (Path(filepath).stem + '_GAPS_scores.json')
            scores_file.write_text(json.dumps({
                'schema_version': 1,
                'residues': score_records(structure, raw_scores),
            }), encoding='utf-8')
            structure = encode_bfactor(structure, raw_scores)
            out_file = out_dir / (Path(filepath).stem + '_GAPS_output.pdb')
            save_pdb(split_by_chain(structure), str(out_file))


def main():
    parser = argparse.ArgumentParser(description='GAPS binding-site prediction')
    parser.add_argument('-p', '--pdb', required=True, type=Path)
    parser.add_argument('-o', '--out-dir', required=True, type=Path)
    parser.add_argument('-d', '--device', default='cuda')
    args = parser.parse_args()
    run_prediction(args.pdb, args.out_dir, args.device)


if __name__ == '__main__':
    main()
