#!/usr/bin/env python3
"""Run the GAPS model inside the dedicated ``gaps`` environment."""

import argparse
from pathlib import Path

import torch as pt
from tqdm import tqdm

from GAPS.model.model import Model
from GAPS.src.config import config_model
from GAPS.src.data_encoding import en_features, en_structure, ext_topology
from GAPS.src.dataset import StructuresDataset, col_batch
from GAPS.src.structure import (concatenate_chains, encode_bfactor,
                                split_by_chain)
from GAPS.src.structure_io import save_pdb


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
            structure = encode_bfactor(structure, p.cpu().numpy())
            out_file = out_dir / (Path(filepath).stem + '_GAPS_output.pdb')
            save_pdb(split_by_chain(structure), str(out_file))


def main():
    parser = argparse.ArgumentParser(description='GAPS binding-site prediction')
    parser.add_argument('--pdb', required=True, type=Path)
    parser.add_argument('--out-dir', required=True, type=Path)
    parser.add_argument('--device', default='cuda')
    args = parser.parse_args()
    run_prediction(args.pdb, args.out_dir, args.device)


if __name__ == '__main__':
    main()
