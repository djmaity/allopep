""" Run GAPS to predict peptide binding sites """

import os
from glob import glob

import torch as pt
from tqdm import tqdm

from GAPS.model.model import Model
from GAPS.src.config import config_model
from GAPS.src.data_encoding import en_features, en_structure, ext_topology
from GAPS.src.dataset import StructuresDataset, col_batch
from GAPS.src.structure import (concatenate_chains, encode_bfactor,
                                split_by_chain)
from GAPS.src.structure_io import save_pdb


def run_gaps(pdb_filepaths, out_dir='.', pytorch_device='cuda'):
    """ Get peptide binging site prediction using GAPS """

    model_weights = 'GAPS/checkpoints/model.pt'
    device = pt.device(pytorch_device)

    model = Model(config_model)
    model.load_state_dict(pt.load(model_weights,
                                  map_location=pt.device(pytorch_device),
                                  weights_only=True)
    )
    model = model.eval().to(device)

    if not isinstance(pdb_filepaths, list):
        pdb_filepaths = [pdb_filepaths]
    dataset = StructuresDataset(pdb_filepaths, with_preprocessing=True)

    os.makedirs(out_dir, exist_ok=True)
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

            filename = os.path.splitext(os.path.basename(filepath))[0]
            out_file = os.path.join(out_dir, filename + '_GAPS_output.pdb')
            save_pdb(split_by_chain(structure), out_file)


if __name__ == '__main__':
    # TODO: add argument parser
    run_gaps(['input/2VH7.pdb', 'input/5lvp_holo.pdb'], out_dir='output')
