"""Run with `conda run -n gaps python tests/gaps_scores_case.py`."""

import json
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from GAPS.src.structure import cle_structure, encode_bfactor, split_by_chain
from GAPS.src.structure_io import save_pdb
from gaps_predict import score_records


class RawGapsScoresTests(unittest.TestCase):
    def test_original_ids_and_full_precision_survive_preprocessing(self):
        source = {
            'name': np.array(['CA', 'N', 'CA', 'N', 'CA', 'N']),
            'element': np.array(['C', 'N', 'C', 'N', 'C', 'N']),
            'het_flag': np.array(['A'] * 6),
            'resname': np.array(['ALA'] * 6),
            'resid': np.array([42, 42, 42, 42, 7, 7]),
            'icode': np.array(['', '', 'A', 'A', '', '']),
            'chain_name': np.array(['A:0'] * 4 + ['B:0'] * 2),
            'xyz': np.array([
                [0, 0, 0], [1, 0, 0], [4, 0, 0], [5, 0, 0],
                [20, 0, 0], [21, 0, 0],
            ], dtype=np.float32),
        }
        structure = cle_structure(source)
        probabilities = np.array([0.873241, 0.712345, 0.123456], dtype=np.float32)
        records = score_records(structure, probabilities)
        self.assertEqual([(r['chain'], r['number'], r['icode']) for r in records], [
            ('A', 42, ' '), ('A', 42, 'A'), ('B', 7, ' '),
        ])
        serialized = json.loads(json.dumps(records))
        self.assertGreater(abs(serialized[0]['score'] - 0.87), 0.003)
        encoded = encode_bfactor(structure, probabilities)
        self.assertAlmostEqual(float(encoded['bfactor'][0]), serialized[0]['score'])
        with tempfile.TemporaryDirectory() as directory:
            pdb = Path(directory) / 'scored.pdb'
            save_pdb(split_by_chain(encoded), str(pdb))
            first_atom = pdb.read_text().splitlines()[0]
            self.assertEqual(float(first_atom[60:66]), 0.87)


if __name__ == '__main__':
    unittest.main()
