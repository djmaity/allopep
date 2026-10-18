"""Run with `conda run -n PepGLAD python tests/pepglad_guidance_case.py`."""

import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import torch
import ray

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'PepGLAD'))

from api.run import DesignDataset, peptide_coordinates_supported, relax_designs


class PepgladPocketTests(unittest.TestCase):
    def test_individual_relaxation_failure_is_quarantined_and_other_samples_finish(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            paths = [str(root / f'target_{index}.pdb') for index in range(2)]
            for path in paths:
                Path(path).write_text('generated PDB\n')
            (root / 'summary.jsonl').write_text(''.join(json.dumps({
                'id': f'target_{index}', 'pep_seq': 'LHLKA',
                'pep_chain': 'B', 'rec_chains': ['A'],
            }) + '\n' for index in range(2)))
            failure = ray.exceptions.RayTaskError('openmm_relax', 'test traceback',
                                                   ValueError('bad sample'))
            def get_result(future):
                if future == paths[0]:
                    raise failure
                return future
            with patch('api.run.ray.init'), \
                 patch('api.run.ray.shutdown') as shutdown, \
                 patch('api.run.openmm_relax.remote', side_effect=lambda path: path), \
                 patch('api.run.ray.wait', side_effect=lambda futures, **kw: (futures[:1], futures[1:])), \
                 patch('api.run.ray.get', side_effect=get_result) as get:
                relax_designs(paths, directory)
            self.assertEqual(get.call_count, 2)
            shutdown.assert_called_once()
            self.assertFalse(Path(paths[0]).exists())
            self.assertTrue((root / 'invalid_designs/target_0.pdb').is_file())
            self.assertTrue(Path(paths[1]).is_file())
            records = [json.loads(line) for line in (root / 'summary.jsonl').read_text().splitlines()]
            self.assertEqual(records[0]['status'], 'error')
            self.assertNotIn('status', records[1])
            self.assertTrue((root / 'target_0.pdb.relax_error.json').is_file())

    def test_divergent_coordinates_are_rejected_before_pdb_write(self):
        atom = SimpleNamespace(coordinate=torch.tensor([-30.0, 20.0, 10.0]))
        self.assertTrue(peptide_coordinates_supported([[atom]]))
        atom.coordinate = torch.tensor([5256507.0, 58924.664, -4087467.25])
        self.assertFalse(peptide_coordinates_supported([[atom]]))
        atom.coordinate = torch.tensor([float('nan'), 20.0, 10.0])
        self.assertFalse(peptide_coordinates_supported([[atom]]))
        atom.coordinate = torch.tensor([-999.9998, 20.0, 10.0])
        self.assertFalse(peptide_coordinates_supported([[atom]]))

    def test_two_residue_pocket_prepares_finite_covariance(self):
        pdb = ROOT / 'PepGLAD/assets/1ssc_A_B.pdb'
        entries = json.loads((ROOT / 'PepGLAD/assets/1ssc_A_pocket.json').read_text())[:2]
        with tempfile.TemporaryDirectory() as directory:
            pocket = Path(directory) / 'pocket.json'
            pocket.write_text(json.dumps(entries))
            dataset = DesignDataset([str(pdb)], [str(pocket)],
                                    lengths_range=[(8, 9)])
            sample = dataset[0]
            self.assertEqual((~sample['mask']).sum().item(), 2)
            self.assertEqual(tuple(sample['L'].shape), (1, 3, 3))
            self.assertTrue(torch.isfinite(sample['L']).all().item())
            self.assertNotIn('site_scores', sample)

    def test_two_field_pocket_with_insertion_code_and_no_scores(self):
        pdb = ROOT / 'PepGLAD/assets/1ssc_A_B.pdb'
        original = ROOT / 'PepGLAD/assets/1ssc_A_pocket.json'
        entries = json.loads(original.read_text())[:3]
        with tempfile.TemporaryDirectory() as directory:
            lines = pdb.read_text().splitlines()
            lines = [line[:26] + 'A' + line[27:]
                     if line.startswith('ATOM') and line[21] == 'A'
                     and line[22:26].strip() == '3' and line[26] == ' '
                     else line for line in lines]
            modified = Path(directory) / 'receptor.pdb'
            modified.write_text('\n'.join(lines) + '\n')
            entries[0][1][1] = 'A'
            pocket = Path(directory) / 'pocket.json'
            pocket.write_text(json.dumps(entries[::-1]))
            dataset = DesignDataset([str(modified)], [str(pocket)],
                                    lengths_range=[(8, 9)])
            residues, positions, _ = dataset.get_epitope(0)
            self.assertEqual([block.id for block in residues],
                             [(3, 'A'), (4, ' '), (5, ' ')])
            self.assertEqual(positions, [3, 4, 5])
            self.assertNotIn('site_scores', dataset[0])

            scored = Path(directory) / 'scored.json'
            scored.write_text(json.dumps([entries[0] + [0.9]]))
            with self.assertRaises(ValueError):
                DesignDataset([str(modified)], [str(scored)],
                              lengths_range=[(8, 9)]).get_epitope(0)


if __name__ == '__main__':
    unittest.main()
