"""Run with `conda run -n PepGLAD python tests/pepglad_guidance_case.py`."""

import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'PepGLAD'))

from api.run import DesignDataset, combined_guidance_energy, score_proximity_energy
from models.LDM.diffusion.transition import ContinuousTransition


class PepgladGuidanceTests(unittest.TestCase):
    def test_scored_and_unscored_pocket_formats(self):
        pdb = ROOT / 'PepGLAD/assets/1ssc_A_B.pdb'
        original = ROOT / 'PepGLAD/assets/1ssc_A_pocket.json'
        entries = json.loads(original.read_text())[:3]
        with tempfile.TemporaryDirectory() as directory:
            scored = Path(directory) / 'scored.json'
            scored.write_text(json.dumps([
                [entries[2][0], entries[2][1], 0.3],
                [entries[0][0], entries[0][1], 0.9],
                [entries[1][0], entries[1][1], 0.6],
            ]))
            dataset = DesignDataset([str(pdb)], [str(scored)], lengths_range=[(8, 9)])
            item = dataset[0]
            for actual, expected in zip(item['site_scores'][:3].tolist(),
                                        [0.9, 0.6, 0.3]):
                self.assertAlmostEqual(actual, expected)
            self.assertTrue(torch.all(item['site_scores'][3:] == 0))

            legacy = DesignDataset([str(pdb)], [str(original)], lengths_range=[(8, 9)])
            self.assertTrue(torch.all(legacy[0]['site_scores'] == 0))

    def test_score_gradient_and_zero_strength(self):
        X = torch.zeros((3, 14, 3), dtype=torch.double, requires_grad=True)
        with torch.no_grad():
            X[0, 1, 0] = 0.0
            X[1, 1, 0] = 10.0
            X[2, 0, 0] = 6.0
        mask = torch.tensor([False, False, True])
        batch = torch.zeros(3, dtype=torch.long)
        scores = torch.tensor([0.9, 0.1, 0.0])
        energy = score_proximity_energy(X, mask, batch, scores)
        energy.backward()
        self.assertTrue(torch.isfinite(X.grad).all())
        self.assertGreater(X.grad[2, 0, 0].item(), 0)

        model = SimpleNamespace(latent_geometry_guidance=lambda X, m, b: X[m, 0, 0].sum())
        baseline = model.latent_geometry_guidance(X, mask, batch)
        zero = combined_guidance_energy(model, X, mask, batch, scores, 0.0)
        self.assertEqual(zero.item(), baseline.item())

    def test_guidance_weight_is_applied(self):
        transition = ContinuousTransition(num_steps=3)
        X = torch.zeros((1, 14, 3))
        mask = torch.tensor([True])
        batch = torch.tensor([0])
        t = torch.tensor([2])
        guidance = torch.ones_like(X)
        torch.manual_seed(3)
        baseline = transition.denoise(X, X, mask, batch, t)
        torch.manual_seed(3)
        zero = transition.denoise(X, X, mask, batch, t,
                                  guidance=guidance, guidance_weight=0.0)
        torch.manual_seed(3)
        weighted = transition.denoise(X, X, mask, batch, t,
                                      guidance=guidance, guidance_weight=1.0)
        self.assertTrue(torch.equal(baseline, zero))
        self.assertFalse(torch.equal(baseline, weighted))


if __name__ == '__main__':
    unittest.main()
