"""Tests for GAPS score provenance and connected receptor pockets."""

import json
import math
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from allopep.gaps import (ScoredResidue, load_scored_residues, run_gaps,
                      select_contexts, select_guidance_contexts)


class GapsGuidanceTests(unittest.TestCase):
    def test_prediction_uses_gaps_environment(self):
        pdb = Path(__file__).parent / 'fixtures' / 'sample_GAPS_output.pdb'
        with patch('allopep.gaps.subprocess.run') as subprocess_run:
            run_gaps(pdb, out_dir='output/test', pytorch_device='cpu')
        command = subprocess_run.call_args.args[0]
        self.assertEqual(command[:6], [
            'conda', 'run', '--no-capture-output', '-n', 'gaps', 'python',
        ])
        self.assertEqual(command[command.index('--device') + 1], 'cpu')

    def test_original_number_and_insertion_code_are_checked(self):
        fixture = Path(__file__).parent / 'fixtures' / 'sample_GAPS_output.pdb'
        lines = fixture.read_text().splitlines()
        positions = [('A', 42, ' '), ('A', 42, 'A'), ('A', 43, ' ')]
        for index, (chain, number, icode) in enumerate(positions):
            line = lines[index]
            lines[index] = line[:21] + chain + f'{number:4d}' + icode + line[27:]
        records = []
        for index, (chain, number, icode) in enumerate(positions):
            records.append({
                'chain': chain, 'number': number, 'icode': icode,
                'resname': 'ALA', 'score': [0.873241, 0.772346, 0.123456][index],
                'ca': [[0.0, 0.0, 0.0], [4.0, 0.0, 0.0], [8.0, 0.0, 0.0]][index],
            })
        areas = {'A': {
            '42': SimpleNamespace(total=12.5, relativeTotal=0.0),
            '42A': SimpleNamespace(total=0.25, relativeTotal=float('nan')),
            '43': SimpleNamespace(total=0.0, relativeTotal=0.9),
        }}
        mock_sasa = SimpleNamespace(residueAreas=lambda: areas)
        with tempfile.TemporaryDirectory() as directory:
            pdb = Path(directory) / 'target.pdb'
            scores = Path(directory) / 'target_GAPS_scores.json'
            pdb.write_text('\n'.join(lines) + '\n')
            scores.write_text(json.dumps({'schema_version': 1, 'residues': records}))
            with patch('allopep.gaps.freesasa.calcBioPDB', return_value=(mock_sasa, {})):
                scored = load_scored_residues(pdb, scores)
                self.assertEqual([r.key for r in scored], positions)
                self.assertEqual(scored[0].score, 0.873241)
                self.assertEqual([r.sasa for r in scored], [12.5, 0.25, 0.0])
                pocket, = select_contexts(scored)
                self.assertEqual([r.key for r in pocket.residues], positions[:2])
                records[1]['icode'] = 'B'
                scores.write_text(json.dumps({'schema_version': 1, 'residues': records}))
                with self.assertRaisesRegex(ValueError, 'missing from input PDB'):
                    load_scored_residues(pdb, scores)
                records[1]['icode'] = ' '
                scores.write_text(json.dumps({'schema_version': 1, 'residues': records}))
                with self.assertRaisesRegex(ValueError, 'Duplicate GAPS residue'):
                    load_scored_residues(pdb, scores)

    def test_score_cutoff_and_positive_absolute_sasa_are_strict(self):
        base = [
            ScoredResidue('A', 1, ' ', 0.8, (0.0, 0.0, 0.0), 0.5),
            ScoredResidue('A', 2, ' ', 0.7, (1.0, 0.0, 0.0), 0.5),
        ]
        candidate = ScoredResidue('A', 3, ' ', 0.9, (2.0, 0.0, 0.0), 0.5)
        for score, sasa, included in (
                (0.0, 0.5, False), (0.499999, 0.5, False),
                (0.5, 0.5, False), (0.500001, 0.5, True),
                (0.9, -0.1, False), (0.9, 0.0, False),
                (0.9, 0.000001, True), (0.9, 0.05, True),
                (0.9, 0.1, True), (0.9, 0.25, True)):
            with self.subTest(score=score, sasa=sasa):
                residue = replace(candidate, score=score, sasa=sasa)
                pocket, = select_contexts([*base, residue])
                self.assertEqual(pocket.residues,
                                 tuple([*base, residue] if included else base))

    def test_connection_distance_boundary_and_pocket_separation(self):
        for distance, expected_sizes in (
                (7.999999, [4]), (8.0, [4]), (8.000001, [2, 2])):
            with self.subTest(distance=distance):
                residues = [
                    ScoredResidue('A', 1, ' ', 0.8, (0.0, 0.0, -1.0), 0.5),
                    ScoredResidue('A', 2, ' ', 0.8, (0.0, 0.0, 0.0), 0.5),
                    ScoredResidue('A', 3, ' ', 0.8, (0.0, 0.0, distance), 0.5),
                    ScoredResidue('A', 4, ' ', 0.8, (0.0, 0.0, distance + 1), 0.5),
                ]
                pockets = select_contexts(residues)
                self.assertEqual([len(p.residues) for p in pockets], expected_sizes)
                keys = [r.key for p in pockets for r in p.residues]
                self.assertEqual(len(keys), len(set(keys)))
                if len(pockets) == 2:
                    self.assertTrue(all(
                        math.dist(a.ca, b.ca) > 8.0
                        for a in pockets[0].residues for b in pockets[1].residues))

    def test_connections_are_transitive_without_a_center_radius(self):
        residues = [
            ScoredResidue('A', i + 1, ' ', 0.9 if i == 0 else 0.8,
                          (float(i * 7), 0.0, 0.0), 0.5)
            for i in range(6)
        ]
        pocket, = select_contexts(residues)
        self.assertEqual(pocket.residues, tuple(residues))
        self.assertEqual(pocket.seed, residues[0])
        self.assertGreater(math.dist(pocket.seed.ca, residues[-1].ca), 16.0)

    def test_ineligible_residues_cannot_bridge_pockets(self):
        residues = [
            ScoredResidue('A', 1, ' ', 0.8, (-1.0, 0.0, 0.0), 0.5),
            ScoredResidue('A', 2, ' ', 0.8, (0.0, 0.0, 0.0), 0.5),
            ScoredResidue('A', 4, ' ', 0.8, (16.0, 0.0, 0.0), 0.5),
            ScoredResidue('A', 5, ' ', 0.8, (17.0, 0.0, 0.0), 0.5),
        ]
        for score, sasa in ((0.5, 0.5), (0.0, 0.5), (0.9, -0.1), (0.9, 0.0)):
            with self.subTest(score=score, sasa=sasa):
                bridge = ScoredResidue('A', 3, ' ', score, (8.0, 0.0, 0.0), sasa)
                pockets = select_contexts([*residues, bridge])
                self.assertEqual([p.residues for p in pockets],
                                 [tuple(residues[:2]), tuple(residues[2:])])
        bridge = ScoredResidue('A', 3, ' ', 0.500001, (8.0, 0.0, 0.0), 0.000001)
        pocket, = select_contexts([*residues, bridge])
        self.assertEqual([r.number for r in pocket.residues], [1, 2, 3, 4, 5])

    def test_singletons_are_discarded_and_pairs_are_kept(self):
        residues = [
            ScoredResidue('A', 1, ' ', 1.0, (0.0, 0.0, 0.0), 0.5),
            ScoredResidue('A', 2, ' ', 0.8, (20.0, 0.0, 0.0), 0.5),
            ScoredResidue('A', 3, ' ', 0.9, (24.0, 0.0, 0.0), 0.5),
        ]
        pocket, = select_contexts(residues)
        self.assertEqual(pocket.residues, tuple(residues[1:]))
        self.assertEqual(pocket.seed, residues[2])

    def test_connected_pockets_have_no_size_cap(self):
        residues = [
            ScoredResidue('A', i + 1, ' ', 0.8, (float(i * 4), 0.0, 0.0), 0.5)
            for i in range(60)
        ]
        pocket, = select_contexts(residues)
        self.assertEqual(pocket.residues, tuple(residues))

    def test_pockets_rank_by_mean_instead_of_maximum_or_total_score(self):
        residues = [
            ScoredResidue('A', 1, ' ', 0.99, (0.0, 0.0, 0.0), 0.5),
            ScoredResidue('A', 2, ' ', 0.51, (1.0, 0.0, 0.0), 0.5),
            ScoredResidue('A', 3, ' ', 0.8, (20.0, 0.0, 0.0), 0.5),
            ScoredResidue('A', 4, ' ', 0.8, (21.0, 0.0, 0.0), 0.5),
            ScoredResidue('A', 5, ' ', 0.6, (40.0, 0.0, 0.0), 0.5),
            ScoredResidue('A', 6, ' ', 0.6, (41.0, 0.0, 0.0), 0.5),
            ScoredResidue('A', 7, ' ', 0.6, (42.0, 0.0, 0.0), 0.5),
        ]
        pockets = select_contexts(residues)
        self.assertEqual([p.residues for p in pockets], [
            tuple(residues[2:4]), tuple(residues[:2]), tuple(residues[4:]),
        ])

    def test_default_selects_all_and_rank_limits_apply_after_ranking(self):
        residues = [
            ScoredResidue('A', site * 10 + offset + 1, ' ',
                          0.9 - site * 0.1 - offset * 0.05,
                          (site * 20.0 + offset * 2.0, 0.0, 0.0), 0.5)
            for site in range(4) for offset in range(2)
        ]
        with patch('allopep.gaps.load_scored_residues', return_value=residues):
            default = select_guidance_contexts('unused', 'unused')
            explicit_all = select_guidance_contexts('unused', 'unused', max_rank='all')
            top = select_guidance_contexts('unused', 'unused', max_rank=2)
        self.assertEqual(len(default), 4)
        self.assertEqual(default, explicit_all)
        self.assertEqual(default, select_contexts(residues))
        self.assertEqual(top, default[:2])
        self.assertEqual(select_contexts(residues, max_rank=1), default[:1])
        self.assertEqual(select_contexts(residues, max_rank=10), default)
        self.assertEqual(default[0].as_pepglad_pocket(),
                         [['A', [1, ' ']], ['A', [2, ' ']]])

    def test_identity_ties_cross_chain_pockets_and_input_order(self):
        residues = [
            ScoredResidue('B', 1, ' ', 0.875, (6.0, 0.0, 0.0), 0.5),
            ScoredResidue('A', 42, 'A', 0.75, (4.0, 0.0, 0.0), 0.5),
            ScoredResidue('A', 42, ' ', 0.625, (0.0, 0.0, 0.0), 0.5),
            ScoredResidue('A', 100, ' ', 0.75, (30.0, 0.0, 0.0), 0.5),
            ScoredResidue('B', 2, ' ', 0.75, (34.0, 0.0, 0.0), 0.5),
        ]
        expected = select_contexts(residues)
        self.assertEqual([p.as_pepglad_pocket() for p in expected], [
            [['A', [42, ' ']], ['A', [42, 'A']], ['B', [1, ' ']]],
            [['A', [100, ' ']], ['B', [2, ' ']]],
        ])
        self.assertEqual([p.seed.key for p in expected],
                         [('B', 1, ' '), ('A', 100, ' ')])
        for inputs in (list(reversed(residues)), residues[2:] + residues[:2]):
            self.assertEqual(select_contexts(inputs), expected)

    def test_no_usable_pockets_preserves_selection_error(self):
        for residues in (
                [],
                [ScoredResidue('A', 1, ' ', 0.9, (0.0, 0.0, 0.0), 0.5)],
                [ScoredResidue('A', 1, ' ', 0.5, (0.0, 0.0, 0.0), 0.5)],
                [ScoredResidue('A', i + 1, ' ', 0.9, (float(i * 20), 0.0, 0.0), 0.5)
                 for i in range(2)]):
            with self.subTest(residues=residues), self.assertRaisesRegex(
                    ValueError, 'No usable GAPS-guided receptor contexts found'):
                select_contexts(residues)

    def test_invalid_rank_is_rejected(self):
        for rank in (0, -1, True, 1.0, '2', None, float('nan')):
            with self.subTest(rank=rank), self.assertRaisesRegex(ValueError, 'max_rank'):
                select_contexts([], max_rank=rank)

    def test_invalid_sidecar_schema_is_rejected_before_pdb_loading(self):
        with tempfile.TemporaryDirectory() as directory:
            scores = Path(directory) / 'scores.json'
            for payload in ([], {'schema_version': 2, 'residues': []},
                            {'schema_version': 1, 'residues': {}}):
                scores.write_text(json.dumps(payload))
                with self.assertRaisesRegex(ValueError, 'Invalid GAPS score sidecar'):
                    load_scored_residues('unused.pdb', scores)


if __name__ == '__main__':
    unittest.main()
