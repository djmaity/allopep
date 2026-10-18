"""Tests for GAPS score provenance and local receptor contexts."""

import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from run_gaps import (GuidanceSelection, ScoredResidue, load_scored_residues,
                      run_gaps, select_contexts, select_guidance_contexts)


class GapsGuidanceTests(unittest.TestCase):
    def test_prediction_uses_gaps_environment(self):
        pdb = Path(__file__).parent / 'fixtures' / 'sample_GAPS_output.pdb'
        with patch('run_gaps.subprocess.run') as subprocess_run:
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
            '42': SimpleNamespace(relativeTotal=0.5),
            '42A': SimpleNamespace(relativeTotal=0.5),
            '43': SimpleNamespace(relativeTotal=0.5),
        }}
        mock_sasa = SimpleNamespace(residueAreas=lambda: areas)
        with tempfile.TemporaryDirectory() as directory:
            pdb = Path(directory) / 'target.pdb'
            scores = Path(directory) / 'target_GAPS_scores.json'
            pdb.write_text('\n'.join(lines) + '\n')
            scores.write_text(json.dumps({'schema_version': 1, 'residues': records}))
            with patch('run_gaps.freesasa.calcBioPDB', return_value=(mock_sasa, {})):
                scored = load_scored_residues(pdb, scores)
                self.assertEqual([r.key for r in scored], positions)
                self.assertEqual(scored[0].score, 0.873241)
                records[1]['icode'] = 'B'
                scores.write_text(json.dumps({'schema_version': 1, 'residues': records}))
                with self.assertRaisesRegex(ValueError, 'missing from input PDB'):
                    load_scored_residues(pdb, scores)
                records[1]['icode'] = ' '
                scores.write_text(json.dumps({'schema_version': 1, 'residues': records}))
                with self.assertRaisesRegex(ValueError, 'Duplicate GAPS residue'):
                    load_scored_residues(pdb, scores)

    def test_separated_contexts_keep_continuous_scores(self):
        residues = []
        for site in range(4):
            for offset in range(3):
                residues.append(ScoredResidue(
                    'A', site * 10 + offset + 1, ' ',
                    0.9 - site * 0.1 if offset == 0 else 0.123456,
                    (site * 20.0 + offset * 2.0, 0.0, 0.0), 0.5,
                ))
        with patch('run_gaps.load_scored_residues', return_value=residues):
            top = select_guidance_contexts('unused', 'unused')
            all_sites = select_guidance_contexts('unused', 'unused', max_rank='all')
        self.assertEqual(len(top), 3)
        self.assertEqual(len(all_sites), 4)
        self.assertEqual(len(top[0].residues), 3)
        self.assertEqual(top[0].as_pepglad_pocket()[1], ['A', [2, ' ']])
        self.assertEqual(top[0].residues[1].score, 0.123456)
        self.assertEqual(top, select_contexts(residues, GuidanceSelection()))

    def test_invalid_rank_is_rejected(self):
        with self.assertRaisesRegex(ValueError, 'max_rank'):
            select_guidance_contexts('unused', 'unused', max_rank=0)
        for value in (float('nan'), float('inf'), -1):
            with self.subTest(value=value), \
                 self.assertRaisesRegex(ValueError, 'spatial context'):
                GuidanceSelection(context_radius=value)

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
