"""Regression test for the Rosetta interface-score calculation."""

import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from allopep import prepare_rosetta_inputs, read_interface_scores, run_pepglad
from run_gaps import SpatialPatch


class InterfaceScoreTests(unittest.TestCase):
    def test_pepglad_inputs_and_output_paths(self):
        args = SimpleNamespace(
            max_rank=1, length_min=5, length_max=10, num_samples=1,
        )
        patch_score = SpatialPatch((('A', 1, ' '),), 0.8, 0.9)
        with tempfile.TemporaryDirectory(dir='.') as directory:
            out_dir = Path(directory) / 'output'
            output_pdb = (
                out_dir / 'target_PepGLAD_outputs'
                / 'target_allosteric_site_1_codesign' / 'target_0.pdb'
            )
            output_pdb.parent.mkdir(parents=True)
            output_pdb.touch()

            with patch('allopep.subprocess.run') as subprocess_run:
                generated = run_pepglad(
                    args, '/input/target.pdb', out_dir, 'target', [patch_score]
                )

            pocket_file = (
                out_dir / 'target_allosteric_sites'
                / 'target_allosteric_site_1.json'
            )
            self.assertEqual(pocket_file.read_text(), '[["A", [1, " "]]]')
            self.assertEqual(generated, [(output_pdb, 1, 0)])
            command = subprocess_run.call_args.args[0]
            self.assertEqual(command[command.index('--pocket') + 1], f'../{pocket_file}')
            self.assertEqual(command[command.index('--out_dir') + 1], f'../{output_pdb.parent}')

    def test_rosetta_input_lists_only_nonclashing_structures(self):
        fixture = Path(__file__).parent / 'fixtures' / 'sample_GAPS_output.pdb'
        with tempfile.TemporaryDirectory() as directory:
            score_dir = Path(directory) / 'scores'
            with patch('allopep.count_clashes.count_clashes', return_value=0):
                list_file = prepare_rosetta_inputs(
                    [(fixture, 1, 0)], score_dir, 'target', 'B'
                )
            listed = [Path(line) for line in list_file.read_text().splitlines()]
            self.assertEqual(len(listed), 3)
            self.assertTrue(all(path.is_file() for path in listed))
            self.assertEqual(
                [path.stem.rsplit('_', 1)[-1] for path in listed],
                ['complex', 'protein', 'peptide'],
            )

    def test_combines_matching_positive_structure_scores(self):
        with tempfile.TemporaryDirectory() as directory:
            score_file = Path(directory) / 'scores.tsv'
            score_file.write_text(
                'SEQUENCE:\n'
                'SCORE: total_score description\n'
                'SCORE: 10 target_allosteric_site_1_sample_0_complex_0001\n'
                'SCORE: 8 target_allosteric_site_1_sample_0_protein_0001\n'
                'SCORE: 5 target_allosteric_site_1_sample_0_peptide_0001\n'
                'SCORE: -1 target_allosteric_site_1_sample_1_complex_0001\n',
                encoding='utf-8',
            )
            scores = read_interface_scores(score_file)

        self.assertEqual(len(scores), 1)
        self.assertEqual(scores.loc[('target', '1', '0'), 'total_score'], 3)


if __name__ == '__main__':
    unittest.main()
