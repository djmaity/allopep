"""Regression test for the Rosetta interface-score calculation."""

import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from allopep import (find_existing_pepglad_outputs, prepare_input_pdb,
                     prepare_rosetta_inputs, read_interface_scores,
                     run_openstructure_validation, run_pepglad)
from run_gaps import SpatialPatch


class InterfaceScoreTests(unittest.TestCase):
    def test_validation_reuses_saved_input_pdb(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / 'target.pdb'
            source.write_bytes(b'original input')
            out_dir = root / 'output' / 'target_output'

            saved = Path(prepare_input_pdb(source, out_dir))
            self.assertEqual(saved, (out_dir / 'target.pdb').resolve())
            self.assertEqual(saved.read_bytes(), source.read_bytes())

            source.write_bytes(b'changed input')
            self.assertEqual(
                prepare_input_pdb(source, out_dir, validate_existing=True),
                str(saved),
            )
            self.assertEqual(saved.read_bytes(), b'original input')
            with self.assertRaisesRegex(ValueError, 'different input PDB'):
                prepare_input_pdb(source, out_dir)
            self.assertEqual(saved.read_bytes(), b'original input')

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
            self.assertEqual(command[:7], [
                'conda', 'run', '--no-capture-output', '-n', 'PepGLAD',
                'python', '-m',
            ])
            self.assertEqual(command[command.index('--pocket') + 1], f'../{pocket_file}')
            self.assertEqual(command[command.index('--out_dir') + 1], f'../{output_pdb.parent}')

    def test_rosetta_input_lists_only_validated_structures(self):
        fixture = Path(__file__).parent / 'fixtures' / 'sample_GAPS_output.pdb'
        with tempfile.TemporaryDirectory() as directory:
            score_dir = Path(directory) / 'scores'
            with patch('allopep.run_openstructure_validation',
                       return_value={str(fixture.resolve())}):
                list_file = prepare_rosetta_inputs(
                    [(fixture, 1, 0)], score_dir, 'target', 'B', fixture
                )
            listed = [Path(line) for line in list_file.read_text().splitlines()]
            self.assertEqual(len(listed), 3)
            self.assertTrue(all(path.is_file() for path in listed))
            self.assertEqual(
                [path.stem.rsplit('_', 1)[-1] for path in listed],
                ['complex', 'protein', 'peptide'],
            )

    def test_invalid_structure_is_excluded_from_rosetta(self):
        fixture = Path(__file__).parent / 'fixtures' / 'sample_GAPS_output.pdb'
        with tempfile.TemporaryDirectory() as directory:
            score_dir = Path(directory) / 'scores'
            with patch('allopep.run_openstructure_validation', return_value=set()):
                list_file = prepare_rosetta_inputs(
                    [(fixture, 1, 0)], score_dir, 'target', 'B', fixture
                )
            self.assertEqual(list_file.read_text(), '')
            self.assertFalse(list(score_dir.glob('*.pdb')))

    def test_validation_manifest_uses_pepglad_summary_and_input(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            site = root / 'codesign'
            site.mkdir()
            output = site / 'target_0.pdb'
            output.touch()
            (site / 'summary.jsonl').write_text(
                '{"id":"target_0","pep_seq":"LHLKA","pep_chain":"C",'
                '"rec_chains":["A","B"]}\n', encoding='utf-8')
            input_pdb = root / 'input.pdb'
            input_pdb.touch()
            score_dir = root / 'scores'
            score_dir.mkdir()
            def write_result(command, check):
                csv_path = Path(command[command.index('--csv') + 1])
                csv_path.write_text('path,status\n' + str(output.resolve()) + ',valid\n')

            with patch('allopep.shutil.which', return_value='/usr/bin/ost'), \
                 patch('allopep.subprocess.run', side_effect=write_result):
                valid = run_openstructure_validation(
                    [(output, 1, 0)], score_dir, 'target', None, input_pdb)
            manifest = json.loads((score_dir / 'target_pepglad_validation_manifest.json').read_text())
            self.assertEqual(valid, {str(output.resolve())})
            self.assertEqual(manifest[0]['pep_seq'], 'LHLKA')
            self.assertEqual(manifest[0]['rec_chains'], ['A', 'B'])
            self.assertEqual(manifest[0]['input_pdb'], str(input_pdb.resolve()))

    def test_existing_outputs_reject_missing_samples(self):
        with tempfile.TemporaryDirectory() as directory:
            site = (Path(directory) / 'target_PepGLAD_outputs' /
                    'target_allosteric_site_1_codesign')
            site.mkdir(parents=True)
            (site / 'target_0.pdb').touch()
            (site / 'summary.jsonl').write_text(
                '{"id":"target_0"}\n{"id":"target_1"}\n', encoding='utf-8'
            )
            with self.assertRaisesRegex(ValueError, 'missing=.*target_1'):
                find_existing_pepglad_outputs(Path(directory), 'target')

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
