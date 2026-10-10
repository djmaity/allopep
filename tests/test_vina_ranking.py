"""Tests for PepGLAD validation and Vina ranking."""

import csv
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from allopep import (find_existing_pepglad_outputs, prepare_input_pdb,
                     rank_pepglad_outputs, run_openstructure_validation,
                     run_pepglad, write_pymol_script)
from run_gaps import SpatialPatch


class VinaRankingTests(unittest.TestCase):
    def test_pymol_script_aligns_every_complex_to_input_receptor(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / 'target output'
            root.mkdir()
            input_pdb = root / 'target.pdb'
            input_pdb.touch()
            generated = []
            for site in (1, 2):
                site_dir = root / f'target_allosteric_site_{site}_codesign'
                site_dir.mkdir()
                model = site_dir / 'target_0.pdb'
                model.touch()
                (site_dir / 'summary.jsonl').write_text(json.dumps({
                    'id': 'target_0', 'rec_chains': ['A', 'B'],
                    'pep_chain': 'C',
                }) + '\n', encoding='utf-8')
                generated.append((model, site, 0))

            script = write_pymol_script(generated, input_pdb, root, 'target')
            content = script.read_text(encoding='utf-8')
            self.assertEqual(script.name, 'target_pepglad_aligned.pml')
            self.assertIn(f'load {json.dumps(str(input_pdb.resolve()))}, '
                          'allopep_target_input', content)
            for model, site, _ in generated:
                name = f'allopep_target_site_{site}_sample_0'
                self.assertIn(f'load {json.dumps(str(model.resolve()))}, {name}',
                              content)
                self.assertIn(
                    f'align ({name} and name CA and (chain A or chain B)), '
                    '(allopep_target_input and name CA and (chain A or chain B)), '
                    'cycles=0', content)
                self.assertIn(f'show sticks, allopep_target_site_{site}_* '
                              'and chain C', content)
            self.assertEqual(content.count('\nalign '), 2)
            self.assertNotIn('\npython\n', content)

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

    def test_ranks_only_validated_structures_by_lower_vina_score(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            structures = [root / f'target_{i}.pdb' for i in range(4)]
            generated = [(structures[0], 1, 0), (structures[1], 1, 1),
                         (structures[2], 2, 0), (structures[3], 2, 1)]
            valid = {str(structures[0].resolve()): 'B',
                     str(structures[2].resolve()): 'C',
                     str(structures[3].resolve()): 'C'}

            def write_results(command, check):
                self.assertTrue(check)
                self.assertIn('vina_scoring/vina_score.py',
                              command[1].replace('\\', '/'))
                manifest = Path(command[command.index('--manifest') + 1])
                with manifest.open(newline='') as handle:
                    jobs = list(csv.DictReader(handle))
                self.assertEqual([job['id'] for job in jobs],
                                 ['site_1_sample_0', 'site_2_sample_0',
                                  'site_2_sample_1'])
                output = Path(command[command.index('--output') + 1])
                with output.open('w', newline='') as handle:
                    writer = csv.DictWriter(handle, fieldnames=(
                        'id', 'status', 'score_total', 'score_ligand_inter',
                        'vina_version', 'meeko_version', 'error'))
                    writer.writeheader()
                    writer.writerows([
                        {'id': 'site_1_sample_0', 'status': 'ok',
                         'score_total': -4.0, 'score_ligand_inter': -4.2,
                         'vina_version': '1.2.7', 'meeko_version': '0.8.0'},
                        {'id': 'site_2_sample_0', 'status': 'ok',
                         'score_total': -8.0, 'score_ligand_inter': -8.1,
                         'vina_version': '1.2.7', 'meeko_version': '0.8.0'},
                        {'id': 'site_2_sample_1', 'status': 'error',
                         'error': 'Meeko preparation failed'},
                    ])

            with patch('allopep.subprocess.run', side_effect=write_results) as run:
                ranked = rank_pepglad_outputs(
                    generated, valid, root, 'target', box_padding=6,
                    hydrogen_policy='rebuild')
            self.assertEqual([(r['rank'], r['site'], r['sample']) for r in ranked],
                             [(1, 2, 0), (2, 1, 0)])
            self.assertEqual([r['vina_score_kcal_mol'] for r in ranked],
                             [-8.0, -4.0])
            self.assertEqual(run.call_count, 1)
            self.assertIn('rebuild', run.call_args.args[0])
            with (root / 'target_vina_ranked.csv').open(newline='') as handle:
                csv_rows = list(csv.DictReader(handle))
            self.assertEqual([row['sample'] for row in csv_rows], ['0', '0'])
            self.assertEqual([row['site'] for row in csv_rows], ['2', '1'])
            self.assertEqual(csv_rows[0]['vina_hydrogen_policy'], 'rebuild')
            with (root / 'target_vina_results.csv').open(newline='') as handle:
                results = list(csv.DictReader(handle))
            self.assertEqual([row['status'] for row in results],
                             ['ok', 'ok', 'error'])

    def test_all_vina_failures_keep_diagnostics_and_stop_ranking(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            pdb = root / 'target_0.pdb'

            def write_failure(command, check):
                output = Path(command[command.index('--output') + 1])
                output.write_text(
                    'id,status,error\n'
                    'site_1_sample_0,error,Meeko preparation failed\n',
                    encoding='utf-8')

            with patch('allopep.subprocess.run', side_effect=write_failure), \
                 self.assertRaisesRegex(RuntimeError, 'All 1 Vina scoring jobs failed'):
                rank_pepglad_outputs(
                    [(pdb, 1, 0)], {str(pdb.resolve()): 'B'}, root, 'target')
            self.assertTrue((root / 'target_vina_results.csv').is_file())
            self.assertFalse((root / 'target_vina_ranked.csv').exists())

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
                csv_path.write_text('path,status,peptide_chain\n' +
                                    str(output.resolve()) + ',valid,C\n')

            with patch('allopep.shutil.which', return_value='/usr/bin/ost'), \
                 patch('allopep.subprocess.run', side_effect=write_result):
                valid = run_openstructure_validation(
                    [(output, 1, 0)], score_dir, 'target', None, input_pdb)
            manifest = json.loads((score_dir / 'target_pepglad_validation_manifest.json').read_text())
            self.assertEqual(valid, {str(output.resolve()): 'C'})
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


if __name__ == '__main__':
    unittest.main()
