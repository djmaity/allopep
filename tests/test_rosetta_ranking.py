"""Tests for PepGLAD validation and Rosetta ranking."""

import csv
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from allopep import (find_existing_pepglad_outputs, prepare_input_pdb,
                     main, rank_pepglad_outputs, run_openstructure_validation,
                     run_pepglad, write_pymol_session)
from pepglad_outputs import DesignOutput, load_design_outputs
from rosetta_scoring import parse_scorefile, resolve_rosetta_binary
from run_gaps import GuidanceContext, ScoredResidue


def design(path, site=1, sample=0, sequence='LHLKA', peptide_chain='B',
           receptor_chains=('A',)):
    return DesignOutput(path, site, sample, sequence,
                        peptide_chain, receptor_chains)


class RosettaRankingTests(unittest.TestCase):
    def test_pymol_session_shows_gaps_receptor_and_pepglad_peptides(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / 'target output'
            root.mkdir()
            gaps_pdb = root / 'target_GAPS_output.pdb'
            gaps_pdb.touch()
            generated = []
            for site in (1, 2, 3, 4):
                site_dir = root / f'target_allosteric_site_{site}_codesign'
                site_dir.mkdir()
                pocket_dir = root / 'target_allosteric_sites'
                pocket_dir.mkdir(exist_ok=True)
                (pocket_dir / f'target_allosteric_site_{site}.json').write_text(
                    json.dumps([['A', [10, ' '], 0.8],
                                ['A', [11, 'A'], 0.7],
                                ['B', [20, ' '], 0.6]]), encoding='utf-8')
                model = site_dir / 'target_0.pdb'
                model.touch()
                (site_dir / 'summary.jsonl').write_text(json.dumps({
                    'id': 'target_0', 'rec_chains': ['A', 'B'],
                    'pep_chain': 'C',
                }) + '\n', encoding='utf-8')
                generated.append(design(model, site, peptide_chain='C',
                                        receptor_chains=('A', 'B')))
            second_model = root / 'target_allosteric_site_1_codesign' / 'target_1.pdb'
            second_model.touch()
            generated.append(design(second_model, 1, sample=1,
                                    peptide_chain='C', receptor_chains=('A', 'B')))

            legacy_script = root / 'target_pepglad_aligned.pml'
            legacy_script.touch()
            scripts = []

            def build_session(command, **kwargs):
                self.assertEqual(command[:2], ['pymol', '-cqQ'])
                self.assertTrue(kwargs['capture_output'])
                content = (kwargs['cwd'] / command[2]).read_text(encoding='utf-8')
                scripts.append(content)
                self.assertEqual(content.splitlines()[-1], 'save session.pse')
                session = kwargs['cwd'] / 'session.pse'
                session.write_bytes(b'PyMOL session')
                return SimpleNamespace(returncode=0, stderr='', stdout='')

            with patch('allopep.subprocess.run', side_effect=build_session):
                session = write_pymol_session(generated, root, 'target')
            content = scripts[0]
            self.assertEqual(session.name, 'target_pepglad_aligned.pse')
            self.assertEqual(session.read_bytes(), b'PyMOL session')
            self.assertFalse(legacy_script.exists())
            self.assertFalse(list(root.glob('**/*.pml')))
            self.assertIn(f'load {json.dumps(str(gaps_pdb.resolve()))}, '
                          'target', content)
            self.assertIn('hide everything, all', content)
            self.assertIn('show cartoon, target', content)
            self.assertIn('spectrum b, blue_white_red, target '
                          'and polymer.protein, minimum=0, maximum=1', content)
            site_colors = {1: 'green', 2: 'cyan',
                           3: 'lightmagenta', 4: 'yellow'}
            for output in generated:
                name = f'site_{output.site}_peptide_{output.sample + 1}'
                self.assertIn(f'load {json.dumps(output.resolved_path)}, {name}',
                              content)
                self.assertIn(
                    f'align ({name} and name CA and (chain A or chain B)), '
                    '(target and name CA and (chain A or chain B)), '
                    'cycles=0', content)
                self.assertIn(f'remove {name} and not chain C', content)
                self.assertIn(f'show sticks, site_{output.site}_peptide_* '
                              'and chain C', content)
                self.assertIn(
                    f'color {site_colors[output.site]}, '
                    f'site_{output.site}_peptide_* and chain C', content)
                self.assertIn(
                    f'color atomic, site_{output.site}_peptide_* '
                    'and chain C and not elem C', content)
                surface = f'site_{output.site}'
                self.assertIn(f'create {surface}, target '
                              'and polymer.protein', content)
                self.assertIn(
                    f'show surface, {surface} and '
                    '((chain A and resi 10+11A) or (chain B and resi 20))',
                    content)
                self.assertIn(f'set transparency, 0.4, {surface}', content)
                self.assertIn(f'color {site_colors[output.site]}, {surface}',
                              content)
                self.assertIn(f'group sites, {surface}', content)
                self.assertNotIn(f'site_{output.site}_surface', content)
                self.assertIn(
                    f'group site_{output.site}_peptides, '
                    f'site_{output.site}_peptide_*',
                    content)
                self.assertNotIn(
                    f'group sites, site_{output.site}_peptides\n', content)
                self.assertIn(f'group peptides, site_{output.site}_peptides',
                              content)
            self.assertNotIn('allopep_', content)
            self.assertEqual(content.count('\nalign '), 5)
            self.assertNotIn('\npython\n', content)

    def test_pymol_failure_preserves_existing_session(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / 'target_GAPS_output.pdb').touch()
            session = root / 'target_pepglad_aligned.pse'
            session.write_bytes(b'previous session')
            with patch('allopep.subprocess.run', return_value=SimpleNamespace(
                    returncode=1, stderr='PyMOL failed', stdout='')):
                with self.assertRaisesRegex(RuntimeError, 'PyMOL failed'):
                    write_pymol_session([], root, 'target')
            self.assertEqual(session.read_bytes(), b'previous session')
            self.assertFalse(list(root.glob('.allopep-pymol-*')))

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
        residue = ScoredResidue('A', 1, ' ', 0.873241, (0.0, 0.0, 0.0), 0.5)
        context = GuidanceContext(residue, (residue,))
        with tempfile.TemporaryDirectory(dir='.') as directory:
            out_dir = Path(directory) / 'output'
            output_pdb = (
                out_dir / 'target_PepGLAD_outputs'
                / 'target_allosteric_site_1_codesign' / 'target_0.pdb'
            )
            output_pdb.parent.mkdir(parents=True)
            output_pdb.touch()
            (output_pdb.parent / 'summary.jsonl').write_text(json.dumps({
                'id': 'target_0', 'pep_seq': 'LHLKA',
                'pep_chain': 'B', 'rec_chains': ['A'],
            }) + '\n', encoding='utf-8')

            with patch('allopep.subprocess.run') as subprocess_run:
                generated = run_pepglad(
                    args, '/input/target.pdb', out_dir, 'target', [context]
                )

            pocket_file = (
                out_dir / 'target_allosteric_sites'
                / 'target_allosteric_site_1.json'
            )
            self.assertEqual(pocket_file.read_text(), '[["A", [1, " "]]]')
            self.assertEqual(generated, [design(output_pdb)])
            command = subprocess_run.call_args.args[0]
            self.assertEqual(command[:7], [
                'conda', 'run', '--no-capture-output', '-n', 'PepGLAD',
                'python', '-m',
            ])
            self.assertEqual(command[command.index('--pocket') + 1], f'../{pocket_file}')
            self.assertEqual(command[command.index('--out_dir') + 1], f'../{output_pdb.parent}')
            self.assertNotIn('--score-guidance-strength', command)

    def test_resume_reuses_only_relaxed_matching_site(self):
        residue = ScoredResidue('A', 1, ' ', 0.873241, (0.0, 0.0, 0.0), 0.5)
        context = GuidanceContext(residue, (residue,))
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            pocket = root / 'target_allosteric_sites' / 'target_allosteric_site_1.json'
            pocket.parent.mkdir()
            pocket.write_text(json.dumps(context.as_pepglad_pocket()))
            site = root / 'target_PepGLAD_outputs' / 'target_allosteric_site_1_codesign'
            site.mkdir(parents=True)
            pdb = site / 'target_0.pdb'
            pdb.write_text('REMARK   1 CREATED WITH OPENMM 8.0\n')
            (site / 'summary.jsonl').write_text(json.dumps({
                'id': 'target_0', 'pep_seq': 'LHLKA',
                'pep_chain': 'B', 'rec_chains': ['A'],
            }) + '\n')
            args = SimpleNamespace(resume=True, num_samples=1)
            with patch('allopep.subprocess.run') as subprocess_run:
                self.assertEqual(run_pepglad(
                    args, '/input/target.pdb', root, 'target', [context]),
                    [design(pdb)])
                subprocess_run.assert_not_called()
            pdb.write_text('ATOM      1\n')
            with self.assertRaisesRegex(RuntimeError, 'unrelaxed samples'):
                run_pepglad(args, '/input/target.pdb', root, 'target', [context])

    def test_ranks_only_validated_structures_by_lower_rosetta_score(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            structures = [root / f'target_{i}.pdb' for i in range(4)]
            valid = [design(structures[0]),
                     design(structures[2], 2, 0, peptide_chain='C',
                            receptor_chains=('A', 'B')),
                     design(structures[3], 2, 1, peptide_chain='C',
                            receptor_chains=('A', 'B'))]
            for output in valid:
                output.path.write_text(''.join(
                    f'ATOM      1  CA  ALA {chain}   1      0.000   0.000   0.000\n'
                    for chain in (*output.receptor_chains, output.peptide_chain)))

            def write_results(command, **kwargs):
                self.assertEqual(command[0], '/fake/InterfaceAnalyzer')
                self.assertEqual(command[command.index('-score:weights') + 1],
                                 str(Path('weights/optimized_rosetta_26_terms.wts').resolve()))
                self.assertEqual(command[command.index('-pack_input') + 1], 'false')
                self.assertEqual(command[command.index('-pack_separated') + 1], 'false')
                self.assertEqual(kwargs['timeout'], 1800)
                pdb = Path(command[command.index('-s') + 1])
                expected_interface = 'A_B' if pdb == structures[0] else 'AB_C'
                self.assertEqual(command[command.index('-interface') + 1],
                                 expected_interface)
                if pdb == structures[3]:
                    return SimpleNamespace(returncode=1, stdout='', stderr='bad residue')
                score = -8.0 if pdb == structures[2] else -4.0
                output = Path(command[command.index('-out:file:score_only') + 1])
                output.write_text('SCORE: dG_separated dSASA_int description\n'
                                  f'SCORE: {score} 120.0 test\n')
                return SimpleNamespace(returncode=0, stdout='ok', stderr='')

            with patch('rosetta_scoring.subprocess.run', side_effect=write_results) as run:
                ranked = rank_pepglad_outputs(
                    valid, root, 'target', '/fake/InterfaceAnalyzer')
            self.assertEqual([(r['rank'], r['site'], r['sample']) for r in ranked],
                             [(1, 2, 0), (2, 1, 0)])
            self.assertEqual([r['rosetta_dG_separated'] for r in ranked],
                             [-8.0, -4.0])
            self.assertEqual(run.call_count, 3)
            with (root / 'target_rosetta_manifest.csv').open(newline='') as handle:
                jobs = list(csv.DictReader(handle))
            self.assertEqual([job['id'] for job in jobs],
                             ['site_1_sample_0', 'site_2_sample_0',
                              'site_2_sample_1'])
            with (root / 'target_rosetta_ranked.csv').open(newline='') as handle:
                csv_rows = list(csv.DictReader(handle))
            self.assertEqual([row['sample'] for row in csv_rows], ['0', '0'])
            self.assertEqual([row['site'] for row in csv_rows], ['2', '1'])
            with (root / 'target_rosetta_results.csv').open(newline='') as handle:
                results = list(csv.DictReader(handle))
            self.assertEqual([row['status'] for row in results],
                             ['ok', 'ok', 'error'])
            self.assertIn('bad residue', results[2]['error'])
            self.assertEqual(results[1]['rosetta_dSASA_int'], '120.0')
            self.assertEqual(results[1]['rosetta_binary'], '/fake/InterfaceAnalyzer')
            self.assertEqual(results[1]['weights_sha256'],
                             'e193a7acc0cfe0918a15662cac9d021f232797da81d57aa89e6a1c6005983044')
            self.assertTrue(Path(results[1]['scorefile']).is_file())
            self.assertTrue(Path(results[1]['log_file']).is_file())

    def test_all_rosetta_failures_keep_diagnostics_and_stop_ranking(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            pdb = root / 'target_0.pdb'
            pdb.write_text('ATOM      1  CA  ALA A   1      0.000   0.000   0.000\n'
                           'ATOM      2  CA  ALA B   1      0.000   0.000   0.000\n')

            with patch('rosetta_scoring.subprocess.run', return_value=SimpleNamespace(
                    returncode=1, stdout='', stderr='Rosetta failed')), \
                 self.assertRaisesRegex(RuntimeError, 'All 1 Rosetta scoring jobs failed'):
                rank_pepglad_outputs(
                    [design(pdb)], root, 'target', '/fake/InterfaceAnalyzer')
            self.assertTrue((root / 'target_rosetta_results.csv').is_file())
            self.assertFalse((root / 'target_rosetta_ranked.csv').exists())

    def test_stale_scorefile_is_not_reused(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            pdb = root / 'target_0.pdb'
            pdb.write_text('ATOM      1  CA  ALA A   1      0.000   0.000   0.000\n'
                           'ATOM      2  CA  ALA B   1      0.000   0.000   0.000\n')
            run_dir = root / 'rosetta_runs' / 'site_1_sample_0'
            run_dir.mkdir(parents=True)
            (run_dir / 'interface.sc').write_text(
                'SCORE: dG_separated description\nSCORE: -5.64 old\n')
            (root / 'target_rosetta_ranked.csv').write_text('old rank\n')

            with patch('rosetta_scoring.subprocess.run', return_value=SimpleNamespace(
                    returncode=0, stdout='', stderr='')), \
                 self.assertRaisesRegex(RuntimeError, 'All 1 Rosetta scoring jobs failed'):
                rank_pepglad_outputs([design(pdb)], root, 'target',
                                     '/fake/InterfaceAnalyzer')
            self.assertFalse((run_dir / 'interface.sc').exists())
            self.assertFalse((root / 'target_rosetta_ranked.csv').exists())
            with (root / 'target_rosetta_results.csv').open(newline='') as handle:
                result = next(csv.DictReader(handle))
            self.assertIn('scorefile missing', result['error'])

    def test_scorefile_rejects_nonfinite_interface_energy(self):
        with tempfile.TemporaryDirectory() as directory:
            scorefile = Path(directory) / 'interface.sc'
            scorefile.write_text('SCORE: dG_separated description\n'
                                 'SCORE: nan test\n')
            with self.assertRaisesRegex(ValueError, 'Non-finite dG_separated'):
                parse_scorefile(scorefile)

    def test_missing_rosetta_binary_is_reported(self):
        with patch('rosetta_scoring.shutil.which', return_value=None), \
             self.assertRaisesRegex(FileNotFoundError, 'InterfaceAnalyzer is not on PATH'):
            resolve_rosetta_binary()

    def test_main_checks_rosetta_before_preparing_input(self):
        with patch('allopep.resolve_rosetta_binary', side_effect=FileNotFoundError(
                'InterfaceAnalyzer missing')), \
             patch('allopep.prepare_input_pdb') as prepare, \
             self.assertRaisesRegex(FileNotFoundError, 'InterfaceAnalyzer missing'):
            main(['target.pdb'])
        prepare.assert_not_called()

    def test_chain_mismatch_records_error_without_running_rosetta(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            pdb = root / 'target_0.pdb'
            pdb.write_text('ATOM      1  CA  ALA A   1      0.000   0.000   0.000\n')
            with patch('rosetta_scoring.subprocess.run') as run, \
                 self.assertRaisesRegex(RuntimeError, 'All 1 Rosetta scoring jobs failed'):
                rank_pepglad_outputs([design(pdb)], root, 'target',
                                     '/fake/InterfaceAnalyzer')
            run.assert_not_called()
            with (root / 'target_rosetta_results.csv').open(newline='') as handle:
                result = next(csv.DictReader(handle))
            self.assertIn('PDB chains', result['error'])
            self.assertIn('PDB chains', Path(result['log_file']).read_text())

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
                generated = load_design_outputs([(output, 1, 0)])
                valid = run_openstructure_validation(
                    generated, score_dir, 'target', None, input_pdb)
            manifest = json.loads((score_dir / 'target_pepglad_validation_manifest.json').read_text())
            self.assertEqual(valid, generated)
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
            (site / 'target_1.pdb').touch()
            (site / 'summary.jsonl').write_text(''.join(json.dumps({
                'id': f'target_{sample}', 'pep_seq': 'LHLKA',
                'pep_chain': 'B', 'rec_chains': ['A'],
            }) + '\n' for sample in (0, 1)))
            self.assertEqual(find_existing_pepglad_outputs(
                Path(directory), 'target'), [
                    design(site / 'target_0.pdb'),
                    design(site / 'target_1.pdb', sample=1),
                ])

    def test_existing_outputs_accept_reported_invalid_sample(self):
        with tempfile.TemporaryDirectory() as directory:
            site = (Path(directory) / 'target_PepGLAD_outputs' /
                    'target_allosteric_site_1_codesign')
            site.mkdir(parents=True)
            (site / 'target_0.pdb').touch()
            (site / 'summary.jsonl').write_text(''.join(json.dumps({
                'id': f'target_{sample}', 'pep_seq': 'LHLKA',
                'pep_chain': 'B', 'rec_chains': ['A'],
                **({'status': 'invalid', 'error': 'bad coordinates'}
                   if sample == 1 else {}),
            }) + '\n' for sample in (0, 1)))
            self.assertEqual(find_existing_pepglad_outputs(
                Path(directory), 'target'), [design(site / 'target_0.pdb')])
            (site / 'target_1.pdb').touch()
            with self.assertRaisesRegex(ValueError, 'extra=.*target_1'):
                find_existing_pepglad_outputs(Path(directory), 'target')

    def test_summary_loader_rejects_duplicate_and_missing_ids(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            output = root / 'target_0.pdb'
            output.touch()
            summary = root / 'summary.jsonl'
            record = {'id': 'target_0', 'pep_seq': 'LHLKA',
                      'pep_chain': 'B', 'rec_chains': ['A']}
            summary.write_text(json.dumps(record) + '\n' + json.dumps(record) + '\n')
            with self.assertRaisesRegex(ValueError, 'Duplicate PepGLAD IDs'):
                load_design_outputs([(output, 1, 0)])
            record['id'] = 'target_1'
            summary.write_text(json.dumps(record) + '\n')
            with self.assertRaisesRegex(ValueError, 'target_0 is missing'):
                load_design_outputs([(output, 1, 0)])

    def test_invalid_cli_options_fail_before_pipeline_stages(self):
        for option in ('--length_min 0', '--length_min 8 --length_max 7',
                       '--num_samples 0',
                       '--peptide-chain AA'):
            with self.subTest(option=option), \
                 patch('allopep.prepare_input_pdb') as prepare:
                with self.assertRaises(ValueError):
                    main(['target.pdb', *option.split()])
                prepare.assert_not_called()


if __name__ == '__main__':
    unittest.main()
