"""Tests for PepGLAD validation and Vina ranking."""

import csv
import json
import os
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from allopep import (find_existing_pepglad_outputs, parse_args, prepare_input_pdb,
                     main, rank_pepglad_outputs, run_openstructure_validation,
                     run_pepglad_round, write_pymol_session)
from pepglad_outputs import DesignOutput, load_design_outputs
from vina_ranking import check_vina_dependencies
from vina_scoring.models import WorkflowConfig
from vina_scoring.results import ScoreSuccess, failure
from vina_scoring.scoring import PoseMetrics, SCORE_TERMS
from run_gaps import ScoredResidue, select_contexts, select_guidance_contexts


def design(path, site=1, sample=0, sequence='LHLKA', peptide_chain='B',
           receptor_chains=('A',)):
    return DesignOutput(path, site, sample, sequence,
                        peptide_chain, receptor_chains)


class VinaRankingTests(unittest.TestCase):
    def test_pymol_session_shows_original_receptor_and_pepglad_peptides(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / 'target output'
            root.mkdir()
            input_pdb = root / 'target.pdb'
            input_pdb.touch()
            generated = []
            for site in (1, 2, 3, 4):
                site_dir = root / f'target_allosteric_site_{site}_codesign'
                site_dir.mkdir()
                pocket_dir = root / 'target_allosteric_sites'
                pocket_dir.mkdir(exist_ok=True)
                entries = [['A', [10, ' ']], ['A', [11, 'A']], ['B', [20, ' ']]]
                if site == 4:  # Read-only viewing of an older scored pocket.
                    entries = [entry + [0.8] for entry in entries]
                (pocket_dir / f'target_allosteric_site_{site}.json').write_text(
                    json.dumps(entries), encoding='utf-8')
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

            with patch('allopep._pymol_score_commands', return_value=[]) as scores, \
                 patch('allopep.subprocess.run', side_effect=build_session):
                session = write_pymol_session(generated, root, 'target')
            scores.assert_called_once_with(root, 'target', 'target')
            content = scripts[0]
            self.assertEqual(session.name, 'target_pepglad_aligned.pse')
            self.assertEqual(session.read_bytes(), b'PyMOL session')
            self.assertFalse(legacy_script.exists())
            self.assertFalse(list(root.glob('**/*.pml')))
            self.assertIn(f'load {json.dumps(str(input_pdb.resolve()))}, '
                          'target', content)
            self.assertIn('hide everything, all', content)
            self.assertIn('show cartoon, target', content)
            self.assertIn('spectrum b, blue_white_red, target '
                          'and polymer.protein and b > -0.5, minimum=0, maximum=1', content)
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
            (root / 'target.pdb').touch()
            session = root / 'target_pepglad_aligned.pse'
            session.write_bytes(b'previous session')
            with patch('allopep._pymol_score_commands', return_value=[]), \
                 patch('allopep.subprocess.run', return_value=SimpleNamespace(
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

    def test_pepglad_round_uses_absolute_paths_and_an_isolated_log(self):
        args = SimpleNamespace(length_min=5, length_max=10)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            pocket = root / 'pocket.json'
            round_dir = root / 'round_000001'
            with patch('allopep.run_logged') as run:
                run_pepglad_round(args, '/input/target.pdb', pocket, round_dir, 16,
                                  deadline=1234)
            command = run.call_args.args[0]
            self.assertEqual(command[:7], [
                'conda', 'run', '--no-capture-output', '-n', 'PepGLAD',
                'python', '-m'])
            self.assertEqual(command[command.index('--pocket') + 1], str(pocket.resolve()))
            self.assertEqual(command[command.index('--out_dir') + 1], str(round_dir.resolve()))
            self.assertEqual(command[command.index('--n_samples') + 1], '16')
            self.assertEqual(run.call_args.kwargs['deadline'], 1234)
            self.assertEqual(run.call_args.kwargs['log_file'], round_dir / 'pepglad.log')
            self.assertNotIn('--score-guidance-strength', command)

    @staticmethod
    def score(job, total):
        """Build a scorer outcome while keeping adapter tests deterministic."""
        receptor = {'chains': list(job.source.receptor_chains),
                    'warnings': [], 'source_sha256': 'input-digest'}
        peptide = {'warnings': [], 'source_sha256': 'input-digest'}
        metrics = PoseMetrics(2, 2, 1, 0, ((0.0, 0.0, 0.0),))
        scores = {term: total if term in ('total', 'ligand_inter') else 0.0
                  for term in SCORE_TERMS}
        return ScoreSuccess(job, WorkflowConfig(), receptor, peptide,
                            metrics, scores, (0, 0, 0), (10, 10, 10),
                            'union_of_peptide_bounds', '1.2.7', '0.8.0')

    def test_ranks_selected_designs_by_total_with_stable_ties(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            valid = [design(root / 'pose_0.pdb', 2, 1),
                     design(root / 'pose_1.pdb', 1, 1),
                     design(root / 'pose_2.pdb', 1, 0),
                     design(root / 'pose_3.pdb', 3, 0, peptide_chain='C',
                            receptor_chains=('A', 'B')),
                     design(root / 'pose_4.pdb', 3, 1)]
            for output in valid:
                output.path.write_text(''.join(
                    f'ATOM      1  CA  ALA {chain}   1      0.000   0.000   0.000\n'
                    for chain in (*output.receptor_chains,
                                  output.peptide_chain)))
            (root / 'unselected.pdb').touch()

            def score_jobs(jobs):
                self.assertEqual([job.job_id for job in jobs],
                                 [output.job_id for output in valid])
                self.assertEqual(jobs[0].source.receptor_chains, ('A',))
                self.assertEqual(jobs[3].source.receptor_chains, ('A', 'B'))
                self.assertEqual(jobs[3].source.peptide_chain, 'C')
                return [self.score(job, -8.0 if index == 3 else -4.0)
                        if index < 4 else failure(
                            job, WorkflowConfig(), ValueError('bad residue'))
                        for index, job in enumerate(jobs)]

            with patch('vina_scoring.batch.BatchRunner.run',
                       side_effect=score_jobs):
                ranked = rank_pepglad_outputs(valid, root, 'target')
            self.assertEqual([(r['site'], r['sample']) for r in ranked],
                             [(3, 0), (1, 0), (1, 1), (2, 1)])
            self.assertEqual([r['rank'] for r in ranked], [1, 2, 3, 4])
            with (root / 'target_vina_manifest.csv').open(newline='') as handle:
                manifest = list(csv.DictReader(handle))
            self.assertEqual(len(manifest), len(valid))
            self.assertEqual(manifest[0]['box_padding'], '5.0')
            self.assertEqual(manifest[0]['hydrogen_policy'], 'keep')
            self.assertEqual(manifest[0]['spacing'], '0.375')
            self.assertEqual(manifest[3]['receptor_chains'], 'A;B')
            with (root / 'target_vina_results.csv').open(newline='') as handle:
                results = list(csv.DictReader(handle))
            self.assertEqual([r['status'] for r in results],
                             ['ok', 'ok', 'ok', 'ok', 'error'])
            self.assertIn('bad residue', results[-1]['error'])
            self.assertEqual(results[0]['protocol_id'],
                             'rigid_meeko_vina_residue_sum_v1')
            self.assertEqual(results[0]['vina_version'], '1.2.7')
            with (root / 'target_vina_ranked.csv').open(newline='') as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual([r['id'] for r in rows],
                             [r['id'] for r in ranked])

    def test_all_failures_keep_diagnostics_and_clear_stale_reports(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            pdb = root / 'target_0.pdb'
            pdb.write_text('ATOM      1  CA  ALA A   1      0.000   0.000   0.000\n'
                           'ATOM      2  CA  ALA B   1      0.000   0.000   0.000\n')
            for suffix in ('manifest', 'results', 'ranked'):
                (root / f'target_vina_{suffix}.csv').write_text('old report\n')

            def fail_jobs(jobs):
                self.assertFalse((root / 'target_vina_ranked.csv').exists())
                self.assertFalse((root / 'target_vina_results.csv').exists())
                return [failure(job, WorkflowConfig(), RuntimeError('failed'))
                        for job in jobs]

            with patch('vina_scoring.batch.BatchRunner.run',
                       side_effect=fail_jobs), self.assertRaisesRegex(
                           RuntimeError, 'All 1 Vina scoring jobs failed'):
                rank_pepglad_outputs([design(pdb)], root, 'target')
            self.assertFalse((root / 'target_vina_ranked.csv').exists())
            with (root / 'target_vina_results.csv').open(newline='') as handle:
                result = next(csv.DictReader(handle))
            self.assertIn('failed', result['error'])

    def test_nonfinite_total_is_a_scoring_failure(self):
        for total in (float('nan'), float('inf'), -float('inf')):
            with self.subTest(total=total), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                pdb = root / 'target.pdb'
                pdb.write_text('ATOM      1  CA  ALA A   1      0.000   0.000   0.000\n'
                               'ATOM      2  CA  ALA B   1      0.000   0.000   0.000\n')
                with patch('vina_scoring.batch.BatchRunner.run', side_effect=(
                        lambda jobs: [self.score(job, total) for job in jobs])), \
                     self.assertRaisesRegex(RuntimeError, 'All 1 Vina'):
                    rank_pepglad_outputs([design(pdb)], root, 'target')
                with (root / 'target_vina_results.csv').open(newline='') as handle:
                    result = next(csv.DictReader(handle))
                self.assertEqual(result['status'], 'error')
                self.assertIn('non-finite', result['error'])

    def test_missing_vina_dependency_is_actionable(self):
        with patch('vina_ranking.importlib.import_module',
                   side_effect=ImportError('missing vina')), \
             self.assertRaisesRegex(RuntimeError, 'conda env update'):
            check_vina_dependencies()

    def test_main_checks_vina_before_preparing_input(self):
        with patch('allopep.check_vina_dependencies',
                   side_effect=RuntimeError('Vina missing')), \
             patch('allopep.prepare_input_pdb') as prepare, \
             self.assertRaisesRegex(RuntimeError, 'Vina missing'):
            main(['target.pdb'])
        prepare.assert_not_called()

    def test_validation_only_skips_scoring_dependency_check(self):
        with patch('allopep.check_vina_dependencies') as check, \
             patch('allopep.prepare_input_pdb', side_effect=FileNotFoundError), \
             self.assertRaises(FileNotFoundError):
            main(['target.pdb', '--validate-existing'])
        check.assert_not_called()

    def test_max_rank_cli_defaults_to_all_and_accepts_a_limit(self):
        self.assertEqual(parse_args(['target.pdb']).max_rank, 'all')
        self.assertEqual(parse_args(['target.pdb', '--max_rank', 'all']).max_rank, 'all')
        self.assertEqual(parse_args(['target.pdb', '--max_rank', '2']).max_rank, 2)

    def test_new_design_uses_effective_gaps_selector_defaults(self):
        residues = [
            ScoredResidue('A', site * 10 + offset + 1, ' ', 0.9 - site * 0.1,
                          (site * 20.0 + offset * 2.0, 0.0, 0.0), 0.5)
            for site in range(4) for offset in range(2)
        ]
        expected = select_contexts(residues)
        with tempfile.TemporaryDirectory() as directory:
            previous_dir = os.getcwd()
            try:
                os.chdir(directory)
                with patch('allopep.check_vina_dependencies', return_value='/fake/ia'), \
                     patch('allopep.prepare_input_pdb', return_value='/fake/target.pdb'), \
                     patch('allopep.run_gaps.run_gaps'), \
                     patch('allopep.run_gaps.load_scored_residues',
                           return_value=residues), \
                     patch('allopep.run_gaps.select_guidance_contexts',
                           wraps=select_guidance_contexts) as selector, \
                     patch('allopep.run_pepglad', return_value=SimpleNamespace(
                         selected=[], shortfall=True, summary_path='summary.json')) as pepglad, \
                     patch('allopep.write_pymol_session', return_value='session.pse'), \
                     patch('allopep.run_openstructure_validation', return_value=[]):
                    main(['target.pdb'])
                self.assertEqual(selector.call_args.kwargs, {'max_rank': 'all'})
                self.assertEqual(pepglad.call_args.args[4], expected)
                pockets = pepglad.call_args.args[4]
                self.assertEqual(len(pockets), 4)
                self.assertEqual(pockets[0].seed, residues[0])
                self.assertEqual(pockets[0].residues, tuple(residues[:2]))
            finally:
                os.chdir(previous_dir)

    def test_chain_mismatch_records_error_without_scoring(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            pdb = root / 'target_0.pdb'
            pdb.write_text('ATOM      1  CA  ALA A   1      0.000   0.000   0.000\n')
            with patch('vina_scoring.batch.BatchRunner.run') as run, \
                 self.assertRaisesRegex(RuntimeError, 'All 1 Vina scoring jobs failed'):
                rank_pepglad_outputs([design(pdb)], root, 'target')
            run.assert_not_called()
            with (root / 'target_vina_results.csv').open(newline='') as handle:
                result = next(csv.DictReader(handle))
            self.assertIn('PDB chains', result['error'])

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
                       '--num_samples 0', '--max_num_peptides 0',
                       '--attempt-multiplier 0', '--attempt-multiplier -1',
                       '-n 5 --max_num_peptides 4', '--round-timeout-seconds 0',
                       '--peptide-chain AA'):
            with self.subTest(option=option), \
                 patch('allopep.prepare_input_pdb') as prepare:
                with self.assertRaises(ValueError):
                    main(['target.pdb', *option.split()])
                prepare.assert_not_called()


if __name__ == '__main__':
    unittest.main()
