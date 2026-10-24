"""Behavioral tests for bounded multi-site sampling and saved-round recovery."""

import csv
import json
import os
import tempfile
import time
import unittest
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from allopep.cli import (find_existing_pepglad_outputs, main, parse_args,
                         sampling_config)
from allopep.fileio import atomic_json, sha256
from allopep.layout import OutputLayout
from allopep.sampling import (SamplingConfig, SamplingController,
                              next_round_size, sampling_lock)
from allopep.recover_relaxation import inspect_site


FIELDS = ('path', 'site', 'sample', 'status', 'reasons', 'peptide_chain',
          'sha256', 'input_sha256')


class SamplingTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.pdb = self.root / 'target.pdb'
        self.pdb.write_text('input receptor\n')
        self.out = self.root / 'output'
        self.layout = OutputLayout.at(self.out, 'target')
        self.config = SamplingConfig(num_samples=1, max_attempts=32,
                                     round_timeout_seconds=60, length_min=5,
                                     length_max=10, peptides_per_site=True)
        self.contexts = [SimpleNamespace(as_pepglad_pocket=lambda: [['A', [1, ' ']]])]
        self.events = []
        self.validation_calls = []
        self.outcome = lambda site, sample: 'invalid'
        self.generation_failures = {}

    def generate(self, config, pdb, pocket, directory, count, *, deadline):
        metadata = json.loads((directory / 'round.json').read_text())
        site, start = metadata['site'], metadata['start']
        self.events.append((site, start, count))
        self.assertLess(time.monotonic(), deadline)
        records = []
        for index in range(count):
            record = {'id': f'target_{index}', 'pep_seq': 'LHLKA',
                      'pep_chain': 'B', 'rec_chains': ['A']}
            failure = self.generation_failures.get((site, start + index))
            if failure:
                record.update(status=failure, error='candidate failure')
                if failure == 'error':
                    quarantine = directory / 'invalid_designs'
                    quarantine.mkdir(exist_ok=True)
                    (quarantine / f'target_{index}.pdb').write_text('failed minimization\n')
                    atomic_json(directory / f'target_{index}.pdb.relax_error.json',
                                {'error': 'candidate failure'})
            else:
                (directory / f'target_{index}.pdb').write_text(
                    f'REMARK   1 CREATED WITH OPENMM 8.0\n{site}/{start + index}\n')
            records.append(record)
        (directory / 'summary.jsonl').write_text(
            ''.join(json.dumps(record) + '\n' for record in records))

    def validate(self, designs, directory, layout, peptide_chain, input_pdb, **kwargs):
        self.validation_calls.append([(design.site, design.sample) for design in designs])
        rows = []
        for design in designs:
            status = self.outcome(design.site, design.sample)
            rows.append(dict(path=design.resolved_path, site=design.site, sample=design.sample,
                             status=status, reasons='' if status == 'valid' else 'test failure',
                             peptide_chain=design.peptide_chain, sha256=sha256(design.path),
                             input_sha256=sha256(input_pdb)))
        csv_path, jsonl_path, manifest_path = layout.validation_reports(directory)
        with csv_path.open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=FIELDS)
            writer.writeheader()
            writer.writerows(rows)
        jsonl_path.write_text(
            ''.join(json.dumps({'result': row, 'details': {}}) + '\n' for row in rows))
        atomic_json(manifest_path,
                    [{'path': design.resolved_path} for design in designs])
        return [design for design, row in zip(designs, rows) if row['status'] == 'valid']

    def controller(self, generate=None, validate=None):
        return SamplingController(self.config, str(self.pdb), self.layout,
                                  self.contexts, generate or self.generate,
                                  validate or self.validate)

    def state(self):
        return json.loads((self.out / 'target_sampling_summary.json').read_text())

    def test_cli_defaults_and_large_default_cap(self):
        parsed = parse_args(['target.pdb'])
        self.assertEqual(parsed.num_samples, 10)
        self.assertEqual(parsed.attempt_multiplier, 1000)
        self.assertIsNone(parsed.max_num_peptides)
        self.assertFalse(sampling_config(parsed).peptides_per_site)
        self.assertTrue(sampling_config(parse_args([
            'target.pdb', '--peptides-per-site'])).peptides_per_site)
        with patch('allopep.cli.check_vina_dependencies', side_effect=FileNotFoundError) as resolve:
            with self.assertRaises(FileNotFoundError):
                main(['target.pdb', '-n', '3'])
        resolve.assert_called_once()

    def test_cli_multiplier_changes_budget_and_absolute_cap_takes_precedence(self):
        cases = ((['--attempt-multiplier', '500'], 1500),
                 (['--attempt-multiplier', '2000'], 6000),
                 (['--attempt-multiplier', '1'], 3),
                 (['--attempt-multiplier', '500', '--max_num_peptides', '20'], 20))
        for options, expected in cases:
            with self.subTest(options=options):
                parsed = parse_args(['target.pdb', '-n', '3', *options])
                self.assertEqual(sampling_config(parsed).max_attempts, expected)

    def test_adaptive_sizes_and_clipping(self):
        self.assertEqual(next_round_size(10, 0, 0, 8, 100), 8)
        self.assertEqual(next_round_size(10, 8, 0, 8, 100), 16)
        self.assertEqual(next_round_size(10, 24, 0, 16, 100), 32)
        self.assertEqual(next_round_size(10, 56, 0, 32, 1000), 64)
        self.assertEqual(next_round_size(10, 120, 0, 64, 1000), 64)
        self.assertEqual(next_round_size(10, 120, 9, 64, 1000), 16)
        self.assertEqual(next_round_size(10, 24, 9, 16, 1000), 8)
        self.assertEqual(next_round_size(10, 24, 0, 16, 29), 5)
        self.assertEqual(next_round_size(10, 29, 0, 5, 29), 0)
        self.assertEqual(next_round_size(10, 40, 10, 16, 100), 0)

    def test_success_after_several_rounds_and_exact_target_at_cap(self):
        self.config = replace(self.config, num_samples=4)
        self.config = replace(self.config, max_attempts=40)
        self.outcome = lambda site, sample: 'valid' if sample in (0, 8, 24, 25) else 'invalid'
        result = self.controller().run()
        self.assertFalse(result.shortfall)
        self.assertEqual([design.sample for design in result.selected], [0, 8, 24, 25])
        self.assertEqual(self.events, [(1, 0, 8), (1, 8, 16), (1, 24, 16)])
        site = self.state()['sites'][0]
        self.assertEqual(site['attempted'], 40)
        self.assertEqual(site['completion_reason'], 'target_reached')
        csv_path = self.out / 'target_vina_score/target_pepglad_validation.csv'
        with csv_path.open() as handle:
            self.assertEqual(len(list(csv.DictReader(handle))), 40)

    def test_zero_acceptance_uses_leftover_budget_below_eight(self):
        self.config = replace(self.config, max_attempts=29)
        result = self.controller().run()
        self.assertTrue(result.shortfall)
        self.assertEqual(result.selected, [])
        self.assertEqual(self.events, [(1, 0, 8), (1, 8, 16), (1, 24, 5)])
        self.assertEqual(self.state()['sites'][0]['failures']['validation_invalid'], 29)

    def test_round_robin_and_independent_site_completion(self):
        self.contexts *= 3
        self.config = replace(self.config, max_attempts=29)
        self.outcome = lambda site, sample: 'valid' if (site, sample) in ((1, 0), (3, 8)) else 'invalid'
        result = self.controller().run()
        self.assertEqual(self.events, [(1, 0, 8), (2, 0, 8), (3, 0, 8),
                                       (2, 8, 16), (3, 8, 16), (2, 24, 5)])
        self.assertEqual([(d.site, d.sample) for d in result.selected], [(1, 0), (3, 8)])
        self.assertTrue(result.shortfall)

    def test_total_mode_pools_sites_and_stops_at_target(self):
        self.contexts *= 3
        self.config = replace(self.config, num_samples=3,
                              peptides_per_site=False)
        self.outcome = lambda site, sample: (
            'valid' if sample == 0 else 'invalid')
        result = self.controller().run()
        self.assertFalse(result.shortfall)
        self.assertEqual(self.events, [(1, 0, 8), (2, 0, 8), (3, 0, 8)])
        self.assertEqual([(d.site, d.sample) for d in result.selected],
                         [(1, 0), (2, 0), (3, 0)])
        self.assertEqual(self.state()['attempted'], 24)
        self.assertFalse(self.state()['peptides_per_site'])

    def test_total_budget_is_shared_and_preserves_partial_results(self):
        self.contexts *= 3
        self.config = replace(self.config, num_samples=2, max_attempts=29,
                              peptides_per_site=False)
        self.outcome = lambda site, sample: (
            'valid' if (site, sample) == (2, 0) else 'invalid')
        result = self.controller().run()
        self.assertTrue(result.shortfall)
        self.assertEqual(self.events, [(1, 0, 8), (2, 0, 8), (3, 0, 8),
                                       (1, 8, 5)])
        self.assertEqual([(d.site, d.sample) for d in result.selected],
                         [(2, 0)])
        self.assertEqual(self.state()['attempted'], 29)
        self.config = replace(self.config, resume=True, max_attempts=28)
        with self.assertRaisesRegex(ValueError, 'charged 29 total attempts'):
            self.controller()

    def test_total_mode_keeps_extras_without_starting_another_site(self):
        self.contexts *= 2
        self.config = replace(self.config, num_samples=2,
                              peptides_per_site=False)
        self.outcome = lambda site, sample: 'valid'
        result = self.controller().run()
        self.assertEqual(self.events, [(1, 0, 8)])
        self.assertEqual([d.sample for d in result.selected], [0, 1])
        self.assertEqual(self.state()['valid_count'], 8)
        self.assertEqual(self.state()['selected_count'], 2)
        self.assertEqual(len(find_existing_pepglad_outputs(self.layout)), 8)

    def test_overshoot_keeps_first_valid_structures_and_repeated_sequences(self):
        self.config = replace(self.config, num_samples=2)
        self.outcome = lambda site, sample: 'valid'
        result = self.controller().run()
        self.assertEqual([design.sample for design in result.selected], [0, 1])
        self.assertEqual(self.state()['sites'][0]['valid_count'], 8)
        self.assertEqual(len(find_existing_pepglad_outputs(self.layout)), 8)

    def test_generation_and_validation_failures_consume_attempts(self):
        self.generation_failures = {(1, 0): 'invalid', (1, 1): 'error'}
        self.outcome = lambda site, sample: {2: 'error', 3: 'valid'}.get(sample, 'invalid')
        result = self.controller().run()
        self.assertEqual([design.sample for design in result.selected], [3])
        site = self.state()['sites'][0]
        self.assertEqual(site['attempted'], 8)
        self.assertEqual(site['failures'], dict(generation_invalid=1, relaxation_error=1,
                                               validation_error=1, validation_invalid=4,
                                               unreported=0))
        self.assertEqual(len(list(self.out.glob('**/*.pdb'))), 7)
        self.assertEqual(len(find_existing_pepglad_outputs(self.layout)), 6)
        self.assertTrue(next(self.out.glob('**/invalid_designs/target_1.pdb'))
                        .is_file())
        self.assertTrue(next(self.out.glob('**/target_1.pdb.relax_error.json')).is_file())

    def test_rejected_structures_survive_the_next_site_and_round(self):
        self.contexts *= 2
        self.config = replace(self.config, max_attempts=16)
        self.outcome = lambda site, sample: 'valid' if (site, sample) == (1, 8) else 'invalid'
        saved = {}
        def check_then_generate(*args, **kwargs):
            for path, content in saved.items():
                self.assertEqual(path.read_bytes(), content)
            self.generate(*args, **kwargs)
            for pdb in self.out.glob('**/*.pdb'):
                saved[pdb] = pdb.read_bytes()
        result = self.controller(generate=check_then_generate).run()
        self.assertTrue(result.shortfall)
        self.assertEqual([(d.site, d.sample) for d in result.selected], [(1, 8)])
        self.assertEqual(len(list(self.out.glob('**/*.pdb'))), 32)
        for path, content in saved.items():
            self.assertEqual(path.read_bytes(), content)
        self.assertFalse(list(self.out.glob('**/removed_structures.json')))

    def test_all_invalid_structures_and_reports_survive_resume(self):
        self.config = replace(self.config, max_attempts=8)
        self.controller().run()
        before = {path: path.read_bytes()
                  for path in self.out.glob('**/*.pdb')}
        self.assertEqual(len(before), 8)
        self.assertFalse(list(self.out.glob('**/removed_structures.json')))
        self.assertEqual(len(find_existing_pepglad_outputs(self.layout)), 8)
        self.config = replace(self.config, resume=True)
        self.validation_calls.clear()
        self.assertTrue(self.controller().run().shortfall)
        self.assertEqual(self.validation_calls, [])
        self.assertEqual(self.state()['sites'][0]['attempted'], 8)
        for path, content in before.items():
            self.assertEqual(path.read_bytes(), content)
        with (self.out / 'target_vina_score/target_pepglad_validation.csv').open() as handle:
            self.assertEqual(len(list(csv.DictReader(handle))), 8)

    def test_validated_rounds_resume_without_deleting_any_pdb(self):
        self.outcome = lambda site, sample: 'valid' if sample == 0 else 'invalid'
        self.controller().run()
        self.config = replace(self.config, resume=True)
        self.validation_calls.clear()
        with patch.object(Path, 'unlink') as unlink:
            result = self.controller().run()
        unlink.assert_not_called()
        self.assertEqual(self.validation_calls, [])
        self.assertEqual(self.events, [(1, 0, 8)])
        self.assertEqual(len(list(self.out.glob('**/*.pdb'))), 8)
        _, _, relaxed, pending, invalid = inspect_site(result.selected[0].path.parent)
        self.assertEqual(len(relaxed), 8)
        self.assertEqual(pending, [])
        self.assertEqual(invalid, [])

    def test_missing_valid_pdb_is_not_treated_as_intentionally_removed(self):
        self.outcome = lambda site, sample: 'valid' if sample == 0 else 'invalid'
        result = self.controller().run()
        result.selected[0].path.unlink()
        self.config = replace(self.config, resume=True)
        with self.assertRaisesRegex(ValueError, 'missing=.*target_0'):
            self.controller().run()

    def test_changed_rejected_pdb_is_detected_on_resume(self):
        self.outcome = lambda site, sample: 'valid' if sample == 0 else 'invalid'
        result = self.controller().run()
        pdb = result.selected[0].path.with_name('target_1.pdb')
        pdb.write_text(pdb.read_text() + 'changed\n')
        self.config = replace(self.config, resume=True)
        with self.assertRaisesRegex(ValueError, 'artifact changed'):
            self.controller().run()
        self.assertTrue(pdb.is_file())

    def test_resume_reuses_validated_rounds_and_does_not_reset_budget(self):
        self.config = replace(self.config, max_attempts=8)
        self.assertTrue(self.controller().run().shortfall)
        # Simulate cumulative reports saved before the scorer migration.
        old_reports = self.out / 'target_rosetta_score'
        (self.out / 'target_vina_score').rename(old_reports)
        historical_csv = old_reports / 'target_pepglad_validation.csv'
        historical_bytes = historical_csv.read_bytes()
        self.config = replace(self.config, resume=True)
        self.config = replace(self.config, max_attempts=16)
        self.outcome = lambda site, sample: 'valid' if sample == 8 else 'invalid'
        self.validation_calls.clear()
        result = self.controller().run()
        self.assertEqual(self.events, [(1, 0, 8), (1, 8, 8)])
        self.assertEqual(self.validation_calls, [[(1, sample) for sample in range(8, 16)]])
        self.assertEqual(result.selected[0].sample, 8)
        self.assertEqual(self.state()['sites'][0]['attempted'], 16)
        self.assertEqual(historical_csv.read_bytes(), historical_bytes)
        with (self.out / 'target_vina_score/target_pepglad_validation.csv').open() as handle:
            self.assertEqual(len(list(csv.DictReader(handle))), 16)

    def test_interruption_before_outputs_keeps_reserved_attempts_charged(self):
        def interrupted(*args, **kwargs):
            raise RuntimeError('interrupted generation')
        with self.assertRaisesRegex(RuntimeError, 'interrupted'):
            self.controller(generate=interrupted).run()
        self.assertEqual(self.state()['sites'][0]['attempted'], 8)
        self.config = replace(self.config, resume=True)
        self.config = replace(self.config, max_attempts=16)
        self.outcome = lambda site, sample: 'valid' if sample == 8 else 'invalid'
        result = self.controller().run()
        self.assertEqual(result.selected[0].sample, 8)
        self.assertEqual(self.events, [(1, 8, 8)])
        self.assertEqual(self.state()['sites'][0]['failures']['unreported'], 8)

    def test_partial_generation_can_resume_without_resampling_reserved_ids(self):
        def interrupted(args, pdb, pocket, directory, count, **kwargs):
            self.generate(args, pdb, pocket, directory, 2, **kwargs)
            raise RuntimeError('interrupted generation')
        with self.assertRaises(RuntimeError):
            self.controller(generate=interrupted).run()
        self.config = replace(self.config, resume=True)
        self.outcome = lambda site, sample: 'valid' if sample == 0 else 'invalid'
        result = self.controller().run()
        self.assertEqual(result.selected[0].sample, 0)
        self.assertEqual(self.state()['sites'][0]['attempted'], 8)
        self.assertEqual(self.state()['sites'][0]['failures']['unreported'], 6)

    def test_unrelaxed_outputs_require_recovery_before_resume(self):
        def unrelaxed(*args, **kwargs):
            self.generate(*args, **kwargs)
            (args[3] / 'target_0.pdb').write_text('unrelaxed coordinates\n')
        with self.assertRaisesRegex(RuntimeError, 'allopep.recover_relaxation'):
            self.controller(generate=unrelaxed).run()
        self.config = replace(self.config, resume=True)
        with self.assertRaisesRegex(RuntimeError, 'unrelaxed'):
            self.controller().run()
        self.assertEqual(self.state()['sites'][0]['attempted'], 8)
        self.assertEqual(len(self.events), 1)

    def test_resume_rejects_changed_input_settings_pockets_and_artifacts(self):
        self.outcome = lambda site, sample: 'valid'
        self.controller().run()
        self.config = replace(self.config, resume=True)
        self.config = replace(self.config, peptides_per_site=False)
        with self.assertRaisesRegex(ValueError, 'settings'):
            self.controller()
        self.config = replace(self.config, peptides_per_site=True)
        self.config = replace(self.config, length_max=self.config.length_max + 1)
        with self.assertRaisesRegex(ValueError, 'settings'):
            self.controller()
        self.config = replace(self.config, length_max=self.config.length_max - 1)
        self.pdb.write_text('different receptor')
        with self.assertRaisesRegex(ValueError, 'input'):
            self.controller()
        self.pdb.write_text('input receptor\n')
        pocket = self.out / 'target_allosteric_sites/target_allosteric_site_1.json'
        atomic_json(pocket, [['A', [2, ' ']]])
        with self.assertRaisesRegex(ValueError, 'pocket differs'):
            self.controller()
        atomic_json(pocket, [['A', [1, ' ']]])
        model = next(self.out.glob('**/target_0.pdb'))
        model.write_text(model.read_text() + 'changed\n')
        with self.assertRaisesRegex(ValueError, 'artifact changed'):
            self.controller().run()

    def test_fresh_run_and_smaller_resume_cap_are_rejected(self):
        self.config = replace(self.config, max_attempts=8)
        self.controller().run()
        with self.assertRaisesRegex(ValueError, 'require --resume'):
            self.controller()
        self.config = replace(self.config, resume=True)
        self.config = replace(self.config, max_attempts=7)
        with self.assertRaisesRegex(ValueError, 'already charged 8'):
            self.controller()

    def test_pre_round_site_layout_is_rejected(self):
        """Designs directly in a site directory are the removed legacy shape."""
        self.config = replace(self.config, resume=True)
        directory = self.out / 'target_PepGLAD_outputs/target_allosteric_site_1_codesign'
        directory.mkdir(parents=True)
        atomic_json(directory / 'round.json', {'site': 1, 'start': 0, 'count': 8})
        self.generate(self.config, str(self.pdb), None, directory, 8,
                      deadline=time.monotonic() + 60)
        (directory / 'round.json').unlink()
        atomic_json(self.out / 'target_allosteric_sites/target_allosteric_site_1.json',
                    [['A', [1, ' ']]])
        with self.assertRaisesRegex(ValueError, 'round_NNNNNN subdirectories'):
            self.controller()

    def test_timeout_preserves_budget_and_error_summary(self):
        def timeout(*args, **kwargs):
            raise TimeoutError('deadline expired')
        with self.assertRaises(TimeoutError):
            self.controller(validate=timeout).run()
        state = self.state()
        self.assertEqual(state['status'], 'error')
        self.assertEqual(state['sites'][0]['attempted'], 8)
        self.assertIn('TimeoutError', state['error'])
        self.assertEqual(len(list(self.out.glob('**/*.pdb'))), 8)
        self.assertFalse(list(self.out.glob('**/removed_structures.json')))

    def test_invalid_validation_status_is_an_operational_failure(self):
        self.outcome = lambda site, sample: 'unknown'
        with self.assertRaisesRegex(ValueError, 'validation status'):
            self.controller().run()
        self.assertEqual(self.state()['status'], 'error')
        self.assertEqual(self.state()['sites'][0]['attempted'], 8)

    def test_missing_or_duplicate_sample_metadata_is_rejected(self):
        def corrupt(args, pdb, pocket, directory, count, **kwargs):
            self.generate(args, pdb, pocket, directory, count, **kwargs)
            summary = directory / 'summary.jsonl'
            summary.write_text(summary.read_text() + summary.read_text().splitlines()[0] + '\n')
        with self.assertRaisesRegex(ValueError, 'Duplicate PepGLAD IDs'):
            self.controller(generate=corrupt).run()
        self.assertEqual(self.state()['sites'][0]['attempted'], 8)

    def test_sampling_lock_prevents_concurrent_budget_use(self):
        with sampling_lock(self.layout):
            with self.assertRaisesRegex(RuntimeError, 'Another sampling'):
                with sampling_lock(self.layout):
                    self.fail('second lock acquired')

    def test_main_returns_shortfall_and_ranks_only_selected_designs(self):
        old_dir = os.getcwd()
        self.addCleanup(os.chdir, old_dir)
        os.chdir(self.root)
        selected = [SimpleNamespace(sample=0)]
        for shortfall, expected, cap in ((False, 0, None), (True, 2, None),
                                         (False, 0, 20), (False, 0, 50000)):
            result = SimpleNamespace(selected=selected, shortfall=shortfall, summary_path='summary')
            with patch('allopep.cli.check_vina_dependencies', return_value='/fake/ia'), \
                 patch('allopep.cli.prepare_input_pdb', return_value=str(self.pdb)), \
                 patch('allopep.cli.gaps.run_gaps'), \
                 patch('allopep.cli.gaps.select_guidance_contexts', return_value=self.contexts), \
                 patch('allopep.cli.run_pepglad', return_value=result) as sample, \
                 patch('allopep.cli.write_pymol_session', return_value='session') as pymol, \
                 patch('allopep.cli.rank_pepglad_outputs', return_value=[{}]) as rank, \
                 patch('allopep.cli.export_peptides'), \
                 patch('allopep.cli.run_openstructure_validation') as validate:
                arguments = ['target.pdb', '-n', '3']
                if cap is not None:
                    arguments.extend(['--max_num_peptides', str(cap)])
                self.assertEqual(main(arguments), expected)
                self.assertEqual(sample.call_args.args[0].max_attempts,
                                 3000 if cap is None else cap)
                self.assertEqual(rank.call_args.args[0], selected)
                self.assertEqual(pymol.call_args.args[0], selected)
                validate.assert_not_called()

    def test_existing_modes_validate_without_generation(self):
        old_dir = os.getcwd()
        self.addCleanup(os.chdir, old_dir)
        os.chdir(self.root)
        designs = [SimpleNamespace(site=site, sample=sample)
                   for site in (1, 2) for sample in range(3)]
        cases = [(option, per_site) for option in (
            '--validate-existing', '--rank-existing')
                 for per_site in (False, True)]
        for option, per_site in cases:
            with patch('allopep.cli.check_vina_dependencies', return_value='/fake/ia'), \
                 patch('allopep.cli.prepare_input_pdb', return_value=str(self.pdb)), \
                 patch('allopep.cli.find_existing_pepglad_outputs', return_value=designs), \
                 patch('allopep.cli.write_allosteric_site_visualization',
                       return_value=Path('target_allosteric_site.pse')) as site_view, \
                 patch('allopep.cli.run_openstructure_validation', return_value=designs), \
                 patch('allopep.cli.write_pymol_session', return_value='session'), \
                 patch('allopep.cli.rank_pepglad_outputs', return_value=[{}]) as rank, \
                 patch('allopep.cli.export_peptides') as export, \
                 patch('allopep.cli.run_pepglad') as sample, \
                 patch('allopep.cli.gaps.run_gaps') as gaps:
                arguments = ['target.pdb', option, '-n', '2']
                if per_site:
                    arguments.append('--peptides-per-site')
                self.assertEqual(main(arguments), 0)
                sample.assert_not_called()
                gaps.assert_not_called()
                site_view.assert_called_once_with(
                    OutputLayout.at(Path('output/target_output'), 'target'))
                self.assertEqual(rank.call_count, int(option == '--rank-existing'))
                expected = [designs[0], designs[1]]
                if per_site:
                    expected.extend(designs[3:5])
                self.assertEqual(export.call_args.args[0], expected)
                if option == '--rank-existing':
                    self.assertEqual(rank.call_args.args[0], expected)


if __name__ == '__main__':
    unittest.main()
