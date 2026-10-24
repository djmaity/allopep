"""Check selected-complex exports and their joined pocket/score table."""

import csv
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from allopep.cli import main
from allopep.fileio import atomic_csv, atomic_json
from allopep.layout import OutputLayout
from allopep.outputs import DesignOutput
from allopep.peptide_exports import FIELDS, export_peptides
from allopep.scoring.models import WorkflowConfig
from allopep.scoring.records import failure


class PeptideExportTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.layout = OutputLayout.at(self.root, 'target')
        fixture = Path(__file__).parent / 'fixtures/sample_GAPS_output.pdb'
        lines = fixture.read_text().splitlines()
        identities = [('A', 'ARG', 42, ' '), ('A', 'GLY', 43, ' '),
                      ('B', 'TYR', 105, 'A')]
        records = []
        for index, (chain, name, number, icode) in enumerate(identities):
            line = lines[index]
            lines[index] = (line[:17] + name + line[20:21] + chain
                            + f'{number:4d}' + icode + line[27:])
            records.append({
                'chain': chain, 'resname': name, 'number': number,
                'icode': icode, 'score': [0.873241, 0.772346, 0.649125][index],
                'ca': [index * 4.0, 0.0, 0.0],
            })
        self.layout.input_pdb.write_text('\n'.join(lines) + '\n')
        atomic_json(self.layout.gaps_scores,
                    {'schema_version': 1, 'residues': records})
        atomic_json(self.layout.pocket_file(1),
                    [['A', [42, ' ']], ['A', [43, ' ']]])
        atomic_json(self.layout.pocket_file(2),
                    [['A', [43, ' ']], ['B', [105, 'A']]])
        areas = {'A': {'42': SimpleNamespace(total=12.5),
                       '43': SimpleNamespace(total=0.25)},
                 'B': {'105A': SimpleNamespace(total=7.0)}}
        sasa = SimpleNamespace(residueAreas=lambda: areas)
        mock = patch('allopep.gaps.freesasa.calcBioPDB',
                     return_value=(sasa, {}))
        self.sasa = mock.start()
        self.addCleanup(mock.stop)
        self.designs = []
        for site, sample, round_index in ((1, 0, 1), (1, 8, 2), (2, 0, 1)):
            path = self.layout.round_dir(site, round_index) / 'target_0.pdb'
            path.parent.mkdir(parents=True, exist_ok=True)
            peptide = lines[0][:21] + 'C' + lines[0][22:]
            path.write_text('\n'.join([*lines, peptide]) + '\n')
            self.designs.append(DesignOutput(
                path, site, sample, 'LHLKA', 'C', ('A', 'B')))
        self.rejected = self.designs[0].path.with_name('target_1.pdb')
        self.rejected.write_text('rejected complex retained\n')
        validation = [{
            'path': d.resolved_path, 'status': 'valid',
            'buried_sasa_a2': 150.5, 'receptor_contact_residues_5a': 4,
            'minimum_interface_distance_a': 2.5,
        } for d in self.designs]
        validation.append({'path': str(self.rejected.resolve()),
                           'status': 'invalid'})
        atomic_csv(self.layout.validation_reports(self.layout.score_dir)[0],
                   tuple(validation[0]), validation)

    def rows(self):
        with self.layout.peptides_csv.open(newline='') as handle:
            return list(csv.DictReader(handle))

    def test_copies_complexes_and_reports_residue_names_scores_and_sasa(self):
        export_peptides(list(reversed(self.designs)), self.layout,
                        self.layout.input_pdb)
        rows = self.rows()
        self.assertEqual([(r['site'], r['sample']) for r in rows],
                         [('1', '0'), ('1', '8'), ('2', '0')])
        self.assertEqual(rows[0]['site_residues'], 'A-ARG-42;A-GLY-43')
        self.assertEqual(rows[2]['site_residues'], 'A-GLY-43;B-TYR-105A')
        self.assertEqual(float(rows[0]['pocket_sasa_a2']), 12.75)
        self.assertEqual(float(rows[2]['pocket_sasa_a2']), 7.25)
        self.assertAlmostEqual(float(rows[0]['gaps_pocket_score']),
                               (0.873241 + 0.772346) / 2)
        self.assertEqual(rows[0]['peptide_sequence'], 'LHLKA')
        self.assertEqual(rows[0]['peptide_length'], '5')
        self.assertEqual(rows[0]['buried_sasa_a2'], '150.5')
        self.assertEqual(rows[0]['vina_status'], 'not_run')
        self.assertEqual(rows[0]['vina_score'], '')
        for row, design in zip(rows, self.designs):
            copied = self.root / row['complex_path']
            self.assertEqual(copied.name,
                             f'target_site_{design.site}_'
                             f'peptide_{design.sample}.pdb')
            self.assertEqual(copied.read_bytes(), design.path.read_bytes())
            self.assertEqual(self.root / row['source_path'], design.path)
            self.assertFalse(Path(row['complex_path']).is_absolute())
        self.assertTrue(self.rejected.is_file())

    def test_refresh_only_copies_and_keep_header_for_empty_selection(self):
        export_peptides(self.designs, self.layout, self.layout.input_pdb)
        export_peptides(self.designs[2:], self.layout, self.layout.input_pdb)
        self.assertEqual(list(self.layout.valid_complexes_dir.glob(
            'site_*/*.pdb')), [self.layout.valid_complex_file(2, 0)])
        export_peptides([], self.layout, self.layout.input_pdb)
        self.assertEqual(self.rows(), [])
        with self.layout.peptides_csv.open(newline='') as handle:
            self.assertEqual(csv.DictReader(handle).fieldnames, list(FIELDS))
        self.assertFalse(list(self.layout.valid_complexes_dir.glob(
            '**/*.pdb')))
        self.assertTrue(all(d.path.is_file() for d in self.designs))
        self.assertTrue(self.rejected.is_file())

    def test_csv_joins_scores_and_keeps_failed_or_missing_outcomes(self):
        scores = [
            {'id': self.designs[0].job_id, 'status': 'ok', 'score_total': -7.5,
             'score_units': 'kcal/mol', 'protocol_id': 'test', 'error': ''},
            {'id': self.designs[1].job_id, 'status': 'error',
             'score_total': 'nan', 'error': 'bad residue'},
        ]
        atomic_csv(self.layout.vina_results, tuple(scores[0]), scores)
        atomic_csv(self.layout.vina_ranked, ('id', 'rank'),
                   [{'id': self.designs[0].job_id, 'rank': 1},
                    {'id': self.designs[1].job_id, 'rank': 2}])
        export_peptides(self.designs, self.layout, self.layout.input_pdb,
                        scored=True, scoring_error='worker stopped')
        first, failed, missing = self.rows()
        self.assertEqual((first['vina_score'], first['vina_rank']),
                         ('-7.5', '1'))
        self.assertEqual(first['vina_units'], 'kcal/mol')
        self.assertEqual(failed['vina_status'], 'error')
        self.assertEqual(failed['vina_error'], 'bad residue')
        self.assertEqual((failed['vina_score'], failed['vina_rank']), ('', ''))
        self.assertEqual(missing['vina_error'], 'worker stopped')
        self.assertEqual(len(list(self.layout.valid_complexes_dir.glob(
            '**/*.pdb'))), 3)

    def existing_main(self, arguments):
        """Run the real selection/export flow around saved validation data."""
        with patch('allopep.cli.OutputLayout.for_input',
                   return_value=self.layout), \
             patch('allopep.cli.find_existing_pepglad_outputs',
                   return_value=self.designs), \
             patch('allopep.cli.run_openstructure_validation',
                   return_value=self.designs), \
             patch('allopep.cli.write_allosteric_site_visualization',
                   return_value=self.root / 'sites.pse'), \
             patch('allopep.cli.write_pymol_session',
                   return_value='session'), \
             patch('allopep.cli.check_vina_dependencies'):
            return main([str(self.layout.input_pdb), *arguments])

    def test_existing_exports_partial_results_and_skips_stale_scores(self):
        atomic_csv(self.layout.vina_results, ('id', 'status', 'score_total'),
                   [{'id': self.designs[0].job_id, 'status': 'ok',
                     'score_total': -99}])
        self.assertEqual(self.existing_main([
            '--validate-existing', '-n', '2']), 0)
        self.assertEqual(len(self.rows()), 2)
        self.assertEqual(self.rows()[0]['vina_status'], 'not_run')
        self.assertEqual(self.existing_main([
            '--validate-existing', '-n', '2', '--peptides-per-site']), 2)
        self.assertEqual(len(self.rows()), 3)

    def test_all_vina_failures_still_export_selected_complexes_and_csv(self):
        def fail_jobs(jobs):
            return [failure(job, WorkflowConfig(), ValueError('bad residue'))
                    for job in jobs]

        with patch('allopep.scoring.service.BatchRunner.run',
                   side_effect=fail_jobs), \
             self.assertRaisesRegex(RuntimeError, 'All 2 Vina scoring jobs '
                                    'failed'):
            self.existing_main(['--rank-existing', '-n', '2'])
        rows = self.rows()
        self.assertEqual(len(rows), 2)
        self.assertTrue(all(r['vina_status'] == 'error' for r in rows))
        self.assertTrue(all('bad residue' in r['vina_error'] for r in rows))
        self.assertTrue(all(not r['vina_score'] and not r['vina_rank']
                            for r in rows))
        self.assertEqual(len(list(self.layout.valid_complexes_dir.glob(
            '**/*.pdb'))), 2)
        self.assertTrue(self.rejected.is_file())


if __name__ == '__main__':
    unittest.main()
