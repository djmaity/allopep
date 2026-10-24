"""Allosteric site views must be available independently of peptide success."""

import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from allopep.cli import run_pepglad
from allopep.layout import OutputLayout
from allopep.pymol_session import (_pymol_commands,
                                   write_allosteric_site_visualization)


class AllostericVisualizationTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.out = Path(self.temporary.name) / 'target output'
        self.out.mkdir()
        self.pdb = self.out / 'target.pdb'
        self.scores = self.out / 'target_GAPS_scores.json'
        identities = [('A', 10, ' '), ('A', 11, 'A'), ('B', 20, ' '),
                      ('B', 12, ' '), (' ', 30, ' ')]
        self.records = []
        lines = []
        for index, ((chain, number, icode), score) in enumerate(zip(
                identities, (0.123456789, 0.8, 0.9, 0.5, 0.6)), 1):
            x = index * 4.0
            lines.append(
                f'ATOM  {index:5d}  CA  ALA {chain}{number:4d}{icode}   '
                f'{x:8.3f}{0:8.3f}{0:8.3f}{1:6.2f}{25:6.2f}           C  ')
            self.records.append({'chain': chain, 'number': number, 'icode': icode,
                                 'resname': 'ALA', 'score': score, 'ca': [x, 0.0, 0.0]})
        self.pdb.write_text('\n'.join(lines) + '\nEND\n')
        self.scores.write_text(json.dumps({'schema_version': 1, 'residues': self.records}))
        (self.out / 'target_GAPS_output.pdb').write_text('unused renumbered export\n')
        self.sites = self.out / 'target_allosteric_sites'
        self.sites.mkdir()
        self.layout = OutputLayout.at(self.out, 'target')

    def pocket(self, site, residues):
        (self.sites / f'target_allosteric_site_{site}.json').write_text(
            json.dumps(residues), encoding='utf-8')

    def test_script_and_session_include_all_sites_without_peptides(self):
        self.pocket(1, [['A', [10, ' ']], ['A', [11, 'A']], ['B', [20, ' ']]])
        self.pocket(2, [['B', [12, ' ']]])
        self.pocket(10, [[' ', [30, ' ']]])
        (self.sites / 'unrelated.json').write_text('null')

        def build(command, **kwargs):
            script = (kwargs['cwd'] / command[2]).read_text()
            self.assertEqual(script.splitlines()[-1], 'save session.pse')
            (kwargs['cwd'] / 'session.pse').write_bytes(b'new session')
            return SimpleNamespace(returncode=0, stdout='', stderr='')

        with patch('allopep.pymol_session.subprocess.run', side_effect=build):
            session = write_allosteric_site_visualization(self.layout)
        self.assertEqual(session, self.out / 'target_allosteric_site.pse')
        self.assertEqual(session.read_bytes(), b'new session')
        content = session.with_suffix('.pml').read_text()
        self.assertIn('spectrum b, blue_white_red, target and polymer.protein and b > -0.5, '
                      'minimum=0, maximum=1', content)
        self.assertIn('show surface, site_1 and ((chain A and resi 10+11A) '
                      'or (chain B and resi 20))', content)
        self.assertIn('show surface, site_2 and ((chain B and resi 12))', content)
        self.assertIn('show surface, site_10 and ((chain "" and resi 30))', content)
        self.assertLess(content.index('group sites, site_2'),
                        content.index('group sites, site_10'))
        self.assertNotIn('peptide', content)
        self.assertNotIn('\npython\n', content)
        self.assertNotIn('\nsave ', content)
        self.assertFalse(list(self.out.glob('.allopep-pymol-*')))

    def test_pymol_error_preserves_both_previous_files(self):
        self.pocket(1, [['A', [10, ' ']]])
        session = self.out / 'target_allosteric_site.pse'
        script = session.with_suffix('.pml')
        session.write_bytes(b'previous session')
        script.write_text('previous script')

        def fail(command, **kwargs):
            (kwargs['cwd'] / 'session.pse').write_bytes(b'broken session')
            return SimpleNamespace(returncode=0, stdout='Error: invalid selection',
                                   stderr='')

        with patch('allopep.pymol_session.subprocess.run', side_effect=fail):
            with self.assertRaisesRegex(RuntimeError, 'invalid selection'):
                write_allosteric_site_visualization(self.layout)
        self.assertEqual(session.read_bytes(), b'previous session')
        self.assertEqual(script.read_text(), 'previous script')
        self.assertFalse(list(self.out.glob('.allopep-pymol-*')))

    def test_full_precision_scores_and_sites_use_original_residue_ids(self):
        original = self.pdb.read_text()
        self.pocket(1, [['A', [10, ' ']], ['A', [11, 'A']], ['B', [20, ' ']]])
        content = '\n'.join(_pymol_commands([], self.layout, include_all_sites=True))
        self.assertIn(f'load {json.dumps(str(self.pdb.resolve()))}, target', content)
        self.assertNotIn('target_GAPS_output.pdb', content)
        self.assertIn('alter target and polymer.protein and chain "A" and resi 10, '
                      'b=0.123456789', content)
        self.assertIn('alter target and polymer.protein and chain "A" and resi 11A, '
                      'b=0.8', content)
        self.assertIn('alter target and polymer.protein and chain "B" and resi 20, '
                      'b=0.9', content)
        self.assertIn('show surface, site_1 and ((chain A and resi 10+11A) '
                      'or (chain B and resi 20))', content)
        self.assertEqual(self.pdb.read_text(), original)

    def test_scores_from_different_coordinates_are_rejected(self):
        self.records[1]['ca'][0] += 1.0
        self.scores.write_text(json.dumps({'schema_version': 1, 'residues': self.records}))
        with self.assertRaisesRegex(ValueError, 'coordinates differ from original PDB'):
            _pymol_commands([], self.layout, include_all_sites=True)

    def test_missing_original_pdb_does_not_use_renumbered_export(self):
        self.pdb.unlink()
        with self.assertRaisesRegex(FileNotFoundError, 'Original receptor structure'):
            write_allosteric_site_visualization(self.layout)

    def test_missing_scores_do_not_use_original_temperature_factors(self):
        self.scores.unlink()
        with self.assertRaises(FileNotFoundError):
            write_allosteric_site_visualization(self.layout)

    def test_native_pymol_session_preserves_ids_scores_and_unscored_color(self):
        import pymol2

        original = self.pdb.read_text().replace(
            'END\n',
            'ATOM      6  CA  ALA C  99      24.000   0.000   0.000  1.00 25.00           C  \nEND\n')
        self.pdb.write_text(original)
        self.pocket(1, [['A', [11, 'A']], ['B', [20, ' ']]])
        session = write_allosteric_site_visualization(self.layout)
        with pymol2.PyMOL() as pymol:
            cmd = pymol.cmd
            cmd.load(str(session))
            atoms = cmd.get_model('target').atom
            actual = {(atom.chain or ' ', atom.resi): atom.b for atom in atoms}
            for record in self.records:
                identity = record['chain'], f"{record['number']}{record['icode'].strip()}"
                self.assertAlmostEqual(actual[identity], record['score'], places=7)
            self.assertEqual(actual[('C', '99')], -1)
            colors = []
            cmd.iterate('target and chain C', 'colors.append(color)', space={'colors': colors})
            self.assertEqual(colors, [cmd.get_color_index('gray70')])
            surface = {(atom.chain, atom.resi)
                       for atom in cmd.get_model('site_1 and rep surface').atom}
            self.assertEqual(surface, {('A', '11A'), ('B', '20')})
        self.assertEqual(self.pdb.read_text(), original)

    def test_missing_sites_are_reported_before_launching_pymol(self):
        with patch('allopep.pymol_session.subprocess.run') as launch:
            with self.assertRaisesRegex(FileNotFoundError, 'No allosteric site'):
                write_allosteric_site_visualization(self.layout)
        launch.assert_not_called()

    def test_site_view_is_created_before_sampling_starts(self):
        events = []
        result = object()

        def start_sampling():
            self.assertTrue((self.out / 'target_allosteric_site.pse').is_file())
            events.append('sampling')
            return result

        def create_view(layout):
            session = layout.root / f'{layout.prefix}_allosteric_site.pse'
            session.touch()
            events.append('visualization')
            return session

        with patch('allopep.cli.SamplingController') as controller, \
             patch('allopep.cli.write_allosteric_site_visualization',
                   side_effect=create_view):
            controller.return_value.run.side_effect = start_sampling
            self.assertIs(run_pepglad(None, 'target.pdb', self.layout, []),
                          result)
        self.assertEqual(events, ['visualization', 'sampling'])


if __name__ == '__main__':
    unittest.main()
