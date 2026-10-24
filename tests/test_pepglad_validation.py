"""Checks for the physical-validity policy and its report formats."""

import csv
import json
import os
import tempfile
import time
import unittest
from pathlib import Path
from unittest.mock import patch

from allopep.validation import (FIELDS, validate_manifest, validation_reasons,
                                write_reports)


def _probe_validation(entry, peptide_chain):
    time.sleep(0.05)
    if entry.get('fail'):
        raise ValueError('sample failure')
    return {'path': entry['path'], 'status': 'valid', 'pid': os.getpid()}, {}


class ValidationPolicyTests(unittest.TestCase):
    def test_parallel_validation_runs_in_multiple_processes_and_keeps_order(self):
        entries = [{'path': f'sample_{index}.pdb', 'site': 1, 'sample': index,
                    'fail': index == 2} for index in range(4)]
        with patch('allopep.validation.validate_structure', _probe_validation):
            reports = list(validate_manifest(entries, workers=2))
        rows = [row for row, _ in reports]
        self.assertEqual([row['path'] for row in rows],
                         [entry['path'] for entry in entries])
        self.assertEqual([row['status'] for row in rows],
                         ['valid', 'valid', 'error', 'valid'])
        self.assertIn('ValueError: sample failure', rows[2]['reasons'])
        self.assertEqual(len({row['pid'] for row in rows if row['status'] == 'valid'}), 2)

    def test_rejection_reasons_preserve_policy_and_order(self):
        metrics = {field: 0 for field in FIELDS}
        metrics['contact_residues_5a'] = 1
        self.assertEqual(validation_reasons(metrics), [])
        metrics['inverted_ca'] = 2
        metrics['rama_outliers'] = 1
        metrics['contact_residues_5a'] = 0
        self.assertEqual(validation_reasons(metrics), [
            'inverted_ca', 'rama_outliers', 'no_receptor_contact_5a'])
        metrics['receptor_clashes'] = 3
        metrics['cis_peptide_bonds'] = 1
        self.assertEqual(validation_reasons(metrics), [
            'inverted_ca', 'rama_outliers', 'no_receptor_contact_5a'])

    def test_report_files_keep_csv_columns_and_jsonl_details(self):
        reports = [({'path': 'sample.pdb', 'site': 1, 'sample': 0,
                     'status': 'invalid', 'reasons': 'inverted_ca'},
                    {'inverted_ca': [{'residue': 'B.1.ALA'}]})]
        with tempfile.TemporaryDirectory() as directory:
            csv_path = Path(directory) / 'validation.csv'
            jsonl_path = Path(directory) / 'validation.jsonl'
            write_reports(reports, csv_path, jsonl_path)
            with csv_path.open(newline='', encoding='utf-8') as handle:
                reader = csv.DictReader(handle)
                self.assertEqual(tuple(reader.fieldnames), FIELDS)
                self.assertEqual(next(reader)['reasons'], 'inverted_ca')
            detail = json.loads(jsonl_path.read_text(encoding='utf-8'))
            self.assertEqual(detail['details']['inverted_ca'][0]['residue'], 'B.1.ALA')


if __name__ == '__main__':
    unittest.main()
