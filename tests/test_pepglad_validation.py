"""Checks for the physical-validity policy and its report formats."""

import csv
import json
import tempfile
import unittest
from pathlib import Path

from pepglad_validation import FIELDS, validation_reasons, write_reports


class ValidationPolicyTests(unittest.TestCase):
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
