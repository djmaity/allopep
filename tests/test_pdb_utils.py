"""Tests for the shared PDB chain splitter."""

import tempfile
import unittest
from pathlib import Path

from Bio.PDB import PDBParser

from pdb_utils import split_protein_peptide


class PDBUtilsTests(unittest.TestCase):
    def setUp(self):
        self.fixture = (
            Path(__file__).parent / 'fixtures' / 'sample_GAPS_output.pdb'
        )

    def test_shortest_chain_is_saved_as_peptide(self):
        with tempfile.TemporaryDirectory() as directory:
            prefix = Path(directory) / 'sample'
            structure = split_protein_peptide(self.fixture, output_prefix=prefix)
            parser = PDBParser(QUIET=True)
            peptide = parser.get_structure('peptide', f'{prefix}_peptide.pdb')
            protein = parser.get_structure('protein', f'{prefix}_protein.pdb')

            self.assertEqual([chain.id for chain in structure[0]], ['A', 'B'])
            self.assertEqual([chain.id for chain in peptide[0]], ['B'])
            self.assertEqual([chain.id for chain in protein[0]], ['A'])

    def test_known_peptide_chain_overrides_detection(self):
        with tempfile.TemporaryDirectory() as directory:
            prefix = Path(directory) / 'sample'
            split_protein_peptide(
                self.fixture, peptide_chain='A', output_prefix=prefix
            )
            peptide = PDBParser(QUIET=True).get_structure(
                'peptide', f'{prefix}_peptide.pdb'
            )
            self.assertEqual([chain.id for chain in peptide[0]], ['A'])


if __name__ == '__main__':
    unittest.main()
