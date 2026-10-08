"""Regression test for atom-pair clash counting."""

import unittest
from pathlib import Path

from count_clashes import count_clashes


class ClashCountingTests(unittest.TestCase):
    def test_counts_each_clash_once(self):
        pdb_path = Path(__file__).parent / "fixtures" / "sample_GAPS_output.pdb"
        self.assertEqual(count_clashes(pdb_path), 2)


if __name__ == "__main__":
    unittest.main()
