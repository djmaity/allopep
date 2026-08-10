"""Tests for spatial selection of GAPS-scored residues."""

import unittest
from pathlib import Path

from select_gaps_residues import (
    find_high_bfactor_spatial_patches,
    select_contiguous_high_bfactor_residues,
)


class SpatialResidueSelectionTests(unittest.TestCase):
    def setUp(self):
        self.pdb_path = (
            Path(__file__).parent / "fixtures" / "sample_GAPS_output.pdb"
        )

    def test_finds_cross_chain_connected_patch(self):
        patches = find_high_bfactor_spatial_patches(
            self.pdb_path, bfactor_threshold=0.5,
            distance_cutoff=5.0, min_patch_size=2,
        )

        self.assertEqual(len(patches), 1)
        self.assertEqual(
            patches[0].residues,
            (("A", 1, " "), ("A", 2, " "),
             ("A", 3, " "), ("B", 1, " ")),
        )
        self.assertAlmostEqual(patches[0].mean_bfactor, 0.75)

    def test_selector_returns_pepglad_format(self):
        selected = select_contiguous_high_bfactor_residues(
            self.pdb_path, bfactor_threshold=0.5,
            distance_cutoff=5.0, min_patch_size=2,
        )

        self.assertEqual(selected, [
            ["A", [1, " "]], ["A", [2, " "]],
            ["A", [3, " "]], ["B", [1, " "]],
        ])

    def test_returns_empty_list_when_no_residue_passes(self):
        selected = select_contiguous_high_bfactor_residues(
            self.pdb_path, bfactor_threshold=1.0,
        )

        self.assertEqual(selected, [])

    def test_rejects_invalid_parameters(self):
        with self.assertRaises(ValueError):
            find_high_bfactor_spatial_patches(
                self.pdb_path, distance_cutoff=0,
            )


if __name__ == "__main__":
    unittest.main()
