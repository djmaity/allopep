import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import ModuleType
from typing import ClassVar
from unittest import mock

import batch
import preparation
import scoring
import vina_score
from models import ScoreJob, WorkflowConfig, group_jobs


def pdbqt_atom_line(serial: int, atom: str, residue: str, chain: str,
    residue_number: int, x: float, y: float, z: float, atom_type: str,) -> str:
    return (
        f"ATOM  {serial:5d} {atom:<4} {residue:>3} "
        f"{chain}{residue_number:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00    -0.100 {atom_type:>2}\n")


def pdb_atom_line(serial: int, chain: str, x: float, y: float, z: float,
    element: str = "C", *, record: str = "ATOM", residue: str = "ALA",) -> str:
    return (
        f"{record:<6}{serial:5d}  CA  {residue:>3} {chain}{serial:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00          {element:>2}\n")


class ScoringTests(unittest.TestCase):
    def test_box_padding_is_applied_to_each_side(self):
        center, size = scoring.box_from_coordinates(
            [[(0.0, 2.0, 4.0), (10.0, 6.0, 8.0)]], padding=5.0)
        self.assertEqual(center, [5.0, 4.0, 6.0])
        self.assertEqual(size, [20.0, 14.0, 14.0])

    def test_box_validation_rejects_atom_outside(self):
        with self.assertRaisesRegex(ValueError, "outside"):
            scoring.validate_box([0, 0, 0], [2, 2, 2], [[(2, 0, 0)]])

    def test_named_energy_components_and_normalization(self):
        values = scoring.named_energies(range(8))
        self.assertEqual(values["total"], 0.0)
        self.assertEqual(values["torsional"], 6.0)
        normalized = scoring.size_normalized_diagnostics(
            {"total": -10.0, "ligand_inter": -12.0},
            heavy_atoms=5, residues=2,)
        self.assertEqual(normalized["total_per_heavy_atom"], -2.0)
        self.assertEqual(normalized["interaction_per_residue"], -6.0)

    def test_parse_pdbqt_metrics(self):
        content = "".join(
            ["ROOT\n", pdbqt_atom_line(1, "C", "ALA", "B", 1, 1, 2, 3, "C"),
                pdbqt_atom_line(2, "H", "ALA", "B", 1, 2, 3, 4, "HD"),
                pdbqt_atom_line(3, "N", "GLY", "B", 2, 3, 4, 5, "N"),
                "ENDROOT\n", "TORSDOF 0\n",])
        with tempfile.TemporaryDirectory() as temporary:
            filename = Path(temporary) / "pose.pdbqt"
            filename.write_text(content)
            metrics = scoring.parse_pdbqt(filename)
        self.assertEqual(metrics["atom_count"], 3)
        self.assertEqual(metrics["heavy_atom_count"], 2)
        self.assertEqual(metrics["residue_count"], 2)
        self.assertEqual(metrics["active_torsions"], 0)

    def test_scorer_sums_raw_residue_scores_after_computing_maps(self):
        fake_vina_module = ModuleType("vina")

        class FakeVina:
            instances: ClassVar[list["FakeVina"]] = []

            def __init__(self, **kwargs):
                self.events = []
                self.ligands = []
                self._vina = self
                self.__class__.instances.append(self)

            def set_receptor(self, receptor):
                self.receptor = receptor

            def compute_vina_maps(self, **kwargs):
                self.events.append("maps")

            def set_ligand_from_string(self, ligand):
                self.events.append("ligand")
                self.ligands.append(ligand)

            def score(self):
                self.events.append("score")
                value = -0.4444 if len(self.ligands) == 1 else -0.5555
                return [value, value, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

        fake_vina_module.Vina = FakeVina
        receptor_pose = (
            "ROOT\n" + pdbqt_atom_line(1, "C", "ALA", "B", 1, 0, 0, 0, "C")
            + "ENDROOT\nTORSDOF 0\n")
        peptide_pose = (
            "ROOT\n"
            + pdbqt_atom_line(1, "C", "ALA", "B", 1, 0, 0, 0, "C")
            + pdbqt_atom_line(2, "N", "GLY", "B", 2, 1, 0, 0, "N")
            + "ENDROOT\nTORSDOF 0\n")
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            receptor = directory / "receptor.pdbqt"
            peptide = directory / "peptide.pdbqt"
            receptor.write_text(receptor_pose)
            peptide.write_text(peptide_pose)
            metrics = scoring.parse_pdbqt(peptide)
            with mock.patch.dict(sys.modules, {"vina": fake_vina_module}):
                scorer = scoring.VinaPoseScorer(
                    receptor, center=[0, 0, 0], box_size=[10, 10, 10])
                result = scorer.score_pose(peptide, metrics)

        instance = FakeVina.instances[0]
        self.assertEqual(instance.events,
            ["maps", "ligand", "score", "ligand", "score"])
        self.assertEqual(result["scores"]["total"], -1.0)
        self.assertEqual(len(instance.ligands), 2)
        self.assertNotIn("GLY", instance.ligands[0])
        self.assertNotIn("ALA", instance.ligands[1])


class PreparationTests(unittest.TestCase):
    def test_combined_pdb_chain_selection(self):
        content = (
            pdb_atom_line(1, "A", 0, 0, 0)
            + pdb_atom_line(2, "B", 5, 6, 7) + "END\n")
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "complex.pdb"
            source.write_text(content)
            preparer = preparation.PdbPreparer(hydrogen_policy="keep")
            coordinates = preparer.selected_coordinates(source, ("B",))
        self.assertEqual(coordinates, [(5.0, 6.0, 7.0)])

    def test_meeko_creates_a_rigid_ligand_record_without_docking(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = directory / "peptide.pdb"
            destination = directory / "peptide.pdbqt"
            source.write_text(pdb_atom_line(1, "B", 0, 0, 0))

            def fake_run(command, **kwargs):
                destination.write_text("REMARK generated\n"
                    + pdbqt_atom_line(7, "C", "ALA", "B", 1, 0, 0, 0, "C"))
                return subprocess.CompletedProcess(command, 0, "", "")

            preparer = preparation.PdbPreparer(hydrogen_policy="rebuild")
            with (
                mock.patch.object(preparation, "find_meeko_receptor_preparer",
                    return_value=Path("mk_prepare_receptor.py"),
                ), mock.patch.object(
                    preparation.subprocess, "run", side_effect=fake_run
                ) as run,):
                preparer._convert_to_rigid_pdbqt(
                    source, destination, as_ligand=True)

            command = run.call_args.args[0]
            output = destination.read_text()
        self.assertIn("--read_pdb", command)
        self.assertIn("--write_pdbqt", command)
        self.assertIn("ROOT\n", output)
        self.assertIn("TORSDOF 0\n", output)

    def test_peptide_keeps_modified_residues_but_excludes_water(self):
        content = (
            pdb_atom_line(1, "B", 0, 0, 0)
            + pdb_atom_line(2, "B", 1, 0, 0, record="HETATM", residue="MLY")
            + pdb_atom_line(
                3, "B", 2, 0, 0, element="O", record="HETATM", residue="HOH"
            ) + "END\n")
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "peptide.pdb"
            source.write_text(content)
            preparer = preparation.PdbPreparer(hydrogen_policy="keep")
            receptor_coordinates = preparer.selected_coordinates(source)
            peptide_coordinates = preparer.selected_coordinates(
                source, as_ligand=True)

        self.assertEqual(receptor_coordinates, [(0.0, 0.0, 0.0)])
        self.assertEqual(
            peptide_coordinates, [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)],)

    def test_mmcif_input_is_supported(self):
        from Bio.PDB import MMCIFIO, PDBParser

        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            pdb = directory / "input.pdb"
            mmcif = directory / "input.cif"
            pdb.write_text(pdb_atom_line(1, "B", 1, 2, 3) + "END\n")
            structure = PDBParser(QUIET=True).get_structure("input", pdb)
            writer = MMCIFIO()
            writer.set_structure(structure)
            writer.save(str(mmcif))

            coordinates = preparation.PdbPreparer(
                hydrogen_policy="keep").selected_coordinates(mmcif, ("B",))

        self.assertEqual(coordinates, [(1.0, 2.0, 3.0)])

    def test_unsupported_modified_residue_fails_explicitly(self):
        content = pdb_atom_line(
            1, "B", 0, 0, 0, record="HETATM", residue="ZZZ")
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "peptide.pdb"
            selected = Path(temporary) / "selected.pdb"
            source.write_text(content + "END\n")
            preparer = preparation.PdbPreparer(hydrogen_policy="keep")
            preparer._write_selected_structure(
                source, selected, None, as_ligand=True)
            with self.assertRaisesRegex(
                ValueError, "Unsupported Meeko residue template.*ZZZ"):
                preparer._selected_atom_metrics(selected)

    def test_meeko_warning_is_fatal(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = directory / "peptide.pdb"
            destination = directory / "peptide.pdbqt"
            source.write_text(pdb_atom_line(1, "B", 0, 0, 0))

            def fake_run(command, **kwargs):
                destination.write_text(
                    pdbqt_atom_line(1, "C", "ALA", "B", 1, 0, 0, 0, "C"))
                return subprocess.CompletedProcess(
                    command, 0, "", "Warning: unsafe residue chemistry")

            preparer = preparation.PdbPreparer(hydrogen_policy="keep")
            with (
                mock.patch.object(preparation, "find_meeko_receptor_preparer",
                    return_value=Path("mk_prepare_receptor.py"),
                ), mock.patch.object(
                    preparation.subprocess, "run", side_effect=fake_run
                ), self.assertRaisesRegex(RuntimeError, "unsafe chemistry"),):
                preparer._convert_to_rigid_pdbqt(
                    source, destination, as_ligand=True)


class BatchTests(unittest.TestCase):
    def test_manifest_paths_are_relative_and_jobs_are_isolated(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            receptor = directory / "receptor.pdb"
            peptide_a = directory / "a.pdb"
            peptide_b = directory / "b.pdb"
            for filename in (receptor, peptide_a, peptide_b):
                filename.write_text(pdb_atom_line(1, "A", 0, 0, 0))
            manifest = directory / "jobs.csv"
            manifest.write_text("id,receptor_pdb,peptide_pdb\n"
                "a,receptor.pdb,a.pdb\n" "b,receptor.pdb,b.pdb\n")
            jobs = batch.read_manifest(manifest)
            groups = group_jobs(jobs)
        self.assertEqual(len(groups), 2)
        self.assertTrue(all(len(group.jobs) == 1 for group in groups))

    def test_memory_budget_rejects_a_group_before_scoring(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            receptor = directory / "receptor.pdb"
            peptide = directory / "peptide.pdb"
            receptor.write_text(pdb_atom_line(1, "A", 0, 0, 0))
            peptide.write_text(pdb_atom_line(1, "B", 0, 0, 0))
            job = ScoreJob(index=0,
                job_id="limited", receptor_pdb=receptor, peptide_pdb=peptide,)
            results = batch.MemoryAwareBatchRunner(
                WorkflowConfig(), max_workers=2, memory_budget_bytes=1
            ).run([job])
        self.assertEqual(results[0]["status"], "error")
        self.assertIn("exceeds the batch budget", results[0]["error"])

    def test_cli_defaults_to_keeping_hydrogens(self):
        args = vina_score.argument_parser([])
        self.assertEqual(args.hydrogen_policy, "keep")
        self.assertEqual(vina_score.argument_parser(
                ["--hydrogen-policy", "rebuild"]).hydrogen_policy,
            "rebuild",)


if __name__ == "__main__":
    unittest.main()
