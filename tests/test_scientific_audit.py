"""Synthetic regressions for scientific identity and control accounting."""
# ruff: noqa: E402 -- Exercise the source checkout in public and private layouts.
from copy import deepcopy
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "skinpipe" if (ROOT / "skinpipe/skinpipe").is_dir() else ROOT
sys.path.insert(0, str(SOURCE))

from rdkit import Chem
from rdkit.Chem import AllChem
from skinpipe.dock import poses
from skinpipe.platform import algorithm_chemistry, algorithm_evaluation
from skinpipe.platform.algorithm_profiles import select_bands
from skinpipe.platform.confirmation import confirmation_policy, select_confirmation
from skinpipe.platform.evaluation import evaluate_ranked_candidates
from skinpipe.platform.planning import plan_target
from skinpipe.platform.references import measurement


def molecule(smiles="C[C@H](O)F"):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    if AllChem.EmbedMolecule(mol, randomSeed=42) != 0:
        raise AssertionError("synthetic molecule did not embed")
    return mol


def candidate(identifier, smiles, **extra):
    return {"candidate_id": identifier, "target_id": "synthetic-target",
            "smiles": smiles, "priority_score_v2": 80.0, "purchasable": True,
            **extra}


class PoseIdentityTests(unittest.TestCase):
    def test_reflected_defined_stereocentre_is_rejected(self):
        mol = molecule()
        for index in range(mol.GetNumAtoms()):
            point = mol.GetConformer().GetAtomPosition(index)
            point.x = -point.x
            mol.GetConformer().SetAtomPosition(index, point)
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "reflected.sdf"
            self.assertIsNone(poses.write_pose_sdf(mol, 0, path))
            self.assertFalse(path.exists())

    def test_matching_defined_stereochemistry_roundtrips(self):
        mol = molecule("[13CH3][C@H](O)F")
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "correct.sdf"
            self.assertEqual(poses.write_pose_sdf(mol, 0, path), path)
            self.assertEqual(Chem.MolToSmiles(Chem.RemoveHs(poses.load_single_mol(path))),
                             Chem.MolToSmiles(Chem.RemoveHs(mol)))

    def test_unassigned_source_stereo_can_gain_coordinate_assignment(self):
        # This low-level writer supports historical inspection. Preparation
        # separately rejects unspecified persistent source stereochemistry.
        mol = molecule("CC(O)F")
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "unassigned.sdf"
            self.assertEqual(poses.write_pose_sdf(mol, 0, path), path)

    def test_multiple_sdf_records_cannot_be_one_pose(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "multiple.sdf"
            with Chem.SDWriter(str(path)) as writer:
                writer.write(molecule("CCO"))
                writer.write(molecule("CCN"))
            self.assertFalse(poses.sdf_roundtrips(path))

    def test_infinite_pose_coordinate_is_rejected(self):
        mol = molecule("CCO")
        point = mol.GetConformer().GetAtomPosition(0)
        point.x = float("inf")
        mol.GetConformer().SetAtomPosition(0, point)
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "nonfinite.sdf"
            self.assertIsNone(poses.write_pose_sdf(mol, 0, path))
            self.assertFalse(path.exists())

    def test_malformed_second_record_cannot_be_discarded(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "malformed.sdf"
            with Chem.SDWriter(str(path)) as writer:
                writer.write(molecule("CCO"))
            with path.open("a") as stream:
                stream.write("broken\nrecord\n\nthis is not a counts line\n$$$$\n")
            self.assertFalse(poses.sdf_roundtrips(path))

    def test_multiple_rebuilt_molecules_cannot_be_first_molecule_only(self):
        for extra in (molecule("CCN"), None):
            with self.subTest(second_record_valid=extra is not None), patch(
                    "meeko.PDBQTMolecule"), patch(
                    "meeko.RDKitMolCreate.from_pdbqt_mol",
                    return_value=[molecule("CCO"), extra]):
                self.assertIsNone(poses.mols_from_docked_pdbqt("synthetic"))


class ControlAccountingTests(unittest.TestCase):
    def setUp(self):
        self.rows = [candidate("reference", "CCN", role="reference_control"),
                     candidate("alias", "CC[NH3+].[Cl-]"),
                     candidate("new", "CCO")]

    def test_reference_control_does_not_count_as_new_hit(self):
        reference = {**self.rows[0], "score": 99.0, "label": 1, "label_measured": True}
        novel = {**self.rows[2], "score": 1.0, "label": 0, "label_measured": True}
        output = evaluate_ranked_candidates([reference, novel], k=1)
        self.assertEqual(output["control_count"], 1)
        self.assertEqual(output["candidate_count"], 1)
        self.assertEqual(output["hit_at_k"], 0)
        self.assertEqual(output["pr_auc"], None)

    def test_control_parent_alias_cannot_inflate_benchmark(self):
        for role in ("positive_control", "reference_control"):
            with self.subTest(role=role):
                rows = [{**self.rows[0], "role": role, "parent_id": "same-parent", "score": 99.0},
                        {**self.rows[1], "parent_id": "same-parent", "score": 90.0}]
                with self.assertRaisesRegex(ValueError, "control parent"):
                    evaluate_ranked_candidates(rows, k=1)

    def test_learning_identity_gate_recognizes_reference_control_alias(self):
        self.assertTrue(algorithm_evaluation._control(self.rows[0]))
        with self.assertRaisesRegex(ValueError, "control parent"):
            algorithm_evaluation._identities(self.rows)

    def test_confirmation_excludes_reference_and_its_parent_alias(self):
        results = [{"candidate_id": row["candidate_id"], "status": "complete"}
                   for row in self.rows]
        self.assertEqual(select_confirmation(self.rows, results, confirmation_policy({})), ["new"])

    def test_algorithm_bands_exclude_reference_and_its_parent_alias(self):
        selected = select_bands(self.rows, "synthetic-target", "dimorphite")
        self.assertEqual([row["candidate_id"] for row in selected["candidates"]], ["new"])

    def test_plan_keeps_reference_as_a_separate_confirmation_control(self):
        context = {"target_id": "synthetic-target", "target_identity_valid": True,
                   "actual_protein_id": "synthetic-protein", "receptor_id": "synthetic-receptor",
                   "pocket_id": "synthetic-pocket"}
        result = plan_target(self.rows, "synthetic-target", snapshot_id="synthetic-snapshot",
                             target_context=context)
        self.assertEqual(result["population_count"], 3)
        self.assertEqual(result["stages"]["broad_docking"]["candidate_ids"], ["new"])
        self.assertEqual(result["stages"]["confirmation"]["control_ids"], ["reference"])

    def test_chemistry_keeps_reference_and_alias_outside_candidate_compute(self):
        original = deepcopy(self.rows)
        with tempfile.TemporaryDirectory() as temporary, patch.object(
                algorithm_chemistry, "_dimorphite", return_value=({"synthetic": True}, [])):
            output = algorithm_chemistry.execute_chemistry(
                "dimorphite", self.rows, {}, workdir=Path(temporary))
        self.assertEqual(self.rows, original)
        reasons = {row["candidate_id"]: row.get("reason") for row in output["candidates"]}
        self.assertEqual(reasons["reference"], "controls_separate")
        self.assertEqual(reasons["alias"], "controls_separate")
        self.assertEqual(output["measurements"]["selected_count"], 1)


class MeasuredReferenceTests(unittest.TestCase):
    @staticmethod
    def record(**overrides):
        raw = {"target_chembl_id": "SYNTHETIC_TARGET", "canonical_smiles": "CCO",
               "standard_value": "1", "standard_units": "nM", "standard_type": "Ki",
               "standard_relation": "=", "assay_chembl_id": "SYNTHETIC_ASSAY",
               **overrides}
        descriptor = {"kind": "protein", "target_id": "synthetic-target",
                      "chembl_target_id": "SYNTHETIC_TARGET", "organism": "Homo sapiens"}
        assay = {"confidence_score": "9", "assay_organism": "Homo sapiens", "assay_type": "B"}
        return raw, descriptor, assay

    def test_boolean_concentration_is_not_a_one_nanomolar_reference(self):
        raw, descriptor, assay = self.record(standard_value=True)
        output = measurement(raw, descriptor, assay)
        self.assertIsNone(output["value"])
        self.assertFalse(output["is_measured"])
        self.assertFalse(output["reference_eligible"])
        self.assertFalse(output["ranking_eligible"])
        self.assertEqual(output["original"], raw)

    def test_explicit_unmeasured_comments_override_numeric_placeholders(self):
        for comment in ("Not Evaluated", "Not Determined", "not tested",
                        "Not Assayed", "Not Measured", "Not Available"):
            with self.subTest(comment=comment):
                raw, descriptor, assay = self.record(activity_comment=comment)
                output = measurement(raw, descriptor, assay)
                self.assertEqual(output["value"], 1.0)
                self.assertFalse(output["is_measured"])
                self.assertFalse(output["reference_eligible"])
                self.assertFalse(output["ranking_eligible"])
                self.assertIn("measurement_not_evaluated", output["reference_hold_reasons"])
                self.assertEqual(output["original"], raw)

    def test_source_decimal_string_remains_an_exact_reference(self):
        raw, descriptor, assay = self.record()
        output = measurement(raw, descriptor, assay)
        self.assertEqual(output["value"], 1.0)
        self.assertTrue(output["reference_eligible"])
        self.assertTrue(output["ranking_eligible"])

    def test_fractional_confidence_cannot_be_truncated_into_direct_mapping(self):
        raw, descriptor, assay = self.record()
        for value in (9.5, True):
            with self.subTest(confidence=value):
                output = measurement(raw, descriptor, {**assay, "confidence_score": value})
                self.assertFalse(output["reference_eligible"])
                self.assertIn("target_mapping_not_direct", output["reference_hold_reasons"])
        # The complete raw assay context must also be canonical JSON before
        # it can supply a comparison-cohort digest.
        with self.assertRaises(ValueError):
            measurement(raw, descriptor, {**assay, "confidence_score": float("inf")})


if __name__ == "__main__":
    unittest.main()
