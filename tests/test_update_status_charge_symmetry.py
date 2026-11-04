from __future__ import annotations

import importlib.util
import sys
import unittest
from copy import deepcopy
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "research" / "update_status_charge_symmetry.py"
spec = importlib.util.spec_from_file_location("status_updater", MODULE_PATH)
assert spec and spec.loader
module = importlib.util.module_from_spec(spec)
sys.modules["status_updater"] = module
spec.loader.exec_module(module)


def status_fixture() -> dict:
    return {
        "schema": "old",
        "updated": "old",
        "working_drug_discovered": True,
        "measured_dia_health_one_hit": True,
        "primary_candidate": {
            "id": "MC-Q1",
            "evidence_level": "E4",
            "candidate_ranking_authorized": True,
        },
        "public_benchmark_calibration": {
            "pag_calibration_completed": True,
            "candidate_ranking_authorized": True,
            "existing_field": "preserved",
        },
        "aws_execution": {
            "candidate_screen_authorized": True,
            "current_gate": "old",
            "next_action": "old",
        },
        "next_confidence_changes": ["existing"],
    }


def audit_fixture() -> dict:
    return {
        "cross_system": {
            "host_charge_vectors_identical": True,
            "host_charge_vector_sha256": "abc",
        },
        "graph_symmetry": {
            "automorphism_count": 6,
            "heavy_orbit_count": 28,
            "nontrivial_heavy_orbit_count": 28,
            "heavy_orbit_size_distribution": {"3": 28},
        },
        "summary": {
            "maximum_nontrivial_heavy_orbit_charge_range_e": 0.0,
            "maximum_nontrivial_hydrogen_group_charge_range_e": 0.0,
            "symmetry_averaging_adjustment": {
                "rms_per_atom_adjustment_e": 0.0,
            },
        },
        "decision": {
            "compiled_charges_are_exactly_graph_symmetry_averaged": True,
            "large_graph_equivalent_charge_asymmetry_detected": False,
            "automatic_charge_symmetrisation_authorized": False,
            "direct_transfer_to_pag_or_mcq1_authorized": False,
            "source_resp_derivation_reconstructable": False,
            "candidate_ranking_authorized": False,
        },
    }


class StatusUpdaterTests(unittest.TestCase):
    def test_updates_bounded_fields_and_preserves_existing_data(self) -> None:
        status = status_fixture()
        updated = module.update_status(status, audit_fixture())
        benchmark = updated["public_benchmark_calibration"]

        self.assertEqual(updated["schema"], module.STATUS_SCHEMA)
        self.assertEqual(benchmark["existing_field"], "preserved")
        self.assertEqual(
            benchmark["parameter_labeled_heavy_graph_automorphism_count"], 6
        )
        self.assertEqual(
            benchmark["parameter_labeled_heavy_graph_orbit_size_distribution"],
            {"3": 28},
        )
        self.assertTrue(
            benchmark[
                "compiled_charges_exactly_equal_within_preserved_graph_symmetries"
            ]
        )
        self.assertFalse(
            benchmark["large_graph_equivalent_charge_asymmetry_detected"]
        )
        self.assertFalse(
            benchmark["direct_compiled_charge_transfer_to_pag_or_mcq1_authorized"]
        )

    def test_restores_core_evidence_boundaries(self) -> None:
        updated = module.update_status(status_fixture(), audit_fixture())
        self.assertFalse(updated["working_drug_discovered"])
        self.assertFalse(updated["measured_dia_health_one_hit"])
        self.assertEqual(updated["primary_candidate"]["evidence_level"], "E2")
        self.assertFalse(
            updated["primary_candidate"]["candidate_ranking_authorized"]
        )
        self.assertFalse(
            updated["aws_execution"]["candidate_screen_authorized"]
        )
        self.assertFalse(
            updated["public_benchmark_calibration"]["pag_calibration_completed"]
        )
        self.assertFalse(
            updated["public_benchmark_calibration"]["candidate_ranking_authorized"]
        )

    def test_next_change_is_inserted_once(self) -> None:
        status = status_fixture()
        first = module.update_status(status, audit_fixture())
        second = module.update_status(first, audit_fixture())
        matches = [
            value
            for value in second["next_confidence_changes"]
            if value.startswith("obtain source RESP conformers")
        ]
        self.assertEqual(len(matches), 1)

    def test_rejects_authorized_parameter_transfer(self) -> None:
        audit = deepcopy(audit_fixture())
        audit["decision"]["direct_transfer_to_pag_or_mcq1_authorized"] = True
        with self.assertRaises(module.StatusUpdateError):
            module.update_status(status_fixture(), audit)

    def test_rejects_nonidentical_charge_vectors(self) -> None:
        audit = deepcopy(audit_fixture())
        audit["cross_system"]["host_charge_vectors_identical"] = False
        with self.assertRaises(module.StatusUpdateError):
            module.update_status(status_fixture(), audit)


if __name__ == "__main__":
    unittest.main()
