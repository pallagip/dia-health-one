from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "research" / "audit_gluhut_binpos_coordinates.py"
spec = importlib.util.spec_from_file_location("binpos_audit", MODULE_PATH)
assert spec and spec.loader
module = importlib.util.module_from_spec(spec)
sys.modules["binpos_audit"] = module
spec.loader.exec_module(module)


class BinposCoordinateAuditTests(unittest.TestCase):
    def test_mass_weighted_com(self) -> None:
        xyz = np.asarray(
            [
                [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
                [[1.0, 0.0, 0.0], [3.0, 0.0, 0.0]],
            ]
        )
        com = module.mass_weighted_com(xyz, [1, 2], [1.0, 3.0])
        np.testing.assert_allclose(com[:, 0], [1.5, 2.5])
        np.testing.assert_allclose(com[:, 1:], 0.0)

    def test_reaction_coordinate(self) -> None:
        xyz = np.asarray(
            [
                [[0.0, 0.0, 0.0], [3.0, 4.0, 0.0]],
                [[0.0, 0.0, 0.0], [6.0, 8.0, 0.0]],
            ]
        )
        coordinate = module.reaction_coordinate_A(
            xyz, [1], [2], [12.0, 16.0]
        )
        np.testing.assert_allclose(coordinate, [5.0, 10.0])

    def test_trajectory_summary(self) -> None:
        coordinate = np.asarray([1.0, 2.0, 1.5, 4.0])
        summary = module.trajectory_summary(coordinate)
        self.assertEqual(summary["frame_count"], 4)
        self.assertEqual(summary["first_A"], 1.0)
        self.assertEqual(summary["last_A"], 4.0)
        self.assertEqual(summary["net_change_A"], 3.0)
        self.assertAlmostEqual(summary["fraction_positive_frame_steps"], 2 / 3)

    def test_delta_summary(self) -> None:
        deposited = np.asarray([1.0, 2.0, 3.0])
        alternative = np.asarray([1.5, 1.5, 4.0])
        summary = module.delta_summary(alternative, deposited)
        self.assertAlmostEqual(summary["first_delta_A"], 0.5)
        self.assertAlmostEqual(summary["last_delta_A"], 1.0)
        self.assertAlmostEqual(summary["mean_signed_delta_A"], 1.0 / 3.0)
        self.assertAlmostEqual(summary["mean_absolute_delta_A"], 2.0 / 3.0)
        self.assertAlmostEqual(
            summary["rms_delta_A"], np.sqrt((0.25 + 0.25 + 1.0) / 3.0)
        )

    def test_first_frame_bond_statistics(self) -> None:
        xyz = np.asarray(
            [
                [
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [1.0, 1.5, 0.0],
                ]
            ]
        )
        result = module.first_frame_bond_statistics(xyz, [(1, 2), (2, 3)])
        self.assertEqual(result["minimum_A"], 1.0)
        self.assertEqual(result["maximum_A"], 1.5)
        self.assertEqual(result["median_A"], 1.25)

    def test_alternative_map(self) -> None:
        record = {
            "viable_all_heavy_connected_hypotheses": [
                {
                    "label": "five_ring",
                    "atom_indices_1_based": [1, 2, 3, 4, 5],
                },
                {
                    "label": "ring_plus_C6",
                    "atom_indices_1_based": [1, 2, 3, 4, 5, 6],
                },
            ]
        }
        self.assertEqual(
            module.alternative_map(record),
            {
                "five_ring": [1, 2, 3, 4, 5],
                "ring_plus_C6": [1, 2, 3, 4, 5, 6],
            },
        )

    def test_invalid_mass_rejected(self) -> None:
        xyz = np.zeros((1, 2, 3))
        with self.assertRaises(module.CoordinateAuditError):
            module.mass_weighted_com(xyz, [1, 2], [12.0, 0.0])

    def test_safe_member_rejects_traversal(self) -> None:
        for value in ("../x", "/x", "a\\b"):
            with self.assertRaises(module.CoordinateAuditError):
                module.safe_member(value)


if __name__ == "__main__":
    unittest.main()
