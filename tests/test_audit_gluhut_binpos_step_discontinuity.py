from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "research" / "audit_gluhut_binpos_step_discontinuity.py"
spec = importlib.util.spec_from_file_location("step_audit", MODULE_PATH)
assert spec and spec.loader
module = importlib.util.module_from_spec(spec)
sys.modules["step_audit"] = module
spec.loader.exec_module(module)


class StepDiscontinuityTests(unittest.TestCase):
    def test_smooth_coordinate_has_no_pbc_like_jump(self) -> None:
        coordinate = np.linspace(0.0, 15.0, 1000)
        result = module.step_summary(coordinate)
        self.assertFalse(result["pbc_like_jump_detected"])
        self.assertLess(result["maximum_absolute_step_A"], 0.1)
        self.assertEqual(result["counts_above_threshold_A"]["5.0"], 0)

    def test_single_large_jump_detected(self) -> None:
        coordinate = np.concatenate(
            [np.linspace(0.0, 5.0, 500), np.linspace(15.0, 20.0, 500)]
        )
        result = module.step_summary(coordinate)
        self.assertTrue(result["pbc_like_jump_detected"])
        self.assertGreater(result["maximum_absolute_step_A"], 9.0)
        self.assertEqual(result["maximum_step_from_frame_index_0_based"], 499)
        self.assertEqual(result["counts_above_threshold_A"]["5.0"], 1)

    def test_negative_jump_detected(self) -> None:
        coordinate = np.asarray([0.0, 0.1, 0.2, -10.0, -9.9])
        result = module.step_summary(coordinate)
        self.assertTrue(result["pbc_like_jump_detected"])
        self.assertLess(result["maximum_step_signed_A"], -10.0)

    def test_percentile_summary(self) -> None:
        values = np.asarray([1.0, 2.0, 3.0, 4.0])
        result = module.percentile_summary(values)
        self.assertEqual(result["minimum"], 1.0)
        self.assertEqual(result["maximum"], 4.0)
        self.assertEqual(result["median"], 2.5)
        self.assertEqual(result["mean"], 2.5)

    def test_safe_member_rejects_traversal(self) -> None:
        for value in ("../x", "/x", "a\\b"):
            with self.assertRaises(module.StepAuditError):
                module.safe_member(value)


if __name__ == "__main__":
    unittest.main()
