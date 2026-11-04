from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "research" / "audit_gluhut_short_window_sensitivity.py"
spec = importlib.util.spec_from_file_location("short_window", MODULE_PATH)
assert spec and spec.loader
module = importlib.util.module_from_spec(spec)
sys.modules["short_window"] = module
spec.loader.exec_module(module)


class ShortWindowSensitivityTests(unittest.TestCase):
    def test_center_parser(self):
        self.assertEqual(module.center_from_token("12-6"), 12.6)
        self.assertEqual(module.center_from_token("2-25"), 2.25)
        with self.assertRaises(module.SensitivityError):
            module.center_from_token("12.6")

    def test_payload_split_and_numeric_column(self):
        headers, rows = module.split_payload(
            b"# header\n@ metadata\n0 12.5\n1 12.6 extra\n"
        )
        self.assertEqual(len(headers), 2)
        self.assertEqual(len(rows), 2)
        np.testing.assert_allclose(module.last_numeric_column(rows), [12.5, 12.6])

    def test_identical_histograms_have_full_overlap(self):
        values = np.linspace(11.0, 13.0, 1000)
        self.assertAlmostEqual(module.bhattacharyya(values, values), 1.0, places=12)

    def test_separated_histograms_have_zero_overlap(self):
        left = np.linspace(1.0, 2.0, 1000)
        right = np.linspace(10.0, 11.0, 1000)
        self.assertAlmostEqual(module.bhattacharyya(left, right), 0.0, places=12)

    def test_binding_energy_sign_and_bulk_reference(self):
        profile = [(0.0, 0.0), (1.0, 2.0), (13.0, 8.0), (14.0, 10.0), (15.0, 9.0)]
        self.assertAlmostEqual(module.binding_energy(profile), -9.0)


if __name__ == "__main__":
    unittest.main()
