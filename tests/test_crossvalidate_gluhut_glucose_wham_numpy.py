from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "research" / "crossvalidate_gluhut_glucose_wham_numpy.py"
spec = importlib.util.spec_from_file_location("numpy_wham", MODULE_PATH)
assert spec and spec.loader
module = importlib.util.module_from_spec(spec)
sys.modules["numpy_wham"] = module
spec.loader.exec_module(module)


class NumPyWhamCrosscheckTests(unittest.TestCase):
    def test_center_parser(self):
        self.assertEqual(module.center_from_token("0-0"), 0.0)
        self.assertEqual(module.center_from_token("12-6"), 12.6)
        self.assertEqual(module.center_from_token("2-25"), 2.25)
        with self.assertRaises(module.CrosscheckError):
            module.center_from_token("2.5")

    def test_coordinate_parser_uses_last_numeric_column(self):
        payload = b"# comment\n@ metadata\n0 1.25\n1 2.50 other\n"
        values = module.parse_coordinate_payload(payload)
        np.testing.assert_allclose(values, [1.25, 2.5])

    def test_unsafe_member_rejected(self):
        for name in ("../escape.dat", "/absolute.dat", "a\\b.dat"):
            with self.assertRaises(module.CrosscheckError):
                module.safe_member(name)

    def test_bhattacharyya_is_one_for_identical_histograms(self):
        counts = np.asarray([[1.0, 2.0, 3.0], [1.0, 2.0, 3.0]])
        coefficient = module.adjacent_bhattacharyya(counts)
        self.assertAlmostEqual(float(coefficient[0]), 1.0, places=12)

    def test_solver_recovers_known_harmonic_potential(self):
        rng = np.random.default_rng(2215)
        underlying_k = 0.8
        umbrella_k = 8.0
        beta = 1.0 / (module.K_B_KCAL_MOL_K * module.TEMPERATURE_K)
        windows = []
        for index, center in enumerate(np.linspace(-3.0, 3.0, 13)):
            mean = umbrella_k * center / (underlying_k + umbrella_k)
            sd = (beta * (underlying_k + umbrella_k)) ** -0.5
            values = rng.normal(mean, sd, 40000)
            windows.append((float(center), values, f"window-{index}.dat"))

        original_minimum = module.MINIMUM_A
        original_maximum = module.MAXIMUM_A
        original_bins = module.BINS
        original_force = module.FORCE_CONSTANT
        original_bulk_min = module.BULK_MINIMUM_A
        original_bulk_max = module.BULK_MAXIMUM_A
        try:
            module.MINIMUM_A = -4.0
            module.MAXIMUM_A = 4.0
            module.BINS = 320
            module.FORCE_CONSTANT = umbrella_k
            module.BULK_MINIMUM_A = 3.0
            module.BULK_MAXIMUM_A = 4.0
            result = module.solve_route(windows)
        finally:
            module.MINIMUM_A = original_minimum
            module.MAXIMUM_A = original_maximum
            module.BINS = original_bins
            module.FORCE_CONSTANT = original_force
            module.BULK_MINIMUM_A = original_bulk_min
            module.BULK_MAXIMUM_A = original_bulk_max

        profile = result["pmf_profile"]
        x = np.asarray([row["coordinate_A"] for row in profile], dtype=float)
        y = np.asarray(
            [
                np.nan
                if row["relative_free_energy_kcal_mol"] is None
                else row["relative_free_energy_kcal_mol"]
                for row in profile
            ],
            dtype=float,
        )
        finite = np.isfinite(y) & (np.abs(x) <= 2.2)
        coefficient, _ = np.polyfit(x[finite] ** 2, y[finite], 1)
        self.assertAlmostEqual(coefficient, 0.5 * underlying_k, delta=0.05)
        self.assertGreater(result["minimum_adjacent_bhattacharyya"], 0.01)


if __name__ == "__main__":
    unittest.main()
