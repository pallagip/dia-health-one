from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "research" / "reproduce_gluhut_all_sugars.py"
spec = importlib.util.spec_from_file_location("all_sugars", MODULE_PATH)
assert spec and spec.loader
module = importlib.util.module_from_spec(spec)
sys.modules["all_sugars"] = module
spec.loader.exec_module(module)


class AllSugarToolTests(unittest.TestCase):
    def test_center_parser(self):
        self.assertEqual(module.center_from_token("0-0"), 0.0)
        self.assertEqual(module.center_from_token("12-6"), 12.6)
        self.assertEqual(module.center_from_token("2-25"), 2.25)
        with self.assertRaises(module.AnalysisError):
            module.center_from_token("2.5")

    def test_route_group_mapping(self):
        self.assertEqual(module.route_group("glucose", "US-1"), "glucose")
        self.assertEqual(module.route_group("galactose", "US-2"), "galactose-pathB")
        self.assertEqual(module.route_group("galactose", "US-routeA"), "galactose-pathA")

    def test_binding_energy_sign_and_bulk(self):
        profile = [(0.0, 0.0), (1.0, 1.0), (13.0, 8.0), (14.0, 10.0), (15.0, 9.0)]
        self.assertAlmostEqual(module.binding_energy(profile), -9.0)

    def test_expected_ranking_is_published_order(self):
        self.assertEqual(
            module.EXPECTED_RANKING,
            [
                "glucose",
                "galactose-pathB",
                "fructofuranose",
                "galactose-pathA",
                "fructopyranose",
            ],
        )

    def test_frame_count_ignores_comments_and_empty_lines(self):
        payload = b"# comment\n@ metadata\n0 1.0\n\n1 2.0\n"
        self.assertEqual(module.frame_count(payload), 2)


if __name__ == "__main__":
    unittest.main()
