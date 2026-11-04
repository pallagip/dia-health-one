from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "research" / "audit_gluhut_force_convention.py"
spec = importlib.util.spec_from_file_location("force_convention", MODULE_PATH)
assert spec and spec.loader
module = importlib.util.module_from_spec(spec)
sys.modules["force_convention"] = module
spec.loader.exec_module(module)


class GluHUTForceConventionTests(unittest.TestCase):
    def test_energy_conversion(self) -> None:
        rows = module.energy_penalties(50.0)
        by_displacement = {row["displacement_A"]: row for row in rows}
        self.assertAlmostEqual(
            by_displacement[1.0]["amber_energy_kcal_mol"], 50.0
        )
        self.assertAlmostEqual(
            by_displacement[1.0]["standard_harmonic_energy_kcal_mol"],
            50.0,
        )
        self.assertAlmostEqual(
            by_displacement[0.5]["amber_energy_kcal_mol"], 12.5
        )

    def test_collect_smd_rk2(self) -> None:
        mapping = {
            "mappings": [
                {"system_name": "glucose", "parameters": {"rk2": 200.0}},
                {"system_name": "galactose", "parameters": {"rk2": 200.0}},
            ]
        }
        self.assertEqual(module.collect_smd_rk2(mapping), [200.0, 200.0])

    def test_missing_rk2_fails(self) -> None:
        with self.assertRaises(module.ForceConventionError):
            module.collect_smd_rk2(
                {"mappings": [{"system_name": "glucose", "parameters": {}}]}
            )

    def test_full_audit(self) -> None:
        mapping = {
            "mappings": [
                {
                    "system_name": name,
                    "parameters": {"rk2": 200.0},
                }
                for name in (
                    "glucose",
                    "galactose",
                    "fructopyranose",
                    "fructofuranose",
                )
            ]
        }
        wham = {
            "implementation": {
                "force_constant_kcal_mol_A2": 100.0,
            }
        }
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            mapping_path = root / "mapping.json"
            wham_path = root / "wham.json"
            output_path = root / "output.json"
            mapping_path.write_text(json.dumps(mapping), encoding="utf-8")
            wham_path.write_text(json.dumps(wham), encoding="utf-8")
            report = module.audit(mapping_path, wham_path, output_path)

            self.assertTrue(
                report["decision"]["force_constant_scalar_convention_validated"]
            )
            self.assertEqual(
                report["smd"]["standard_harmonic_K_kcal_mol_A2"], 400.0
            )
            self.assertEqual(
                report["umbrella"][
                    "expected_grossfield_wham_K_kcal_mol_A2"
                ],
                100.0,
            )
            self.assertFalse(
                report["decision"]["new_smd_or_umbrella_generation_authorized"]
            )
            self.assertFalse(
                report["decision"]["parent_pag_calibration_authorized"]
            )
            self.assertEqual(
                json.loads(output_path.read_text(encoding="utf-8")), report
            )

    def test_multiple_smd_coefficients_fail(self) -> None:
        mapping = {
            "mappings": [
                {"system_name": "glucose", "parameters": {"rk2": 200.0}},
                {"system_name": "galactose", "parameters": {"rk2": 50.0}},
            ]
        }
        wham = {"implementation": {"force_constant_kcal_mol_A2": 100.0}}
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            mapping_path = root / "mapping.json"
            wham_path = root / "wham.json"
            output_path = root / "output.json"
            mapping_path.write_text(json.dumps(mapping), encoding="utf-8")
            wham_path.write_text(json.dumps(wham), encoding="utf-8")
            with self.assertRaises(module.ForceConventionError):
                module.audit(mapping_path, wham_path, output_path)


if __name__ == "__main__":
    unittest.main()
