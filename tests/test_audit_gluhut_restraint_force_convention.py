from __future__ import annotations

import importlib.util
import json
import math
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "research" / "audit_gluhut_restraint_force_convention.py"
spec = importlib.util.spec_from_file_location("force_convention", MODULE_PATH)
assert spec and spec.loader
module = importlib.util.module_from_spec(spec)
sys.modules["force_convention"] = module
spec.loader.exec_module(module)


class ForceConventionTests(unittest.TestCase):
    def test_round_trip_conversion(self) -> None:
        for rk2 in (0.0, 50.0, 200.0, 123.456):
            conventional = module.conventional_k_from_amber_rk2(rk2)
            self.assertAlmostEqual(
                module.amber_rk2_from_conventional_k(conventional),
                rk2,
                places=12,
            )

    def test_energy_forms_are_equivalent(self) -> None:
        for rk2 in (50.0, 200.0):
            conventional = module.conventional_k_from_amber_rk2(rk2)
            for displacement in (0.0, 0.1, 0.5, 1.0, 2.5):
                self.assertAlmostEqual(
                    module.amber_energy(rk2, displacement),
                    module.conventional_energy(conventional, displacement),
                    places=12,
                )

    def test_unique_numeric_values(self) -> None:
        mappings = [
            {"parameters": {"rk2": 200.0, "outxyz": None}},
            {"parameters": {"rk2": 200.0, "outxyz": 1.0}},
            {"parameters": {"rk2": 200.0, "outxyz": None}},
        ]
        self.assertEqual(module.unique_numeric_values(mappings, "rk2"), [200.0])
        self.assertEqual(
            module.unique_numeric_values(mappings, "outxyz"),
            [1.0, None],
        )

    def test_full_report(self) -> None:
        mappings = [
            {
                "system_name": system,
                "parameters": {
                    "rk2": 200.0,
                    "r2": 0.0,
                    "r2a": 15.0,
                },
            }
            for system in (
                "glucose",
                "galactose",
                "fructopyranose",
                "fructofuranose",
            )
        ]
        mapping_record = {"mappings": mappings}
        wham_record = {
            "analysis": {
                "parameters": {
                    "wham_metadata_spring_kcal_mol_angstrom2": 100.0,
                    "article_reported_amber_spring_kcal_mol_angstrom2": 50.0,
                }
            }
        }

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            mapping_path = root / "mapping.json"
            wham_path = root / "wham.json"
            output_path = root / "output.json"
            mapping_path.write_text(json.dumps(mapping_record), encoding="utf-8")
            wham_path.write_text(json.dumps(wham_record), encoding="utf-8")

            report = module.run(mapping_path, wham_path, output_path)

            self.assertEqual(report["smd"]["amber_rk2_kcal_mol_A2"], 200.0)
            self.assertEqual(
                report["smd"]["conventional_half_K_kcal_mol_A2"],
                400.0,
            )
            self.assertEqual(
                report["umbrella"]["deposited_wham_metadata_K_kcal_mol_A2"],
                100.0,
            )
            self.assertEqual(
                report["umbrella"]["amber_rk2_recovered_from_wham_K_kcal_mol_A2"],
                50.0,
            )
            self.assertTrue(
                report["decision"]["force_constant_convention_resolved"]
            )
            self.assertFalse(
                report["decision"]["guest_atom_selection_validated"]
            )
            self.assertFalse(
                report["decision"]["new_smd_or_umbrella_generation_authorized"]
            )
            self.assertEqual(json.loads(output_path.read_text()), report)

    def test_unexpected_smd_value_fails(self) -> None:
        mapping_record = {
            "mappings": [
                {
                    "parameters": {
                        "rk2": 100.0,
                        "r2": 0.0,
                        "r2a": 15.0,
                    }
                }
                for _ in range(4)
            ]
        }
        wham_record = {
            "analysis": {
                "parameters": {
                    "wham_metadata_spring_kcal_mol_angstrom2": 100.0,
                    "article_reported_amber_spring_kcal_mol_angstrom2": 50.0,
                }
            }
        }
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            mapping_path = root / "mapping.json"
            wham_path = root / "wham.json"
            output_path = root / "output.json"
            mapping_path.write_text(json.dumps(mapping_record), encoding="utf-8")
            wham_path.write_text(json.dumps(wham_record), encoding="utf-8")
            with self.assertRaises(module.ConventionError):
                module.run(mapping_path, wham_path, output_path)


if __name__ == "__main__":
    unittest.main()
