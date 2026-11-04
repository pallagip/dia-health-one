from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "research" / "audit_gluhut_smd_restraints.py"
spec = importlib.util.spec_from_file_location("restraint_audit", MODULE_PATH)
assert spec and spec.loader
module = importlib.util.module_from_spec(spec)
sys.modules["restraint_audit"] = module
spec.loader.exec_module(module)


class RestraintAuditTests(unittest.TestCase):
    def test_assignment_parser(self) -> None:
        text = """
&rst
 iat=-1,-1,
 r2=0.0,
 rk2=200.0,
 igr1=1,2,3,4,
 igr2=5,6,7,
/
"""
        assignments = module.parse_namelist_assignments(text)
        self.assertEqual(module.as_int_list(assignments, "iat"), [-1, -1])
        self.assertEqual(module.as_int_list(assignments, "igr1"), [1, 2, 3, 4])
        self.assertEqual(module.as_int_list(assignments, "igr2"), [5, 6, 7])
        self.assertEqual(module.as_float(assignments, "rk2"), 200.0)
        self.assertIsNone(module.as_float(assignments, "r3"))

    def test_comments_and_fortran_floats(self) -> None:
        text = """
&rst
 rk2=2.0D+02, ! comment
 r2a=1.5d+01,
 igr1=1,2,
 igr2=3,4,
 iat=-1,-1,
/
"""
        assignments = module.parse_namelist_assignments(text)
        self.assertEqual(module.as_float(assignments, "rk2"), 200.0)
        self.assertEqual(module.as_float(assignments, "r2a"), 15.0)

    def test_induced_group(self) -> None:
        topology = {
            "natom": 5,
            "atoms": {
                1: {"atom_1_based": 1, "atom_name": "C1", "amber_atom_type": "ca", "atomic_number": 6, "element": "C", "charge_e": 0.0, "residue_index_1_based": 1, "residue_label": "MOL"},
                2: {"atom_1_based": 2, "atom_name": "C2", "amber_atom_type": "ca", "atomic_number": 6, "element": "C", "charge_e": 0.0, "residue_index_1_based": 1, "residue_label": "MOL"},
                3: {"atom_1_based": 3, "atom_name": "O1", "amber_atom_type": "Os", "atomic_number": 8, "element": "O", "charge_e": 0.0, "residue_index_1_based": 2, "residue_label": "0GB"},
                4: {"atom_1_based": 4, "atom_name": "C3", "amber_atom_type": "Cg", "atomic_number": 6, "element": "C", "charge_e": 0.0, "residue_index_1_based": 2, "residue_label": "0GB"},
                5: {"atom_1_based": 5, "atom_name": "H1", "amber_atom_type": "H1", "atomic_number": 1, "element": "H", "charge_e": 0.0, "residue_index_1_based": 2, "residue_label": "0GB"},
            },
            "bonds": [(1, 2), (3, 4), (4, 5)],
        }
        group = module.induced_group([3, 4], topology)
        self.assertEqual(group["element_counts"], {"C": 1, "O": 1})
        self.assertEqual(group["residue_counts"], {"0GB": 2})
        self.assertEqual(group["induced_edge_count"], 1)
        self.assertEqual(group["connected_component_count"], 1)
        self.assertEqual(group["induced_cycle_rank"], 0)

    def test_duplicate_group_indices_rejected(self) -> None:
        topology = {
            "natom": 2,
            "atoms": {
                1: {"atom_1_based": 1, "atom_name": "C1", "amber_atom_type": "ca", "atomic_number": 6, "element": "C", "charge_e": 0.0, "residue_index_1_based": 1, "residue_label": "MOL"},
                2: {"atom_1_based": 2, "atom_name": "C2", "amber_atom_type": "ca", "atomic_number": 6, "element": "C", "charge_e": 0.0, "residue_index_1_based": 1, "residue_label": "MOL"},
            },
            "bonds": [(1, 2)],
        }
        with self.assertRaises(module.RestraintAuditError):
            module.induced_group([1, 1], topology)

    def test_restraint_filename_mapping(self) -> None:
        payload = b"&rst\n iat=-1,-1,\n igr1=1,2,\n igr2=3,4,\n/\n"
        parsed = module.parse_restraint(
            "SMD_trajectories/Input-files/Glucose-and-galactose-SMD.RST",
            payload,
        )
        self.assertEqual(parsed["applies_to_systems"], ["glucose", "galactose"])
        self.assertEqual(parsed["iat"], [-1, -1])

    def test_safe_member_rejects_traversal(self) -> None:
        for value in ("../x", "/x", "a\\b"):
            with self.assertRaises(module.RestraintAuditError):
                module.safe_member(value)


if __name__ == "__main__":
    unittest.main()
