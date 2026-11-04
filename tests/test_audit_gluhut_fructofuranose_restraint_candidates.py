from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "research" / "audit_gluhut_fructofuranose_restraint_candidates.py"
spec = importlib.util.spec_from_file_location("candidate_audit", MODULE_PATH)
assert spec and spec.loader
module = importlib.util.module_from_spec(spec)
sys.modules["candidate_audit"] = module
spec.loader.exec_module(module)


def atom(index: int, name: str, element: str) -> dict:
    return {
        "atom_1_based": index,
        "atom_name": name,
        "amber_atom_type": "H1" if element == "H" else "Cg",
        "atomic_number": {"H": 1, "C": 6, "O": 8}[element],
        "element": element,
        "charge_e": 0.0,
        "residue_index_1_based": 2,
        "residue_label": "0CU",
    }


def topology() -> dict:
    # Five-member ring 1-2-3-4-5, exocyclic C1 atom 6 on ring atom 1,
    # exocyclic C6 atom 7 on ring atom 3, and H62 atom 8 attached to C6.
    atoms = {
        1: atom(1, "C2", "C"),
        2: atom(2, "O5", "O"),
        3: atom(3, "C5", "C"),
        4: atom(4, "C4", "C"),
        5: atom(5, "C3", "C"),
        6: atom(6, "C1", "C"),
        7: atom(7, "C6", "C"),
        8: atom(8, "H62", "H"),
    }
    bonds = [(1, 2), (2, 3), (3, 4), (4, 5), (5, 1), (1, 6), (3, 7), (7, 8)]
    return {"natom": 8, "atoms": atoms, "bonds": bonds}


class FructofuranoseCandidateTests(unittest.TestCase):
    def test_adjacency_finds_hydrogen_parent(self) -> None:
        graph = module.adjacency(topology())
        self.assertEqual(graph[8], {7})

    def test_replacement_summary_for_five_ring(self) -> None:
        result = module.replacement_summary(
            "omit_H62_use_five_ring_atoms", [1, 2, 3, 4, 5], topology()
        )
        self.assertTrue(result["all_selected_atoms_are_heavy"])
        self.assertTrue(result["single_connected_component"])
        self.assertTrue(result["retains_five_member_cyclic_core"])
        self.assertEqual(result["selection_class"], "cyclic_heavy_atom_core_only")

    def test_replacement_with_exocyclic_carbon(self) -> None:
        result = module.replacement_summary(
            "replace_H62_with_C6", [1, 2, 3, 4, 5, 7], topology()
        )
        self.assertTrue(result["all_selected_atoms_are_heavy"])
        self.assertTrue(result["single_connected_component"])
        self.assertTrue(result["retains_five_member_cyclic_core"])
        self.assertEqual(
            result["selection_class"],
            "cyclic_core_plus_exocyclic_heavy_atom",
        )
        self.assertEqual(result["noncore_selected_atoms"][0]["atom_name"], "C6")

    def test_selected_hydrogen_group_is_disconnected(self) -> None:
        result = module.replacement_summary(
            "deposited", [1, 2, 3, 4, 5, 8], topology()
        )
        self.assertFalse(result["all_selected_atoms_are_heavy"])
        self.assertFalse(result["single_connected_component"])
        self.assertTrue(result["retains_five_member_cyclic_core"])

    def test_safe_member_rejects_traversal(self) -> None:
        for value in ("../x", "/x", "a\\b"):
            with self.assertRaises(module.CandidateAuditError):
                module.safe_member(value)


if __name__ == "__main__":
    unittest.main()
