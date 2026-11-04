from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "research" / "audit_gluhut_graph_chemistry.py"
spec = importlib.util.spec_from_file_location("graph_chemistry", MODULE_PATH)
assert spec and spec.loader
module = importlib.util.module_from_spec(spec)
sys.modules["graph_chemistry"] = module
spec.loader.exec_module(module)


class GraphChemistryTests(unittest.TestCase):
    def test_hill_formula_order(self) -> None:
        formula, counts = module.formula_from_atomic_numbers(
            [8] * 6 + [1] * 12 + [6] * 6
        )
        self.assertEqual(formula, "C6H12O6")
        self.assertEqual(counts, {"C": 6, "H": 12, "O": 6})

    def test_noncarbon_formula_is_alphabetical(self) -> None:
        formula, counts = module.formula_from_atomic_numbers([1, 1, 8])
        self.assertEqual(formula, "H2O")
        self.assertEqual(counts, {"H": 2, "O": 1})

    def test_components_respect_atom_subset(self) -> None:
        adjacency = [set() for _ in range(7)]
        for left, right in ((1, 2), (2, 3), (3, 4), (5, 6)):
            adjacency[left].add(right)
            adjacency[right].add(left)
        self.assertEqual(module.components({1, 2, 3, 4, 5, 6}, adjacency), [[1, 2, 3, 4], [5, 6]])
        self.assertEqual(module.components({2, 3, 5}, adjacency), [[2, 3], [5]])

    def test_urea_and_amide_motif_classification(self) -> None:
        # Atom 1: urea carbonyl C linked to O2, N3 and N4.
        # Atom 5: amide carbonyl C linked to O6, N7 and C8.
        atomic_numbers = [0, 6, 8, 7, 7, 6, 8, 7, 6]
        names = ["", "CU", "OU", "NU1", "NU2", "CA", "OA", "NA", "CX"]
        types = ["", "c", "o", "n", "n", "c", "o", "n", "c3"]
        adjacency = [set() for _ in range(9)]
        bonds = {}
        for left, right, length in (
            (1, 2, 1.23), (1, 3, 1.36), (1, 4, 1.36),
            (5, 6, 1.23), (5, 7, 1.36), (5, 8, 1.52),
        ):
            adjacency[left].add(right)
            adjacency[right].add(left)
            bonds[tuple(sorted((left, right)))] = {
                "atoms_1_based": list(sorted((left, right))),
                "equilibrium_A": length,
                "force_constant": 1.0,
            }
        motifs = module.classify_carbonyl_motifs(
            set(range(1, 9)), adjacency, bonds, atomic_numbers, names, types
        )
        by_atom = {row["carbon_atom_1_based"]: row for row in motifs}
        self.assertEqual(by_atom[1]["classification"], "urea_carbonyl")
        self.assertEqual(by_atom[5]["classification"], "amide_carbonyl")

    def test_methyl_group_detection(self) -> None:
        # Carbon 1 has H2-H4 plus heavy atom 5.
        atomic_numbers = [0, 6, 1, 1, 1, 6]
        names = ["", "CM", "H1", "H2", "H3", "C2"]
        types = ["", "c3", "hc", "hc", "hc", "c3"]
        adjacency = [set() for _ in range(6)]
        for neighbour in (2, 3, 4, 5):
            adjacency[1].add(neighbour)
            adjacency[neighbour].add(1)
        rows = module.methyl_groups(
            set(range(1, 6)), adjacency, atomic_numbers, names, types
        )
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["carbon_atom_1_based"], 1)
        self.assertEqual(rows[0]["attached_heavy_atom_1_based"], 5)

    def test_safe_member_rejects_traversal(self) -> None:
        for value in ("../x", "/x", "a\\b"):
            with self.assertRaises(module.ChemistryAuditError):
                module.safe_member(value)


if __name__ == "__main__":
    unittest.main()
