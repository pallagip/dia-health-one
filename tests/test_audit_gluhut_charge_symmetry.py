from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path

import networkx as nx

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "research" / "audit_gluhut_charge_symmetry.py"
spec = importlib.util.spec_from_file_location("charge_symmetry", MODULE_PATH)
assert spec and spec.loader
module = importlib.util.module_from_spec(spec)
sys.modules["charge_symmetry"] = module
spec.loader.exec_module(module)


class ChargeSymmetryTests(unittest.TestCase):
    def test_union_find(self) -> None:
        union = module.UnionFind([1, 2, 3, 4])
        union.union(1, 2)
        union.union(3, 4)
        union.union(2, 3)
        self.assertEqual(union.groups(), [[1, 2, 3, 4]])

    def test_statistics(self) -> None:
        result = module.statistics([-0.2, 0.0, 0.2])
        self.assertAlmostEqual(result["mean"], 0.0)
        self.assertAlmostEqual(result["range"], 0.4)
        self.assertAlmostEqual(
            result["population_standard_deviation"],
            (0.08 / 3.0) ** 0.5,
        )

    def test_automorphism_orbits_on_labeled_symmetric_graph(self) -> None:
        graph = nx.Graph()
        for index in (1, 2):
            graph.add_node(
                index,
                element="C",
                amber_atom_type="c3",
                full_degree=1,
                heavy_degree=1,
                hydrogen_count=3,
            )
        graph.add_edge(
            1,
            2,
            force_constant=310.0,
            equilibrium_A=1.526,
        )
        count, orbits = module.automorphism_orbits(graph)
        self.assertEqual(count, 2)
        self.assertEqual(orbits, [[1, 2]])

    def test_edge_labels_break_false_symmetry(self) -> None:
        graph = nx.Graph()
        for index in (1, 2, 3):
            graph.add_node(
                index,
                element="C",
                amber_atom_type="c3",
                full_degree=2 if index == 2 else 1,
                heavy_degree=2 if index == 2 else 1,
                hydrogen_count=0,
            )
        graph.add_edge(1, 2, force_constant=100.0, equilibrium_A=1.5)
        graph.add_edge(2, 3, force_constant=200.0, equilibrium_A=1.5)
        count, orbits = module.automorphism_orbits(graph)
        self.assertEqual(count, 1)
        self.assertEqual(sorted(orbits), [[1], [2], [3]])

    def test_orbit_and_hydrogen_charge_records(self) -> None:
        atoms = {
            1: {"atom_name": "C1", "element": "C", "amber_atom_type": "c3", "charge_e": 0.1},
            2: {"atom_name": "C2", "element": "C", "amber_atom_type": "c3", "charge_e": 0.2},
            3: {"atom_name": "H11", "element": "H", "amber_atom_type": "hc", "charge_e": 0.01},
            4: {"atom_name": "H12", "element": "H", "amber_atom_type": "hc", "charge_e": 0.02},
            5: {"atom_name": "H21", "element": "H", "amber_atom_type": "hc", "charge_e": 0.03},
            6: {"atom_name": "H22", "element": "H", "amber_atom_type": "hc", "charge_e": 0.04},
        }
        adjacency = {
            1: {2, 3, 4},
            2: {1, 5, 6},
            3: {1},
            4: {1},
            5: {2},
            6: {2},
        }
        host = {
            "atoms": atoms,
            "adjacency": adjacency,
            "host_atom_indices": [1, 2, 3, 4, 5, 6],
            "host_charge_e": 0.4,
        }
        heavy = module.orbit_records([[1, 2]], host)
        hydrogens = module.hydrogen_group_records([[1, 2]], host)
        self.assertEqual(heavy[0]["charge_statistics_e"]["range"], 0.1)
        self.assertEqual(len(hydrogens), 1)
        self.assertEqual(hydrogens[0]["hydrogen_count"], 4)
        self.assertAlmostEqual(
            hydrogens[0]["charge_statistics_e"]["range"], 0.03
        )
        adjustment = module.symmetry_averaging_adjustment(
            heavy, hydrogens, host
        )
        self.assertAlmostEqual(adjustment["total_charge_difference_e"], 0.0)
        self.assertGreater(adjustment["rms_per_atom_adjustment_e"], 0.0)

    def test_empty_statistics_rejected(self) -> None:
        with self.assertRaises(module.ChargeSymmetryError):
            module.statistics([])

    def test_safe_member_rejects_traversal(self) -> None:
        for value in ("../x", "/x", "a\\b"):
            with self.assertRaises(module.ChargeSymmetryError):
                module.safe_member(value)


if __name__ == "__main__":
    unittest.main()
