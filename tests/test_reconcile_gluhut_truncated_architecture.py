from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "research" / "reconcile_gluhut_truncated_architecture.py"
spec = importlib.util.spec_from_file_location("architecture", MODULE_PATH)
assert spec and spec.loader
module = importlib.util.module_from_spec(spec)
sys.modules["architecture"] = module
spec.loader.exec_module(module)


def atom(index: int, name: str, element: str, atom_type: str) -> dict:
    return {
        "atom_1_based": index,
        "atom_name": name,
        "element": element,
        "atomic_number": {"H": 1, "C": 6, "N": 7, "O": 8}[element],
        "amber_atom_type": atom_type,
        "charge_e": 0.0,
        "residue_index_1_based": 1,
        "residue_label": "MOL",
    }


def connect(topology: dict, left: int, right: int) -> None:
    topology["bonds"].append(tuple(sorted((left, right))))


def synthetic_ring(classification: str) -> tuple[dict, set[int], list[int]]:
    topology = {"atoms": {}, "bonds": []}
    ring = list(range(1, 7))
    for index in ring:
        topology["atoms"][index] = atom(index, f"R{index}", "C", "ca")
    for left, right in ((1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 1)):
        connect(topology, left, right)
    next_index = 7
    if classification == "cap":
        for ring_atom in ring:
            topology["atoms"][next_index] = atom(next_index, f"B{next_index}", "C", "c3")
            connect(topology, ring_atom, next_index)
            next_index += 1
    elif classification == "pillar":
        for ring_atom, element, atom_type in ((1, "N", "n"), (3, "N", "n"), (5, "C", "c")):
            topology["atoms"][next_index] = atom(next_index, f"B{next_index}", element, atom_type)
            connect(topology, ring_atom, next_index)
            next_index += 1
        for ring_atom in (2, 4, 6):
            topology["atoms"][next_index] = atom(next_index, f"H{next_index}", "H", "ha")
            connect(topology, ring_atom, next_index)
            next_index += 1
    else:
        raise ValueError(classification)
    return topology, set(topology["atoms"]), ring


class ArchitectureTests(unittest.TestCase):
    def test_module_formula_matches_compiled_formula(self) -> None:
        total = module.module_formula_total(module.EXPECTED_MODULES)
        self.assertEqual(total, module.EXPECTED_HOST_FORMULA)
        self.assertEqual(module.formula_string(total), "C60H75N15O9")

    def test_formula_arithmetic(self) -> None:
        total: dict[str, int] = {}
        module.add_formula(total, {"C": 15, "H": 24, "N": 3}, 2)
        module.add_formula(total, {"C": 8, "H": 9, "N": 3, "O": 1}, 3)
        module.add_formula(total, {"C": 1, "O": 1}, 6)
        self.assertEqual(dict(sorted(total.items())), module.EXPECTED_HOST_FORMULA)

    def test_classifies_cap_ring(self) -> None:
        topology, host, ring = synthetic_ring("cap")
        graph = module.adjacency(topology)
        result = module.classify_ring("cap", ring, host, topology, graph)
        self.assertEqual(result["classification"], "triethyl_tris_methylene_cap")
        self.assertEqual(result["external_heavy_elements"], {"C": 6})
        self.assertEqual(result["ring_hydrogen_count"], 0)

    def test_classifies_pillar_ring(self) -> None:
        topology, host, ring = synthetic_ring("pillar")
        graph = module.adjacency(topology)
        result = module.classify_ring("pillar", ring, host, topology, graph)
        self.assertEqual(result["classification"], "diaminobenzamide_pillar")
        self.assertEqual(result["external_heavy_elements"], {"C": 1, "N": 2})
        self.assertEqual(result["ring_hydrogen_count"], 3)

    def test_aromatic_ring_component_requires_six_cycle(self) -> None:
        topology, host, ring = synthetic_ring("cap")
        graph = module.adjacency(topology)
        rings = module.aromatic_ring_components(host, topology, graph)
        self.assertEqual(rings, [ring])
        topology["bonds"].remove((1, 2))
        graph = module.adjacency(topology)
        with self.assertRaises(module.ReconciliationError):
            module.aromatic_ring_components(host, topology, graph)

    def test_formula_string_hill_order(self) -> None:
        self.assertEqual(
            module.formula_string({"O": 9, "N": 15, "H": 75, "C": 60}),
            "C60H75N15O9",
        )

    def test_safe_member_rejects_traversal(self) -> None:
        for value in ("../x", "/x", "a\\b"):
            with self.assertRaises(module.ReconciliationError):
                module.safe_member(value)


if __name__ == "__main__":
    unittest.main()
