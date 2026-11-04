from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "research" / "summarize_gluhut_smd_restraint_anomaly.py"
spec = importlib.util.spec_from_file_location("restraint_summary", MODULE_PATH)
assert spec and spec.loader
module = importlib.util.module_from_spec(spec)
sys.modules["restraint_summary"] = module
spec.loader.exec_module(module)


def atom(index: int, name: str, element: str, residue: str = "SUG") -> dict:
    return {
        "atom_1_based": index,
        "atom_name": name,
        "element": element,
        "amber_atom_type": "ca" if element == "C" and residue == "MOL" else "X",
        "residue_label": residue,
    }


def group(indices: list[int], atoms: list[dict], edges: list[list[int]]) -> dict:
    components = []
    adjacency = {index: set() for index in indices}
    for left, right in edges:
        adjacency[left].add(right)
        adjacency[right].add(left)
    remaining = set(indices)
    while remaining:
        start = min(remaining)
        remaining.remove(start)
        stack = [start]
        component = []
        while stack:
            current = stack.pop()
            component.append(current)
            for neighbour in adjacency[current]:
                if neighbour in remaining:
                    remaining.remove(neighbour)
                    stack.append(neighbour)
        components.append(sorted(component))
    components.sort(key=lambda row: (len(row), row), reverse=True)
    element_counts = {}
    type_counts = {}
    residue_counts = {}
    for item in atoms:
        element_counts[item["element"]] = element_counts.get(item["element"], 0) + 1
        type_counts[item["amber_atom_type"]] = type_counts.get(item["amber_atom_type"], 0) + 1
        residue_counts[item["residue_label"]] = residue_counts.get(item["residue_label"], 0) + 1
    return {
        "atom_indices_1_based": indices,
        "atom_count": len(indices),
        "atoms": atoms,
        "element_counts": element_counts,
        "amber_atom_type_counts": type_counts,
        "residue_counts": residue_counts,
        "induced_edges": edges,
        "connected_component_count": len(components),
        "component_atom_indices": components,
        "induced_cycle_rank": len(edges) - len(indices) + len(components),
    }


def host_group() -> dict:
    indices = list(range(1, 13))
    edges = [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 1],
             [7, 8], [8, 9], [9, 10], [10, 11], [11, 12], [12, 7]]
    atoms = [atom(index, f"C{index}", "C", "MOL") for index in indices]
    return group(indices, atoms, edges)


def guest_ring_plus_c6(system: str) -> dict:
    indices = [20, 21, 22, 23, 24, 25, 26]
    names = ["C1", "O5", "C5", "C4", "C3", "C2", "C6"]
    atoms = [atom(index, name, "O" if name.startswith("O") else "C", system) for index, name in zip(indices, names)]
    edges = [[20, 21], [21, 22], [22, 23], [23, 24], [24, 25], [25, 20], [22, 26]]
    return group(indices, atoms, edges)


def six_ring(system: str) -> dict:
    indices = [30, 31, 32, 33, 34, 35]
    names = ["C2", "O6", "C6", "C5", "C4", "C3"]
    atoms = [atom(index, name, "O" if name.startswith("O") else "C", system) for index, name in zip(indices, names)]
    edges = [[30, 31], [31, 32], [32, 33], [33, 34], [34, 35], [35, 30]]
    return group(indices, atoms, edges)


def five_ring_plus_h() -> dict:
    indices = [40, 41, 42, 43, 44, 45]
    names = ["C2", "O5", "C5", "H62", "C4", "C3"]
    elements = ["C", "O", "C", "H", "C", "C"]
    atoms = [atom(index, name, element, "0CU") for index, name, element in zip(indices, names, elements)]
    edges = [[40, 41], [41, 42], [42, 44], [44, 45], [45, 40]]
    return group(indices, atoms, edges)


class RestraintConsistencySummaryTests(unittest.TestCase):
    def test_two_core_removes_leaf_and_isolated_atom(self) -> None:
        glucose = guest_ring_plus_c6("0GB")
        self.assertEqual(module.two_core(glucose), [20, 21, 22, 23, 24, 25])
        furanose = five_ring_plus_h()
        self.assertEqual(module.two_core(furanose), [40, 41, 42, 44, 45])

    def test_group_classification(self) -> None:
        glucose = module.summarize_group("glucose", guest_ring_plus_c6("0GB"))
        self.assertEqual(glucose["selection_class"], "cyclic_core_plus_exocyclic_heavy_atom")
        self.assertEqual(glucose["noncore_selected_atoms"][0]["atom_name"], "C6")

        pyranose = module.summarize_group("fructopyranose", six_ring("0CB"))
        self.assertEqual(pyranose["selection_class"], "cyclic_heavy_atom_core_only")

        furanose = module.summarize_group("fructofuranose", five_ring_plus_h())
        self.assertEqual(furanose["selection_class"], "cyclic_core_plus_selected_hydrogen_or_mixed_atoms")
        self.assertEqual(furanose["selected_hydrogen_atoms"][0]["atom_name"], "H62")
        self.assertEqual(furanose["isolated_selected_atoms"][0]["atom_name"], "H62")

    def test_full_summary_preserves_failed_consistency_gate(self) -> None:
        host = host_group()
        rows = [
            ("glucose", guest_ring_plus_c6("0GB")),
            ("galactose", guest_ring_plus_c6("0LB")),
            ("fructopyranose", six_ring("0CB")),
            ("fructofuranose", five_ring_plus_h()),
        ]
        mapping = {
            "archive": {"md5": "x", "sha256": "y"},
            "cross_system_summary": {"host_group_mapped_identity_identical": True},
            "mappings": [
                {
                    "system_name": system,
                    "group1": host,
                    "group2": guest,
                    "parameters": {
                        "iat": [-1, -1],
                        "r2": 0.0,
                        "r2a": 15.0,
                        "rk2": 200.0,
                        "iresid": 0.0,
                        "outxyz": 1.0 if system == "fructopyranose" else None,
                    },
                }
                for system, guest in rows
            ],
        }
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "mapping.json"
            output = root / "summary.json"
            source.write_text(json.dumps(mapping), encoding="utf-8")
            report = module.summarize(source, output)
            self.assertFalse(
                report["decision"]
                ["guest_reaction_coordinate_definition_is_uniform_across_sugars"]
            )
            self.assertTrue(
                report["decision"]
                ["fructofuranose_group_contains_unexplained_hydrogen_H62"]
            )
            self.assertTrue(
                report["decision"]
                ["host_reference_group_matches_two_aromatic_six_carbon_rings"]
            )
            self.assertFalse(
                report["decision"]["new_smd_or_umbrella_generation_authorized"]
            )
            self.assertTrue(
                report["parameter_consistency"]["core_parameters_identical"]
            )
            self.assertFalse(
                report["parameter_consistency"]
                ["all_recorded_parameters_identical"]
            )
            self.assertEqual(json.loads(output.read_text()), report)


if __name__ == "__main__":
    unittest.main()
