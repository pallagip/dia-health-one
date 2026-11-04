from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "research" / "export_gluhut_host_graph.py"
spec = importlib.util.spec_from_file_location("host_export", MODULE_PATH)
assert spec and spec.loader
module = importlib.util.module_from_spec(spec)
sys.modules["host_export"] = module
spec.loader.exec_module(module)


class GluHUTHostGraphExportTests(unittest.TestCase):
    def setUp(self) -> None:
        self.nodes = [
            {
                "compiled_atom_index_1_based": 1,
                "host_atom_index_1_based": 1,
                "atom_name": "C1",
                "element": "C",
                "atomic_number": 6,
                "amber_atom_type": "ca",
                "charge_e": 0.1,
                "mass": 12.01,
                "total_degree": 3,
                "heavy_degree": 2,
                "heavy_neighbour_host_indices_1_based": [2, 3],
            },
            {
                "compiled_atom_index_1_based": 2,
                "host_atom_index_1_based": 2,
                "atom_name": "N1",
                "element": "N",
                "atomic_number": 7,
                "amber_atom_type": "n",
                "charge_e": -0.1,
                "mass": 14.01,
                "total_degree": 2,
                "heavy_degree": 1,
                "heavy_neighbour_host_indices_1_based": [1],
            },
        ]
        self.edges = [
            {
                "begin_compiled_atom_index_1_based": 1,
                "end_compiled_atom_index_1_based": 2,
                "begin_host_atom_index_1_based": 1,
                "end_host_atom_index_1_based": 2,
                "force_constant": 300.0,
                "equilibrium_A": 1.4,
                "source_flag": "BONDS_WITHOUT_HYDROGEN",
                "chemical_bond_order_known": False,
            }
        ]

    def test_structural_payload_excludes_compiled_parameters(self) -> None:
        payload = module.structural_node_payload(self.nodes)
        self.assertEqual(payload[0]["element"], "C")
        self.assertNotIn("charge_e", payload[0])
        self.assertNotIn("amber_atom_type", payload[0])
        self.assertEqual(
            payload[0]["heavy_neighbour_host_indices_1_based"],
            [2, 3],
        )

    def test_compiled_payload_retains_charge_type_and_terms(self) -> None:
        node_payload = module.compiled_node_payload(self.nodes)
        edge_payload = module.compiled_edge_payload(self.edges)
        self.assertEqual(node_payload[0]["charge_e"], 0.1)
        self.assertEqual(node_payload[0]["amber_atom_type"], "ca")
        self.assertEqual(edge_payload[0]["force_constant"], 300.0)
        self.assertEqual(edge_payload[0]["equilibrium_A"], 1.4)

    def test_structural_fingerprint_ignores_charge_changes(self) -> None:
        changed = [dict(row) for row in self.nodes]
        changed[0]["charge_e"] = 0.9
        first = module.canonical_hash(
            {
                "nodes": module.structural_node_payload(self.nodes),
                "edges": module.structural_edge_payload(self.edges),
            }
        )
        second = module.canonical_hash(
            {
                "nodes": module.structural_node_payload(changed),
                "edges": module.structural_edge_payload(self.edges),
            }
        )
        self.assertEqual(first, second)
        compiled_first = module.canonical_hash(
            {
                "nodes": module.compiled_node_payload(self.nodes),
                "edges": module.compiled_edge_payload(self.edges),
            }
        )
        compiled_second = module.canonical_hash(
            {
                "nodes": module.compiled_node_payload(changed),
                "edges": module.compiled_edge_payload(self.edges),
            }
        )
        self.assertNotEqual(compiled_first, compiled_second)

    def test_safe_member_rejects_traversal(self) -> None:
        for value in ("../x", "/x", "a\\b"):
            with self.assertRaises(module.HostGraphExportError):
                module.safe_member(value)

    def test_system_inference(self) -> None:
        self.assertEqual(
            module.infer_system("SMD/Glucose/strip.glucose.parm7"),
            "glucose",
        )
        self.assertEqual(
            module.infer_system("SMD/Fructofuranose/host.parm7"),
            "fructofuranose",
        )


if __name__ == "__main__":
    unittest.main()
