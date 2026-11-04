from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "research" / "audit_gluhut_compiled_systems.py"
spec = importlib.util.spec_from_file_location("compiled_audit", MODULE_PATH)
assert spec and spec.loader
module = importlib.util.module_from_spec(spec)
sys.modules["compiled_audit"] = module
spec.loader.exec_module(module)


def section(flag: str, kind: str, width: int, count: int, values: list) -> str:
    lines = [f"%FLAG {flag}", f"%FORMAT({count}{kind}{width})"]
    if kind.upper() == "A":
        encoded = [f"{str(value):<{width}}"[:width] for value in values]
    elif kind.upper() == "I":
        encoded = [f"{int(value):>{width}d}" for value in values]
    else:
        encoded = [f"{float(value):>{width}.8E}" for value in values]
    for start in range(0, len(encoded), count):
        lines.append("".join(encoded[start : start + count]))
    return "\n".join(lines)


def minimal_parm7() -> bytes:
    pointers = [0] * len(module.POINTER_NAMES)
    pointer = {name: i for i, name in enumerate(module.POINTER_NAMES)}
    pointers[pointer["NATOM"]] = 3
    pointers[pointer["NTYPES"]] = 2
    pointers[pointer["NRES"]] = 2
    payload = [
        "%VERSION VERSION_STAMP = V0001.000",
        section("POINTERS", "I", 8, 10, pointers),
        section("ATOM_NAME", "A", 4, 20, ["H1", "H2", "G1"]),
        section("CHARGE", "E", 16, 5, [18.2223, -18.2223, 0.0]),
        section("MASS", "E", 16, 5, [12.0, 16.0, 12.0]),
        section("ATOM_TYPE_INDEX", "I", 8, 10, [1, 2, 1]),
        section("NUMBER_EXCLUDED_ATOMS", "I", 8, 10, [0, 0, 0]),
        section("NONBONDED_PARM_INDEX", "I", 8, 10, [1, 2, 2, 3]),
        section("RESIDUE_LABEL", "A", 4, 20, ["HST", "GST"]),
        section("RESIDUE_POINTER", "I", 8, 10, [1, 3]),
        section("LENNARD_JONES_ACOEF", "E", 16, 5, [1.0, 2.0, 3.0]),
        section("LENNARD_JONES_BCOEF", "E", 16, 5, [0.1, 0.2, 0.3]),
        section("AMBER_ATOM_TYPE", "A", 4, 20, ["c3", "oh", "c3"]),
        section("TREE_CHAIN_CLASSIFICATION", "A", 4, 20, ["M", "M", "M"]),
        section("JOIN_ARRAY", "I", 8, 10, [0, 0, 0]),
        section("IROTAT", "I", 8, 10, [0, 0, 0]),
    ]
    return ("\n".join(payload) + "\n").encode()


def host_summary(*, type_indices: list[int], charge_shift: float = 0.0) -> dict:
    physical_props = {
        "atom_names": ["C1", "O1"],
        "charges_e": [0.2 + charge_shift, -0.2],
        "masses": [12.01, 16.0],
        "atomic_numbers": [6, 8],
        "amber_atom_types": ["c3", "oh"],
    }
    indexing = {
        "atom_type_indices": type_indices,
        "tree_chain_classification": ["M", "M"],
        "join_array": [0, 0],
        "irotat": [0, 0],
    }
    physical_components = {
        "atom_properties": module.canonical_hash(physical_props),
        "bonds": module.canonical_hash(
            [{"atoms": [1, 2], "force_constant": 1.0, "equilibrium_A": 1.4}]
        ),
        "angles": module.canonical_hash([]),
        "dihedrals": module.canonical_hash([]),
        "exclusions": module.canonical_hash([[2], []]),
        "host_lj_pairs": module.canonical_hash(
            [{"types": ["c3", "oh"], "acoef": 2.0, "bcoef": 0.2}]
        ),
    }
    indexed_lj = [
        {
            "types": ["c3", "oh"],
            "type_indices": type_indices,
            "nonbonded_parameter_index": 1,
            "acoef": 2.0,
            "bcoef": 0.2,
        }
    ]
    indexing_components = {
        "compiled_indexing_properties": module.canonical_hash(indexing),
        "indexed_host_lj_pairs": module.canonical_hash(indexed_lj),
    }
    return {
        "path": "x",
        "host": {
            "atom_count": 2,
            "label": "HST",
            "charge_e": sum(physical_props["charges_e"]),
            "physical_atom_properties": physical_props,
            "compiled_indexing_properties": indexing,
            "physical_component_hashes": physical_components,
            "indexing_component_hashes": indexing_components,
            "physical_parameter_fingerprint_sha256": module.canonical_hash(
                physical_components
            ),
            "compiled_indexing_fingerprint_sha256": module.canonical_hash(
                indexing_components
            ),
        },
    }


class CompiledSystemAuditTests(unittest.TestCase):
    def test_fortran_format_parser(self):
        self.assertEqual(module.parse_fortran_format("%FORMAT(10I8)"), (10, "I", 8))
        self.assertEqual(module.parse_fortran_format("%FORMAT(5E16.8)"), (5, "E", 16))
        self.assertEqual(module.parse_fortran_format("%FORMAT(20a4)"), (20, "A", 4))

    def test_minimal_topology_summary(self):
        summary = module.topology_summary(
            "fixture/strip.glucose.parm7", minimal_parm7()
        )
        self.assertEqual(summary["natom"], 3)
        self.assertEqual(summary["nres"], 2)
        self.assertEqual(summary["host"]["label"], "HST")
        self.assertEqual(summary["host"]["atom_count"], 2)
        self.assertAlmostEqual(summary["host"]["charge_e"], 0.0, places=12)
        self.assertEqual(summary["guest"]["label"], "GST")
        self.assertEqual(summary["host"]["interaction_counts"]["lj_type_pairs"], 3)

    def test_physical_identity_ignores_global_type_numbering(self):
        left = host_summary(type_indices=[1, 2])
        right = host_summary(type_indices=[7, 9])
        left["path"] = "left"
        right["path"] = "right"
        result = module.compare_hosts([left, right])
        self.assertTrue(
            result["physical_host_parameters_identical_across_all_topologies"]
        )
        self.assertFalse(
            result["compiled_indexing_identical_across_all_topologies"]
        )
        self.assertTrue(result["pairwise"][0]["physical_host_parameters_equal"])

    def test_charge_difference_breaks_physical_identity(self):
        left = host_summary(type_indices=[1, 2])
        right = host_summary(type_indices=[1, 2], charge_shift=1e-4)
        left["path"] = "left"
        right["path"] = "right"
        result = module.compare_hosts([left, right])
        self.assertFalse(
            result["physical_host_parameters_identical_across_all_topologies"]
        )
        self.assertGreater(
            result["pairwise"][0]["maximum_absolute_charge_difference_e"], 0.0
        )

    def test_restart_parser(self):
        payload = b"title\n       2  1.2500000\n 0 0 0 1 1 1\n 10 10 10\n"
        result = module.parse_restart("fixture.rst7", payload)
        self.assertEqual(result["natom"], 2)
        self.assertFalse(result["velocities_present"])
        self.assertEqual(result["box_value_count"], 3)
        self.assertEqual(result["time"], 1.25)

    def test_ambiguous_rst_is_classified_by_content(self):
        restraint = b"&rst\n iat=-1,-1,\n r1=0.0, r2=1.0, r3=2.0, r4=3.0,\n/\n"
        result = module.classify_ambiguous_rst("distance.rst", restraint)
        self.assertEqual(result["kind"], "restraint_input")
        self.assertTrue(result["contains_iat"])
        self.assertTrue(result["contains_rst_namelist"])

        restart = b"title\n       2\n 0 0 0 1 1 1\n"
        result = module.classify_ambiguous_rst("coordinate.rst", restart)
        self.assertEqual(result["kind"], "restart_coordinate")
        self.assertEqual(result["restart"]["natom"], 2)

    def test_unsafe_paths(self):
        for value in ("../x", "/x", "a\\b"):
            with self.assertRaises(module.AuditError):
                module.safe_member(value)


if __name__ == "__main__":
    unittest.main()
