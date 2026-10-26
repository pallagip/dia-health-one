#!/usr/bin/env python3
"""Enumerate explicit group-level charge states for full-conjugate model prep.

The audit starts from the neutral constitutional graph. It applies declared
proton-count assumptions at the functional-group level; it does not predict
pKa values, microstate populations, atom-localized partial charges, or force
field parameters.
"""

from __future__ import annotations

import hashlib
import itertools
import json
from collections import Counter
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
INPUT_PATH = ROOT / "research/data/nnc2215-full-conjugate-graph.json"
OUTPUT_PATH = ROOT / "research/data/nnc2215-physiological-charge-states.json"
EXPECTED_NEUTRAL_FORMULA = "C325H461N79O97S6"
EXPECTED_CANONICAL_FORMULA = "C325H454N79O97S6"
MACROCYCLE_GROUP_MAPS = ([52, 53, 54], [60, 61, 62], [68, 69, 70])
HISTIDINE_RING_ATOMS = ("CG", "ND1", "CD2", "CE1", "NE2")
EXPECTED_CHAIN_A = "GIVEQCCTSICSLYQLENYCN"
EXPECTED_CHAIN_B = "FVNQHLCGSHLVEALYLVCGERGFFYTPK"


def fail(message: str) -> None:
    raise SystemExit(f"charge-state audit failed: {message}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        fail(message)


def atom_provenance(atom: dict[str, Any]) -> dict[str, Any]:
    value = atom.get("provenance")
    require(isinstance(value, dict), f"{atom.get('node_id')} has no provenance object")
    return value


def template_atoms(
    atoms: list[dict[str, Any]], chain: str, residue_number: int, residue_name: str
) -> list[dict[str, Any]]:
    return [
        atom
        for atom in atoms
        if atom_provenance(atom).get("origin") == "insulin_chain_template"
        and atom_provenance(atom).get("chain") == chain
        and atom_provenance(atom).get("residue_number") == residue_number
        and atom_provenance(atom).get("residue_name") == residue_name
    ]


def template_atom(
    atoms: list[dict[str, Any]],
    chain: str,
    residue_number: int,
    residue_name: str,
    atom_name: str,
) -> dict[str, Any]:
    matches = [
        atom
        for atom in template_atoms(atoms, chain, residue_number, residue_name)
        if atom_provenance(atom).get("atom_name") == atom_name
    ]
    require(
        len(matches) == 1,
        f"expected one {chain}:{residue_number} {residue_name} {atom_name}, found {len(matches)}",
    )
    return matches[0]


def group_record(
    site_id: str,
    label: str,
    site_type: str,
    atoms: list[dict[str, Any]],
    *,
    names: list[str] | None = None,
    neutral_hydrogen_count: int,
    source_map_ids: list[int] | None = None,
    states: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    return {
        "site_id": site_id,
        "label": label,
        "site_type": site_type,
        "atom_node_ids": [atom["node_id"] for atom in atoms],
        "atom_names": names,
        "source_map_ids": source_map_ids,
        "neutral_graph_site_hydrogen_count": neutral_hydrogen_count,
        "state_options": states,
    }


def formula_with_hydrogen_count(base_counts: dict[str, int], hydrogen_count: int) -> str:
    order = ("C", "H", "N", "O", "S", "P", "F", "Cl", "Br", "I")
    pieces: list[str] = []
    for element in order:
        count = hydrogen_count if element == "H" else base_counts.get(element, 0)
        if count:
            pieces.append(element)
            if count != 1:
                pieces.append(str(count))
    for element in sorted(set(base_counts) - set(order)):
        count = base_counts[element]
        if count:
            pieces.extend((element, str(count) if count != 1 else ""))
    return "".join(pieces)


def bond_index(bonds: list[dict[str, Any]], atom_by_id: dict[str, dict[str, Any]]) -> dict[frozenset[str], dict[str, Any]]:
    indexed: dict[frozenset[str], dict[str, Any]] = {}
    for bond in bonds:
        endpoints = bond.get("atom_node_ids", [])
        require(len(endpoints) == 2 and all(node_id in atom_by_id for node_id in endpoints),
                f"invalid endpoints for {bond.get('edge_id', 'unnamed bond')}")
        key = frozenset(endpoints)
        require(key not in indexed, f"duplicate bond between {sorted(key)}")
        indexed[key] = bond
    return indexed


def require_bond(
    indexed_bonds: dict[frozenset[str], dict[str, Any]],
    left: dict[str, Any],
    right: dict[str, Any],
    expected_order: float,
    context: str,
) -> dict[str, Any]:
    bond = indexed_bonds.get(frozenset((left["node_id"], right["node_id"])))
    require(bond is not None, f"{context} bond is missing")
    require(bond.get("bond_order") == expected_order, f"{context} bond order changed")
    return bond


def heavy_neighbors(
    atom: dict[str, Any], bonds: list[dict[str, Any]], atom_by_id: dict[str, dict[str, Any]]
) -> list[dict[str, Any]]:
    neighbors: list[dict[str, Any]] = []
    for bond in bonds:
        endpoints = bond["atom_node_ids"]
        if atom["node_id"] in endpoints:
            other_id = endpoints[0] if endpoints[1] == atom["node_id"] else endpoints[1]
            if atom_by_id[other_id]["element"] != "H":
                neighbors.append(atom_by_id[other_id])
    return neighbors


def require_neutral_carboxyl(
    members: list[dict[str, Any]],
    names: list[str],
    indexed_bonds: dict[frozenset[str], dict[str, Any]],
    label: str,
    *,
    anchor: dict[str, Any] | None = None,
) -> None:
    require([atom["element"] for atom in members] == ["C", "O", "O"], f"{label} elements changed")
    require(sum(atom["formal_charge"] for atom in members) == 0, f"{label} neutral graph charge changed")
    require(sum(atom["hydrogen_count_not_as_nodes"] for atom in members) == 1,
            f"{label} must contain one neutral-acid proton")
    carbon, oxygen_1, oxygen_2 = members
    bond_1 = indexed_bonds.get(frozenset((carbon["node_id"], oxygen_1["node_id"])))
    bond_2 = indexed_bonds.get(frozenset((carbon["node_id"], oxygen_2["node_id"])))
    require(bond_1 is not None and bond_2 is not None, f"{label} C–O bond is missing")
    order_1 = bond_1["bond_order"]
    order_2 = bond_2["bond_order"]
    require(sorted((order_1, order_2)) == [1.0, 2.0], f"{label} C–O bond orders changed")
    h_1 = oxygen_1["hydrogen_count_not_as_nodes"]
    h_2 = oxygen_2["hydrogen_count_not_as_nodes"]
    require((h_1, h_2) in ((1, 0), (0, 1)), f"{label} proton is not assigned to one oxygen")
    require((order_1 if h_1 else order_2) == 1.0,
            f"{label} hydroxyl oxygen is not attached by a single bond")
    carbon_bonds = [bond for bond in indexed_bonds.values() if carbon["node_id"] in bond["atom_node_ids"]]
    require(len(carbon_bonds) == 3, f"{label} carboxyl carbon has unexpected heavy-atom neighbors")
    for oxygen in (oxygen_1, oxygen_2):
        oxygen_bonds = [bond for bond in indexed_bonds.values() if oxygen["node_id"] in bond["atom_node_ids"]]
        require(len(oxygen_bonds) == 1, f"{label} oxygen has an unexpected heavy-atom neighbor")
    if anchor is not None:
        require_bond(indexed_bonds, anchor, carbon, 1.0, f"{label} residue anchor")


def mapped_source_atom(
    atoms: list[dict[str, Any]], component: str, map_id: int
) -> dict[str, Any]:
    matches = [
        atom for atom in atoms
        if atom_provenance(atom).get("origin") == "mapped_source_fragment"
        and atom_provenance(atom).get("component") == component
        and atom_provenance(atom).get("source_map_id") == map_id
    ]
    require(len(matches) == 1, f"expected one {component} map {map_id}, found {len(matches)}")
    return matches[0]


def main() -> None:
    report = json.loads(INPUT_PATH.read_text(encoding="utf-8"))
    neutral_formula = report.get("formula", {}).get("atom_ledger_formula")
    require(neutral_formula == EXPECTED_NEUTRAL_FORMULA, "neutral input formula changed")
    require(report.get("schema") == "dia-health-one.nnc2215-full-conjugate-graph.v1", "unexpected input graph schema")

    atoms: list[dict[str, Any]] = report.get("atoms", [])
    require(len(atoms) == 507, f"expected 507 graph atoms, found {len(atoms)}")
    atom_by_id = {atom["node_id"]: atom for atom in atoms}
    require(len(atom_by_id) == len(atoms), "duplicate graph node ID")
    bonds: list[dict[str, Any]] = report.get("bonds", [])
    indexed_bonds = bond_index(bonds, atom_by_id)
    require(len(bonds) == 527, f"expected 527 graph bonds, found {len(bonds)}")
    actual_element_counts = Counter(atom["element"] for atom in atoms)
    recorded_element_counts = report.get("formula", {}).get("element_counts", {})
    require(actual_element_counts.get("H", 0) == 0,
            "neutral graph unexpectedly contains explicit hydrogen nodes")
    require(
        {element: count for element, count in actual_element_counts.items() if element != "H"}
        == {element: count for element, count in recorded_element_counts.items() if element != "H"},
        "neutral graph heavy-element counts do not match the formula ledger",
    )
    require(report.get("sequences", {}).get("chain_A") == EXPECTED_CHAIN_A,
            "A-chain insulin sequence changed")
    require(report.get("sequences", {}).get("chain_B_desB30") == EXPECTED_CHAIN_B,
            "B-chain desB30 insulin sequence changed")
    require(report.get("sequences", {}).get("sequence_audit_match") is True,
            "upstream sequence audit is not passing")
    neutral_h = report["formula"]["element_counts"]["H"]
    require(
        sum(atom.get("hydrogen_count_not_as_nodes", 0) for atom in atoms) == neutral_h,
        "neutral graph per-atom hydrogen ledger does not match its formula",
    )
    require(
        sum(atom.get("formal_charge", 0) for atom in atoms) == 0,
        "neutral graph per-atom formal charges are not neutral",
    )
    require(all(atom.get("formal_charge", 0) == 0 for atom in atoms),
            "neutral graph contains a nonzero atom formal charge")

    fixed_groups: list[dict[str, Any]] = []
    variable_groups: list[dict[str, Any]] = []

    # Four glutamate sidechains are represented as neutral carboxylic acids in
    # the source graph; the pH 7.4 baseline treats each as a carboxylate.
    glutamates = [("A", 4), ("A", 17), ("B", 13), ("B", 21)]
    for chain, residue_number in glutamates:
        residue = template_atoms(atoms, chain, residue_number, "GLU")
        atom_names = ["CD", "OE1", "OE2"]
        members = [
            template_atom(atoms, chain, residue_number, "GLU", name)
            for name in atom_names
        ]
        require(len(residue) >= len(members), f"incomplete glutamate {chain}:{residue_number}")
        anchor = template_atom(atoms, chain, residue_number, "GLU", "CG")
        require_neutral_carboxyl(
            members, atom_names, indexed_bonds, f"Glu {chain}{residue_number}", anchor=anchor
        )
        fixed_groups.append(
            group_record(
                f"glu_{chain}{residue_number}",
                f"Glu {chain}{residue_number} sidechain carboxyl",
                "fixed_carboxylate_assumption",
                members,
                names=atom_names,
                neutral_hydrogen_count=1,
                states=[{"state": "deprotonated", "group_formal_charge": -1, "hydrogen_delta": -1}],
            )
        )

    # The two peptide-chain C-terminal acids remain separate from the three
    # macrocycle acids and each loses one proton in the conventional baseline.
    termini = [("A", 21, "ASN"), ("B", 29, "LYS")]
    for chain, residue_number, residue_name in termini:
        atom_names = ["C", "O", "OXT"]
        members = [
            template_atom(atoms, chain, residue_number, residue_name, name)
            for name in atom_names
        ]
        anchor = template_atom(atoms, chain, residue_number, residue_name, "CA")
        require_neutral_carboxyl(
            members, atom_names, indexed_bonds,
            f"C-terminal group {chain}:{residue_number}", anchor=anchor,
        )
        fixed_groups.append(
            group_record(
                f"cterm_{chain}{residue_number}",
                f"Chain {chain} C-terminal carboxyl",
                "fixed_carboxylate_assumption",
                members,
                names=atom_names,
                neutral_hydrogen_count=1,
                states=[{"state": "deprotonated", "group_formal_charge": -1, "hydrogen_delta": -1}],
            )
        )

    # Each mapped macrocycle acid is recorded as one three-map functional group.
    # The ledger deliberately avoids naming a particular oxygen as the anion.
    for index, map_ids in enumerate(MACROCYCLE_GROUP_MAPS, start=1):
        members = [
            atom
            for atom in atoms
            if atom_provenance(atom).get("origin") == "mapped_source_fragment"
            and atom_provenance(atom).get("component") == "B29"
            and atom_provenance(atom).get("source_map_id") in map_ids
        ]
        found_maps = sorted(atom_provenance(atom)["source_map_id"] for atom in members)
        require(len(members) == 3 and found_maps == sorted(map_ids),
                f"macrocycle acid {index} map membership changed: {found_maps}")
        mapped_members = [mapped_source_atom(atoms, "B29", map_id) for map_id in map_ids]
        require_neutral_carboxyl(
            mapped_members, [f"map_{map_id}" for map_id in map_ids], indexed_bonds,
            f"B29 macrocycle carboxyl {index}",
        )
        variable_groups.append(
            group_record(
                f"macrocycle_carboxyl_{index}",
                f"B29 macrocycle carboxyl {index}",
                "variable_carboxyl_group",
                mapped_members,
                source_map_ids=list(map_ids),
                neutral_hydrogen_count=1,
                states=[
                    {"state": "deprotonated", "group_formal_charge": -1, "hydrogen_delta": -1},
                    {"state": "protonated", "group_formal_charge": 0, "hydrogen_delta": 0},
                ],
            )
        )

    # One variable N-terminal group occurs at GlyA1. PheB1 is conjugated and is
    # not treated as a free N terminus.
    a1_n = template_atom(atoms, "A", 1, "GLY", "N")
    require(
        a1_n["element"] == "N"
        and a1_n["hydrogen_count_not_as_nodes"] == 2
        and a1_n["formal_charge"] == 0,
        "neutral GlyA1 N-terminus representation changed",
    )
    a1_ca = template_atom(atoms, "A", 1, "GLY", "CA")
    require_bond(indexed_bonds, a1_n, a1_ca, 1.0, "GlyA1 N–CA")
    require(len(heavy_neighbors(a1_n, bonds, atom_by_id)) == 1,
            "GlyA1 N terminus has an unexpected heavy-atom attachment")
    variable_groups.append(
        group_record(
            "nterm_A1",
            "A-chain GlyA1 N terminus",
            "variable_ammonium_group",
            [a1_n],
            names=["N"],
            neutral_hydrogen_count=a1_n["hydrogen_count_not_as_nodes"],
            states=[
                {"state": "protonated", "group_formal_charge": 1, "hydrogen_delta": 1},
                {"state": "neutral", "group_formal_charge": 0, "hydrogen_delta": 0},
            ],
        )
    )

    # The standard pH 7.4 baseline treats ArgB22 guanidinium as protonated.
    arg_residue = template_atoms(atoms, "B", 22, "ARG")
    require(arg_residue, "ArgB22 residue not found")
    arg_names = [
        atom_provenance(atom).get("atom_name")
        for atom in arg_residue
        if atom_provenance(atom).get("atom_name") in {"NE", "CZ", "NH1", "NH2"}
    ]
    require(set(arg_names) == {"NE", "CZ", "NH1", "NH2"}, "ArgB22 guanidinium atom identity changed")
    arg_members = [
        template_atom(atoms, "B", 22, "ARG", name)
        for name in ("NE", "CZ", "NH1", "NH2")
    ]
    require([atom["element"] for atom in arg_members] == ["N", "C", "N", "N"],
            "ArgB22 guanidinium elements changed")
    require(sum(atom["formal_charge"] for atom in arg_members) == 0,
            "neutral ArgB22 graph group charge changed")
    require(sum(atom["hydrogen_count_not_as_nodes"] for atom in arg_members) == 4,
            "neutral ArgB22 graph hydrogen count changed")
    for nitrogen_name in ("NE", "NH1", "NH2"):
        nitrogen = template_atom(atoms, "B", 22, "ARG", nitrogen_name)
        cz = template_atom(atoms, "B", 22, "ARG", "CZ")
        require(indexed_bonds.get(frozenset((nitrogen["node_id"], cz["node_id"]))) is not None,
                f"ArgB22 {nitrogen_name}–CZ bond is missing")
    require(sorted(
        indexed_bonds[frozenset((template_atom(atoms, "B", 22, "ARG", name)["node_id"],
                                 template_atom(atoms, "B", 22, "ARG", "CZ")["node_id"]))]["bond_order"]
        for name in ("NE", "NH1", "NH2")
    ) == [1.0, 1.0, 2.0], "ArgB22 guanidinium bond-order pattern changed")
    cz = template_atom(atoms, "B", 22, "ARG", "CZ")
    require({atom["node_id"] for atom in heavy_neighbors(cz, bonds, atom_by_id)}
            == {atom["node_id"] for atom in (arg_members[0], arg_members[2], arg_members[3])},
            "ArgB22 guanidinium carbon has unexpected heavy-atom neighbors")
    fixed_groups.append(
        group_record(
            "arg_B22",
            "ArgB22 guanidinium",
            "fixed_protonated_base_assumption",
            arg_members,
            names=["NE", "CZ", "NH1", "NH2"],
            neutral_hydrogen_count=sum(atom["hydrogen_count_not_as_nodes"] for atom in arg_members),
            states=[{"state": "protonated", "group_formal_charge": 1, "hydrogen_delta": 1}],
        )
    )

    # Histidine neutral tautomers and protonated states are enumerated without
    # assigning populations or selecting a tautomer from coordinates.
    histidine_states: dict[str, list[str]] = {}
    for residue_number in (5, 10):
        residue = template_atoms(atoms, "B", residue_number, "HIS")
        require(residue, f"HisB{residue_number} residue not found")
        ring_members = [
            template_atom(atoms, "B", residue_number, "HIS", name)
            for name in HISTIDINE_RING_ATOMS
        ]
        nd1 = template_atom(atoms, "B", residue_number, "HIS", "ND1")
        ne2 = template_atom(atoms, "B", residue_number, "HIS", "NE2")
        require([atom["element"] for atom in ring_members] == ["C", "N", "C", "C", "N"],
                f"HisB{residue_number} ring elements changed")
        nd1_h = nd1["hydrogen_count_not_as_nodes"]
        ne2_h = ne2["hydrogen_count_not_as_nodes"]
        require(nd1_h + ne2_h == 1, f"HisB{residue_number} neutral template tautomer is not unambiguous")
        require(
            sum(atom["formal_charge"] for atom in ring_members) == 0,
            f"HisB{residue_number} neutral template ring has nonzero formal charge",
        )
        require(sum(atom["hydrogen_count_not_as_nodes"] for atom in ring_members) == 3,
                f"HisB{residue_number} neutral ring hydrogen count changed")
        ring_bond_names = (("CG", "ND1"), ("ND1", "CE1"), ("CE1", "NE2"),
                           ("NE2", "CD2"), ("CD2", "CG"))
        ring_edges = [
            indexed_bonds.get(frozenset((
                template_atom(atoms, "B", residue_number, "HIS", left)["node_id"],
                template_atom(atoms, "B", residue_number, "HIS", right)["node_id"],
            )))
            for left, right in ring_bond_names
        ]
        require(all(edge is not None for edge in ring_edges),
                f"HisB{residue_number} imidazole ring connectivity changed")
        require(all(edge["bond_type"] == "AROMATIC" and edge["bond_order"] == 1.5
                    and edge["is_aromatic"] is True for edge in ring_edges),
                f"HisB{residue_number} imidazole ring bond labels changed")
        require(len(heavy_neighbors(nd1, bonds, atom_by_id)) == 2
                and len(heavy_neighbors(ne2, bonds, atom_by_id)) == 2,
                f"HisB{residue_number} ring nitrogen attachment changed")
        template_tautomer = "delta1" if nd1_h == 1 else "epsilon2"
        alternate_tautomer = "epsilon2" if template_tautomer == "delta1" else "delta1"
        histidine_states[f"his_B{residue_number}"] = [
            "neutral_template_tautomer",
            "neutral_alternate_tautomer",
            "protonated",
        ]
        variable_groups.append(
            group_record(
                f"his_B{residue_number}",
                f"HisB{residue_number} imidazole",
                "variable_histidine_group",
                ring_members,
                names=list(HISTIDINE_RING_ATOMS),
                neutral_hydrogen_count=sum(atom["hydrogen_count_not_as_nodes"] for atom in ring_members),
                states=[
                    {
                        "state": "neutral_template_tautomer",
                        "tautomer": template_tautomer,
                        "group_formal_charge": 0,
                        "hydrogen_delta": 0,
                    },
                    {
                        "state": "neutral_alternate_tautomer",
                        "tautomer": alternate_tautomer,
                        "group_formal_charge": 0,
                        "hydrogen_delta": 0,
                    },
                    {"state": "protonated", "tautomer": "both_ring_nitrogens_protonated", "group_formal_charge": 1, "hydrogen_delta": 1},
                ],
            )
        )

    # Record high-pKa tyrosine phenols as excluded from this baseline rather
    # than silently overlooking them.
    excluded_sites: list[dict[str, Any]] = []
    for chain, residue_number in (("A", 14), ("A", 19), ("B", 16), ("B", 26)):
        oh = template_atom(atoms, chain, residue_number, "TYR", "OH")
        require(oh["element"] == "O" and oh["formal_charge"] == 0
                and oh["hydrogen_count_not_as_nodes"] == 1,
                f"Tyr {chain}{residue_number} phenol neutral graph state changed")
        excluded_sites.append(
            {
                "site_id": f"tyr_{chain}{residue_number}_phenol",
                "label": f"Tyr {chain}{residue_number} phenol",
                "atom_node_id": oh["node_id"],
                "assumed_state": "neutral",
                "reason": "Excluded from the conventional pH 7.4 baseline; no local pKa prediction or measurement is available.",
            }
        )

    # Modified termini are recorded explicitly as neutral amide/carbamate groups.
    # Validate both the serialized join operation and its actual graph edges.
    expected_joins = {
        "B1": {
            "removed_cap_map_id": 1,
            "source_endpoint_map_id": 2,
            "protein_endpoint": {"chain": "B", "residue_number": 1, "residue_name": "PHE", "atom_name": "N"},
            "carbonyl_oxygen_map_id": 3,
            "carbonyl_oxygen_element": "O",
            "other_substituent_map_id": 4,
            "other_substituent_element": "C",
            "edge_provenance": "B1 protein attachment boundary; source cap replaced at map 2",
        },
        "B29": {
            "removed_cap_map_id": 92,
            "source_endpoint_map_id": 84,
            "protein_endpoint": {"chain": "B", "residue_number": 29, "residue_name": "LYS", "atom_name": "NZ"},
            "carbonyl_oxygen_map_id": 90,
            "carbonyl_oxygen_element": "O",
            "other_substituent_map_id": 91,
            "other_substituent_element": "O",
            "edge_provenance": "B29 protein attachment boundary; source cap replaced at map 84",
        },
    }
    join_operations = report.get("assembly", {}).get("join_operations", [])
    require(len(join_operations) == len(expected_joins), "expected the B1 and B29 conjugate join operations")
    joins_by_component = {join.get("component"): join for join in join_operations}
    require(set(joins_by_component) == set(expected_joins), "B1/B29 join component identities changed")
    endpoint_atoms: dict[str, dict[str, Any]] = {}
    for component, expectation in expected_joins.items():
        join = joins_by_component[component]
        require(join.get("removed_cap_map_id") == expectation["removed_cap_map_id"],
                f"{component} removed-cap map changed")
        require(join.get("source_endpoint_map_id") == expectation["source_endpoint_map_id"],
                f"{component} source endpoint map changed")
        require(not any(
            atom_provenance(atom).get("origin") == "mapped_source_fragment"
            and atom_provenance(atom).get("component") == component
            and atom_provenance(atom).get("source_map_id") == expectation["removed_cap_map_id"]
            for atom in atoms
        ), f"{component} removed cap map is still present in the full graph")
        require(join.get("protein_endpoint") == expectation["protein_endpoint"],
                f"{component} protein attachment identity changed")
        require(join.get("bond_order") == "single", f"{component} join bond order changed")
        endpoint = template_atom(
            atoms, expectation["protein_endpoint"]["chain"],
            expectation["protein_endpoint"]["residue_number"],
            expectation["protein_endpoint"]["residue_name"],
            expectation["protein_endpoint"]["atom_name"],
        )
        source = mapped_source_atom(atoms, component, expectation["source_endpoint_map_id"])
        carbonyl_oxygen = mapped_source_atom(atoms, component, expectation["carbonyl_oxygen_map_id"])
        substituent = mapped_source_atom(atoms, component, expectation["other_substituent_map_id"])
        require(
            endpoint["element"] == "N"
            and source["element"] == "C"
            and carbonyl_oxygen["element"] == expectation["carbonyl_oxygen_element"]
            and substituent["element"] == expectation["other_substituent_element"],
            f"{component} attachment atom elements changed",
        )
        require(endpoint["formal_charge"] == 0 and endpoint["hydrogen_count_not_as_nodes"] == 1,
                f"{component} conjugated nitrogen neutral amide/carbamate representation changed")
        require_bond(indexed_bonds, endpoint, source, 1.0, f"{component} protein–conjugate join")
        require_bond(
            indexed_bonds, endpoint,
            template_atom(atoms, "B", 1 if component == "B1" else 29,
                         "PHE" if component == "B1" else "LYS",
                         "CA" if component == "B1" else "CE"),
            1.0, f"{component} native residue attachment",
        )
        boundary_bond = indexed_bonds[frozenset((endpoint["node_id"], source["node_id"]))]
        require(boundary_bond.get("provenance") == expectation["edge_provenance"],
                f"{component} graph attachment-edge provenance changed")
        require_bond(indexed_bonds, source, carbonyl_oxygen, 2.0, f"{component} amide/carbamate carbonyl")
        require_bond(indexed_bonds, source, substituent, 1.0, f"{component} carbonyl substituent")
        expected_neighbors = {endpoint["node_id"], carbonyl_oxygen["node_id"], substituent["node_id"]}
        require({atom["node_id"] for atom in heavy_neighbors(source, bonds, atom_by_id)} == expected_neighbors,
                f"{component} attachment carbonyl has unexpected heavy-atom neighbors")
        require({atom["node_id"] for atom in heavy_neighbors(endpoint, bonds, atom_by_id)} == {
            source["node_id"],
            template_atom(atoms, "B", 1 if component == "B1" else 29,
                          "PHE" if component == "B1" else "LYS",
                          "CA" if component == "B1" else "CE")["node_id"],
        }, f"{component} modified nitrogen has unexpected heavy-atom neighbors")
        endpoint_atoms[component] = endpoint
    phe_b1_n = endpoint_atoms["B1"]
    lys_b29_nz = endpoint_atoms["B29"]
    endpoint_group_notes = [
        {
            "site_id": "conjugated_PheB1_N",
            "atom_node_id": phe_b1_n["node_id"],
            "graph_formal_charge": 0,
            "interpretation": "Acylated B-chain N terminus; not a free amine in this graph.",
            "hydrogen_ownership": "Unresolved in the candidate PDB records; this ledger does not change it.",
        },
        {
            "site_id": "conjugated_LysB29_NZ",
            "atom_node_id": lys_b29_nz["node_id"],
            "graph_formal_charge": 0,
            "interpretation": "Carbamate-bound Lys epsilon nitrogen; not a free lysine amine in this graph.",
        },
    ]

    require(len(fixed_groups) == 7, f"expected six fixed carboxylates and ArgB22, found {len(fixed_groups)} fixed groups")
    require(len(variable_groups) == 6, f"expected three macrocycle acids, A1 amine, and two histidines, found {len(variable_groups)} variable groups")

    # The scope excludes free thiols only if all six cysteines are paired in
    # the declared insulin disulfides and those bonds exist in the graph.
    expected_disulfides = {
        frozenset((('A', 6, 'SG'), ('A', 11, 'SG'))),
        frozenset((('A', 7, 'SG'), ('B', 7, 'SG'))),
        frozenset((('A', 20, 'SG'), ('B', 19, 'SG'))),
    }
    disulfide_records = report.get("assembly", {}).get("disulfides", [])
    observed_disulfides = {
        frozenset((endpoint.get("chain"), endpoint.get("residue_number"), endpoint.get("atom_name"))
                  for endpoint in record.get("atom_identities", []))
        for record in disulfide_records
    }
    require(observed_disulfides == expected_disulfides and len(disulfide_records) == 3,
            "declared insulin disulfide identities changed")
    for pair in expected_disulfides:
        pair_atoms = [template_atom(atoms, chain, residue_number, "CYS", atom_name)
                      for chain, residue_number, atom_name in sorted(pair)]
        require(all(atom["element"] == "S" and atom["formal_charge"] == 0
                    and atom["hydrogen_count_not_as_nodes"] == 0 for atom in pair_atoms),
                "disulfide cysteine sulfur identity or neutral graph state changed")
        require_bond(indexed_bonds, pair_atoms[0], pair_atoms[1], 1.0,
                     "declared cysteine disulfide")

    # Group-level charge/hydrogen deltas relative to the source-neutral graph.
    fixed_charge_delta = sum(group["state_options"][0]["group_formal_charge"] for group in fixed_groups)
    fixed_hydrogen_delta = sum(group["state_options"][0]["hydrogen_delta"] for group in fixed_groups)
    require(fixed_charge_delta == -5 and fixed_hydrogen_delta == -5, "fixed-group charge/hydrogen arithmetic changed")

    macro_ids = [f"macrocycle_carboxyl_{index}" for index in range(1, 4)]
    his_ids = ["his_B5", "his_B10"]
    variable_by_id = {group["site_id"]: group for group in variable_groups}
    a1_states = variable_by_id["nterm_A1"]["state_options"]
    macro_state_options = [
        variable_by_id[site_id]["state_options"] for site_id in macro_ids
    ]
    histidine_options = [
        variable_by_id[site_id]["state_options"] for site_id in his_ids
    ]

    base_counts = report["formula"]["element_counts"]
    enumerated_states: list[dict[str, Any]] = []
    for index, combination in enumerate(
        itertools.product(*macro_state_options, a1_states, *histidine_options), start=1
    ):
        macro_states = combination[:3]
        a1_state = combination[3]
        his_states = combination[4:]
        state_charge = fixed_charge_delta + sum(
            state["group_formal_charge"]
            for state in (*macro_states, a1_state, *his_states)
        )
        hydrogen_delta = fixed_hydrogen_delta + sum(
            state["hydrogen_delta"] for state in (*macro_states, a1_state, *his_states)
        )
        hydrogen_count = neutral_h + hydrogen_delta
        state_id = "_".join(
            [
                "mc" + "".join("d" if state["state"] == "deprotonated" else "p" for state in macro_states),
                "a1p" if a1_state["state"] == "protonated" else "a1n",
                "h5" + {"neutral_template_tautomer": "t", "neutral_alternate_tautomer": "a", "protonated": "p"}[his_states[0]["state"]],
                "h10" + {"neutral_template_tautomer": "t", "neutral_alternate_tautomer": "a", "protonated": "p"}[his_states[1]["state"]],
            ]
        )
        enumerated_states.append(
            {
                "state_index": index,
                "state_id": state_id,
                "variable_site_states": {
                    **{site_id: state["state"] for site_id, state in zip(macro_ids, macro_states)},
                    "nterm_A1": a1_state["state"],
                    **{site_id: state["state"] for site_id, state in zip(his_ids, his_states)},
                },
                "histidine_tautomers": {
                    site_id: state.get("tautomer")
                    for site_id, state in zip(his_ids, his_states)
                },
                "net_formal_charge": state_charge,
                "hydrogen_delta_from_neutral_graph": hydrogen_delta,
                "formula": formula_with_hydrogen_count(base_counts, hydrogen_count),
                "hydrogen_count": hydrogen_count,
                "probability_assigned": False,
            }
        )

    require(len(enumerated_states) == 144, f"expected 144 combinations, found {len(enumerated_states)}")
    state_ids = [state["state_id"] for state in enumerated_states]
    require(len(set(state_ids)) == len(state_ids), "duplicate enumerated state ID")
    charge_counts = Counter(state["net_formal_charge"] for state in enumerated_states)
    hydrogen_counts = [state["hydrogen_count"] for state in enumerated_states]
    canonical = next(
        state
        for state in enumerated_states
        if all(state["variable_site_states"][site_id] == "deprotonated" for site_id in macro_ids)
        and state["variable_site_states"]["nterm_A1"] == "protonated"
        and all(
            state["variable_site_states"][site_id] == "neutral_template_tautomer"
            for site_id in his_ids
        )
    )
    require(canonical["net_formal_charge"] == -7, "canonical candidate net charge is not −7")
    require(canonical["formula"] == EXPECTED_CANONICAL_FORMULA, "canonical candidate formula mismatch")
    require(min(charge_counts) == -8 and max(charge_counts) == -2, "enumerated net charge range changed")
    require(min(hydrogen_counts) == 453 and max(hydrogen_counts) == 459, "enumerated hydrogen range changed")

    result = {
        "schema": "dia-health-one.nnc2215-physiological-charge-states.v1",
        "run": {
            "audit_script": "research/audit_nnc2215_physiological_charge_states.py",
            "input_graph": "research/data/nnc2215-full-conjugate-graph.json",
            "input_graph_sha256": sha256(INPUT_PATH),
            "pH_reference": 7.4,
            "method": "Group-level formal-charge and proton-count bookkeeping from the neutral graph; no coordinates or pKa model used.",
        },
        "representation_boundaries": {
            "scope": "Conditional pH 7.4 candidate charge-state ledger for the source-derived full desB30 insulin–NNC2215 neutral constitutional graph.",
            "not_a_prediction": "No pKa values, Boltzmann populations, or experimentally measured microstate are assigned.",
            "not_force_field_ready": "Formal charges are summed at functional-group level. No atom-localized partial charges, atom types, Lennard-Jones parameters, counterions, or solvent are provided.",
            "coordinate_usage": "No coordinate files were read or edited; coordinates do not select any atom map or site state.",
            "oxygen_labeling": "Carboxylate charge belongs to the functional group; no carbonyl versus anionic oxygen label is selected.",
        },
        "neutral_input": {
            "schema": report["schema"],
            "formula": neutral_formula,
            "total_formal_charge": report["formula"].get("total_formal_charge"),
            "hydrogen_count": neutral_h,
        },
        "fixed_conventional_assumptions": {
            "description": "These conventional ionization assignments are fixed for this bookkeeping envelope; local protein and macrocycle pKa shifts are not modeled.",
            "groups": fixed_groups,
            "group_formal_charge_delta": fixed_charge_delta,
            "hydrogen_delta_from_neutral_graph": fixed_hydrogen_delta,
        },
        "enumerated_variable_groups": variable_groups,
        "excluded_sites": {
            "tyrosine_phenols": excluded_sites,
            "disulfide_cysteines": "All six cysteines are in the three declared disulfides; no free thiol state is included.",
            "modified_conjugate_endpoints": endpoint_group_notes,
            "other_groups": "All remaining neutral amides, alcohols, ethers, triazole, and aromatic groups are held in their neutral graph states for this envelope.",
        },
        "canonical_candidate": {
            "state_id": canonical["state_id"],
            "pH_reference": 7.4,
            "site_states": canonical["variable_site_states"],
            "fixed_group_charge_delta": fixed_charge_delta,
            "net_formal_charge": canonical["net_formal_charge"],
            "formula": canonical["formula"],
            "hydrogen_count": canonical["hydrogen_count"],
            "hydrogen_delta_from_neutral_graph": canonical["hydrogen_delta_from_neutral_graph"],
            "interpretation": "Convenient conventional starting candidate under declared assumptions; it is neither a measured state nor ranked by a pKa or population calculation.",
        },
        "enumeration": {
            "state_count": len(enumerated_states),
            "factorization": {
                "macrocycle_carboxyl_groups": {"count": 3, "states_each": 2},
                "A_chain_N_terminus": {"count": 1, "states": ["protonated", "neutral"]},
                "histidines": {"count": 2, "states_each": ["neutral_template_tautomer", "neutral_alternate_tautomer", "protonated"]},
            },
            "net_formal_charge_counts": {str(key): charge_counts[key] for key in sorted(charge_counts)},
            "net_formal_charge_range": [min(charge_counts), max(charge_counts)],
            "hydrogen_count_range": [min(hydrogen_counts), max(hydrogen_counts)],
            "probabilities_assigned": False,
            "states": enumerated_states,
        },
        "checks": {
            "input_graph_hash_recorded": True,
            "neutral_input_formula_matches_expected": neutral_formula == EXPECTED_NEUTRAL_FORMULA,
            "neutral_atom_element_counts_match_formula_ledger": True,
            "input_chain_sequences_match_expected": True,
            "all_neutral_input_atom_formal_charges_are_zero": True,
            "carboxyl_groups_have_neutral_input_charge_hydrogens_and_topology": True,
            "conjugate_join_records_and_graph_edges_match_declared_sites": True,
            "declared_disulfides_present_with_neutral_cysteine_sulfur_states": True,
            "four_Glu_sidechains_found": len(glutamates) == 4,
            "two_chain_C_termini_found": len(termini) == 2,
            "three_mapped_macrocycle_acids_found": len(MACROCYCLE_GROUP_MAPS) == 3,
            "A1_amine_and_ArgB22_found": True,
            "HisB5_and_HisB10_template_tautomers_found": set(histidine_states) == {"his_B5", "his_B10"},
            "histidine_site_hydrogen_counts_include_all_five_ring_atoms": all(
                group["neutral_graph_site_hydrogen_count"] == 3
                for group in variable_groups if group["site_type"] == "variable_histidine_group"
            ),
            "modified_PheB1_and_LysB29_nitrogens_remain_neutral": all(
                item["graph_formal_charge"] == 0 for item in endpoint_group_notes
            ),
            "canonical_candidate_is_charge_minus_7_and_formula_C325H454N79O97S6": canonical["net_formal_charge"] == -7 and canonical["formula"] == EXPECTED_CANONICAL_FORMULA,
            "144_group_state_combinations_enumerated_without_probabilities": len(enumerated_states) == 144 and all(not state["probability_assigned"] for state in enumerated_states),
            "charge_and_hydrogen_ranges_match_assumptions": min(charge_counts) == -8 and max(charge_counts) == -2 and min(hydrogen_counts) == 453 and max(hydrogen_counts) == 459,
            "no_atom_localized_partial_charges_or_force_field_claimed": True,
        },
        "open_gates": [
            "local pKa shifts and measured microstate populations",
            "atom-localized force-field charges and atom types for the NNC2215 conjugate",
            "histidine tautomer energetics and site-specific protonation",
            "tyrosine phenol sensitivity",
            "PheB1 hydrogen ownership and unique B29 coordinate map",
            "roof stereochemistry, solvent, counterions, and MD readiness",
        ],
    }
    OUTPUT_PATH.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote research/data/nnc2215-physiological-charge-states.json")
    print(f"Neutral graph: {neutral_formula}, net formal charge {report['formula'].get('total_formal_charge')}")
    print(f"Canonical pH 7.4 candidate: {canonical['formula']}, net formal charge {canonical['net_formal_charge']}")
    print(f"Enumerated states: {len(enumerated_states)}; charge range {min(charge_counts)} to {max(charge_counts)}")
    print(f"Hydrogen-count range: {min(hydrogen_counts)} to {max(hydrogen_counts)}")
    print("pKa/population model: none; atom-localized force-field charges: none")
    print("Coordinates read or modified: no")


if __name__ == "__main__":
    main()
