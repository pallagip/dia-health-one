#!/usr/bin/env python3
"""Audit NNC2215 B1/B29 attachment boundaries in corrected candidates.

This read-only audit joins source-map edges to exact PDB atom identities and
PDB LINK/CONECT records. Distances are reported only after those identities
and adjacency records are resolved; they never choose a mapping or assign a
hydrogen. The only file written is the JSON report.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import platform
from collections import Counter, deque
from typing import Any

from rdkit import Chem, rdBase


REPO_ROOT = Path(__file__).resolve().parents[1]
RESEARCH = REPO_ROOT / "research"
DATA = RESEARCH / "data"
STRUCTURES = RESEARCH / "structures"

PATHS = {
    "coordinate_reconciliation": DATA / "nnc2215-carbon-reconciliation-candidate.json",
    "b1_source_graph": STRUCTURES / "nnc2215-b1-glucoside-attachment.smi",
    "b1_map_legend": DATA / "nnc2215-b1-glucoside-atom-map.json",
    "b1_topology_audit": DATA / "nnc2215-b1-topology-audit.json",
    "b29_source_graph": STRUCTURES / "nnc2215-b29-linker-roof-product.smi",
    "b29_map_legend": DATA / "nnc2215-b29-linker-roof-product-map.json",
    "b29_topology_audit": DATA / "nnc2215-b29-topology-audit.json",
    "roof_source_graph": STRUCTURES / "nnc2215-roof-propyl-azide.smi",
    "roof_map_legend": DATA / "nnc2215-roof-macrocycle-atom-map.json",
    "active_open_candidate": STRUCTURES / "nnc2215-open-source-corrected-candidate.pdb",
    "inactive_closed_candidate": STRUCTURES / "nnc2215-closed-source-corrected-candidate.pdb",
}
OUTPUT_PATH = DATA / "nnc2215-attachment-boundary-audit.json"

STATE_ORDER = ("active_open", "inactive_closed")
PDB_PATH_KEY = {
    "active_open": "active_open_candidate",
    "inactive_closed": "inactive_closed_candidate",
}


@dataclass(frozen=True)
class Atom:
    model_number: int
    record_type: str
    serial: int
    name: str
    altloc: str
    residue: str
    chain: str
    residue_number: str
    insertion_code: str
    element: str
    xyz: tuple[float, float, float]

    def semantic_key(self) -> tuple[Any, ...]:
        return (
            self.model_number,
            self.record_type,
            self.chain,
            self.residue_number,
            self.insertion_code,
            self.residue,
            self.altloc,
            self.name,
            self.element,
        )

    def identity(self, state: str) -> dict[str, Any]:
        return {
            "state": state,
            "model_number": self.model_number,
            "record_type": self.record_type,
            "chain": self.chain,
            "residue_number": self.residue_number,
            "insertion_code": self.insertion_code,
            "residue_name": self.residue,
            "atom_name": self.name,
            "altloc": self.altloc,
            "element": self.element,
            "serial": self.serial,
            "coordinates_angstrom": [round(value, 3) for value in self.xyz],
        }


@dataclass(frozen=True)
class AtomSelector:
    record_type: str
    chain: str
    residue_number: str
    insertion_code: str
    residue: str
    atom_name: str
    element: str
    altloc: str = ""
    model_number: int = 1

    def matches(self, atom: Atom) -> bool:
        return (
            atom.model_number == self.model_number
            and atom.record_type == self.record_type
            and atom.chain == self.chain
            and atom.residue_number == self.residue_number
            and atom.insertion_code == self.insertion_code
            and atom.residue == self.residue
            and atom.name == self.atom_name
            and atom.element == self.element
            and atom.altloc == self.altloc
        )

    def link_matches(self, endpoint: dict[str, str], atom: Atom) -> bool:
        return (
            atom.model_number == self.model_number
            and atom.name == endpoint["atom_name"]
            and atom.altloc == endpoint["altloc"]
            and atom.residue == endpoint["residue"]
            and atom.chain == endpoint["chain"]
            and atom.residue_number == endpoint["residue_number"]
            and atom.insertion_code == endpoint["insertion_code"]
        )


@dataclass(frozen=True)
class LinkRecord:
    first_serial: int
    second_serial: int
    line_number: int
    raw_line: str


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT).as_posix()


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object in {relative(path)}")
    return value


def mapped_graph(path: Path) -> tuple[dict[int, str], dict[tuple[int, int], str]]:
    tokens = path.read_text(encoding="utf-8").split()
    if not tokens:
        raise ValueError(f"Empty mapped SMILES: {relative(path)}")
    molecule = Chem.MolFromSmiles(tokens[0])
    if molecule is None:
        raise ValueError(f"RDKit could not parse mapped SMILES: {relative(path)}")
    Chem.SanitizeMol(molecule)
    elements: dict[int, str] = {}
    for atom in molecule.GetAtoms():
        map_id = atom.GetAtomMapNum()
        if map_id <= 0 or map_id in elements:
            raise ValueError(f"Non-unique/invalid atom map ID {map_id} in {relative(path)}")
        elements[map_id] = atom.GetSymbol()
    edges: dict[tuple[int, int], str] = {}
    for bond in molecule.GetBonds():
        left = molecule.GetAtomWithIdx(bond.GetBeginAtomIdx()).GetAtomMapNum()
        right = molecule.GetAtomWithIdx(bond.GetEndAtomIdx()).GetAtomMapNum()
        pair = tuple(sorted((left, right)))
        if pair in edges:
            raise ValueError(f"Duplicate mapped source edge {pair} in {relative(path)}")
        edges[pair] = "aromatic" if bond.GetIsAromatic() else bond.GetBondType().name.lower()
    return elements, edges


def edge_record(edges: dict[tuple[int, int], str], left: int, right: int, label: str) -> dict[str, Any]:
    pair = tuple(sorted((left, right)))
    if pair not in edges:
        raise ValueError(f"Required source graph edge {left}-{right} is absent")
    if label and label not in edges[pair]:
        raise ValueError(f"Source edge {left}-{right} conflicts with its map legend")
    return {
        "map_ids": [left, right],
        "source_bond_class": edges[pair],
        "source_graph_adjacency": True,
        "provenance": "mapped source SMILES plus map legend",
    }


def graph_summary(elements: dict[int, str], edges: dict[tuple[int, int], str], target_elements: dict[int, str] | None = None) -> dict[str, Any]:
    inventory = target_elements if target_elements is not None else elements
    adjacency: dict[int, set[int]] = {map_id: set() for map_id in elements}
    for left, right in edges:
        adjacency[left].add(right)
        adjacency[right].add(left)
    remaining = set(adjacency)
    component_count = 0
    while remaining:
        component_count += 1
        seed = remaining.pop()
        queue = deque([seed])
        while queue:
            current = queue.popleft()
            for neighbor in adjacency[current] & remaining:
                remaining.remove(neighbor)
                queue.append(neighbor)
    return {
        "node_count": len(elements),
        "element_counts": dict(sorted(Counter(inventory.values()).items())),
        "edge_count": len(edges),
        "connected_component_count": component_count,
        "cycle_rank_E_minus_V_plus_components": len(edges) - len(elements) + component_count,
    }


def parse_link_endpoint(line: str, start: int) -> dict[str, str]:
    return {
        "atom_name": line[start : start + 4].strip(),
        "altloc": line[start + 4 : start + 5].strip(),
        "residue": line[start + 5 : start + 8].strip(),
        "chain": line[start + 9 : start + 10].strip(),
        "residue_number": line[start + 10 : start + 14].strip(),
        "insertion_code": line[start + 14 : start + 15].strip(),
    }


def parse_pdb(path: Path) -> tuple[list[Atom], list[LinkRecord], dict[int, set[int]], int]:
    lines = path.read_text(encoding="utf-8").splitlines()
    atoms: list[Atom] = []
    links_unresolved: list[tuple[dict[str, str], dict[str, str], int, str, int]] = []
    conect_pairs: set[tuple[int, int]] = set()
    model_number = 1
    saw_model = False
    for line_number, line in enumerate(lines, start=1):
        record = line[:6].strip()
        if record == "MODEL":
            saw_model = True
            raw_model = line[10:14].strip()
            model_number = int(raw_model) if raw_model else model_number + 1
        elif record == "ENDMDL":
            model_number += 1
        elif record in {"ATOM", "HETATM"}:
            if len(line) < 78:
                raise ValueError(f"Short PDB atom record at {relative(path)}:{line_number}")
            atom = Atom(
                model_number=model_number,
                record_type=record,
                serial=int(line[6:11]),
                name=line[12:16].strip(),
                altloc=line[16:17].strip(),
                residue=line[17:20].strip(),
                chain=line[21:22].strip(),
                residue_number=line[22:26].strip(),
                insertion_code=line[26:27].strip(),
                xyz=(float(line[30:38]), float(line[38:46]), float(line[46:54])),
                element=line[76:78].strip().upper(),
            )
            if not atom.element:
                raise ValueError(f"Missing element in PDB atom record at {relative(path)}:{line_number}")
            atoms.append(atom)
        elif record == "LINK":
            links_unresolved.append(
                (parse_link_endpoint(line, 12), parse_link_endpoint(line, 42), line_number, line, model_number)
            )
        elif record == "CONECT":
            if len(line) < 11 or not line[6:11].strip():
                raise ValueError(f"Malformed CONECT record at {relative(path)}:{line_number}")
            source = int(line[6:11])
            targets = [
                int(line[offset : offset + 5])
                for offset in range(11, len(line), 5)
                if line[offset : offset + 5].strip()
            ]
            for target in targets:
                if source != target:
                    conect_pairs.add(tuple(sorted((source, target))))

    if not atoms:
        raise ValueError(f"No ATOM/HETATM records in {relative(path)}")
    serials: dict[int, Atom] = {}
    for atom in atoms:
        if atom.serial in serials:
            raise ValueError(f"Duplicate PDB serial {atom.serial} in {relative(path)}")
        serials[atom.serial] = atom
    for pair in conect_pairs:
        if pair[0] not in serials or pair[1] not in serials:
            raise ValueError(f"CONECT references absent PDB atom serials {pair} in {relative(path)}")

    links: list[LinkRecord] = []
    for first, second, line_number, raw_line, link_model_number in links_unresolved:
        first_matches = [
            atom for atom in atoms
            if atom.model_number == link_model_number
            and atom.name == first["atom_name"]
            and atom.altloc == first["altloc"]
            and atom.residue == first["residue"]
            and atom.chain == first["chain"]
            and atom.residue_number == first["residue_number"]
            and atom.insertion_code == first["insertion_code"]
        ]
        second_matches = [
            atom for atom in atoms
            if atom.model_number == link_model_number
            and atom.name == second["atom_name"]
            and atom.altloc == second["altloc"]
            and atom.residue == second["residue"]
            and atom.chain == second["chain"]
            and atom.residue_number == second["residue_number"]
            and atom.insertion_code == second["insertion_code"]
        ]
        if len(first_matches) != 1 or len(second_matches) != 1:
            raise ValueError(
                f"LINK endpoints must resolve uniquely at {relative(path)}:{line_number}; "
                f"found {len(first_matches)} and {len(second_matches)}"
            )
        links.append(LinkRecord(first_matches[0].serial, second_matches[0].serial, line_number, raw_line))

    # PDB files without MODEL records have one implicit coordinate model.
    if not saw_model and {atom.model_number for atom in atoms} != {1}:
        raise ValueError(f"Unexpected model numbering in {relative(path)}")
    conect_adjacency: dict[int, set[int]] = {serial: set() for serial in serials}
    for left, right in conect_pairs:
        conect_adjacency[left].add(right)
        conect_adjacency[right].add(left)
    return atoms, links, conect_adjacency, len(lines)


def require_one_atom(atoms: list[Atom], selector: AtomSelector, description: str) -> Atom:
    matches = [atom for atom in atoms if selector.matches(atom)]
    if len(matches) != 1:
        raise ValueError(f"{description} must resolve to one full PDB identity; found {len(matches)}")
    return matches[0]


def pair_distance(left: Atom, right: Atom) -> float:
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(left.xyz, right.xyz)))


def atom_sort_key(atom: Atom) -> tuple[Any, ...]:
    return atom.model_number, atom.chain, atom.residue_number, atom.insertion_code, atom.residue, atom.name, atom.serial


def neighbor_ledger(
    state: str,
    atom: Atom,
    atoms_by_serial: dict[int, Atom],
    links: list[LinkRecord],
    conect_adjacency: dict[int, set[int]],
) -> dict[str, Any]:
    conect_atoms = [atoms_by_serial[serial] for serial in conect_adjacency.get(atom.serial, set())]
    link_serials: set[int] = set()
    for link in links:
        if link.first_serial == atom.serial:
            link_serials.add(link.second_serial)
        elif link.second_serial == atom.serial:
            link_serials.add(link.first_serial)
    link_atoms = [atoms_by_serial[serial] for serial in link_serials]
    combined: dict[int, Atom] = {neighbor.serial: neighbor for neighbor in [*conect_atoms, *link_atoms]}
    heavy = sorted((neighbor for neighbor in combined.values() if neighbor.element != "H"), key=atom_sort_key)
    hydrogens = sorted((neighbor for neighbor in conect_atoms if neighbor.element == "H"), key=atom_sort_key)
    return {
        "center_atom": atom.identity(state),
        "conect_neighbors": {
            "heavy_atom_identities": [n.identity(state) for n in sorted((x for x in conect_atoms if x.element != "H"), key=atom_sort_key)],
            "hydrogen_atom_identities": [n.identity(state) for n in hydrogens],
        },
        "link_neighbors": [n.identity(state) for n in sorted(link_atoms, key=atom_sort_key)],
        "adjacent_heavy_atom_identities_from_conect_or_link": [n.identity(state) for n in heavy],
        "adjacency_sources_by_neighbor_serial": {
            str(neighbor.serial): [
                source for source, exists in (
                    ("CONECT", neighbor.serial in conect_adjacency.get(atom.serial, set())),
                    ("LINK", any(
                        {link.first_serial, link.second_serial} == {atom.serial, neighbor.serial}
                        for link in links
                    )),
                ) if exists
            ]
            for neighbor in sorted(combined.values(), key=atom_sort_key)
        },
    }


def observed_pair_sources(
    left: Atom,
    right: Atom,
    links: list[LinkRecord],
    conect_adjacency: dict[int, set[int]],
) -> list[str]:
    sources: list[str] = []
    if right.serial in conect_adjacency.get(left.serial, set()):
        sources.append("CONECT")
    if any({item.first_serial, item.second_serial} == {left.serial, right.serial} for item in links):
        sources.append("LINK")
    return sources


def unique_link_records(left: Atom, right: Atom, links: list[LinkRecord]) -> list[LinkRecord]:
    return [item for item in links if {item.first_serial, item.second_serial} == {left.serial, right.serial}]


def serialize_join(
    state: str,
    join_id: str,
    join_type: str,
    left: Atom,
    right: Atom,
    adjacency_sources: list[str],
    source_edge: dict[str, Any] | None = None,
    link_records: list[LinkRecord] | None = None,
) -> dict[str, Any]:
    record: dict[str, Any] = {
        "join_id": join_id,
        "join_type": join_type,
        "state": state,
        "atom_identities": [left.identity(state), right.identity(state)],
        "adjacency_record_sources": adjacency_sources,
        "distance_angstrom": round(pair_distance(left, right), 6),
        "distance_interpretation": "coordinate observation only; not used for identity, mapping, bond-order, or hydrogen assignment",
    }
    if source_edge is not None:
        record["source_graph_edge"] = source_edge
    if link_records is not None:
        record["link_record_count_for_resolved_pair"] = len(link_records)
        record["link_records"] = [
            {"line_number": item.line_number, "record": item.raw_line}
            for item in link_records
        ]
    return record


def recorded_input(
    path: Path,
    observed_hash: str,
    expected_hash: str | None = None,
    expected_source: str | None = None,
) -> dict[str, Any]:
    return {
        "path": relative(path),
        "sha256_observed": observed_hash,
        "sha256_expected": expected_hash,
        "expected_hash_source": expected_source,
        "hash_check": (
            "pass" if expected_hash is not None and observed_hash == expected_hash
            else "recorded_only_no_independent_expected_hash" if expected_hash is None
            else "fail"
        ),
    }


def main() -> None:
    initial_hashes = {key: sha256(path) for key, path in PATHS.items()}
    script_hash = sha256(Path(__file__))
    reconciliation = load_json(PATHS["coordinate_reconciliation"])
    b1_legend = load_json(PATHS["b1_map_legend"])
    b1_audit = load_json(PATHS["b1_topology_audit"])
    b29_legend = load_json(PATHS["b29_map_legend"])
    b29_audit = load_json(PATHS["b29_topology_audit"])
    roof_legend = load_json(PATHS["roof_map_legend"])

    b1_elements, b1_edges = mapped_graph(PATHS["b1_source_graph"])
    b29_elements, b29_edges = mapped_graph(PATHS["b29_source_graph"])
    roof_elements, roof_edges = mapped_graph(PATHS["roof_source_graph"])

    b1_expected_keys = {str(i) for i in range(1, 20)}
    b29_expected_keys = {str(i) for i in range(1, 93)}
    roof_expected_keys = {str(i) for i in range(1, 84)}
    if set(b1_elements) != set(range(1, 20)) or set(b1_legend.get("map_id_legend", {})) != b1_expected_keys:
        raise ValueError("B1 source graph/map legend must cover map IDs 1 through 19")
    if set(b29_elements) != set(range(1, 93)):
        raise ValueError("B29 source graph must cover map IDs 1 through 92")
    if set(roof_elements) != set(range(1, 84)) or set(roof_legend.get("map_id_legend", {}).get("roof_tail_source_C21_C22_C23", {})) != {"C21", "C22", "C23"}:
        raise ValueError("Roof source graph/map legend must identify C21-C23 at maps 78-80")
    if b1_elements[1] != "*" or b1_legend["map_id_legend"]["1"].get("model_match_element") != "N":
        raise ValueError("B1 map 1 must be the dummy placeholder for PheB1 N")
    if b1_elements[2] != "C" or "carbonyl" not in b1_legend["map_id_legend"]["2"].get("role", "").lower():
        raise ValueError("B1 map 2 must be the linker carbonyl carbon")
    if b29_legend.get("map_id_legend", {}).get("added_product_atoms_84_to_92", {}).get("84") != "source linker C3, carbamate carbonyl carbon":
        raise ValueError("B29 map 84 legend no longer identifies the source C3 carbamate carbonyl carbon")
    if b29_legend.get("map_id_legend", {}).get("added_product_atoms_84_to_92", {}).get("92") != "LysB29 NZ attachment endpoint":
        raise ValueError("B29 map 92 legend no longer identifies the LysB29 NZ endpoint")
    if b29_elements[84] != "C" or b29_elements[92] != "N":
        raise ValueError("B29 map 84/92 source elements must be carbonyl carbon and nitrogen")
    if roof_elements[78] != "C" or roof_elements[79] != "C":
        raise ValueError("Roof maps 78/79 must be carbon atoms")
    roof_tail_map = roof_legend["map_id_legend"]["roof_tail_source_C21_C22_C23"]
    if roof_tail_map != {"C21": 78, "C22": 79, "C23": 80}:
        raise ValueError("Roof map legend C21-C23 assignments differ from the approved boundary map IDs")
    b29_roof_tail_map = b29_legend["map_id_legend"]["inherited_roof_map_ids_1_to_83"]["roof_tail_source_C21_C22_C23"]
    if b29_roof_tail_map != roof_tail_map:
        raise ValueError("B29 product map roof-tail assignments conflict with the roof map legend")

    b1_legend_bonds: dict[tuple[int, int], str] = {}
    for item in b1_legend.get("key_connectivity", []):
        pair = tuple(sorted(int(value) for value in item["atoms"]))
        if pair in b1_legend_bonds:
            raise ValueError(f"Duplicate edge in B1 map legend: {pair}")
        b1_legend_bonds[pair] = str(item["bond"]).lower()
    b1_graph_bonds = {pair: bond for pair, bond in b1_edges.items()}
    if b1_graph_bonds != b1_legend_bonds:
        raise ValueError("B1 mapped SMILES edge set differs from its complete map-legend connectivity")

    b1_edge = edge_record(b1_edges, 1, 2, b1_legend_bonds.get((1, 2), ""))
    b29_edge = edge_record(b29_edges, 84, 92, "")
    roof_edge = edge_record(roof_edges, 78, 79, "")
    b29_roof_edge = edge_record(b29_edges, 78, 79, "")
    if roof_edge["source_bond_class"] != b29_roof_edge["source_bond_class"]:
        raise ValueError("Roof source edge 78-79 has inconsistent bond classes in the roof and B29 product graphs")

    # Compare source graph summaries and the exact-match counts already
    # reported by the component audits before inspecting candidate geometry.
    b1_target_elements = dict(b1_elements)
    b1_target_elements[1] = b1_legend["map_id_legend"]["1"]["model_match_element"]
    b1_expected_summary = graph_summary(b1_elements, b1_edges, b1_target_elements)
    b1_states_by_name = {state["state"]: state for state in b1_audit["candidate_states"]}
    b29_states_by_name = {state["state"]: state for state in b29_audit["states"]}
    if set(b1_states_by_name) != set(STATE_ORDER) or set(b29_states_by_name) != set(STATE_ORDER):
        raise ValueError("Existing topology audits do not contain both expected candidate states")
    if b1_audit.get("overall_correspondence_status") != "both_candidates_exact":
        raise ValueError("Existing B1 topology audit does not report exact correspondence in both states")
    b1_mapping_cross_checks: dict[str, dict[str, Any]] = {}
    b29_mapping_cross_checks: dict[str, dict[str, Any]] = {}
    for state in STATE_ORDER:
        b1_state = b1_states_by_name[state]
        b1_source_summary = b1_state.get("source_graph", {})
        b1_target_summary = b1_state.get("target_graph_summary", {})
        if b1_state.get("exact_mapping_count") != 1:
            raise ValueError(f"Existing B1 audit mapping count is not exactly one in {state}")
        if b1_source_summary != b1_expected_summary or b1_target_summary != b1_expected_summary:
            raise ValueError(f"B1 topology-audit source/target graph summary changed in {state}")
        exact_b1 = b1_state["all_exact_mappings"]
        if len(exact_b1) != 1:
            raise ValueError(f"B1 exact mapping ledger has an unexpected size in {state}")
        b1_map_to_target = exact_b1[0]["map_id_to_target"]
        b1_mapping_cross_checks[state] = {
            "exact_mapping_count": 1,
            "source_graph_summary": b1_source_summary,
            "map_1_target_atom": b1_map_to_target["1"]["target_atom"],
            "map_2_target_atom": b1_map_to_target["2"]["target_atom"],
        }

        b29_state = b29_states_by_name[state]
        comparison = b29_state["complete_product_exact_graph_comparison"]
        if comparison.get("exact_element_labeled_graph_isomorphism_found") is not True:
            raise ValueError(f"Existing B29 topology audit lacks exact graph correspondence in {state}")
        if comparison.get("exact_mapping_count") != 16 or comparison.get("matches_independent_feasibility_review_count") is not True:
            raise ValueError(f"Existing B29 exact mapping count is not 16 in {state}")
        if comparison.get("source_node_count") != len(b29_elements) or comparison.get("source_edge_count") != len(b29_edges):
            raise ValueError(f"B29 graph summary in the topology audit differs from the source mapped SMILES in {state}")
        mappings = comparison.get("exact_mappings", [])
        if len(mappings) != 16:
            raise ValueError(f"B29 exact mapping ledger does not retain all 16 mappings in {state}")
        anchored_serials: dict[str, set[int]] = {str(map_id): set() for map_id in (78, 79, 84, 92)}
        for mapping in mappings:
            source_to_serial = mapping.get("map_id_to_pdb_serial", {})
            for map_id in anchored_serials:
                if map_id not in source_to_serial:
                    raise ValueError(f"B29 mapping {mapping.get('mapping_id')} lacks source map {map_id} in {state}")
                anchored_serials[map_id].add(int(source_to_serial[map_id]))
        if any(len(values) != 1 for values in anchored_serials.values()):
            raise ValueError(f"B29 boundary source maps 78/79/84/92 are not invariant across the 16 mappings in {state}")
        b29_mapping_cross_checks[state] = {
            "exact_mapping_count": len(mappings),
            "full_ambiguity_preserved": True,
            "source_map_serials_invariant_across_all_mappings": {
                map_id: next(iter(values)) for map_id, values in anchored_serials.items()
            },
            "source_graph_summary": {
                "node_count": comparison.get("source_node_count"),
                "edge_count": comparison.get("source_edge_count"),
                "connected_component_count": comparison.get("source_connected_component_count"),
                "cycle_rank_E_minus_V_plus_components": comparison.get("source_cycle_rank_E_minus_V_plus_components"),
            },
        }

    reconciliation_models = {item["state"]: item for item in reconciliation["models"]}
    if set(reconciliation_models) != set(STATE_ORDER):
        raise ValueError("Coordinate reconciliation record does not contain both candidate states")

    b1_hashes = b1_audit["input_hashes_sha256"]
    expected_hashes: dict[str, tuple[str, str]] = {
        "coordinate_reconciliation": (b1_hashes["coordinate_reconciliation"]["sha256"], "B1 topology-audit input hash"),
        "b1_source_graph": (b1_hashes["source_graph"]["sha256"], "B1 topology-audit input hash"),
        "b1_map_legend": (b1_hashes["map_legend"]["sha256"], "B1 topology-audit input hash"),
        "b29_source_graph": (b29_legend["graph_sha256"], "B29 product map legend graph_sha256"),
        "b29_map_legend": (b29_audit["source_product_graph"]["map_legend_sha256"], "B29 topology-audit source_product_graph hash"),
        "roof_source_graph": (b29_legend["input_audit_files"]["roof_graph"]["sha256"], "B29 product map legend input hash"),
        "roof_map_legend": (b29_legend["input_audit_files"]["roof_atom_map_legend"]["sha256"], "B29 product map legend input hash"),
        "active_open_candidate": (reconciliation_models["active_open"]["output_sha256"], "coordinate reconciliation candidate output hash"),
        "inactive_closed_candidate": (reconciliation_models["inactive_closed"]["output_sha256"], "coordinate reconciliation candidate output hash"),
    }
    # Use the cross-state reconciliation file's own hashes as the direct
    # candidate expectations; the input record above is checked separately.
    input_hash_records: dict[str, dict[str, Any]] = {}
    for key, path in PATHS.items():
        expected_hash: str | None = None
        expected_source: str | None = None
        if key in expected_hashes:
            expected_hash, expected_source = expected_hashes[key]
        input_hash_records[key] = recorded_input(path, initial_hashes[key], expected_hash, expected_source)
        if input_hash_records[key]["hash_check"] == "fail":
            raise ValueError(f"Input SHA-256 mismatch for {relative(path)}")
    input_hash_records["audit_script"] = recorded_input(Path(__file__), script_hash)

    # Validate both corrected PDB hashes against the reconciliation record and
    # existing topology audits before parsing any coordinates.
    candidate_integrity: dict[str, dict[str, Any]] = {}
    for state in STATE_ORDER:
        path = PATHS[PDB_PATH_KEY[state]]
        observed = initial_hashes[PDB_PATH_KEY[state]]
        expected = reconciliation_models[state]["output_sha256"]
        if observed != expected:
            raise ValueError(f"Candidate PDB SHA-256 differs from reconciliation record for {state}")
        b1_state = b1_states_by_name[state]
        if b1_state.get("candidate_pdb_sha256") != observed or b1_state.get("candidate_pdb_hash_matches_reconciliation") is not True:
            raise ValueError(f"B1 topology audit candidate hash cross-check failed for {state}")
        b29_state = b29_states_by_name[state]
        if b29_state.get("coordinate_candidate_sha256_before_audit") != observed or b29_state.get("coordinate_candidate_sha256_after_audit") != observed:
            raise ValueError(f"B29 topology audit candidate hash cross-check failed for {state}")
        candidate_integrity[state] = {
            "path": relative(path),
            "sha256_before_parsing": observed,
            "sha256_expected_by_reconciliation": expected,
            "sha256_expected_by_existing_topology_audits": {
                "b1_topology_audit": b1_state.get("candidate_pdb_sha256"),
                "b29_topology_audit_before": b29_state.get("coordinate_candidate_sha256_before_audit"),
                "b29_topology_audit_after": b29_state.get("coordinate_candidate_sha256_after_audit"),
            },
        }

    state_reports: list[dict[str, Any]] = []
    attachment_link_total = 0
    native_context_total = 0
    atom_selectors_by_state: dict[str, dict[str, AtomSelector]] = {}
    for state in STATE_ORDER:
        roof_name = "CD6" if state == "active_open" else "C"
        atom_selectors_by_state[state] = {
            "phe_n": AtomSelector("ATOM", "B", "1", "", "PHE", "N", "N"),
            "phe_ca": AtomSelector("ATOM", "B", "1", "", "PHE", "CA", "C"),
            "b1_carbonyl": AtomSelector("HETATM", "B", "0", "", "PHX", "C2", "C"),
            "lys_nz": AtomSelector("ATOM", "B", "29", "", "LYS", "NZ", "N"),
            "lys_ce": AtomSelector("ATOM", "B", "29", "", "LYS", "CE", "C"),
            "b29_carbonyl": AtomSelector("HETATM", "B", "30", "", "LYX", "CD1", "C"),
            "roof_c21": AtomSelector("HETATM", "C", "1", "", "MC", "C31", "C"),
            "roof_c22": AtomSelector("HETATM", "B", "30", "", "LYX", roof_name, "C"),
            "closed_nearby_h2a": AtomSelector("HETATM", "B", "0", "", "PHX", "H2A", "H"),
        }
        path = PATHS[PDB_PATH_KEY[state]]
        atoms, links, conect_adjacency, _line_count = parse_pdb(path)
        atoms_by_serial = {atom.serial: atom for atom in atoms}
        selector = atom_selectors_by_state[state]
        resolved = {
            key: require_one_atom(atoms, value, f"{state} {key}")
            for key, value in selector.items()
            if key != "closed_nearby_h2a"
        }
        # Candidate coordinates may contain H2A in either state; this is only
        # a distance observation and never creates a connectivity edge.
        h2a_matches = [atom for atom in atoms if selector["closed_nearby_h2a"].matches(atom)]
        if len(h2a_matches) > 1:
            raise ValueError(f"{state} has duplicate PHX H2A identities")
        h2a = h2a_matches[0] if h2a_matches else None

        # Existing topology-audit correspondences are checked against the
        # uniquely resolved atoms; they never select an atom by geometry.
        b1_map1 = b1_mapping_cross_checks[state]["map_1_target_atom"]
        b1_map2 = b1_mapping_cross_checks[state]["map_2_target_atom"]
        for mapped, observed, description in (
            (b1_map1, resolved["phe_n"], "B1 map 1 / PheB1 N"),
            (b1_map2, resolved["b1_carbonyl"], "B1 map 2 / PHX C2"),
        ):
            expected_identity = {
                "serial": observed.serial,
                "record_type": observed.record_type,
                "chain": observed.chain,
                "residue_number": observed.residue_number,
                "insertion_code": observed.insertion_code,
                "residue": observed.residue,
                "altloc": observed.altloc,
                "name": observed.name,
                "element": observed.element,
            }
            if any(mapped.get(field) != value for field, value in expected_identity.items()):
                raise ValueError(f"Existing B1 audit {description} full identity differs in {state}")
        for source_map, key in (("78", "roof_c21"), ("79", "roof_c22"), ("84", "b29_carbonyl"), ("92", "lys_nz")):
            mapped_serial = b29_mapping_cross_checks[state]["source_map_serials_invariant_across_all_mappings"][source_map]
            if mapped_serial != resolved[key].serial:
                raise ValueError(f"B29 source map {source_map} does not agree with the resolved {key} identity in {state}")

        attachment_specs = [
            ("b1_glucoside_attachment", "b1_component_boundary", resolved["phe_n"], resolved["b1_carbonyl"], b1_edge),
            ("b29_lys_attachment", "b29_component_boundary", resolved["b29_carbonyl"], resolved["lys_nz"], b29_edge),
            ("roof_tail_attachment", "roof_tail_inter_residue_boundary", resolved["roof_c21"], resolved["roof_c22"], roof_edge),
        ]
        attachment_joins: list[dict[str, Any]] = []
        for join_id, join_type, left, right, source_edge in attachment_specs:
            matching_links = unique_link_records(left, right, links)
            if len(matching_links) != 1:
                raise ValueError(f"{state} {join_id} must have exactly one matching LINK record; found {len(matching_links)}")
            sources = observed_pair_sources(left, right, links, conect_adjacency)
            if "LINK" not in sources:
                raise ValueError(f"{state} {join_id} lacks a LINK adjacency record")
            attachment_joins.append(serialize_join(
                state, join_id, join_type, left, right, sources, source_edge, matching_links
            ))
        attachment_link_total += len(attachment_joins)

        native_specs = [
            ("pheb1_n_ca_context", resolved["phe_n"], resolved["phe_ca"]),
            ("lysb29_ce_nz_context", resolved["lys_ce"], resolved["lys_nz"]),
        ]
        native_joins: list[dict[str, Any]] = []
        for join_id, left, right in native_specs:
            sources = observed_pair_sources(left, right, links, conect_adjacency)
            if not sources:
                raise ValueError(f"{state} native context pair {join_id} has no CONECT/LINK adjacency")
            native_joins.append(serialize_join(state, join_id, "native_residue_context", left, right, sources))
        native_context_total += len(native_joins)

        phe_neighbors = neighbor_ledger(state, resolved["phe_n"], atoms_by_serial, links, conect_adjacency)
        lys_neighbors = neighbor_ledger(state, resolved["lys_nz"], atoms_by_serial, links, conect_adjacency)
        phe_heavy = {
            item["serial"] for item in phe_neighbors["adjacent_heavy_atom_identities_from_conect_or_link"]
        }
        lys_heavy = {
            item["serial"] for item in lys_neighbors["adjacent_heavy_atom_identities_from_conect_or_link"]
        }
        phe_expected_heavy = {resolved["phe_ca"].serial, resolved["b1_carbonyl"].serial}
        lys_expected_heavy = {resolved["lys_ce"].serial, resolved["b29_carbonyl"].serial}
        if phe_heavy != phe_expected_heavy:
            raise ValueError(f"PheB1 N heavy-neighbor pattern differs from the expected local boundary in {state}")
        if lys_heavy != lys_expected_heavy:
            raise ValueError(f"LysB29 NZ heavy-neighbor pattern differs from the expected local boundary in {state}")

        phe_hydrogens = phe_neighbors["conect_neighbors"]["hydrogen_atom_identities"]
        lys_hydrogens = lys_neighbors["conect_neighbors"]["hydrogen_atom_identities"]
        if state == "active_open":
            if len(phe_hydrogens) != 1 or phe_hydrogens[0]["atom_name"] != "H2":
                raise ValueError("Open PheB1 N must retain its observed single CONECT hydrogen neighbor H2")
        elif phe_hydrogens:
            raise ValueError("Closed PheB1 N unexpectedly has a CONECT hydrogen neighbor")
        if len(lys_hydrogens) != 1 or lys_hydrogens[0]["atom_name"] != "HZ2":
            raise ValueError(f"LysB29 NZ must have the recorded HZ2 CONECT hydrogen neighbor in {state}")

        h2a_observation = None
        if h2a is not None:
            h2a_distance = pair_distance(resolved["phe_n"], h2a)
            h2a_connected = h2a.serial in conect_adjacency.get(resolved["phe_n"].serial, set())
            h2a_observation = {
                "atom_identity": h2a.identity(state),
                "distance_from_pheB1_N_angstrom": round(h2a_distance, 6),
                "has_pheB1_N_conect_edge": h2a_connected,
                "interpretation": "distance observation only; no hydrogen ownership is assigned",
            }

        expected_neutral_patterns = {
            "pheB1_N": {
                "expected_heavy_neighbors": [resolved["phe_ca"].identity(state), resolved["b1_carbonyl"].identity(state)],
                "expected_conect_hydrogen_neighbor_count": 1,
                "observed_conect_hydrogen_neighbors": phe_hydrogens,
                "record_pattern_status": "matches_expected" if len(phe_hydrogens) == 1 else "hydrogen_record_discrepancy",
                "chemistry_validation": "not_established by coordinate-file connectivity",
            },
            "lysB29_NZ": {
                "expected_heavy_neighbors": [resolved["lys_ce"].identity(state), resolved["b29_carbonyl"].identity(state)],
                "expected_conect_hydrogen_neighbor_count": 1,
                "observed_conect_hydrogen_neighbors": lys_hydrogens,
                "record_pattern_status": "matches_expected_record_count",
                "chemistry_validation": "not_established by coordinate-file connectivity",
            },
        }
        state_reports.append({
            "state": state,
            "candidate_pdb": relative(path),
            "mapping_cross_checks": {
                "b1_topology_audit": b1_mapping_cross_checks[state],
                "b29_topology_audit": b29_mapping_cross_checks[state],
            },
            "attachment_boundary_joins": attachment_joins,
            "native_context_joins": native_joins,
            "local_attachment_neighborhoods": {
                "pheB1_N": phe_neighbors,
                "lysB29_NZ": lys_neighbors,
            },
            "expected_neutral_local_boundary_record_patterns": expected_neutral_patterns,
            "pheB1_H2A_proximity_observation": h2a_observation,
        })

    # Preserve the precise, approved open/closed discrepancy without treating
    # the closed PHX H2A distance as a bond or changing either coordinate file.
    open_report, closed_report = state_reports
    open_h = open_report["local_attachment_neighborhoods"]["pheB1_N"]["conect_neighbors"]["hydrogen_atom_identities"]
    closed_h = closed_report["local_attachment_neighborhoods"]["pheB1_N"]["conect_neighbors"]["hydrogen_atom_identities"]
    closed_h2a = closed_report["pheB1_H2A_proximity_observation"]
    if len(open_h) != 1 or open_h[0]["atom_name"] != "H2" or closed_h:
        raise ValueError("The approved open/closed PheB1 N hydrogen-record discrepancy did not resolve as expected")
    if closed_h2a is None or closed_h2a["has_pheB1_N_conect_edge"] is not False or not math.isclose(
        closed_h2a["distance_from_pheB1_N_angstrom"], 1.006, abs_tol=0.001
    ):
        raise ValueError("Closed PHX H2A must remain a 1.006 Å proximity observation without PheB1 N CONECT")

    final_hashes = {key: sha256(path) for key, path in PATHS.items()}
    if final_hashes != initial_hashes:
        changed = [key for key in PATHS if initial_hashes[key] != final_hashes[key]]
        raise ValueError(f"One or more audit input files changed during the read-only audit: {changed}")
    for state in STATE_ORDER:
        key = PDB_PATH_KEY[state]
        candidate_integrity[state]["sha256_after_parsing"] = final_hashes[key]
        candidate_integrity[state]["unchanged_during_audit"] = final_hashes[key] == initial_hashes[key]

    output: dict[str, Any] = {
        "schema": "dia-health-one.nnc2215-attachment-boundary-audit.v1",
        "run": {
            "generated_at_utc": datetime.now(timezone.utc).isoformat(),
            "script": relative(Path(__file__)),
            "script_sha256": script_hash,
            "python_version": platform.python_version(),
            "rdkit_version": rdBase.rdkitVersion,
            "coordinate_files_modified": False,
            "hydrogen_or_coordinate_records_modified": False,
        },
        "inputs": {
            "sha256_file_integrity_records": input_hash_records,
            "candidate_pdb_before_after_hashes": candidate_integrity,
        },
        "source_graph_boundary_edges": {
            "b1_maps_1_2": {
                "graph_path": relative(PATHS["b1_source_graph"]),
                "legend_path": relative(PATHS["b1_map_legend"]),
                **b1_edge,
                "map_legend_roles": {
                    "1": b1_legend["map_id_legend"]["1"]["role"],
                    "2": b1_legend["map_id_legend"]["2"]["role"],
                },
            },
            "b29_maps_84_92": {
                "graph_path": relative(PATHS["b29_source_graph"]),
                "legend_path": relative(PATHS["b29_map_legend"]),
                **b29_edge,
                "map_legend_roles": {
                    "84": b29_legend["map_id_legend"]["added_product_atoms_84_to_92"]["84"],
                    "92": b29_legend["map_id_legend"]["added_product_atoms_84_to_92"]["92"],
                },
            },
            "roof_tail_maps_78_79": {
                "graph_path": relative(PATHS["roof_source_graph"]),
                "legend_path": relative(PATHS["roof_map_legend"]),
                **roof_edge,
                "also_present_in_b29_product_graph": b29_roof_edge["source_graph_adjacency"],
                "map_legend_roles": {"78": "source macrocycle C21", "79": "source macrocycle C22"},
            },
        },
        "states": state_reports,
        "discrepancies": [{
            "id": "pheB1_N_hydrogen_connectivity_open_closed_difference",
            "state": "inactive_closed",
            "expected_neutral_local_pattern": "PheB1 N has one hydrogen neighbor in addition to CA and the B1 carbonyl carbon",
            "observed_pdb_connectivity": "PheB1 N has no CONECT hydrogen neighbor in the closed candidate; the open candidate records H2 as a CONECT neighbor",
            "nearby_observation": closed_h2a,
            "interpretation": "PHX H2A is 1.006 Å from PheB1 N but has no CONECT edge to it. The coordinate file does not establish hydrogen ownership; no reassignment or repair was made.",
            "heavy_atom_boundary_closure_affected": False,
        }],
        "assessment": {
            "heavy_atom_boundary_closure": "pass" if attachment_link_total == 6 and native_context_total == 4 else "fail",
            "unique_attachment_LINK_count": attachment_link_total,
            "unique_native_context_join_count": native_context_total,
            "all_source_boundary_edges_cross_checked": True,
            "b1_exact_mappings_per_state_preserved": 1,
            "b29_exact_mappings_per_state_preserved": 16,
            "b29_ambiguity_resolved_by_coordinates": False,
            "hydrogen_record_chemistry_resolved": False,
            "full_conjugate_graph_bond_orders_stereochemistry_resolved": False,
            "whole_insulin_protonation_or_charge_resolved": False,
            "force_field_validated": False,
            "target_binding_md_run": False,
            "clinical_readiness_established": False,
            "interpretation": "The three named heavy-atom attachment boundaries and two native residue context joins per state resolve by source-map and PDB record evidence in corrected model-built candidates. Distances are observations only. The PheB1 hydrogen-record discrepancy and broader chemistry/model gates remain open.",
        },
        "limits": [
            "Corrected candidate PDBs are model-built derivatives, not experimentally solved NNC2215 structures.",
            "PDB LINK and CONECT records provide adjacency evidence only; they do not validate bond order, formal charge, stereochemistry, protonation, or hydrogen ownership.",
            "The 16 exact B29 source-to-candidate mappings per state are preserved; this audit does not choose a unique full B29 map or resolve carboxyl oxygen-label ambiguity.",
            "This local boundary audit does not establish the complete conjugate graph, full insulin charge/protonation, force-field parameters, target-binding MD readiness, or clinical readiness.",
        ],
    }
    OUTPUT_PATH.write_text(json.dumps(output, indent=2, sort_keys=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": relative(OUTPUT_PATH),
        "schema": output["schema"],
        "heavy_atom_boundary_closure": output["assessment"]["heavy_atom_boundary_closure"],
        "unique_attachment_LINK_count": attachment_link_total,
        "unique_native_context_join_count": native_context_total,
        "pheB1_hydrogen_discrepancy": output["discrepancies"][0]["id"],
        "closed_PHX_H2A_distance_angstrom": closed_h2a["distance_from_pheB1_N_angstrom"],
        "candidate_pdbs_unchanged": all(item["unchanged_during_audit"] for item in candidate_integrity.values()),
    }, indent=2))


if __name__ == "__main__":
    main()
