#!/usr/bin/env python3
"""Map the source-derived NNC2215 heavy-atom graph onto both model candidates.

Insulin atoms are matched by exact PDB identity. Existing B1/B29 audits supply
the conjugate mappings; all 16 B29 assignments are retained per state.
Coordinates and distances never select or rank mappings. Candidate PDBs are
read-only inputs.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
RESEARCH = ROOT / "research"
DATA = RESEARCH / "data"
STRUCTURES = RESEARCH / "structures"
OUTPUT_PATH = DATA / "nnc2215-full-coordinate-mapping.json"
PATHS = {
    "full_graph": DATA / "nnc2215-full-conjugate-graph.json",
    "reconciliation": DATA / "nnc2215-carbon-reconciliation-candidate.json",
    "b1_topology": DATA / "nnc2215-b1-topology-audit.json",
    "b29_topology": DATA / "nnc2215-b29-topology-audit.json",
    "boundary_audit": DATA / "nnc2215-attachment-boundary-audit.json",
    "charge_ledger": DATA / "nnc2215-physiological-charge-states.json",
    "hybrid_audit": DATA / "nnc2215-hybrid-feasibility.json",
    "active_open_pdb": STRUCTURES / "nnc2215-open-source-corrected-candidate.pdb",
    "inactive_closed_pdb": STRUCTURES / "nnc2215-closed-source-corrected-candidate.pdb",
}
STATE_ORDER = ("active_open", "inactive_closed")
STATE_PDB_KEY = {"active_open": "active_open_pdb", "inactive_closed": "inactive_closed_pdb"}
EXPECTED_ATOMS = 507
EXPECTED_BONDS = 527
EXPECTED_PROTEIN_ATOMS = 398
EXPECTED_B1_ATOMS = 18
EXPECTED_B29_ATOMS = 91
EXPECTED_BGL_ATOMS = 12
EXPECTED_B29_MAPPINGS = 16


class AuditError(ValueError):
    """Raised when a provenance or coordinate-mapping requirement fails."""


@dataclass(frozen=True)
class PdbAtom:
    serial: int
    record_type: str
    chain: str
    residue_number: int
    insertion_code: str
    residue_name: str
    atom_name: str
    altloc: str
    element: str
    xyz: tuple[float, float, float]

    def identity(self) -> tuple[str, int, str, str, str, str]:
        return (self.chain, self.residue_number, self.insertion_code,
                self.residue_name, self.atom_name, self.altloc)

    def protein_identity(self) -> tuple[str, int, str, str, str, str]:
        return (self.chain, self.residue_number, self.insertion_code,
                self.residue_name, self.atom_name, self.element)

    def to_json(self) -> dict[str, Any]:
        return {
            "serial": self.serial,
            "record_type": self.record_type,
            "chain": self.chain,
            "residue_number": self.residue_number,
            "insertion_code": self.insertion_code,
            "residue_name": self.residue_name,
            "atom_name": self.atom_name,
            "altloc": self.altloc,
            "element": self.element,
            "coordinates_angstrom": list(self.xyz),
        }


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def relpath(path: Path) -> str:
    return path.resolve().relative_to(ROOT).as_posix()


def research_path(path: Path) -> str:
    return path.resolve().relative_to(RESEARCH).as_posix()


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise AuditError(f"Expected JSON object in {relpath(path)}")
    return value


def link_endpoint(line: str, second: bool = False) -> dict[str, Any]:
    start = 42 if second else 12
    try:
        return {
            "atom_name": line[start:start + 4].strip(),
            "altloc": line[start + 4].strip(),
            "residue_name": line[start + 5:start + 8].strip(),
            "chain": line[start + 9].strip(),
            "residue_number": int(line[start + 10:start + 14]),
            "insertion_code": line[start + 14].strip(),
        }
    except (IndexError, ValueError) as exc:
        raise AuditError(f"Malformed PDB LINK endpoint: {line.rstrip()}") from exc


def endpoint_identity(endpoint: dict[str, Any]) -> tuple[str, int, str, str, str, str]:
    return (
        str(endpoint["chain"]), int(endpoint["residue_number"]),
        str(endpoint["insertion_code"]), str(endpoint["residue_name"]),
        str(endpoint["atom_name"]), str(endpoint["altloc"]),
    )


def parse_pdb(path: Path) -> dict[str, Any]:
    atoms: dict[int, PdbAtom] = {}
    conect_pairs: list[tuple[int, int]] = []
    link_records: list[tuple[str, dict[str, Any], dict[str, Any]]] = []
    try:
        lines = path.read_text(encoding="ascii").splitlines()
    except (OSError, UnicodeError) as exc:
        raise AuditError(f"Cannot read {relpath(path)}: {exc}") from exc

    for line_no, line in enumerate(lines, start=1):
        record = line[:6].strip()
        if record in {"ATOM", "HETATM"}:
            try:
                atom = PdbAtom(
                    serial=int(line[6:11]),
                    record_type=record,
                    chain=line[21:22].strip(),
                    residue_number=int(line[22:26]),
                    insertion_code=line[26:27].strip(),
                    residue_name=line[17:20].strip(),
                    atom_name=line[12:16].strip(),
                    altloc=line[16:17].strip(),
                    element=line[76:78].strip().upper(),
                    xyz=(float(line[30:38]), float(line[38:46]), float(line[46:54])),
                )
            except (IndexError, ValueError) as exc:
                raise AuditError(f"Malformed atom record at {relpath(path)}:{line_no}") from exc
            if atom.serial in atoms:
                raise AuditError(f"Duplicate PDB serial {atom.serial} in {relpath(path)}")
            if not atom.atom_name or not atom.element:
                raise AuditError(f"Missing atom name/element at {relpath(path)}:{line_no}")
            if not all(math.isfinite(value) for value in atom.xyz):
                raise AuditError(f"Non-finite coordinates at {relpath(path)}:{line_no}")
            atoms[atom.serial] = atom
        elif line.startswith("CONECT"):
            fields = [line[index:index + 5] for index in range(6, len(line), 5)]
            try:
                serials = [int(field) for field in fields if field.strip()]
            except ValueError as exc:
                raise AuditError(f"Malformed CONECT at {relpath(path)}:{line_no}") from exc
            if serials:
                conect_pairs.extend((serials[0], neighbor) for neighbor in serials[1:])
        elif line.startswith("LINK  "):
            link_records.append((line, link_endpoint(line), link_endpoint(line, second=True)))

    edges: set[tuple[int, int]] = set()
    for left, right in conect_pairs:
        if left not in atoms or right not in atoms:
            raise AuditError(f"CONECT references absent serial {left}-{right} in {relpath(path)}")
        if left != right:
            edges.add(tuple(sorted((left, right))))

    atoms_by_identity: dict[tuple[str, int, str, str, str, str], list[int]] = defaultdict(list)
    for serial, atom in atoms.items():
        atoms_by_identity[atom.identity()].append(serial)
    resolved_links = []
    for raw, left, right in link_records:
        resolved = []
        for endpoint in (left, right):
            matches = atoms_by_identity.get(endpoint_identity(endpoint), [])
            if len(matches) != 1:
                raise AuditError(f"LINK endpoint resolves to {len(matches)} atoms in {relpath(path)}")
            resolved.append(matches[0])
        if resolved[0] == resolved[1]:
            raise AuditError(f"LINK resolves both ends to one serial in {relpath(path)}: {raw.rstrip()}")
        pair = tuple(sorted(resolved))
        edges.add(pair)
        resolved_links.append({"record": raw.rstrip(), "serials": list(pair)})
    return {"atoms_by_serial": atoms, "adjacency_edges": edges, "resolved_links": resolved_links}


def check_digest(path: Path, expected: str, label: str) -> dict[str, str]:
    observed = sha256_file(path)
    if observed != expected:
        raise AuditError(f"SHA-256 mismatch for {label}: expected {expected}; observed {observed}")
    return {"path": relpath(path), "sha256": observed, "status": "pass"}


def walk_hash_records(value: Any):
    if isinstance(value, dict):
        if value.get("path") and value.get("observed_sha256"):
            yield value["path"], value["observed_sha256"]
        for child in value.values():
            yield from walk_hash_records(child)
    elif isinstance(value, list):
        for child in value:
            yield from walk_hash_records(child)


def load_inputs() -> tuple[dict[str, dict[str, Any]], list[dict[str, str]]]:
    loaded = {key: load_json(path) for key, path in PATHS.items() if path.suffix == ".json"}
    full_graph = loaded["full_graph"]
    graph_hash = sha256_file(PATHS["full_graph"])
    charge_hash = loaded["charge_ledger"].get("run", {}).get("input_graph_sha256")
    hybrid_hash = loaded["hybrid_audit"].get("inputs", {}).get("full_conjugate_graph", {}).get("sha256")
    if not charge_hash or charge_hash != hybrid_hash or graph_hash != charge_hash:
        raise AuditError("Full graph hash disagrees with the independent charge-ledger/hybrid records")
    if len(full_graph.get("atoms", [])) != EXPECTED_ATOMS or len(full_graph.get("bonds", [])) != EXPECTED_BONDS:
        raise AuditError("Full graph must contain 507 heavy atoms and 527 heavy-atom bonds")
    if full_graph.get("graph", {}).get("heavy_atom_count") != EXPECTED_ATOMS:
        raise AuditError("Full graph heavy-atom summary is not 507")

    hash_checks = [check_digest(PATHS["full_graph"], graph_hash, "full graph JSON")]
    for source_path, expected in sorted(set(walk_hash_records(full_graph.get("input_hashes", {})))):
        path = ROOT / source_path
        if not path.is_file():
            raise AuditError(f"Recorded full-graph input is missing: {source_path}")
        hash_checks.append(check_digest(path, expected, f"full-graph input {source_path}"))
    for key in ("reconciliation", "b1_topology", "b29_topology", "boundary_audit", "charge_ledger", "hybrid_audit"):
        hash_checks.append({"path": relpath(PATHS[key]), "sha256": sha256_file(PATHS[key]), "status": "observed"})
    return loaded, hash_checks


def records_by_state(items: list[dict[str, Any]], label: str) -> dict[str, dict[str, Any]]:
    if len(items) != len(STATE_ORDER):
        raise AuditError(f"{label} must contain exactly two state records")
    result = {}
    for item in items:
        state = item.get("state")
        if state not in STATE_ORDER or state in result:
            raise AuditError(f"{label} has an unknown or duplicate state record: {state!r}")
        result[state] = item
    if set(result) != set(STATE_ORDER):
        raise AuditError(f"{label} must contain active_open and inactive_closed")
    return result


def extract_component_mappings(audit: dict[str, Any], component: str) -> dict[str, list[dict[str, Any]]]:
    result = {}
    if component == "B1":
        states = records_by_state(audit.get("candidate_states", []), "B1 candidate states")
        for state in STATE_ORDER:
            mappings = states[state].get("all_exact_mappings", [])
            if states[state].get("exact_mapping_count") != 1 or len(mappings) != 1:
                raise AuditError(f"B1 must have one exact mapping in {state}")
            result[state] = mappings
    else:
        states = records_by_state(audit.get("states", []), "B29 candidate states")
        for state in STATE_ORDER:
            comparison = states[state].get("complete_product_exact_graph_comparison", {})
            mappings = comparison.get("exact_mappings", [])
            if comparison.get("exact_mapping_count") != EXPECTED_B29_MAPPINGS or len(mappings) != EXPECTED_B29_MAPPINGS:
                raise AuditError(f"B29 must have 16 exact mappings in {state}")
            result[state] = mappings
    for state, mappings in result.items():
        mapping_ids = [mapping.get("mapping_id") for mapping in mappings]
        if any(not isinstance(mapping_id, str) or not mapping_id for mapping_id in mapping_ids):
            raise AuditError(f"{component} mappings in {state} must all have non-empty string IDs")
        if len(set(mapping_ids)) != len(mapping_ids):
            raise AuditError(f"{component} mappings in {state} contain duplicate mapping IDs")
    if set(result) != set(STATE_ORDER):
        raise AuditError(f"{component} audit does not contain both candidate states")
    return result


def graph_indexes(graph_atoms: list[dict[str, Any]]) -> tuple[dict[str, dict[str, Any]], dict[tuple[str, int], str]]:
    nodes = {atom["node_id"]: atom for atom in graph_atoms}
    if len(nodes) != len(graph_atoms):
        raise AuditError("Duplicate full-graph node IDs")
    indices = [atom.get("atom_index_zero_based") for atom in graph_atoms]
    if any(type(index) is not int for index in indices) or sorted(indices) != list(range(len(graph_atoms))):
        raise AuditError("Full-graph atom indices must be unique, contiguous, and zero-based")
    source_nodes = {}
    for atom in graph_atoms:
        provenance = atom.get("provenance", {})
        if provenance.get("origin") == "mapped_source_fragment":
            key = (str(provenance["component"]), int(provenance["source_map_id"]))
            if key in source_nodes:
                raise AuditError(f"Duplicate full-graph source map {key}")
            source_nodes[key] = atom["node_id"]
    expected = {("B1", map_id) for map_id in range(2, 20)} | {
        ("B29", map_id) for map_id in range(1, 92)
    }
    if set(source_nodes) != expected:
        raise AuditError("Full graph mapped-source node inventory differs from B1 maps 2-19 and B29 maps 1-91")
    return nodes, source_nodes


def map_protein_atoms(
    graph_atoms: list[dict[str, Any]],
    heavy_atoms: list[PdbAtom],
) -> dict[str, PdbAtom]:
    by_identity: dict[tuple[str, int, str, str, str, str], list[PdbAtom]] = defaultdict(list)
    for atom in heavy_atoms:
        by_identity[atom.protein_identity()].append(atom)
    mapping = {}
    for graph_atom in graph_atoms:
        provenance = graph_atom.get("provenance", {})
        if provenance.get("origin") != "insulin_chain_template":
            continue
        key = (
            str(provenance["chain"]), int(provenance["residue_number"]), "",
            str(provenance["residue_name"]), str(provenance["atom_name"]),
            str(graph_atom["element"]).upper(),
        )
        matches = by_identity.get(key, [])
        if len(matches) != 1:
            raise AuditError(
                f"Protein graph node {graph_atom['node_id']} identity {key} resolves to {len(matches)} PDB atoms"
            )
        mapping[graph_atom["node_id"]] = matches[0]
    if len(mapping) != EXPECTED_PROTEIN_ATOMS:
        raise AuditError(f"Mapped {len(mapping)} insulin-template atoms; expected {EXPECTED_PROTEIN_ATOMS}")
    return mapping


def upstream_map_serials(component: str, mapping: dict[str, Any]) -> dict[int, int]:
    if component == "B1":
        records = mapping.get("map_id_to_target", {})
        serials = {
            int(map_id): int(record["target_atom"]["serial"])
            for map_id, record in records.items()
        }
        expected_maps = set(range(1, 20))
    else:
        records = mapping.get("map_id_to_pdb_serial", {})
        serials = {int(map_id): int(serial) for map_id, serial in records.items()}
        expected_maps = set(range(1, 93))
    if set(records) != {str(map_id) for map_id in expected_maps} or set(serials) != expected_maps:
        raise AuditError(f"{component} upstream mapping {mapping.get('mapping_id')} has an incomplete map-ID set")
    return serials


def compose_component_mapping(
    component: str,
    upstream_mapping: dict[str, Any],
    graph_nodes: dict[str, dict[str, Any]],
    source_nodes: dict[tuple[str, int], str],
    parsed: dict[str, Any],
    node_to_atom: dict[str, PdbAtom],
) -> dict[str, int]:
    serials = upstream_map_serials(component, upstream_mapping)
    endpoint_map = 1 if component == "B1" else 92
    endpoint_identity = (
        ("B", 1, "", "PHE", "N", "N")
        if component == "B1"
        else ("B", 29, "", "LYS", "NZ", "N")
    )
    for map_id, serial in sorted(serials.items()):
        try:
            atom = parsed["atoms_by_serial"][serial]
        except KeyError as exc:
            raise AuditError(f"{component} map {map_id} refers to absent PDB serial {serial}") from exc
        if map_id == endpoint_map:
            if atom.protein_identity() != endpoint_identity:
                raise AuditError(f"{component} endpoint map {map_id} is not the expected insulin atom")
            continue
        node_id = source_nodes.get((component, map_id))
        if node_id is None:
            raise AuditError(f"No full-graph node for {component} source map {map_id}")
        graph_atom = graph_nodes[node_id]
        if atom.element != graph_atom["element"]:
            raise AuditError(
                f"{component} map {map_id} element {atom.element} differs from graph element {graph_atom['element']}"
            )
        if component == "B1":
            source_element = upstream_mapping["map_id_to_target"][str(map_id)]["source_element"]
            if source_element != atom.element:
                raise AuditError(f"B1 map {map_id} element differs from its topology-audit record")
        if node_id in node_to_atom:
            raise AuditError(f"Full-graph node {node_id} received more than one coordinate atom")
        node_to_atom[node_id] = atom
    return {str(map_id): serial for map_id, serial in sorted(serials.items())}


def validate_candidate_inventory(
    state: str,
    parsed: dict[str, Any],
    mapped_serials: set[int],
) -> dict[str, Any]:
    heavy_atoms = [atom for atom in parsed["atoms_by_serial"].values() if atom.element not in {"H", "D"}]
    bgl_atoms = [atom for atom in heavy_atoms if atom.residue_name == "BGL"]
    bgl_residues = {
        (atom.chain, atom.residue_number, atom.insertion_code, atom.residue_name)
        for atom in bgl_atoms
    }
    if state == "active_open":
        if len(bgl_atoms) != EXPECTED_BGL_ATOMS or len(bgl_residues) != 1:
            raise AuditError("Open candidate must contain one separate 12-heavy-atom BGL residue")
    elif bgl_atoms:
        raise AuditError("Closed candidate unexpectedly contains BGL")
    bgl_serials = {atom.serial for atom in bgl_atoms}
    for edge in parsed["adjacency_edges"]:
        if (edge[0] in bgl_serials and edge[1] in mapped_serials) or (
            edge[1] in bgl_serials and edge[0] in mapped_serials
        ):
            raise AuditError("BGL has a CONECT/LINK edge to the insulin–conjugate graph")

    candidate_serials = {atom.serial for atom in heavy_atoms}
    extra_serials = candidate_serials - mapped_serials
    missing_serials = mapped_serials - candidate_serials
    if extra_serials != bgl_serials or missing_serials:
        raise AuditError(
            f"Candidate inventory mismatch: {len(extra_serials)} extra heavy atoms and "
            f"{len(missing_serials)} mapped serials absent from the heavy-atom inventory"
        )
    mapped_elements = Counter(atom.element for atom in heavy_atoms if atom.serial in mapped_serials)
    return {
        "candidate_heavy_atom_count": len(heavy_atoms),
        "mapped_graph_heavy_atom_count": len(mapped_serials),
        "mapped_element_counts": dict(sorted(mapped_elements.items())),
        "unlinked_extra_component": {
            "residue_name": "BGL",
            "heavy_atom_count": len(bgl_atoms),
            "residue_identities": [list(item) for item in sorted(bgl_residues)],
            "atom_serials": sorted(bgl_serials),
            "adjacent_to_mapped_graph": False,
        },
        "other_extra_heavy_atoms": [],
        "missing_graph_heavy_atoms": [],
    }


def geometry_records(
    bonds: list[dict[str, Any]],
    graph_nodes: dict[str, dict[str, Any]],
    node_to_atom: dict[str, PdbAtom],
    geometry_method: dict[str, Any],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], int]:
    radii = geometry_method["basis"]["reported_covalent_radii_angstrom_used"]
    offsets = geometry_method["heuristic_bond_class_offsets_from_radii_sum_angstrom"]
    observations = []
    flags = []
    screened_count = 0
    for bond in bonds:
        left_id, right_id = bond["atom_node_ids"]
        left_graph, right_graph = graph_nodes[left_id], graph_nodes[right_id]
        left_atom, right_atom = node_to_atom[left_id], node_to_atom[right_id]
        distance = math.sqrt(sum((left_atom.xyz[i] - right_atom.xyz[i]) ** 2 for i in range(3)))
        if not math.isfinite(distance):
            raise AuditError(f"Non-finite mapped bond distance for {bond['edge_id']}")
        bond_class = "aromatic" if bond.get("is_aromatic") else str(bond["bond_type"]).lower()
        interval = None
        screen = "unscreened_sulfur_edge"
        if left_graph["element"] in radii and right_graph["element"] in radii:
            if bond_class not in offsets:
                raise AuditError(f"No disclosed geometry offset for bond class {bond_class}")
            base = radii[left_graph["element"]] + radii[right_graph["element"]]
            interval = [base + offsets[bond_class][0], base + offsets[bond_class][1]]
            screen = "in_range" if interval[0] <= distance <= interval[1] else "out_of_range"
            screened_count += 1
        elif "S" not in {left_graph["element"], right_graph["element"]}:
            raise AuditError(f"No disclosed geometry screen for source bond {bond['edge_id']}")
        record = {
            "edge_id": bond["edge_id"],
            "graph_node_ids": [left_id, right_id],
            "source_bond_type": bond["bond_type"],
            "source_bond_order": bond["bond_order"],
            "pdb_serials": [left_atom.serial, right_atom.serial],
            "element_pair": [left_graph["element"], right_graph["element"]],
            "distance_angstrom": round(distance, 4),
            "broad_heuristic_interval_angstrom": [round(value, 4) for value in interval] if interval else None,
            "geometry_screen_result": screen,
        }
        observations.append(record)
        if screen == "out_of_range":
            flags.append(record)
    return observations, flags, screened_count


def build_source_bond_ledger(
    bonds: list[dict[str, Any]],
    graph_nodes: dict[str, dict[str, Any]],
    geometry_method: dict[str, Any],
) -> list[dict[str, Any]]:
    radii = geometry_method["basis"]["reported_covalent_radii_angstrom_used"]
    offsets = geometry_method["heuristic_bond_class_offsets_from_radii_sum_angstrom"]
    edge_ids = [bond.get("edge_id") for bond in bonds]
    if any(not isinstance(edge_id, str) or not edge_id for edge_id in edge_ids):
        raise AuditError("Full-graph source bonds must have non-empty string edge IDs")
    if len(set(edge_ids)) != len(edge_ids):
        raise AuditError("Full graph contains duplicate source-bond edge IDs")
    ledger = []
    for bond in bonds:
        left_id, right_id = bond["atom_node_ids"]
        left, right = graph_nodes[left_id], graph_nodes[right_id]
        bond_class = "aromatic" if bond.get("is_aromatic") else str(bond["bond_type"]).lower()
        interval = None
        if left["element"] in radii and right["element"] in radii:
            if bond_class not in offsets:
                raise AuditError(f"No disclosed geometry offset for bond class {bond_class}")
            base = radii[left["element"]] + radii[right["element"]]
            interval = [base + offsets[bond_class][0], base + offsets[bond_class][1]]
        elif "S" not in {left["element"], right["element"]}:
            raise AuditError(f"No disclosed geometry screen for source bond {bond['edge_id']}")
        ledger.append({
            "edge_id": bond["edge_id"],
            "graph_node_ids": [left_id, right_id],
            "source_bond_type": bond["bond_type"],
            "source_bond_order": bond["bond_order"],
            "element_pair": [left["element"], right["element"]],
            "broad_heuristic_interval_angstrom": [round(value, 4) for value in interval] if interval else None,
            "geometry_screen_applicable": interval is not None,
            "geometry_screen_status": "heuristic_screen" if interval is not None else "unscreened_sulfur_edge",
        })
    return ledger


def compose_state(
    state: str,
    graph: dict[str, Any],
    parsed: dict[str, Any],
    b1_mappings: list[dict[str, Any]],
    b29_mappings: list[dict[str, Any]],
    geometry_method: dict[str, Any],
) -> dict[str, Any]:
    if len(b1_mappings) != 1 or len(b29_mappings) != EXPECTED_B29_MAPPINGS:
        raise AuditError(f"Unexpected upstream mapping counts for {state}")
    graph_nodes, source_nodes = graph_indexes(graph["atoms"])
    heavy_atoms = [atom for atom in parsed["atoms_by_serial"].values() if atom.element not in {"H", "D"}]
    protein_mapping = map_protein_atoms(graph["atoms"], heavy_atoms)
    b1_anchor = next(
        atom["node_id"] for atom in graph["atoms"]
        if atom["provenance"].get("origin") == "insulin_chain_template"
        and atom["provenance"].get("chain") == "B"
        and atom["provenance"].get("residue_number") == 1
        and atom["provenance"].get("atom_name") == "N"
    )
    b29_anchor = next(
        atom["node_id"] for atom in graph["atoms"]
        if atom["provenance"].get("origin") == "insulin_chain_template"
        and atom["provenance"].get("chain") == "B"
        and atom["provenance"].get("residue_number") == 29
        and atom["provenance"].get("atom_name") == "NZ"
    )

    full_mappings = []
    graph_node_order = sorted(graph_nodes, key=lambda key: graph_nodes[key]["atom_index_zero_based"])
    common_serials: set[int] | None = None
    for b1 in b1_mappings:
        for b29 in b29_mappings:
            node_to_atom = dict(protein_mapping)
            b1_serials = compose_component_mapping("B1", b1, graph_nodes, source_nodes, parsed, node_to_atom)
            b29_serials = compose_component_mapping("B29", b29, graph_nodes, source_nodes, parsed, node_to_atom)
            if node_to_atom[b1_anchor].serial != int(b1_serials["1"]):
                raise AuditError("B1 source map 1 does not resolve to the PheB1 N graph node")
            if node_to_atom[b29_anchor].serial != int(b29_serials["92"]):
                raise AuditError("B29 source map 92 does not resolve to the LysB29 NZ graph node")
            if set(node_to_atom) != set(graph_nodes):
                raise AuditError(f"{state} mapping does not cover all 507 full-graph nodes")
            serials = [atom.serial for atom in node_to_atom.values()]
            if len(serials) != EXPECTED_ATOMS or len(set(serials)) != EXPECTED_ATOMS:
                raise AuditError(f"{state} mapping does not use 507 distinct candidate serials")
            serial_set = set(serials)
            if common_serials is None:
                common_serials = serial_set
            elif common_serials != serial_set:
                raise AuditError(f"{state} mappings do not cover the same target-heavy-atom set")

            bonds, flags, screened_count = geometry_records(
                graph["bonds"], graph_nodes, node_to_atom, geometry_method
            )
            if screened_count != 518:
                raise AuditError(f"Expected 518 C/N/O geometry observations; found {screened_count}")
            full_mappings.append({
                "mapping_id": f"{b1['mapping_id']}+{b29['mapping_id']}",
                "b1_mapping_id": b1["mapping_id"],
                "b29_mapping_id": b29["mapping_id"],
                "b1_source_map_to_pdb_serial": b1_serials,
                "b29_source_map_to_pdb_serial": b29_serials,
                "coverage": {
                    "graph_node_count": len(graph_node_order),
                    "unique_pdb_serial_count": len(serial_set),
                    "all_graph_nodes_mapped_once": True,
                    "all_coordinates_finite": True,
                },
                # These vectors align to top-level graph_atoms and source_bonds.
                # Candidate atom identity and XYZ are stored once by serial in the
                # state ledger so all 16 mappings remain complete without repeating
                # immutable PDB records thousands of times.
                "pdb_serial_by_graph_atom_index": [node_to_atom[node_id].serial for node_id in graph_node_order],
                "bond_distance_observations_angstrom": [item["distance_angstrom"] for item in bonds],
                "geometry_screen": {
                    "screened_CNO_edge_count": screened_count,
                    "unscreened_sulfur_edge_count": EXPECTED_BONDS - screened_count,
                    "out_of_range_edge_count": len(flags),
                    "out_of_range_edge_ids": [item["edge_id"] for item in flags],
                    "distances_used_to_choose_or_rank_mapping": False,
                },
            })

    inventory = validate_candidate_inventory(state, parsed, common_serials or set())
    expected_elements = Counter(atom["element"] for atom in graph["atoms"])
    if inventory["mapped_element_counts"] != dict(sorted(expected_elements.items())):
        raise AuditError(f"{state} element inventory differs from full graph")
    return {
        "state": state,
        "coordinate_candidate": relpath(PATHS[STATE_PDB_KEY[state]]),
        "mapping_count": len(full_mappings),
        "protein_template_atom_count": EXPECTED_PROTEIN_ATOMS,
        "b1_conjugate_atom_count": EXPECTED_B1_ATOMS,
        "b29_conjugate_atom_count": EXPECTED_B29_ATOMS,
        "coordinate_inventory": inventory,
        "coordinate_atoms": [atom.to_json() for atom in sorted(heavy_atoms, key=lambda atom: atom.serial)],
        "mappings": full_mappings,
        "hydrogen_record_context": {
            "source_audit": relpath(PATHS["boundary_audit"]),
            "pdb_hydrogens_mapped": False,
            "pheB1_hydrogen_ownership_resolved": False,
            "lyx_hydrogen_neighbor_differences_resolved": False,
        },
    }


def build_report(preflight: bool = False) -> dict[str, Any]:
    report: dict[str, Any] = {
        "schema": "dia-health-one.nnc2215-full-coordinate-mapping.v2",
        "status": "running",
        "scope": "heavy-atom graph-to-coordinate mapping ledger for model-built candidates",
        "script": relpath(Path(__file__)),
        "script_sha256": sha256_file(Path(__file__)),
        "inputs": {},
        "graph_atoms": [],
        "source_bonds": [],
        "mapping_record_layout": {
            "pdb_serial_by_graph_atom_index": "Aligns with graph_atoms ordered by graph_atom_index_zero_based; look up serial in that state's coordinate_atoms for PDB identity and XYZ.",
            "bond_distance_observations_angstrom": "Aligns by list position with source_bonds, which preserves the full graph's bond-array order; each value is in angstroms.",
        },
        "states": [],
        "diagnostics": [],
        "assessment": {},
        "limits": [
            "Both coordinate candidates are model-built derivatives, not experimentally resolved NNC2215 structures.",
            "Sixteen exact B29 graph-to-coordinate mappings remain per state; the audit does not select or rank them.",
            "PDB adjacency and coordinates do not independently establish bond orders, stereochemistry, protonation, or hydrogen ownership.",
            "The geometry intervals are broad heuristic screens, not evidence of bond order or mapping validity.",
            "No complete force field, OpenMM target system, target MD, affinity result, or clinical conclusion is established.",
        ],
    }
    before_hashes = {}
    for state, path_key in STATE_PDB_KEY.items():
        try:
            before_hashes[state] = sha256_file(PATHS[path_key])
        except OSError:
            before_hashes[state] = None
    try:
        loaded, hash_checks = load_inputs()
        reconciliation = records_by_state(loaded["reconciliation"].get("models", []), "reconciliation models")
        b1_by_state = extract_component_mappings(loaded["b1_topology"], "B1")
        b29_by_state = extract_component_mappings(loaded["b29_topology"], "B29")
        graph = loaded["full_graph"]
        report["inputs"] = {
            "full_graph": {"path": relpath(PATHS["full_graph"]), "sha256": sha256_file(PATHS["full_graph"]),
                           "heavy_atom_count": len(graph["atoms"]), "heavy_bond_count": len(graph["bonds"])},
            "audits": {
                key: {"path": relpath(PATHS[key]), "sha256": sha256_file(PATHS[key])}
                for key in ("reconciliation", "b1_topology", "b29_topology", "boundary_audit", "charge_ledger", "hybrid_audit")
            },
            "hash_checks": hash_checks,
        }
        if preflight:
            candidate_checks = {}
            for state, path_key in STATE_PDB_KEY.items():
                expected = reconciliation[state]["output_sha256"]
                candidate_checks[state] = {"observed": before_hashes[state], "expected": expected,
                                           "matches": before_hashes[state] == expected}
                if before_hashes[state] != expected:
                    raise AuditError(f"{state} candidate hash differs from reconciliation")
            report["preflight"] = {
                "full_graph_heavy_atom_count": len(graph["atoms"]),
                "full_graph_bond_count": len(graph["bonds"]),
                "b1_mapping_count_by_state": {state: len(b1_by_state[state]) for state in STATE_ORDER},
                "b29_mapping_count_by_state": {state: len(b29_by_state[state]) for state in STATE_ORDER},
                "candidate_hash_checks": candidate_checks,
            }
            report["status"] = "preflight_passed"
            return report

        report["graph_atoms"] = [
            {
                "graph_node_id": atom["node_id"],
                "graph_atom_index_zero_based": atom["atom_index_zero_based"],
                "element": atom["element"],
                "provenance": atom["provenance"],
            }
            for atom in sorted(graph["atoms"], key=lambda item: item["atom_index_zero_based"])
        ]

        boundary_pdbs = loaded["boundary_audit"].get("inputs", {}).get("candidate_pdb_before_after_hashes", {})
        b1_states = records_by_state(loaded["b1_topology"].get("candidate_states", []), "B1 candidate states")
        b29_states = records_by_state(loaded["b29_topology"].get("states", []), "B29 candidate states")
        geometry_method = loaded["b29_topology"].get("geometry_screen_method")
        if not geometry_method:
            raise AuditError("B29 geometry-screen method is missing")
        graph_nodes, _ = graph_indexes(graph["atoms"])
        report["source_bonds"] = build_source_bond_ledger(graph["bonds"], graph_nodes, geometry_method)

        for state in STATE_ORDER:
            path = PATHS[STATE_PDB_KEY[state]]
            model = reconciliation[state]
            if model.get("output_file") != research_path(path):
                raise AuditError(f"Reconciliation names the wrong candidate path for {state}")
            expected_hashes = {
                "reconciliation": model.get("output_sha256"),
                "B1 topology audit": b1_states[state].get("candidate_pdb_sha256"),
                "B29 topology audit": b29_states[state].get("coordinate_candidate_sha256_expected_by_reconciliation"),
            }
            boundary_hash = boundary_pdbs.get(state, {}).get("sha256_expected_by_reconciliation")
            if boundary_hash:
                expected_hashes["attachment-boundary audit"] = boundary_hash
            if any(not value or value != before_hashes[state] for value in expected_hashes.values()):
                raise AuditError(f"Candidate PDB hash disagrees with upstream records for {state}")
            parsed = parse_pdb(path)
            completed = compose_state(
                state, graph, parsed, b1_by_state[state], b29_by_state[state], geometry_method
            )
            completed["coordinate_candidate_sha256"] = before_hashes[state]
            report["states"].append(completed)
            if sha256_file(path) != before_hashes[state]:
                raise AuditError(f"Candidate PDB changed during audit for {state}")

        report["status"] = "coordinate_mapping_coverage_passed"
        report["assessment"] = {
            "both_candidate_states_mapped": len(report["states"]) == 2,
            "full_mapping_count_by_state": {item["state"]: item["mapping_count"] for item in report["states"]},
            "each_full_mapping_covers_507_unique_graph_nodes_and_pdb_serials": all(
                mapping["coverage"]["graph_node_count"] == EXPECTED_ATOMS
                and mapping["coverage"]["unique_pdb_serial_count"] == EXPECTED_ATOMS
                for item in report["states"] for mapping in item["mappings"]
            ),
            "b29_ambiguity_preserved": all(item["mapping_count"] == EXPECTED_B29_MAPPINGS for item in report["states"]),
            "coordinates_used_to_select_or_rank_mappings": False,
            "whole_system_microstate_assigned": False,
            "complete_force_field_validated": False,
            "openmm_system_created": False,
            "target_md_run": False,
            "clinical_readiness_established": False,
        }
    except Exception as exc:
        report["status"] = "partial" if report["states"] else "failed"
        report["diagnostics"].append({"type": type(exc).__name__, "message": str(exc)})
        report["assessment"] = {
            "complete_coordinate_mapping_established": False,
            "incomplete_assignments_published_as_complete": False,
            "whole_system_microstate_assigned": False,
            "complete_force_field_validated": False,
            "openmm_system_created": False,
            "target_md_run": False,
            "clinical_readiness_established": False,
        }
    finally:
        after_hashes = {}
        for state, path_key in STATE_PDB_KEY.items():
            try:
                after = sha256_file(PATHS[path_key])
                unchanged = before_hashes[state] == after
                after_hashes[state] = {"before": before_hashes[state], "after": after, "unchanged": unchanged}
                if not unchanged:
                    report["diagnostics"].append({"type": "CandidateMutationError", "message": f"{state} candidate changed during audit"})
                    report["status"] = "partial" if report["states"] else "failed"
                    report["assessment"]["complete_coordinate_mapping_established"] = False
            except OSError as exc:
                after_hashes[state] = {"before": before_hashes[state], "after": None, "unchanged": False}
                report["diagnostics"].append({"type": "HashReadError", "message": str(exc)})
                report["status"] = "partial" if report["states"] else "failed"
        report["candidate_pdb_before_after_hashes"] = after_hashes
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preflight", action="store_true", help="validate graph inputs, hashes, and mapping counts only")
    args = parser.parse_args()
    report = build_report(preflight=args.preflight)
    if args.preflight:
        print(json.dumps(report.get("preflight", report.get("diagnostics", [])), indent=2, sort_keys=True))
        return 0 if report["status"] == "preflight_passed" else 1
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_PATH.write_text(json.dumps(report, indent=1, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    print(f"Status: {report['status']}")
    if report["diagnostics"]:
        print(f"Diagnostics: {len(report['diagnostics'])}; see {relpath(OUTPUT_PATH)}")
    else:
        print(f"Wrote {relpath(OUTPUT_PATH)} with 16 mappings per state")
    return 0 if report["status"] == "coordinate_mapping_coverage_passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
