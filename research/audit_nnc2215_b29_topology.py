#!/usr/bin/env python3
"""Audit the complete source-mapped B29 product against corrected NNC2215 PDBs.

PDB CONECT/LINK records are treated as adjacency-only evidence. Source bond
labels are transferred after exact graph mapping and are used only to annotate
a broad coordinate geometry screen; distances are not mapping or bond-order
evidence.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import platform
import re
from typing import Any

import networkx as nx
from rdkit import Chem, rdBase
from rdkit.Chem import rdMolDescriptors


ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
PRODUCT_MAP_PATH = DATA / "nnc2215-b29-linker-roof-product-map.json"
RECONCILIATION_PATH = DATA / "nnc2215-carbon-reconciliation-candidate.json"
PRIOR_ROOF_MAPPING_PATH = DATA / "nnc2215-roof-graph-coordinate-mapping-audit.json"
LINKER_MAPPING_PATH = DATA / "nnc2215-linker-mapping-audit.json"
ROOF_GRAPH_AUDIT_PATH = DATA / "nnc2215-roof-macrocycle-graph-audit.json"
ROOF_MAP_PATH = DATA / "nnc2215-roof-macrocycle-atom-map.json"
SOURCE_CHEMISTRY_PATH = DATA / "nnc2215-chemistry-source-audit.json"
OUTPUT_PATH = DATA / "nnc2215-b29-topology-audit.json"

EXPECTED_PRODUCT_ELEMENTS = Counter({"C": 62, "N": 16, "O": 14})
EXPECTED_PRODUCT_NODES = 92
EXPECTED_PRODUCT_EDGES = 99
EXPECTED_PRODUCT_COMPONENTS = 1
EXPECTED_PRODUCT_CYCLE_RANK = 8
EXPECTED_ROOF_ELEMENTS = Counter({"C": 56, "N": 15, "O": 12})
EXPECTED_ROOF_NODES = 83
EXPECTED_ROOF_EDGES = 89
EXPECTED_ROOF_EMBEDDINGS = 16
CARBOXYLATE_MAP_GROUPS = ((52, 53, 54), (60, 61, 62), (68, 69, 70))
EXPECTED_CARBOXYLATE_MAP_GROUPS = ((52, 53, 54), (60, 61, 62), (68, 69, 70))
EXPECTED_PRODUCT_FORMULA = "C62H72N16O14"
EXPECTED_DEPROTONATED_FRAGMENT_FORMULA = "C62H69N16O14"
ASSUMED_CARBOXYLATE_PROTON_REMOVALS = 3
ASSUMED_CARBOXYLATE_LOCAL_CHARGE_DELTA = -3
STATE_ORDER = ("active_open", "inactive_closed")

# Literature radii are rounded to the values specified in the audit request.
COVALENT_RADII_ANGSTROM = {"C": 0.76, "N": 0.71, "O": 0.66}
# Heuristic offsets from the sum of single-bond covalent radii. The source
# article reports radii, not these broad per-bond screening intervals.
BOND_CLASS_OFFSETS_ANGSTROM = {
    "single": [-0.18, 0.25],
    "double": [-0.35, -0.06],
    "triple": [-0.45, -0.18],
    "aromatic": [-0.20, 0.05],
}

ROOF_ANCHOR_MAP_KEYS = {
    78: "source_macrocycle_c21_candidate",
    79: "source_macrocycle_c22_candidate",
    80: "source_macrocycle_c23_candidate",
    81: "source_azide_n_attached_to_c23_candidate",
    82: "source_triazole_n2_candidate",
    83: "source_triazole_n1_candidate",
}

HYDROGEN_OBSERVATIONS = {
    78: {"label": "roof C21", "open": ["H1", "H2"], "closed": ["H1", "H2"], "relation": "same"},
    79: {"label": "roof C22", "open": ["HX22", "HX31"], "closed": ["HCB", "HCC"], "relation": "different"},
    80: {"label": "roof C23", "open": ["HC2", "HC3"], "closed": ["HCA", "HCD"], "relation": "different"},
    81: {"label": "triazole N1", "open": [], "closed": [], "relation": "same"},
    82: {"label": "triazole N2", "open": [], "closed": [], "relation": "same"},
    83: {"label": "triazole N3", "open": [], "closed": [], "relation": "same"},
    85: {"label": "linker C4", "open": ["HC1", "HXTA"], "closed": ["HXTA"], "relation": "different"},
    86: {"label": "linker C5", "open": ["H1", "H2A"], "closed": ["H1", "H2"], "relation": "different"},
    87: {"label": "linker C6", "open": ["H2B", "H3"], "closed": ["HD3", "HD4"], "relation": "different"},
    88: {"label": "linker C7 / triazole C4", "open": [], "closed": [], "relation": "same"},
    89: {"label": "linker C8 / triazole C5", "open": ["HXTB"], "closed": ["HXTB"], "relation": "same"},
    92: {"label": "LysB29 NZ", "open": ["HZ2"], "closed": ["HZ2"], "relation": "same"},
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def subtract_formula_hydrogens(formula: str, hydrogen_count: int) -> str:
    tokens = list(re.finditer(r"([A-Z][a-z]?)(\d*)", formula))
    if not tokens or "".join(token.group(0) for token in tokens) != formula:
        raise ValueError(f"Invalid molecular formula: {formula!r}")
    element_counts = {
        token.group(1): int(token.group(2) or "1")
        for token in tokens
    }
    if element_counts.get("H", 0) < hydrogen_count:
        raise ValueError(
            f"Cannot remove {hydrogen_count} hydrogens from formula {formula}"
        )
    element_counts["H"] = element_counts.get("H", 0) - hydrogen_count
    if element_counts["H"] == 0:
        del element_counts["H"]
    element_order = ["C", "H"] + sorted(
        element for element in element_counts if element not in {"C", "H"}
    )
    return "".join(
        f"{element}{element_counts[element] if element_counts[element] != 1 else ''}"
        for element in element_order
        if element in element_counts
    )


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def project_path(relative_path: str) -> Path:
    path = (ROOT / relative_path).resolve()
    if not path.is_relative_to(ROOT):
        raise ValueError(f"Input path escapes research directory: {relative_path}")
    return path


def check_hash(path: Path, expected: str, description: str) -> str:
    observed = sha256(path)
    if observed != expected:
        raise ValueError(
            f"SHA-256 mismatch for {description}: expected {expected}, observed {observed}"
        )
    return observed


def parse_molecule(path: Path, description: str) -> Chem.Mol:
    tokens = path.read_text(encoding="utf-8").split()
    if not tokens:
        raise ValueError(f"Empty mapped SMILES file: {path}")
    molecule = Chem.MolFromSmiles(tokens[0])
    if molecule is None:
        raise ValueError(f"RDKit could not parse {description}: {path}")
    Chem.SanitizeMol(molecule)
    return molecule


def bond_label(bond: Chem.Bond) -> tuple[str, float]:
    if bond.GetIsAromatic():
        return "aromatic", 1.5
    order = bond.GetBondTypeAsDouble()
    if math.isclose(order, 1.0):
        return "single", 1.0
    if math.isclose(order, 2.0):
        return "double", 2.0
    if math.isclose(order, 3.0):
        return "triple", 3.0
    raise ValueError(f"Unsupported source bond order: {bond.GetBondType()}")


def molecule_graph(molecule: Chem.Mol, description: str) -> nx.Graph:
    graph = nx.Graph()
    for atom in molecule.GetAtoms():
        if atom.GetAtomicNum() == 1:
            continue
        map_id = atom.GetAtomMapNum()
        if map_id <= 0 or map_id in graph:
            raise ValueError(f"{description} requires unique positive heavy-atom map IDs")
        graph.add_node(map_id, map_id=map_id, element=atom.GetSymbol())
    for bond in molecule.GetBonds():
        left = molecule.GetAtomWithIdx(bond.GetBeginAtomIdx())
        right = molecule.GetAtomWithIdx(bond.GetEndAtomIdx())
        if left.GetAtomicNum() == 1 or right.GetAtomicNum() == 1:
            continue
        label, order = bond_label(bond)
        graph.add_edge(
            left.GetAtomMapNum(),
            right.GetAtomMapNum(),
            bond_class=label,
            bond_order=order,
        )
    return graph


def element_counts(graph: nx.Graph) -> dict[str, int]:
    return dict(sorted(Counter(data["element"] for _, data in graph.nodes(data=True)).items()))


def graph_summary(graph: nx.Graph) -> dict[str, Any]:
    components = nx.number_connected_components(graph)
    return {
        "heavy_atom_count": graph.number_of_nodes(),
        "heavy_element_counts": element_counts(graph),
        "adjacency_edge_count": graph.number_of_edges(),
        "connected_component_count": components,
        "cycle_rank_E_minus_V_plus_components": graph.number_of_edges()
        - graph.number_of_nodes()
        + components,
    }


def node_element_match(left: dict[str, Any], right: dict[str, Any]) -> bool:
    return left["element"] == right["element"]


def source_bond_match(left: dict[str, Any], right: dict[str, Any]) -> bool:
    return left["bond_class"] == right["bond_class"] and left["bond_order"] == right["bond_order"]


def graph_automorphisms(graph: nx.Graph, preserve_bond_labels: bool) -> list[dict[int, int]]:
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        graph,
        graph,
        node_match=node_element_match,
        edge_match=source_bond_match if preserve_bond_labels else None,
    )
    return [dict(mapping) for mapping in matcher.isomorphisms_iter()]


def carboxylate_resonance_graph(source_graph: nx.Graph) -> tuple[nx.Graph, list[dict[str, Any]]]:
    """Copy the neutral graph and normalize only the three declared carboxylate C–O pairs."""
    if CARBOXYLATE_MAP_GROUPS != EXPECTED_CARBOXYLATE_MAP_GROUPS:
        raise ValueError(
            "Carboxylate map groups differ from the approved groups: "
            f"{CARBOXYLATE_MAP_GROUPS}"
        )
    if len(CARBOXYLATE_MAP_GROUPS) != 3:
        raise ValueError(
            f"Expected exactly three carboxylate map groups; found {len(CARBOXYLATE_MAP_GROUPS)}"
        )

    normalized_graph = source_graph.copy()
    expected_edges: set[tuple[int, int]] = set()
    group_records: list[dict[str, Any]] = []
    for carbon_map_id, first_oxygen_map_id, second_oxygen_map_id in CARBOXYLATE_MAP_GROUPS:
        oxygen_map_ids = (first_oxygen_map_id, second_oxygen_map_id)
        if len(set((carbon_map_id, *oxygen_map_ids))) != 3:
            raise ValueError(f"Carboxylate group contains duplicate map IDs: {carbon_map_id, *oxygen_map_ids}")
        if carbon_map_id not in source_graph or source_graph.nodes[carbon_map_id]["element"] != "C":
            raise ValueError(f"Carboxylate map {carbon_map_id} must identify a source carbon")
        if any(
            oxygen_map_id not in source_graph
            or source_graph.nodes[oxygen_map_id]["element"] != "O"
            for oxygen_map_id in oxygen_map_ids
        ):
            raise ValueError(
                f"Carboxylate carbon map {carbon_map_id} must be paired with the declared oxygen maps {oxygen_map_ids}"
            )
        oxygen_neighbors = {
            neighbor
            for neighbor in source_graph.neighbors(carbon_map_id)
            if source_graph.nodes[neighbor]["element"] == "O"
        }
        if oxygen_neighbors != set(oxygen_map_ids) or any(
            not source_graph.has_edge(carbon_map_id, oxygen_map_id)
            for oxygen_map_id in oxygen_map_ids
        ):
            raise ValueError(
                f"Carboxylate carbon map {carbon_map_id} must be bonded to exactly declared oxygen maps {oxygen_map_ids}; "
                f"observed oxygen neighbors {sorted(oxygen_neighbors)}"
            )

        neutral_bond_labels = {
            (
                source_graph.edges[carbon_map_id, oxygen_map_id]["bond_class"],
                source_graph.edges[carbon_map_id, oxygen_map_id]["bond_order"],
            )
            for oxygen_map_id in oxygen_map_ids
        }
        if neutral_bond_labels != {("single", 1.0), ("double", 2.0)}:
            raise ValueError(
                f"Neutral carboxylate group {carbon_map_id, *oxygen_map_ids} must contain one single and one double C–O edge; "
                f"observed {sorted(neutral_bond_labels)}"
            )

        for oxygen_map_id in oxygen_map_ids:
            normalized_graph.edges[carbon_map_id, oxygen_map_id]["bond_class"] = "carboxylate_resonance"
            normalized_graph.edges[carbon_map_id, oxygen_map_id]["bond_order"] = 1.5
            expected_edges.add(tuple(sorted((carbon_map_id, oxygen_map_id))))
        group_records.append(
            {
                "carbon_map_id": carbon_map_id,
                "oxygen_map_ids": list(oxygen_map_ids),
                "equivalent_edge_map_ids": [
                    [carbon_map_id, first_oxygen_map_id],
                    [carbon_map_id, second_oxygen_map_id],
                ],
            }
        )

    changed_edges = {
        tuple(sorted((left, right)))
        for left, right in source_graph.edges
        if source_graph.edges[left, right] != normalized_graph.edges[left, right]
    }
    if changed_edges != expected_edges:
        raise ValueError(
            "Resonance normalization must change only the six declared C–O edges; "
            f"observed changed edges {sorted(changed_edges)}"
        )
    return normalized_graph, group_records


def parse_link_endpoint(line: str, first: bool) -> dict[str, str]:
    if first:
        fields = (line[12:16], line[16:17], line[17:20], line[21:22], line[22:26], line[26:27])
    else:
        fields = (line[42:46], line[46:47], line[47:50], line[51:52], line[52:56], line[56:57])
    atom_name, altloc, residue, chain, residue_number, insertion_code = (
        value.strip() for value in fields
    )
    return {
        "name": atom_name,
        "altloc": altloc,
        "residue": residue,
        "chain": chain,
        "residue_number": residue_number,
        "insertion_code": insertion_code,
    }


def link_match_key(atom: dict[str, Any]) -> tuple[str, ...]:
    return tuple(
        atom[key]
        for key in ("chain", "residue_number", "insertion_code", "residue", "altloc", "name")
    )


def endpoint_match_key(endpoint: dict[str, str]) -> tuple[str, ...]:
    return tuple(
        endpoint[key]
        for key in ("chain", "residue_number", "insertion_code", "residue", "altloc", "name")
    )


def parse_pdb(path: Path) -> dict[str, Any]:
    atoms: dict[int, dict[str, Any]] = {}
    raw_conect: list[tuple[int, int]] = []
    link_records: list[tuple[str, dict[str, str], dict[str, str]]] = []

    for line_number, line in enumerate(path.read_text(encoding="ascii").splitlines(), start=1):
        if line.startswith(("ATOM  ", "HETATM")):
            serial = int(line[6:11])
            if serial in atoms:
                raise ValueError(f"Duplicate PDB atom serial {serial} in {path}")
            element = line[76:78].strip().upper()
            if not element:
                raise ValueError(f"Missing element at PDB line {line_number} in {path}")
            atom = {
                "serial": serial,
                "record_type": line[:6].strip(),
                "name": line[12:16].strip(),
                "altloc": line[16:17].strip(),
                "residue": line[17:20].strip(),
                "chain": line[21:22].strip(),
                "residue_number": line[22:26].strip(),
                "insertion_code": line[26:27].strip(),
                "element": element,
                "coordinates_angstrom": [
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54]),
                ],
            }
            atoms[serial] = atom
        elif line.startswith("CONECT"):
            serials = [
                int(line[index : index + 5])
                for index in range(6, len(line), 5)
                if line[index : index + 5].strip()
            ]
            if serials:
                raw_conect.extend((serials[0], neighbor) for neighbor in serials[1:])
        elif line.startswith("LINK  "):
            link_records.append(
                (line, parse_link_endpoint(line, True), parse_link_endpoint(line, False))
            )

    directed_conect: set[tuple[int, int]] = set()
    for left, right in raw_conect:
        if left not in atoms or right not in atoms:
            raise ValueError(f"CONECT references absent atom serial {left}-{right} in {path}")
        if left == right:
            continue
        directed_conect.add((left, right))

    edge_sources: dict[tuple[int, int], set[str]] = defaultdict(set)
    for left, right in directed_conect:
        edge_sources[tuple(sorted((left, right)))].add("CONECT")
    conect_only_undirected_edge_count = len(edge_sources)

    resolved_links: list[dict[str, Any]] = []
    all_atoms_by_identity: dict[tuple[str, ...], list[int]] = defaultdict(list)
    for serial, atom in atoms.items():
        all_atoms_by_identity[link_match_key(atom)].append(serial)
    for raw_line, left_endpoint, right_endpoint in link_records:
        resolved_serials: list[int] = []
        for endpoint in (left_endpoint, right_endpoint):
            matches = all_atoms_by_identity.get(endpoint_match_key(endpoint), [])
            if len(matches) != 1:
                raise ValueError(
                    f"LINK endpoint {endpoint} resolves to {matches} in {path}; expected exactly one"
                )
            resolved_serials.append(matches[0])
        if resolved_serials[0] == resolved_serials[1]:
            raise ValueError(f"LINK resolves both endpoints to one atom in {path}: {raw_line}")
        pair = tuple(sorted(resolved_serials))
        edge_sources[pair].add("LINK")
        resolved_links.append(
            {
                "record": raw_line,
                "first_endpoint": left_endpoint,
                "second_endpoint": right_endpoint,
                "resolved_serials": resolved_serials,
                "edge_pair": list(pair),
            }
        )

    return {
        "atoms": atoms,
        "directed_conect": directed_conect,
        "conect_edges": {pair for pair, sources in edge_sources.items() if "CONECT" in sources},
        "edge_sources": {pair: set(sources) for pair, sources in edge_sources.items()},
        "conect_link_overlap_pair_count": sum(
            sources == {"CONECT", "LINK"} for sources in edge_sources.values()
        ),
        "resolved_links": resolved_links,
        "raw_conect_reference_count": len(raw_conect),
        "unique_directed_conect_reference_count": len(directed_conect),
        "conect_only_undirected_edge_count": conect_only_undirected_edge_count,
    }


def pdb_atom_record(atom: dict[str, Any], include_coordinates: bool = False) -> dict[str, Any]:
    record = {
        "serial": atom["serial"],
        "record_type": atom["record_type"],
        "chain": atom["chain"],
        "residue_number": atom["residue_number"],
        "insertion_code": atom["insertion_code"],
        "residue": atom["residue"],
        "altloc": atom["altloc"],
        "name": atom["name"],
        "element": atom["element"],
    }
    if include_coordinates:
        record["coordinates_angstrom"] = atom["coordinates_angstrom"]
    return record


def selected_target_graph(
    parsed: dict[str, Any], include_nz: bool
) -> tuple[nx.Graph, dict[int, dict[str, Any]], dict[str, Any]]:
    atoms = parsed["atoms"]
    selected: dict[int, dict[str, Any]] = {}
    for serial, atom in atoms.items():
        in_lyx_or_mc = (
            atom["chain"],
            atom["residue_number"],
            atom["insertion_code"],
            atom["residue"],
        ) in {
            ("B", "30", "", "LYX"),
            ("C", "1", "", "MC"),
        }
        is_nz = (
            atom["chain"] == "B"
            and atom["residue_number"] == "29"
            and atom["insertion_code"] == ""
            and atom["residue"] == "LYS"
            and atom["name"] == "NZ"
        )
        if (in_lyx_or_mc and atom["element"] not in {"H", "D"}) or (include_nz and is_nz):
            selected[serial] = atom

    if include_nz:
        nz_atoms = [
            atom
            for atom in selected.values()
            if atom["chain"] == "B"
            and atom["residue_number"] == "29"
            and atom["insertion_code"] == ""
            and atom["residue"] == "LYS"
            and atom["name"] == "NZ"
        ]
        if len(nz_atoms) != 1 or nz_atoms[0]["element"] != "N":
            raise ValueError(f"Expected exactly one heavy LYS B29 NZ; found {len(nz_atoms)}")

    graph = nx.Graph()
    for serial, atom in selected.items():
        graph.add_node(serial, serial=serial, element=atom["element"], atom=atom)
    for pair in parsed["edge_sources"]:
        left, right = pair
        if left in selected and right in selected:
            graph.add_edge(
                left,
                right,
                pdb_adjacency_sources=sorted(parsed["edge_sources"][pair]),
            )

    node_set = set(selected)
    boundary_edges = []
    for pair, sources in parsed["edge_sources"].items():
        left, right = pair
        if (left in node_set) != (right in node_set):
            selected_serial = left if left in node_set else right
            outside_serial = right if left in node_set else left
            if selected[selected_serial]["name"] == "NZ":
                outside_atom = atoms[outside_serial]
                boundary_edges.append(
                    {
                        "selected_atom": pdb_atom_record(selected[selected_serial]),
                        "outside_scope_atom": pdb_atom_record(outside_atom),
                        "adjacency_sources": sorted(sources),
                    }
                )

    component_counts = Counter(
        (atom["chain"], atom["residue_number"], atom["residue"], atom["element"])
        for atom in selected.values()
    )
    inventory = {
        f"{chain}:{residue_number}:{residue}:{element}": count
        for (chain, residue_number, residue, element), count in sorted(component_counts.items())
    }
    metadata = {
        "selected_atom_count": len(selected),
        "selected_heavy_element_counts": dict(sorted(Counter(a["element"] for a in selected.values()).items())),
        "selected_residue_element_inventory": inventory,
        "selected_connected_component_count": nx.number_connected_components(graph) if graph else 0,
        "selected_adjacency_edge_count": graph.number_of_edges(),
        "selected_cycle_rank_E_minus_V_plus_components": graph.number_of_edges()
        - graph.number_of_nodes()
        + (nx.number_connected_components(graph) if graph else 0),
        "excluded_b29_nz_boundary_edges": sorted(
            boundary_edges,
            key=lambda item: item["outside_scope_atom"]["serial"],
        ),
    }
    return graph, selected, metadata


def resolve_anchor(
    selected_atoms: dict[int, dict[str, Any]],
    source_record: dict[str, Any],
    map_id: int,
    expected_element: str,
) -> tuple[int, dict[str, Any]]:
    pdb_atom = source_record["pdb_atom"]
    expected = {
        "chain": str(pdb_atom["chain"]).strip(),
        "residue_number": str(pdb_atom["residue_number"]).strip(),
        "insertion_code": "",
        "residue": str(pdb_atom["residue"]).strip(),
        "name": str(pdb_atom["name"]).strip(),
        "altloc": "",
    }
    if str(pdb_atom["element"]).upper() != expected_element:
        raise ValueError(
            f"Candidate anchor map {map_id} element {pdb_atom['element']} != {expected_element}"
        )
    matches = [
        serial
        for serial, atom in selected_atoms.items()
        if all(atom[key] == value for key, value in expected.items())
        and atom["element"] == expected_element
    ]
    if len(matches) != 1:
        raise ValueError(f"Candidate anchor {map_id} {expected} resolves to {matches}; expected one")
    return matches[0], {
        "map_id": map_id,
        "source_anchor_key": source_record.get("source_key"),
        "candidate_source_model_serial": source_record["serial"],
        "candidate_coordinate_atom": pdb_atom_record(selected_atoms[matches[0]]),
        "resolution_rule": "Full chain/residue/insertion/altloc/atom-name/element identity; unique match required.",
    }


def roof_only_embeddings(
    roof_graph: nx.Graph,
    target: nx.Graph,
    selected_atoms: dict[int, dict[str, Any]],
    linker_state: dict[str, Any],
) -> tuple[list[dict[int, int]], dict[str, Any]]:
    anchors: dict[int, int] = {}
    anchor_details: dict[str, Any] = {}
    source_map = linker_state["candidate_source_to_model_atom_map"]
    expected_elements = {map_id: roof_graph.nodes[map_id]["element"] for map_id in ROOF_ANCHOR_MAP_KEYS}
    for map_id, source_key in ROOF_ANCHOR_MAP_KEYS.items():
        source_record = source_map[source_key]
        serial, detail = resolve_anchor(
            selected_atoms,
            {**source_record, "source_key": source_key},
            map_id,
            expected_elements[map_id],
        )
        anchors[map_id] = serial
        detail["source_graph_label"] = {
            78: "roof tail C21",
            79: "roof tail C22",
            80: "roof tail C23",
            81: "azide-derived triazole N1 attached to C23",
            82: "triazole N2",
            83: "triazole N3",
        }[map_id]
        anchor_details[str(map_id)] = detail

    def constrained_node_match(target_node: dict[str, Any], source_node: dict[str, Any]) -> bool:
        if target_node["element"] != source_node["element"]:
            return False
        fixed_serial = anchors.get(source_node["map_id"])
        return fixed_serial is None or target_node["serial"] == fixed_serial

    matcher = nx.algorithms.isomorphism.GraphMatcher(
        target,
        roof_graph,
        node_match=constrained_node_match,
    )
    embeddings: list[dict[int, int]] = []
    for target_to_source in matcher.subgraph_monomorphisms_iter():
        source_to_target = {
            source_map_id: target_serial
            for target_serial, source_map_id in target_to_source.items()
        }
        if len(source_to_target) != roof_graph.number_of_nodes():
            continue
        if any(source_to_target.get(map_id) != serial for map_id, serial in anchors.items()):
            continue
        embeddings.append(source_to_target)
    embeddings.sort(key=lambda mapping: tuple(mapping[map_id] for map_id in sorted(mapping)))
    return embeddings, {
        "fixed_candidate_anchors_by_map_id": anchor_details,
        "anchored_embedding_count": len(embeddings),
        "embeddings": [
            {
                "embedding_id": f"roof_embedding_{index:02d}",
                "map_id_to_pdb_serial": {
                    str(map_id): serial for map_id, serial in sorted(mapping.items())
                },
            }
            for index, mapping in enumerate(embeddings, start=1)
        ],
        "search_scope": "Element-labeled, edge-preserving subgraph monomorphism; extra target edges are allowed. This is a precursor check, not a complete-product assignment.",
    }


def interval_for_edge(left_element: str, right_element: str, bond_class: str) -> list[float]:
    radii_sum = COVALENT_RADII_ANGSTROM[left_element] + COVALENT_RADII_ANGSTROM[right_element]
    low_offset, high_offset = BOND_CLASS_OFFSETS_ANGSTROM[bond_class]
    return [round(radii_sum + low_offset, 2), round(radii_sum + high_offset, 2)]


def distance_angstrom(left: list[float], right: list[float]) -> float:
    return math.sqrt(sum((left[index] - right[index]) ** 2 for index in range(3)))


def transfer_and_measure_edges(
    source_graph: nx.Graph,
    target_graph: nx.Graph,
    selected_atoms: dict[int, dict[str, Any]],
    mapping: dict[int, int],
) -> list[dict[str, Any]]:
    edge_records: list[dict[str, Any]] = []
    for source_left, source_right, source_edge in sorted(
        source_graph.edges(data=True), key=lambda edge: tuple(sorted(edge[:2]))
    ):
        target_left = mapping[source_left]
        target_right = mapping[source_right]
        pair = tuple(sorted((target_left, target_right)))
        if not target_graph.has_edge(*pair):
            raise ValueError(
                f"Mapped source edge {source_left}-{source_right} has no PDB adjacency edge {pair}"
            )
        left_atom = selected_atoms[target_left]
        right_atom = selected_atoms[target_right]
        bond_class = source_edge["bond_class"]
        interval = interval_for_edge(
            source_graph.nodes[source_left]["element"],
            source_graph.nodes[source_right]["element"],
            bond_class,
        )
        distance = distance_angstrom(
            left_atom["coordinates_angstrom"], right_atom["coordinates_angstrom"]
        )
        edge_records.append(
            {
                "source_map_ids": [source_left, source_right],
                "pdb_serials": [target_left, target_right],
                "pdb_atom_names_in_serial_order": [left_atom["name"], right_atom["name"]],
                "coordinates_angstrom_in_serial_order": [
                    left_atom["coordinates_angstrom"],
                    right_atom["coordinates_angstrom"],
                ],
                "element_pair": [
                    source_graph.nodes[source_left]["element"],
                    source_graph.nodes[source_right]["element"],
                ],
                "source_bond_order_class": bond_class,
                "source_bond_order_numeric": source_edge["bond_order"],
                "pdb_adjacency_sources": target_graph.edges[pair]["pdb_adjacency_sources"],
                "distance_angstrom": round(distance, 4),
                "broad_plausible_interval_angstrom": interval,
                "geometry_screen_result": "in_range" if interval[0] <= distance <= interval[1] else "out_of_range",
            }
        )
    return edge_records


def mapping_symmetry_classes(
    embeddings: list[dict[int, int]],
    automorphisms: list[dict[int, int]],
    map_ids: list[int],
    class_id_prefix: str = "source_bond_label_class",
) -> list[dict[str, Any]]:
    groups: dict[tuple[int, ...], list[int]] = defaultdict(list)
    embedding_index: dict[int, dict[int, int]] = {}
    for index, mapping in enumerate(embeddings, start=1):
        embedding_index[index] = mapping
        signatures = [
            tuple(mapping[automorphism[map_id]] for map_id in map_ids)
            for automorphism in automorphisms
        ]
        groups[min(signatures)].append(index)
    result = []
    for group_index, (signature, members) in enumerate(sorted(groups.items()), start=1):
        result.append(
            {
                "symmetry_class_id": f"{class_id_prefix}_{group_index:02d}",
                "embedding_ids": [f"exact_mapping_{index:02d}" for index in members],
                "embedding_count": len(members),
                "representative_embedding_id": f"exact_mapping_{members[0]:02d}",
                "representative_map_id_to_pdb_serial": {
                    str(map_id): embedding_index[members[0]][map_id] for map_id in map_ids
                },
                "class_signature_serials": list(signature),
            }
        )
    return result


def hydrogen_neighbors_for_selected(
    selected: dict[int, dict[str, Any]],
    parsed: dict[str, Any],
) -> dict[int, list[dict[str, Any]]]:
    result: dict[int, list[dict[str, Any]]] = {}
    adjacency: dict[int, set[int]] = defaultdict(set)
    for left, right in parsed["conect_edges"]:
        adjacency[left].add(right)
        adjacency[right].add(left)
    for serial, atom in selected.items():
        neighbors = [
            parsed["atoms"][neighbor]
            for neighbor in adjacency[serial]
            if parsed["atoms"][neighbor]["element"] in {"H", "D"}
        ]
        result[serial] = [
            {
                "serial": neighbor["serial"],
                "name": neighbor["name"],
                "element": neighbor["element"],
                "connectivity_source": "CONECT",
            }
            for neighbor in sorted(neighbors, key=lambda item: (item["name"], item["serial"]))
        ]
    return result


def expected_anchor_consistency(
    embeddings: list[dict[int, int]],
    linker_state: dict[str, Any],
    selected_atoms: dict[int, dict[str, Any]],
) -> dict[str, Any]:
    linker_map = linker_state["candidate_source_to_model_atom_map"]
    checks: dict[str, Any] = {}
    all_ok = True
    for map_id, source_key in ROOF_ANCHOR_MAP_KEYS.items():
        source_record = linker_map[source_key]
        candidate_serials = [
            serial
            for serial, atom in selected_atoms.items()
            if atom["chain"] == source_record["pdb_atom"]["chain"]
            and atom["residue_number"] == str(source_record["pdb_atom"]["residue_number"])
            and atom["insertion_code"] == ""
            and atom["residue"] == source_record["pdb_atom"]["residue"]
            and atom["name"] == source_record["pdb_atom"]["name"]
            and atom["altloc"] == ""
            and atom["element"] == source_record["pdb_atom"]["element"].upper()
        ]
        candidate_serial = candidate_serials[0] if len(candidate_serials) == 1 else None
        observed_serials = sorted({mapping[map_id] for mapping in embeddings})
        consistent = candidate_serial is not None and observed_serials == [candidate_serial]
        all_ok = all_ok and consistent
        checks[str(map_id)] = {
            "candidate_coordinate_serial": candidate_serial,
            "exact_mapping_candidate_serials": observed_serials,
            "unique_and_consistent": consistent,
        }
    return {"all_tail_azide_anchor_maps_consistent": all_ok, "by_map_id": checks}


def state_audit(
    model: dict[str, Any],
    linker_state: dict[str, Any],
    prior_roof_state: dict[str, Any],
    source_graph: nx.Graph,
    normalized_source_graph: nx.Graph,
    normalized_source_automorphisms: list[dict[int, int]],
    normalized_carboxylate_groups: list[dict[str, Any]],
    roof_graph: nx.Graph,
    expected_elements: Counter[str],
) -> dict[str, Any]:
    coordinate_path = project_path(model["output_file"])
    before_hash = check_hash(
        coordinate_path,
        model["output_sha256"],
        f"corrected candidate {model['state']}",
    )
    parsed = parse_pdb(coordinate_path)
    target, selected_atoms, target_metadata = selected_target_graph(parsed, include_nz=True)
    roof_target, roof_selected_atoms, roof_metadata = selected_target_graph(parsed, include_nz=False)

    roof_embeddings, roof_preliminary = roof_only_embeddings(
        roof_graph,
        roof_target,
        roof_selected_atoms,
        linker_state,
    )

    inventory = Counter(target_metadata["selected_heavy_element_counts"])
    target_has_expected_inventory = inventory == expected_elements
    source_to_target_matcher = nx.algorithms.isomorphism.GraphMatcher(
        source_graph,
        target,
        node_match=node_element_match,
    )
    embeddings = [
        dict(mapping)
        for mapping in source_to_target_matcher.isomorphisms_iter()
        if len(mapping) == source_graph.number_of_nodes()
    ]
    embeddings.sort(key=lambda mapping: tuple(mapping[map_id] for map_id in sorted(mapping)))
    target_node_ids = set(target.nodes)
    target_edge_pairs = {tuple(sorted(edge)) for edge in target.edges}
    unmatched_target_nodes: list[dict[str, Any]] = []
    unmatched_target_edges: list[list[int]] = []
    mapping_records: list[dict[str, Any]] = []
    map_candidates: dict[int, set[int]] = defaultdict(set)
    target_to_map_candidates: dict[int, set[int]] = defaultdict(set)

    full_source_automorphisms = graph_automorphisms(source_graph, preserve_bond_labels=True)
    exact_mapping_symmetry_classes = mapping_symmetry_classes(
        embeddings,
        full_source_automorphisms,
        sorted(source_graph.nodes),
    ) if embeddings else []
    normalized_symmetry_classes = mapping_symmetry_classes(
        embeddings,
        normalized_source_automorphisms,
        sorted(normalized_source_graph.nodes),
        class_id_prefix="carboxylate_equivalence_class",
    ) if embeddings else []
    normalized_map_ambiguity = {
        str(map_id): {
            "candidate_count": len({mapping[map_id] for mapping in embeddings}),
            "unique_across_all_exact_mappings": len({mapping[map_id] for mapping in embeddings}) == 1,
            "candidate_pdb_serials": sorted({mapping[map_id] for mapping in embeddings}),
        }
        for map_id in sorted(normalized_source_graph.nodes)
    }
    normalized_comparison = {
        "comparison_type": "Separate graph-equivalence analysis for the declared fully deprotonated tri-carboxylate scenario.",
        "normalized_graph": {
            "heavy_atom_count": normalized_source_graph.number_of_nodes(),
            "adjacency_edge_count": normalized_source_graph.number_of_edges(),
            "normalized_carboxylate_groups": normalized_carboxylate_groups,
            "normalization_rule": "Only the six declared C–O edges are assigned shared carboxylate_resonance labels with numeric bond order 1.5; every other source edge label remains unchanged.",
        },
        "exact_mapping_count": len(embeddings),
        "exact_mappings": [
            {
                "mapping_id": f"exact_mapping_{index:02d}",
                "map_id_to_pdb_serial": {
                    str(map_id): serial for map_id, serial in sorted(mapping.items())
                },
            }
            for index, mapping in enumerate(embeddings, start=1)
        ],
        "normalized_source_automorphism_count": len(normalized_source_automorphisms),
        "symmetry_equivalence_rule": "Two exact coordinate mappings share a normalized class when composition with an element- and normalized-bond-label-preserving automorphism of the copied graph yields the other.",
        "symmetry_classes": normalized_symmetry_classes,
        "distinct_chemically_equivalent_class_count": len(normalized_symmetry_classes),
        "map_id_ambiguity": normalized_map_ambiguity,
        "interpretation_limits": [
            "This graph normalization represents equivalence of the two C–O edges in each assumed carboxylate; it does not assign PDB oxygen names.",
            "The exact coordinate mappings are reused without ranking or selection by coordinates or distances.",
            "The normalized graph is not a pKa or protonation validation and does not establish partial charges or force-field parameters.",
        ],
    }

    if embeddings:
        used_nodes = set(embeddings[0].values())
        unmatched_target_nodes = [
            pdb_atom_record(selected_atoms[serial])
            for serial in sorted(target_node_ids - used_nodes)
        ]
        mapped_edges_first = {
            tuple(sorted((embeddings[0][left], embeddings[0][right])))
            for left, right in source_graph.edges
        }
        unmatched_target_edges = [list(edge) for edge in sorted(target_edge_pairs - mapped_edges_first)]
    else:
        # Without an exact bijection there is no representative assignment.
        # Report every target atom/edge as unassigned to the source graph.
        unmatched_target_nodes = [
            pdb_atom_record(selected_atoms[serial]) for serial in sorted(target_node_ids)
        ]
        unmatched_target_edges = [list(edge) for edge in sorted(target_edge_pairs)]
        expected_degrees = Counter(
            (source_graph.nodes[node]["element"], source_graph.degree[node])
            for node in source_graph.nodes
        )
        observed_degrees = Counter(
            (target.nodes[node]["element"], target.degree[node]) for node in target.nodes
        )
        target_metadata["degree_element_signature_mismatch"] = {
            "source_minus_target": {
                f"{element}:degree_{degree}": count - observed_degrees[(element, degree)]
                for (element, degree), count in expected_degrees.items()
                if count != observed_degrees[(element, degree)]
            },
            "target_minus_source": {
                f"{element}:degree_{degree}": count - expected_degrees[(element, degree)]
                for (element, degree), count in observed_degrees.items()
                if count != expected_degrees[(element, degree)]
            },
        }

    hydrogen_by_serial = hydrogen_neighbors_for_selected(selected_atoms, parsed)
    for index, mapping in enumerate(embeddings, start=1):
        for map_id, serial in mapping.items():
            map_candidates[map_id].add(serial)
            target_to_map_candidates[serial].add(map_id)
        edge_records = transfer_and_measure_edges(source_graph, target, selected_atoms, mapping)
        mapping_records.append(
            {
                "mapping_id": f"exact_mapping_{index:02d}",
                "map_id_to_pdb_serial": {
                    str(map_id): serial for map_id, serial in sorted(mapping.items())
                },
                "source_bond_labels_transferred_after_mapping": True,
                "mapped_heavy_heavy_edge_count": len(edge_records),
                "mapped_heavy_heavy_edges": edge_records,
            }
        )

    map_ambiguity = {
        str(map_id): {
            "candidate_count": len(map_candidates.get(map_id, set())),
            "unique_across_all_exact_mappings": len(map_candidates.get(map_id, set())) == 1,
            "candidate_pdb_serials": sorted(map_candidates.get(map_id, set())),
            "candidate_pdb_atoms": [
                pdb_atom_record(selected_atoms[serial])
                for serial in sorted(map_candidates.get(map_id, set()))
            ],
        }
        for map_id in sorted(source_graph.nodes)
    }
    target_atom_inventory = {
        str(serial): {
            **pdb_atom_record(atom, include_coordinates=True),
            "hydrogen_neighbors_from_conect": hydrogen_by_serial[serial],
            "candidate_source_map_ids_across_exact_mappings": sorted(target_to_map_candidates.get(serial, set())),
        }
        for serial, atom in sorted(selected_atoms.items())
    }
    anchor_consistency = expected_anchor_consistency(embeddings, linker_state, selected_atoms)
    after_hash = sha256(coordinate_path)

    same_name_groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for serial, atom in sorted(selected_atoms.items()):
        same_name_groups[atom["name"]].append(pdb_atom_record(atom))
    duplicate_selected_names = {
        name: records for name, records in sorted(same_name_groups.items()) if len(records) > 1
    }
    unique_maps = {
        map_id: next(iter(serials))
        for map_id, serials in map_candidates.items()
        if len(serials) == 1
    }
    explicit_candidate_aliases = {}
    for map_id, label, expected_name in (
        (79, "roof tail C22", "CD6" if model["state"] == "active_open" else "C"),
        (90, "linker carbonyl oxygen", "OD" if model["state"] == "active_open" else "O"),
    ):
        serial = unique_maps.get(map_id)
        atom = selected_atoms.get(serial) if serial is not None else None
        explicit_candidate_aliases[str(map_id)] = {
            "source_label": label,
            "expected_state_specific_pdb_name": expected_name,
            "observed_target_atom": pdb_atom_record(atom) if atom is not None else None,
            "unique_and_as_expected": atom is not None and atom["name"] == expected_name,
        }

    source_edges = source_graph.number_of_edges()
    exact_edge_counts_match = source_edges == target.number_of_edges()
    exact_node_counts_match = source_graph.number_of_nodes() == target.number_of_nodes()
    full_isomorphism_found = bool(embeddings)
    all_exact_embeddings_cover_all_target_nodes = bool(embeddings) and all(
        set(mapping.values()) == target_node_ids for mapping in embeddings
    )
    maps_78_92_unique = all(
        map_ambiguity[str(map_id)]["unique_across_all_exact_mappings"] for map_id in range(78, 93)
    ) if embeddings else False

    return {
        "state": model["state"],
        "coordinate_candidate_file": model["output_file"],
        "coordinate_candidate_sha256_expected_by_reconciliation": model["output_sha256"],
        "coordinate_candidate_sha256_before_audit": before_hash,
        "coordinate_candidate_sha256_after_audit": after_hash,
        "coordinate_candidate_unchanged_during_audit": before_hash == after_hash,
        "source_pdb_file_before_reconciliation": model["source_file"],
        "source_pdb_sha256_before_reconciliation": model["source_sha256"],
        "selected_target_scope": {
            "heavy_atoms": "all heavy atoms from LYX B30 and MC C1 plus exactly LYS B29 NZ",
            "lys_b29_nz_selection": "chain B, residue 29, insertion code blank, residue LYS, atom NZ; exactly one match required",
            "excluded_atoms": "LYS B29 CE and every other atom outside LYX B30, MC C1, and the single selected NZ are excluded",
            "hydrogen_and_deuterium_excluded_from_graph": True,
            "selected_node_induction": "Only CONECT/LINK edges whose two endpoints are selected nodes enter the target graph.",
        },
        "target_inventory": target_metadata,
        "roof_only_preliminary_result": {
            "source_roof_graph": graph_summary(roof_graph),
            "candidate_anchor_embedding": roof_preliminary,
            "expected_anchored_embeddings": EXPECTED_ROOF_EMBEDDINGS,
            "prior_roof_audit_anchored_embedding_count": prior_roof_state["anchored_embedding_count"],
            "matches_prior_roof_audit_count": roof_preliminary["anchored_embedding_count"]
            == prior_roof_state["anchored_embedding_count"],
            "preliminary_only": True,
            "interpretation": "The anchored roof-only embeddings reproduce the precursor search. They are not complete-product assignments.",
        },
        "resolved_pdb_link_records": parsed["resolved_links"],
        "pdb_adjacency_record_statistics": {
            "raw_directed_conect_reference_count": parsed["raw_conect_reference_count"],
            "unique_directed_conect_reference_count": parsed["unique_directed_conect_reference_count"],
            "directional_or_reciprocal_conect_reference_duplicates_collapsed": parsed["raw_conect_reference_count"]
            - parsed["unique_directed_conect_reference_count"],
            "unique_undirected_conect_edge_count": parsed["conect_only_undirected_edge_count"],
            "conect_link_overlapping_undirected_edge_count_deduplicated": parsed["conect_link_overlap_pair_count"],
            "resolved_link_record_count": len(parsed["resolved_links"]),
            "link_endpoints_require_unique_full_identity_match": True,
            "full_link_identity_fields": [
                "chain",
                "residue_number",
                "insertion_code",
                "residue",
                "altloc",
                "atom_name",
            ],
            "reciprocal_conect_and_conect_link_overlap_deduplicated": True,
            "selected_graph_edges_include_link_only_edges": True,
        },
        "complete_product_exact_graph_comparison": {
            "comparison_type": "Exact element-labeled graph isomorphism; no extra target nodes or edges and no monomorphism.",
            "source_node_count": source_graph.number_of_nodes(),
            "source_edge_count": source_graph.number_of_edges(),
            "source_connected_component_count": nx.number_connected_components(source_graph),
            "source_cycle_rank_E_minus_V_plus_components": source_graph.number_of_edges()
            - source_graph.number_of_nodes()
            + nx.number_connected_components(source_graph),
            "target_node_count": target.number_of_nodes(),
            "target_edge_count": target.number_of_edges(),
            "target_connected_component_count": nx.number_connected_components(target) if target else 0,
            "target_cycle_rank_E_minus_V_plus_components": target.number_of_edges()
            - target.number_of_nodes()
            + (nx.number_connected_components(target) if target else 0),
            "source_and_target_node_counts_equal": exact_node_counts_match,
            "source_and_target_edge_counts_equal": exact_edge_counts_match,
            "source_and_target_component_counts_equal": nx.number_connected_components(source_graph)
            == (nx.number_connected_components(target) if target else 0),
            "source_and_target_cycle_ranks_equal": source_graph.number_of_edges()
            - source_graph.number_of_nodes()
            + nx.number_connected_components(source_graph)
            == target.number_of_edges()
            - target.number_of_nodes()
            + (nx.number_connected_components(target) if target else 0),
            "source_and_target_element_counts_equal": inventory == Counter(element_counts(source_graph)),
            "exact_element_labeled_graph_isomorphism_found": full_isomorphism_found,
            "exact_mapping_count": len(embeddings),
            "expected_feasibility_review_mapping_count": 16,
            "matches_independent_feasibility_review_count": len(embeddings) == 16,
            "every_exact_mapping_uses_all_target_nodes": all_exact_embeddings_cover_all_target_nodes,
            "unmatched_target_atoms": unmatched_target_nodes,
            "unmatched_target_edges_after_representative_mapping": unmatched_target_edges,
            "unmapped_source_map_ids_when_no_exact_isomorphism": sorted(source_graph.nodes)
            if not embeddings
            else [],
            "unmapped_source_edges_when_no_exact_isomorphism": [
                list(sorted(edge)) for edge in sorted(source_graph.edges)
            ]
            if not embeddings
            else [],
            "full_source_bond_label_preserving_automorphism_count": len(full_source_automorphisms),
            "full_source_bond_label_preserving_automorphisms": [
                {str(source): target for source, target in sorted(automorphism.items())}
                for automorphism in full_source_automorphisms
            ],
            "symmetry_equivalence_rule": "Two exact source-to-PDB mappings share a class when composing one with an element- and source bond-label-preserving automorphism of the complete source graph yields the other.",
            "symmetry_classes": exact_mapping_symmetry_classes,
            "distinct_bond_label_transfer_pattern_count": len(exact_mapping_symmetry_classes),
            "multiple_bond_label_transfer_patterns_remain": len(exact_mapping_symmetry_classes) > 1,
            "topology_assignment_status": (
                "exact_element_adjacency_match_but_multiple_bond_label_mapping_classes_remain"
                if full_isomorphism_found and len(exact_mapping_symmetry_classes) > 1
                else "exact_element_adjacency_match_with_one_bond_label_mapping_class"
                if full_isomorphism_found
                else "exact_isomorphism_not_found_mismatch_reported"
            ),
            "unresolved_scaffold_mapping_and_bond_label_transfer_report": {
                "exact_unlabeled_element_adjacency_match": full_isomorphism_found,
                "multiple_bond_label_mapping_classes_remain": len(exact_mapping_symmetry_classes) > 1,
                "distinct_source_bond_label_transfer_patterns": len(exact_mapping_symmetry_classes),
                "unique_topology_or_bond_label_assignment_claimed": False,
                "details": "The source drawing supports the reported 1,4-disubstituted click regiochemistry. The remaining symmetry classes concern scaffold-to-coordinate mapping and source bond-label transfer; they do not represent alternative click regioisomer claims. PDB coordinates and adjacency are not used to choose among these assignments.",
            },
            "map_id_ambiguity": map_ambiguity,
            "map_ids_78_to_92_unique_across_all_exact_mappings": maps_78_92_unique,
            "map_ids_1_to_77_with_multiple_pdb_candidates": [
                map_id
                for map_id in range(1, 78)
                if map_ambiguity[str(map_id)]["candidate_count"] > 1
            ],
            "target_atom_inventory_by_pdb_serial": target_atom_inventory,
            "exact_mappings": mapping_records,
            "interpretation": "Source bond labels are attached to target adjacency edges only after each exact graph mapping has been enumerated. PDB adjacency does not carry bond order, charge, or stereochemistry.",
        },
        "conditional_carboxylate_resonance_normalized_comparison": normalized_comparison,
        "candidate_linker_tail_anchor_consistency": anchor_consistency,
        "candidate_coordinate_aliases_and_name_collisions": {
            "state_specific_map_id_aliases": explicit_candidate_aliases,
            "same_name_collisions_in_selected_graph_resolved_by_full_identity_and_serial": duplicate_selected_names,
            "lyx_c1_and_mc_c1_are_distinguished_by_chain_and_residue": [
                pdb_atom_record(atom)
                for atom in sorted(selected_atoms.values(), key=lambda item: item["serial"])
                if atom["name"] == "C1"
            ],
        },
        "hydrogen_neighbor_connectivity_from_conect": {
            "data_location": "complete_product_exact_graph_comparison.target_atom_inventory_by_pdb_serial[*].hydrogen_neighbors_from_conect",
            "hydrogen_edges_source": "CONECT records only",
            "hydrogen_names_are_model_records_not_validated_chemistry": True,
        },
        "input_pdb_bytes_unchanged": before_hash == after_hash,
    }


def input_hash_record(path: Path, expected: str | None, expected_from: str | None) -> dict[str, Any]:
    observed = sha256(path)
    if expected is not None and observed != expected:
        raise ValueError(f"Hash check failed for {path.relative_to(ROOT)}: {observed} != {expected}")
    return {
        "path": path.relative_to(ROOT).as_posix(),
        "sha256_observed": observed,
        "sha256_expected": expected,
        "expected_hash_source": expected_from,
        "hash_check": "pass" if expected is not None else "recorded_only_no_independent_expected_hash",
    }


def prepare_inputs() -> dict[str, Any]:
    product_map = load_json(PRODUCT_MAP_PATH)
    reconciliation = load_json(RECONCILIATION_PATH)
    prior_roof = load_json(PRIOR_ROOF_MAPPING_PATH)
    linker = load_json(LINKER_MAPPING_PATH)
    roof_graph_audit = load_json(ROOF_GRAPH_AUDIT_PATH)
    roof_map = load_json(ROOF_MAP_PATH)
    source_chemistry = load_json(SOURCE_CHEMISTRY_PATH)

    product_graph_path = project_path(product_map["graph_file"])
    roof_graph_path = project_path(roof_graph_audit["input"]["graph_file"])

    product_graph_hash = check_hash(product_graph_path, product_map["graph_sha256"], "product SMILES")
    roof_graph_hash = check_hash(
        roof_graph_path,
        product_map["input_audit_files"]["roof_graph"]["sha256"],
        "source roof SMILES in product map legend",
    )
    check_hash(ROOF_MAP_PATH, product_map["input_audit_files"]["roof_atom_map_legend"]["sha256"], "roof atom-map legend")
    check_hash(LINKER_MAPPING_PATH, product_map["input_audit_files"]["linker_mapping_audit"]["sha256"], "linker mapping audit")
    check_hash(SOURCE_CHEMISTRY_PATH, product_map["input_audit_files"]["source_chemistry_audit"]["sha256"], "source chemistry audit")

    check_hash(roof_graph_path, roof_graph_audit["input"]["graph_sha256"], "roof graph validation input")
    check_hash(ROOF_MAP_PATH, roof_graph_audit["input"]["atom_map_sha256"], "roof graph validation map legend")
    check_hash(roof_graph_path, prior_roof["inputs"]["mapped_graph_sha256"], "prior roof mapping graph input")
    check_hash(ROOF_MAP_PATH, prior_roof["inputs"]["atom_map_legend_sha256"], "prior roof mapping atom-map legend")
    check_hash(ROOF_GRAPH_AUDIT_PATH, prior_roof["inputs"]["graph_validation_audit_sha256"], "prior roof graph validation audit")
    check_hash(RECONCILIATION_PATH, prior_roof["inputs"]["coordinate_reconciliation_audit_sha256"], "coordinate reconciliation audit")
    check_hash(LINKER_MAPPING_PATH, prior_roof["inputs"]["anchor_mapping_audit_sha256"], "prior roof mapping linker audit")

    product_inputs = product_map["input_audit_files"]
    if product_inputs["roof_graph"]["path"] != roof_graph_audit["input"]["graph_file"]:
        raise ValueError("Product legend and roof validation audit refer to different source roof SMILES")
    if product_inputs["roof_atom_map_legend"]["path"] != roof_graph_audit["input"]["atom_map_file"]:
        raise ValueError("Product legend and roof validation audit refer to different atom-map legends")
    if product_map["graph_file"] != "structures/nnc2215-b29-linker-roof-product.smi":
        raise ValueError(f"Unexpected complete product graph path: {product_map['graph_file']}")

    reconciliation_sha = sha256(RECONCILIATION_PATH)
    roof_mapping_sha = sha256(PRIOR_ROOF_MAPPING_PATH)
    linker_sha = sha256(LINKER_MAPPING_PATH)
    reconciliation_models = {item["state"]: item for item in reconciliation["models"]}
    linker_models = {item["state"]: item for item in linker["states"]}
    prior_roof_models = {item["state"]: item for item in prior_roof["states"]}
    if (
        set(reconciliation_models) != set(STATE_ORDER)
        or set(linker_models) != set(STATE_ORDER)
        or set(prior_roof_models) != set(STATE_ORDER)
    ):
        raise ValueError("Reconciliation/linker/roof records must contain active_open and inactive_closed")

    pdb_hashes: dict[str, dict[str, Any]] = {}
    for state in STATE_ORDER:
        model = reconciliation_models[state]
        linker_model = linker_models[state]
        if model["source_file"] != linker_model["model_file"]:
            raise ValueError(f"Source PDB path mismatch for {state} in linker audit")
        if model["source_sha256"] != linker_model["sha256"]:
            raise ValueError(f"Source PDB hash mismatch for {state} in linker audit")
        source_path = project_path(model["source_file"])
        candidate_path = project_path(model["output_file"])
        source_hash = check_hash(source_path, model["source_sha256"], f"pre-correction PDB {state}")
        candidate_hash = check_hash(candidate_path, model["output_sha256"], f"corrected candidate PDB {state}")
        pdb_hashes[state] = {
            "pre_reconciliation_source_file": model["source_file"],
            "pre_reconciliation_source_sha256": source_hash,
            "corrected_candidate_file": model["output_file"],
            "corrected_candidate_sha256": candidate_hash,
            "corrected_candidate_expected_by_reconciliation": model["output_sha256"],
        }

    product_molecule = parse_molecule(product_graph_path, "complete product graph")
    roof_molecule = parse_molecule(roof_graph_path, "source roof graph")
    product_graph = molecule_graph(product_molecule, "complete product graph")
    roof_graph = molecule_graph(roof_molecule, "source roof graph")
    product_map_ids = sorted(product_graph.nodes)
    roof_map_ids = sorted(roof_graph.nodes)
    product_summary = graph_summary(product_graph)
    roof_summary = graph_summary(roof_graph)
    product_formula = rdMolDescriptors.CalcMolFormula(product_molecule)
    roof_formula = rdMolDescriptors.CalcMolFormula(roof_molecule)

    expected_product_validation = product_map["validation"]
    checks = {
        "product_map_ids_complete_1_to_92": product_map_ids == list(range(1, 93)),
        "product_inventory_matches_legend": Counter(product_summary["heavy_element_counts"]) == EXPECTED_PRODUCT_ELEMENTS,
        "product_node_edge_component_cycle_counts_match_legend": (
            product_summary["heavy_atom_count"] == EXPECTED_PRODUCT_NODES
            and product_summary["adjacency_edge_count"] == EXPECTED_PRODUCT_EDGES
            and product_summary["connected_component_count"] == EXPECTED_PRODUCT_COMPONENTS
            and product_summary["cycle_rank_E_minus_V_plus_components"] == EXPECTED_PRODUCT_CYCLE_RANK
        ),
        "product_graph_matches_committed_validation_counts": (
            expected_product_validation["heavy_atom_count"] == product_summary["heavy_atom_count"]
            and expected_product_validation["adjacency_edge_count"] == product_summary["adjacency_edge_count"]
            and expected_product_validation["cycle_rank_E_minus_V_plus_components"]
            == product_summary["cycle_rank_E_minus_V_plus_components"]
            and Counter(expected_product_validation["heavy_element_inventory"])
            == Counter(product_summary["heavy_element_counts"])
        ),
        "roof_map_ids_complete_1_to_83": roof_map_ids == list(range(1, 84)),
        "roof_inventory_matches_source_audit": Counter(roof_summary["heavy_element_counts"]) == EXPECTED_ROOF_ELEMENTS,
        "roof_counts_match_source_audit": (
            roof_summary["heavy_atom_count"] == EXPECTED_ROOF_NODES
            and roof_summary["adjacency_edge_count"] == EXPECTED_ROOF_EDGES
            and roof_summary["connected_component_count"] == 1
            and roof_summary["cycle_rank_E_minus_V_plus_components"] == 7
        ),
        "roof_legend_matches_product_inherited_map_ids": roof_map_ids
        == product_map["source_transcription"]["roof_azide"]["preserved_roof_map_ids"],
    }
    failed = [name for name, passed in checks.items() if not passed]
    if failed:
        raise ValueError(f"Product/source graph and legend integrity checks failed: {failed}")

    file_hashes = {
        "product_smiles": input_hash_record(product_graph_path, product_map["graph_sha256"], "product map legend graph_sha256"),
        "product_map_legend": input_hash_record(PRODUCT_MAP_PATH, None, None),
        "coordinate_reconciliation_record": input_hash_record(RECONCILIATION_PATH, reconciliation_sha, "prior roof mapping audit input hash"),
        "prior_roof_coordinate_mapping_audit": input_hash_record(PRIOR_ROOF_MAPPING_PATH, None, None),
        "linker_mapping_audit": input_hash_record(LINKER_MAPPING_PATH, product_map["input_audit_files"]["linker_mapping_audit"]["sha256"], "product map legend input hash"),
        "source_roof_smiles": input_hash_record(roof_graph_path, roof_graph_hash, "product map legend and source graph audit"),
        "source_roof_atom_map_legend": input_hash_record(ROOF_MAP_PATH, product_map["input_audit_files"]["roof_atom_map_legend"]["sha256"], "product map legend input hash"),
        "source_roof_graph_validation_audit": input_hash_record(ROOF_GRAPH_AUDIT_PATH, prior_roof["inputs"]["graph_validation_audit_sha256"], "prior roof mapping audit input hash"),
        "source_chemistry_audit": input_hash_record(SOURCE_CHEMISTRY_PATH, product_map["input_audit_files"]["source_chemistry_audit"]["sha256"], "product map legend input hash"),
    }
    for state, record in pdb_hashes.items():
        file_hashes[f"{state}_pre_reconciliation_pdb"] = {
            "path": record["pre_reconciliation_source_file"],
            "sha256_observed": record["pre_reconciliation_source_sha256"],
            "sha256_expected": record["pre_reconciliation_source_sha256"],
            "expected_hash_source": "coordinate reconciliation record and linker mapping audit",
            "hash_check": "pass",
        }
        file_hashes[f"{state}_corrected_candidate_pdb"] = {
            "path": record["corrected_candidate_file"],
            "sha256_observed": record["corrected_candidate_sha256"],
            "sha256_expected": record["corrected_candidate_expected_by_reconciliation"],
            "expected_hash_source": "coordinate reconciliation record output_sha256",
            "hash_check": "pass",
        }

    return {
        "product_map": product_map,
        "reconciliation": reconciliation,
        "prior_roof_mapping": prior_roof,
        "linker_mapping": linker,
        "roof_graph_audit": roof_graph_audit,
        "roof_map": roof_map,
        "source_chemistry": source_chemistry,
        "product_graph_path": product_graph_path,
        "roof_graph_path": roof_graph_path,
        "product_molecule": product_molecule,
        "roof_molecule": roof_molecule,
        "product_graph": product_graph,
        "roof_graph": roof_graph,
        "product_formula": product_formula,
        "roof_formula": roof_formula,
        "product_summary": product_summary,
        "roof_summary": roof_summary,
        "product_map_ids": product_map_ids,
        "roof_map_ids": roof_map_ids,
        "integrity_checks": checks,
        "input_hashes": file_hashes,
        "pdb_hashes": pdb_hashes,
        "reconciliation_models": reconciliation_models,
        "linker_models": linker_models,
        "prior_roof_models": prior_roof_models,
        "prior_roof_sha256": roof_mapping_sha,
        "linker_sha256": linker_sha,
    }


def hydrogen_comparison(states: list[dict[str, Any]]) -> dict[str, Any]:
    by_state = {state["state"]: state for state in states}
    comparison: dict[str, Any] = {}
    for map_id, expected in HYDROGEN_OBSERVATIONS.items():
        per_state: dict[str, list[dict[str, Any]]] = {}
        names: dict[str, list[str]] = {}
        unique: dict[str, bool] = {}
        for state_name in STATE_ORDER:
            complete = by_state[state_name]["complete_product_exact_graph_comparison"]
            candidates = complete["map_id_ambiguity"][str(map_id)]["candidate_pdb_serials"]
            target_inventory = complete["target_atom_inventory_by_pdb_serial"]
            unique[state_name] = len(candidates) == 1
            neighbor_records: list[dict[str, Any]] = []
            neighbor_names: set[str] = set()
            for serial in candidates:
                atom_info = target_inventory[str(serial)]
                for neighbor in atom_info["hydrogen_neighbors_from_conect"]:
                    neighbor_names.add(neighbor["name"])
                    neighbor_records.append(
                        {
                            "parent_pdb_serial": serial,
                            "parent_pdb_atom": {
                                "chain": atom_info["chain"],
                                "residue_number": atom_info["residue_number"],
                                "residue": atom_info["residue"],
                                "altloc": atom_info["altloc"],
                                "name": atom_info["name"],
                            },
                            "hydrogen_neighbor": neighbor,
                        }
                    )
            per_state[state_name] = sorted(
                neighbor_records,
                key=lambda item: (
                    item["parent_pdb_serial"],
                    item["hydrogen_neighbor"]["name"],
                    item["hydrogen_neighbor"]["serial"],
                ),
            )
            names[state_name] = sorted(neighbor_names)
        observed_relation = "same" if names["active_open"] == names["inactive_closed"] else "different"
        comparison[str(map_id)] = {
            "source_map_id": map_id,
            "source_label": expected["label"],
            "candidate_map_is_unique_in_each_state": unique,
            "hydrogen_neighbors_by_state": per_state,
            "hydrogen_neighbor_names_by_state": names,
            "observed_relation": observed_relation,
            "expected_relation": expected["relation"],
            "observed_open_names_match_expected": names["active_open"] == expected["open"],
            "observed_closed_names_match_expected": names["inactive_closed"] == expected["closed"],
            "observed_as_expected": (
                names["active_open"] == expected["open"]
                and names["inactive_closed"] == expected["closed"]
                and observed_relation == expected["relation"]
            ),
            "hydrogen_labels_and_counts_are_model_records_not_validated_chemistry": True,
        }
    return {
        "source": "PDB CONECT records only; hydrogen names/counts are not treated as chemically validated.",
        "by_map_id": comparison,
        "all_required_observations_surface_as_expected": all(
            item["observed_as_expected"] for item in comparison.values()
        ),
        "known_differences": [
            {
                "map_id": map_id,
                "source_label": HYDROGEN_OBSERVATIONS[map_id]["label"],
                "open_names": comparison[str(map_id)]["hydrogen_neighbor_names_by_state"]["active_open"],
                "closed_names": comparison[str(map_id)]["hydrogen_neighbor_names_by_state"]["inactive_closed"],
                "open_count": len(comparison[str(map_id)]["hydrogen_neighbor_names_by_state"]["active_open"]),
                "closed_count": len(comparison[str(map_id)]["hydrogen_neighbor_names_by_state"]["inactive_closed"]),
                "open_minus_closed_names": sorted(
                    set(comparison[str(map_id)]["hydrogen_neighbor_names_by_state"]["active_open"])
                    - set(comparison[str(map_id)]["hydrogen_neighbor_names_by_state"]["inactive_closed"])
                ),
                "closed_minus_open_names": sorted(
                    set(comparison[str(map_id)]["hydrogen_neighbor_names_by_state"]["inactive_closed"])
                    - set(comparison[str(map_id)]["hydrogen_neighbor_names_by_state"]["active_open"])
                ),
                "open_count_minus_closed_count": len(
                    comparison[str(map_id)]["hydrogen_neighbor_names_by_state"]["active_open"]
                ) - len(comparison[str(map_id)]["hydrogen_neighbor_names_by_state"]["inactive_closed"]),
                "closed_has_one_fewer_recorded_hydrogen": len(
                    comparison[str(map_id)]["hydrogen_neighbor_names_by_state"]["inactive_closed"]
                )
                < len(comparison[str(map_id)]["hydrogen_neighbor_names_by_state"]["active_open"]),
            }
            for map_id in (79, 80, 85, 86, 87)
        ],
    }


def require_carboxylate_audit_invariants(
    product_formula: str,
    states: list[dict[str, Any]],
    candidate_pdb_hashes_unchanged: bool,
) -> None:
    if product_formula != EXPECTED_PRODUCT_FORMULA:
        raise ValueError(
            "Neutral NZ-NH2-capped fragment formula must be "
            f"{EXPECTED_PRODUCT_FORMULA}; observed {product_formula}"
        )
    if ASSUMED_CARBOXYLATE_PROTON_REMOVALS != 3:
        raise ValueError(
            "The declared scenario requires exactly three assumed proton removals; "
            f"observed {ASSUMED_CARBOXYLATE_PROTON_REMOVALS}"
        )
    if ASSUMED_CARBOXYLATE_LOCAL_CHARGE_DELTA != -3:
        raise ValueError(
            "The declared scenario requires a local charge delta of -3; "
            f"observed {ASSUMED_CARBOXYLATE_LOCAL_CHARGE_DELTA}"
        )
    derived_deprotonated_formula = subtract_formula_hydrogens(
        product_formula,
        ASSUMED_CARBOXYLATE_PROTON_REMOVALS,
    )
    if EXPECTED_DEPROTONATED_FRAGMENT_FORMULA != derived_deprotonated_formula:
        raise ValueError(
            "The declared proton removals do not match the expected deprotonated formula: "
            f"derived {derived_deprotonated_formula} from {product_formula} minus "
            f"{ASSUMED_CARBOXYLATE_PROTON_REMOVALS} hydrogens, configured "
            f"{EXPECTED_DEPROTONATED_FRAGMENT_FORMULA}"
        )
    if not candidate_pdb_hashes_unchanged:
        raise ValueError("Corrected candidate PDB hashes changed during the audit")

    for state in states:
        neutral = state["complete_product_exact_graph_comparison"]
        normalized = state["conditional_carboxylate_resonance_normalized_comparison"]
        if neutral["exact_mapping_count"] != 16:
            raise ValueError(
                f"{state['state']} neutral exact mapping count must remain 16; "
                f"observed {neutral['exact_mapping_count']}"
            )
        if neutral["distinct_bond_label_transfer_pattern_count"] != 8:
            raise ValueError(
                f"{state['state']} neutral bond-label class count must remain 8; "
                f"observed {neutral['distinct_bond_label_transfer_pattern_count']}"
            )
        if normalized["exact_mapping_count"] != 16:
            raise ValueError(
                f"{state['state']} normalized exact mapping count must be 16; "
                f"observed {normalized['exact_mapping_count']}"
            )
        if normalized["distinct_chemically_equivalent_class_count"] != 1:
            raise ValueError(
                f"{state['state']} carboxylate-normalized class count must be 1; "
                f"observed {normalized['distinct_chemically_equivalent_class_count']} "
                f"for {normalized['exact_mapping_count']} exact mappings"
            )


def main() -> None:
    inputs = prepare_inputs()
    product_graph: nx.Graph = inputs["product_graph"]
    roof_graph: nx.Graph = inputs["roof_graph"]
    if inputs["product_formula"] != EXPECTED_PRODUCT_FORMULA:
        raise ValueError(
            "Neutral NZ-NH2-capped fragment formula must be "
            f"{EXPECTED_PRODUCT_FORMULA}; observed {inputs['product_formula']}"
        )
    normalized_product_graph, normalized_carboxylate_groups = carboxylate_resonance_graph(product_graph)
    normalized_product_autos = graph_automorphisms(
        normalized_product_graph,
        preserve_bond_labels=True,
    )
    product_autos = graph_automorphisms(product_graph, preserve_bond_labels=True)
    roof_autos = graph_automorphisms(roof_graph, preserve_bond_labels=True)
    if len(roof_autos) != 2:
        raise ValueError(f"Expected two bond-label-preserving roof graph automorphisms; found {len(roof_autos)}")

    state_results: list[dict[str, Any]] = []
    for state_name in STATE_ORDER:
        state_results.append(
            state_audit(
                inputs["reconciliation_models"][state_name],
                inputs["linker_models"][state_name],
                inputs["prior_roof_models"][state_name],
                product_graph,
                normalized_product_graph,
                normalized_product_autos,
                normalized_carboxylate_groups,
                roof_graph,
                EXPECTED_PRODUCT_ELEMENTS,
            )
        )
    hydrogen_results = hydrogen_comparison(state_results)

    candidate_hashes_post = {
        model["output_file"]: sha256(project_path(model["output_file"]))
        for model in inputs["reconciliation_models"].values()
    }
    pdb_unchanged_all = all(
        candidate_hashes_post[record["corrected_candidate_file"]]
        == record["corrected_candidate_sha256"]
        for record in inputs["pdb_hashes"].values()
    ) and all(state["coordinate_candidate_unchanged_during_audit"] for state in state_results)
    require_carboxylate_audit_invariants(
        inputs["product_formula"],
        state_results,
        pdb_unchanged_all,
    )

    state_isomorphisms = [
        state["complete_product_exact_graph_comparison"]["exact_element_labeled_graph_isomorphism_found"]
        for state in state_results
    ]
    map_ids_78_92_unique = all(
        state["complete_product_exact_graph_comparison"]["map_ids_78_to_92_unique_across_all_exact_mappings"]
        for state in state_results
    )
    expected_exact_count = all(
        state["complete_product_exact_graph_comparison"]["exact_mapping_count"] == 16
        for state in state_results
    )
    all_exact = all(state_isomorphisms)
    all_edges_measured = all(
        state["complete_product_exact_graph_comparison"]["exact_mapping_count"] > 0
        and all(
            len(mapping["mapped_heavy_heavy_edges"]) == EXPECTED_PRODUCT_EDGES
            for mapping in state["complete_product_exact_graph_comparison"]["exact_mappings"]
        )
        for state in state_results
    )

    product_map = inputs["product_map"]
    result = {
        "schema": "dia-health-one.nnc2215-b29-topology-audit.v1",
        "run": {
            "generated_at_utc": datetime.now(timezone.utc).isoformat(),
            "script": "research/audit_nnc2215_b29_topology.py",
            "script_sha256": sha256(Path(__file__)),
            "python_version": platform.python_version(),
            "rdkit_version": rdBase.rdkitVersion,
            "networkx_version": nx.__version__,
            "gpu_used": False,
            "coordinate_files_modified": False,
        },
        "inputs": {
            "sha256_file_integrity_records": inputs["input_hashes"],
            "hash_validation_boundary": [
                "The product map legend and prior roof coordinate-mapping audit have recorded SHA-256 values but no independent expected self-checksum in the committed input set; their bytes are recorded, and their graph, legend, reconciliation, linker, and source-audit references are cross-checked against the corresponding files.",
                "All artifacts with an upstream expected checksum are verified against that checksum before graph matching; both pre-reconciliation source PDBs and corrected candidate PDBs are checked against reconciliation/linker records.",
            ],
            "input_integrity_checks": inputs["integrity_checks"],
            "corrected_candidate_pdbs": inputs["pdb_hashes"],
            "all_corrected_candidate_pdb_hashes_unchanged_after_audit": pdb_unchanged_all,
            "candidate_hashes_after_audit": candidate_hashes_post,
        },
        "ambiguity_and_mapping_assumptions": [
            "The six map IDs 78-83 used in the roof-only precursor search are fixed to candidate tail/azide anchors from the prior linker mapping audit; this does not authorize those anchors as source chemistry evidence.",
            "Complete-product exact mappings are enumerated without fixing the prior candidate labels. The graph matcher compares element labels and simple adjacency only.",
            "PDB CONECT resolves by atom serial. PDB LINK resolves only by a unique full chain, residue, insertion code, altloc, and atom-name identity.",
            "Source bond-label-preserving automorphisms define mapping equivalence classes. Distinct classes remain separate because they can transfer different source bond labels onto the unlabeled PDB adjacency.",
            "PDB coordinates do not choose among exact mappings or bond-label transfer classes; geometry is reported after mapping as a screen only.",
            "PDB atom names, including the open/closed aliases for map 79 and map 90, are candidate coordinate labels and are not source-authoritative atom names.",
        ],
        "source_product_graph": {
            "graph_file": product_map["graph_file"],
            "graph_sha256": sha256(inputs["product_graph_path"]),
            "map_legend_file": PRODUCT_MAP_PATH.relative_to(ROOT).as_posix(),
            "map_legend_sha256": sha256(PRODUCT_MAP_PATH),
            "source_product_drawing_regiochemistry": product_map["source_transcription"]["click_product"]["regiochemistry_status"],
            "source_assignment_basis": product_map["source_transcription"]["click_product"]["evidence"],
            "product_map_id_legend": product_map["map_id_legend"],
            "formula_with_nz_nh2_cap": inputs["product_formula"],
            "formula_is_actual_conjugate_formula": product_map["validation"]["computed_formula_is_actual_conjugate_formula"],
            **inputs["product_summary"],
            "graph_integrity_checks": inputs["integrity_checks"],
            "bond_order_representation": product_map["bond_order_representation"],
            "representation_boundaries": product_map["representation_boundaries"],
        },
        "conditional_carboxylate_state_scenario": {
            "declared_pH": 7.4,
            "protonation_assumption": "All three declared macrocycle carboxyl groups are fully deprotonated in this pH 7.4 modeling scenario; site-specific pKa evidence is not available in the source record.",
            "carboxylate_map_groups": normalized_carboxylate_groups,
            "neutral_reference": {
                "formula_with_lys_b29_nz_nh2_boundary_cap": EXPECTED_PRODUCT_FORMULA,
                "formal_charge_delta_reference": "neutral source fragment before the assumed three deprotonations",
            },
            "conditional_fragment_delta": {
                "assumed_removed_protons": ASSUMED_CARBOXYLATE_PROTON_REMOVALS,
                "assumed_local_charge_delta": ASSUMED_CARBOXYLATE_LOCAL_CHARGE_DELTA,
                "formula_with_lys_b29_nz_nh2_boundary_cap": EXPECTED_DEPROTONATED_FRAGMENT_FORMULA,
                "formula_derivation": "C62H72N16O14 minus three assumed acid protons gives C62H69N16O14.",
                "charge_scope": "The -3 value applies to the isolated NZ-NH2-capped B29 fragment; it is not the full conjugate charge.",
            },
            "graph_equivalence_rule": "For mapping equivalence only, the two C–O edges in each declared carboxylate group share bond class carboxylate_resonance and numeric order 1.5. All other source bond labels are preserved.",
            "state_class_counts": {
                state["state"]: state["conditional_carboxylate_resonance_normalized_comparison"]["distinct_chemically_equivalent_class_count"]
                for state in state_results
            },
            "assumptions_and_limits": [
                "pH 7.4 and full deprotonation are modeling assumptions, not measured site-specific pKa evidence.",
                "The normalization is a graph-equivalence representation, not a validated protonation calculation.",
                "No PDB oxygen atom names are assigned, and coordinates or distances do not select among mappings.",
                "No partial charges or force-field parameters are created or validated; complete-conjugate protonation and charge remain open.",
            ],
        },
        "source_roof_graph_precursor": {
            "graph_file": inputs["roof_graph_path"].relative_to(ROOT).as_posix(),
            "graph_sha256": sha256(inputs["roof_graph_path"]),
            "map_legend_file": ROOF_MAP_PATH.relative_to(ROOT).as_posix(),
            "map_legend_sha256": sha256(ROOF_MAP_PATH),
            "formula": inputs["roof_formula"],
            **inputs["roof_summary"],
            "bond_label_preserving_automorphism_count": len(roof_autos),
            "bond_label_preserving_automorphisms": [
                {str(source): target for source, target in sorted(automorphism.items())}
                for automorphism in roof_autos
            ],
            "interpretation": "The two source roof automorphisms are a precursor symmetry result; the complete clicked product is audited separately.",
        },
        "geometry_screen_method": {
            "purpose": "Conservative broad plausibility screen for mapped heavy-heavy coordinate distances.",
            "basis": {
                "citation": "B. Cordero et al., Covalent radii revisited, Dalton Trans. 2008, 2832-2838, DOI: 10.1039/B801115J.",
                "doi_url": "https://doi.org/10.1039/B801115J",
                "reported_covalent_radii_angstrom_used": COVALENT_RADII_ANGSTROM,
                "values_are_rounded_from_reported_reference_radii": True,
            },
            "heuristic_bond_class_offsets_from_radii_sum_angstrom": BOND_CLASS_OFFSETS_ANGSTROM,
            "interval_definition": "For each source bond class, add the tabulated low and high heuristic offsets to the sum of the two element covalent radii.",
            "custom_offsets_are_heuristic_not_supplied_by_cordero_paper": True,
            "interpretation_limits": [
                "These numeric intervals are a conservative geometry screen, not bond-order evidence.",
                "Distances are not used to select, rank, or reject graph mappings.",
                "All mapped edges are retained, including those outside the screen.",
                "The PDB adjacency records do not encode bond order, charge, or stereochemistry.",
            ],
        },
        "states": state_results,
        "open_closed_hydrogen_neighbor_comparison": hydrogen_results,
        "assessment": {
            "both_states_have_exact_full_element_labeled_adjacency_isomorphism": all_exact,
            "both_states_match_92_nodes_99_edges_and_inventory": all(
                state["complete_product_exact_graph_comparison"]["source_and_target_node_counts_equal"]
                and state["complete_product_exact_graph_comparison"]["source_and_target_edge_counts_equal"]
                and state["complete_product_exact_graph_comparison"]["source_and_target_element_counts_equal"]
                for state in state_results
            ),
            "exact_mappings_per_state_match_feasibility_review_16": expected_exact_count,
            "map_ids_78_to_92_unique_across_exact_mappings_in_both_states": map_ids_78_92_unique,
            "multiple_bond_label_transfer_classes_remain_in_both_states": all(
                state["complete_product_exact_graph_comparison"]["multiple_bond_label_transfer_patterns_remain"]
                for state in state_results
            ),
            "bond_label_transfer_class_count_by_state": {
                state["state"]: state["complete_product_exact_graph_comparison"]["distinct_bond_label_transfer_pattern_count"]
                for state in state_results
            },
            "carboxylate_normalized_class_count_by_state": {
                state["state"]: state["conditional_carboxylate_resonance_normalized_comparison"]["distinct_chemically_equivalent_class_count"]
                for state in state_results
            },
            "all_states_have_one_carboxylate_normalized_class": all(
                state["conditional_carboxylate_resonance_normalized_comparison"]["distinct_chemically_equivalent_class_count"] == 1
                for state in state_results
            ),
            "state_specific_tail_and_carbonyl_aliases_resolved_as_expected": all(
                all(alias["unique_and_as_expected"] for alias in state["candidate_coordinate_aliases_and_name_collisions"]["state_specific_map_id_aliases"].values())
                for state in state_results
            ),
            "all_mapped_source_bonds_have_coordinate_screen_records": all_edges_measured,
            "required_hydrogen_neighbor_observations_match": hydrogen_results["all_required_observations_surface_as_expected"],
            "corrected_candidate_pdbs_unchanged": pdb_unchanged_all,
            "complete_source_to_pdb_atom_map_unique": all(
                state["complete_product_exact_graph_comparison"]["exact_mapping_count"] == 1
                for state in state_results
            ),
            "unique_topology_claimed": False,
            "interpretation": (
                "The complete 92-node, 99-edge source product graph has an exact element-labeled adjacency isomorphism to each selected corrected PDB target. "
                "Each state has 16 exact map-to-coordinate assignments; map IDs 78-92 are consistent and unique across them, while roof/scaffold map IDs retain ambiguity. "
                "The two bond-label-preserving source roof automorphisms are reported as a precursor result. Source bond labels are transferred only after mapping and the PDB records remain adjacency-only. "
                "Coordinates are screened against disclosed broad heuristic intervals and do not validate bond order or mapping. The source drawing supports 1,4-disubstituted click regiochemistry. In the neutral graph, the 8 bond-label symmetry classes per state arise from the three independent carboxyl carbonyl-O/OH label swaps; they are not alternative click regioisomer claims. Under the separate conditional tri-carboxylate graph, those paired C–O edges are resonance-equivalent and yield one chemistry-equivalent label class per state. The 16 exact coordinate mappings and source-to-PDB atom-label ambiguity remain unresolved."
            ),
            "mismatch_handling": "If any exact isomorphism is absent, or unmatched target atoms/edges remain, the corresponding state reports the unmatched inventory and cannot support a unique topology assignment. Here, multiple source bond-label mapping classes are reported as unresolved scaffold-to-coordinate mapping and bond-label transfer; they do not challenge the source-drawing 1,4 click regiochemistry.",
            "bond_order_or_regiochemistry_validated_from_pdb_coordinates": False,
            "pdb_hydrogen_records_validated_as_chemistry": False,
            "simulation_readiness_established": False,
            "limitations": [
                "Corrected candidate PDBs are model-coordinate derivatives, not experimentally resolved NNC2215 structures.",
                "PDB CONECT and LINK records are adjacency-only; bond orders, formal charges, stereochemistry, and protonation are not validated here.",
                "The source product graph represents the source-supported click product connectivity and does not establish simulation-pH protonation or force-field parameters.",
                "The geometry intervals are intentionally broad heuristic screens and do not establish chemical bonds, mapping, or click regiochemistry.",
                "Hydrogen neighbor labels and counts reproduce model CONECT records and are not validated chemistry.",
            ],
        },
    }

    # Take the final no-write check immediately before replacing only the JSON
    # audit output. No coordinate or simulation input file is written.
    final_pdb_hashes = {
        model["output_file"]: sha256(project_path(model["output_file"]))
        for model in inputs["reconciliation_models"].values()
    }
    if any(
        final_pdb_hashes[model["output_file"]] != model["output_sha256"]
        for model in inputs["reconciliation_models"].values()
    ):
        raise ValueError("A corrected candidate PDB changed during the audit")
    result["inputs"]["candidate_hashes_after_audit"] = final_pdb_hashes
    result["inputs"]["all_corrected_candidate_pdb_hashes_unchanged_after_audit"] = True
    result["run"]["coordinate_files_modified"] = False

    OUTPUT_PATH.write_text(json.dumps(result, indent=2, sort_keys=False) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "output": OUTPUT_PATH.relative_to(ROOT).as_posix(),
                "sha256": sha256(OUTPUT_PATH),
                "states": {
                    state["state"]: {
                        "exact_mappings": state["complete_product_exact_graph_comparison"]["exact_mapping_count"],
                        "neutral_bond_label_classes": state["complete_product_exact_graph_comparison"]["distinct_bond_label_transfer_pattern_count"],
                        "carboxylate_normalized_classes": state["conditional_carboxylate_resonance_normalized_comparison"]["distinct_chemically_equivalent_class_count"],
                        "roof_embeddings": state["roof_only_preliminary_result"]["candidate_anchor_embedding"]["anchored_embedding_count"],
                        "nodes": state["complete_product_exact_graph_comparison"]["target_node_count"],
                        "edges": state["complete_product_exact_graph_comparison"]["target_edge_count"],
                    }
                    for state in state_results
                },
                "assessment": result["assessment"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
