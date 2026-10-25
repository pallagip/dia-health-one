#!/usr/bin/env python3
"""Audit candidate embeddings of the source roof graph in NNC2215 PDBs.

The comparison is element-labeled adjacency only. PDB CONECT/LINK records do
not establish bond order, protonation, stereochemistry, or force-field chemistry.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
from typing import Any

import networkx as nx
from rdkit import Chem, rdBase
from rdkit.Chem import rdMolDescriptors


ROOT = Path(__file__).resolve().parent
GRAPH_PATH = ROOT / "structures" / "nnc2215-roof-propyl-azide.smi"
GRAPH_MAP_PATH = ROOT / "data" / "nnc2215-roof-macrocycle-atom-map.json"
GRAPH_AUDIT_PATH = ROOT / "data" / "nnc2215-roof-macrocycle-graph-audit.json"
COORDINATE_AUDIT_PATH = ROOT / "data" / "nnc2215-carbon-reconciliation-candidate.json"
ANCHOR_AUDIT_PATH = ROOT / "data" / "nnc2215-linker-mapping-audit.json"
OUTPUT_PATH = ROOT / "data" / "nnc2215-roof-graph-coordinate-mapping-audit.json"

EXPECTED_FORMULA = "C56H63N15O12"
EXPECTED_QUERY_HEAVY_ATOMS = 83
EXPECTED_QUERY_EDGES = 89
EXPECTED_TARGET_ELEMENTS = Counter({"C": 62, "N": 15, "O": 14})
EXPECTED_TARGET_HEAVY_ATOMS = 91
EXPECTED_TARGET_EDGES = 98
EXPECTED_MATCHES_PER_STATE = 16

ANCHOR_MAP_KEYS = {
    78: "source_macrocycle_c21_candidate",
    79: "source_macrocycle_c22_candidate",
    80: "source_macrocycle_c23_candidate",
    81: "source_azide_n_attached_to_c23_candidate",
    82: "source_triazole_n2_candidate",
    83: "source_triazole_n1_candidate",
}

EXPECTED_UNMATCHED_NAMES = {
    "CD1",
    "CD2",
    "C1",
    "C2",
    "CD3",
    "CD4",
    "O1",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def pdb_address(atom: dict[str, Any]) -> dict[str, str]:
    return {
        "chain": atom["chain"],
        "residue_number": atom["residue_number"],
        "insertion_code": atom["insertion_code"],
        "residue": atom["residue"],
        "name": atom["name"],
        "element": atom["element"],
    }


def address_key(address: dict[str, str]) -> tuple[str, ...]:
    return tuple(
        address[key]
        for key in (
            "chain",
            "residue_number",
            "insertion_code",
            "residue",
            "name",
            "element",
        )
    )


def parse_link_endpoint(line: str, first: bool) -> dict[str, str]:
    if first:
        fields = (line[12:16], line[17:20], line[21:22], line[22:26], line[26:27])
    else:
        fields = (line[42:46], line[47:50], line[51:52], line[52:56], line[56:57])
    atom_name, residue, chain, residue_number, insertion_code = (
        value.strip() for value in fields
    )
    return {
        "name": atom_name,
        "residue": residue,
        "chain": chain,
        "residue_number": residue_number,
        "insertion_code": insertion_code,
    }


def parse_pdb(path: Path) -> tuple[dict[int, dict[str, Any]], set[tuple[int, int]]]:
    atoms: dict[int, dict[str, Any]] = {}
    raw_conect: list[tuple[int, int]] = []
    link_records: list[tuple[dict[str, str], dict[str, str]]] = []

    for line in path.read_text(encoding="ascii").splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            serial = int(line[6:11])
            if serial in atoms:
                raise ValueError(f"Duplicate PDB atom serial {serial} in {path}")
            atom = {
                "serial": serial,
                "name": line[12:16].strip(),
                "altloc": line[16:17].strip(),
                "residue": line[17:20].strip(),
                "chain": line[21:22].strip(),
                "residue_number": line[22:26].strip(),
                "insertion_code": line[26:27].strip(),
                "element": line[76:78].strip().upper(),
            }
            if not atom["element"]:
                raise ValueError(f"Missing element for atom serial {serial} in {path}")
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
            link_records.append((parse_link_endpoint(line, True), parse_link_endpoint(line, False)))

    edges: set[tuple[int, int]] = set()
    for left, right in raw_conect:
        if left not in atoms or right not in atoms:
            raise ValueError(f"CONECT references absent atom serial {left}-{right} in {path}")
        if left != right:
            edges.add(tuple(sorted((left, right))))

    for left_address, right_address in link_records:
        resolved: list[int] = []
        for endpoint in (left_address, right_address):
            matches = [
                serial
                for serial, atom in atoms.items()
                if all(atom[key] == value for key, value in endpoint.items())
            ]
            if len(matches) != 1:
                raise ValueError(
                    f"LINK endpoint {endpoint} resolves to {matches} in {path}"
                )
            resolved.append(matches[0])
        if resolved[0] != resolved[1]:
            edges.add(tuple(sorted((resolved[0], resolved[1]))))

    return atoms, edges


def build_source_graph(mol: Chem.Mol) -> nx.Graph:
    graph = nx.Graph()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1:
            continue
        map_id = atom.GetAtomMapNum()
        if map_id <= 0:
            raise ValueError("Every source graph heavy atom must have a positive map ID")
        graph.add_node(map_id, map_id=map_id, element=atom.GetSymbol())
    for bond in mol.GetBonds():
        left = mol.GetAtomWithIdx(bond.GetBeginAtomIdx())
        right = mol.GetAtomWithIdx(bond.GetEndAtomIdx())
        if left.GetAtomicNum() != 1 and right.GetAtomicNum() != 1:
            graph.add_edge(left.GetAtomMapNum(), right.GetAtomMapNum())
    return graph


def selected_target_graph(
    path: Path,
) -> tuple[nx.Graph, dict[int, dict[str, Any]], dict[str, int]]:
    atoms, all_edges = parse_pdb(path)
    selected = {
        serial: atom
        for serial, atom in atoms.items()
        if (atom["chain"], atom["residue_number"], atom["residue"])
        in {("B", "30", "LYX"), ("C", "1", "MC")}
        and atom["element"] not in {"H", "D"}
    }
    target = nx.Graph()
    for serial, atom in selected.items():
        target.add_node(serial, serial=serial, element=atom["element"], atom=atom)
    for left, right in all_edges:
        if left in selected and right in selected:
            target.add_edge(left, right)

    counts = Counter(atom["element"] for atom in selected.values())
    return target, selected, dict(sorted(counts.items()))


def find_anchor_serial(
    selected_atoms: dict[int, dict[str, Any]],
    source_atom: dict[str, Any],
    map_id: int,
    expected_element: str,
) -> int:
    source_record = source_atom["pdb_atom"]
    expected = {
        "chain": source_record["chain"],
        "residue_number": str(source_record["residue_number"]),
        "insertion_code": "",
        "residue": source_record["residue"],
        "name": source_record["name"],
        "element": source_record["element"].upper(),
    }
    if expected["element"] != expected_element:
        raise ValueError(
            f"Anchor map {map_id} has element {expected['element']}, expected {expected_element}"
        )
    matches = [
        serial
        for serial, atom in selected_atoms.items()
        if address_key(pdb_address(atom)) == address_key(expected)
    ]
    if len(matches) != 1:
        raise ValueError(f"Anchor map {map_id} address {expected} resolves to {matches}")
    return matches[0]


def atom_record(atom: dict[str, Any]) -> dict[str, Any]:
    return {"serial": atom["serial"], **pdb_address(atom)}


def audit_state(
    model: dict[str, Any],
    anchor_state: dict[str, Any],
    source_graph: nx.Graph,
    expected_coordinate_record: dict[str, Any],
) -> dict[str, Any]:
    relative_path = model["output_file"]
    pdb_path = ROOT / relative_path
    observed_hash = sha256(pdb_path)
    if observed_hash != model["output_sha256"]:
        raise ValueError(f"Coordinate-reconciliation checksum mismatch for {relative_path}")

    target, selected_atoms, element_counts = selected_target_graph(pdb_path)
    observed_counts = Counter(element_counts)
    if len(target) != EXPECTED_TARGET_HEAVY_ATOMS:
        raise ValueError(f"Expected 91 selected heavy atoms in {relative_path}; found {len(target)}")
    if observed_counts != EXPECTED_TARGET_ELEMENTS:
        raise ValueError(
            f"Unexpected selected heavy-element inventory in {relative_path}: {dict(observed_counts)}"
        )
    if target.number_of_edges() != EXPECTED_TARGET_EDGES:
        raise ValueError(
            f"Expected 98 selected adjacency edges in {relative_path}; found {target.number_of_edges()}"
        )

    map_to_expected_atom = {
        atom.GetAtomMapNum(): atom.GetSymbol()
        for atom in Chem.MolFromSmiles(GRAPH_PATH.read_text(encoding="utf-8").split()[0]).GetAtoms()
        if atom.GetAtomicNum() != 1
    }
    anchors: dict[int, int] = {}
    anchor_details: dict[str, Any] = {}
    source_map = anchor_state["candidate_source_to_model_atom_map"]
    for map_id, source_key in ANCHOR_MAP_KEYS.items():
        source_record = source_map[source_key]
        serial = find_anchor_serial(
            selected_atoms,
            source_record,
            map_id,
            map_to_expected_atom[map_id],
        )
        anchors[map_id] = serial
        anchor_details[str(map_id)] = {
            "source_graph_label": {
                78: "C21",
                79: "C22",
                80: "C23",
                81: "azide N attached to C23",
                82: "azide N 2",
                83: "azide N 3",
            }[map_id],
            "source_anchor_key": source_key,
            "source_model_atom": source_record["pdb_atom"],
            "source_model_serial": source_record["serial"],
            "corrected_candidate_atom": atom_record(selected_atoms[serial]),
            "anchor_basis": "Previously documented candidate linker/tail map from data/nnc2215-linker-mapping-audit.json; fixed for this subgraph search, not independently established here.",
        }

    def node_match(target_node: dict[str, Any], query_node: dict[str, Any]) -> bool:
        if target_node["element"] != query_node["element"]:
            return False
        fixed_serial = anchors.get(query_node["map_id"])
        return fixed_serial is None or target_node["serial"] == fixed_serial

    matcher = nx.algorithms.isomorphism.GraphMatcher(target, source_graph, node_match=node_match)
    embeddings: list[dict[int, int]] = []
    for target_to_query in matcher.subgraph_monomorphisms_iter():
        query_to_target = {query_id: target_id for target_id, query_id in target_to_query.items()}
        if len(query_to_target) != len(source_graph):
            continue
        if any(query_to_target.get(map_id) != serial for map_id, serial in anchors.items()):
            continue
        embeddings.append(query_to_target)
    embeddings.sort(
        key=lambda mapping: tuple(mapping[map_id] for map_id in sorted(mapping))
    )
    if len(embeddings) != EXPECTED_MATCHES_PER_STATE:
        raise ValueError(
            f"Expected {EXPECTED_MATCHES_PER_STATE} anchored embeddings for {model['state']}; found {len(embeddings)}"
        )

    match_records: list[dict[str, Any]] = []
    candidate_serials: dict[int, set[int]] = defaultdict(set)
    unmatched_by_match: list[list[dict[str, Any]]] = []
    for index, mapping in enumerate(embeddings, start=1):
        for map_id, serial in mapping.items():
            candidate_serials[map_id].add(serial)
        used = set(mapping.values())
        unmatched = sorted(
            (atom_record(atom) for serial, atom in selected_atoms.items() if serial not in used),
            key=lambda atom: (atom["chain"], atom["residue_number"], atom["name"], atom["serial"]),
        )
        unmatched_by_match.append(unmatched)
        match_records.append(
            {
                "embedding_id": f"embedding_{index:02d}",
                "map_id_to_pdb_serial": {
                    str(map_id): serial
                    for map_id, serial in sorted(mapping.items())
                },
            }
        )

    unmatched_signature = {
        tuple(sorted((atom["chain"], atom["residue_number"], atom["residue"], atom["name"]) for atom in items))
        for items in unmatched_by_match
    }
    if len(unmatched_signature) != 1:
        raise ValueError(f"Unmatched target atoms vary across embeddings in {model['state']}")
    unmatched_atoms = unmatched_by_match[0]
    unmatched_names = {atom["name"] for atom in unmatched_atoms}
    expected_names = EXPECTED_UNMATCHED_NAMES | ({"OD"} if model["state"] == "active_open" else {"O"})
    if unmatched_names != expected_names:
        raise ValueError(
            f"Unexpected unmatched fragment atoms in {model['state']}: {sorted(unmatched_names)}"
        )

    ambiguity = {
        str(map_id): {
            "candidate_count": len(candidate_serials[map_id]),
            "unique_across_all_embeddings": len(candidate_serials[map_id]) == 1,
            "candidate_pdb_serials": sorted(candidate_serials[map_id]),
        }
        for map_id in sorted(source_graph.nodes)
    }

    post_hash = sha256(pdb_path)
    if post_hash != observed_hash or post_hash != expected_coordinate_record["output_sha256"]:
        raise ValueError(f"PDB input changed during audit: {relative_path}")

    return {
        "state": model["state"],
        "coordinate_candidate_file": relative_path,
        "coordinate_candidate_sha256": observed_hash,
        "coordinate_candidate_sha256_expected_by_reconciliation_audit": expected_coordinate_record[
            "output_sha256"
        ],
        "selected_residues": ["LYX B30", "MC C1"],
        "selected_heavy_atom_count": len(target),
        "selected_heavy_element_counts": element_counts,
        "selected_adjacency_edge_count": target.number_of_edges(),
        "selected_connected_component_count": nx.number_connected_components(target),
        "selected_target_atom_inventory": {
            str(serial): atom_record(atom)
            for serial, atom in sorted(selected_atoms.items())
        },
        "source_graph_heavy_atom_count": source_graph.number_of_nodes(),
        "source_graph_adjacency_edge_count": source_graph.number_of_edges(),
        "anchor_provenance": {
            "audit_file": "data/nnc2215-linker-mapping-audit.json",
            "state": anchor_state["state"],
            "anchors_by_map_id": anchor_details,
        },
        "anchored_embedding_count": len(match_records),
        "embeddings": match_records,
        "embedding_map_encoding": "Each embedding stores map_id_to_pdb_serial; resolve serials through selected_target_atom_inventory.",
        "map_id_ambiguity": ambiguity,
        "ambiguous_graph_map_id_count": sum(
            not item["unique_across_all_embeddings"] for item in ambiguity.values()
        ),
        "unmatched_target_fragment_atoms": unmatched_atoms,
        "unmatched_fragment_interpretation": "Eight LYX B30 atoms attributed to the previously mapped six-carbon B29 linker (C6O2) under the candidate linker map; this audit does not independently prove that chemical interpretation.",
        "graph_match_scope": "Element-labeled graph monomorphism using CONECT/LINK adjacency only; extra target edges are allowed and bond orders are ignored.",
    }


def main() -> None:
    graph_text = GRAPH_PATH.read_text(encoding="utf-8").split()
    if not graph_text:
        raise ValueError(f"Empty mapped graph file: {GRAPH_PATH}")
    molecule = Chem.MolFromSmiles(graph_text[0])
    if molecule is None:
        raise ValueError(f"RDKit could not parse mapped graph: {GRAPH_PATH}")
    Chem.SanitizeMol(molecule)
    formula = rdMolDescriptors.CalcMolFormula(molecule)
    source_graph = build_source_graph(molecule)
    map_ids = sorted(source_graph.nodes)
    if formula != EXPECTED_FORMULA:
        raise ValueError(f"Expected formula {EXPECTED_FORMULA}; got {formula}")
    if len(source_graph) != EXPECTED_QUERY_HEAVY_ATOMS:
        raise ValueError(f"Expected 83 source graph heavy atoms; got {len(source_graph)}")
    if source_graph.number_of_edges() != EXPECTED_QUERY_EDGES:
        raise ValueError(f"Expected 89 source graph adjacency edges; got {source_graph.number_of_edges()}")
    if map_ids != list(range(1, EXPECTED_QUERY_HEAVY_ATOMS + 1)):
        raise ValueError("Source graph map IDs must be unique and complete from 1 through 83")

    graph_map = load_json(GRAPH_MAP_PATH)
    graph_audit = load_json(GRAPH_AUDIT_PATH)
    coordinate_audit = load_json(COORDINATE_AUDIT_PATH)
    anchor_audit = load_json(ANCHOR_AUDIT_PATH)
    graph_hash = sha256(GRAPH_PATH)
    if graph_audit["input"]["graph_sha256"] != graph_hash:
        raise ValueError("Mapped SMILES checksum does not match its graph audit")

    models_by_state = {model["state"]: model for model in coordinate_audit["models"]}
    anchor_states = {state["state"]: state for state in anchor_audit["states"]}
    if set(models_by_state) != {"active_open", "inactive_closed"}:
        raise ValueError(f"Unexpected coordinate candidate states: {sorted(models_by_state)}")
    if set(anchor_states) != set(models_by_state):
        raise ValueError("Coordinate candidates and anchor audit states do not match")

    states = [
        audit_state(
            models_by_state[state],
            anchor_states[state],
            source_graph,
            models_by_state[state],
        )
        for state in ("active_open", "inactive_closed")
    ]

    pdb_hashes_after = {
        model["output_file"]: sha256(ROOT / model["output_file"])
        for model in models_by_state.values()
    }
    if any(
        pdb_hashes_after[model["output_file"]] != model["output_sha256"]
        for model in models_by_state.values()
    ):
        raise ValueError("A coordinate candidate changed while the audit ran")

    result = {
        "schema": "dia-health-one.nnc2215-roof-graph-coordinate-mapping-audit.v1",
        "run": {
            "generated_at_utc": datetime.now(timezone.utc).isoformat(),
            "script": "research/audit_nnc2215_roof_graph_coordinate_mapping.py",
            "script_sha256": sha256(Path(__file__)),
            "python_version": platform.python_version(),
            "rdkit_version": rdBase.rdkitVersion,
            "networkx_version": nx.__version__,
            "gpu_used": False,
            "coordinate_files_modified": False,
        },
        "inputs": {
            "mapped_graph_file": "structures/nnc2215-roof-propyl-azide.smi",
            "mapped_graph_sha256": graph_hash,
            "atom_map_legend_file": "data/nnc2215-roof-macrocycle-atom-map.json",
            "atom_map_legend_sha256": sha256(GRAPH_MAP_PATH),
            "graph_validation_audit_file": "data/nnc2215-roof-macrocycle-graph-audit.json",
            "graph_validation_audit_sha256": sha256(GRAPH_AUDIT_PATH),
            "coordinate_reconciliation_audit_file": "data/nnc2215-carbon-reconciliation-candidate.json",
            "coordinate_reconciliation_audit_sha256": sha256(COORDINATE_AUDIT_PATH),
            "anchor_mapping_audit_file": "data/nnc2215-linker-mapping-audit.json",
            "anchor_mapping_audit_sha256": sha256(ANCHOR_AUDIT_PATH),
        },
        "source_graph": {
            "article_doi": graph_map["source"]["article_doi"],
            "figure": graph_map["source"]["figure"],
            "formula": formula,
            "heavy_atom_count": source_graph.number_of_nodes(),
            "adjacency_edge_count": source_graph.number_of_edges(),
            "connected_component_count": nx.number_connected_components(source_graph),
            "cycle_rank": source_graph.number_of_edges()
            - source_graph.number_of_nodes()
            + nx.number_connected_components(source_graph),
            "map_id_count": len(map_ids),
            "tail_source_numbering": graph_map["source"]["tail_source_numbering"],
        },
        "states": states,
        "assessment": {
            "both_coordinate_candidates_contain_anchored_roof_graph": all(
                state["anchored_embedding_count"] == EXPECTED_MATCHES_PER_STATE for state in states
            ),
            "one_authoritative_full_graph_to_pdb_map_established": False,
            "unmatched_b29_linker_atoms_per_state": 8,
            "interpretation": "The isolated 83-heavy-atom source roof graph embeds in each corrected coordinate candidate with 16 element-matched, edge-preserving mappings after fixing the six previously documented C21-C23/azide-tail anchors. The six anchors are assumptions imported from the candidate linker map. Many remaining graph map IDs have multiple PDB atom candidates; no unique full graph-to-PDB atom map follows.",
            "limitations": [
                "The check compares element labels and adjacency from PDB CONECT/LINK records only; PDB records do not encode authoritative bond orders or formal charges.",
                "Subgraph monomorphism allows additional target edges, including adjacency introduced when the source azide is clicked into the B29 alkyne linker.",
                "The six fixed tail anchors come from the previously documented candidate map and are not independently proven by this search.",
                "No coordinate fitting, geometry validation, stereochemical assignment, protonation assignment, partial-charge assignment, or force-field validation was performed.",
                "The source graph is an isolated roof macrocycle component, not the complete NNC2215 insulin conjugate.",
                "The eight unmatched LYX atoms are identified as the B29 linker only under the existing candidate partition; this audit does not independently prove that partition.",
                "The corrected PDB files remain candidate derivatives of published model coordinates; this audit does not promote them to experimental structures or production simulation inputs.",
            ],
        },
    }
    OUTPUT_PATH.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "audit": str(OUTPUT_PATH),
                "formula": formula,
                "source_heavy_atoms": source_graph.number_of_nodes(),
                "state_embedding_counts": {
                    state["state"]: state["anchored_embedding_count"] for state in states
                },
                "unique_full_atom_map": False,
                "networkx": nx.__version__,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
