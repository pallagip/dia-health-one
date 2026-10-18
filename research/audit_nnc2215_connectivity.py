#!/usr/bin/env python3
"""Audit atom-serial connectivity in publisher-supplied NNC2215 model PDBs.

This records the adjacency encoded by PDB CONECT/LINK records. It deliberately
does not infer bond orders, formal charges, protonation, or force-field terms.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent
MANIFEST_PATH = ROOT / "data" / "structure-manifest.json"
OUTPUT_PATH = ROOT / "data" / "nnc2215-connectivity-audit.json"
WWPDB_CONNECTIVITY_URL = (
    "https://www.wwpdb.org/documentation/file-format-content/format23/sect10.html"
)
def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_atom(line: str) -> dict[str, object]:
    return {
        "serial": int(line[6:11]),
        "record": line[:6].strip(),
        "name": line[12:16].strip(),
        "altloc": line[16].strip(),
        "residue": line[17:20].strip(),
        "chain": line[21].strip(),
        "residue_number": line[22:26].strip(),
        "insertion_code": line[26].strip(),
        "element": line[76:78].strip().upper(),
        "xyz_angstrom": [float(line[i : i + 8]) for i in (30, 38, 46)],
    }


def residue_key(atom: dict[str, object]) -> tuple[str, str, str, str]:
    return (
        str(atom["chain"]),
        str(atom["residue_number"]),
        str(atom["insertion_code"]),
        str(atom["residue"]),
    )


def parse_link_endpoint(line: str, first: bool) -> tuple[str, str, str, str, str, str]:
    if first:
        atom_name, altloc = line[12:16].strip(), line[16].strip()
        residue, chain = line[17:20].strip(), line[21].strip()
        residue_number, insertion = line[22:26].strip(), line[26].strip()
    else:
        atom_name, altloc = line[42:46].strip(), line[46].strip()
        residue, chain = line[47:50].strip(), line[51].strip()
        residue_number, insertion = line[52:56].strip(), line[56].strip()
    return chain, residue_number, insertion, residue, atom_name, altloc


def endpoint_matches(
    atoms: dict[int, dict[str, object]],
    endpoint: tuple[str, str, str, str, str, str],
) -> list[int]:
    chain, residue_number, insertion, residue, atom_name, altloc = endpoint
    return [
        serial
        for serial, atom in atoms.items()
        if atom["chain"] == chain
        and atom["residue_number"] == residue_number
        and atom["insertion_code"] == insertion
        and atom["residue"] == residue
        and atom["name"] == atom_name
        and atom["altloc"] == altloc
    ]


def audit_model(model: dict[str, object]) -> dict[str, object]:
    pdb_path = ROOT / str(model["path"])
    observed_sha = sha256(pdb_path)
    if observed_sha != model["sha256"]:
        raise ValueError(f"Source PDB checksum mismatch for {pdb_path}")

    atoms: dict[int, dict[str, object]] = {}
    directed_conect: dict[int, set[int]] = defaultdict(set)
    raw_conect_reference_count = 0
    link_records: list[str] = []
    with pdb_path.open(encoding="ascii") as stream:
        for line in stream:
            if line.startswith(("ATOM  ", "HETATM")):
                atom = parse_atom(line)
                serial = int(atom["serial"])
                if serial in atoms:
                    raise ValueError(f"Duplicate atom serial {serial} in {pdb_path}")
                atoms[serial] = atom
            elif line.startswith("CONECT"):
                values = [
                    int(line[offset : offset + 5])
                    for offset in range(6, len(line), 5)
                    if line[offset : offset + 5].strip()
                ]
                if not values:
                    continue
                for target in values[1:]:
                    directed_conect[values[0]].add(target)
                    raw_conect_reference_count += 1
            elif line.startswith("LINK  "):
                link_records.append(line.rstrip("\n"))

    missing_conect_serials = sorted(
        {source for source in directed_conect if source not in atoms}
        | {
            target
            for targets in directed_conect.values()
            for target in targets
            if target not in atoms
        }
    )
    if missing_conect_serials:
        raise ValueError(f"CONECT references absent serials in {pdb_path}: {missing_conect_serials}")

    # CONECT is redundant by design. Collapse directed references into unique
    # undirected adjacency pairs; do not treat repeated fields as bond order.
    conect_pairs = {
        tuple(sorted((source, target)))
        for source, targets in directed_conect.items()
        for target in targets
    }
    hetero_serials = {
        serial
        for serial, atom in atoms.items()
        if atom["record"] == "HETATM"
    }
    hetero_neighbors: dict[int, set[int]] = defaultdict(set)
    for left, right in conect_pairs:
        if left in hetero_serials:
            hetero_neighbors[left].add(right)
        if right in hetero_serials:
            hetero_neighbors[right].add(left)

    by_name: dict[tuple[tuple[str, str, str, str], str], list[int]] = defaultdict(list)
    for serial in hetero_serials:
        atom = atoms[serial]
        by_name[(residue_key(atom), str(atom["name"]))].append(serial)
    duplicates = []
    for (reskey, atom_name), serials in sorted(by_name.items()):
        if len(serials) < 2:
            continue
        entries = []
        for serial in sorted(serials):
            neighbors = []
            for neighbor in sorted(hetero_neighbors.get(serial, set())):
                neighbor_atom = atoms[neighbor]
                neighbors.append(
                    {
                        "serial": neighbor,
                        "atom_name": neighbor_atom["name"],
                        "residue": neighbor_atom["residue"],
                    }
                )
            entries.append({"serial": serial, "neighbors": neighbors})
        duplicates.append(
            {
                "residue_key": list(reskey),
                "atom_name": atom_name,
                "serial_instances": entries,
                "serial_connectivity_disambiguates_instances": all(
                    len(entry["neighbors"]) > 0 for entry in entries
                ),
            }
        )

    resolved_links = []
    ambiguous_links = []
    for line in link_records:
        first = parse_link_endpoint(line, True)
        second = parse_link_endpoint(line, False)
        first_matches = endpoint_matches(atoms, first)
        second_matches = endpoint_matches(atoms, second)
        entry: dict[str, object] = {
            "record": line,
            "first_endpoint": list(first),
            "second_endpoint": list(second),
            "first_serial_matches": first_matches,
            "second_serial_matches": second_matches,
        }
        if len(first_matches) == 1 and len(second_matches) == 1:
            a, b = atoms[first_matches[0]], atoms[second_matches[0]]
            actual_distance = math.dist(a["xyz_angstrom"], b["xyz_angstrom"])
            annotated = line[73:78].strip()
            entry.update(
                {
                    "endpoint_resolution": "unique",
                    "actual_coordinate_distance_angstrom": round(actual_distance, 3),
                    "pdb_annotated_distance_angstrom": (
                        float(annotated) if annotated else None
                    ),
                    "both_endpoints_nonstandard": (
                        a["record"] == "HETATM" and b["record"] == "HETATM"
                    ),
                }
            )
            resolved_links.append(entry)
        else:
            entry["endpoint_resolution"] = "ambiguous_or_missing"
            ambiguous_links.append(entry)

    hetero_by_residue: dict[str, Counter[str]] = defaultdict(Counter)
    for serial in hetero_serials:
        atom = atoms[serial]
        hetero_by_residue[str(atom["residue"])][str(atom["element"])] += 1

    orphan_hetero_serials = sorted(
        serial for serial in hetero_serials if not hetero_neighbors.get(serial)
    )
    unresolved_hetero_atoms = [
        {
            "serial": serial,
            "atom_name": atoms[serial]["name"],
            "residue": atoms[serial]["residue"],
            "chain": atoms[serial]["chain"],
            "residue_number": atoms[serial]["residue_number"],
            "element": atoms[serial]["element"],
        }
        for serial in orphan_hetero_serials
    ]

    duplicate_label_instances = sum(len(x["serial_instances"]) for x in duplicates)
    return {
        "state": model["state"],
        "file": model["path"],
        "sha256": observed_sha,
        "atom_serials_unique": True,
        "coordinate_atom_count": len(atoms),
        "hetero_atom_count": len(hetero_serials),
        "hetero_elements_by_residue": {
            residue: dict(sorted(counts.items()))
            for residue, counts in sorted(hetero_by_residue.items())
        },
        "conect": {
            "directed_reference_count": raw_conect_reference_count,
            "unique_undirected_pair_count": len(conect_pairs),
            "hetero_atoms_with_at_least_one_conect_neighbor": sum(
                bool(hetero_neighbors.get(serial)) for serial in hetero_serials
            ),
            "hetero_atoms_with_no_conect_neighbor": unresolved_hetero_atoms,
            "interpretation": (
                "CONECT identifies bonded atom serials. Repeated atom labels do not "
                "erase this serial-level adjacency, but CONECT does not supply bond "
                "orders or formal charges."
            ),
        },
        "duplicate_atom_labels": {
            "duplicate_name_group_count": len(duplicates),
            "instances_in_duplicate_groups": duplicate_label_instances,
            "all_instances_have_serial_connectivity": all(
                x["serial_connectivity_disambiguates_instances"] for x in duplicates
            ),
            "groups": duplicates,
        },
        "link_records": {
            "count": len(link_records),
            "resolved_unique_endpoint_count": len(resolved_links),
            "ambiguous_or_missing_endpoint_count": len(ambiguous_links),
            "resolved": resolved_links,
            "ambiguous_or_missing": ambiguous_links,
        },
    }


def load_serial_graph(model: dict[str, object]) -> tuple[dict[int, dict[str, object]], set[tuple[int, int]], set[tuple[int, int]]]:
    """Load atom records, CONECT pairs, and LINK-augmented adjacency for comparison."""
    pdb_path = ROOT / str(model["path"])
    atoms: dict[int, dict[str, object]] = {}
    conect_pairs: set[tuple[int, int]] = set()
    link_records: list[str] = []
    with pdb_path.open(encoding="ascii") as stream:
        for line in stream:
            if line.startswith(("ATOM  ", "HETATM")):
                atom = parse_atom(line)
                atoms[int(atom["serial"])] = atom
            elif line.startswith("CONECT"):
                values = [
                    int(line[offset : offset + 5])
                    for offset in range(6, len(line), 5)
                    if line[offset : offset + 5].strip()
                ]
                if values:
                    conect_pairs.update(
                        tuple(sorted((values[0], target)))
                        for target in values[1:]
                    )
            elif line.startswith("LINK  "):
                link_records.append(line.rstrip("\n"))

    all_pairs = set(conect_pairs)
    for line in link_records:
        first = endpoint_matches(atoms, parse_link_endpoint(line, True))
        second = endpoint_matches(atoms, parse_link_endpoint(line, False))
        if len(first) == 1 and len(second) == 1:
            all_pairs.add(tuple(sorted((first[0], second[0]))))
    return atoms, conect_pairs, all_pairs


def named_atom(
    atoms: dict[int, dict[str, object]],
    chain: str,
    residue_number: str,
    residue: str,
    atom_name: str,
) -> dict[str, object]:
    matches = [
        atom for atom in atoms.values()
        if atom["chain"] == chain
        and atom["residue_number"] == residue_number
        and atom["residue"] == residue
        and atom["name"] == atom_name
    ]
    if len(matches) != 1:
        raise ValueError(
            f"Expected one {chain}:{residue_number} {residue} {atom_name}; found {len(matches)}."
        )
    return matches[0]


def cross_state_atom_record_comparison(
    models: list[dict[str, object]],
) -> dict[str, object]:
    """Compare the linker roof and N-terminal hydrogen records across model states."""
    by_state = {str(model["state"]): model for model in models}
    open_atoms, open_conect, open_graph = load_serial_graph(by_state["active_open"])
    closed_atoms, closed_conect, closed_graph = load_serial_graph(by_state["inactive_closed"])

    def neighborhood(
        atoms: dict[int, dict[str, object]],
        graph: set[tuple[int, int]],
        target: dict[str, object],
    ) -> list[dict[str, object]]:
        serial = int(target["serial"])
        neighbors = []
        for left, right in graph:
            neighbor_serial = right if left == serial else left if right == serial else None
            if neighbor_serial is None:
                continue
            atom = atoms[neighbor_serial]
            neighbors.append(
                {
                    "serial": neighbor_serial,
                    "atom_name": atom["name"],
                    "residue": atom["residue"],
                    "element": atom["element"],
                }
            )
        return sorted(neighbors, key=lambda item: (str(item["residue"]), str(item["atom_name"]), int(item["serial"])))

    open_roof = named_atom(open_atoms, "B", "30", "LYX", "CD6")
    closed_roof = named_atom(closed_atoms, "B", "30", "LYX", "C")
    open_roof_neighbors = neighborhood(open_atoms, open_graph, open_roof)
    closed_roof_neighbors = neighborhood(closed_atoms, closed_graph, closed_roof)
    open_roof_elements = sorted(str(item["element"]) for item in open_roof_neighbors)
    closed_roof_elements = sorted(str(item["element"]) for item in closed_roof_neighbors)

    def nearest_hydrogen_to_phenyl_n(
        atoms: dict[int, dict[str, object]], conect: set[tuple[int, int]]
    ) -> dict[str, object]:
        nitrogen = named_atom(atoms, "B", "1", "PHE", "N")
        hydrogen = min(
            (atom for atom in atoms.values() if atom["element"] == "H"),
            key=lambda atom: math.dist(nitrogen["xyz_angstrom"], atom["xyz_angstrom"]),
        )
        pair = tuple(sorted((int(nitrogen["serial"]), int(hydrogen["serial"]))))
        return {
            "nitrogen": {
                "serial": nitrogen["serial"],
                "residue": nitrogen["residue"],
                "atom_name": nitrogen["name"],
            },
            "nearest_hydrogen": {
                "serial": hydrogen["serial"],
                "record": hydrogen["record"],
                "residue": hydrogen["residue"],
                "residue_number": hydrogen["residue_number"],
                "atom_name": hydrogen["name"],
            },
            "distance_angstrom": round(
                math.dist(nitrogen["xyz_angstrom"], hydrogen["xyz_angstrom"]), 3
            ),
            "covalent_adjacency_recorded_in_conect": pair in conect,
        }

    def lyx_cd2_hydrogens(
        atoms: dict[int, dict[str, object]], conect: set[tuple[int, int]]
    ) -> list[dict[str, object]]:
        carbon = named_atom(atoms, "B", "30", "LYX", "CD2")
        serial = int(carbon["serial"])
        return sorted(
            [
                {
                    "serial": neighbor_serial,
                    "atom_name": atoms[neighbor_serial]["name"],
                }
                for left, right in conect
                if (neighbor_serial := (right if left == serial else left if right == serial else None)) is not None
                and atoms[neighbor_serial]["element"] == "H"
            ],
            key=lambda item: (str(item["atom_name"]), int(item["serial"])),
        )

    open_hydrogens = lyx_cd2_hydrogens(open_atoms, open_conect)
    closed_hydrogens = lyx_cd2_hydrogens(closed_atoms, closed_conect)
    open_n_terminal_h = nearest_hydrogen_to_phenyl_n(open_atoms, open_conect)
    closed_n_terminal_h = nearest_hydrogen_to_phenyl_n(closed_atoms, closed_conect)

    return {
        "macrocycle_roof_local_neighborhood": {
            "open": {
                "atom": {"serial": open_roof["serial"], "name": open_roof["name"], "element": open_roof["element"]},
                "neighbors_including_link_records": open_roof_neighbors,
            },
            "closed": {
                "atom": {"serial": closed_roof["serial"], "name": closed_roof["name"], "element": closed_roof["element"]},
                "neighbors_including_link_records": closed_roof_neighbors,
            },
            "neighbor_element_multisets_match": open_roof_elements == closed_roof_elements,
            "interpretation": (
                "The open CD6 and closed C endpoint labels have the same local element-neighbor pattern, "
                "including the MC C31 LINK. This supports a local atom-label correspondence but does not "
                "prove full conjugate graph equivalence."
            ),
        },
        "phenyl_b1_n_terminal_hydrogen_record": {
            "open": open_n_terminal_h,
            "closed": closed_n_terminal_h,
            "interpretation": (
                "The nearest H to PheB1 N occupies the corresponding local position in both models, "
                "but the closed PDB assigns it to PHX and omits a CONECT edge while the open PDB "
                "records it as PHE B1 H2 with an N-H CONECT edge. Treat this as a model-record "
                "discrepancy pending source-level reconciliation, not as validated topology."
            ),
        },
        "lyx_cd2_hydrogen_neighbors": {
            "open": open_hydrogens,
            "closed": closed_hydrogens,
            "counts_match": len(open_hydrogens) == len(closed_hydrogens),
            "interpretation": (
                "The open model records two CONECT hydrogen neighbors for LYX CD2; the closed model "
                "records one. The source coordinates do not establish whether the missing neighbor "
                "was omitted or assigned elsewhere."
            ),
        },
        "full_model_topology_equivalence_established": False,
    }


def selected_heavy_atom_adjacency_comparison(
    models: list[dict[str, object]],
) -> dict[str, object]:
    """Compare CONECT/LINK heavy-atom adjacency with explicit LYX label aliases.

    The comparison covers PHX, LYX, and MC heavy atoms plus the PheB1 nitrogen
    and LysB29 NZ attachment atoms. It excludes hydrogens and free glucose.
    """
    by_state = {str(model["state"]): model for model in models}
    aliases = {
        "active_open": {
            ("LYX", "CD6"): "roof_carbon_linked_to_MC_C31",
            ("LYX", "OD"): "terminal_oxygen_with_CD1_neighbor",
        },
        "inactive_closed": {
            ("LYX", "C"): "roof_carbon_linked_to_MC_C31",
            ("LYX", "O"): "terminal_oxygen_with_CD1_neighbor",
        },
    }

    def graph_for_state(state: str) -> tuple[set[tuple[str, ...]], set[tuple[tuple[str, ...], tuple[str, ...]]]]:
        atoms, _, graph = load_serial_graph(by_state[state])
        serial_to_key: dict[int, tuple[str, ...]] = {}
        for serial, atom in atoms.items():
            residue = str(atom["residue"])
            chain = str(atom["chain"])
            residue_number = str(atom["residue_number"])
            element = str(atom["element"])
            is_modified_heavy_atom = residue in {"PHX", "LYX", "MC"} and element != "H"
            is_b1_anchor = residue == "PHE" and chain == "B" and residue_number == "1" and atom["name"] == "N"
            is_b29_anchor = residue == "LYS" and chain == "B" and residue_number == "29" and atom["name"] == "NZ"
            if not (is_modified_heavy_atom or is_b1_anchor or is_b29_anchor):
                continue
            atom_name = str(atom["name"])
            atom_name = aliases[state].get((residue, atom_name), atom_name)
            key = (chain, residue_number, residue, atom_name, element)
            if key in serial_to_key.values():
                raise ValueError(f"Duplicate normalized heavy-atom key in {state}: {key}")
            serial_to_key[serial] = key

        node_keys = set(serial_to_key.values())
        edge_keys = {
            tuple(sorted((serial_to_key[left], serial_to_key[right])))
            for left, right in graph
            if left in serial_to_key and right in serial_to_key
        }
        return node_keys, edge_keys

    open_nodes, open_edges = graph_for_state("active_open")
    closed_nodes, closed_edges = graph_for_state("inactive_closed")
    nodes_match = open_nodes == closed_nodes
    edges_match = open_edges == closed_edges
    return {
        "scope": (
            "PHX, LYX, and MC heavy atoms plus PheB1 N and LysB29 NZ; excludes all hydrogens, "
            "free glucose, and the remainder of the insulin protein."
        ),
        "explicit_atom_name_aliases": [
            {"open": "LYX B30 CD6", "closed": "LYX B30 C", "basis": "same element-neighbor pattern and MC C31 LINK"},
            {"open": "LYX B30 OD", "closed": "LYX B30 O", "basis": "same element and sole heavy-atom neighbor LYX CD1"},
        ],
        "open_node_count": len(open_nodes),
        "closed_node_count": len(closed_nodes),
        "open_adjacency_edge_count": len(open_edges),
        "closed_adjacency_edge_count": len(closed_edges),
        "normalized_node_sets_match": nodes_match,
        "normalized_adjacency_edge_sets_match": edges_match,
        "selected_heavy_atom_adjacency_matches": nodes_match and edges_match,
        "interpretation": (
            "The selected 113-node heavy-atom CONECT/LINK adjacency graph matches after the two "
            "explicit LYX atom-name aliases. This is an adjacency comparison only; it does not "
            "validate bond orders, stereochemistry, protonation, formal charges, or force-field parameters."
            if nodes_match and edges_match
            else "The selected heavy-atom node or adjacency sets differ after the explicit LYX atom-name aliases."
        ),
    }


def main() -> None:
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    model_records = manifest["files"]
    results = [audit_model(model) for model in model_records]

    links_by_state = {
        model["state"]: [
            (
                link["first_endpoint"][3],
                link["first_endpoint"][4],
                link["second_endpoint"][3],
                link["second_endpoint"][4],
            )
            for link in model["link_records"]["resolved"]
        ]
        for model in results
    }
    output = {
        "schema": "dia-health-one.nnc2215-connectivity-audit.v1",
        "run": {
            "run_id": "DHO-NNC2215-CONNECTIVITY-"
            + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ"),
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
            "script": "research/audit_nnc2215_connectivity.py",
            "script_sha256": sha256(Path(__file__)),
            "runtime_python": __import__("sys").version.split()[0],
            "nvidia_gpu_used": False,
        },
        "source": {
            "manifest": "research/data/structure-manifest.json",
            "dataset_doi": manifest["dataset_doi"],
            "dataset_url": manifest["dataset_url"],
            "wwpdb_pdb_connectivity_spec": WWPDB_CONNECTIVITY_URL,
        },
        "bond_order_and_charge_evidence": {
            "conect_encodes_bond_order": False,
            "conect_encodes_formal_charge": False,
            "source_basis": (
                "The wwPDB PDB v2.3 connectivity specification defines CONECT as "
                "atom-serial adjacency and states that it does not distinguish "
                "atoms with excess positive or negative charge."
            ),
            "production_bond_order_inference_performed": False,
        },
        "open_closed_link_comparison": {
            "open_links": links_by_state.get("active_open", []),
            "closed_links": links_by_state.get("inactive_closed", []),
            "macrocycle_to_lyx_link_differs": (
                links_by_state.get("active_open", [])[0:1]
                != links_by_state.get("inactive_closed", [])[0:1]
            ),
            "interpretation": (
                "The PDB LINK records identify endpoints within each file. The macrocycle-to-LYX "
                "atom label differs, but the selected heavy-atom CONECT/LINK adjacency graphs "
                "match after the documented local LYX label aliases. This does not establish "
                "bond orders or a force-field-ready chemical graph."
            ),
        },
        "cross_state_atom_record_comparison": cross_state_atom_record_comparison(model_records),
        "selected_heavy_atom_adjacency_comparison": selected_heavy_atom_adjacency_comparison(model_records),
        "models": results,
        "assessment": (
            "Atom-serial CONECT records and LINK endpoints recover a traceable "
            "selected heavy-atom adjacency graph for the modified components despite repeated atom labels. "
            "After explicit normalization of two LYX atom-name differences, the open and closed "
            "selected heavy-atom node and edge sets match. Hydrogen-level correspondence remains "
            "unresolved: the closed model assigns the PheB1 N-terminal H to PHX without a CONECT "
            "edge and records one fewer hydrogen neighbor on LYX CD2. Full topology equivalence is not established. "
            "Bond orders, protonation, formal charges, stereochemical validation, and "
            "force-field parameters remain unresolved; the models are not production-MD-ready."
        ),
        "target_binding_md_started": False,
        "clinical_use": False,
    }
    OUTPUT_PATH.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "output": str(OUTPUT_PATH),
                "assessment": output["assessment"],
                "models": [
                    {
                        "state": model["state"],
                        "coordinate_atom_count": model["coordinate_atom_count"],
                        "hetero_atom_count": model["hetero_atom_count"],
                        "duplicate_atom_name_groups": model["duplicate_atom_labels"][
                            "duplicate_name_group_count"
                        ],
                        "duplicate_instances_with_unique_conect_neighbors": model[
                            "duplicate_atom_labels"
                        ]["all_instances_have_serial_connectivity"],
                        "hetero_atoms_without_conect_neighbor": model["conect"][
                            "hetero_atoms_with_no_conect_neighbor"
                        ],
                        "resolved_link_records": model["link_records"][
                            "resolved_unique_endpoint_count"
                        ],
                        "ambiguous_link_records": model["link_records"][
                            "ambiguous_or_missing_endpoint_count"
                        ],
                    }
                    for model in results
                ],
                "macrocycle_attachment_labels_differ": output[
                    "open_closed_link_comparison"
                ]["macrocycle_to_lyx_link_differs"],
                "macrocycle_attachment_local_neighbor_pattern_matches": output[
                    "cross_state_atom_record_comparison"
                ]["macrocycle_roof_local_neighborhood"]["neighbor_element_multisets_match"],
                "selected_heavy_atom_adjacency": output[
                    "selected_heavy_atom_adjacency_comparison"
                ],
                "open_closed_phe_b1_nearest_hydrogen": output[
                    "cross_state_atom_record_comparison"
                ]["phenyl_b1_n_terminal_hydrogen_record"],
                "lyx_cd2_hydrogen_neighbor_counts": {
                    "open": len(output["cross_state_atom_record_comparison"]["lyx_cd2_hydrogen_neighbors"]["open"]),
                    "closed": len(output["cross_state_atom_record_comparison"]["lyx_cd2_hydrogen_neighbors"]["closed"]),
                },
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
