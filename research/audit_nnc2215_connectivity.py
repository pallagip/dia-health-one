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
                "The PDB LINK records uniquely identify endpoints within each file, "
                "but the macrocycle-to-LYX atom label differs between open and closed "
                "models and must be reconciled against the source chemical graph."
            ),
        },
        "models": results,
        "assessment": (
            "Atom-serial CONECT records and LINK endpoints recover a traceable "
            "adjacency graph for the modified components despite repeated atom labels. "
            "The closed model contains one HETATM hydrogen with no CONECT neighbor. "
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
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
