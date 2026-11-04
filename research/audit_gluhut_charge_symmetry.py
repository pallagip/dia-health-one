#!/usr/bin/env python3
"""Quantify symmetry breaking in the deposited compiled GluHUT partial charges.

The public methyl-truncated host has a graph-reconciled two-cap/three-pillar
architecture. This audit enumerates automorphisms of its labelled heavy-atom
graph, derives exact graph-symmetry orbits, and measures charge variation within
those orbits. Attached hydrogens are compared by parent heavy-atom orbit and
AMBER hydrogen type.

The result describes the compiled RESP/GAFF model. It does not reconstruct RESP
inputs, establish that asymmetry is wrong, symmetrise production parameters, or
authorize transfer to PaG, MC-Q1, or complete NNC2215.
"""
from __future__ import annotations

import argparse
import collections
import hashlib
import json
import math
import stat
import zipfile
from pathlib import Path, PurePosixPath
from typing import Any

import networkx as nx

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import audit_gluhut_compiled_systems as amber  # noqa: E402
import audit_gluhut_graph_chemistry as chemistry  # noqa: E402
import audit_gluhut_smd_restraints as topology_reader  # noqa: E402

ARCHIVE_MD5 = "9f3c6733dc608d2bd37c488c560a5443"
ARCHIVE_SIZE_BYTES = 37931874
CHUNK_SIZE = 1024 * 1024
AUTOMORPHISM_LIMIT = 100_000
LARGE_CHARGE_RANGE_E = 0.05
EXACT_CHARGE_TOLERANCE_E = 1.0e-8


class ChargeSymmetryError(RuntimeError):
    pass


class UnionFind:
    def __init__(self, values: list[int]) -> None:
        self.parent = {value: value for value in values}

    def find(self, value: int) -> int:
        root = value
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[value] != value:
            next_value = self.parent[value]
            self.parent[value] = root
            value = next_value
        return root

    def union(self, left: int, right: int) -> None:
        left_root = self.find(left)
        right_root = self.find(right)
        if left_root != right_root:
            self.parent[max(left_root, right_root)] = min(left_root, right_root)

    def groups(self) -> list[list[int]]:
        grouped: dict[int, list[int]] = collections.defaultdict(list)
        for value in sorted(self.parent):
            grouped[self.find(value)].append(value)
        return sorted(grouped.values(), key=lambda group: (len(group), group), reverse=True)


def digest_file(path: Path, algorithm: str) -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(CHUNK_SIZE), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def canonical_hash(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("utf-8")
    ).hexdigest()


def safe_member(name: str) -> PurePosixPath:
    if "\\" in name:
        raise ChargeSymmetryError(f"backslash in member path: {name!r}")
    path = PurePosixPath(name)
    if path.is_absolute() or any(part in {"", ".", ".."} for part in path.parts):
        raise ChargeSymmetryError(f"unsafe member path: {name!r}")
    return path


def member_kind(info: zipfile.ZipInfo) -> str:
    if info.is_dir():
        return "directory"
    mode = (info.external_attr >> 16) & 0xFFFF
    if stat.S_ISLNK(mode):
        return "symlink"
    file_type = stat.S_IFMT(mode)
    if file_type not in {0, stat.S_IFREG}:
        return "special"
    return "file"


def infer_system(path: str) -> str:
    lower = path.lower()
    for name in ("fructofuranose", "fructopyranose", "galactose", "glucose"):
        if name in lower:
            return name
    raise ChargeSymmetryError(f"cannot infer system from {path}")


def parse_host(path: str, payload: bytes) -> dict[str, Any]:
    topology = topology_reader.topology_atoms(path, payload)
    sections = amber.parse_parm7_text(payload.decode("utf-8"))
    pointers = amber.pointer_map(sections)
    natom = pointers["NATOM"]
    raw_charges = [float(value) for value in amber.required(sections, "CHARGE", natom)[:natom]]
    charges = [value / amber.AMBER_CHARGE_SCALE for value in raw_charges]
    masses = [float(value) for value in amber.required(sections, "MASS", natom)[:natom]]
    bonds = chemistry.all_bonds(sections, pointers)
    host_atoms = sorted(
        index for index, atom in topology["atoms"].items() if atom["residue_label"] == "MOL"
    )
    host_set = set(host_atoms)
    adjacency = {index: set() for index in host_atoms}
    host_bonds = []
    for bond in bonds:
        left, right = bond["atoms_1_based"]
        if left in host_set and right in host_set:
            adjacency[left].add(right)
            adjacency[right].add(left)
            host_bonds.append(bond)
    atoms = {}
    for index in host_atoms:
        atom = dict(topology["atoms"][index])
        atom["charge_e"] = charges[index - 1]
        atom["mass"] = masses[index - 1]
        atom["full_degree"] = len(adjacency[index])
        atom["hydrogen_neighbour_count"] = sum(
            topology["atoms"][neighbor]["element"] == "H" for neighbor in adjacency[index]
        )
        atom["heavy_neighbour_count"] = sum(
            topology["atoms"][neighbor]["element"] != "H" for neighbor in adjacency[index]
        )
        atoms[index] = atom
    return {
        "path": path,
        "sha256": sha256_bytes(payload),
        "system_name": infer_system(path),
        "atoms": atoms,
        "bonds": host_bonds,
        "adjacency": adjacency,
        "host_atom_indices": host_atoms,
        "host_charge_e": sum(atoms[index]["charge_e"] for index in host_atoms),
        "charge_vector_sha256": canonical_hash([atoms[index]["charge_e"] for index in host_atoms]),
    }


def build_heavy_graph(host: dict[str, Any]) -> nx.Graph:
    graph = nx.Graph()
    heavy = {
        index for index in host["host_atom_indices"] if host["atoms"][index]["element"] != "H"
    }
    for index in sorted(heavy):
        atom = host["atoms"][index]
        graph.add_node(
            index,
            element=atom["element"],
            amber_atom_type=atom["amber_atom_type"],
            full_degree=atom["full_degree"],
            heavy_degree=atom["heavy_neighbour_count"],
            hydrogen_count=atom["hydrogen_neighbour_count"],
        )
    for bond in host["bonds"]:
        left, right = bond["atoms_1_based"]
        if left not in heavy or right not in heavy:
            continue
        graph.add_edge(
            left,
            right,
            force_constant=round(float(bond["force_constant"]), 8),
            equilibrium_A=round(float(bond["equilibrium_A"]), 8),
        )
    if not nx.is_connected(graph):
        raise ChargeSymmetryError("heavy-atom host graph is disconnected")
    return graph


def automorphism_orbits(graph: nx.Graph) -> tuple[int, list[list[int]]]:
    node_match = nx.algorithms.isomorphism.categorical_node_match(
        ["element", "amber_atom_type", "full_degree", "heavy_degree", "hydrogen_count"],
        [None, None, None, None, None],
    )
    edge_match = nx.algorithms.isomorphism.categorical_edge_match(
        ["force_constant", "equilibrium_A"], [None, None]
    )
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        graph, graph, node_match=node_match, edge_match=edge_match
    )
    union = UnionFind(list(graph.nodes))
    count = 0
    for mapping in matcher.isomorphisms_iter():
        count += 1
        if count > AUTOMORPHISM_LIMIT:
            raise ChargeSymmetryError(
                f"automorphism count exceeded safety limit {AUTOMORPHISM_LIMIT}"
            )
        for source, target in mapping.items():
            union.union(int(source), int(target))
    if count == 0:
        raise ChargeSymmetryError("no graph automorphism found")
    return count, union.groups()


def statistics(values: list[float]) -> dict[str, float]:
    if not values:
        raise ChargeSymmetryError("cannot summarize empty values")
    mean = sum(values) / len(values)
    variance = sum((value - mean) ** 2 for value in values) / len(values)
    return {
        "minimum": min(values),
        "maximum": max(values),
        "mean": mean,
        "population_standard_deviation": math.sqrt(variance),
        "range": max(values) - min(values),
        "maximum_absolute_deviation_from_mean": max(abs(value - mean) for value in values),
    }


def orbit_records(orbits: list[list[int]], host: dict[str, Any]) -> list[dict[str, Any]]:
    records = []
    for orbit_id, indices in enumerate(orbits, start=1):
        atoms = [host["atoms"][index] for index in indices]
        charges = [float(atom["charge_e"]) for atom in atoms]
        records.append(
            {
                "orbit_id": orbit_id,
                "orbit_size": len(indices),
                "atom_indices_1_based": indices,
                "atom_names": [atom["atom_name"] for atom in atoms],
                "element": atoms[0]["element"],
                "amber_atom_type": atoms[0]["amber_atom_type"],
                "charges_e": charges,
                "charge_statistics_e": statistics(charges),
                "exactly_equal_within_tolerance": (
                    max(charges) - min(charges) <= EXACT_CHARGE_TOLERANCE_E
                ),
                "range_exceeds_0_05e": (
                    max(charges) - min(charges) > LARGE_CHARGE_RANGE_E
                ),
            }
        )
    return sorted(
        records,
        key=lambda row: (
            row["charge_statistics_e"]["range"], row["orbit_size"], row["atom_indices_1_based"]
        ),
        reverse=True,
    )


def hydrogen_group_records(
    orbits: list[list[int]], host: dict[str, Any]
) -> list[dict[str, Any]]:
    orbit_by_atom = {
        atom: orbit_id
        for orbit_id, orbit in enumerate(orbits, start=1)
        for atom in orbit
    }
    groups: dict[tuple[int, str], list[int]] = collections.defaultdict(list)
    for parent in sorted(orbit_by_atom):
        for neighbor in sorted(host["adjacency"][parent]):
            atom = host["atoms"][neighbor]
            if atom["element"] == "H":
                groups[(orbit_by_atom[parent], atom["amber_atom_type"])].append(neighbor)
    records = []
    for (parent_orbit_id, atom_type), indices in sorted(groups.items()):
        unique_indices = sorted(set(indices))
        atoms = [host["atoms"][index] for index in unique_indices]
        charges = [float(atom["charge_e"]) for atom in atoms]
        parent_indices = sorted(
            {
                neighbor
                for index in unique_indices
                for neighbor in host["adjacency"][index]
                if neighbor in orbit_by_atom
            }
        )
        records.append(
            {
                "parent_heavy_orbit_id": parent_orbit_id,
                "hydrogen_amber_atom_type": atom_type,
                "hydrogen_count": len(unique_indices),
                "parent_heavy_atom_indices_1_based": parent_indices,
                "hydrogen_atom_indices_1_based": unique_indices,
                "hydrogen_atom_names": [atom["atom_name"] for atom in atoms],
                "charges_e": charges,
                "charge_statistics_e": statistics(charges),
                "exactly_equal_within_tolerance": (
                    max(charges) - min(charges) <= EXACT_CHARGE_TOLERANCE_E
                ),
                "range_exceeds_0_05e": (
                    max(charges) - min(charges) > LARGE_CHARGE_RANGE_E
                ),
            }
        )
    return sorted(
        records,
        key=lambda row: (
            row["charge_statistics_e"]["range"], row["hydrogen_count"]
        ),
        reverse=True,
    )


def symmetry_averaging_adjustment(
    heavy_records: list[dict[str, Any]],
    hydrogen_records: list[dict[str, Any]],
    host: dict[str, Any],
) -> dict[str, Any]:
    averaged = {
        index: float(host["atoms"][index]["charge_e"])
        for index in host["host_atom_indices"]
    }
    grouped_indices: set[int] = set()
    for record in heavy_records:
        mean = record["charge_statistics_e"]["mean"]
        for index in record["atom_indices_1_based"]:
            averaged[index] = mean
            grouped_indices.add(index)
    for record in hydrogen_records:
        mean = record["charge_statistics_e"]["mean"]
        for index in record["hydrogen_atom_indices_1_based"]:
            averaged[index] = mean
            grouped_indices.add(index)
    if grouped_indices != set(host["host_atom_indices"]):
        missing = sorted(set(host["host_atom_indices"]) - grouped_indices)
        raise ChargeSymmetryError(f"symmetry groups do not cover host atoms: {missing}")
    deltas = [
        averaged[index] - float(host["atoms"][index]["charge_e"])
        for index in host["host_atom_indices"]
    ]
    original_total = host["host_charge_e"]
    averaged_total = sum(averaged.values())
    return {
        "original_total_charge_e": original_total,
        "averaged_total_charge_e": averaged_total,
        "total_charge_difference_e": averaged_total - original_total,
        "rms_per_atom_adjustment_e": math.sqrt(
            sum(delta * delta for delta in deltas) / len(deltas)
        ),
        "mean_absolute_per_atom_adjustment_e": sum(abs(delta) for delta in deltas)
        / len(deltas),
        "maximum_absolute_per_atom_adjustment_e": max(abs(delta) for delta in deltas),
        "atoms_with_absolute_adjustment_above_0_05e": sum(
            abs(delta) > LARGE_CHARGE_RANGE_E for delta in deltas
        ),
    }


def audit(archive_path: Path, output_path: Path) -> dict[str, Any]:
    if archive_path.stat().st_size != ARCHIVE_SIZE_BYTES:
        raise ChargeSymmetryError(
            f"size mismatch: expected {ARCHIVE_SIZE_BYTES}, got {archive_path.stat().st_size}"
        )
    md5 = digest_file(archive_path, "md5")
    if md5 != ARCHIVE_MD5:
        raise ChargeSymmetryError(f"MD5 mismatch: expected {ARCHIVE_MD5}, got {md5}")
    archive_sha256 = digest_file(archive_path, "sha256")

    hosts: list[dict[str, Any]] = []
    seen: set[str] = set()
    with zipfile.ZipFile(archive_path) as archive:
        bad = archive.testzip()
        if bad:
            raise ChargeSymmetryError(f"ZIP CRC failure: {bad}")
        for info in archive.infolist():
            safe_member(info.filename)
            if info.filename in seen:
                raise ChargeSymmetryError(f"duplicate member: {info.filename}")
            seen.add(info.filename)
            if info.flag_bits & 0x1:
                raise ChargeSymmetryError(f"encrypted member: {info.filename}")
            kind = member_kind(info)
            if kind in {"symlink", "special"}:
                raise ChargeSymmetryError(f"unsupported {kind}: {info.filename}")
            if kind != "file" or PurePosixPath(info.filename).suffix.lower() not in {
                ".parm7",
                ".prmtop",
            }:
                continue
            hosts.append(parse_host(info.filename, archive.read(info)))
    hosts.sort(key=lambda row: row["system_name"])
    if len(hosts) != 4:
        raise ChargeSymmetryError(f"expected four host topologies, found {len(hosts)}")
    charge_hashes = {host["charge_vector_sha256"] for host in hosts}
    if len(charge_hashes) != 1:
        raise ChargeSymmetryError("host charge vectors differ across sugar systems")

    reference = hosts[0]
    graph = build_heavy_graph(reference)
    automorphism_count, orbits = automorphism_orbits(graph)
    heavy_records = orbit_records(orbits, reference)
    hydrogen_records = hydrogen_group_records(orbits, reference)
    averaging = symmetry_averaging_adjustment(
        heavy_records, hydrogen_records, reference
    )
    nontrivial_heavy = [row for row in heavy_records if row["orbit_size"] > 1]
    nontrivial_hydrogen = [
        row for row in hydrogen_records if row["hydrogen_count"] > 1
    ]
    max_heavy_range = max(
        row["charge_statistics_e"]["range"] for row in nontrivial_heavy
    )
    max_hydrogen_range = max(
        row["charge_statistics_e"]["range"] for row in nontrivial_hydrogen
    )

    report = {
        "schema": "dia-health-one.gluhut-charge-symmetry.v1",
        "status": "compiled_charge_symmetry_quantified_source_resp_derivation_missing",
        "scope": "graph-equivalent atom charge variation in the public methyl-truncated GluHUT host",
        "not_target_claim": (
            "not NNC2215, PaG, MC-Q1, binding affinity, efficacy, safety, or a drug result"
        ),
        "archive": {
            "filename": archive_path.name,
            "size_bytes": archive_path.stat().st_size,
            "md5": md5,
            "sha256": archive_sha256,
            "crc_verified": True,
        },
        "software": {
            "python": sys.version.split()[0],
            "networkx": nx.__version__,
        },
        "cross_system": {
            "topology_count": len(hosts),
            "systems": [host["system_name"] for host in hosts],
            "host_charge_vector_sha256": next(iter(charge_hashes)),
            "host_charge_vectors_identical": True,
            "host_total_charges_e": {
                host["system_name"]: host["host_charge_e"] for host in hosts
            },
        },
        "graph_symmetry": {
            "heavy_atom_count": graph.number_of_nodes(),
            "heavy_bond_count": graph.number_of_edges(),
            "automorphism_count": automorphism_count,
            "heavy_orbit_count": len(orbits),
            "nontrivial_heavy_orbit_count": len(nontrivial_heavy),
            "heavy_orbit_size_distribution": dict(
                sorted(
                    collections.Counter(len(orbit) for orbit in orbits).items()
                )
            ),
            "node_labels": [
                "element",
                "AMBER atom type",
                "full degree",
                "heavy degree",
                "attached hydrogen count",
            ],
            "edge_labels": ["bond force constant", "equilibrium distance"],
        },
        "heavy_atom_orbits": heavy_records,
        "hydrogen_groups": hydrogen_records,
        "summary": {
            "maximum_nontrivial_heavy_orbit_charge_range_e": max_heavy_range,
            "maximum_nontrivial_hydrogen_group_charge_range_e": max_hydrogen_range,
            "heavy_orbits_with_range_above_0_05e": sum(
                row["range_exceeds_0_05e"] for row in nontrivial_heavy
            ),
            "hydrogen_groups_with_range_above_0_05e": sum(
                row["range_exceeds_0_05e"] for row in nontrivial_hydrogen
            ),
            "all_nontrivial_heavy_orbits_exactly_equal": all(
                row["exactly_equal_within_tolerance"] for row in nontrivial_heavy
            ),
            "all_nontrivial_hydrogen_groups_exactly_equal": all(
                row["exactly_equal_within_tolerance"] for row in nontrivial_hydrogen
            ),
            "symmetry_averaging_adjustment": averaging,
        },
        "decision": {
            "heavy_graph_automorphisms_enumerated": True,
            "nontrivial_graph_symmetry_present": bool(nontrivial_heavy),
            "compiled_charges_are_exactly_graph_symmetry_averaged": (
                all(row["exactly_equal_within_tolerance"] for row in nontrivial_heavy)
                and all(
                    row["exactly_equal_within_tolerance"]
                    for row in nontrivial_hydrogen
                )
            ),
            "large_graph_equivalent_charge_asymmetry_detected": (
                max(max_heavy_range, max_hydrogen_range) > LARGE_CHARGE_RANGE_E
            ),
            "source_resp_derivation_reconstructable": False,
            "automatic_charge_symmetrisation_authorized": False,
            "direct_transfer_to_pag_or_mcq1_authorized": False,
            "parent_pag_calibration_authorized": False,
            "candidate_ranking_authorized": False,
            "next_action": (
                "request source RESP conformers, ESP grids, equivalencing constraints and final mol2; "
                "use compiled charges only for exact public benchmark reproduction"
            ),
        },
        "limitations": [
            "Graph-equivalent atoms can receive different RESP charges when charges are fit to an asymmetric conformation without equivalencing constraints.",
            "Charge asymmetry is descriptive and is not by itself evidence that the published parameterisation is invalid.",
            "Heavy-atom automorphisms use compiled atom types and bonded parameters as graph labels; alternative source atom typing could define different equivalence classes.",
            "Hydrogen groups are defined from heavy-parent graph orbits and AMBER hydrogen type rather than full all-hydrogen automorphism enumeration.",
            "No production parameter file is modified or symmetrised by this audit.",
        ],
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return report


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument(
        "--output-json",
        type=Path,
        default=Path("research/data/gluhut-charge-symmetry.json"),
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    report = audit(args.archive, args.output_json)
    print(
        json.dumps(
            {
                "cross_system": report["cross_system"],
                "graph_symmetry": report["graph_symmetry"],
                "summary": report["summary"],
                "decision": report["decision"],
                "largest_heavy_orbits": report["heavy_atom_orbits"][:10],
                "largest_hydrogen_groups": report["hydrogen_groups"][:10],
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
