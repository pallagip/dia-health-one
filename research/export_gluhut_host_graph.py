#!/usr/bin/env python3
"""Export the exact compiled heavy-atom graph of the public GluHUT host.

The output is derived from the checksum-verified AMBER parm7 systems in the
public SMD archive. It records element/connectivity identity separately from
compiled charges and force-field terms so private downstream work can compare
constitutional graphs without exposing private candidate structures.

AMBER topology files do not encode chemical bond order. The exported edges are
therefore connectivity plus compiled equilibrium length/force constant, not a
bond-order reconstruction. This is not an affinity or candidate result.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import stat
import zipfile
from pathlib import Path, PurePosixPath
from typing import Any

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import audit_gluhut_compiled_systems as amber  # noqa: E402
import audit_gluhut_graph_chemistry as chemistry  # noqa: E402

ARCHIVE_MD5 = "9f3c6733dc608d2bd37c488c560a5443"
ARCHIVE_SIZE_BYTES = 37931874
CHUNK_SIZE = 1024 * 1024


class HostGraphExportError(RuntimeError):
    pass


def digest_file(path: Path, algorithm: str) -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(CHUNK_SIZE), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def canonical_hash(value: Any) -> str:
    return sha256_bytes(
        json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
        ).encode("utf-8")
    )


def safe_member(name: str) -> PurePosixPath:
    if "\\" in name:
        raise HostGraphExportError(f"backslash in member path: {name!r}")
    path = PurePosixPath(name)
    if path.is_absolute() or any(part in {"", ".", ".."} for part in path.parts):
        raise HostGraphExportError(f"unsafe member path: {name!r}")
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
    return PurePosixPath(path).parent.name.lower()


def structural_node_payload(nodes: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        {
            "host_atom_index_1_based": row["host_atom_index_1_based"],
            "element": row["element"],
            "atomic_number": row["atomic_number"],
            "heavy_neighbour_host_indices_1_based": row[
                "heavy_neighbour_host_indices_1_based"
            ],
        }
        for row in nodes
    ]


def structural_edge_payload(edges: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        {
            "begin_host_atom_index_1_based": row[
                "begin_host_atom_index_1_based"
            ],
            "end_host_atom_index_1_based": row["end_host_atom_index_1_based"],
        }
        for row in edges
    ]


def compiled_node_payload(nodes: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        {
            "host_atom_index_1_based": row["host_atom_index_1_based"],
            "atom_name": row["atom_name"],
            "element": row["element"],
            "amber_atom_type": row["amber_atom_type"],
            "charge_e": row["charge_e"],
            "mass": row["mass"],
            "heavy_neighbour_host_indices_1_based": row[
                "heavy_neighbour_host_indices_1_based"
            ],
        }
        for row in nodes
    ]


def compiled_edge_payload(edges: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        {
            "begin_host_atom_index_1_based": row[
                "begin_host_atom_index_1_based"
            ],
            "end_host_atom_index_1_based": row["end_host_atom_index_1_based"],
            "force_constant": row["force_constant"],
            "equilibrium_A": row["equilibrium_A"],
        }
        for row in edges
    ]


def export_topology(path: str, payload: bytes) -> dict[str, Any]:
    try:
        text = payload.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise HostGraphExportError(f"{path} is not UTF-8 AMBER topology text") from exc

    sections = amber.parse_parm7_text(text)
    pointers = amber.pointer_map(sections)
    natom = pointers["NATOM"]
    atom_names = [
        str(value) for value in amber.required(sections, "ATOM_NAME", natom)[:natom]
    ]
    atom_types = [
        str(value)
        for value in amber.required(sections, "AMBER_ATOM_TYPE", natom)[:natom]
    ]
    atomic_numbers = [
        int(value)
        for value in amber.required(sections, "ATOMIC_NUMBER", natom)[:natom]
    ]
    masses = [
        float(value) for value in amber.required(sections, "MASS", natom)[:natom]
    ]
    raw_charges = [
        float(value) for value in amber.required(sections, "CHARGE", natom)[:natom]
    ]
    charges_e = [value / amber.AMBER_CHARGE_SCALE for value in raw_charges]

    residues = chemistry.residue_spans(sections, pointers)
    host = max(residues, key=lambda row: row["atom_count"])
    start = host["first_atom_1_based"]
    end = host["last_atom_1_based"]
    host_atoms = set(range(start, end + 1))
    heavy_atoms = {
        index
        for index in host_atoms
        if atomic_numbers[index - 1] != 1
    }

    bonds = chemistry.all_bonds(sections, pointers)
    adjacency: dict[int, set[int]] = {index: set() for index in host_atoms}
    for bond in bonds:
        left, right = bond["atoms_1_based"]
        if left in host_atoms and right in host_atoms:
            adjacency[left].add(right)
            adjacency[right].add(left)

    nodes: list[dict[str, Any]] = []
    for compiled_index in sorted(heavy_atoms):
        host_index = compiled_index - start + 1
        heavy_neighbours = sorted(
            neighbour - start + 1
            for neighbour in adjacency[compiled_index]
            if neighbour in heavy_atoms
        )
        nodes.append(
            {
                "compiled_atom_index_1_based": compiled_index,
                "host_atom_index_1_based": host_index,
                "atom_name": atom_names[compiled_index - 1],
                "element": chemistry.ELEMENT_SYMBOLS[atomic_numbers[compiled_index - 1]],
                "atomic_number": atomic_numbers[compiled_index - 1],
                "amber_atom_type": atom_types[compiled_index - 1],
                "charge_e": charges_e[compiled_index - 1],
                "mass": masses[compiled_index - 1],
                "total_degree": len(adjacency[compiled_index]),
                "heavy_degree": len(heavy_neighbours),
                "heavy_neighbour_host_indices_1_based": heavy_neighbours,
            }
        )

    edges: list[dict[str, Any]] = []
    for bond in bonds:
        left, right = bond["atoms_1_based"]
        if left not in heavy_atoms or right not in heavy_atoms:
            continue
        begin, finish = sorted((left, right))
        edges.append(
            {
                "begin_compiled_atom_index_1_based": begin,
                "end_compiled_atom_index_1_based": finish,
                "begin_host_atom_index_1_based": begin - start + 1,
                "end_host_atom_index_1_based": finish - start + 1,
                "force_constant": bond["force_constant"],
                "equilibrium_A": bond["equilibrium_A"],
                "source_flag": bond["source_flag"],
                "chemical_bond_order_known": False,
            }
        )
    edges.sort(
        key=lambda row: (
            row["begin_host_atom_index_1_based"],
            row["end_host_atom_index_1_based"],
        )
    )

    host_atomic_numbers = [atomic_numbers[index - 1] for index in sorted(host_atoms)]
    formula, element_counts = chemistry.formula_from_atomic_numbers(
        host_atomic_numbers
    )
    structural_payload = {
        "nodes": structural_node_payload(nodes),
        "edges": structural_edge_payload(edges),
    }
    compiled_payload = {
        "nodes": compiled_node_payload(nodes),
        "edges": compiled_edge_payload(edges),
    }

    return {
        "system_name": infer_system(path),
        "topology_path": path,
        "topology_sha256": sha256_bytes(payload),
        "host_residue": {
            "index_1_based": host["index_1_based"],
            "label": host["label"],
            "first_compiled_atom_index_1_based": start,
            "last_compiled_atom_index_1_based": end,
            "atom_count": len(host_atoms),
            "heavy_atom_count": len(heavy_atoms),
            "formula": formula,
            "element_counts": element_counts,
            "compiled_charge_e": sum(charges_e[index - 1] for index in host_atoms),
        },
        "heavy_graph": {
            "node_count": len(nodes),
            "edge_count": len(edges),
            "nodes": nodes,
            "edges": edges,
            "chemical_bond_order_known": False,
            "structural_fingerprint_sha256": canonical_hash(structural_payload),
            "compiled_parameter_fingerprint_sha256": canonical_hash(
                compiled_payload
            ),
        },
    }


def audit(archive_path: Path, output_json: Path) -> dict[str, Any]:
    if archive_path.stat().st_size != ARCHIVE_SIZE_BYTES:
        raise HostGraphExportError(
            f"size mismatch: expected {ARCHIVE_SIZE_BYTES}, "
            f"got {archive_path.stat().st_size}"
        )
    md5 = digest_file(archive_path, "md5")
    if md5 != ARCHIVE_MD5:
        raise HostGraphExportError(f"MD5 mismatch: expected {ARCHIVE_MD5}, got {md5}")
    archive_sha256 = digest_file(archive_path, "sha256")

    systems: list[dict[str, Any]] = []
    seen: set[str] = set()
    with zipfile.ZipFile(archive_path) as archive:
        bad = archive.testzip()
        if bad:
            raise HostGraphExportError(f"ZIP CRC failure: {bad}")
        for info in archive.infolist():
            safe_member(info.filename)
            if info.filename in seen:
                raise HostGraphExportError(f"duplicate archive member: {info.filename}")
            seen.add(info.filename)
            if info.flag_bits & 0x1:
                raise HostGraphExportError(f"encrypted member: {info.filename}")
            kind = member_kind(info)
            if kind in {"symlink", "special"}:
                raise HostGraphExportError(f"unsupported {kind}: {info.filename}")
            if kind != "file":
                continue
            suffix = PurePosixPath(info.filename).suffix.lower()
            if suffix in {".parm7", ".prmtop"}:
                systems.append(export_topology(info.filename, archive.read(info)))

    systems.sort(key=lambda row: row["system_name"])
    if {row["system_name"] for row in systems} != {
        "glucose",
        "galactose",
        "fructopyranose",
        "fructofuranose",
    }:
        raise HostGraphExportError(
            f"unexpected system set: {[row['system_name'] for row in systems]}"
        )

    structural_hashes = {
        row["heavy_graph"]["structural_fingerprint_sha256"] for row in systems
    }
    parameter_hashes = {
        row["heavy_graph"]["compiled_parameter_fingerprint_sha256"]
        for row in systems
    }
    formulas = {row["host_residue"]["formula"] for row in systems}
    node_counts = {row["heavy_graph"]["node_count"] for row in systems}
    edge_counts = {row["heavy_graph"]["edge_count"] for row in systems}
    if len(structural_hashes) != 1:
        raise HostGraphExportError("public host heavy-graph structure differs across sugars")
    if len(parameter_hashes) != 1:
        raise HostGraphExportError("public host compiled parameters differ across sugars")

    canonical = next(row for row in systems if row["system_name"] == "glucose")
    report = {
        "schema": "dia-health-one.gluhut-public-host-heavy-graph.v1",
        "status": "public_compiled_host_graph_exported_bond_orders_unknown",
        "scope": "public truncated GluHUT benchmark host only",
        "not_target_claim": (
            "not MC-P0, NNC2215, PaG, MC-Q1, affinity, efficacy, safety, or a drug result"
        ),
        "archive": {
            "filename": archive_path.name,
            "size_bytes": archive_path.stat().st_size,
            "md5": md5,
            "sha256": archive_sha256,
            "crc_verified": True,
        },
        "cross_system": {
            "system_count": len(systems),
            "systems": [row["system_name"] for row in systems],
            "host_formula_identical": len(formulas) == 1,
            "host_heavy_node_count_identical": len(node_counts) == 1,
            "host_heavy_edge_count_identical": len(edge_counts) == 1,
            "host_structural_graph_identical": len(structural_hashes) == 1,
            "host_compiled_parameters_identical": len(parameter_hashes) == 1,
            "structural_fingerprint_sha256": next(iter(structural_hashes)),
            "compiled_parameter_fingerprint_sha256": next(iter(parameter_hashes)),
        },
        "canonical_host": {
            "source_system": canonical["system_name"],
            "source_topology_path": canonical["topology_path"],
            "source_topology_sha256": canonical["topology_sha256"],
            "residue": canonical["host_residue"],
            "heavy_graph": canonical["heavy_graph"],
        },
        "system_sources": [
            {
                "system_name": row["system_name"],
                "topology_path": row["topology_path"],
                "topology_sha256": row["topology_sha256"],
                "structural_fingerprint_sha256": row["heavy_graph"][
                    "structural_fingerprint_sha256"
                ],
                "compiled_parameter_fingerprint_sha256": row["heavy_graph"][
                    "compiled_parameter_fingerprint_sha256"
                ],
            }
            for row in systems
        ],
        "decision": {
            "safe_for_public_constitutional_connectivity_comparison": True,
            "safe_for_parameter_transfer_to_private_candidates": False,
            "chemical_bond_orders_reconstructed": False,
            "candidate_ranking_authorized": False,
            "next_action": (
                "compare this element-connectivity graph privately with the frozen MC-P0 "
                "graph; treat every non-common atom, edge, charge, type and term as "
                "non-transferable until independently parameterized"
            ),
        },
        "limitations": [
            "AMBER parm7 connectivity does not preserve formal chemical bond order.",
            "Compiled charges and terms are internally identical but their RESP/GAFF source derivation is not deposited.",
            "The exported host is the public truncated benchmark receptor, not MC-P0 or complete NNC2215.",
            "A common element-connectivity subgraph does not authorize charge or force-field transfer.",
        ],
    }
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(
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
        default=Path("research/data/gluhut-public-host-heavy-graph.json"),
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    report = audit(args.archive, args.output_json)
    print(
        json.dumps(
            {
                "archive": report["archive"],
                "cross_system": report["cross_system"],
                "canonical_residue": report["canonical_host"]["residue"],
                "decision": report["decision"],
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
