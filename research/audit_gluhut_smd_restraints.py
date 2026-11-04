#!/usr/bin/env python3
"""Map deposited GluHUT SMD restraint groups onto compiled AMBER atoms.

The audit verifies the public SMD archive, parses the three AMBER NMR-style
restraint files, maps every igr1/igr2 atom to the matching topology, and records
the induced graph and restraint parameters. It does not convert the restraint
into a binding free energy, validate the force-constant convention, or authorize
PaG/MC-Q1 ranking.
"""
from __future__ import annotations

import argparse
import collections
import hashlib
import json
import re
import stat
import zipfile
from pathlib import Path, PurePosixPath
from typing import Any

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import audit_gluhut_compiled_systems as amber  # noqa: E402

ARCHIVE_MD5 = "9f3c6733dc608d2bd37c488c560a5443"
ARCHIVE_SIZE_BYTES = 37931874
CHUNK_SIZE = 1024 * 1024

ELEMENT_SYMBOLS = {
    1: "H", 5: "B", 6: "C", 7: "N", 8: "O", 9: "F",
    15: "P", 16: "S", 17: "Cl", 35: "Br", 53: "I",
}

RESTRAINT_TO_SYSTEMS = {
    "glucose-and-galactose-smd.rst": ["glucose", "galactose"],
    "fructopyranose-smd.rst": ["fructopyranose"],
    "fructofuranose-smd.rst": ["fructofuranose"],
}


class RestraintAuditError(RuntimeError):
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
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def safe_member(name: str) -> PurePosixPath:
    if "\\" in name:
        raise RestraintAuditError(f"backslash in member path: {name!r}")
    path = PurePosixPath(name)
    if path.is_absolute() or any(part in {"", ".", ".."} for part in path.parts):
        raise RestraintAuditError(f"unsafe member path: {name!r}")
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


def strip_comments(text: str) -> str:
    lines: list[str] = []
    for raw in text.splitlines():
        line = raw.split("!", 1)[0].strip()
        if not line or line.startswith("&") or line == "/":
            continue
        lines.append(line)
    return " ".join(lines)


def parse_namelist_assignments(text: str) -> dict[str, list[str]]:
    """Parse the simple scalar/list assignments used in the deposited files."""
    body = strip_comments(text)
    matches = list(re.finditer(r"(?i)\b([a-z][a-z0-9_]*)\s*=", body))
    if not matches:
        raise RestraintAuditError("no namelist assignments found")
    result: dict[str, list[str]] = {}
    for index, match in enumerate(matches):
        key = match.group(1).lower()
        start = match.end()
        end = matches[index + 1].start() if index + 1 < len(matches) else len(body)
        raw_value = body[start:end].strip().strip(",").strip()
        tokens = [token.strip().strip(",") for token in raw_value.split(",")]
        tokens = [token for token in tokens if token]
        if key in result:
            raise RestraintAuditError(f"duplicate assignment {key}")
        result[key] = tokens
    return result


def as_int_list(assignments: dict[str, list[str]], key: str) -> list[int]:
    if key not in assignments:
        raise RestraintAuditError(f"missing {key}")
    try:
        return [int(token) for token in assignments[key]]
    except ValueError as exc:
        raise RestraintAuditError(f"non-integer value in {key}: {assignments[key]}") from exc


def as_float(assignments: dict[str, list[str]], key: str) -> float | None:
    if key not in assignments:
        return None
    if len(assignments[key]) != 1:
        raise RestraintAuditError(f"expected scalar {key}, got {assignments[key]}")
    return float(assignments[key][0].replace("D", "E").replace("d", "e"))


def infer_system(path: str) -> str:
    lower = path.lower()
    for name in ("fructofuranose", "fructopyranose", "galactose", "glucose"):
        if name in lower:
            return name
    return PurePosixPath(path).parent.name.lower()


def topology_atoms(path: str, payload: bytes) -> dict[str, Any]:
    sections = amber.parse_parm7_text(payload.decode("utf-8"))
    pointers = amber.pointer_map(sections)
    natom = pointers["NATOM"]
    nres = pointers["NRES"]
    names = [str(v) for v in amber.required(sections, "ATOM_NAME", natom)[:natom]]
    types = [str(v) for v in amber.required(sections, "AMBER_ATOM_TYPE", natom)[:natom]]
    atomic_numbers = [
        int(v) for v in amber.required(sections, "ATOMIC_NUMBER", natom)[:natom]
    ]
    raw_charges = [float(v) for v in amber.required(sections, "CHARGE", natom)[:natom]]
    charges = [value / amber.AMBER_CHARGE_SCALE for value in raw_charges]
    residue_labels = [
        str(v) for v in amber.required(sections, "RESIDUE_LABEL", nres)[:nres]
    ]
    residue_starts = [
        int(v) for v in amber.required(sections, "RESIDUE_POINTER", nres)[:nres]
    ]
    atom_residue = [0] * (natom + 1)
    for idx, start in enumerate(residue_starts, start=1):
        end = residue_starts[idx] - 1 if idx < len(residue_starts) else natom
        for atom in range(start, end + 1):
            atom_residue[atom] = idx

    bonds: list[tuple[int, int]] = []
    for flag in ("BONDS_INC_HYDROGEN", "BONDS_WITHOUT_HYDROGEN"):
        values = [int(v) for v in sections.get(flag, [])]
        if len(values) % 3:
            raise RestraintAuditError(f"{path}: malformed {flag}")
        for offset in range(0, len(values), 3):
            left = abs(values[offset]) // 3 + 1
            right = abs(values[offset + 1]) // 3 + 1
            bonds.append(tuple(sorted((left, right))))
    bonds = sorted(set(bonds))

    atoms = {}
    for atom in range(1, natom + 1):
        residue_index = atom_residue[atom]
        atoms[atom] = {
            "atom_1_based": atom,
            "atom_name": names[atom - 1],
            "amber_atom_type": types[atom - 1],
            "atomic_number": atomic_numbers[atom - 1],
            "element": ELEMENT_SYMBOLS.get(atomic_numbers[atom - 1], str(atomic_numbers[atom - 1])),
            "charge_e": charges[atom - 1],
            "residue_index_1_based": residue_index,
            "residue_label": residue_labels[residue_index - 1],
        }
    return {
        "path": path,
        "sha256": sha256_bytes(payload),
        "system_name": infer_system(path),
        "natom": natom,
        "atoms": atoms,
        "bonds": bonds,
    }


def induced_group(indices: list[int], topology: dict[str, Any]) -> dict[str, Any]:
    if len(indices) != len(set(indices)):
        raise RestraintAuditError(f"duplicate atom index in group: {indices}")
    natom = topology["natom"]
    if any(index < 1 or index > natom for index in indices):
        raise RestraintAuditError(f"group atom outside 1..{natom}: {indices}")
    selected = set(indices)
    atoms = [topology["atoms"][index] for index in indices]
    edges = [list(pair) for pair in topology["bonds"] if set(pair).issubset(selected)]
    adjacency = {index: set() for index in indices}
    for left, right in edges:
        adjacency[left].add(right)
        adjacency[right].add(left)
    components: list[list[int]] = []
    remaining = set(indices)
    while remaining:
        start = min(remaining)
        remaining.remove(start)
        stack = [start]
        component = []
        while stack:
            atom = stack.pop()
            component.append(atom)
            for neighbour in sorted(adjacency[atom], reverse=True):
                if neighbour in remaining:
                    remaining.remove(neighbour)
                    stack.append(neighbour)
        components.append(sorted(component))
    components.sort(key=lambda row: (len(row), row), reverse=True)
    element_counts = collections.Counter(atom["element"] for atom in atoms)
    residue_counts = collections.Counter(atom["residue_label"] for atom in atoms)
    type_counts = collections.Counter(atom["amber_atom_type"] for atom in atoms)
    degree_counts = collections.Counter(len(adjacency[index]) for index in indices)
    cycle_rank = len(edges) - len(indices) + len(components)
    return {
        "atom_indices_1_based": indices,
        "atom_count": len(indices),
        "atoms": atoms,
        "element_counts": dict(sorted(element_counts.items())),
        "residue_counts": dict(sorted(residue_counts.items())),
        "amber_atom_type_counts": dict(sorted(type_counts.items())),
        "induced_edges": edges,
        "induced_edge_count": len(edges),
        "connected_component_count": len(components),
        "component_atom_indices": components,
        "induced_degree_counts": {str(k): v for k, v in sorted(degree_counts.items())},
        "induced_cycle_rank": cycle_rank,
        "atom_index_set_sha256": canonical_hash(sorted(indices)),
        "mapped_atom_identity_sha256": canonical_hash(
            [
                {
                    "atom_1_based": atom["atom_1_based"],
                    "atom_name": atom["atom_name"],
                    "element": atom["element"],
                    "amber_atom_type": atom["amber_atom_type"],
                    "residue_label": atom["residue_label"],
                }
                for atom in atoms
            ]
        ),
    }


def parse_restraint(path: str, payload: bytes) -> dict[str, Any]:
    text = payload.decode("utf-8")
    assignments = parse_namelist_assignments(text)
    iat = as_int_list(assignments, "iat")
    if iat != [-1, -1]:
        raise RestraintAuditError(f"{path}: expected group-COM iat=-1,-1, got {iat}")
    igr1 = as_int_list(assignments, "igr1")
    igr2 = as_int_list(assignments, "igr2")
    basename = PurePosixPath(path).name.lower()
    if basename not in RESTRAINT_TO_SYSTEMS:
        raise RestraintAuditError(f"unrecognized restraint filename: {path}")
    return {
        "path": path,
        "sha256": sha256_bytes(payload),
        "applies_to_systems": RESTRAINT_TO_SYSTEMS[basename],
        "iat": iat,
        "igr1": igr1,
        "igr2": igr2,
        "r1": as_float(assignments, "r1"),
        "r2": as_float(assignments, "r2"),
        "r3": as_float(assignments, "r3"),
        "r4": as_float(assignments, "r4"),
        "rk2": as_float(assignments, "rk2"),
        "rk3": as_float(assignments, "rk3"),
        "r1a": as_float(assignments, "r1a"),
        "r2a": as_float(assignments, "r2a"),
        "r3a": as_float(assignments, "r3a"),
        "r4a": as_float(assignments, "r4a"),
        "rk2a": as_float(assignments, "rk2a"),
        "rk3a": as_float(assignments, "rk3a"),
        "nstep1": as_float(assignments, "nstep1"),
        "nstep2": as_float(assignments, "nstep2"),
        "iresid": as_float(assignments, "iresid"),
        "outxyz": as_float(assignments, "outxyz"),
        "raw_assignments": assignments,
    }


def audit(archive_path: Path, output_json: Path) -> dict[str, Any]:
    if archive_path.stat().st_size != ARCHIVE_SIZE_BYTES:
        raise RestraintAuditError(
            f"size mismatch: expected {ARCHIVE_SIZE_BYTES}, got {archive_path.stat().st_size}"
        )
    md5 = digest_file(archive_path, "md5")
    if md5 != ARCHIVE_MD5:
        raise RestraintAuditError(f"MD5 mismatch: expected {ARCHIVE_MD5}, got {md5}")
    archive_sha256 = digest_file(archive_path, "sha256")

    topologies: dict[str, dict[str, Any]] = {}
    restraints: list[dict[str, Any]] = []
    seen: set[str] = set()
    with zipfile.ZipFile(archive_path) as archive:
        bad = archive.testzip()
        if bad:
            raise RestraintAuditError(f"ZIP CRC failure: {bad}")
        for info in archive.infolist():
            safe_member(info.filename)
            if info.filename in seen:
                raise RestraintAuditError(f"duplicate member: {info.filename}")
            seen.add(info.filename)
            if info.flag_bits & 0x1:
                raise RestraintAuditError(f"encrypted member: {info.filename}")
            kind = member_kind(info)
            if kind in {"symlink", "special"}:
                raise RestraintAuditError(f"unsupported {kind}: {info.filename}")
            if kind != "file":
                continue
            suffix = PurePosixPath(info.filename).suffix.lower()
            payload = archive.read(info)
            if suffix in {".parm7", ".prmtop"}:
                topology = topology_atoms(info.filename, payload)
                topologies[topology["system_name"]] = topology
            elif suffix == ".rst":
                restraints.append(parse_restraint(info.filename, payload))

    if set(topologies) != {"glucose", "galactose", "fructopyranose", "fructofuranose"}:
        raise RestraintAuditError(f"unexpected topology systems: {sorted(topologies)}")
    if len(restraints) != 3:
        raise RestraintAuditError(f"expected 3 restraint files, found {len(restraints)}")

    mappings: list[dict[str, Any]] = []
    for restraint in sorted(restraints, key=lambda row: row["path"]):
        for system_name in restraint["applies_to_systems"]:
            topology = topologies[system_name]
            group1 = induced_group(restraint["igr1"], topology)
            group2 = induced_group(restraint["igr2"], topology)
            mappings.append(
                {
                    "restraint_path": restraint["path"],
                    "restraint_sha256": restraint["sha256"],
                    "system_name": system_name,
                    "topology_path": topology["path"],
                    "topology_sha256": topology["sha256"],
                    "parameters": {
                        key: restraint[key]
                        for key in (
                            "iat", "r1", "r2", "r3", "r4", "rk2", "rk3",
                            "r1a", "r2a", "r3a", "r4a", "rk2a", "rk3a",
                            "nstep1", "nstep2", "iresid", "outxyz",
                        )
                    },
                    "group1": group1,
                    "group2": group2,
                }
            )

    host_group_hashes = {row["group1"]["atom_index_set_sha256"] for row in mappings}
    host_identity_hashes = {row["group1"]["mapped_atom_identity_sha256"] for row in mappings}
    host_groups_all_carbon = all(row["group1"]["element_counts"] == {"C": 12} for row in mappings)
    host_groups_all_mol = all(row["group1"]["residue_counts"] == {"MOL": 12} for row in mappings)
    guest_groups_all_guest = all("MOL" not in row["group2"]["residue_counts"] for row in mappings)
    guest_groups_heavy_only = all("H" not in row["group2"]["element_counts"] for row in mappings)
    parameter_signatures = {
        canonical_hash(row["parameters"]) for row in mappings
    }

    group2_summary = [
        {
            "system_name": row["system_name"],
            "atom_indices_1_based": row["group2"]["atom_indices_1_based"],
            "atom_count": row["group2"]["atom_count"],
            "element_counts": row["group2"]["element_counts"],
            "residue_counts": row["group2"]["residue_counts"],
            "amber_atom_type_counts": row["group2"]["amber_atom_type_counts"],
            "induced_edge_count": row["group2"]["induced_edge_count"],
            "connected_component_count": row["group2"]["connected_component_count"],
            "induced_cycle_rank": row["group2"]["induced_cycle_rank"],
            "atom_names": [atom["atom_name"] for atom in row["group2"]["atoms"]],
        }
        for row in mappings
    ]

    report = {
        "schema": "dia-health-one.gluhut-smd-restraint-mapping.v1",
        "status": "restraint_atom_mapping_complete_force_convention_not_validated",
        "scope": "public GluHUT SMD group-COM restraint definitions",
        "not_target_claim": "not NNC2215, PaG, MC-Q1, affinity, efficacy, safety, or a drug result",
        "archive": {
            "filename": archive_path.name,
            "size_bytes": archive_path.stat().st_size,
            "md5": md5,
            "sha256": archive_sha256,
            "crc_verified": True,
        },
        "restraint_files": sorted(restraints, key=lambda row: row["path"]),
        "mappings": mappings,
        "cross_system_summary": {
            "mapping_count": len(mappings),
            "host_group_index_set_identical": len(host_group_hashes) == 1,
            "host_group_mapped_identity_identical": len(host_identity_hashes) == 1,
            "host_group_all_carbon": host_groups_all_carbon,
            "host_group_all_MOL_residue": host_groups_all_mol,
            "host_group_atom_indices_1_based": mappings[0]["group1"]["atom_indices_1_based"],
            "host_group_atom_names": [atom["atom_name"] for atom in mappings[0]["group1"]["atoms"]],
            "host_group_amber_types": [atom["amber_atom_type"] for atom in mappings[0]["group1"]["atoms"]],
            "guest_groups_exclude_host_residue": guest_groups_all_guest,
            "guest_groups_heavy_atoms_only": guest_groups_heavy_only,
            "restraint_parameters_identical_across_mappings": len(parameter_signatures) == 1,
            "guest_group_summary": group2_summary,
        },
        "decision": {
            "group1_is_same_12_host_carbons_across_all_systems": (
                len(host_group_hashes) == 1 and len(host_identity_hashes) == 1
                and host_groups_all_carbon and host_groups_all_mol
            ),
            "group2_is_guest_only_heavy_atom_selection": (
                guest_groups_all_guest and guest_groups_heavy_only
            ),
            "force_constant_convention_validated": False,
            "reaction_coordinate_matches_paper_description": None,
            "candidate_ranking_authorized": False,
            "pag_parameterization_authorized": False,
            "mcq1_parameter_transfer_authorized": False,
            "next_action": (
                "compare the exact mapped group atoms and rk2/r2a convention with the paper and "
                "author response; establish whether group2 is a ring COM or includes an exocyclic "
                "atom before reproducing SMD or umbrella generation"
            ),
        },
        "limitations": [
            "The audit maps atom indices and compiled connectivity but does not reconstruct coordinates or trajectories.",
            "Group membership does not by itself establish mass weighting, force convention, or the exact paper reporting convention.",
            "The fructopyranose group has a different atom count from the other guest groups and requires chemical interpretation.",
            "No SMD work value, PMF, affinity, kinetics, or candidate result is calculated.",
        ],
    }
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument(
        "--output-json",
        type=Path,
        default=Path("research/data/gluhut-smd-restraint-mapping.json"),
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    report = audit(args.archive, args.output_json)
    print(json.dumps({
        "cross_system_summary": report["cross_system_summary"],
        "decision": report["decision"],
        "mappings": report["mappings"],
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
