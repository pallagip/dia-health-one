#!/usr/bin/env python3
"""Recover and audit recorded coordinates from the public GluHUT BINPOS files.

CPPTRAJ's Scripps BINPOS implementation defines a four-byte ``fxyz`` header,
followed by frame records containing one native-endian integer atom count and
``3*natom`` native-endian single-precision coordinates. This audit verifies and
parses every deposited trajectory, associates it with the matching AMBER
compiled topology and deposited group-COM restraint, and measures the recorded
host–guest reaction-coordinate path.

The first trajectory frame is a recorded simulation frame, not necessarily the
pre-production starting coordinate or the exact source restart. This audit does
not calculate work, PMFs, affinity, or candidate performance.
"""
from __future__ import annotations

import argparse
import collections
import hashlib
import json
import math
import stat
import struct
import sys
import zipfile
from pathlib import Path, PurePosixPath
from typing import Any, Iterable

sys.path.insert(0, str(Path(__file__).resolve().parent))
import audit_gluhut_compiled_systems as amber  # noqa: E402
import audit_gluhut_smd_restraints as restraints  # noqa: E402

ARCHIVE_MD5 = "9f3c6733dc608d2bd37c488c560a5443"
ARCHIVE_SIZE_BYTES = 37931874
CHUNK_SIZE = 1024 * 1024
BINPOS_MAGIC = b"fxyz"
EXPECTED_SYSTEMS = {
    "glucose",
    "galactose",
    "fructopyranose",
    "fructofuranose",
}


class BinposAuditError(RuntimeError):
    pass


def digest_file(path: Path, algorithm: str) -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(CHUNK_SIZE), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_bytes(payload: bytes | memoryview) -> str:
    return hashlib.sha256(payload).hexdigest()


def canonical_hash(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
        ).encode("utf-8")
    ).hexdigest()


def safe_member(name: str) -> PurePosixPath:
    if "\\" in name:
        raise BinposAuditError(f"backslash in member path: {name!r}")
    path = PurePosixPath(name)
    if path.is_absolute() or any(part in {"", ".", ".."} for part in path.parts):
        raise BinposAuditError(f"unsafe member path: {name!r}")
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
    for name in (
        "fructofuranose",
        "fructopyranose",
        "galactose",
        "glucose",
    ):
        if name in lower:
            return name
    raise BinposAuditError(f"could not infer sugar system from path: {path}")


def topology_record(path: str, payload: bytes) -> dict[str, Any]:
    sections = amber.parse_parm7_text(payload.decode("utf-8"))
    pointers = amber.pointer_map(sections)
    natom = pointers["NATOM"]
    nres = pointers["NRES"]
    atom_names = [
        str(value)
        for value in amber.required(sections, "ATOM_NAME", natom)[:natom]
    ]
    atom_types = [
        str(value)
        for value in amber.required(sections, "AMBER_ATOM_TYPE", natom)[:natom]
    ]
    masses = [
        float(value)
        for value in amber.required(sections, "MASS", natom)[:natom]
    ]
    atomic_numbers = [
        int(value)
        for value in amber.required(sections, "ATOMIC_NUMBER", natom)[:natom]
    ]
    residue_labels = [
        str(value)
        for value in amber.required(sections, "RESIDUE_LABEL", nres)[:nres]
    ]
    residue_starts = [
        int(value)
        for value in amber.required(sections, "RESIDUE_POINTER", nres)[:nres]
    ]
    atom_residue = [0] * (natom + 1)
    for index, start in enumerate(residue_starts, start=1):
        end = residue_starts[index] - 1 if index < len(residue_starts) else natom
        for atom_index in range(start, end + 1):
            atom_residue[atom_index] = index
    if any(index == 0 for index in atom_residue[1:]):
        raise BinposAuditError(f"unassigned residue atom in {path}")

    atoms = {
        atom_index: {
            "atom_1_based": atom_index,
            "atom_name": atom_names[atom_index - 1],
            "amber_atom_type": atom_types[atom_index - 1],
            "atomic_number": atomic_numbers[atom_index - 1],
            "mass": masses[atom_index - 1],
            "residue_index_1_based": atom_residue[atom_index],
            "residue_label": residue_labels[atom_residue[atom_index] - 1],
        }
        for atom_index in range(1, natom + 1)
    }
    return {
        "path": path,
        "sha256": sha256_bytes(payload),
        "system_name": infer_system(path),
        "natom": natom,
        "atoms": atoms,
        "masses": masses,
        "total_mass": sum(masses),
    }


def detect_endian(payload: bytes, expected_natom: int) -> str:
    if len(payload) < 8:
        raise BinposAuditError("BINPOS payload is too short")
    if payload[:4] != BINPOS_MAGIC:
        raise BinposAuditError(
            f"invalid BINPOS magic: expected {BINPOS_MAGIC!r}, got {payload[:4]!r}"
        )
    little = struct.unpack_from("<i", payload, 4)[0]
    big = struct.unpack_from(">i", payload, 4)[0]
    matches = []
    if little == expected_natom:
        matches.append("<")
    if big == expected_natom:
        matches.append(">")
    if len(matches) != 1:
        raise BinposAuditError(
            f"cannot determine BINPOS endian for expected natom {expected_natom}: "
            f"little={little}, big={big}"
        )
    return matches[0]


def binpos_layout(payload: bytes, expected_natom: int) -> dict[str, Any]:
    endian = detect_endian(payload, expected_natom)
    coordinate_count = 3 * expected_natom
    coordinate_bytes = coordinate_count * 4
    record_bytes = 4 + coordinate_bytes
    data_bytes = len(payload) - 4
    if data_bytes <= 0 or data_bytes % record_bytes != 0:
        raise BinposAuditError(
            f"BINPOS size mismatch: data_bytes={data_bytes}, "
            f"record_bytes={record_bytes}, remainder={data_bytes % record_bytes}"
        )
    frame_count = data_bytes // record_bytes
    if frame_count < 1:
        raise BinposAuditError("BINPOS contains no frames")
    return {
        "endian": endian,
        "natom": expected_natom,
        "coordinate_count_per_frame": coordinate_count,
        "coordinate_bytes_per_frame": coordinate_bytes,
        "record_bytes": record_bytes,
        "frame_count": frame_count,
    }


def frame_raw_view(
    payload: bytes,
    frame_index: int,
    layout: dict[str, Any],
) -> memoryview:
    if not 0 <= frame_index < layout["frame_count"]:
        raise IndexError(frame_index)
    record_offset = 4 + frame_index * layout["record_bytes"]
    natom = struct.unpack_from(layout["endian"] + "i", payload, record_offset)[0]
    if natom != layout["natom"]:
        raise BinposAuditError(
            f"frame {frame_index}: natom {natom} != expected {layout['natom']}"
        )
    start = record_offset + 4
    end = start + layout["coordinate_bytes_per_frame"]
    return memoryview(payload)[start:end]


def frame_coordinates(
    payload: bytes,
    frame_index: int,
    layout: dict[str, Any],
) -> tuple[float, ...]:
    raw = frame_raw_view(payload, frame_index, layout)
    return struct.unpack(
        layout["endian"]
        + f"{layout['coordinate_count_per_frame']}f",
        raw,
    )


def validate_finite_coordinates(
    coordinates: Iterable[float],
    path: str,
    frame_index: int,
) -> None:
    for coordinate_index, value in enumerate(coordinates):
        if not math.isfinite(value):
            raise BinposAuditError(
                f"{path} frame {frame_index} coordinate {coordinate_index} "
                f"is not finite: {value}"
            )


def centre_of_mass(
    coordinates: tuple[float, ...],
    atom_indices_1_based: list[int],
    masses: list[float],
) -> tuple[float, float, float]:
    total_mass = 0.0
    x = 0.0
    y = 0.0
    z = 0.0
    for atom_index in atom_indices_1_based:
        if not 1 <= atom_index <= len(masses):
            raise BinposAuditError(f"COM atom outside topology: {atom_index}")
        mass = masses[atom_index - 1]
        offset = 3 * (atom_index - 1)
        total_mass += mass
        x += mass * coordinates[offset]
        y += mass * coordinates[offset + 1]
        z += mass * coordinates[offset + 2]
    if total_mass <= 0.0:
        raise BinposAuditError("COM group has non-positive total mass")
    return x / total_mass, y / total_mass, z / total_mass


def distance(
    left: tuple[float, float, float],
    right: tuple[float, float, float],
) -> float:
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(left, right)))


def pearson_index_correlation(values: list[float]) -> float | None:
    count = len(values)
    if count < 2:
        return None
    mean_x = (count - 1) / 2.0
    mean_y = sum(values) / count
    var_x = sum((index - mean_x) ** 2 for index in range(count))
    var_y = sum((value - mean_y) ** 2 for value in values)
    if var_x <= 0.0 or var_y <= 0.0:
        return None
    covariance = sum(
        (index - mean_x) * (value - mean_y)
        for index, value in enumerate(values)
    )
    return covariance / math.sqrt(var_x * var_y)


def linear_schedule_rmse(
    values: list[float],
    start: float,
    end: float,
) -> float:
    if not values:
        raise ValueError("empty values")
    if len(values) == 1:
        expected = [start]
    else:
        expected = [
            start + (end - start) * index / (len(values) - 1)
            for index in range(len(values))
        ]
    return math.sqrt(
        sum((value - target) ** 2 for value, target in zip(values, expected))
        / len(values)
    )


def coordinate_bounds(coordinates: tuple[float, ...]) -> dict[str, list[float]]:
    xs = coordinates[0::3]
    ys = coordinates[1::3]
    zs = coordinates[2::3]
    return {
        "minimum_xyz_A": [min(xs), min(ys), min(zs)],
        "maximum_xyz_A": [max(xs), max(ys), max(zs)],
    }


def audit_trajectory(
    path: str,
    payload: bytes,
    topology: dict[str, Any],
    restraint: dict[str, Any],
) -> dict[str, Any]:
    layout = binpos_layout(payload, topology["natom"])
    group1 = restraint["igr1"]
    group2 = restraint["igr2"]
    distances: list[float] = []
    first_coordinates: tuple[float, ...] | None = None
    last_coordinates: tuple[float, ...] | None = None

    for frame_index in range(layout["frame_count"]):
        coordinates = frame_coordinates(payload, frame_index, layout)
        validate_finite_coordinates(coordinates, path, frame_index)
        if frame_index == 0:
            first_coordinates = coordinates
        if frame_index == layout["frame_count"] - 1:
            last_coordinates = coordinates
        host_com = centre_of_mass(coordinates, group1, topology["masses"])
        guest_com = centre_of_mass(coordinates, group2, topology["masses"])
        distances.append(distance(host_com, guest_com))

    assert first_coordinates is not None
    assert last_coordinates is not None
    deltas = [right - left for left, right in zip(distances, distances[1:])]
    monotonic_tolerance_A = 0.05
    nondecreasing_fraction = (
        sum(delta >= -monotonic_tolerance_A for delta in deltas) / len(deltas)
        if deltas
        else None
    )
    correlation = pearson_index_correlation(distances)
    first_raw = frame_raw_view(payload, 0, layout)
    last_raw = frame_raw_view(payload, layout["frame_count"] - 1, layout)

    expected_start = float(restraint["r2"])
    expected_end = float(restraint["r2a"])
    net_change = distances[-1] - distances[0]
    smd_like = (
        net_change > 5.0
        and correlation is not None
        and correlation > 0.8
    )

    return {
        "path": path,
        "sha256": sha256_bytes(payload),
        "size_bytes": len(payload),
        "system_name": topology["system_name"],
        "topology_path": topology["path"],
        "topology_sha256": topology["sha256"],
        "restraint_path": restraint["path"],
        "restraint_sha256": restraint["sha256"],
        "layout": layout,
        "first_frame": {
            "coordinate_bytes_sha256": sha256_bytes(first_raw),
            "reaction_coordinate_A": distances[0],
            **coordinate_bounds(first_coordinates),
        },
        "last_frame": {
            "coordinate_bytes_sha256": sha256_bytes(last_raw),
            "reaction_coordinate_A": distances[-1],
            **coordinate_bounds(last_coordinates),
        },
        "reaction_coordinate": {
            "group1_atom_indices_1_based": group1,
            "group2_atom_indices_1_based": group2,
            "group1_total_mass_Da": sum(
                topology["masses"][index - 1] for index in group1
            ),
            "group2_total_mass_Da": sum(
                topology["masses"][index - 1] for index in group2
            ),
            "minimum_A": min(distances),
            "minimum_frame_index_zero_based": distances.index(min(distances)),
            "maximum_A": max(distances),
            "maximum_frame_index_zero_based": distances.index(max(distances)),
            "first_A": distances[0],
            "last_A": distances[-1],
            "net_change_A": net_change,
            "mean_A": sum(distances) / len(distances),
            "pearson_correlation_with_frame_index": correlation,
            "fraction_nondecreasing_with_0_05A_tolerance": nondecreasing_fraction,
            "rmse_from_linear_0_to_15A_schedule_A": linear_schedule_rmse(
                distances,
                expected_start,
                expected_end,
            ),
            "deposited_target_start_A": expected_start,
            "deposited_target_end_A": expected_end,
            "smd_like_forward_unbinding_path": smd_like,
            "distance_series_sha256": canonical_hash(distances),
        },
    }


def audit(archive_path: Path, output_json: Path) -> dict[str, Any]:
    if archive_path.stat().st_size != ARCHIVE_SIZE_BYTES:
        raise BinposAuditError(
            f"size mismatch: expected {ARCHIVE_SIZE_BYTES}, "
            f"got {archive_path.stat().st_size}"
        )
    archive_md5 = digest_file(archive_path, "md5")
    if archive_md5 != ARCHIVE_MD5:
        raise BinposAuditError(
            f"MD5 mismatch: expected {ARCHIVE_MD5}, got {archive_md5}"
        )
    archive_sha256 = digest_file(archive_path, "sha256")

    topologies: dict[str, dict[str, Any]] = {}
    restraint_records: list[dict[str, Any]] = []
    trajectory_payloads: list[tuple[str, bytes]] = []
    seen: set[str] = set()

    with zipfile.ZipFile(archive_path) as archive:
        bad = archive.testzip()
        if bad:
            raise BinposAuditError(f"ZIP CRC failure: {bad}")
        for info in archive.infolist():
            safe_member(info.filename)
            if info.filename in seen:
                raise BinposAuditError(f"duplicate member: {info.filename}")
            seen.add(info.filename)
            if info.flag_bits & 0x1:
                raise BinposAuditError(f"encrypted member: {info.filename}")
            kind = member_kind(info)
            if kind in {"symlink", "special"}:
                raise BinposAuditError(f"unsupported {kind}: {info.filename}")
            if kind != "file":
                continue
            suffix = PurePosixPath(info.filename).suffix.lower()
            payload = archive.read(info)
            if suffix in {".parm7", ".prmtop"}:
                topology = topology_record(info.filename, payload)
                if topology["system_name"] in topologies:
                    raise BinposAuditError(
                        f"duplicate topology for {topology['system_name']}"
                    )
                topologies[topology["system_name"]] = topology
            elif suffix == ".rst":
                restraint_records.append(
                    restraints.parse_restraint(info.filename, payload)
                )
            elif suffix == ".binpos":
                trajectory_payloads.append((info.filename, payload))

    if set(topologies) != EXPECTED_SYSTEMS:
        raise BinposAuditError(
            f"unexpected topology systems: {sorted(topologies)}"
        )
    restraint_by_system: dict[str, dict[str, Any]] = {}
    for restraint in restraint_records:
        for system_name in restraint["applies_to_systems"]:
            if system_name in restraint_by_system:
                raise BinposAuditError(f"duplicate restraint for {system_name}")
            restraint_by_system[system_name] = restraint
    if set(restraint_by_system) != EXPECTED_SYSTEMS:
        raise BinposAuditError(
            f"unexpected restraint systems: {sorted(restraint_by_system)}"
        )
    if len(trajectory_payloads) != 20:
        raise BinposAuditError(
            f"expected 20 BINPOS trajectories, found {len(trajectory_payloads)}"
        )

    trajectory_rows = []
    for path, payload in sorted(trajectory_payloads):
        system_name = infer_system(path)
        trajectory_rows.append(
            audit_trajectory(
                path,
                payload,
                topologies[system_name],
                restraint_by_system[system_name],
            )
        )

    by_system: dict[str, list[dict[str, Any]]] = collections.defaultdict(list)
    for row in trajectory_rows:
        by_system[row["system_name"]].append(row)
    system_summaries = []
    for system_name in sorted(by_system):
        rows = by_system[system_name]
        first_values = [
            row["reaction_coordinate"]["first_A"] for row in rows
        ]
        last_values = [
            row["reaction_coordinate"]["last_A"] for row in rows
        ]
        frame_counts = [row["layout"]["frame_count"] for row in rows]
        system_summaries.append(
            {
                "system_name": system_name,
                "trajectory_count": len(rows),
                "trajectory_paths": [row["path"] for row in rows],
                "frame_counts": frame_counts,
                "total_recorded_frames": sum(frame_counts),
                "first_reaction_coordinate_range_A": [
                    min(first_values),
                    max(first_values),
                ],
                "last_reaction_coordinate_range_A": [
                    min(last_values),
                    max(last_values),
                ],
                "smd_like_forward_path_count": sum(
                    row["reaction_coordinate"][
                        "smd_like_forward_unbinding_path"
                    ]
                    for row in rows
                ),
                "unique_first_frame_coordinate_hash_count": len(
                    {
                        row["first_frame"]["coordinate_bytes_sha256"]
                        for row in rows
                    }
                ),
                "unique_last_frame_coordinate_hash_count": len(
                    {
                        row["last_frame"]["coordinate_bytes_sha256"]
                        for row in rows
                    }
                ),
            }
        )

    all_headers_valid = all(
        row["layout"]["natom"] == topologies[row["system_name"]]["natom"]
        for row in trajectory_rows
    )
    all_little_endian = all(
        row["layout"]["endian"] == "<" for row in trajectory_rows
    )
    all_smd_like = all(
        row["reaction_coordinate"]["smd_like_forward_unbinding_path"]
        for row in trajectory_rows
    )
    total_frames = sum(row["layout"]["frame_count"] for row in trajectory_rows)

    report = {
        "schema": "dia-health-one.gluhut-binpos-coordinate-recovery.v1",
        "status": "recorded_trajectory_coordinates_recovered_not_original_restart",
        "scope": "twenty public GluHUT Scripps BINPOS trajectories",
        "not_target_claim": (
            "not source restart recovery, not SMD work or PMF, not NNC2215, PaG, "
            "MC-Q1, affinity, efficacy, safety, or a drug result"
        ),
        "format_source": {
            "repository": "Amber-MD/cpptraj",
            "path": "src/Traj_Binpos.cpp",
            "format": (
                "four-byte fxyz header; each frame is one native-endian int natom "
                "followed by 3*natom native-endian float coordinates"
            ),
        },
        "archive": {
            "filename": archive_path.name,
            "size_bytes": archive_path.stat().st_size,
            "md5": archive_md5,
            "sha256": archive_sha256,
            "crc_verified": True,
        },
        "topologies": {
            system_name: {
                "path": topology["path"],
                "sha256": topology["sha256"],
                "natom": topology["natom"],
                "total_mass": topology["total_mass"],
            }
            for system_name, topology in sorted(topologies.items())
        },
        "trajectory_count": len(trajectory_rows),
        "total_recorded_frames": total_frames,
        "trajectories": trajectory_rows,
        "system_summaries": system_summaries,
        "decision": {
            "all_binpos_headers_and_record_sizes_valid": all_headers_valid,
            "all_files_little_endian": all_little_endian,
            "recorded_coordinate_frames_recoverable": True,
            "recorded_first_and_last_frames_hashed": True,
            "all_paths_match_forward_unbinding_shape": all_smd_like,
            "original_pretrajectory_start_coordinates_recovered": False,
            "recorded_frames_can_support_coordinate_level_provenance_audit": True,
            "recorded_frames_can_replace_missing_source_restart": False,
            "new_smd_or_umbrella_generation_authorized": False,
            "parent_pag_calibration_authorized": False,
            "mcq1_parameter_transfer_authorized": False,
            "candidate_ranking_authorized": False,
            "next_action": (
                "use recorded first/last frames to audit deposited paths and reconstruct "
                "representative coordinate snapshots; still request exact source restarts "
                "and corrected guest groups before prospective generation"
            ),
        },
        "limitations": [
            "The first recorded frame may follow minimisation, equilibration, or earlier unrecorded integration steps.",
            "BINPOS does not store periodic box information, velocities, forces, or simulation time.",
            "Reaction coordinates use the deposited atom groups, including the unresolved fructofuranose H62 selection.",
            "A monotonic recorded COM path does not validate the work protocol, force field, or free-energy calculation.",
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
        default=Path("research/data/gluhut-binpos-coordinate-recovery.json"),
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    report = audit(args.archive, args.output_json)
    print(
        json.dumps(
            {
                "trajectory_count": report["trajectory_count"],
                "total_recorded_frames": report["total_recorded_frames"],
                "system_summaries": report["system_summaries"],
                "decision": report["decision"],
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
