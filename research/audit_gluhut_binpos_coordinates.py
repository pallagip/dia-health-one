#!/usr/bin/env python3
"""Reconstruct deposited GluHUT SMD group-COM coordinates from BINPOS files.

The audit reads all public SMD trajectories with MDTraj's low-level BINPOS
reader, maps the deposited AMBER group restraints onto each topology, and
recalculates the mass-weighted host-to-guest centre-of-mass distance.

For fructofuranose it also evaluates the bounded graph-coherent alternatives to
selected hydrogen H62 on the already deposited coordinates. This does not
repair the original restraint, regenerate trajectories, calculate work or PMFs,
or authorize PaG/MC-Q1 ranking.
"""
from __future__ import annotations

import argparse
import collections
import hashlib
import json
import math
import shutil
import stat
import tempfile
import zipfile
from pathlib import Path, PurePosixPath
from typing import Any

import numpy as np

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import audit_gluhut_compiled_systems as amber  # noqa: E402
import audit_gluhut_smd_restraints as restraints  # noqa: E402

ARCHIVE_MD5 = "9f3c6733dc608d2bd37c488c560a5443"
ARCHIVE_SIZE_BYTES = 37931874
CHUNK_SIZE = 1024 * 1024
EXPECTED_NATOM = 183
EXPECTED_TRAJECTORY_COUNT = 20
TARGET_FINAL_A = 15.0


class CoordinateAuditError(RuntimeError):
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
        json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
        ).encode("utf-8")
    ).hexdigest()


def safe_member(name: str) -> PurePosixPath:
    if "\\" in name:
        raise CoordinateAuditError(f"backslash in member path: {name!r}")
    path = PurePosixPath(name)
    if path.is_absolute() or any(part in {"", ".", ".."} for part in path.parts):
        raise CoordinateAuditError(f"unsafe member path: {name!r}")
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
    raise CoordinateAuditError(f"cannot infer sugar system from path: {path}")


def topology_payload(path: str, payload: bytes) -> dict[str, Any]:
    parsed = restraints.topology_atoms(path, payload)
    sections = amber.parse_parm7_text(payload.decode("utf-8"))
    pointers = amber.pointer_map(sections)
    natom = pointers["NATOM"]
    masses = [float(v) for v in amber.required(sections, "MASS", natom)[:natom]]
    if natom != EXPECTED_NATOM:
        raise CoordinateAuditError(f"{path}: expected {EXPECTED_NATOM} atoms, found {natom}")
    parsed["masses"] = masses
    return parsed


def mass_weighted_com(
    xyz_A: np.ndarray,
    atom_indices_1_based: list[int],
    masses: list[float],
) -> np.ndarray:
    indices = np.asarray(atom_indices_1_based, dtype=int) - 1
    if indices.size == 0:
        raise CoordinateAuditError("empty COM atom selection")
    if np.any(indices < 0) or np.any(indices >= xyz_A.shape[1]):
        raise CoordinateAuditError("COM atom selection outside coordinate range")
    weights = np.asarray([masses[index] for index in indices], dtype=float)
    if np.any(weights <= 0.0) or not np.isfinite(weights).all():
        raise CoordinateAuditError("invalid masses in COM selection")
    return np.einsum("fai,a->fi", xyz_A[:, indices, :], weights) / weights.sum()


def reaction_coordinate_A(
    xyz_A: np.ndarray,
    host_indices: list[int],
    guest_indices: list[int],
    masses: list[float],
) -> np.ndarray:
    host = mass_weighted_com(xyz_A, host_indices, masses)
    guest = mass_weighted_com(xyz_A, guest_indices, masses)
    return np.linalg.norm(guest - host, axis=1)


def percentile_summary(values: np.ndarray) -> dict[str, float]:
    if values.ndim != 1 or values.size == 0:
        raise CoordinateAuditError("expected non-empty one-dimensional values")
    q = np.quantile(values, [0.0, 0.05, 0.25, 0.5, 0.75, 0.95, 1.0])
    return {
        "minimum": float(q[0]),
        "p05": float(q[1]),
        "p25": float(q[2]),
        "median": float(q[3]),
        "p75": float(q[4]),
        "p95": float(q[5]),
        "maximum": float(q[6]),
        "mean": float(values.mean()),
        "standard_deviation": float(values.std(ddof=0)),
    }


def trajectory_summary(coordinate: np.ndarray) -> dict[str, Any]:
    differences = np.diff(coordinate)
    return {
        "frame_count": int(coordinate.size),
        "first_A": float(coordinate[0]),
        "last_A": float(coordinate[-1]),
        "net_change_A": float(coordinate[-1] - coordinate[0]),
        "last_minus_target_15A": float(coordinate[-1] - TARGET_FINAL_A),
        "fraction_positive_frame_steps": (
            float(np.mean(differences > 0.0)) if differences.size else None
        ),
        "fraction_nonnegative_frame_steps": (
            float(np.mean(differences >= 0.0)) if differences.size else None
        ),
        "summary_A": percentile_summary(coordinate),
        "coordinate_sha256_float64": hashlib.sha256(
            np.asarray(coordinate, dtype="<f8").tobytes()
        ).hexdigest(),
    }


def first_frame_bond_statistics(
    xyz_A: np.ndarray,
    bonds: list[tuple[int, int]],
) -> dict[str, float]:
    first = xyz_A[0]
    distances = np.asarray(
        [
            np.linalg.norm(first[left - 1] - first[right - 1])
            for left, right in bonds
        ],
        dtype=float,
    )
    if distances.size == 0:
        raise CoordinateAuditError("topology has no bonds")
    return {
        "minimum_A": float(distances.min()),
        "median_A": float(np.median(distances)),
        "p95_A": float(np.quantile(distances, 0.95)),
        "maximum_A": float(distances.max()),
    }


def delta_summary(alternative: np.ndarray, deposited: np.ndarray) -> dict[str, Any]:
    if alternative.shape != deposited.shape:
        raise CoordinateAuditError("coordinate arrays have different shapes")
    delta = alternative - deposited
    return {
        "first_delta_A": float(delta[0]),
        "last_delta_A": float(delta[-1]),
        "mean_signed_delta_A": float(delta.mean()),
        "mean_absolute_delta_A": float(np.abs(delta).mean()),
        "rms_delta_A": float(np.sqrt(np.mean(delta**2))),
        "maximum_absolute_delta_A": float(np.abs(delta).max()),
        "delta_summary_A": percentile_summary(delta),
    }


def read_binpos(path: Path) -> tuple[np.ndarray, str]:
    try:
        import mdtraj
        from mdtraj.formats import BINPOSTrajectoryFile
    except ImportError as exc:
        raise CoordinateAuditError(
            "MDTraj is required; install the pinned workflow dependency"
        ) from exc
    with BINPOSTrajectoryFile(str(path), mode="r") as handle:
        xyz = handle.read()
    xyz = np.asarray(xyz, dtype=float)
    if xyz.ndim != 3 or xyz.shape[2] != 3:
        raise CoordinateAuditError(f"unexpected BINPOS shape {xyz.shape} for {path}")
    return xyz, mdtraj.__version__


def choose_restraint_for_system(
    system_name: str,
    parsed_restraints: list[dict[str, Any]],
) -> dict[str, Any]:
    matches = [
        row for row in parsed_restraints if system_name in row["applies_to_systems"]
    ]
    if len(matches) != 1:
        raise CoordinateAuditError(
            f"{system_name}: expected one restraint, found {len(matches)}"
        )
    return matches[0]


def alternative_map(candidate_record: dict[str, Any]) -> dict[str, list[int]]:
    result: dict[str, list[int]] = {}
    for row in candidate_record.get("viable_all_heavy_connected_hypotheses", []):
        label = str(row["label"])
        indices = [int(value) for value in row["atom_indices_1_based"]]
        if label in result:
            raise CoordinateAuditError(f"duplicate alternative label: {label}")
        result[label] = indices
    if not result:
        raise CoordinateAuditError("no viable fructofuranose alternatives found")
    return result


def audit(
    archive_path: Path,
    candidates_path: Path,
    output_path: Path,
) -> dict[str, Any]:
    if archive_path.stat().st_size != ARCHIVE_SIZE_BYTES:
        raise CoordinateAuditError(
            f"size mismatch: expected {ARCHIVE_SIZE_BYTES}, got {archive_path.stat().st_size}"
        )
    md5 = digest_file(archive_path, "md5")
    if md5 != ARCHIVE_MD5:
        raise CoordinateAuditError(f"MD5 mismatch: expected {ARCHIVE_MD5}, got {md5}")
    archive_sha256 = digest_file(archive_path, "sha256")
    candidate_record = json.loads(candidates_path.read_text(encoding="utf-8"))
    furanose_alternatives = alternative_map(candidate_record)

    topologies: dict[str, dict[str, Any]] = {}
    parsed_restraints: list[dict[str, Any]] = []
    trajectory_members: list[tuple[str, bytes]] = []
    seen: set[str] = set()

    with zipfile.ZipFile(archive_path) as archive:
        bad = archive.testzip()
        if bad:
            raise CoordinateAuditError(f"ZIP CRC failure: {bad}")
        for info in archive.infolist():
            safe_member(info.filename)
            if info.filename in seen:
                raise CoordinateAuditError(f"duplicate member: {info.filename}")
            seen.add(info.filename)
            if info.flag_bits & 0x1:
                raise CoordinateAuditError(f"encrypted member: {info.filename}")
            kind = member_kind(info)
            if kind in {"symlink", "special"}:
                raise CoordinateAuditError(f"unsupported {kind}: {info.filename}")
            if kind != "file":
                continue
            suffix = PurePosixPath(info.filename).suffix.lower()
            payload = archive.read(info)
            if suffix in {".parm7", ".prmtop"}:
                topology = topology_payload(info.filename, payload)
                topologies[topology["system_name"]] = topology
            elif suffix == ".rst":
                parsed_restraints.append(restraints.parse_restraint(info.filename, payload))
            elif suffix == ".binpos":
                trajectory_members.append((info.filename, payload))

    if set(topologies) != {"glucose", "galactose", "fructopyranose", "fructofuranose"}:
        raise CoordinateAuditError(f"unexpected topology systems: {sorted(topologies)}")
    if len(trajectory_members) != EXPECTED_TRAJECTORY_COUNT:
        raise CoordinateAuditError(
            f"expected {EXPECTED_TRAJECTORY_COUNT} BINPOS files, found {len(trajectory_members)}"
        )

    trajectory_rows: list[dict[str, Any]] = []
    mdtraj_versions: set[str] = set()
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        for member_path, payload in sorted(trajectory_members):
            system_name = infer_system(member_path)
            topology = topologies[system_name]
            restraint = choose_restraint_for_system(system_name, parsed_restraints)
            local = root / PurePosixPath(member_path).name
            local.write_bytes(payload)
            xyz_A, mdtraj_version = read_binpos(local)
            mdtraj_versions.add(mdtraj_version)
            if xyz_A.shape[1] != topology["natom"]:
                raise CoordinateAuditError(
                    f"{member_path}: trajectory has {xyz_A.shape[1]} atoms, "
                    f"topology has {topology['natom']}"
                )
            if not np.isfinite(xyz_A).all():
                raise CoordinateAuditError(f"{member_path}: non-finite coordinates")

            deposited_coordinate = reaction_coordinate_A(
                xyz_A,
                restraint["igr1"],
                restraint["igr2"],
                topology["masses"],
            )
            row: dict[str, Any] = {
                "path": member_path,
                "sha256": sha256_bytes(payload),
                "size_bytes": len(payload),
                "system_name": system_name,
                "topology_path": topology["path"],
                "topology_sha256": topology["sha256"],
                "coordinate_shape": list(xyz_A.shape),
                "coordinate_dtype_after_read": str(xyz_A.dtype),
                "first_frame_bond_statistics": first_frame_bond_statistics(
                    xyz_A, topology["bonds"]
                ),
                "deposited_reaction_coordinate": trajectory_summary(
                    deposited_coordinate
                ),
                "fructofuranose_alternatives": {},
            }

            if system_name == "fructofuranose":
                alternatives = {}
                for label, indices in sorted(furanose_alternatives.items()):
                    coordinate = reaction_coordinate_A(
                        xyz_A,
                        restraint["igr1"],
                        indices,
                        topology["masses"],
                    )
                    alternatives[label] = {
                        "atom_indices_1_based": indices,
                        "coordinate": trajectory_summary(coordinate),
                        "difference_from_deposited_H62_coordinate": delta_summary(
                            coordinate, deposited_coordinate
                        ),
                    }
                row["fructofuranose_alternatives"] = alternatives
            trajectory_rows.append(row)

    if len(mdtraj_versions) != 1:
        raise CoordinateAuditError(f"multiple MDTraj versions observed: {mdtraj_versions}")

    by_system: dict[str, list[dict[str, Any]]] = collections.defaultdict(list)
    for row in trajectory_rows:
        by_system[row["system_name"]].append(row)

    system_summaries: dict[str, Any] = {}
    for system_name, rows in sorted(by_system.items()):
        first_values = np.asarray(
            [row["deposited_reaction_coordinate"]["first_A"] for row in rows]
        )
        last_values = np.asarray(
            [row["deposited_reaction_coordinate"]["last_A"] for row in rows]
        )
        net_values = np.asarray(
            [row["deposited_reaction_coordinate"]["net_change_A"] for row in rows]
        )
        frame_counts = [
            row["deposited_reaction_coordinate"]["frame_count"] for row in rows
        ]
        system_summaries[system_name] = {
            "trajectory_count": len(rows),
            "frame_counts": frame_counts,
            "first_coordinate_A": percentile_summary(first_values),
            "last_coordinate_A": percentile_summary(last_values),
            "net_change_A": percentile_summary(net_values),
            "all_net_changes_positive": bool(np.all(net_values > 0.0)),
            "all_first_frame_median_bonds_in_angstrom_range": all(
                0.8
                <= row["first_frame_bond_statistics"]["median_A"]
                <= 1.8
                for row in rows
            ),
        }

    furanose_rows = by_system["fructofuranose"]
    alternative_cross_run: dict[str, Any] = {}
    for label in sorted(furanose_alternatives):
        mean_abs = np.asarray(
            [
                row["fructofuranose_alternatives"][label][
                    "difference_from_deposited_H62_coordinate"
                ]["mean_absolute_delta_A"]
                for row in furanose_rows
            ]
        )
        max_abs = np.asarray(
            [
                row["fructofuranose_alternatives"][label][
                    "difference_from_deposited_H62_coordinate"
                ]["maximum_absolute_delta_A"]
                for row in furanose_rows
            ]
        )
        last_delta = np.asarray(
            [
                row["fructofuranose_alternatives"][label][
                    "difference_from_deposited_H62_coordinate"
                ]["last_delta_A"]
                for row in furanose_rows
            ]
        )
        alternative_cross_run[label] = {
            "run_count": len(furanose_rows),
            "mean_absolute_delta_A_across_runs": percentile_summary(mean_abs),
            "maximum_absolute_delta_A_across_runs": percentile_summary(max_abs),
            "last_frame_delta_A_across_runs": percentile_summary(last_delta),
        }

    all_bond_ranges_ok = all(
        row["first_frame_bond_statistics"]["median_A"] >= 0.8
        and row["first_frame_bond_statistics"]["median_A"] <= 1.8
        for row in trajectory_rows
    )
    all_shapes_ok = all(
        row["coordinate_shape"][1:] == [EXPECTED_NATOM, 3]
        for row in trajectory_rows
    )

    report = {
        "schema": "dia-health-one.gluhut-binpos-coordinate-audit.v1",
        "status": "deposited_smd_coordinates_reconstructed_no_corrected_restraint_selected",
        "scope": "all public GluHUT SMD BINPOS trajectories and deposited COM groups",
        "not_target_claim": (
            "not NNC2215, PaG, MC-Q1, work, PMF, binding affinity, efficacy, safety, or a drug result"
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
            "numpy": np.__version__,
            "mdtraj": next(iter(mdtraj_versions)),
            "binpos_low_level_return_unit": "angstrom",
        },
        "candidate_record": {
            "path": str(candidates_path),
            "sha256": digest_file(candidates_path, "sha256"),
            "alternative_count": len(furanose_alternatives),
        },
        "trajectory_count": len(trajectory_rows),
        "total_frame_count": sum(
            row["deposited_reaction_coordinate"]["frame_count"]
            for row in trajectory_rows
        ),
        "system_summaries": system_summaries,
        "fructofuranose_alternative_cross_run_summary": alternative_cross_run,
        "trajectories": trajectory_rows,
        "decision": {
            "all_twenty_binpos_files_read": len(trajectory_rows)
            == EXPECTED_TRAJECTORY_COUNT,
            "all_trajectory_atom_shapes_match_topology": all_shapes_ok,
            "coordinate_unit_self_check_passed": all_bond_ranges_ok,
            "deposited_group_coordinate_reconstructed": True,
            "fructofuranose_alternative_coordinate_effect_quantified": True,
            "automatic_corrected_restraint_selected": False,
            "deposited_smd_work_reproduced": False,
            "new_smd_or_umbrella_generation_authorized": False,
            "parent_pag_calibration_authorized": False,
            "candidate_ranking_authorized": False,
            "next_action": (
                "use the coordinate-delta results to prioritize the author clarification "
                "and any preregistered multi-definition trajectory sensitivity; do not "
                "reinterpret deposited work or PMFs as corrected"
            ),
        },
        "limitations": [
            "BINPOS stores coordinates only; frame time and thermodynamic work are not reconstructed here.",
            "The recalculated COM distance uses compiled AMBER masses and deposited group membership.",
            "Alternative fructofuranose coordinates are post hoc measurements on trajectories generated with the deposited H62 group; they are not corrected trajectories.",
            "No periodic box is stored in BINPOS and no imaging correction is applied; the deposited visualisation trajectories are assumed to keep the complex contiguous.",
            "A small coordinate delta does not prove that an alternative restraint would generate the same path or work distribution.",
        ],
    }

    if not all(
        report["decision"][key]
        for key in (
            "all_twenty_binpos_files_read",
            "all_trajectory_atom_shapes_match_topology",
            "coordinate_unit_self_check_passed",
            "deposited_group_coordinate_reconstructed",
            "fructofuranose_alternative_coordinate_effect_quantified",
        )
    ):
        raise CoordinateAuditError("one or more coordinate-audit integrity checks failed")

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
        "--candidates-json",
        type=Path,
        default=Path(
            "research/data/gluhut-fructofuranose-restraint-candidates.json"
        ),
    )
    parser.add_argument(
        "--output-json",
        type=Path,
        default=Path("research/data/gluhut-binpos-coordinate-audit.json"),
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    report = audit(args.archive, args.candidates_json, args.output_json)
    print(
        json.dumps(
            {
                "software": report["software"],
                "trajectory_count": report["trajectory_count"],
                "total_frame_count": report["total_frame_count"],
                "system_summaries": report["system_summaries"],
                "fructofuranose_alternatives": report[
                    "fructofuranose_alternative_cross_run_summary"
                ],
                "decision": report["decision"],
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
