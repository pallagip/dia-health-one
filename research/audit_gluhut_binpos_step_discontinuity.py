#!/usr/bin/env python3
"""Test deposited GluHUT SMD coordinates for wrapping or jump artifacts.

The script reconstructs the deposited mass-weighted COM coordinate for all 20
BINPOS trajectories and measures every frame-to-frame increment.  A very large
single-frame jump would support a periodic-imaging or coordinate-discontinuity
explanation for an anomalous endpoint; a smooth progression to the anomalous
endpoint would instead indicate a different deposited pulling path or file.

This audit does not unwrap coordinates, alter trajectories, calculate work or
PMFs, or authorize PaG/MC-Q1 calculations.
"""
from __future__ import annotations

import argparse
import collections
import hashlib
import json
import stat
import tempfile
import zipfile
from pathlib import Path, PurePosixPath
from typing import Any

import numpy as np

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import audit_gluhut_binpos_coordinates as core  # noqa: E402
import audit_gluhut_smd_restraints as restraints  # noqa: E402
from binpos_native_reader import READER_VERSION, read_binpos  # noqa: E402

ARCHIVE_MD5 = "9f3c6733dc608d2bd37c488c560a5443"
ARCHIVE_SIZE_BYTES = 37931874
CHUNK_SIZE = 1024 * 1024
EXPECTED_TRAJECTORY_COUNT = 20
STEP_THRESHOLDS_A = (0.25, 0.5, 1.0, 2.0, 5.0)
PBC_LIKE_JUMP_THRESHOLD_A = 5.0


class StepAuditError(RuntimeError):
    pass


def digest_file(path: Path, algorithm: str) -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(CHUNK_SIZE), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def safe_member(name: str) -> PurePosixPath:
    if "\\" in name:
        raise StepAuditError(f"backslash in member path: {name!r}")
    path = PurePosixPath(name)
    if path.is_absolute() or any(part in {"", ".", ".."} for part in path.parts):
        raise StepAuditError(f"unsafe member path: {name!r}")
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


def percentile_summary(values: np.ndarray) -> dict[str, float]:
    if values.ndim != 1 or values.size == 0:
        raise StepAuditError("expected non-empty one-dimensional values")
    q = np.quantile(values, [0.0, 0.5, 0.9, 0.95, 0.99, 1.0])
    return {
        "minimum": float(q[0]),
        "median": float(q[1]),
        "p90": float(q[2]),
        "p95": float(q[3]),
        "p99": float(q[4]),
        "maximum": float(q[5]),
        "mean": float(values.mean()),
        "standard_deviation": float(values.std(ddof=0)),
    }


def step_summary(coordinate: np.ndarray) -> dict[str, Any]:
    if coordinate.ndim != 1 or coordinate.size < 2:
        raise StepAuditError("coordinate must contain at least two frames")
    signed = np.diff(coordinate)
    absolute = np.abs(signed)
    maximum_index = int(np.argmax(absolute))
    counts = {
        str(threshold): int(np.sum(absolute > threshold))
        for threshold in STEP_THRESHOLDS_A
    }
    first_quarter = max(1, coordinate.size // 4)
    last_quarter_start = coordinate.size - first_quarter
    return {
        "frame_step_count": int(signed.size),
        "signed_step_summary_A": percentile_summary(signed),
        "absolute_step_summary_A": percentile_summary(absolute),
        "counts_above_threshold_A": counts,
        "maximum_absolute_step_A": float(absolute[maximum_index]),
        "maximum_step_signed_A": float(signed[maximum_index]),
        "maximum_step_from_frame_index_0_based": maximum_index,
        "maximum_step_to_frame_index_0_based": maximum_index + 1,
        "coordinate_before_maximum_step_A": float(coordinate[maximum_index]),
        "coordinate_after_maximum_step_A": float(coordinate[maximum_index + 1]),
        "pbc_like_jump_detected": bool(
            absolute[maximum_index] > PBC_LIKE_JUMP_THRESHOLD_A
        ),
        "first_quarter_coordinate_mean_A": float(
            coordinate[:first_quarter].mean()
        ),
        "last_quarter_coordinate_mean_A": float(
            coordinate[last_quarter_start:].mean()
        ),
        "fraction_positive_steps": float(np.mean(signed > 0.0)),
        "fraction_negative_steps": float(np.mean(signed < 0.0)),
    }


def audit(archive_path: Path, output_path: Path) -> dict[str, Any]:
    if archive_path.stat().st_size != ARCHIVE_SIZE_BYTES:
        raise StepAuditError(
            f"size mismatch: expected {ARCHIVE_SIZE_BYTES}, got {archive_path.stat().st_size}"
        )
    md5 = digest_file(archive_path, "md5")
    if md5 != ARCHIVE_MD5:
        raise StepAuditError(f"MD5 mismatch: expected {ARCHIVE_MD5}, got {md5}")
    archive_sha256 = digest_file(archive_path, "sha256")

    topologies: dict[str, dict[str, Any]] = {}
    parsed_restraints: list[dict[str, Any]] = []
    trajectories: list[tuple[str, bytes]] = []
    seen: set[str] = set()
    with zipfile.ZipFile(archive_path) as archive:
        bad = archive.testzip()
        if bad:
            raise StepAuditError(f"ZIP CRC failure: {bad}")
        for info in archive.infolist():
            safe_member(info.filename)
            if info.filename in seen:
                raise StepAuditError(f"duplicate member: {info.filename}")
            seen.add(info.filename)
            if info.flag_bits & 0x1:
                raise StepAuditError(f"encrypted member: {info.filename}")
            kind = member_kind(info)
            if kind in {"symlink", "special"}:
                raise StepAuditError(f"unsupported {kind}: {info.filename}")
            if kind != "file":
                continue
            suffix = PurePosixPath(info.filename).suffix.lower()
            payload = archive.read(info)
            if suffix in {".parm7", ".prmtop"}:
                topology = core.topology_payload(info.filename, payload)
                topologies[topology["system_name"]] = topology
            elif suffix == ".rst":
                parsed_restraints.append(
                    restraints.parse_restraint(info.filename, payload)
                )
            elif suffix == ".binpos":
                trajectories.append((info.filename, payload))

    if len(trajectories) != EXPECTED_TRAJECTORY_COUNT:
        raise StepAuditError(
            f"expected {EXPECTED_TRAJECTORY_COUNT} trajectories, found {len(trajectories)}"
        )

    rows: list[dict[str, Any]] = []
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        for member_path, payload in sorted(trajectories):
            system_name = core.infer_system(member_path)
            topology = topologies[system_name]
            restraint = core.choose_restraint_for_system(
                system_name, parsed_restraints
            )
            local = root / PurePosixPath(member_path).name
            local.write_bytes(payload)
            xyz_A, reader_version = read_binpos(
                local, expected_natoms=topology["natom"]
            )
            if reader_version != READER_VERSION:
                raise StepAuditError(
                    f"unexpected reader version {reader_version}"
                )
            coordinate = core.reaction_coordinate_A(
                xyz_A,
                restraint["igr1"],
                restraint["igr2"],
                topology["masses"],
            )
            summary = step_summary(coordinate)
            rows.append(
                {
                    "path": member_path,
                    "file_sha256": sha256_bytes(payload),
                    "system_name": system_name,
                    "frame_count": int(coordinate.size),
                    "first_A": float(coordinate[0]),
                    "last_A": float(coordinate[-1]),
                    "net_change_A": float(coordinate[-1] - coordinate[0]),
                    "step_analysis": summary,
                }
            )

    by_system: dict[str, list[dict[str, Any]]] = collections.defaultdict(list)
    for row in rows:
        by_system[row["system_name"]].append(row)

    system_summaries: dict[str, Any] = {}
    for system_name, system_rows in sorted(by_system.items()):
        max_steps = np.asarray(
            [
                row["step_analysis"]["maximum_absolute_step_A"]
                for row in system_rows
            ],
            dtype=float,
        )
        system_summaries[system_name] = {
            "trajectory_count": len(system_rows),
            "pbc_like_jump_trajectory_count": sum(
                bool(row["step_analysis"]["pbc_like_jump_detected"])
                for row in system_rows
            ),
            "maximum_absolute_step_A_across_trajectories": (
                percentile_summary(max_steps)
            ),
            "paths_with_pbc_like_jump": [
                row["path"]
                for row in system_rows
                if row["step_analysis"]["pbc_like_jump_detected"]
            ],
        }

    galactose_A = [
        row
        for row in rows
        if row["system_name"] == "galactose"
        and PurePosixPath(row["path"]).name.endswith("-A.binpos")
    ]
    galactose_B = [
        row
        for row in rows
        if row["system_name"] == "galactose"
        and PurePosixPath(row["path"]).name.endswith("-B.binpos")
    ]
    if len(galactose_A) != 2 or len(galactose_B) != 3:
        raise StepAuditError(
            f"expected 2 galactose A and 3 B files; got {len(galactose_A)} and {len(galactose_B)}"
        )

    galactose_A_has_large_jump = any(
        row["step_analysis"]["pbc_like_jump_detected"]
        for row in galactose_A
    )
    galactose_A_smooth_to_anomalous_endpoint = all(
        not row["step_analysis"]["pbc_like_jump_detected"]
        and row["last_A"] > 20.0
        for row in galactose_A
    )

    report = {
        "schema": "dia-health-one.gluhut-binpos-step-discontinuity.v1",
        "status": (
            "galactose_A_large_jump_detected"
            if galactose_A_has_large_jump
            else "galactose_A_no_large_single_frame_jump"
        ),
        "scope": "frame-to-frame continuity of all public GluHUT SMD COM coordinates",
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
        "settings": {
            "step_thresholds_A": list(STEP_THRESHOLDS_A),
            "pbc_like_jump_threshold_A": PBC_LIKE_JUMP_THRESHOLD_A,
            "binpos_reader": READER_VERSION,
        },
        "trajectory_count": len(rows),
        "total_frame_count": sum(row["frame_count"] for row in rows),
        "system_summaries": system_summaries,
        "galactose_path_A": galactose_A,
        "galactose_path_B": galactose_B,
        "trajectories": rows,
        "decision": {
            "all_twenty_trajectories_analysed": len(rows)
            == EXPECTED_TRAJECTORY_COUNT,
            "galactose_A_contains_pbc_like_single_frame_jump": (
                galactose_A_has_large_jump
            ),
            "galactose_A_progresses_smoothly_to_endpoint_above_20A": (
                galactose_A_smooth_to_anomalous_endpoint
            ),
            "periodic_imaging_explanation_supported": (
                galactose_A_has_large_jump
            ),
            "different_or_mislabeled_deposited_path_supported": (
                galactose_A_smooth_to_anomalous_endpoint
            ),
            "automatic_unwrapping_or_correction_authorized": False,
            "new_smd_or_umbrella_generation_authorized": False,
            "parent_pag_calibration_authorized": False,
            "candidate_ranking_authorized": False,
            "next_action": (
                "request author clarification and the missing intended galactose A replicate; "
                "do not treat the duplicated A file or mixed endpoints as five independent runs"
            ),
        },
        "limitations": [
            "The 5 A threshold is a conservative discontinuity screen, not a formal periodic-boundary reconstruction.",
            "BINPOS contains no box vectors, so minimum-image unwrapping cannot be performed from this archive alone.",
            "Smooth coordinate progression does not distinguish an intentional longer target from a mislabeled or wrong trajectory.",
            "The audit does not reconstruct applied target versus time or nonequilibrium work.",
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
        default=Path(
            "research/data/gluhut-binpos-step-discontinuity.json"
        ),
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    report = audit(args.archive, args.output_json)
    print(
        json.dumps(
            {
                "system_summaries": report["system_summaries"],
                "galactose_A": report["galactose_path_A"],
                "galactose_B": report["galactose_path_B"],
                "decision": report["decision"],
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
