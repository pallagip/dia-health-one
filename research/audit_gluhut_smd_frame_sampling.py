#!/usr/bin/env python3
"""Audit frame-count and target-progress provenance in the public GluHUT SMD deposit.

The deposited AMBER control requests a 30 ns trajectory and writes coordinates
every 1000 integration steps.  The public BINPOS files contain 1000 frames each
and are often named ``visualisation``.  This audit compares the requested output
cadence, deposited frame count, and reconstructed COM progress without assuming
that omitted frames were selected uniformly.

It does not reconstruct frame times, target histories, nonequilibrium work or
PMFs, and it does not authorize PaG or MC-Q1 calculations.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from pathlib import Path
from typing import Any

EXPECTED_SMD_TARGET_CHANGE_A = 15.0
NORMAL_NET_CHANGE_TOLERANCE_A = 0.5


class FrameSamplingError(RuntimeError):
    pass


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_assignments(lines: list[str]) -> dict[str, float]:
    text = " ".join(lines)
    result: dict[str, float] = {}
    for key in ("nstlim", "dt", "ntwx", "ntwr", "ntpr"):
        match = re.search(
            rf"(?i)(?:^|[\s,]){re.escape(key)}\s*=\s*([0-9.+\-EeDd]+)",
            text,
        )
        if not match:
            raise FrameSamplingError(f"SMD control preview is missing {key}")
        result[key] = float(match.group(1).replace("D", "E").replace("d", "e"))
    return result


def find_smd_preview(inventory: dict[str, Any]) -> tuple[str, list[str]]:
    matches = [
        row
        for row in inventory.get("members", [])
        if row.get("path", "").endswith("/SMD.in")
    ]
    if len(matches) != 1:
        raise FrameSamplingError(f"expected one SMD.in, found {len(matches)}")
    preview = matches[0].get("preview")
    if not preview or preview.get("encoding") != "utf-8":
        raise FrameSamplingError("SMD.in has no UTF-8 preview")
    return str(matches[0]["path"]), list(preview.get("head", []))


def classify_trajectory(row: dict[str, Any]) -> dict[str, Any]:
    coordinate = row["deposited_reaction_coordinate"]
    frame_count = int(coordinate["frame_count"])
    if frame_count < 2:
        raise FrameSamplingError(f"trajectory has fewer than two frames: {row['path']}")
    net_change = float(coordinate["net_change_A"])
    mean_signed_step = net_change / (frame_count - 1)
    deviation = net_change - EXPECTED_SMD_TARGET_CHANGE_A
    return {
        "system_name": row["system_name"],
        "path": row["path"],
        "file_sha256": row["sha256"],
        "frame_count": frame_count,
        "first_A": float(coordinate["first_A"]),
        "last_A": float(coordinate["last_A"]),
        "net_change_A": net_change,
        "mean_signed_coordinate_increment_per_deposited_frame_A": mean_signed_step,
        "deviation_from_reported_15A_target_change_A": deviation,
        "consistent_with_reported_15A_change_within_tolerance": (
            abs(deviation) <= NORMAL_NET_CHANGE_TOLERANCE_A
        ),
        "filename_contains_visualisation": "visualisation" in row["path"].lower(),
    }


def summarize(
    inventory_path: Path,
    coordinate_path: Path,
    output_path: Path,
) -> dict[str, Any]:
    inventory = json.loads(inventory_path.read_text(encoding="utf-8"))
    coordinates = json.loads(coordinate_path.read_text(encoding="utf-8"))

    control_path, control_lines = find_smd_preview(inventory)
    control = parse_assignments(control_lines)
    nstlim = int(control["nstlim"])
    ntwx = int(control["ntwx"])
    dt_ps = float(control["dt"])
    if nstlim <= 0 or ntwx <= 0 or nstlim % ntwx:
        raise FrameSamplingError(
            f"invalid nstlim/ntwx combination: {nstlim}/{ntwx}"
        )
    expected_written_frames = nstlim // ntwx
    simulation_time_ps = nstlim * dt_ps
    simulation_time_ns = simulation_time_ps / 1000.0

    rows = [classify_trajectory(row) for row in coordinates["trajectories"]]
    if len(rows) != 20:
        raise FrameSamplingError(f"expected 20 trajectories, found {len(rows)}")
    deposited_frame_counts = sorted({row["frame_count"] for row in rows})
    if len(deposited_frame_counts) != 1:
        raise FrameSamplingError(
            f"deposited trajectories have multiple frame counts: {deposited_frame_counts}"
        )
    deposited_frames = deposited_frame_counts[0]
    output_to_deposit_ratio = expected_written_frames / deposited_frames

    by_system: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        by_system.setdefault(row["system_name"], []).append(row)

    system_summaries: dict[str, Any] = {}
    for system_name, system_rows in sorted(by_system.items()):
        normal = [
            row
            for row in system_rows
            if row["consistent_with_reported_15A_change_within_tolerance"]
        ]
        anomalous = [
            row
            for row in system_rows
            if not row["consistent_with_reported_15A_change_within_tolerance"]
        ]
        system_summaries[system_name] = {
            "trajectory_count": len(system_rows),
            "reported_15A_consistent_count": len(normal),
            "reported_15A_inconsistent_count": len(anomalous),
            "consistent_paths": [row["path"] for row in normal],
            "inconsistent_paths": [row["path"] for row in anomalous],
            "mean_net_change_A": sum(row["net_change_A"] for row in system_rows)
            / len(system_rows),
            "mean_signed_increment_per_deposited_frame_A": sum(
                row["mean_signed_coordinate_increment_per_deposited_frame_A"]
                for row in system_rows
            )
            / len(system_rows),
        }

    expected_raw_output_target_increment_A = (
        EXPECTED_SMD_TARGET_CHANGE_A / expected_written_frames
    )
    expected_deposited_increment_if_uniform_15x_subsample_A = (
        EXPECTED_SMD_TARGET_CHANGE_A / (deposited_frames - 1)
    )
    normal_rows = [
        row
        for row in rows
        if row["consistent_with_reported_15A_change_within_tolerance"]
    ]
    anomalous_rows = [
        row
        for row in rows
        if not row["consistent_with_reported_15A_change_within_tolerance"]
    ]

    report = {
        "schema": "dia-health-one.gluhut-smd-frame-sampling.v1",
        "status": "visualisation_files_are_downsampled_or_incomplete_time_mapping_unavailable",
        "scope": "requested AMBER output cadence versus deposited BINPOS frame counts and COM progress",
        "not_target_claim": (
            "not NNC2215, PaG, MC-Q1, SMD work, PMF, binding affinity, efficacy, safety, or a drug result"
        ),
        "inputs": {
            "archive_inventory_path": str(inventory_path),
            "archive_inventory_sha256": sha256(inventory_path),
            "coordinate_audit_path": str(coordinate_path),
            "coordinate_audit_sha256": sha256(coordinate_path),
            "smd_control_path": control_path,
        },
        "smd_control": {
            "nstlim": nstlim,
            "dt_ps": dt_ps,
            "simulation_time_ps": simulation_time_ps,
            "simulation_time_ns": simulation_time_ns,
            "ntwx_steps": ntwx,
            "ntwr_steps": int(control["ntwr"]),
            "ntpr_steps": int(control["ntpr"]),
            "expected_coordinate_outputs_from_nstlim_div_ntwx": expected_written_frames,
            "reported_target_change_A": EXPECTED_SMD_TARGET_CHANGE_A,
            "nominal_target_increment_per_requested_output_A": (
                expected_raw_output_target_increment_A
            ),
        },
        "deposit": {
            "trajectory_count": len(rows),
            "deposited_frames_per_trajectory": deposited_frames,
            "total_deposited_frames": sum(row["frame_count"] for row in rows),
            "requested_output_to_deposited_frame_ratio": output_to_deposit_ratio,
            "filename_contains_visualisation_count": sum(
                row["filename_contains_visualisation"] for row in rows
            ),
            "nominal_increment_if_uniformly_reduced_to_deposited_frame_count_A": (
                expected_deposited_increment_if_uniform_15x_subsample_A
            ),
            "uniform_downsampling_provenance_available": False,
            "frame_time_mapping_available": False,
        },
        "trajectory_classification": {
            "reported_15A_consistent_count": len(normal_rows),
            "reported_15A_inconsistent_count": len(anomalous_rows),
            "reported_15A_consistent_paths": [row["path"] for row in normal_rows],
            "reported_15A_inconsistent_paths": [row["path"] for row in anomalous_rows],
            "normal_change_tolerance_A": NORMAL_NET_CHANGE_TOLERANCE_A,
        },
        "system_summaries": system_summaries,
        "trajectories": rows,
        "decision": {
            "control_requests_30ns": math.isclose(simulation_time_ns, 30.0),
            "control_would_write_15000_coordinate_frames": (
                expected_written_frames == 15000
            ),
            "every_deposited_file_contains_1000_frames": deposited_frames == 1000,
            "requested_to_deposited_frame_ratio_is_15": math.isclose(
                output_to_deposit_ratio, 15.0
            ),
            "deposit_proves_uniform_every_15th_frame_selection": False,
            "thirteen_files_match_reported_15A_change_within_0_5A": (
                len(normal_rows) == 13
            ),
            "seven_files_do_not_match_reported_15A_change": (
                len(anomalous_rows) == 7
            ),
            "deposited_files_are_sufficient_for_time_resolved_work_reconstruction": False,
            "new_smd_or_umbrella_generation_authorized": False,
            "parent_pag_calibration_authorized": False,
            "candidate_ranking_authorized": False,
            "next_action": (
                "request the original full trajectories or exact frame-selection script and "
                "time mapping; request target histories for fructopyranose and galactose A"
            ),
        },
        "limitations": [
            "The deposited files contain coordinates but no per-frame time or target value.",
            "A 15-fold frame-count ratio is consistent with uniform downsampling but does not prove which source frames were retained.",
            "The first deposited frame may not correspond exactly to integration step zero.",
            "COM coordinate progress is not the applied target history or nonequilibrium work.",
            "Files named visualisation may have been transformed or subsampled for presentation, but no transformation script is deposited.",
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
    parser.add_argument(
        "--inventory-json",
        type=Path,
        default=Path("research/data/gluhut-smd-archive-inventory.json"),
    )
    parser.add_argument(
        "--coordinate-json",
        type=Path,
        default=Path("research/data/gluhut-binpos-coordinate-audit.json"),
    )
    parser.add_argument(
        "--output-json",
        type=Path,
        default=Path("research/data/gluhut-smd-frame-sampling.json"),
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    report = summarize(
        args.inventory_json,
        args.coordinate_json,
        args.output_json,
    )
    print(
        json.dumps(
            {
                "smd_control": report["smd_control"],
                "deposit": report["deposit"],
                "trajectory_classification": report[
                    "trajectory_classification"
                ],
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
