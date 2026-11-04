#!/usr/bin/env python3
"""Summarize endpoint provenance of the deposited GluHUT SMD trajectories.

The input is the deterministic coordinate audit for all twenty BINPOS files.
This script preserves file-level path, content hash and coordinate hash,
clusters final centre-of-mass coordinates, detects duplicate trajectories, and
flags systems whose files do not form one endpoint protocol.

It does not infer why an endpoint differs, modify trajectories, calculate work
or PMFs, or authorize prospective PaG/MC-Q1 calculations.
"""
from __future__ import annotations

import argparse
import collections
import hashlib
import json
from pathlib import Path
from typing import Any

ENDPOINT_CLUSTER_TOLERANCE_A = 0.25
EXPECTED_TARGETS_A = (15.0, 20.0)
TARGET_MATCH_TOLERANCE_A = 0.25


class EndpointAuditError(RuntimeError):
    pass


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def cluster_rows(
    rows: list[dict[str, Any]],
    tolerance_A: float,
) -> list[dict[str, Any]]:
    """Greedily cluster sorted endpoints within a bounded diameter."""
    sorted_rows = sorted(
        rows,
        key=lambda row: (
            float(row["last_A"]),
            str(row["path"]),
        ),
    )
    clusters: list[list[dict[str, Any]]] = []
    for row in sorted_rows:
        if not clusters:
            clusters.append([row])
            continue
        current = clusters[-1]
        candidate_values = [float(item["last_A"]) for item in current] + [
            float(row["last_A"])
        ]
        if max(candidate_values) - min(candidate_values) <= tolerance_A:
            current.append(row)
        else:
            clusters.append([row])

    result = []
    for index, cluster in enumerate(clusters, start=1):
        values = [float(row["last_A"]) for row in cluster]
        result.append(
            {
                "cluster_index": index,
                "trajectory_count": len(cluster),
                "minimum_last_A": min(values),
                "maximum_last_A": max(values),
                "mean_last_A": sum(values) / len(values),
                "paths": [row["path"] for row in cluster],
                "file_sha256": [row["file_sha256"] for row in cluster],
                "coordinate_sha256": [row["coordinate_sha256"] for row in cluster],
            }
        )
    return result


def closest_expected_target(value: float) -> dict[str, Any]:
    distances = [
        (abs(value - target), target) for target in EXPECTED_TARGETS_A
    ]
    distance, target = min(distances)
    return {
        "closest_target_A": target,
        "absolute_difference_A": distance,
        "matches_within_tolerance": distance <= TARGET_MATCH_TOLERANCE_A,
    }


def trajectory_rows(data: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    for trajectory in data.get("trajectories", []):
        deposited = trajectory["deposited_reaction_coordinate"]
        rows.append(
            {
                "system_name": trajectory["system_name"],
                "path": trajectory["path"],
                "file_sha256": trajectory["sha256"],
                "coordinate_sha256": deposited[
                    "coordinate_sha256_float64"
                ],
                "frame_count": deposited["frame_count"],
                "first_A": deposited["first_A"],
                "last_A": deposited["last_A"],
                "net_change_A": deposited["net_change_A"],
                "last_minus_target_15A": deposited[
                    "last_minus_target_15A"
                ],
                "fraction_positive_frame_steps": deposited[
                    "fraction_positive_frame_steps"
                ],
                "closest_expected_target": closest_expected_target(
                    float(deposited["last_A"])
                ),
            }
        )
    if len(rows) != 20:
        raise EndpointAuditError(
            f"expected 20 trajectory rows, found {len(rows)}"
        )
    return sorted(rows, key=lambda row: (row["system_name"], row["path"]))


def duplicate_groups(
    rows: list[dict[str, Any]],
    key: str,
) -> list[dict[str, Any]]:
    grouped: dict[str, list[dict[str, Any]]] = collections.defaultdict(list)
    for row in rows:
        grouped[str(row[key])].append(row)
    return [
        {
            "hash": digest,
            "trajectory_count": len(group),
            "paths": [row["path"] for row in group],
            "systems": sorted({row["system_name"] for row in group}),
        }
        for digest, group in sorted(grouped.items())
        if len(group) > 1
    ]


def summarize(input_path: Path, output_path: Path) -> dict[str, Any]:
    data = json.loads(input_path.read_text(encoding="utf-8"))
    rows = trajectory_rows(data)

    by_system: dict[str, list[dict[str, Any]]] = collections.defaultdict(list)
    for row in rows:
        by_system[row["system_name"]].append(row)

    systems = {}
    heterogeneous_systems = []
    target_mismatch_rows = []
    for system_name, system_rows in sorted(by_system.items()):
        clusters = cluster_rows(
            system_rows, ENDPOINT_CLUSTER_TOLERANCE_A
        )
        heterogeneous = len(clusters) > 1
        if heterogeneous:
            heterogeneous_systems.append(system_name)
        for row in system_rows:
            if not row["closest_expected_target"][
                "matches_within_tolerance"
            ]:
                target_mismatch_rows.append(row)
        systems[system_name] = {
            "trajectory_count": len(system_rows),
            "endpoint_cluster_count": len(clusters),
            "endpoint_protocol_heterogeneous": heterogeneous,
            "endpoint_clusters": clusters,
            "trajectories": system_rows,
        }

    duplicate_files = duplicate_groups(rows, "file_sha256")
    duplicate_coordinates = duplicate_groups(rows, "coordinate_sha256")

    report = {
        "schema": "dia-health-one.gluhut-smd-endpoint-provenance.v1",
        "status": (
            "endpoint_protocol_heterogeneity_detected"
            if heterogeneous_systems
            else "one_endpoint_cluster_per_system"
        ),
        "scope": "file-level endpoint provenance of twenty public GluHUT SMD BINPOS trajectories",
        "not_target_claim": (
            "not NNC2215, PaG, MC-Q1, SMD work, PMF, binding affinity, efficacy, safety, or a drug result"
        ),
        "source": {
            "path": str(input_path),
            "sha256": sha256(input_path),
            "trajectory_count": len(rows),
            "total_frame_count": sum(row["frame_count"] for row in rows),
        },
        "settings": {
            "endpoint_cluster_tolerance_A": ENDPOINT_CLUSTER_TOLERANCE_A,
            "expected_targets_A": list(EXPECTED_TARGETS_A),
            "target_match_tolerance_A": TARGET_MATCH_TOLERANCE_A,
        },
        "systems": systems,
        "duplicate_file_content_groups": duplicate_files,
        "duplicate_coordinate_series_groups": duplicate_coordinates,
        "trajectories_outside_15_or_20A_target_tolerance": (
            target_mismatch_rows
        ),
        "decision": {
            "all_twenty_trajectories_accounted_for": len(rows) == 20,
            "one_endpoint_protocol_per_system": not heterogeneous_systems,
            "heterogeneous_systems": heterogeneous_systems,
            "duplicate_file_content_detected": bool(duplicate_files),
            "duplicate_coordinate_series_detected": bool(
                duplicate_coordinates
            ),
            "new_smd_or_umbrella_generation_authorized": False,
            "parent_pag_calibration_authorized": False,
            "candidate_ranking_authorized": False,
            "next_action": (
                "inspect file naming and restraint provenance for each endpoint cluster; "
                "request author clarification for mixed-target or duplicate trajectories; "
                "do not pool heterogeneous files as one protocol"
                if heterogeneous_systems
                else "retain per-system endpoint definitions and continue source-coordinate provenance audit"
            ),
        },
        "limitations": [
            "Endpoint clusters are based only on the reconstructed deposited COM coordinate and a prespecified 0.25 A tolerance.",
            "A 15 or 20 A endpoint does not by itself establish the intended target or production duration.",
            "File-content or coordinate-series equality detects exact duplicates only, not trajectories sharing a starting structure or random seed.",
            "The audit does not reconstruct time, work, applied target history, or umbrella-window generation.",
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
        "--input-json",
        type=Path,
        default=Path(
            "research/data/gluhut-binpos-coordinate-audit.json"
        ),
    )
    parser.add_argument(
        "--output-json",
        type=Path,
        default=Path(
            "research/data/gluhut-smd-endpoint-provenance.json"
        ),
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    report = summarize(args.input_json, args.output_json)
    print(
        json.dumps(
            {
                "systems": {
                    name: {
                        "endpoint_cluster_count": value[
                            "endpoint_cluster_count"
                        ],
                        "endpoint_clusters": value[
                            "endpoint_clusters"
                        ],
                    }
                    for name, value in report["systems"].items()
                },
                "duplicate_files": report[
                    "duplicate_file_content_groups"
                ],
                "duplicate_coordinates": report[
                    "duplicate_coordinate_series_groups"
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
