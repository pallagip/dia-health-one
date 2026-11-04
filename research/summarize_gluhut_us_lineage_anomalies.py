#!/usr/bin/env python3
"""Summarize material provenance failures in the GluHUT umbrella archive.

The full umbrella inventory and route-lineage records are intentionally verbose.
This module extracts only the confidence-changing anomalies:

- files named ``100-frames`` that contain fewer recorded frames;
- WHAM windows that have no corresponding trajectory file;
- exact duplicate windows within a route;
- exact duplicate complete routes;
- the gap between deposited visualization trajectories and full histogram input.

It does not regenerate histograms or PMFs and does not authorize PaG/MC-Q1 work.
"""
from __future__ import annotations

import argparse
import collections
import hashlib
import json
from pathlib import Path
from typing import Any

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import summarize_gluhut_us_route_lineage as lineage_tools  # noqa: E402


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def round_center(value: float) -> float:
    return round(float(value), 6)


def mode(values: list[int]) -> int | None:
    if not values:
        return None
    counts = collections.Counter(values)
    highest = max(counts.values())
    winners = sorted(value for value, count in counts.items() if count == highest)
    return winners[0] if len(winners) == 1 else None


def group_rows(rows: list[dict[str, Any]]) -> dict[str, list[dict[str, Any]]]:
    grouped: dict[str, list[dict[str, Any]]] = collections.defaultdict(list)
    for row in rows:
        grouped[row["route_key"]].append(row)
    for route_rows in grouped.values():
        route_rows.sort(key=lambda row: row["center_tenths"])
    return dict(sorted(grouped.items()))


def frame_anomalies(
    rows_by_route: dict[str, list[dict[str, Any]]],
    expected_frames: int = 100,
) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    for route_key, rows in rows_by_route.items():
        anomalies = [
            {
                "center_A": row["center_A"],
                "frame_count": row["frame_count"],
                "path": row["path"],
                "sha256": row["sha256"],
            }
            for row in rows
            if row["frame_count"] != expected_frames
        ]
        if not anomalies:
            continue
        all_counts = [row["frame_count"] for row in rows]
        result.append(
            {
                "route_key": route_key,
                "mapped_wham_key": rows[0]["mapped_wham_key"],
                "expected_frames_from_filename": expected_frames,
                "unique_frame_counts": sorted(set(all_counts)),
                "modal_frame_count": mode(all_counts),
                "anomalous_file_count": len(anomalies),
                "total_file_count": len(rows),
                "anomalies": anomalies,
            }
        )
    return result


def center_mismatches(
    route_summaries: list[dict[str, Any]],
    wham_by_key: dict[str, dict[str, Any]],
) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    for route in route_summaries:
        trajectory_centers = [round_center(value) for value in route["window_centres_A"]]
        wham_route = wham_by_key[route["mapped_wham_key"]]
        wham_centers = [round_center(value) for value in wham_route["centres_A"]]
        trajectory_only = sorted(set(trajectory_centers) - set(wham_centers))
        wham_only = sorted(set(wham_centers) - set(trajectory_centers))
        if not trajectory_only and not wham_only:
            continue
        result.append(
            {
                "route_key": route["route_key"],
                "mapped_wham_key": route["mapped_wham_key"],
                "trajectory_window_count": len(trajectory_centers),
                "wham_window_count": len(wham_centers),
                "trajectory_only_centres_A": trajectory_only,
                "wham_only_centres_A": wham_only,
                "trajectory_centres_A": trajectory_centers,
                "wham_centres_A": wham_centers,
            }
        )
    return result


def histogram_visualization_ratios(
    route_summaries: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    result = []
    for route in route_summaries:
        values = route["trajectory_frames_per_window_values"]
        modal = mode([int(value) for value in values])
        result.append(
            {
                "route_key": route["route_key"],
                "mapped_wham_key": route["mapped_wham_key"],
                "trajectory_frame_count_values": values,
                "modal_trajectory_frame_count": modal,
                "wham_minimum_frames": route["wham_frames_minimum_window"],
                "wham_median_frames": route["wham_frames_median_window"],
                "wham_maximum_frames": route["wham_frames_maximum_window"],
                "median_histogram_to_modal_trajectory_ratio": (
                    route["wham_frames_median_window"] / modal
                    if modal
                    else None
                ),
            }
        )
    return result


def run(
    inventory_path: Path,
    wham_path: Path,
    lineage_path: Path,
    output_path: Path,
) -> dict[str, Any]:
    inventory = json.loads(inventory_path.read_text(encoding="utf-8"))
    wham = json.loads(wham_path.read_text(encoding="utf-8"))
    lineage = json.loads(lineage_path.read_text(encoding="utf-8"))

    trajectory_rows = lineage_tools.trajectory_rows(inventory)
    rows_by_route = group_rows(trajectory_rows)
    wham_by_key = lineage_tools.wham_routes(wham)

    frame_failures = frame_anomalies(rows_by_route)
    center_failures = center_mismatches(lineage["route_summaries"], wham_by_key)
    ratio_rows = histogram_visualization_ratios(lineage["route_summaries"])

    frame_failure_keys = {row["route_key"] for row in frame_failures}
    center_failure_keys = {row["route_key"] for row in center_failures}
    all_center_failures_are_pathb_225 = (
        center_failure_keys
        == {"galactose/run1b", "galactose/run2b", "galactose/run3b"}
        and all(
            row["trajectory_only_centres_A"] == []
            and row["wham_only_centres_A"] == [2.25]
            and row["trajectory_window_count"] == 51
            and row["wham_window_count"] == 52
            for row in center_failures
        )
    )

    furanose_run2 = next(
        (row for row in frame_failures if row["route_key"] == "fructofuranose/run2"),
        None,
    )
    furanose_run3 = next(
        (row for row in frame_failures if row["route_key"] == "fructofuranose/run3"),
        None,
    )
    run2_one_frame_centres = (
        [row["center_A"] for row in furanose_run2["anomalies"] if row["frame_count"] == 1]
        if furanose_run2
        else []
    )
    run2_ten_frame_centres = (
        [row["center_A"] for row in furanose_run2["anomalies"] if row["frame_count"] == 10]
        if furanose_run2
        else []
    )

    wham_anomalies = wham.get("anomalies", [])
    same_center_as_wham_short_window = any(
        float(row.get("center_A")) in run2_one_frame_centres
        and "fructofuranose-US-2" in str(row.get("path", ""))
        for row in wham_anomalies
    )

    report = {
        "schema": "dia-health-one.gluhut-us-lineage-anomalies.v1",
        "status": "umbrella_archive_lineage_failures_explicit_histogram_regeneration_blocked",
        "scope": "material route/window/frame anomalies in the public GluHUT umbrella archive",
        "not_target_claim": (
            "not regenerated histograms, not a new PMF or affinity, not NNC2215, "
            "PaG, MC-Q1, efficacy, safety, or a drug result"
        ),
        "inputs": {
            "umbrella_inventory": {
                "path": str(inventory_path),
                "sha256": file_sha256(inventory_path),
                "archive": inventory["archive"],
            },
            "wham_reproduction": {
                "path": str(wham_path),
                "sha256": file_sha256(wham_path),
                "archive": wham["archive"],
            },
            "route_lineage": {
                "path": str(lineage_path),
                "sha256": file_sha256(lineage_path),
            },
        },
        "summary": {
            "trajectory_file_count": lineage["trajectory_file_count"],
            "raw_route_directory_count": lineage["raw_route_directory_count"],
            "unique_complete_route_content_signature_count": lineage[
                "unique_complete_route_content_signature_count"
            ],
            "deposited_wham_route_count": lineage["deposited_wham_route_count"],
            "frame_failure_route_count": len(frame_failures),
            "center_mismatch_route_count": len(center_failures),
            "complete_duplicate_route_group_count": len(
                lineage["duplicate_complete_route_groups"]
            ),
            "within_route_duplicate_route_count": len(
                lineage["within_route_duplicate_groups_by_route"]
            ),
        },
        "frame_count_failures": frame_failures,
        "center_mismatches": center_failures,
        "histogram_to_visualization_ratios": ratio_rows,
        "complete_duplicate_route_groups": lineage[
            "duplicate_complete_route_groups"
        ],
        "within_route_duplicate_groups_by_route": lineage[
            "within_route_duplicate_groups_by_route"
        ],
        "wham_anomalies": wham_anomalies,
        "decision": {
            "all_route_directory_names_cover_0_to_15_by_0_3A": lineage[
                "decision"
            ]["all_routes_have_51_windows_0_to_15_by_0_3A"],
            "all_named_100_frame_files_actually_contain_100_frames": not bool(
                frame_failures
            ),
            "frame_failure_routes_are_exactly_fructofuranose_run2_and_run3": (
                frame_failure_keys
                == {"fructofuranose/run2", "fructofuranose/run3"}
            ),
            "fructofuranose_run2_has_one_and_ten_frame_files": (
                furanose_run2 is not None
                and furanose_run2["unique_frame_counts"] == [1, 10]
            ),
            "fructofuranose_run3_has_only_ten_frame_files": (
                furanose_run3 is not None
                and furanose_run3["unique_frame_counts"] == [10]
            ),
            "fructofuranose_run2_one_frame_centres_A": run2_one_frame_centres,
            "fructofuranose_run2_ten_frame_centres_A": run2_ten_frame_centres,
            "one_frame_visualization_file_matches_wham_short_window_center": (
                same_center_as_wham_short_window
            ),
            "all_center_mismatches_are_galactose_pathB_extra_2_25A": (
                all_center_failures_are_pathb_225
            ),
            "galactose_pathB_trajectory_archive_has_51_windows": all(
                row["trajectory_window_count"] == 51 for row in center_failures
            ),
            "galactose_pathB_wham_has_52_windows_including_2_25A": all(
                row["wham_window_count"] == 52
                and row["wham_only_centres_A"] == [2.25]
                for row in center_failures
            ),
            "galactose_run4a_run5a_complete_routes_are_exact_duplicates": lineage[
                "decision"
            ]["galactose_run4a_and_run5a_are_exact_duplicate_routes"],
            "galactose_run1b_run2b_have_duplicate_0_0_and_0_3_files": (
                set(lineage["within_route_duplicate_groups_by_route"])
                == {"galactose/run1b", "galactose/run2b"}
                and all(
                    groups[0]["centres_A"] == [0.0, 0.3]
                    for groups in lineage[
                        "within_route_duplicate_groups_by_route"
                    ].values()
                )
            ),
            "fourteen_route_directories_reduce_to_thirteen_unique_content_signatures": lineage[
                "decision"
            ][
                "fourteen_route_directories_reduce_to_thirteen_unique_content_signatures"
            ],
            "umbrella_archive_contains_full_wham_histogram_frame_series": False,
            "authoritative_route_lineage_file_present": False,
            "histogram_generation_reproducible_from_us_archive": False,
            "pmf_recalculation_from_us_trajectories_authorized": False,
            "parent_pag_calibration_authorized": False,
            "candidate_ranking_authorized": False,
            "next_action": (
                "request the missing full-resolution trajectories or exact frame-selection map, "
                "the galactose 2.25 A window lineage, corrected fructofuranose route files, "
                "and authoritative route-to-WHAM provenance before histogram regeneration"
            ),
        },
        "limitations": [
            "Frame counts are inferred from valid 183-atom Scripps BINPOS record sizes.",
            "A filename ending in 100-frames is not treated as evidence of the actual record count.",
            "The archive-to-WHAM route mapping remains path-based rather than author-supplied.",
            "The presence of the same anomalous centre in a short visualization file and short WHAM histogram does not prove direct file lineage.",
            "No trajectory coordinate, histogram, PMF, work, affinity, or candidate calculation is performed."
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
        "--umbrella-inventory",
        type=Path,
        default=Path("research/data/gluhut-us-archive-inventory.json"),
    )
    parser.add_argument(
        "--wham-reproduction",
        type=Path,
        default=Path("research/data/gluhut-all-sugars-wham-reproduction.json"),
    )
    parser.add_argument(
        "--route-lineage",
        type=Path,
        default=Path("research/data/gluhut-us-route-lineage.json"),
    )
    parser.add_argument(
        "--output-json",
        type=Path,
        default=Path("research/data/gluhut-us-lineage-anomalies.json"),
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    report = run(
        args.umbrella_inventory,
        args.wham_reproduction,
        args.route_lineage,
        args.output_json,
    )
    print(
        json.dumps(
            {
                "summary": report["summary"],
                "frame_count_failures": report["frame_count_failures"],
                "center_mismatches": report["center_mismatches"],
                "decision": report["decision"],
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
