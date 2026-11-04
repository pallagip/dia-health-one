#!/usr/bin/env python3
"""Reconcile umbrella-trajectory routes/windows with deposited WHAM inputs.

This consumes the deterministic umbrella archive inventory and the independently
reproduced WHAM route record. It checks route/window coverage, frame-count
reduction, exact duplicate trajectory content, and the bounded path-based mapping
between umbrella route directories and deposited WHAM routes.

It does not reconstruct the omitted full-resolution trajectory frames, prove
SMD-to-umbrella lineage, regenerate histograms, calculate a new PMF, or authorize
PaG/MC-Q1 ranking.
"""
from __future__ import annotations

import argparse
import collections
import hashlib
import json
import re
from pathlib import Path
from typing import Any

NATOM = 183
BINPOS_HEADER_BYTES = 4
BINPOS_RECORD_BYTES = 4 + 3 * NATOM * 4
EXPECTED_CENTRES_TENTHS = tuple(range(0, 151, 3))
EXPECTED_CENTRES_A = tuple(value / 10.0 for value in EXPECTED_CENTRES_TENTHS)

TRAJECTORY_RE = re.compile(
    r"^US_Trajectories/"
    r"(?P<system>glucose|galactose|fructopyranose|fructofuranose)/"
    r"run\s*(?P<run>\d+)\s*(?:\((?P<label>[abAB])\))?/"
    r"prod-(?P<whole>\d+)-(?P<fraction>\d+)-100-frames\.binpos$",
    re.IGNORECASE,
)


class LineageError(RuntimeError):
    pass


def canonical_hash(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
        ).encode("utf-8")
    ).hexdigest()


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_center_tenths(whole: str, fraction: str) -> int:
    if len(fraction) != 1:
        raise LineageError(
            f"expected one decimal digit in deposited window token: {whole}-{fraction}"
        )
    return 10 * int(whole) + int(fraction)


def normalize_route_id(run: int, label: str | None) -> str:
    suffix = label.lower() if label else ""
    return f"run{run}{suffix}"


def wham_target(system: str, run: int, label: str | None) -> tuple[str, str]:
    label = label.lower() if label else None
    if system == "galactose":
        if label == "b" and run in {1, 2, 3}:
            return system, f"US-{run}"
        if label == "a" and run in {4, 5}:
            return system, "US-routeA"
        raise LineageError(
            f"unsupported galactose route mapping: run={run}, label={label!r}"
        )
    if label is not None:
        raise LineageError(
            f"unexpected route label for {system}: run={run}, label={label!r}"
        )
    if run not in {1, 2, 3}:
        raise LineageError(f"unsupported {system} route number: {run}")
    return system, f"US-{run}"


def infer_binpos_frames(size_bytes: int) -> int:
    data_bytes = size_bytes - BINPOS_HEADER_BYTES
    if data_bytes <= 0 or data_bytes % BINPOS_RECORD_BYTES:
        raise LineageError(
            f"BINPOS size {size_bytes} does not match {NATOM}-atom records "
            f"of {BINPOS_RECORD_BYTES} bytes"
        )
    return data_bytes // BINPOS_RECORD_BYTES


def trajectory_rows(inventory: dict[str, Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for member in inventory["members"]:
        if member["kind"] != "file" or member["extension"] != ".binpos":
            continue
        match = TRAJECTORY_RE.fullmatch(member["path"])
        if not match:
            raise LineageError(
                f"umbrella trajectory path does not match expected grammar: {member['path']}"
            )
        system = match.group("system").lower()
        run = int(match.group("run"))
        label = match.group("label")
        center_tenths = parse_center_tenths(
            match.group("whole"), match.group("fraction")
        )
        mapped_system, mapped_route = wham_target(system, run, label)
        rows.append(
            {
                "path": member["path"],
                "sha256": member["sha256"],
                "size_bytes": member["uncompressed_size_bytes"],
                "frame_count": infer_binpos_frames(
                    member["uncompressed_size_bytes"]
                ),
                "system": system,
                "run": run,
                "label": label.lower() if label else None,
                "route_id": normalize_route_id(run, label),
                "route_key": f"{system}/{normalize_route_id(run, label)}",
                "center_tenths": center_tenths,
                "center_A": center_tenths / 10.0,
                "mapped_wham_system": mapped_system,
                "mapped_wham_route": mapped_route,
                "mapped_wham_key": f"{mapped_system}/{mapped_route}",
            }
        )
    rows.sort(key=lambda row: (row["system"], row["route_id"], row["center_tenths"]))
    return rows


def wham_routes(wham: dict[str, Any]) -> dict[str, dict[str, Any]]:
    result: dict[str, dict[str, Any]] = {}
    for route in wham["routes"]:
        key = f"{route['system']}/{route['route']}"
        if key in result:
            raise LineageError(f"duplicate WHAM route key: {key}")
        result[key] = route
    return result


def route_summaries(
    rows: list[dict[str, Any]],
    wham_by_key: dict[str, dict[str, Any]],
) -> list[dict[str, Any]]:
    grouped: dict[str, list[dict[str, Any]]] = collections.defaultdict(list)
    for row in rows:
        grouped[row["route_key"]].append(row)

    summaries: list[dict[str, Any]] = []
    for route_key, route_rows in sorted(grouped.items()):
        route_rows.sort(key=lambda row: row["center_tenths"])
        centers = [row["center_tenths"] for row in route_rows]
        if len(centers) != len(set(centers)):
            raise LineageError(f"duplicate center tokens in route {route_key}")
        missing = sorted(set(EXPECTED_CENTRES_TENTHS) - set(centers))
        unexpected = sorted(set(centers) - set(EXPECTED_CENTRES_TENTHS))
        frames = sorted({row["frame_count"] for row in route_rows})
        mapped_keys = sorted({row["mapped_wham_key"] for row in route_rows})
        if len(mapped_keys) != 1:
            raise LineageError(
                f"route {route_key} maps to multiple WHAM routes: {mapped_keys}"
            )
        mapped_key = mapped_keys[0]
        if mapped_key not in wham_by_key:
            raise LineageError(
                f"route {route_key} maps to absent WHAM route {mapped_key}"
            )
        wham_route = wham_by_key[mapped_key]
        wham_centers_tenths = [
            int(round(float(center) * 10))
            for center in wham_route["centres_A"]
        ]
        if len(wham_centers_tenths) != len(set(wham_centers_tenths)):
            raise LineageError(f"duplicate centers in WHAM route {mapped_key}")

        hash_to_rows: dict[str, list[dict[str, Any]]] = collections.defaultdict(list)
        for row in route_rows:
            hash_to_rows[row["sha256"]].append(row)
        within_route_duplicates = [
            {
                "sha256": digest,
                "centres_A": [row["center_A"] for row in duplicate_rows],
                "paths": [row["path"] for row in duplicate_rows],
                "file_count": len(duplicate_rows),
            }
            for digest, duplicate_rows in sorted(hash_to_rows.items())
            if len(duplicate_rows) > 1
        ]

        signature_payload = [
            {
                "center_tenths": row["center_tenths"],
                "sha256": row["sha256"],
                "frame_count": row["frame_count"],
            }
            for row in route_rows
        ]
        median_histogram_frames = float(wham_route["frames_median_window"])
        trajectory_frame_count = frames[0] if len(frames) == 1 else None
        downsampling_ratio = (
            median_histogram_frames / trajectory_frame_count
            if trajectory_frame_count
            else None
        )

        summaries.append(
            {
                "route_key": route_key,
                "system": route_rows[0]["system"],
                "route_id": route_rows[0]["route_id"],
                "label": route_rows[0]["label"],
                "mapped_wham_key": mapped_key,
                "trajectory_file_count": len(route_rows),
                "expected_window_count": len(EXPECTED_CENTRES_TENTHS),
                "window_centres_A": [row["center_A"] for row in route_rows],
                "missing_window_centres_A": [value / 10.0 for value in missing],
                "unexpected_window_centres_A": [
                    value / 10.0 for value in unexpected
                ],
                "complete_0_to_15_by_0_3A_coverage": (
                    not missing
                    and not unexpected
                    and centers == list(EXPECTED_CENTRES_TENTHS)
                ),
                "trajectory_frames_per_window_values": frames,
                "trajectory_frames_total": sum(
                    row["frame_count"] for row in route_rows
                ),
                "unique_trajectory_content_count": len(hash_to_rows),
                "within_route_duplicate_groups": within_route_duplicates,
                "within_route_duplicate_group_count": len(
                    within_route_duplicates
                ),
                "route_content_signature_sha256": canonical_hash(
                    signature_payload
                ),
                "wham_window_count": int(wham_route["window_count"]),
                "wham_centres_match_trajectory_centres": (
                    wham_centers_tenths == centers
                ),
                "wham_frames_minimum_window": int(
                    wham_route["frames_minimum_window"]
                ),
                "wham_frames_median_window": median_histogram_frames,
                "wham_frames_maximum_window": int(
                    wham_route["frames_maximum_window"]
                ),
                "median_histogram_to_trajectory_frame_ratio": (
                    downsampling_ratio
                ),
            }
        )
    return summaries


def route_signature_groups(
    summaries: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    grouped: dict[str, list[dict[str, Any]]] = collections.defaultdict(list)
    for summary in summaries:
        grouped[summary["route_content_signature_sha256"]].append(summary)
    return [
        {
            "route_content_signature_sha256": signature,
            "route_keys": [row["route_key"] for row in rows],
            "mapped_wham_keys": sorted(
                {row["mapped_wham_key"] for row in rows}
            ),
            "route_count": len(rows),
            "trajectory_file_count_per_route": sorted(
                {row["trajectory_file_count"] for row in rows}
            ),
        }
        for signature, rows in sorted(grouped.items())
        if len(rows) > 1
    ]


def duplicate_trajectory_groups(
    rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    grouped: dict[str, list[dict[str, Any]]] = collections.defaultdict(list)
    for row in rows:
        grouped[row["sha256"]].append(row)
    result = []
    for digest, duplicate_rows in sorted(grouped.items()):
        if len(duplicate_rows) < 2:
            continue
        route_keys = sorted({row["route_key"] for row in duplicate_rows})
        centers = sorted({row["center_A"] for row in duplicate_rows})
        result.append(
            {
                "sha256": digest,
                "file_count": len(duplicate_rows),
                "route_keys": route_keys,
                "route_count": len(route_keys),
                "centres_A": centers,
                "paths": [row["path"] for row in duplicate_rows],
                "same_center_across_routes": (
                    len(route_keys) > 1 and len(centers) == 1
                ),
                "different_centers_within_or_across_routes": len(centers) > 1,
            }
        )
    return result


def run(
    inventory_path: Path,
    wham_path: Path,
    output_path: Path,
) -> dict[str, Any]:
    inventory = json.loads(inventory_path.read_text(encoding="utf-8"))
    wham = json.loads(wham_path.read_text(encoding="utf-8"))
    rows = trajectory_rows(inventory)
    wham_by_key = wham_routes(wham)
    summaries = route_summaries(rows, wham_by_key)
    signature_duplicates = route_signature_groups(summaries)
    trajectory_duplicates = duplicate_trajectory_groups(rows)

    raw_route_count = len(summaries)
    unique_route_signature_count = len(
        {row["route_content_signature_sha256"] for row in summaries}
    )
    mapped_wham_keys = {row["mapped_wham_key"] for row in summaries}
    complete_routes = all(
        row["complete_0_to_15_by_0_3A_coverage"] for row in summaries
    )
    all_centers_match = all(
        row["wham_centres_match_trajectory_centres"] for row in summaries
    )
    all_100_frames = all(
        row["trajectory_frames_per_window_values"] == [100]
        for row in summaries
    )

    by_key = {row["route_key"]: row for row in summaries}
    expected_route_keys = {
        "glucose/run1",
        "glucose/run2",
        "glucose/run3",
        "fructopyranose/run1",
        "fructopyranose/run2",
        "fructopyranose/run3",
        "fructofuranose/run1",
        "fructofuranose/run2",
        "fructofuranose/run3",
        "galactose/run1b",
        "galactose/run2b",
        "galactose/run3b",
        "galactose/run4a",
        "galactose/run5a",
    }
    observed_route_keys = set(by_key)

    galactose_a_exact_duplicate = (
        "galactose/run4a" in by_key
        and "galactose/run5a" in by_key
        and by_key["galactose/run4a"]["route_content_signature_sha256"]
        == by_key["galactose/run5a"]["route_content_signature_sha256"]
    )

    within_route_duplicate_routes = {
        row["route_key"]: row["within_route_duplicate_groups"]
        for row in summaries
        if row["within_route_duplicate_groups"]
    }

    effective_wham_mapping = collections.defaultdict(list)
    for summary in summaries:
        effective_wham_mapping[summary["mapped_wham_key"]].append(
            summary["route_key"]
        )

    report = {
        "schema": "dia-health-one.gluhut-us-route-lineage.v1",
        "status": (
            "route_window_coverage_reconciled_full_resolution_histogram_lineage_missing"
        ),
        "scope": "public GluHUT umbrella trajectory archive versus deposited WHAM routes",
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
        },
        "binpos_layout": {
            "natom": NATOM,
            "header_bytes": BINPOS_HEADER_BYTES,
            "record_bytes": BINPOS_RECORD_BYTES,
        },
        "trajectory_file_count": len(rows),
        "raw_route_directory_count": raw_route_count,
        "unique_complete_route_content_signature_count": (
            unique_route_signature_count
        ),
        "deposited_wham_route_count": len(wham_by_key),
        "mapped_wham_route_count": len(mapped_wham_keys),
        "expected_route_keys": sorted(expected_route_keys),
        "observed_route_keys": sorted(observed_route_keys),
        "route_summaries": summaries,
        "duplicate_complete_route_groups": signature_duplicates,
        "duplicate_trajectory_content_groups": trajectory_duplicates,
        "within_route_duplicate_groups_by_route": within_route_duplicate_routes,
        "wham_route_to_umbrella_directories": {
            key: sorted(values)
            for key, values in sorted(effective_wham_mapping.items())
        },
        "decision": {
            "all_714_trajectory_files_parsed": len(rows) == 714,
            "expected_fourteen_route_directories_present": (
                observed_route_keys == expected_route_keys
                and raw_route_count == 14
            ),
            "all_routes_have_51_windows_0_to_15_by_0_3A": complete_routes,
            "all_trajectory_files_have_100_recorded_frames": all_100_frames,
            "all_route_centres_match_mapped_wham_centres": all_centers_match,
            "fourteen_route_directories_reduce_to_thirteen_unique_content_signatures": (
                raw_route_count == 14
                and unique_route_signature_count == 13
            ),
            "thirteen_unique_route_signatures_match_thirteen_wham_routes": (
                unique_route_signature_count == len(wham_by_key) == 13
                and len(mapped_wham_keys) == 13
            ),
            "galactose_run4a_and_run5a_are_exact_duplicate_routes": (
                galactose_a_exact_duplicate
            ),
            "galactose_routeA_has_two_archive_copies_but_one_wham_route": (
                sorted(
                    effective_wham_mapping.get(
                        "galactose/US-routeA", []
                    )
                )
                == ["galactose/run4a", "galactose/run5a"]
            ),
            "within_route_duplicate_window_content_detected": bool(
                within_route_duplicate_routes
            ),
            "umbrella_archive_contains_full_wham_histogram_frame_series": False,
            "median_histogram_to_visualization_frame_ratio_is_approximately_150": all(
                row["median_histogram_to_trajectory_frame_ratio"]
                is not None
                and 149.0
                <= row["median_histogram_to_trajectory_frame_ratio"]
                <= 151.0
                for row in summaries
            ),
            "authoritative_route_lineage_file_present": False,
            "path_based_route_mapping_is_complete_but_not_authoritative": (
                mapped_wham_keys == set(wham_by_key)
                and complete_routes
                and all_centers_match
            ),
            "histogram_generation_reproducible_from_us_archive": False,
            "pmf_recalculation_from_us_trajectories_authorized": False,
            "parent_pag_calibration_authorized": False,
            "candidate_ranking_authorized": False,
            "next_action": (
                "compare duplicate window hashes with deposited histogram content, request "
                "the omitted full-resolution umbrella trajectories or frame-selection map, "
                "and obtain authoritative route-to-WHAM lineage before histogram regeneration"
            ),
        },
        "limitations": [
            "The route-to-WHAM mapping is inferred from system, run number, A/B label and complete centre coverage; the deposit contains no authoritative lineage manifest.",
            "Each umbrella BINPOS contains 100 recorded frames, whereas deposited WHAM histograms usually contain approximately 15,001 samples; no sampling/downselection map is present.",
            "Byte-identical 100-frame visualization files do not prove that the unavailable full-resolution trajectories were identical.",
            "Matching route centres and route counts does not prove that the deposited histograms were generated from these downsampled files.",
            "No coordinate, histogram, PMF, work, affinity or candidate calculation is performed."
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
        default=Path(
            "research/data/gluhut-all-sugars-wham-reproduction.json"
        ),
    )
    parser.add_argument(
        "--output-json",
        type=Path,
        default=Path("research/data/gluhut-us-route-lineage.json"),
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    report = run(
        args.umbrella_inventory,
        args.wham_reproduction,
        args.output_json,
    )
    print(
        json.dumps(
            {
                "trajectory_file_count": report["trajectory_file_count"],
                "raw_route_directory_count": report[
                    "raw_route_directory_count"
                ],
                "unique_complete_route_content_signature_count": report[
                    "unique_complete_route_content_signature_count"
                ],
                "deposited_wham_route_count": report[
                    "deposited_wham_route_count"
                ],
                "duplicate_complete_route_groups": report[
                    "duplicate_complete_route_groups"
                ],
                "within_route_duplicate_groups_by_route": report[
                    "within_route_duplicate_groups_by_route"
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
