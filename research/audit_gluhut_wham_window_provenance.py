#!/usr/bin/env python3
"""Audit GluHUT WHAM reaction-coordinate time-series provenance.

The deposited WHAM archive can be reanalysed deterministically, but its
relationship to the downsampled umbrella trajectory archive is incomplete. This
audit focuses on three material questions:

1. Are the galactose path-B 0.0/0.3 Å WHAM series duplicates, mirroring the
   duplicate visualization trajectories?
2. Are the WHAM-only 2.25 Å galactose windows independent deposited series,
   copies, or simple sample-wise interpolation of 2.1/2.4 Å neighbours?
3. Is the 205-sample fructofuranose US-2 12.6 Å series a clean prefix on the
   ordinary time grid, or a different sampling process?

The audit hashes and summarizes every deposited window but retains no full time
series in the output. It does not regenerate histograms, calculate a new PMF or
affinity, or authorize PaG/MC-Q1 work.
"""
from __future__ import annotations

import argparse
import collections
import hashlib
import json
import math
import re
import statistics
import stat
import struct
import zipfile
from pathlib import Path, PurePosixPath
from typing import Any

ARCHIVE_MD5 = "dbc6723b8adcb85deb8230f63b6a51da"
ARCHIVE_SHA256 = "e0e7b7c04326b8fe86cf71e2193315f34fb865168cf11117b3caf9253119c6eb"
ARCHIVE_SIZE_BYTES = 41560830
CHUNK_SIZE = 1024 * 1024
FLOAT_TOLERANCE = 1.0e-10

WINDOW_RE = re.compile(
    r"^WHAM_PMF/histogram-collection-"
    r"(?P<system>glucose|galactose|fructopyranose|fructofuranose)-"
    r"(?P<route>US-(?:[123]|routeA))/"
    r"(?P<center>\d+-\d+)-test-dist_vs_t\.dat$",
    re.IGNORECASE,
)

TARGET_KEYS = {
    *(f"galactose/US-{route}/{center:.2f}" for route in range(1, 4) for center in (0.0, 0.3, 2.1, 2.25, 2.4)),
    *(f"fructofuranose/US-2/{center:.2f}" for center in (12.3, 12.6, 12.9)),
}


class ProvenanceError(RuntimeError):
    pass


def digest_file(path: Path, algorithm: str) -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(CHUNK_SIZE), b""):
            digest.update(block)
    return digest.hexdigest()


def safe_member(name: str) -> PurePosixPath:
    if "\\" in name:
        raise ProvenanceError(f"backslash in ZIP member: {name!r}")
    path = PurePosixPath(name)
    if path.is_absolute() or any(part in {"", ".", ".."} for part in path.parts):
        raise ProvenanceError(f"unsafe ZIP member: {name!r}")
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


def center_from_token(token: str) -> float:
    if not re.fullmatch(r"\d+-\d+", token):
        raise ProvenanceError(f"invalid centre token: {token!r}")
    whole, fraction = token.split("-", 1)
    return float(f"{int(whole)}.{fraction}")


def window_key(system: str, route: str, center: float) -> str:
    return f"{system.lower()}/{route}/{center:.2f}"


def hash_float_series(values: list[float]) -> str:
    digest = hashlib.sha256()
    for value in values:
        digest.update(struct.pack("!d", value))
    return digest.hexdigest()


def parse_series(payload: bytes, path: str) -> tuple[list[float], list[float]]:
    try:
        text = payload.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise ProvenanceError(f"non-UTF8 WHAM series: {path}") from exc
    times: list[float] = []
    coordinates: list[float] = []
    for line_number, raw in enumerate(text.splitlines(), start=1):
        line = raw.strip()
        if not line or line.startswith(("#", "@")):
            continue
        fields = line.split()
        if len(fields) < 2:
            raise ProvenanceError(
                f"{path}:{line_number}: expected at least two columns"
            )
        try:
            time_value = float(fields[0])
            coordinate = float(fields[1])
        except ValueError as exc:
            raise ProvenanceError(
                f"{path}:{line_number}: non-numeric data {fields[:2]}"
            ) from exc
        if not math.isfinite(time_value) or not math.isfinite(coordinate):
            raise ProvenanceError(
                f"{path}:{line_number}: non-finite value"
            )
        times.append(time_value)
        coordinates.append(coordinate)
    if not times:
        raise ProvenanceError(f"no numeric samples in {path}")
    if any(right <= left for left, right in zip(times, times[1:])):
        raise ProvenanceError(f"time column is not strictly increasing in {path}")
    return times, coordinates


def population_std(values: list[float]) -> float:
    if len(values) < 2:
        return 0.0
    return statistics.pstdev(values)


def step_summary(times: list[float]) -> dict[str, Any]:
    steps = [right - left for left, right in zip(times, times[1:])]
    if not steps:
        return {
            "count": 0,
            "minimum": None,
            "median": None,
            "maximum": None,
            "all_equal_within_tolerance": True,
            "series_sha256": hash_float_series([]),
        }
    median = statistics.median(steps)
    return {
        "count": len(steps),
        "minimum": min(steps),
        "median": median,
        "maximum": max(steps),
        "all_equal_within_tolerance": all(
            abs(step - median) <= FLOAT_TOLERANCE for step in steps
        ),
        "series_sha256": hash_float_series(steps),
    }


def series_summary(
    payload: bytes,
    path: str,
    system: str,
    route: str,
    center: float,
) -> tuple[dict[str, Any], tuple[list[float], list[float]] | None]:
    times, coordinates = parse_series(payload, path)
    steps = step_summary(times)
    summary = {
        "path": path,
        "system": system,
        "route": route,
        "center_A": center,
        "key": window_key(system, route, center),
        "size_bytes": len(payload),
        "file_sha256": hashlib.sha256(payload).hexdigest(),
        "sample_count": len(times),
        "first_time": times[0],
        "last_time": times[-1],
        "duration": times[-1] - times[0],
        "time_series_sha256": hash_float_series(times),
        "coordinate_series_sha256": hash_float_series(coordinates),
        "time_step": steps,
        "coordinate": {
            "first": coordinates[0],
            "last": coordinates[-1],
            "minimum": min(coordinates),
            "maximum": max(coordinates),
            "mean": statistics.fmean(coordinates),
            "population_sd": population_std(coordinates),
        },
    }
    retained = (
        (times, coordinates)
        if summary["key"] in TARGET_KEYS
        else None
    )
    return summary, retained


def pearson(left: list[float], right: list[float]) -> float | None:
    if len(left) != len(right) or len(left) < 2:
        return None
    mean_left = statistics.fmean(left)
    mean_right = statistics.fmean(right)
    dev_left = [value - mean_left for value in left]
    dev_right = [value - mean_right for value in right]
    norm_left = math.sqrt(sum(value * value for value in dev_left))
    norm_right = math.sqrt(sum(value * value for value in dev_right))
    if norm_left <= 0.0 or norm_right <= 0.0:
        return None
    return sum(a * b for a, b in zip(dev_left, dev_right)) / (
        norm_left * norm_right
    )


def compare_series(
    left: tuple[list[float], list[float]],
    right: tuple[list[float], list[float]],
) -> dict[str, Any]:
    left_times, left_coordinates = left
    right_times, right_coordinates = right
    same_length = len(left_times) == len(right_times)
    exact_time_grid = same_length and left_times == right_times
    if not same_length:
        return {
            "same_length": False,
            "left_count": len(left_times),
            "right_count": len(right_times),
            "exact_time_grid": False,
            "exact_coordinate_series": False,
            "coordinate_pearson": None,
            "coordinate_mean_absolute_difference": None,
            "coordinate_rmse": None,
            "coordinate_maximum_absolute_difference": None,
        }
    differences = [
        left - right
        for left, right in zip(left_coordinates, right_coordinates)
    ]
    return {
        "same_length": True,
        "left_count": len(left_times),
        "right_count": len(right_times),
        "exact_time_grid": exact_time_grid,
        "exact_coordinate_series": left_coordinates == right_coordinates,
        "coordinate_pearson": pearson(left_coordinates, right_coordinates),
        "coordinate_mean_absolute_difference": statistics.fmean(
            abs(value) for value in differences
        ),
        "coordinate_rmse": math.sqrt(
            statistics.fmean(value * value for value in differences)
        ),
        "coordinate_maximum_absolute_difference": max(
            abs(value) for value in differences
        ),
    }


def prefix_comparison(
    short: tuple[list[float], list[float]],
    long: tuple[list[float], list[float]],
) -> dict[str, Any]:
    short_times, short_coordinates = short
    long_times, long_coordinates = long
    prefix_available = len(short_times) <= len(long_times)
    time_prefix_equal = (
        prefix_available and short_times == long_times[: len(short_times)]
    )
    coordinate_prefix_equal = (
        prefix_available
        and short_coordinates == long_coordinates[: len(short_coordinates)]
    )
    return {
        "short_count": len(short_times),
        "long_count": len(long_times),
        "prefix_available": prefix_available,
        "time_grid_is_exact_prefix": time_prefix_equal,
        "coordinate_series_is_exact_prefix": coordinate_prefix_equal,
        "sample_fraction": len(short_times) / len(long_times),
        "duration_fraction": (
            (short_times[-1] - short_times[0])
            / (long_times[-1] - long_times[0])
            if len(long_times) > 1 and long_times[-1] != long_times[0]
            else None
        ),
        "time_step_median_short": step_summary(short_times)["median"],
        "time_step_median_long": step_summary(long_times)["median"],
    }


def interpolation_comparison(
    lower: tuple[list[float], list[float]],
    middle: tuple[list[float], list[float]],
    upper: tuple[list[float], list[float]],
) -> dict[str, Any]:
    lower_times, lower_coordinates = lower
    middle_times, middle_coordinates = middle
    upper_times, upper_coordinates = upper
    same_length = (
        len(lower_times) == len(middle_times) == len(upper_times)
    )
    same_grid = same_length and lower_times == middle_times == upper_times
    if not same_grid:
        return {
            "same_length": same_length,
            "same_time_grid": False,
            "samplewise_arithmetic_mean_exact": False,
            "samplewise_arithmetic_mean_mean_absolute_difference": None,
            "samplewise_arithmetic_mean_rmse": None,
            "samplewise_arithmetic_mean_maximum_absolute_difference": None,
        }
    differences = [
        middle_value - 0.5 * (lower_value + upper_value)
        for lower_value, middle_value, upper_value in zip(
            lower_coordinates,
            middle_coordinates,
            upper_coordinates,
        )
    ]
    return {
        "same_length": True,
        "same_time_grid": True,
        "samplewise_arithmetic_mean_exact": all(
            value == 0.0 for value in differences
        ),
        "samplewise_arithmetic_mean_mean_absolute_difference": statistics.fmean(
            abs(value) for value in differences
        ),
        "samplewise_arithmetic_mean_rmse": math.sqrt(
            statistics.fmean(value * value for value in differences)
        ),
        "samplewise_arithmetic_mean_maximum_absolute_difference": max(
            abs(value) for value in differences
        ),
    }


def audit(archive_path: Path, output_path: Path) -> dict[str, Any]:
    if archive_path.stat().st_size != ARCHIVE_SIZE_BYTES:
        raise ProvenanceError(
            f"size mismatch: expected {ARCHIVE_SIZE_BYTES}, "
            f"got {archive_path.stat().st_size}"
        )
    archive_md5 = digest_file(archive_path, "md5")
    archive_sha256 = digest_file(archive_path, "sha256")
    if archive_md5 != ARCHIVE_MD5:
        raise ProvenanceError(
            f"MD5 mismatch: expected {ARCHIVE_MD5}, got {archive_md5}"
        )
    if archive_sha256 != ARCHIVE_SHA256:
        raise ProvenanceError(
            f"SHA-256 mismatch: expected {ARCHIVE_SHA256}, got {archive_sha256}"
        )

    summaries: list[dict[str, Any]] = []
    retained: dict[str, tuple[list[float], list[float]]] = {}
    seen_paths: set[str] = set()

    with zipfile.ZipFile(archive_path) as archive:
        bad = archive.testzip()
        if bad:
            raise ProvenanceError(f"ZIP CRC failure: {bad}")
        for info in archive.infolist():
            safe_member(info.filename)
            if info.filename in seen_paths:
                raise ProvenanceError(f"duplicate ZIP member path: {info.filename}")
            seen_paths.add(info.filename)
            if info.flag_bits & 0x1:
                raise ProvenanceError(f"encrypted ZIP member: {info.filename}")
            kind = member_kind(info)
            if kind in {"symlink", "special"}:
                raise ProvenanceError(
                    f"unsupported {kind}: {info.filename}"
                )
            if kind != "file":
                continue
            match = WINDOW_RE.fullmatch(info.filename)
            if not match:
                continue
            system = match.group("system").lower()
            route = match.group("route")
            center = center_from_token(match.group("center"))
            summary, series = series_summary(
                archive.read(info),
                info.filename,
                system,
                route,
                center,
            )
            summaries.append(summary)
            if series is not None:
                retained[summary["key"]] = series

    summaries.sort(
        key=lambda row: (row["system"], row["route"], row["center_A"])
    )
    if not summaries:
        raise ProvenanceError("no WHAM time-series windows found")
    missing_targets = sorted(TARGET_KEYS - set(retained))
    if missing_targets:
        raise ProvenanceError(f"missing target windows: {missing_targets}")

    hash_groups: dict[str, list[dict[str, Any]]] = collections.defaultdict(list)
    time_grid_groups: dict[tuple[int, str], list[dict[str, Any]]] = (
        collections.defaultdict(list)
    )
    sample_count_distribution: collections.Counter[int] = collections.Counter()
    for summary in summaries:
        hash_groups[summary["file_sha256"]].append(summary)
        time_grid_groups[
            (summary["sample_count"], summary["time_series_sha256"])
        ].append(summary)
        sample_count_distribution[summary["sample_count"]] += 1

    exact_duplicate_groups = [
        {
            "file_sha256": digest,
            "file_count": len(rows),
            "keys": [row["key"] for row in rows],
            "paths": [row["path"] for row in rows],
        }
        for digest, rows in sorted(hash_groups.items())
        if len(rows) > 1
    ]

    route_summaries: list[dict[str, Any]] = []
    routes: dict[tuple[str, str], list[dict[str, Any]]] = collections.defaultdict(list)
    for summary in summaries:
        routes[(summary["system"], summary["route"])].append(summary)
    for (system, route), rows in sorted(routes.items()):
        rows.sort(key=lambda row: row["center_A"])
        route_median = statistics.median(
            row["sample_count"] for row in rows
        )
        anomalies = [
            {
                "key": row["key"],
                "center_A": row["center_A"],
                "sample_count": row["sample_count"],
                "duration": row["duration"],
                "file_sha256": row["file_sha256"],
            }
            for row in rows
            if row["sample_count"] < 0.5 * route_median
        ]
        route_summaries.append(
            {
                "system": system,
                "route": route,
                "window_count": len(rows),
                "sample_count_values": sorted(
                    {row["sample_count"] for row in rows}
                ),
                "median_sample_count": route_median,
                "time_grid_count": len(
                    {
                        (row["sample_count"], row["time_series_sha256"])
                        for row in rows
                    }
                ),
                "short_window_anomalies": anomalies,
            }
        )

    galactose_checks = []
    for route_number in range(1, 4):
        route = f"US-{route_number}"
        key = lambda center: f"galactose/{route}/{center:.2f}"
        zero_vs_point3 = compare_series(
            retained[key(0.0)], retained[key(0.3)]
        )
        lower_vs_middle = compare_series(
            retained[key(2.1)], retained[key(2.25)]
        )
        middle_vs_upper = compare_series(
            retained[key(2.25)], retained[key(2.4)]
        )
        interpolation = interpolation_comparison(
            retained[key(2.1)],
            retained[key(2.25)],
            retained[key(2.4)],
        )
        middle_summary = next(
            row for row in summaries if row["key"] == key(2.25)
        )
        duplicate_of = [
            row["key"]
            for row in summaries
            if row["file_sha256"] == middle_summary["file_sha256"]
            and row["key"] != middle_summary["key"]
        ]
        galactose_checks.append(
            {
                "route": route,
                "zero_vs_0_3": zero_vs_point3,
                "2_1_vs_2_25": lower_vs_middle,
                "2_25_vs_2_4": middle_vs_upper,
                "2_25_vs_samplewise_mean_of_2_1_and_2_4": interpolation,
                "2_25_sample_count": middle_summary["sample_count"],
                "2_25_file_sha256": middle_summary["file_sha256"],
                "2_25_exact_duplicate_of_other_windows": sorted(duplicate_of),
            }
        )

    furanose_short = retained["fructofuranose/US-2/12.60"]
    furanose_lower = retained["fructofuranose/US-2/12.30"]
    furanose_upper = retained["fructofuranose/US-2/12.90"]
    short_summary = next(
        row
        for row in summaries
        if row["key"] == "fructofuranose/US-2/12.60"
    )
    furanose_check = {
        "short_window_key": short_summary["key"],
        "sample_count": short_summary["sample_count"],
        "duration": short_summary["duration"],
        "time_step": short_summary["time_step"],
        "prefix_of_12_3": prefix_comparison(
            furanose_short, furanose_lower
        ),
        "prefix_of_12_9": prefix_comparison(
            furanose_short, furanose_upper
        ),
        "exact_duplicate_of_other_windows": sorted(
            row["key"]
            for row in summaries
            if row["file_sha256"] == short_summary["file_sha256"]
            and row["key"] != short_summary["key"]
        ),
    }

    gal_225_independent = all(
        not check["2_25_exact_duplicate_of_other_windows"]
        and not check["2_25_vs_samplewise_mean_of_2_1_and_2_4"][
            "samplewise_arithmetic_mean_exact"
        ]
        for check in galactose_checks
    )
    gal_00_03_hist_duplicate = any(
        check["zero_vs_0_3"]["exact_coordinate_series"]
        and check["zero_vs_0_3"]["exact_time_grid"]
        for check in galactose_checks
    )
    furanose_time_prefix = (
        furanose_check["prefix_of_12_3"]["time_grid_is_exact_prefix"]
        and furanose_check["prefix_of_12_9"]["time_grid_is_exact_prefix"]
    )

    report = {
        "schema": "dia-health-one.gluhut-wham-window-provenance.v1",
        "status": "wham_window_series_audited_trajectory_lineage_still_unresolved",
        "scope": "deposited public GluHUT WHAM reaction-coordinate time series",
        "not_target_claim": (
            "not regenerated histograms, not a new PMF or affinity, not NNC2215, "
            "PaG, MC-Q1, efficacy, safety, or a drug result"
        ),
        "archive": {
            "filename": archive_path.name,
            "size_bytes": archive_path.stat().st_size,
            "md5": archive_md5,
            "sha256": archive_sha256,
            "crc_verified": True,
        },
        "window_file_count": len(summaries),
        "sample_count_distribution": {
            str(count): files
            for count, files in sorted(sample_count_distribution.items())
        },
        "exact_duplicate_file_groups": exact_duplicate_groups,
        "time_grid_signature_count": len(time_grid_groups),
        "route_summaries": route_summaries,
        "galactose_pathB_checks": galactose_checks,
        "fructofuranose_US2_12_6_check": furanose_check,
        "decision": {
            "all_window_files_parsed": True,
            "all_time_columns_strictly_increasing": True,
            "galactose_0_0_and_0_3_wham_series_are_exact_duplicates_in_any_route": (
                gal_00_03_hist_duplicate
            ),
            "galactose_2_25_windows_are_exact_duplicates_or_exact_samplewise_interpolations": (
                not gal_225_independent
            ),
            "galactose_2_25_windows_are_distinct_deposited_series_under_tested_exact_checks": (
                gal_225_independent
            ),
            "fructofuranose_US2_12_6_has_205_samples": (
                short_summary["sample_count"] == 205
            ),
            "fructofuranose_US2_12_6_time_grid_is_exact_prefix_of_both_neighbours": (
                furanose_time_prefix
            ),
            "fructofuranose_US2_12_6_coordinate_series_is_exact_prefix_of_either_neighbour": (
                furanose_check["prefix_of_12_3"][
                    "coordinate_series_is_exact_prefix"
                ]
                or furanose_check["prefix_of_12_9"][
                    "coordinate_series_is_exact_prefix"
                ]
            ),
            "authoritative_trajectory_to_wham_lineage_established": False,
            "histogram_generation_reproducible_from_public_trajectory_archive": False,
            "pmf_recalculation_from_public_trajectory_archive_authorized": False,
            "parent_pag_calibration_authorized": False,
            "candidate_ranking_authorized": False,
            "next_action": (
                "request the source 2.25 A galactose trajectories and full-resolution "
                "fructofuranose 12.6 A series, plus the exact trajectory-to-WHAM manifest; "
                "do not infer lineage from matching centers or time grids alone"
            ),
        },
        "limitations": [
            "Exact and sample-wise arithmetic-mean tests do not exclude more complex derivation of a series.",
            "A shared time-grid prefix indicates compatible sampling cadence, not shared source trajectory.",
            "The audit does not compare molecular coordinates or regenerate reaction coordinates from full-resolution trajectories.",
            "No WHAM profile or free energy is recalculated."
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
        default=Path("research/data/gluhut-wham-window-provenance.json"),
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    report = audit(args.archive, args.output_json)
    print(
        json.dumps(
            {
                "archive": report["archive"],
                "window_file_count": report["window_file_count"],
                "sample_count_distribution": report[
                    "sample_count_distribution"
                ],
                "exact_duplicate_file_groups": report[
                    "exact_duplicate_file_groups"
                ],
                "galactose_pathB_checks": report[
                    "galactose_pathB_checks"
                ],
                "fructofuranose_US2_12_6_check": report[
                    "fructofuranose_US2_12_6_check"
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
