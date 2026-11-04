#!/usr/bin/env python3
"""Reprocess every deposited GluHUT sugar PMF with pinned Grossfield WHAM.

This is a public-data calibration and hold-out analysis. It does not simulate
NNC2215, PaG, MC-Q1, insulin-receptor signalling, efficacy, safety, or a drug.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import statistics
import subprocess
import tempfile
import zipfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath

ARCHIVE_MD5" "dbc6723b8adcb85deb8230f63b6a51da"
ARCHIVE_SHA256 = "e0e7b7c04326b8fe86cf71e2193315f34fb865168cf11117b3caf9253119c6eb"
WHAM_SOURCE_COMMIT = "43bdd6fcd1916ae2b169cb52ac326b7eea567abb"
WHAM_VERSION = "2.0.11"
WINDOW_RE = re.compile(
    r"^WHAM_PMF/histogram-collection-"
    r"(?P<system>glucose|galactose|fructopyranose|fructofuranose)-"
    r"(?P<route>US-(:[123]|routeA))/"
    r"(?P<center>\d+-\d+)-test-dist_vs_t\.dat$",
    re.IGNORECASE,
)
EXPECTED_GROUPS = {
    "glucose": {"routes": ["US-1", "US-2", "US-3"], "published_mean": -12.5, "published_sd": 0.9},
    "galactose-pathB": {"routes": ["US-1", "US-2", "US-3"], "published_mean": -8.1, "published_sd": 1.9},
    "galactose-pathA": {"routes": ["US-routeA"], "published_mean": -6.9, "published_sd": 0.3},
    "fructopyranose": {"routes": ["US-1", "US-2", "US-3"], "published_mean": -5.0, "published_sd": 1.3},
    "fructofuranose": {"routes": ["US-1", "US-2", "US-3"], "published_mean": -7.7, "published_sd": 0.5},
}
EXPECTED_RANKING = [
    "glucose",
    "galactose-pathB",
    "fructofuranose",
    "galactose-pathA",
    "fructopyranose",
]


class AnalysisError(RuntimeError):
    pass



@dataclass(frozen=True)
class Window:
    archive_path: str
    center: float
    size_bytes: int
    frame_count: int


def digest(path: Path, algorithm: str) -> str:
    value = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def center_from_token(token: str) -> float:
    if not re.fullmatch(r\"\d+-\d+\"", token):
        raise AnalysisError(f"invalid deposited umbrella-centre token: {token!r}")
    whole, fractional = token.split("-", 1)
    return float(f"{int(whole)}.{fractional}")


def safe_member(name: str) -> None:
    if "\\" in name:
        raise AnalysisError(f"backslash in ZIP member: {name!r}")
    path = PurePosixPath(name)
    if path.is_absolute() or any(part in {"", ".", ".."} for part in path.parts):
        raise AnalysisError(f"unsafe ZIP member: {name!r}")


def frame_count(payload: bytes) -> int:
    return sum(
        1
        for line in payload.splitlines()
        if line.strip() and not line.lstrip().startswith((b"#", b"@"))
    )


def archive_inventory(archive: zipfile.ZipFile) -> tuple[dict[tuple[str, str], list[Window]], list[dict]]:
    groups: dict[tuple[str, str], list[Window]] = {}
    unmatched: list[dict] = []
    seen: set[str] = set()
    for info in archive.infolist():
        safe_member(info.filename)
        if info.filename in seen:
            raise AnalysisError(f"duplicate ZIP member: {info.filename}")
        seen.add(info.filename)
        if info.flag_bits & 0x1:
            raise AnalysisError(f"encrypted ZIP member: {info.filename}")
        if info.is_dir():
            continue
        match = WINDOW_RE.match(info.filename)
        if not match:
            unmatched.append({"path": info.filename, "size_bytes": info.file_size})
            continue
        payload = archive.read(info)
        window = Window(
            archive_path=info.filename,
            center=center_from_token(match.group("center")),
            size_bytes=info.file_size,
            frame_count=frame_count(payload),
        )
        key = (match.group("system").lower(), match.group("route"))
        groups.setdefault(key, []).append(window)
    for key, windows in groups.items():
        windows.sort(key=lambda row: row.center)
        centers = [row.center for row in windows]
        if len(centers) != len(set(centers)):
            raise AnalysisError(f"duplicate umbrella centres in {key}")
    return groups, unmatched


def pmf_from_file(path: Path) -> list[tuple[float, float]]:
    values: list[tuple[float, float]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) < 2:
            continue
        coordinate, free_energy = float(fields[0]), float(fields[1])
        if math.isfinite(coordinate) and math.isfinite(free_energy):
            values.append((coordinate, free_energy))
    if len(values) != 400:
        raise AnalysisError(f"expected 400 finite PMF bins in {path}; found {len(values)}")
    minimum = min(value for _, value in values)
    return [(coordinate, value - minimum) for coordinate, value in values]


def binding_energy(profile: list[tuple[float, float]]) -> float:
    bulk = [value for coordinate, value in profile if 13.0 <= coordinate <= 15.0]
    if len(bulk) < 2:
        raise AnalysisError("fewer than two bulk-reference bins")
    return min(value for _, value in profile) - statistics.mean(bulk)


def route_group(system: str, route: str) -> str:
    if system == "galactose":
        return "galactose-pathA" if route == "US-routeA" else "galactose-pathB"
    return system


def run_route(
    archive: zipfile.ZipFile,
    wham_binary: Path,
    system: str,
    route: str,
    windows: list[Window],
    work_root: Path,
) -> dict:
    route_dir = work_root / f"{system}-{route}"
    route_dir.mkdir(parents=True, exist_ok=True)
    metadata_lines: list[str] = []
    frame_counts: list[int] = []
    for window in windows:
        basename = Path(window.archive_path).name
        payload = archive.read(window.archive_path)
        (route_dir / basename).write_bytes(payload)
        metadata_lines.append(f"basename} {window.center:.8f} 100.0\n")
        frame_counts.append(window.frame_count)
    metadata = route_dir / "metadata.dat"
    metadata.write_text("".join(metadata_lines), encoding="utf-8")
    output = route_dir / "pmf.dat"
    process = subprocess.run(
        [
            str(wham_binary),
            "0.0",
            "15.0",
            "400",
            "0.000001",
            "300",
            "0",
            metadata.name,
            output.name,
        ],
        cwd=route_dir,
        capture_output=True,
        text=True,
        timeout=3600,
        check=False,
    )
    if process.returncode != 0:
        raise AnalysisError(
            f"Grossfield WHAM failed for {system} {route}: "
            f"{process.stderr[-2000:]} or {process.stdout[-2000:]}"
        )
    iterations = [
        int(value) for value in re.findall(r"#Iteration\s+(\d+):", process.stdout)
    ]
    profile = pmf_from_file(output)
    sizes = [window.size_bytes for window in windows]
    median_size = statistics.median(sizes)
    anomalies = [
        {
            "path": window.archive_path,
            "center_A": window.center,
            "size_bytes": window.size_bytes,
            "frame_count": window.frame_count,
            "size_fraction_of_route_median": window.size_bytes / median_size,
        }
        for window in windows
        if window.size_bytes < 0.5 * median_size or window.frame_count < 0.5 * statistics.median(frame_counts)
    ]
    return {
        "system": system,
        "route": route,
        "group": route_group(system, route),
        "window_count": len(windows),
        "centres_A": [window.center for window in windows],
        "minimum_spacing_A": min(
            b.center - a.center for a, b in zip(windows, windows[1:])
        )
        if len(windows) > 1
        else None,
        "maximum_spacing_A": max(
            b.center - a.center for a, b in zip(windows, windows[1:])
        )
        if len(windows) > 1
        else None,
        "frames_total": sum(frame_counts),
        "frames_minimum_window": min(frame_counts),
        "frames_median_window": statistics.median(frame_counts),
        "frames_maximum_window": max(frame_counts),
        "size_minimum_bytes": min(sizes),
        "size_median_bytes": median_size,
        "size_maximum_bytes": max(sizes),
        "anomalous_windows": anomalies,
        "wham_convergence_iteration": max(iterations) if iterations else None,
        "binding_free_energy_kcal_mol": binding_energy(profile),
        "pmf_profile": [
            {"coordinate_A": coordinate, "relative_free_energy_kcal_mol": energy}
            for coordinate, energy in profile
        ],
    }


def summarize(routes: list[dict]) -> list[dict]:
    summaries: list[dict] = []
    for group, expected in EXPECTED_GROUPS.items():
        rows = [row for row in routes if row["group"] == group]
        expected_routes = expected["routes"]
        observed_routes = [row["route"] for row in rows]
        if observed_routes != expected_routes:
            raise AnalysisError(
                f"{group}: routes {observed_routes} do not match expected {expected_routes}"
            )
        values = [row["binding_free_energy_kcal_mol"] for row in rows]
        mean = statistics.mean(values)
        sd = statistics.stdev(values) if len(values) > 1 else None
        summaries.append(
            {
                "group": group,
                "route_count": len(rows),
                "route_ids": observed_routes,
                "route_values_kcal_mol": values,
                "mean_kcal_mol": mean,
                "sample_sd_kcal_mol": sd,
                "published_mean_kcal_mol": expected["published_mean"],
                "published_sd_kcal_mol": expected["published_sd"],
                "difference_from_published_mean_kcal_mol": mean
                - expected["published_mean"],
                "mean_within_published_one_sd": abs(mean - expected["published_mean"])
                <= expected["published_sd"],
            }
        )
    return summaries


def write_summary_csv(path: Path, summaries: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "group",
        "route_count",
        "mean_kcal_mol",
        "sample_sd_kcal_mol",
        "published_mean_kcal_mol",
        "published_sd_kcal_mol",
        "difference_from_published_mean_kcal_mol",
        "mean_within_published_one_sd",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in summaries:
            writer.writerow({key: row[key] for key in fields})


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser()
    result.add_argument("--archive", type=Path, required=True)
    result.add_argument("--wham-binary", type=Path, required=True)
    result.add_argument(
        "--output-json",
        type=Path,
        default=Path("research/data/gluhut-all-sugars-wham-reproduction.json"),
    )
    result.add_argument(
        "--output-csv",
        type=Path,
        default=Path("research/data/gluhut-all-sugars-wham-summary.csv"),
    )
    return result


def main() -> int:
    args = parser().parse_args()
    if digest(args.archive, "md5") != ARCHIVE_MD5:
        raise AnalysisError("WHAM_PMF.zip MD5 does not match Zenodo")
    if digest(args.archive, "sha256") != ARCHIVE_SHA256:
        raise AnalysisError("WHAM_PMF.zip SHA-256 does not match the prior verified archive")
    if not args.wham_binary.is_file():
        raise AnalysisError(f"WHAM binary not found: {args.wham_binary}")
    with zipfile.ZipFile(args.archive) as archive:
        bad = archive.testzip()
        if bad:
            raise AnalysisError(f"ZIP CRC failure: {bad}")
        groups, unmatched = archive_inventory(archive)
        with tempfile.TemporaryDirectory(prefix="dho-all-sugars-wham-") as temporary:
            routes = [
                run_route(archive, args.wham_binary, system, route, windows, Path(temporary))
                for (system, route), windows in sorted(groups.items())
            ]
    summaries = summarize(routes)
    ranking = [
        row["group"]
        for row in sorted(summaries, key=lambda row: row["mean_kcal_mol"])
    ]
    report = {
        "schema": "dia-health-one.gluhut-all-sugars-wham-reproduction.v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "analysis_complete_pending_scientific_review",
        "scope": "public GluHUT sugar host-guest methodological calibration only",
        "not_target_claim": "not NNC2215, PaG, MC-Q1, insulin signalling, efficacy, safety, or a drug result",
        "archive": {
            "filename": args.archive.name,
            "size_bytes": args.archive.stat().st_size,
            "md5": ARCHIVE_MD5,
            "sha256": ARCHIVE_SHA256,
        },
        "wham": {
            "name": "Grossfield WHAM",
            "version": WHAM_VERSION,
            "source_commit": WHAM_SOURCE_COMMIT,
            "binary_sha256": digest(args.wham_binary, "sha256"),
            "range_A": [0.0, 15.0],
            "bins": 400,
            "tolerance": 1e-6,
            "temperature_K": 300,
            "periodic": False,
            "metadata_force_constant_kcal_mol_A2": 100.0,
            "bulk_region_A": [13.0, 15.0],
        },
        "routes": routes,
        "group_summaries": summaries,
        "ranking": {
            "observed_strongest_to_weakest": ranking,
            "published_strongest_to_weakest": EXPECTED_RANKING,
            "qualitative_order_matches": ranking == EXPECTED_RANKING,
        },
        "unmatched_archive_members": unmatched,
        "anomalies": [
            anomaly
            for route in routes
            for anomaly in route["anomalous_windows"]
        ],
        "decision": {
            "glucose_reproduction_present": any(row["group"] == "glucose" for row in summaries),
            "holdout_groups_present": all(
                any(row["group"] == group for row in summaries)
                for group in ("galactose-pathB", "fructopyranose", "fructofuranose")
            ),
            "candidate_ranking_authorized": False,
            "reason": "independent NumPy cross-validation and sensitivity review remain required",
        },
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_summary_csv(args.output_csv, summaries)
    print(json.dumps({
        "ranking": report["ranking"],
        "summaries": summaries,
        "anomaly_count": len(report["anomalies"]),
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
