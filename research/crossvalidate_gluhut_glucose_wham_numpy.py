#!/usr/bin/env python3
"""Independently cross-check deposited GluHUT glucose PMFs with NumPy WHAM.

The calculation uses the deposited reaction-coordinate time series and a clean
implementation of the standard one-dimensional WHAM fixed-point equations. It
must agree with the pinned Grossfield WHAM reproduction before this numerical
implementation is used for hold-out sugars or candidate comparisons.

This is a public GluHUT host-glucose calibration. It is not NNC2215, PaG,
MC-Q1, insulin signalling, efficacy, safety, or a drug result.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import statistics
import zipfile
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath

import numpy as np

ARCHIVE_MD5 = "dbc6723b8adcb85deb8230f63b6a51da"
ARCHIVE_SHA256 = "e0e7b7c04326b8fe86cf71e2193315f34fb865168cf11117b3caf9253119c6eb"
K_B_KCAL_MOL_K = 0.00198720425864083
TEMPERATURE_K = 300.0
BETA = 1.0 / (K_B_KCAL_MOL_K * TEMPERATURE_K)
FORCE_CONSTANT = 100.0
MINIMUM_A = 0.0
MAXIMUM_A = 15.0
BINS = 400
TOLERANCE = 1e-6
MAX_ITERATIONS = 200_000
BULK_MINIMUM_A = 13.0
BULK_MAXIMUM_A = 15.0
EXPECTED_GROSSFIELD = {
    "US-1": -13.794472169811321,
    "US-2": -12.79014541509434,
    "US-3": -11.760403301886793,
}
WINDOW_RE = re.compile(
    r"^WHAM_PMF/histogram-collection-glucose-"
    r"(?P<route>US-[123])/"
    r"(?P<center>\d+-\d+)-test-dist_vs_t\.dat$",
    re.IGNORECASE,
)


class CrosscheckError(RuntimeError):
    pass


def digest(path: Path, algorithm: str) -> str:
    value = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def safe_member(name: str) -> None:
    if "\\" in name:
        raise CrosscheckError(f"backslash in ZIP member: {name!r}")
    path = PurePosixPath(name)
    if path.is_absolute() or any(part in {"", ".", ".."} for part in path.parts):
        raise CrosscheckError(f"unsafe ZIP member: {name!r}")


def center_from_token(token: str) -> float:
    if not re.fullmatch(r"\d+-\d+", token):
        raise CrosscheckError(f"invalid centre token: {token!r}")
    whole, fractional = token.split("-", 1)
    return float(f"{int(whole)}.{fractional}")


def parse_coordinate_payload(payload: bytes) -> np.ndarray:
    values: list[float] = []
    for line_number, raw in enumerate(payload.decode("utf-8").splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith(("#", "@")):
            continue
        numeric: list[float] = []
        for token in line.split():
            try:
                numeric.append(float(token))
            except ValueError:
                continue
        if not numeric or not math.isfinite(numeric[-1]):
            raise CrosscheckError(f"invalid coordinate line {line_number}")
        values.append(numeric[-1])
    if not values:
        raise CrosscheckError("coordinate file contains no numeric samples")
    return np.asarray(values, dtype=float)


def load_glucose_routes(archive_path: Path) -> dict[str, list[tuple[float, np.ndarray, str]]]:
    routes: dict[str, list[tuple[float, np.ndarray, str]]] = {}
    seen: set[str] = set()
    with zipfile.ZipFile(archive_path) as archive:
        bad = archive.testzip()
        if bad:
            raise CrosscheckError(f"ZIP CRC failure: {bad}")
        for info in archive.infolist():
            safe_member(info.filename)
            if info.filename in seen:
                raise CrosscheckError(f"duplicate ZIP member: {info.filename}")
            seen.add(info.filename)
            if info.flag_bits & 0x1:
                raise CrosscheckError(f"encrypted ZIP member: {info.filename}")
            if info.is_dir():
                continue
            match = WINDOW_RE.match(info.filename)
            if not match:
                continue
            route = match.group("route")
            routes.setdefault(route, []).append(
                (
                    center_from_token(match.group("center")),
                    parse_coordinate_payload(archive.read(info)),
                    info.filename,
                )
            )
    for route, windows in routes.items():
        windows.sort(key=lambda row: row[0])
        centers = [row[0] for row in windows]
        if len(centers) != 51:
            raise CrosscheckError(f"{route}: expected 51 windows, found {len(centers)}")
        if len(centers) != len(set(centers)):
            raise CrosscheckError(f"{route}: duplicate centres")
    if sorted(routes) != ["US-1", "US-2", "US-3"]:
        raise CrosscheckError(f"unexpected glucose routes: {sorted(routes)}")
    return routes


def logsumexp(values: np.ndarray, axis: int) -> np.ndarray:
    maximum = np.max(values, axis=axis, keepdims=True)
    shifted = values - maximum
    return np.squeeze(maximum + np.log(np.sum(np.exp(shifted), axis=axis, keepdims=True)), axis=axis)


def adjacent_bhattacharyya(counts: np.ndarray) -> np.ndarray:
    probabilities = counts / counts.sum(axis=1, keepdims=True)
    return np.asarray(
        [
            float(np.sum(np.sqrt(probabilities[index] * probabilities[index + 1])))
            for index in range(probabilities.shape[0] - 1)
        ],
        dtype=float,
    )


def solve_route(windows: list[tuple[float, np.ndarray, str]]) -> dict:
    edges = np.linspace(MINIMUM_A, MAXIMUM_A, BINS + 1)
    bin_centers = 0.5 * (edges[:-1] + edges[1:])
    counts = np.asarray(
        [np.histogram(values, bins=edges)[0] for _, values, _ in windows],
        dtype=float,
    )
    samples = counts.sum(axis=1)
    if np.any(samples <= 0):
        raise CrosscheckError("one or more windows have no samples inside the WHAM range")
    umbrella_centers = np.asarray([center for center, _, _ in windows])
    bias = 0.5 * FORCE_CONSTANT * (
        bin_centers[None, :] - umbrella_centers[:, None]
    ) ** 2
    total_counts = counts.sum(axis=0)
    occupied = total_counts > 0
    widths = np.diff(edges)
    offsets = np.zeros(len(windows), dtype=float)
    probability = np.zeros(BINS, dtype=float)
    maximum_change = math.inf
    converged = False

    for iteration in range(1, MAX_ITERATIONS + 1):
        denominator = logsumexp(
            np.log(samples)[:, None] + BETA * offsets[:, None] - BETA * bias,
            axis=0,
        )
        log_probability = np.full(BINS, -np.inf, dtype=float)
        log_probability[occupied] = (
            np.log(total_counts[occupied]) - denominator[occupied]
        )
        normalization = logsumexp(
            (log_probability[occupied] + np.log(widths[occupied]))[None, :],
            axis=1,
        )[0]
        log_probability[occupied] -= normalization
        probability[:] = 0.0
        probability[occupied] = np.exp(log_probability[occupied])
        log_partition = logsumexp(
            log_probability[None, occupied]
            - BETA * bias[:, occupied]
            + np.log(widths[occupied])[None, :],
            axis=1,
        )
        new_offsets = -log_partition / BETA
        new_offsets -= new_offsets[0]
        maximum_change = float(np.max(np.abs(new_offsets - offsets)))
        offsets = new_offsets
        if maximum_change < TOLERANCE:
            converged = True
            break

    if not converged:
        raise CrosscheckError(
            f"NumPy WHAM did not converge after {MAX_ITERATIONS} iterations"
        )
    pmf = np.full(BINS, np.nan)
    positive = probability > 0
    pmf[positive] = -np.log(probability[positive]) / BETA
    pmf[positive] -= float(np.nanmin(pmf[positive]))
    bulk = positive & (bin_centers >= BULK_MINIMUM_A) & (bin_centers <= BULK_MAXIMUM_A)
    if bulk.sum() < 2:
        raise CrosscheckError("insufficient finite bins in bulk region")
    binding = float(np.nanmin(pmf) - np.mean(pmf[bulk]))
    overlap = adjacent_bhattacharyya(counts)
    return {
        "window_count": len(windows),
        "frames_total_raw": int(sum(values.size for _, values, _ in windows)),
        "frames_inside_range": int(samples.sum()),
        "frames_per_window_raw": [int(values.size) for _, values, _ in windows],
        "frames_per_window_inside_range": samples.astype(int).tolist(),
        "centres_A": umbrella_centers.tolist(),
        "iterations": iteration,
        "maximum_offset_change_kcal_mol": maximum_change,
        "binding_free_energy_kcal_mol": binding,
        "minimum_adjacent_bhattacharyya": float(overlap.min()),
        "median_adjacent_bhattacharyya": float(np.median(overlap)),
        "pmf_profile": [
            {
                "coordinate_A": float(coordinate),
                "relative_free_energy_kcal_mol": (
                    float(value) if math.isfinite(value) else None
                ),
            }
            for coordinate, value in zip(bin_centers, pmf)
        ],
        "source_members": [member for _, _, member in windows],
    }


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser()
    result.add_argument("--archive", type=Path, required=True)
    result.add_argument(
        "--output-json",
        type=Path,
        default=Path("research/data/gluhut-glucose-numpy-wham-crosscheck.json"),
    )
    result.add_argument("--maximum-route-difference-kcal", type=float, default=0.10)
    return result


def main() -> int:
    args = parser().parse_args()
    if digest(args.archive, "md5") != ARCHIVE_MD5:
        raise CrosscheckError("archive MD5 mismatch")
    if digest(args.archive, "sha256") != ARCHIVE_SHA256:
        raise CrosscheckError("archive SHA-256 mismatch")
    routes = load_glucose_routes(args.archive)
    results: list[dict] = []
    for route in sorted(routes):
        result = solve_route(routes[route])
        expected = EXPECTED_GROSSFIELD[route]
        difference = result["binding_free_energy_kcal_mol"] - expected
        result.update(
            {
                "route": route,
                "grossfield_reference_kcal_mol": expected,
                "difference_from_grossfield_kcal_mol": difference,
                "within_route_tolerance": abs(difference)
                <= args.maximum_route_difference_kcal,
            }
        )
        results.append(result)
    numpy_values = [row["binding_free_energy_kcal_mol"] for row in results]
    grossfield_values = [row["grossfield_reference_kcal_mol"] for row in results]
    report = {
        "schema": "dia-health-one.gluhut-glucose-numpy-wham-crosscheck.v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "crosscheck_complete_pending_scientific_review",
        "scope": "public GluHUT host-D-glucose methodological calibration only",
        "not_target_claim": "not NNC2215, PaG, MC-Q1, signalling, efficacy, safety, or a drug result",
        "archive": {
            "filename": args.archive.name,
            "size_bytes": args.archive.stat().st_size,
            "md5": ARCHIVE_MD5,
            "sha256": ARCHIVE_SHA256,
        },
        "implementation": {
            "language": "Python/NumPy",
            "numpy_version": np.__version__,
            "temperature_K": TEMPERATURE_K,
            "force_constant_kcal_mol_A2": FORCE_CONSTANT,
            "range_A": [MINIMUM_A, MAXIMUM_A],
            "bins": BINS,
            "tolerance_kcal_mol": TOLERANCE,
            "maximum_iterations": MAX_ITERATIONS,
            "bulk_region_A": [BULK_MINIMUM_A, BULK_MAXIMUM_A],
        },
        "routes": results,
        "summary": {
            "numpy_mean_kcal_mol": statistics.mean(numpy_values),
            "numpy_sample_sd_kcal_mol": statistics.stdev(numpy_values),
            "grossfield_mean_kcal_mol": statistics.mean(grossfield_values),
            "grossfield_sample_sd_kcal_mol": statistics.stdev(grossfield_values),
            "maximum_absolute_route_difference_kcal_mol": max(
                abs(row["difference_from_grossfield_kcal_mol"])
                for row in results
            ),
            "all_routes_within_tolerance": all(
                row["within_route_tolerance"] for row in results
            ),
        },
        "decision": {
            "independent_wham_crosscheck_pass": all(
                row["within_route_tolerance"] for row in results
            ),
            "candidate_ranking_authorized": False,
            "reason": "hold-out sugar and sensitivity gates remain required",
        },
        "limitations": [
            "Agreement on glucose does not validate PaG or transfer to NNC2215.",
            "Histogram bin-edge conventions can create small implementation differences.",
            "Autocorrelation weighting and sensitivity to trajectory truncation remain separate gates.",
        ],
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({"summary": report["summary"], "decision": report["decision"]}, indent=2, sort_keys=True))
    return 0 if report["decision"]["independent_wham_crosscheck_pass"] else 3


if __name__ == "__main__":
    raise SystemExit(main())
