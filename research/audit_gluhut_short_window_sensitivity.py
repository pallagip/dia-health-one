#!/usr/bin/env python3
"""Sensitivity analysis for the short GluHUT fructofuranose US-2 window.

The deposited 12.6 A reaction-coordinate file contains only 205 frames, versus
roughly 15,001 in neighbouring windows. This script quantifies the effect of
excluding that window, discarding early trajectory fractions, and changing WHAM
bin count. It also records local histogram overlap.

This is public GluHUT host-sugar method calibration. It is not NNC2215, PaG,
MC-Q1, efficacy, safety, or a drug result.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import statistics
import subprocess
import tempfile
import zipfile
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath

import numpy as np

ARCHIVE_MD5 = "dbc6723b8adcb85deb8230f63b6a51da"
ARCHIVE_SHA256 = "e0e7b7c04326b8fe86cf71e2193315f34fb865168cf11117b3caf9253119c6eb"
WHAM_SOURCE_COMMIT = "43bdd6fcd1916ae2b169cb52ac326b7eea567abb"
WHAM_VERSION = "2.0.11"
TARGET_PREFIX = "WHAM_PMF/histogram-collection-fructofuranose-US-2/"
WINDOW_RE = re.compile(
    r"^WHAM_PMF/histogram-collection-fructofuranose-US-2/"
    r"(?P<center>\d+-\d+)-test-dist_vs_t\.dat$"
)
SHORT_CENTER_A = 12.6
BASELINE_REFERENCE_KCAL_MOL = -7.521242245283019
PRE_REGISTERED_MAXIMUM_CHANGE_KCAL_MOL = 0.50


class SensitivityError(RuntimeError):
    pass


def digest(path: Path, algorithm: str) -> str:
    value = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def safe_member(name: str) -> None:
    if "\\" in name:
        raise SensitivityError(f"backslash in ZIP member: {name!r}")
    path = PurePosixPath(name)
    if path.is_absolute() or any(part in {"", ".", ".."} for part in path.parts):
        raise SensitivityError(f"unsafe ZIP member: {name!r}")


def center_from_token(token: str) -> float:
    if not re.fullmatch(r"\d+-\d+", token):
        raise SensitivityError(f"invalid centre token: {token!r}")
    whole, fractional = token.split("-", 1)
    return float(f"{int(whole)}.{fractional}")


def split_payload(payload: bytes) -> tuple[list[bytes], list[bytes]]:
    headers: list[bytes] = []
    rows: list[bytes] = []
    for raw in payload.splitlines(keepends=True):
        stripped = raw.strip()
        if not stripped or stripped.startswith((b"#", b"@")):
            headers.append(raw)
        else:
            rows.append(raw)
    if not rows:
        raise SensitivityError("window contains no numeric rows")
    return headers, rows


def last_numeric_column(rows: list[bytes]) -> np.ndarray:
    values: list[float] = []
    for line_number, raw in enumerate(rows, 1):
        numeric: list[float] = []
        for token in raw.decode("utf-8").split():
            try:
                numeric.append(float(token))
            except ValueError:
                continue
        if not numeric or not math.isfinite(numeric[-1]):
            raise SensitivityError(f"invalid numeric row {line_number}")
        values.append(numeric[-1])
    return np.asarray(values, dtype=float)


def load_windows(archive_path: Path) -> list[dict]:
    windows: list[dict] = []
    seen: set[str] = set()
    with zipfile.ZipFile(archive_path) as archive:
        bad = archive.testzip()
        if bad:
            raise SensitivityError(f"ZIP CRC failure: {bad}")
        for info in archive.infolist():
            safe_member(info.filename)
            if info.filename in seen:
                raise SensitivityError(f"duplicate ZIP member: {info.filename}")
            seen.add(info.filename)
            if info.flag_bits & 0x1:
                raise SensitivityError(f"encrypted ZIP member: {info.filename}")
            if info.is_dir():
                continue
            match = WINDOW_RE.match(info.filename)
            if not match:
                continue
            payload = archive.read(info)
            headers, rows = split_payload(payload)
            windows.append(
                {
                    "archive_path": info.filename,
                    "center_A": center_from_token(match.group("center")),
                    "headers": headers,
                    "rows": rows,
                    "coordinates": last_numeric_column(rows),
                    "size_bytes": info.file_size,
                }
            )
    windows.sort(key=lambda row: row["center_A"])
    if len(windows) != 51:
        raise SensitivityError(f"expected 51 target windows, found {len(windows)}")
    centers = [row["center_A"] for row in windows]
    if len(centers) != len(set(centers)):
        raise SensitivityError("duplicate target centres")
    short = [row for row in windows if math.isclose(row["center_A"], SHORT_CENTER_A)]
    if len(short) != 1 or len(short[0]["rows"]) != 205:
        raise SensitivityError("deposited 12.6 A short-window assumption changed")
    return windows


def bhattacharyya(left: np.ndarray, right: np.ndarray, bins: int = 400) -> float:
    edges = np.linspace(0.0, 15.0, bins + 1)
    left_counts = np.histogram(left, bins=edges)[0].astype(float)
    right_counts = np.histogram(right, bins=edges)[0].astype(float)
    if left_counts.sum() == 0 or right_counts.sum() == 0:
        raise SensitivityError("empty overlap histogram")
    left_probability = left_counts / left_counts.sum()
    right_probability = right_counts / right_counts.sum()
    return float(np.sum(np.sqrt(left_probability * right_probability)))


def pmf_from_file(path: Path, bins: int) -> list[tuple[float, float]]:
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
    if len(values) != bins:
        raise SensitivityError(
            f"expected {bins} finite PMF bins, found {len(values)}"
        )
    minimum = min(value for _, value in values)
    return [(coordinate, value - minimum) for coordinate, value in values]


def binding_energy(profile: list[tuple[float, float]]) -> float:
    bulk = [value for coordinate, value in profile if 13.0 <= coordinate <= 15.0]
    if len(bulk) < 2:
        raise SensitivityError("fewer than two bulk bins")
    return min(value for _, value in profile) - statistics.mean(bulk)


def write_variant_window(window: dict, destination: Path, discard_fraction: float) -> int:
    rows = window["rows"]
    discard_count = int(math.floor(len(rows) * discard_fraction))
    retained = rows[discard_count:]
    if not retained:
        raise SensitivityError("discard removed every frame")
    destination.write_bytes(b"".join(window["headers"] + retained))
    return len(retained)


def run_variant(
    wham_binary: Path,
    windows: list[dict],
    work_root: Path,
    *,
    name: str,
    bins: int,
    discard_fraction: float,
    exclude_short_window: bool,
) -> dict:
    variant_dir = work_root / name
    variant_dir.mkdir(parents=True, exist_ok=True)
    metadata_lines: list[str] = []
    frame_counts: list[int] = []
    selected = [
        row
        for row in windows
        if not (
            exclude_short_window
            and math.isclose(row["center_A"], SHORT_CENTER_A)
        )
    ]
    for window in selected:
        basename = Path(window["archive_path"]).name
        retained = write_variant_window(
            window, variant_dir / basename, discard_fraction
        )
        frame_counts.append(retained)
        metadata_lines.append(f"{basename} {window['center_A']:.8f} 100.0\n")
    metadata = variant_dir / "metadata.dat"
    metadata.write_text("".join(metadata_lines), encoding="utf-8")
    output = variant_dir / "pmf.dat"
    process = subprocess.run(
        [
            str(wham_binary),
            "0.0",
            "15.0",
            str(bins),
            "0.000001",
            "300",
            "0",
            metadata.name,
            output.name,
        ],
        cwd=variant_dir,
        capture_output=True,
        text=True,
        timeout=3600,
        check=False,
    )
    if process.returncode != 0:
        return {
            "name": name,
            "status": "wham_failed",
            "returncode": process.returncode,
            "stderr_tail": process.stderr[-2000:],
            "stdout_tail": process.stdout[-2000:],
            "bins": bins,
            "discard_fraction": discard_fraction,
            "excluded_short_window": exclude_short_window,
            "window_count": len(selected),
        }
    iterations = [
        int(value)
        for value in re.findall(r"#Iteration\s+(\d+):", process.stdout)
    ]
    profile = pmf_from_file(output, bins)
    return {
        "name": name,
        "status": "success",
        "bins": bins,
        "discard_fraction": discard_fraction,
        "excluded_short_window": exclude_short_window,
        "window_count": len(selected),
        "frames_total": sum(frame_counts),
        "frames_minimum_window": min(frame_counts),
        "frames_median_window": statistics.median(frame_counts),
        "frames_maximum_window": max(frame_counts),
        "wham_convergence_iteration": max(iterations) if iterations else None,
        "binding_free_energy_kcal_mol": binding_energy(profile),
    }


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser()
    result.add_argument("--archive", type=Path, required=True)
    result.add_argument("--wham-binary", type=Path, required=True)
    result.add_argument(
        "--output-json",
        type=Path,
        default=Path(
            "research/data/gluhut-fructofuranose-us2-sensitivity.json"
        ),
    )
    return result


def main() -> int:
    args = parser().parse_args()
    if digest(args.archive, "md5") != ARCHIVE_MD5:
        raise SensitivityError("archive MD5 mismatch")
    if digest(args.archive, "sha256") != ARCHIVE_SHA256:
        raise SensitivityError("archive SHA-256 mismatch")
    if not args.wham_binary.is_file():
        raise SensitivityError(f"WHAM binary not found: {args.wham_binary}")
    windows = load_windows(args.archive)
    by_center = {row["center_A"]: row for row in windows}
    local_overlap = {
        "12.3_to_12.6": bhattacharyya(
            by_center[12.3]["coordinates"], by_center[12.6]["coordinates"]
        ),
        "12.6_to_12.9": bhattacharyya(
            by_center[12.6]["coordinates"], by_center[12.9]["coordinates"]
        ),
        "12.3_to_12.9_without_short_window": bhattacharyya(
            by_center[12.3]["coordinates"], by_center[12.9]["coordinates"]
        ),
    }
    specifications = [
        ("baseline", 400, 0.0, False),
        ("exclude_short_12_6", 400, 0.0, True),
        ("discard_first_10_percent", 400, 0.10, False),
        ("discard_first_20_percent", 400, 0.20, False),
        ("bins_200", 200, 0.0, False),
        ("bins_800", 800, 0.0, False),
    ]
    with tempfile.TemporaryDirectory(prefix="dho-fructofuranose-sensitivity-") as temporary:
        variants = [
            run_variant(
                args.wham_binary,
                windows,
                Path(temporary),
                name=name,
                bins=bins,
                discard_fraction=discard,
                exclude_short_window=exclude,
            )
            for name, bins, discard, exclude in specifications
        ]
    baseline = next(row for row in variants if row["name"] == "baseline")
    if baseline["status"] != "success":
        raise SensitivityError("baseline WHAM failed")
    baseline_value = baseline["binding_free_energy_kcal_mol"]
    if abs(baseline_value - BASELINE_REFERENCE_KCAL_MOL) > 1e-9:
        raise SensitivityError(
            "baseline does not reproduce the merged all-sugar record"
        )
    for row in variants:
        if row["status"] == "success":
            row["difference_from_baseline_kcal_mol"] = (
                row["binding_free_energy_kcal_mol"] - baseline_value
            )
            row["within_pre_registered_0_5_kcal_mol"] = (
                abs(row["difference_from_baseline_kcal_mol"])
                <= PRE_REGISTERED_MAXIMUM_CHANGE_KCAL_MOL
            )
        else:
            row["difference_from_baseline_kcal_mol"] = None
            row["within_pre_registered_0_5_kcal_mol"] = False
    nonbaseline = [row for row in variants if row["name"] != "baseline"]
    report = {
        "schema": "dia-health-one.gluhut-fructofuranose-short-window-sensitivity.v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "sensitivity_complete_pending_scientific_review",
        "scope": "public GluHUT fructofuranose US-2 method sensitivity only",
        "not_target_claim": "not NNC2215, PaG, MC-Q1, efficacy, safety, or a drug result",
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
        },
        "short_window": {
            "archive_path": by_center[12.6]["archive_path"],
            "center_A": 12.6,
            "frame_count": len(by_center[12.6]["rows"]),
            "size_bytes": by_center[12.6]["size_bytes"],
            "local_bhattacharyya_overlap": local_overlap,
        },
        "pre_registered_maximum_absolute_change_kcal_mol": PRE_REGISTERED_MAXIMUM_CHANGE_KCAL_MOL,
        "variants": variants,
        "decision": {
            "all_sensitivity_variants_succeeded": all(
                row["status"] == "success" for row in variants
            ),
            "all_nonbaseline_variants_within_0_5_kcal_mol": all(
                row["within_pre_registered_0_5_kcal_mol"]
                for row in nonbaseline
            ),
            "maximum_absolute_change_kcal_mol": max(
                abs(row["difference_from_baseline_kcal_mol"])
                for row in nonbaseline
                if row["difference_from_baseline_kcal_mol"] is not None
            ),
            "candidate_ranking_authorized": False,
            "reason": "independent WHAM and PaG calibration remain required",
        },
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "short_window": report["short_window"],
                "variants": variants,
                "decision": report["decision"],
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
