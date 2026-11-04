#!/usr/bin/env python3
"""Publish umbrella-provenance results into the public status artifacts.

The updater reads the compact deterministic anomaly ledger, modifies only the
relevant machine-readable status fields, and adds one evidence link to the
human-readable status page. It is idempotent and refuses to publish if the
ledger's evidence boundary would permit PMF, PaG, or candidate ranking.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

TARGET_SCHEMA = "dia-health-one.public-confidence-status.v10"
RECORD_RELATIVE_PATH = "../GLUHUT_UMBRELLA_LINEAGE_AUDIT.md"
LINK_HTML = (
    '<a href="../GLUHUT_UMBRELLA_LINEAGE_AUDIT.md">'
    '<span>Umbrella lineage audit</span><span>↗</span></a>'
)


class PublishError(RuntimeError):
    pass


def require_false(decision: dict[str, Any], key: str) -> None:
    if decision.get(key) is not False:
        raise PublishError(f"expected {key}=false, got {decision.get(key)!r}")


def require_true(decision: dict[str, Any], key: str) -> None:
    if decision.get(key) is not True:
        raise PublishError(f"expected {key}=true, got {decision.get(key)!r}")


def update_status(
    status: dict[str, Any],
    anomaly: dict[str, Any],
) -> dict[str, Any]:
    decision = anomaly["decision"]
    summary = anomaly["summary"]

    for key in (
        "galactose_run4a_run5a_complete_routes_are_exact_duplicates",
        "galactose_run1b_run2b_have_duplicate_0_0_and_0_3_files",
        "all_center_mismatches_are_galactose_pathB_extra_2_25A",
        "fructofuranose_run2_has_one_and_ten_frame_files",
        "fructofuranose_run3_has_only_ten_frame_files",
        "one_frame_visualization_file_matches_wham_short_window_center",
        "fourteen_route_directories_reduce_to_thirteen_unique_content_signatures",
    ):
        require_true(decision, key)
    for key in (
        "all_named_100_frame_files_actually_contain_100_frames",
        "umbrella_archive_contains_full_wham_histogram_frame_series",
        "authoritative_route_lineage_file_present",
        "histogram_generation_reproducible_from_us_archive",
        "pmf_recalculation_from_us_trajectories_authorized",
        "parent_pag_calibration_authorized",
        "candidate_ranking_authorized",
    ):
        require_false(decision, key)

    result = json.loads(json.dumps(status))
    result["schema"] = TARGET_SCHEMA
    result["updated"] = "2026-11-04"
    calibration = result.setdefault("public_benchmark_calibration", {})
    calibration.update(
        {
            "umbrella_archive_size_bytes": anomaly["inputs"][
                "umbrella_inventory"
            ]["archive"]["size_bytes"],
            "umbrella_archive_md5": anomaly["inputs"][
                "umbrella_inventory"
            ]["archive"]["md5"],
            "umbrella_archive_sha256": anomaly["inputs"][
                "umbrella_inventory"
            ]["archive"]["sha256"],
            "umbrella_trajectory_file_count": summary[
                "trajectory_file_count"
            ],
            "umbrella_raw_route_directory_count": summary[
                "raw_route_directory_count"
            ],
            "umbrella_unique_complete_route_signature_count": summary[
                "unique_complete_route_content_signature_count"
            ],
            "umbrella_deposited_wham_route_count": summary[
                "deposited_wham_route_count"
            ],
            "umbrella_full_histogram_frame_series_available": False,
            "umbrella_authoritative_route_lineage_available": False,
            "umbrella_histogram_generation_reproducible": False,
            "umbrella_pmf_recalculation_from_public_trajectory_archive_authorized": False,
            "umbrella_galactose_run4a_run5a_exact_duplicate_routes": True,
            "umbrella_galactose_run1b_run2b_duplicate_0_0_0_3_files": True,
            "umbrella_galactose_pathB_trajectory_window_count": 51,
            "umbrella_galactose_pathB_wham_window_count": 52,
            "umbrella_galactose_pathB_wham_only_center_A": 2.25,
            "umbrella_fructofuranose_run2_frame_count_values": [1, 10],
            "umbrella_fructofuranose_run2_one_frame_center_A": 12.6,
            "umbrella_fructofuranose_run3_frame_count_values": [10],
            "umbrella_one_frame_12_6A_file_cooccurs_with_205_sample_wham_anomaly": True,
            "umbrella_direct_file_lineage_to_205_sample_wham_anomaly_established": False,
            "umbrella_lineage_anomaly_record": "../data/gluhut-us-lineage-anomalies.json",
            "umbrella_lineage_record": RECORD_RELATIVE_PATH,
        }
    )

    aws = result.setdefault("aws_execution", {})
    aws["current_gate"] = (
        "AWS-1 public analysis reproduced; host architecture and force convention "
        "reconciled; SMD and umbrella file provenance audited; public umbrella "
        "visualization trajectories cannot regenerate WHAM histograms; source "
        "parameters, full-resolution trajectories, route lineage, start coordinates "
        "and parent PaG calibration remain pending"
    )
    aws["next_action"] = (
        "obtain the full-resolution umbrella trajectories or exact frame-selection "
        "map, the galactose 2.25 A window lineage, corrected fructofuranose route "
        "files, RESP/GAFF/GLYCAM and source start coordinates; then preregister a "
        "source-matched parent glucose-versus-PaG calibration"
    )
    aws["candidate_screen_authorized"] = False

    external = result.setdefault("external_execution", {})
    external["author_umbrella_lineage_clarification_sent"] = True

    result["next_confidence_changes"] = [
        "obtain the original full-resolution umbrella trajectories or exact frame-selection and time map",
        "obtain the missing galactose path-B 2.25 A trajectory and route-to-WHAM lineage",
        "obtain corrected fructofuranose run2/run3 trajectory files and explain the 12.6 A short window",
        "obtain the intended independent galactose path-A route and correct duplicate 0.0/0.3 path-B files",
        "obtain RESP/GAFF/GLYCAM and source starting-coordinate provenance",
        "calibrate a source-matched parent glucose-versus-PaG calculation",
        "obtain provider acceptance and synthesis quotation",
        "measure authenticated free-binder glucose and PaG binding",
    ]
    return result


def update_html(html: str) -> str:
    if LINK_HTML in html:
        return html
    marker = '<a href="status.json"><span>Machine-readable status</span><span>↓</span></a>'
    if marker not in html:
        raise PublishError("could not find machine-readable status link marker")
    return html.replace(marker, LINK_HTML + "\n      " + marker, 1)


def run(
    status_path: Path,
    html_path: Path,
    anomaly_path: Path,
) -> dict[str, Any]:
    status = json.loads(status_path.read_text(encoding="utf-8"))
    anomaly = json.loads(anomaly_path.read_text(encoding="utf-8"))
    updated_status = update_status(status, anomaly)
    updated_html = update_html(html_path.read_text(encoding="utf-8"))
    status_path.write_text(
        json.dumps(updated_status, indent=2, sort_keys=False) + "\n",
        encoding="utf-8",
    )
    html_path.write_text(updated_html, encoding="utf-8")
    return updated_status


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--status-json",
        type=Path,
        default=Path("research/status/status.json"),
    )
    parser.add_argument(
        "--status-html",
        type=Path,
        default=Path("research/status/index.html"),
    )
    parser.add_argument(
        "--anomaly-json",
        type=Path,
        default=Path("research/data/gluhut-us-lineage-anomalies.json"),
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    updated = run(args.status_json, args.status_html, args.anomaly_json)
    print(
        json.dumps(
            {
                "schema": updated["schema"],
                "umbrella_lineage_record": updated[
                    "public_benchmark_calibration"
                ]["umbrella_lineage_record"],
                "aws_execution": updated["aws_execution"],
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
