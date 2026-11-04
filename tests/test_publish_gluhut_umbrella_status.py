from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "research" / "publish_gluhut_umbrella_status.py"
spec = importlib.util.spec_from_file_location("publisher", MODULE_PATH)
assert spec and spec.loader
module = importlib.util.module_from_spec(spec)
sys.modules["publisher"] = module
spec.loader.exec_module(module)


def anomaly_record() -> dict:
    return {
        "inputs": {
            "umbrella_inventory": {
                "archive": {
                    "size_bytes": 123756824,
                    "md5": "c460e6a1ee97ac531fc284d8f07e381e",
                    "sha256": "1" * 64,
                }
            }
        },
        "summary": {
            "trajectory_file_count": 714,
            "raw_route_directory_count": 14,
            "unique_complete_route_content_signature_count": 13,
            "deposited_wham_route_count": 13,
        },
        "decision": {
            "galactose_run4a_run5a_complete_routes_are_exact_duplicates": True,
            "galactose_run1b_run2b_have_duplicate_0_0_and_0_3_files": True,
            "all_center_mismatches_are_galactose_pathB_extra_2_25A": True,
            "fructofuranose_run2_has_one_and_ten_frame_files": True,
            "fructofuranose_run3_has_only_ten_frame_files": True,
            "one_frame_visualization_file_matches_wham_short_window_center": True,
            "fourteen_route_directories_reduce_to_thirteen_unique_content_signatures": True,
            "all_named_100_frame_files_actually_contain_100_frames": False,
            "umbrella_archive_contains_full_wham_histogram_frame_series": False,
            "authoritative_route_lineage_file_present": False,
            "histogram_generation_reproducible_from_us_archive": False,
            "pmf_recalculation_from_us_trajectories_authorized": False,
            "parent_pag_calibration_authorized": False,
            "candidate_ranking_authorized": False,
        },
    }


class UmbrellaStatusPublisherTests(unittest.TestCase):
    def test_update_status_preserves_existing_fields(self) -> None:
        status = {
            "schema": "old",
            "working_drug_discovered": False,
            "primary_candidate": {"id": "MC-Q1", "evidence_level": "E2"},
            "public_benchmark_calibration": {"existing": 1},
            "aws_execution": {"aws_credits_spent_usd": 0},
            "external_execution": {},
            "next_confidence_changes": ["old"],
        }
        updated = module.update_status(status, anomaly_record())
        self.assertEqual(updated["schema"], module.TARGET_SCHEMA)
        self.assertFalse(updated["working_drug_discovered"])
        self.assertEqual(updated["primary_candidate"]["id"], "MC-Q1")
        self.assertEqual(
            updated["public_benchmark_calibration"]["existing"],
            1,
        )
        calibration = updated["public_benchmark_calibration"]
        self.assertEqual(calibration["umbrella_trajectory_file_count"], 714)
        self.assertEqual(
            calibration["umbrella_unique_complete_route_signature_count"],
            13,
        )
        self.assertFalse(
            calibration["umbrella_histogram_generation_reproducible"]
        )
        self.assertTrue(
            calibration[
                "umbrella_one_frame_12_6A_file_cooccurs_with_205_sample_wham_anomaly"
            ]
        )
        self.assertFalse(
            calibration[
                "umbrella_direct_file_lineage_to_205_sample_wham_anomaly_established"
            ]
        )
        self.assertFalse(updated["aws_execution"]["candidate_screen_authorized"])
        self.assertEqual(updated["aws_execution"]["aws_credits_spent_usd"], 0)
        self.assertTrue(
            updated["external_execution"][
                "author_umbrella_lineage_clarification_sent"
            ]
        )

    def test_refuses_overclaiming_anomaly_record(self) -> None:
        status = {"public_benchmark_calibration": {}, "aws_execution": {}}
        anomaly = anomaly_record()
        anomaly["decision"]["candidate_ranking_authorized"] = True
        with self.assertRaises(module.PublishError):
            module.update_status(status, anomaly)

    def test_html_update_is_idempotent(self) -> None:
        html = (
            '<div class="links">\n'
            '  <a href="status.json"><span>Machine-readable status</span><span>↓</span></a>\n'
            '</div>\n'
        )
        first = module.update_html(html)
        second = module.update_html(first)
        self.assertEqual(first, second)
        self.assertEqual(first.count(module.LINK_HTML), 1)

    def test_full_file_update(self) -> None:
        status = {
            "schema": "old",
            "public_benchmark_calibration": {},
            "aws_execution": {},
            "external_execution": {},
        }
        html = (
            '<div class="links">\n'
            '<a href="status.json"><span>Machine-readable status</span><span>↓</span></a>\n'
            '</div>\n'
        )
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            status_path = root / "status.json"
            html_path = root / "index.html"
            anomaly_path = root / "anomaly.json"
            status_path.write_text(json.dumps(status), encoding="utf-8")
            html_path.write_text(html, encoding="utf-8")
            anomaly_path.write_text(
                json.dumps(anomaly_record()), encoding="utf-8"
            )
            result = module.run(status_path, html_path, anomaly_path)
            self.assertEqual(result["schema"], module.TARGET_SCHEMA)
            self.assertEqual(
                json.loads(status_path.read_text()),
                result,
            )
            self.assertIn(module.LINK_HTML, html_path.read_text())


if __name__ == "__main__":
    unittest.main()
