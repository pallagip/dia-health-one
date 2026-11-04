from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "research" / "summarize_gluhut_us_lineage_anomalies.py"
spec = importlib.util.spec_from_file_location("lineage_anomalies", MODULE_PATH)
assert spec and spec.loader
module = importlib.util.module_from_spec(spec)
sys.modules["lineage_anomalies"] = module
spec.loader.exec_module(module)

RECORD_BYTES = 4 + 3 * 183 * 4


def size_for_frames(frames: int) -> int:
    return 4 + frames * RECORD_BYTES


def member(
    path: str,
    sha: str,
    frames: int,
) -> dict:
    return {
        "path": path,
        "kind": "file",
        "extension": ".binpos",
        "sha256": sha,
        "uncompressed_size_bytes": size_for_frames(frames),
    }


def path(system: str, route: str, center: float) -> str:
    whole = int(center)
    fraction = int(round((center - whole) * 10))
    if system == "galactose":
        run = route[3]
        label = route[4]
        directory = f"run {run} ({label})"
    else:
        run = route[3]
        directory = f"run{run}"
    return (
        f"US_Trajectories/{system}/{directory}/"
        f"prod-{whole}-{fraction}-100-frames.binpos"
    )


class LineageAnomalyTests(unittest.TestCase):
    def test_mode_requires_unique_winner(self) -> None:
        self.assertEqual(module.mode([10, 10, 1]), 10)
        self.assertIsNone(module.mode([10, 1]))
        self.assertIsNone(module.mode([]))

    def test_frame_anomalies(self) -> None:
        rows = {
            "fructofuranose/run2": [
                {
                    "route_key": "fructofuranose/run2",
                    "mapped_wham_key": "fructofuranose/US-2",
                    "center_A": 12.6,
                    "frame_count": 1,
                    "path": "a",
                    "sha256": "x",
                },
                {
                    "route_key": "fructofuranose/run2",
                    "mapped_wham_key": "fructofuranose/US-2",
                    "center_A": 12.9,
                    "frame_count": 10,
                    "path": "b",
                    "sha256": "y",
                },
            ],
            "glucose/run1": [
                {
                    "route_key": "glucose/run1",
                    "mapped_wham_key": "glucose/US-1",
                    "center_A": 0.0,
                    "frame_count": 100,
                    "path": "c",
                    "sha256": "z",
                }
            ],
        }
        result = module.frame_anomalies(rows)
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0]["route_key"], "fructofuranose/run2")
        self.assertEqual(result[0]["unique_frame_counts"], [1, 10])
        self.assertIsNone(result[0]["modal_frame_count"])
        self.assertEqual(result[0]["anomalous_file_count"], 2)

    def test_center_mismatch(self) -> None:
        route_summaries = [
            {
                "route_key": "galactose/run1b",
                "mapped_wham_key": "galactose/US-1",
                "window_centres_A": [0.0, 0.3, 0.6],
            },
            {
                "route_key": "glucose/run1",
                "mapped_wham_key": "glucose/US-1",
                "window_centres_A": [0.0, 0.3, 0.6],
            },
        ]
        wham = {
            "galactose/US-1": {"centres_A": [0.0, 0.3, 0.6, 2.25]},
            "glucose/US-1": {"centres_A": [0.0, 0.3, 0.6]},
        }
        result = module.center_mismatches(route_summaries, wham)
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0]["route_key"], "galactose/run1b")
        self.assertEqual(result[0]["trajectory_only_centres_A"], [])
        self.assertEqual(result[0]["wham_only_centres_A"], [2.25])

    def test_full_synthetic_report(self) -> None:
        members = [
            member(path("fructofuranose", "run2", 12.6), "f2-126", 1),
            member(path("fructofuranose", "run2", 12.9), "f2-129", 10),
            member(path("fructofuranose", "run3", 12.6), "f3-126", 10),
        ]
        for run in (1, 2, 3):
            route = f"run{run}b"
            members.extend(
                [
                    member(path("galactose", route, 0.0), f"g{run}-same", 100),
                    member(path("galactose", route, 0.3), f"g{run}-same", 100),
                ]
            )
        inventory = {
            "archive": {
                "filename": "US_Trajectories.zip",
                "size_bytes": 1,
                "md5": "x",
                "sha256": "y",
            },
            "members": members,
        }
        wham_routes = []
        for run in (1, 2, 3):
            wham_routes.append(
                {
                    "system": "galactose",
                    "route": f"US-{run}",
                    "centres_A": [0.0, 0.3, 2.25],
                    "frames_minimum_window": 15001,
                    "frames_median_window": 15001,
                    "frames_maximum_window": 15001,
                }
            )
        wham_routes.extend(
            [
                {
                    "system": "fructofuranose",
                    "route": "US-2",
                    "centres_A": [12.6, 12.9],
                    "frames_minimum_window": 205,
                    "frames_median_window": 15001,
                    "frames_maximum_window": 15001,
                },
                {
                    "system": "fructofuranose",
                    "route": "US-3",
                    "centres_A": [12.6],
                    "frames_minimum_window": 15001,
                    "frames_median_window": 15001,
                    "frames_maximum_window": 15001,
                },
            ]
        )
        wham = {
            "archive": {
                "filename": "WHAM_PMF.zip",
                "size_bytes": 2,
                "md5": "a",
                "sha256": "b",
            },
            "routes": wham_routes,
            "anomalies": [
                {
                    "center_A": 12.6,
                    "frame_count": 205,
                    "path": "WHAM_PMF/histogram-collection-fructofuranose-US-2/12-6-test-dist_vs_t.dat",
                }
            ],
        }
        route_summaries = [
            {
                "route_key": "fructofuranose/run2",
                "mapped_wham_key": "fructofuranose/US-2",
                "window_centres_A": [12.6, 12.9],
                "trajectory_frames_per_window_values": [1, 10],
                "wham_frames_minimum_window": 205,
                "wham_frames_median_window": 15001,
                "wham_frames_maximum_window": 15001,
            },
            {
                "route_key": "fructofuranose/run3",
                "mapped_wham_key": "fructofuranose/US-3",
                "window_centres_A": [12.6],
                "trajectory_frames_per_window_values": [10],
                "wham_frames_minimum_window": 15001,
                "wham_frames_median_window": 15001,
                "wham_frames_maximum_window": 15001,
            },
        ]
        for run in (1, 2, 3):
            route_summaries.append(
                {
                    "route_key": f"galactose/run{run}b",
                    "mapped_wham_key": f"galactose/US-{run}",
                    "window_centres_A": [0.0, 0.3],
                    "trajectory_frames_per_window_values": [100],
                    "wham_frames_minimum_window": 15001,
                    "wham_frames_median_window": 15001,
                    "wham_frames_maximum_window": 15001,
                }
            )
        lineage = {
            "trajectory_file_count": len(members),
            "raw_route_directory_count": 5,
            "unique_complete_route_content_signature_count": 5,
            "deposited_wham_route_count": 5,
            "route_summaries": route_summaries,
            "duplicate_complete_route_groups": [],
            "within_route_duplicate_groups_by_route": {
                "galactose/run1b": [
                    {"centres_A": [0.0, 0.3], "sha256": "g1-same"}
                ],
                "galactose/run2b": [
                    {"centres_A": [0.0, 0.3], "sha256": "g2-same"}
                ],
            },
            "decision": {
                "all_routes_have_51_windows_0_to_15_by_0_3A": False,
                "galactose_run4a_and_run5a_are_exact_duplicate_routes": False,
                "fourteen_route_directories_reduce_to_thirteen_unique_content_signatures": False,
            },
        }

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            inventory_path = root / "inventory.json"
            wham_path = root / "wham.json"
            lineage_path = root / "lineage.json"
            output_path = root / "output.json"
            inventory_path.write_text(json.dumps(inventory), encoding="utf-8")
            wham_path.write_text(json.dumps(wham), encoding="utf-8")
            lineage_path.write_text(json.dumps(lineage), encoding="utf-8")

            report = module.run(
                inventory_path,
                wham_path,
                lineage_path,
                output_path,
            )
            decision = report["decision"]
            self.assertTrue(
                decision[
                    "frame_failure_routes_are_exactly_fructofuranose_run2_and_run3"
                ]
            )
            self.assertTrue(
                decision["fructofuranose_run2_has_one_and_ten_frame_files"]
            )
            self.assertTrue(
                decision["fructofuranose_run3_has_only_ten_frame_files"]
            )
            self.assertTrue(
                decision[
                    "one_frame_visualization_file_matches_wham_short_window_center"
                ]
            )
            self.assertFalse(
                decision[
                    "all_center_mismatches_are_galactose_pathB_extra_2_25A"
                ]
            )
            self.assertFalse(
                decision["histogram_generation_reproducible_from_us_archive"]
            )
            self.assertEqual(json.loads(output_path.read_text()), report)


if __name__ == "__main__":
    unittest.main()
