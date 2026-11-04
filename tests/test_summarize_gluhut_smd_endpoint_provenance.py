from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "research" / "summarize_gluhut_smd_endpoint_provenance.py"
spec = importlib.util.spec_from_file_location("endpoint_audit", MODULE_PATH)
assert spec and spec.loader
module = importlib.util.module_from_spec(spec)
sys.modules["endpoint_audit"] = module
spec.loader.exec_module(module)


def row(path: str, last: float, file_hash: str, coord_hash: str) -> dict:
    return {
        "path": path,
        "last_A": last,
        "file_sha256": file_hash,
        "coordinate_sha256": coord_hash,
    }


class EndpointProvenanceTests(unittest.TestCase):
    def test_clusters_with_bounded_diameter(self) -> None:
        rows = [
            row("a", 14.90, "f1", "c1"),
            row("b", 15.05, "f2", "c2"),
            row("c", 15.20, "f3", "c3"),
            row("d", 20.00, "f4", "c4"),
        ]
        clusters = module.cluster_rows(rows, 0.25)
        self.assertEqual(len(clusters), 3)
        self.assertEqual(clusters[0]["paths"], ["a", "b"])
        self.assertEqual(clusters[1]["paths"], ["c"])
        self.assertEqual(clusters[2]["paths"], ["d"])

    def test_closest_expected_target(self) -> None:
        self.assertEqual(
            module.closest_expected_target(14.9)["closest_target_A"],
            15.0,
        )
        self.assertTrue(
            module.closest_expected_target(20.2)["matches_within_tolerance"]
        )
        self.assertFalse(
            module.closest_expected_target(22.0)["matches_within_tolerance"]
        )

    def test_duplicate_groups(self) -> None:
        rows = [
            {"path": "a", "system_name": "x", "file_sha256": "same"},
            {"path": "b", "system_name": "x", "file_sha256": "same"},
            {"path": "c", "system_name": "y", "file_sha256": "other"},
        ]
        groups = module.duplicate_groups(rows, "file_sha256")
        self.assertEqual(len(groups), 1)
        self.assertEqual(groups[0]["trajectory_count"], 2)
        self.assertEqual(groups[0]["paths"], ["a", "b"])

    def test_full_summary_detects_heterogeneity(self) -> None:
        trajectories = []
        systems = {
            "glucose": [14.95, 15.02, 14.98, 15.01, 15.00],
            "galactose": [15.0, 15.1, 15.05, 22.1, 22.1],
            "fructopyranose": [20.0] * 5,
            "fructofuranose": [15.0] * 5,
        }
        counter = 0
        for system, endpoints in systems.items():
            for endpoint in endpoints:
                counter += 1
                trajectories.append(
                    {
                        "system_name": system,
                        "path": f"{system}-{counter}.binpos",
                        "sha256": f"file-{counter}",
                        "deposited_reaction_coordinate": {
                            "coordinate_sha256_float64": f"coord-{counter}",
                            "frame_count": 1000,
                            "first_A": 0.1,
                            "last_A": endpoint,
                            "net_change_A": endpoint - 0.1,
                            "last_minus_target_15A": endpoint - 15.0,
                            "fraction_positive_frame_steps": 0.8,
                        },
                    }
                )
        source = {"trajectories": trajectories}
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            input_path = root / "input.json"
            output_path = root / "output.json"
            input_path.write_text(json.dumps(source), encoding="utf-8")
            report = module.summarize(input_path, output_path)

        self.assertEqual(
            report["decision"]["heterogeneous_systems"], ["galactose"]
        )
        self.assertFalse(
            report["decision"]["one_endpoint_protocol_per_system"]
        )
        self.assertEqual(
            report["systems"]["galactose"]["endpoint_cluster_count"], 2
        )
        self.assertEqual(
            len(report["trajectories_outside_15_or_20A_target_tolerance"]),
            2,
        )


if __name__ == "__main__":
    unittest.main()
