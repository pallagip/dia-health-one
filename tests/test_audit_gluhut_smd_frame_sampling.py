from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "research" / "audit_gluhut_smd_frame_sampling.py"
spec = importlib.util.spec_from_file_location("frame_sampling", MODULE_PATH)
assert spec and spec.loader
module = importlib.util.module_from_spec(spec)
sys.modules["frame_sampling"] = module
spec.loader.exec_module(module)


class FrameSamplingTests(unittest.TestCase):
    def test_parse_assignments(self) -> None:
        lines = [
            "&cntrl",
            " imin=0, nstlim=15000000, dt=0.002,",
            " ntwx=1000, ntwr=1000, ntpr=1000,",
            "/",
        ]
        parsed = module.parse_assignments(lines)
        self.assertEqual(parsed["nstlim"], 15000000.0)
        self.assertEqual(parsed["dt"], 0.002)
        self.assertEqual(parsed["ntwx"], 1000.0)

    def test_classify_normal_trajectory(self) -> None:
        row = {
            "system_name": "glucose",
            "path": "glucose-pull-1-visualisation.binpos",
            "sha256": "abc",
            "deposited_reaction_coordinate": {
                "frame_count": 1000,
                "first_A": 0.1,
                "last_A": 15.05,
                "net_change_A": 14.95,
            },
        }
        result = module.classify_trajectory(row)
        self.assertTrue(
            result["consistent_with_reported_15A_change_within_tolerance"]
        )
        self.assertTrue(result["filename_contains_visualisation"])

    def test_classify_anomalous_trajectory(self) -> None:
        row = {
            "system_name": "galactose",
            "path": "galactose-pull-4-A.binpos",
            "sha256": "abc",
            "deposited_reaction_coordinate": {
                "frame_count": 1000,
                "first_A": 0.73,
                "last_A": 22.126,
                "net_change_A": 21.396,
            },
        }
        result = module.classify_trajectory(row)
        self.assertFalse(
            result["consistent_with_reported_15A_change_within_tolerance"]
        )
        self.assertFalse(result["filename_contains_visualisation"])

    def test_full_summary(self) -> None:
        inventory = {
            "members": [
                {
                    "path": "SMD_trajectories/Input-files/SMD.in",
                    "preview": {
                        "encoding": "utf-8",
                        "head": [
                            "&cntrl",
                            " nstlim=15000000, dt=0.002,",
                            " ntwx=1000, ntwr=1000, ntpr=1000,",
                        ],
                    },
                }
            ]
        }
        trajectories = []
        systems = {
            "glucose": [15.0] * 5,
            "fructofuranose": [15.0] * 5,
            "fructopyranose": [20.0] * 5,
            "galactose": [15.0, 15.0, 15.0, 22.0, 22.0],
        }
        counter = 0
        for system, net_changes in systems.items():
            for net_change in net_changes:
                counter += 1
                trajectories.append(
                    {
                        "system_name": system,
                        "path": f"{system}-{counter}-visualisation.binpos",
                        "sha256": f"hash-{counter}",
                        "deposited_reaction_coordinate": {
                            "frame_count": 1000,
                            "first_A": 0.0,
                            "last_A": net_change,
                            "net_change_A": net_change,
                        },
                    }
                )
        coordinate = {"trajectories": trajectories}
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            inventory_path = root / "inventory.json"
            coordinate_path = root / "coordinate.json"
            output_path = root / "output.json"
            inventory_path.write_text(json.dumps(inventory), encoding="utf-8")
            coordinate_path.write_text(json.dumps(coordinate), encoding="utf-8")
            report = module.summarize(
                inventory_path, coordinate_path, output_path
            )
        self.assertEqual(
            report["smd_control"][
                "expected_coordinate_outputs_from_nstlim_div_ntwx"
            ],
            15000,
        )
        self.assertEqual(
            report["deposit"]["requested_output_to_deposited_frame_ratio"],
            15.0,
        )
        self.assertEqual(
            report["trajectory_classification"][
                "reported_15A_consistent_count"
            ],
            13,
        )
        self.assertEqual(
            report["trajectory_classification"][
                "reported_15A_inconsistent_count"
            ],
            7,
        )
        self.assertFalse(
            report["decision"]
            ["deposited_files_are_sufficient_for_time_resolved_work_reconstruction"]
        )


if __name__ == "__main__":
    unittest.main()
