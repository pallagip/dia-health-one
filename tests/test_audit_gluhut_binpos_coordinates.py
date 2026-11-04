from __future__ import annotations

import importlib.util
import math
import struct
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "research" / "audit_gluhut_binpos_coordinates.py"
spec = importlib.util.spec_from_file_location("binpos_audit", MODULE_PATH)
assert spec and spec.loader
module = importlib.util.module_from_spec(spec)
sys.modules["binpos_audit"] = module
spec.loader.exec_module(module)


def make_binpos(frames: list[list[float]], natom: int, endian: str = "<") -> bytes:
    payload = bytearray(module.BINPOS_MAGIC)
    for coordinates in frames:
        if len(coordinates) != 3 * natom:
            raise ValueError("wrong coordinate count")
        payload.extend(struct.pack(endian + "i", natom))
        payload.extend(struct.pack(endian + f"{3*natom}f", *coordinates))
    return bytes(payload)


class BinposCoordinateTests(unittest.TestCase):
    def test_little_endian_layout_and_frames(self) -> None:
        frames = [
            [0.0, 0.0, 0.0, 1.0, 2.0, 3.0],
            [0.5, 0.0, 0.0, 2.0, 2.0, 3.0],
        ]
        payload = make_binpos(frames, natom=2, endian="<")
        layout = module.binpos_layout(payload, expected_natom=2)
        self.assertEqual(layout["endian"], "<")
        self.assertEqual(layout["frame_count"], 2)
        self.assertEqual(layout["record_bytes"], 28)
        self.assertEqual(
            module.frame_coordinates(payload, 1, layout),
            tuple(frames[1]),
        )

    def test_big_endian_detection(self) -> None:
        payload = make_binpos(
            [[0.0, 0.0, 0.0, 1.0, 1.0, 1.0]],
            natom=2,
            endian=">",
        )
        layout = module.binpos_layout(payload, expected_natom=2)
        self.assertEqual(layout["endian"], ">")
        self.assertEqual(layout["frame_count"], 1)

    def test_invalid_magic_and_size_fail(self) -> None:
        with self.assertRaises(module.BinposAuditError):
            module.binpos_layout(b"bad!" + b"\x00" * 32, expected_natom=2)
        payload = make_binpos(
            [[0.0, 0.0, 0.0, 1.0, 1.0, 1.0]],
            natom=2,
        ) + b"x"
        with self.assertRaises(module.BinposAuditError):
            module.binpos_layout(payload, expected_natom=2)

    def test_centre_of_mass(self) -> None:
        coordinates = (
            0.0,
            0.0,
            0.0,
            2.0,
            0.0,
            0.0,
            10.0,
            0.0,
            0.0,
        )
        com = module.centre_of_mass(coordinates, [1, 2], [1.0, 3.0, 1.0])
        self.assertEqual(com, (1.5, 0.0, 0.0))
        self.assertEqual(
            module.distance(com, (10.0, 0.0, 0.0)),
            8.5,
        )

    def test_correlation_and_schedule_rmse(self) -> None:
        values = [0.0, 1.0, 2.0, 3.0]
        self.assertAlmostEqual(
            module.pearson_index_correlation(values) or 0.0,
            1.0,
            places=12,
        )
        self.assertAlmostEqual(
            module.linear_schedule_rmse(values, 0.0, 3.0),
            0.0,
            places=12,
        )
        self.assertIsNone(module.pearson_index_correlation([1.0]))

    def test_trajectory_audit(self) -> None:
        # Four atoms: host group 1/2 and guest group 3/4. The guest moves away.
        frames = []
        for step in range(6):
            guest_x = float(step * 2)
            frames.append(
                [
                    0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0,
                    guest_x, 0.0, 0.0,
                    guest_x, 0.0, 0.0,
                ]
            )
        payload = make_binpos(frames, natom=4)
        topology = {
            "system_name": "glucose",
            "path": "glucose.parm7",
            "sha256": "top",
            "natom": 4,
            "masses": [12.0, 12.0, 12.0, 12.0],
        }
        restraint = {
            "path": "glucose.rst",
            "sha256": "rst",
            "igr1": [1, 2],
            "igr2": [3, 4],
            "r2": 0.0,
            "r2a": 10.0,
        }
        row = module.audit_trajectory(
            "glucose.binpos", payload, topology, restraint
        )
        self.assertEqual(row["layout"]["frame_count"], 6)
        self.assertEqual(row["reaction_coordinate"]["first_A"], 0.0)
        self.assertEqual(row["reaction_coordinate"]["last_A"], 10.0)
        self.assertTrue(
            row["reaction_coordinate"]["smd_like_forward_unbinding_path"]
        )
        self.assertAlmostEqual(
            row["reaction_coordinate"]["rmse_from_linear_0_to_15A_schedule_A"],
            0.0,
            places=12,
        )

    def test_nonfinite_coordinate_rejected(self) -> None:
        with self.assertRaises(module.BinposAuditError):
            module.validate_finite_coordinates(
                [0.0, math.inf], "x.binpos", 0
            )

    def test_safe_member_rejects_traversal(self) -> None:
        for value in ("../x", "/x", "a\\b"):
            with self.assertRaises(module.BinposAuditError):
                module.safe_member(value)


if __name__ == "__main__":
    unittest.main()
