from __future__ import annotations

import importlib.util
import struct
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "research" / "binpos_native_reader.py"
spec = importlib.util.spec_from_file_location("binpos_native", MODULE_PATH)
assert spec and spec.loader
module = importlib.util.module_from_spec(spec)
sys.modules["binpos_native"] = module
spec.loader.exec_module(module)


def write_fixture(path: Path, xyz: np.ndarray, endian: str = "<") -> None:
    xyz = np.asarray(xyz, dtype=np.float32)
    with path.open("wb") as handle:
        handle.write(b"fxyz")
        for frame in xyz:
            handle.write(struct.pack(f"{endian}i", frame.shape[0]))
            handle.write(frame.astype(f"{endian}f4").tobytes())


class NativeBinposReaderTests(unittest.TestCase):
    def test_little_endian_round_trip(self) -> None:
        expected = np.arange(2 * 3 * 3, dtype=np.float32).reshape(2, 3, 3)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "fixture.binpos"
            write_fixture(path, expected, "<")
            actual, version = module.read_binpos(path, expected_natoms=3)
        np.testing.assert_allclose(actual, expected)
        self.assertEqual(actual.dtype, np.float64)
        self.assertEqual(version, module.READER_VERSION)

    def test_big_endian_round_trip(self) -> None:
        expected = np.asarray(
            [[[1.5, 2.5, 3.5], [4.5, 5.5, 6.5]]], dtype=np.float32
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "fixture.binpos"
            write_fixture(path, expected, ">")
            actual, _ = module.read_binpos(path, expected_natoms=2)
        np.testing.assert_allclose(actual, expected)

    def test_bad_magic_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "bad.binpos"
            path.write_bytes(b"bad!" + b"\x00" * 20)
            with self.assertRaises(module.BinposFormatError):
                module.read_binpos(path, expected_natoms=1)

    def test_truncated_frame_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "bad.binpos"
            path.write_bytes(b"fxyz" + struct.pack("<i", 2) + b"\x00" * 8)
            with self.assertRaises(module.BinposFormatError):
                module.read_binpos(path, expected_natoms=2)

    def test_frame_atom_mismatch_rejected(self) -> None:
        expected = np.zeros((2, 2, 3), dtype=np.float32)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "bad.binpos"
            write_fixture(path, expected, "<")
            payload = bytearray(path.read_bytes())
            second_header = 4 + (4 + 2 * 3 * 4)
            struct.pack_into("<i", payload, second_header, 3)
            path.write_bytes(payload)
            with self.assertRaises(module.BinposFormatError):
                module.read_binpos(path, expected_natoms=2)

    def test_expected_atom_count_disambiguates_endian(self) -> None:
        expected = np.zeros((1, 3, 3), dtype=np.float32)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "fixture.binpos"
            write_fixture(path, expected, "<")
            with self.assertRaises(module.BinposFormatError):
                module.read_binpos(path, expected_natoms=4)


if __name__ == "__main__":
    unittest.main()
