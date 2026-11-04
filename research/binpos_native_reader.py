#!/usr/bin/env python3
"""Minimal, auditable reader for Scripps/AMBER BINPOS coordinates.

Format, as implemented by the historical VMD/MDTraj ``binposplugin.c``:

- four-byte ASCII magic ``fxyz``;
- for every frame, one int32 atom count;
- immediately followed by ``3 * natom`` float32 coordinates in Å;
- the next frame begins with another int32 atom count.

The reader supports little- and big-endian files, validates every frame atom
count, rejects trailing or truncated bytes, and returns float64 NumPy arrays for
analysis. It reads coordinates only; BINPOS contains no topology, time or box.
"""
from __future__ import annotations

import struct
from pathlib import Path

import numpy as np

READER_VERSION = "diahealthone-native-binpos-v1"
MAGIC = b"fxyz"


class BinposFormatError(ValueError):
    pass


def _choose_endian(payload: bytes, expected_natoms: int | None) -> tuple[str, int]:
    if len(payload) < 8:
        raise BinposFormatError("BINPOS file is too short for magic and atom count")
    candidates = []
    for endian in ("<", ">"):
        natom = struct.unpack_from(f"{endian}i", payload, 4)[0]
        plausible = 0 < natom < 100_000_000
        if expected_natoms is not None:
            plausible = plausible and natom == expected_natoms
        if plausible:
            candidates.append((endian, natom))
    if len(candidates) != 1:
        raise BinposFormatError(
            "could not determine unique BINPOS endian/atom count; "
            f"candidates={candidates}, expected_natoms={expected_natoms}"
        )
    return candidates[0]


def read_binpos(
    path: str | Path,
    expected_natoms: int | None = 183,
) -> tuple[np.ndarray, str]:
    file_path = Path(path)
    payload = file_path.read_bytes()
    if payload[:4] != MAGIC:
        raise BinposFormatError(
            f"{file_path}: expected magic {MAGIC!r}, found {payload[:4]!r}"
        )

    endian, natom = _choose_endian(payload, expected_natoms)
    frame_bytes = 4 + natom * 3 * 4
    data_bytes = len(payload) - 4
    if data_bytes <= 0 or data_bytes % frame_bytes != 0:
        raise BinposFormatError(
            f"{file_path}: {data_bytes} data bytes are not divisible by "
            f"frame size {frame_bytes}"
        )
    frame_count = data_bytes // frame_bytes
    xyz = np.empty((frame_count, natom, 3), dtype=np.float32)
    float_dtype = np.dtype(f"{endian}f4")

    offset = 4
    for frame in range(frame_count):
        frame_natom = struct.unpack_from(f"{endian}i", payload, offset)[0]
        if frame_natom != natom:
            raise BinposFormatError(
                f"{file_path}: frame {frame} atom count {frame_natom} != {natom}"
            )
        offset += 4
        count = natom * 3
        coordinates = np.frombuffer(
            payload,
            dtype=float_dtype,
            count=count,
            offset=offset,
        )
        if coordinates.size != count:
            raise BinposFormatError(
                f"{file_path}: frame {frame} is truncated"
            )
        xyz[frame] = coordinates.reshape(natom, 3)
        offset += count * 4

    if offset != len(payload):
        raise BinposFormatError(
            f"{file_path}: reader consumed {offset} of {len(payload)} bytes"
        )
    if not np.isfinite(xyz).all():
        raise BinposFormatError(f"{file_path}: non-finite coordinate encountered")

    return xyz.astype(np.float64), READER_VERSION
