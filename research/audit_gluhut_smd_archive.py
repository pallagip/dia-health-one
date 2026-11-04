#!/usr/bin/env python3
"""Inventory the public GluHUT steered-MD archive for reusable chemistry inputs.

The audit verifies the deposited archive, rejects unsafe members, hashes every
file, classifies topology/coordinate/parameter/protocol artifacts, and includes
bounded text previews. It does not execute a simulation or establish that any
file is chemically sufficient for PaG or MC-Q1.
"""
from __future__ import annotations

import argparse
import collections
import hashlib
import json
import stat
import zipfile
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath

ARCHIVE_MD5 = "9f3c6733dc608d2bd37c488c560a5443"
ARCHIVE_SIZE_BYTES = 37931874
CHUNK_SIZE = 1024 * 1024
TEXT_PREVIEW_LIMIT = 256 * 1024
TEXT_EXTENSIONS = {
    ".txt", ".md", ".in", ".mdin", ".out", ".log", ".sh", ".py",
    ".pdb", ".mol2", ".frcmod", ".lib", ".off", ".prmtop", ".parm7",
    ".rst7", ".inpcrd", ".crd", ".top", ".dat", ".csv", ".json",
}
TOPOLOGY_EXTENSIONS = {".prmtop", ".parm7", ".top", ".psf"}
COORDINATE_EXTENSIONS = {".pdb", ".rst7", ".inpcrd", ".crd", ".mol2"}
PARAMETER_EXTENSIONS = {".frcmod", ".lib", ".off", ".mol2"}
PROTOCOL_EXTENSIONS = {".in", ".mdin", ".sh", ".py", ".txt", ".md"}
TRAJECTORY_EXTENSIONS = {".nc", ".mdcrd", ".dcd", ".xtc", ".trr"}


class InventoryError(RuntimeError):
    pass


def digest_file(path: Path, algorithm: str) -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(CHUNK_SIZE), b""):
            digest.update(block)
    return digest.hexdigest()


def safe_member(name: str) -> PurePosixPath:
    if "\\" in name:
        raise InventoryError(f"backslash in member path: {name!r}")
    path = PurePosixPath(name)
    if path.is_absolute() or any(part in {"", ".", ".."} for part in path.parts):
        raise InventoryError(f"unsafe member path: {name!r}")
    return path


def member_kind(info: zipfile.ZipInfo) -> str:
    if info.is_dir():
        return "directory"
    mode = (info.external_attr >> 16) & 0xFFFF
    if stat.S_ISLNK(mode):
        return "symlink"
    file_type = stat.S_IFMT(mode)
    if file_type not in {0, stat.S_IFREG}:
        return "special"
    return "file"


def classify(extension: str) -> list[str]:
    classes: list[str] = []
    if extension in TOPOLOGY_EXTENSIONS:
        classes.append("topology")
    if extension in COORDINATE_EXTENSIONS:
        classes.append("coordinate")
    if extension in PARAMETER_EXTENSIONS:
        classes.append("parameter")
    if extension in PROTOCOL_EXTENSIONS:
        classes.append("protocol_or_documentation")
    if extension in TRAJECTORY_EXTENSIONS:
        classes.append("trajectory")
    return classes or ["other"]


def hash_member(archive: zipfile.ZipFile, info: zipfile.ZipInfo) -> str:
    digest = hashlib.sha256()
    with archive.open(info, "r") as handle:
        for block in iter(lambda: handle.read(CHUNK_SIZE), b""):
            digest.update(block)
    return digest.hexdigest()


def bounded_preview(archive: zipfile.ZipFile, info: zipfile.ZipInfo) -> dict | None:
    extension = PurePosixPath(info.filename).suffix.lower()
    if extension not in TEXT_EXTENSIONS or info.file_size > TEXT_PREVIEW_LIMIT:
        return None
    payload = archive.read(info)
    try:
        text = payload.decode("utf-8")
    except UnicodeDecodeError:
        return {"encoding": "binary_or_non_utf8", "line_count": None}
    lines = text.splitlines()
    return {
        "encoding": "utf-8",
        "line_count": len(lines),
        "head": lines[:30],
        "tail": lines[-10:] if len(lines) > 30 else [],
    }


def audit(archive_path: Path, output_json: Path) -> dict:
    if archive_path.stat().st_size != ARCHIVE_SIZE_BYTES:
        raise InventoryError(
            f"size mismatch: expected {ARCHIVE_SIZE_BYTES}, got {archive_path.stat().st_size}"
        )
    md5 = digest_file(archive_path, "md5")
    if md5 != ARCHIVE_MD5:
        raise InventoryError(f"MD5 mismatch: expected {ARCHIVE_MD5}, got {md5}")
    archive_sha256 = digest_file(archive_path, "sha256")

    rows: list[dict] = []
    seen: set[str] = set()
    with zipfile.ZipFile(archive_path) as archive:
        bad = archive.testzip()
        if bad:
            raise InventoryError(f"ZIP CRC failure: {bad}")
        for info in archive.infolist():
            safe_member(info.filename)
            if info.filename in seen:
                raise InventoryError(f"duplicate member: {info.filename}")
            seen.add(info.filename)
            if info.flag_bits & 0x1:
                raise InventoryError(f"encrypted member: {info.filename}")
            kind = member_kind(info)
            if kind in {"symlink", "special"}:
                raise InventoryError(f"unsupported {kind} member: {info.filename}")
            extension = PurePosixPath(info.filename).suffix.lower()
            row = {
                "path": info.filename,
                "kind": kind,
                "extension": extension,
                "classifications": [] if kind == "directory" else classify(extension),
                "compressed_size_bytes": info.compress_size,
                "uncompressed_size_bytes": info.file_size,
                "crc32_hex": f"{info.CRC:08x}",
                "sha256": None if kind == "directory" else hash_member(archive, info),
                "preview": None if kind == "directory" else bounded_preview(archive, info),
            }
            rows.append(row)

    files = [row for row in rows if row["kind"] == "file"]
    type_counts = collections.Counter(row["extension"] or "<none>" for row in files)
    class_counts = collections.Counter(
        classification
        for row in files
        for classification in row["classifications"]
    )
    relevant = [
        row
        for row in files
        if any(
            item in row["classifications"]
            for item in ("topology", "coordinate", "parameter", "protocol_or_documentation")
        )
    ]
    exact_input_components = {
        "topology_present": any("topology" in row["classifications"] for row in files),
        "coordinate_present": any("coordinate" in row["classifications"] for row in files),
        "parameter_present": any("parameter" in row["classifications"] for row in files),
        "protocol_present": any("protocol_or_documentation" in row["classifications"] for row in files),
    }
    report = {
        "schema": "dia-health-one.gluhut-smd-archive-inventory.v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "archive_inventory_complete_not_chemistry_validated",
        "scope": "public GluHUT steered-MD archive provenance and input discovery",
        "not_target_claim": "not NNC2215, PaG, MC-Q1, affinity, efficacy, safety, or a drug result",
        "archive": {
            "filename": archive_path.name,
            "size_bytes": archive_path.stat().st_size,
            "md5": md5,
            "sha256": archive_sha256,
            "member_count": len(rows),
            "file_count": len(files),
            "directory_count": len(rows) - len(files),
            "total_uncompressed_bytes": sum(row["uncompressed_size_bytes"] for row in files),
            "crc_verified": True,
        },
        "extension_counts": dict(sorted(type_counts.items())),
        "classification_counts": dict(sorted(class_counts.items())),
        "exact_input_components": exact_input_components,
        "exact_reusable_input_bundle_present": all(exact_input_components.values()),
        "relevant_members": relevant,
        "members": rows,
        "decision": {
            "can_begin_exact_parent_parameter_audit_from_this_archive_alone": all(
                exact_input_components.values()
            ),
            "candidate_ranking_authorized": False,
            "next_action": (
                "audit identified topology, coordinates, parameters and protocol"
                if all(exact_input_components.values())
                else "inspect other deposited archives or request missing exact inputs from authors"
            ),
        },
        "limitations": [
            "File presence does not establish that topology, charges, protonation or parameters are chemically correct.",
            "SMD inputs may describe a truncated public analogue rather than the complete NNC2215 conjugate.",
            "No file in this audit is assumed transferable to PaG or MC-Q1 without atom-level review.",
        ],
    }
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser()
    result.add_argument("--archive", type=Path, required=True)
    result.add_argument(
        "--output-json",
        type=Path,
        default=Path("research/data/gluhut-smd-archive-inventory.json"),
    )
    return result


def main() -> int:
    args = parser().parse_args()
    report = audit(args.archive, args.output_json)
    print(
        json.dumps(
            {
                "archive": report["archive"],
                "classification_counts": report["classification_counts"],
                "exact_input_components": report["exact_input_components"],
                "decision": report["decision"],
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
