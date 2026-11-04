#!/usr/bin/env python3
"""Inventory the public GluHUT umbrella-sampling trajectory archive.

The audit verifies the checksum and ZIP CRCs, rejects unsafe members, hashes all
files, classifies topology/coordinate/trajectory/protocol content, detects exact
duplicates, and summarizes sugar/route/window path structure. It is a provenance
preflight only: no histogram, PMF, affinity, or candidate result is calculated.
"""
from __future__ import annotations

import argparse
import collections
import hashlib
import json
import re
import stat
import zipfile
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
from typing import Any

ARCHIVE_MD5 = "c460e6a1ee97ac531fc284d8f07e381e"
ARCHIVE_SIZE_BYTES = 123756824
CHUNK_SIZE = 1024 * 1024
TEXT_PREVIEW_LIMIT = 256 * 1024

TOPOLOGY_EXTENSIONS = {".parm7", ".prmtop", ".top", ".psf"}
COORDINATE_EXTENSIONS = {".rst7", ".inpcrd", ".crd", ".pdb", ".mol2"}
TRAJECTORY_EXTENSIONS = {".binpos", ".nc", ".mdcrd", ".dcd", ".xtc", ".trr"}
SOURCE_PARAMETER_EXTENSIONS = {".frcmod", ".lib", ".off", ".mol2"}
PROTOCOL_EXTENSIONS = {
    ".in", ".mdin", ".txt", ".md", ".sh", ".py", ".dat", ".csv", ".json"
}
TEXT_EXTENSIONS = (
    TOPOLOGY_EXTENSIONS
    | COORDINATE_EXTENSIONS
    | SOURCE_PARAMETER_EXTENSIONS
    | PROTOCOL_EXTENSIONS
    | {".rst"}
)

SYSTEM_NAMES = (
    "fructofuranose",
    "fructopyranose",
    "galactose",
    "glucose",
)

ROUTE_PATTERNS = (
    re.compile(r"(?i)(?:US|umbrella|route|path|run)[ _-]*([A-Z]|\d+)"),
    re.compile(r"(?i)\b([ABC])\b"),
)
WINDOW_PATTERNS = (
    re.compile(r"(?i)(?:dist|distance|window|centre|center)[ _=-]*([0-9]+(?:\.[0-9]+)?)"),
    re.compile(r"(?<![A-Za-z0-9])([0-9]+\.[0-9]+)(?![A-Za-z0-9])"),
)


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


def looks_like_amber_restraint(payload: bytes) -> bool:
    try:
        text = payload.decode("utf-8")
    except UnicodeDecodeError:
        return False
    head = "\n".join(
        line.strip().lower()
        for line in text.splitlines()[:80]
        if line.strip()
    )
    return (
        "&rst" in head
        or "iat=" in head
        or "igr1=" in head
        or "igr2=" in head
    )


def detect_binary_format(payload: bytes, extension: str) -> str | None:
    if payload.startswith(b"fxyz"):
        return "scripps_binpos"
    if payload.startswith(b"CDF\x01") or payload.startswith(b"CDF\x02"):
        return "netcdf_classic"
    if payload.startswith(b"\x89HDF\r\n\x1a\n"):
        return "hdf5_or_netcdf4"
    if payload.startswith(b"PK\x03\x04"):
        return "nested_zip"
    if extension == ".dcd" and len(payload) >= 8:
        return "possible_dcd"
    return None


def classify(extension: str, payload: bytes) -> list[str]:
    classes: list[str] = []
    if extension in TOPOLOGY_EXTENSIONS:
        classes.append("topology")
    if extension in COORDINATE_EXTENSIONS:
        classes.append("coordinate")
    if extension == ".rst":
        if looks_like_amber_restraint(payload):
            classes.extend(["restraint_input", "protocol_or_documentation"])
        else:
            classes.append("ambiguous_rst_or_restart")
    if extension in TRAJECTORY_EXTENSIONS:
        classes.append("trajectory")
    if extension in SOURCE_PARAMETER_EXTENSIONS:
        classes.append("source_parameter")
    if extension in PROTOCOL_EXTENSIONS:
        classes.append("protocol_or_documentation")
    binary_format = detect_binary_format(payload, extension)
    if binary_format:
        classes.append(binary_format)
    return classes or ["other"]


def infer_system(path: str) -> str | None:
    lower = path.lower()
    for name in SYSTEM_NAMES:
        if name in lower:
            return name
    return None


def infer_route(path: str) -> str | None:
    for pattern in ROUTE_PATTERNS:
        matches = pattern.findall(path)
        if matches:
            return str(matches[-1]).upper()
    return None


def infer_window_candidates(path: str) -> list[float]:
    candidates: list[float] = []
    for pattern in WINDOW_PATTERNS:
        for match in pattern.findall(path):
            try:
                value = float(match)
            except ValueError:
                continue
            if 0.0 <= value <= 100.0 and value not in candidates:
                candidates.append(value)
    return candidates


def bounded_preview(payload: bytes, extension: str) -> dict[str, Any] | None:
    if extension not in TEXT_EXTENSIONS or len(payload) > TEXT_PREVIEW_LIMIT:
        return None
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


def top_level_directory(path: str) -> str:
    parts = PurePosixPath(path).parts
    return parts[0] if parts else ""


def audit(archive_path: Path, output_json: Path) -> dict[str, Any]:
    if archive_path.stat().st_size != ARCHIVE_SIZE_BYTES:
        raise InventoryError(
            f"size mismatch: expected {ARCHIVE_SIZE_BYTES}, got {archive_path.stat().st_size}"
        )
    md5 = digest_file(archive_path, "md5")
    if md5 != ARCHIVE_MD5:
        raise InventoryError(f"MD5 mismatch: expected {ARCHIVE_MD5}, got {md5}")
    archive_sha256 = digest_file(archive_path, "sha256")

    members: list[dict[str, Any]] = []
    seen: set[str] = set()
    content_hash_to_paths: dict[str, list[str]] = collections.defaultdict(list)

    with zipfile.ZipFile(archive_path) as archive:
        bad = archive.testzip()
        if bad:
            raise InventoryError(f"ZIP CRC failure: {bad}")
        for info in archive.infolist():
            safe_member(info.filename)
            if info.filename in seen:
                raise InventoryError(f"duplicate member path: {info.filename}")
            seen.add(info.filename)
            if info.flag_bits & 0x1:
                raise InventoryError(f"encrypted member: {info.filename}")
            kind = member_kind(info)
            if kind in {"symlink", "special"}:
                raise InventoryError(f"unsupported {kind}: {info.filename}")
            extension = PurePosixPath(info.filename).suffix.lower()
            if kind == "directory":
                payload = b""
                digest = None
                classes: list[str] = []
                preview = None
                binary_format = None
            else:
                payload = archive.read(info)
                digest = hashlib.sha256(payload).hexdigest()
                content_hash_to_paths[digest].append(info.filename)
                classes = classify(extension, payload)
                preview = bounded_preview(payload, extension)
                binary_format = detect_binary_format(payload, extension)
            members.append(
                {
                    "path": info.filename,
                    "kind": kind,
                    "extension": extension,
                    "classifications": classes,
                    "binary_format": binary_format,
                    "compressed_size_bytes": info.compress_size,
                    "uncompressed_size_bytes": info.file_size,
                    "crc32_hex": f"{info.CRC:08x}",
                    "sha256": digest,
                    "system_name": infer_system(info.filename),
                    "route_candidate": infer_route(info.filename),
                    "window_candidates_A": infer_window_candidates(info.filename),
                    "top_level_directory": top_level_directory(info.filename),
                    "preview": preview,
                }
            )

    files = [member for member in members if member["kind"] == "file"]
    extension_counts = collections.Counter(
        member["extension"] or "<none>" for member in files
    )
    classification_counts = collections.Counter(
        classification
        for member in files
        for classification in member["classifications"]
    )
    system_counts = collections.Counter(
        member["system_name"] or "unclassified" for member in files
    )
    top_level_counts = collections.Counter(
        member["top_level_directory"] for member in files
    )
    duplicate_groups = [
        {
            "sha256": digest,
            "paths": sorted(paths),
            "file_count": len(paths),
        }
        for digest, paths in sorted(content_hash_to_paths.items())
        if len(paths) > 1
    ]

    route_summary: dict[str, dict[str, Any]] = {}
    for system in SYSTEM_NAMES:
        system_files = [member for member in files if member["system_name"] == system]
        routes: dict[str, list[dict[str, Any]]] = collections.defaultdict(list)
        for member in system_files:
            routes[member["route_candidate"] or "unclassified"].append(member)
        route_summary[system] = {
            "file_count": len(system_files),
            "route_candidates": {
                route: {
                    "file_count": len(route_files),
                    "extensions": dict(
                        sorted(
                            collections.Counter(
                                item["extension"] or "<none>"
                                for item in route_files
                            ).items()
                        )
                    ),
                    "window_candidates_A": sorted(
                        {
                            candidate
                            for item in route_files
                            for candidate in item["window_candidates_A"]
                        }
                    ),
                    "paths": sorted(item["path"] for item in route_files),
                }
                for route, route_files in sorted(routes.items())
            },
        }

    relevant_classes = {
        "topology",
        "coordinate",
        "trajectory",
        "restraint_input",
        "source_parameter",
        "protocol_or_documentation",
        "scripps_binpos",
        "netcdf_classic",
        "hdf5_or_netcdf4",
    }
    relevant_members = [
        member
        for member in files
        if any(
            classification in relevant_classes
            for classification in member["classifications"]
        )
    ]

    report = {
        "schema": "dia-health-one.gluhut-us-archive-inventory.v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "umbrella_trajectory_archive_inventory_complete_not_lineage_validated",
        "scope": "public GluHUT umbrella-sampling trajectory archive",
        "not_target_claim": (
            "not histogram reproduction, PMF, binding affinity, NNC2215, PaG, "
            "MC-Q1, efficacy, safety, or a drug result"
        ),
        "archive": {
            "filename": archive_path.name,
            "size_bytes": archive_path.stat().st_size,
            "md5": md5,
            "sha256": archive_sha256,
            "member_count": len(members),
            "file_count": len(files),
            "directory_count": len(members) - len(files),
            "total_uncompressed_bytes": sum(
                member["uncompressed_size_bytes"] for member in files
            ),
            "crc_verified": True,
        },
        "extension_counts": dict(sorted(extension_counts.items())),
        "classification_counts": dict(sorted(classification_counts.items())),
        "system_file_counts": dict(sorted(system_counts.items())),
        "top_level_file_counts": dict(sorted(top_level_counts.items())),
        "duplicate_file_content_groups": duplicate_groups,
        "route_summary": route_summary,
        "relevant_members": relevant_members,
        "members": members,
        "decision": {
            "archive_integrity_verified": True,
            "unsafe_or_encrypted_members_detected": False,
            "exact_duplicate_file_content_detected": bool(duplicate_groups),
            "topologies_present": classification_counts.get("topology", 0) > 0,
            "coordinates_present": classification_counts.get("coordinate", 0) > 0,
            "trajectories_present": classification_counts.get("trajectory", 0) > 0,
            "source_parameter_files_present": (
                classification_counts.get("source_parameter", 0) > 0
            ),
            "route_and_window_lineage_validated": False,
            "histogram_generation_reproduced": False,
            "pmf_recalculation_authorized": False,
            "parent_pag_calibration_authorized": False,
            "candidate_ranking_authorized": False,
            "next_action": (
                "inspect archive layout and trajectory formats, match every route/window to "
                "WHAM metadata and histograms, identify duplicates or missing windows, and "
                "only then reconstruct histogram generation"
            ),
        },
        "limitations": [
            "Route and window candidates are inferred from path text and are not treated as authoritative simulation metadata.",
            "Exact duplicate hashes detect byte-identical files only.",
            "File presence does not establish correct topology, restraint, equilibration, sampling, or route lineage.",
            "No trajectory coordinate, histogram, work, PMF, or affinity analysis is performed."
        ],
    }
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return report


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument(
        "--output-json",
        type=Path,
        default=Path("research/data/gluhut-us-archive-inventory.json"),
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    report = audit(args.archive, args.output_json)
    print(
        json.dumps(
            {
                "archive": report["archive"],
                "extension_counts": report["extension_counts"],
                "classification_counts": report["classification_counts"],
                "system_file_counts": report["system_file_counts"],
                "duplicate_groups": report["duplicate_file_content_groups"],
                "decision": report["decision"],
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
