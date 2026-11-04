#!/usr/bin/env python3
"""Audit plausible heavy-atom alternatives to deposited fructofuranose H62.

The script maps the deposited fructofuranose SMD group onto the complete compiled
sugar graph, identifies the heavy atom bonded to H62, enumerates heavy atoms
adjacent to the selected furanose ring, and evaluates bounded replacement or
omission hypotheses. It does not choose a corrected restraint; author
clarification remains required.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import stat
import zipfile
from pathlib import Path, PurePosixPath
from typing import Any

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import audit_gluhut_smd_restraints as mapping_audit  # noqa: E402
import summarize_gluhut_smd_restraint_anomaly as consistency  # noqa: E402

ARCHIVE_MD5 = "9f3c6733dc608d2bd37c488c560a5443"
ARCHIVE_SIZE_BYTES = 37931874
CHUNK_SIZE = 1024 * 1024


class CandidateAuditError(RuntimeError):
    pass


def digest_file(path: Path, algorithm: str) -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(CHUNK_SIZE), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def safe_member(name: str) -> PurePosixPath:
    if "\\" in name:
        raise CandidateAuditError(f"backslash in member path: {name!r}")
    path = PurePosixPath(name)
    if path.is_absolute() or any(part in {"", ".", ".."} for part in path.parts):
        raise CandidateAuditError(f"unsafe member path: {name!r}")
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


def adjacency(topology: dict[str, Any]) -> dict[int, set[int]]:
    result = {index: set() for index in topology["atoms"]}
    for left, right in topology["bonds"]:
        result[left].add(right)
        result[right].add(left)
    return result


def selected_atom_table(indices: list[int], topology: dict[str, Any]) -> list[dict[str, Any]]:
    return [topology["atoms"][index] for index in indices]


def replacement_summary(
    label: str,
    indices: list[int],
    topology: dict[str, Any],
) -> dict[str, Any]:
    group = mapping_audit.induced_group(indices, topology)
    summary = consistency.summarize_group("fructofuranose", group)
    return {
        "label": label,
        "atom_indices_1_based": indices,
        "atom_names": summary["selected_atom_names"],
        "elements": summary["selected_elements"],
        "selected_atom_count": summary["selected_atom_count"],
        "connected_component_count": summary["connected_component_count"],
        "cyclic_core_atom_count": summary["cyclic_core_atom_count"],
        "cyclic_core_atom_names": summary["cyclic_core_atom_names"],
        "noncore_selected_atoms": summary["noncore_selected_atoms"],
        "selected_hydrogen_count": summary["selected_hydrogen_count"],
        "selection_class": summary["selection_class"],
        "all_selected_atoms_are_heavy": summary["selected_hydrogen_count"] == 0,
        "single_connected_component": summary["connected_component_count"] == 1,
        "retains_five_member_cyclic_core": summary["cyclic_core_atom_count"] == 5,
    }


def audit(archive_path: Path, output_json: Path) -> dict[str, Any]:
    if archive_path.stat().st_size != ARCHIVE_SIZE_BYTES:
        raise CandidateAuditError(
            f"size mismatch: expected {ARCHIVE_SIZE_BYTES}, got {archive_path.stat().st_size}"
        )
    md5 = digest_file(archive_path, "md5")
    if md5 != ARCHIVE_MD5:
        raise CandidateAuditError(f"MD5 mismatch: expected {ARCHIVE_MD5}, got {md5}")
    archive_sha256 = digest_file(archive_path, "sha256")

    topologies: dict[str, dict[str, Any]] = {}
    restraints: list[dict[str, Any]] = []
    seen: set[str] = set()
    with zipfile.ZipFile(archive_path) as archive:
        bad = archive.testzip()
        if bad:
            raise CandidateAuditError(f"ZIP CRC failure: {bad}")
        for info in archive.infolist():
            safe_member(info.filename)
            if info.filename in seen:
                raise CandidateAuditError(f"duplicate member: {info.filename}")
            seen.add(info.filename)
            if info.flag_bits & 0x1:
                raise CandidateAuditError(f"encrypted member: {info.filename}")
            kind = member_kind(info)
            if kind in {"symlink", "special"}:
                raise CandidateAuditError(f"unsupported {kind}: {info.filename}")
            if kind != "file":
                continue
            suffix = PurePosixPath(info.filename).suffix.lower()
            payload = archive.read(info)
            if suffix in {".parm7", ".prmtop"}:
                topology = mapping_audit.topology_atoms(info.filename, payload)
                topologies[topology["system_name"]] = topology
            elif suffix == ".rst":
                restraints.append(mapping_audit.parse_restraint(info.filename, payload))

    fructofuranose = topologies.get("fructofuranose")
    fructopyranose = topologies.get("fructopyranose")
    if fructofuranose is None or fructopyranose is None:
        raise CandidateAuditError("missing fructose topologies")
    furanose_restraints = [
        row for row in restraints if row["applies_to_systems"] == ["fructofuranose"]
    ]
    pyranose_restraints = [
        row for row in restraints if row["applies_to_systems"] == ["fructopyranose"]
    ]
    if len(furanose_restraints) != 1 or len(pyranose_restraints) != 1:
        raise CandidateAuditError("unexpected fructose restraint count")
    furanose_restraint = furanose_restraints[0]
    pyranose_restraint = pyranose_restraints[0]

    selected = furanose_restraint["igr2"]
    group = mapping_audit.induced_group(selected, fructofuranose)
    group_summary = consistency.summarize_group("fructofuranose", group)
    hydrogens = group_summary["selected_hydrogen_atoms"]
    if len(hydrogens) != 1 or hydrogens[0]["atom_name"] != "H62":
        raise CandidateAuditError(f"expected one selected H62, found {hydrogens}")
    h62_index = int(hydrogens[0]["atom_1_based"])

    graph = adjacency(fructofuranose)
    h62_neighbours = sorted(graph[h62_index])
    heavy_neighbours = [
        index
        for index in h62_neighbours
        if fructofuranose["atoms"][index]["element"] != "H"
    ]
    if len(heavy_neighbours) != 1:
        raise CandidateAuditError(
            f"expected one heavy atom bonded to H62, found {heavy_neighbours}"
        )
    h62_parent = heavy_neighbours[0]

    ring_core = set(group_summary["cyclic_core_atom_indices_1_based"])
    selected_set = set(selected)
    adjacent_guest_heavy_atoms = sorted(
        {
            neighbour
            for ring_atom in ring_core
            for neighbour in graph[ring_atom]
            if neighbour not in selected_set
            and fructofuranose["atoms"][neighbour]["residue_label"] != "MOL"
            and fructofuranose["atoms"][neighbour]["element"] != "H"
        }
    )

    hypotheses: list[dict[str, Any]] = []
    omit_h62 = [index for index in selected if index != h62_index]
    hypotheses.append(
        replacement_summary("omit_H62_use_five_ring_atoms", omit_h62, fructofuranose)
    )
    candidates = sorted(set(adjacent_guest_heavy_atoms + [h62_parent]))
    for candidate in candidates:
        replacement = [
            candidate if index == h62_index else index for index in selected
        ]
        if len(replacement) != len(set(replacement)):
            continue
        candidate_atom = fructofuranose["atoms"][candidate]
        hypotheses.append(
            replacement_summary(
                f"replace_H62_with_{candidate_atom['atom_name']}_atom_{candidate}",
                replacement,
                fructofuranose,
            )
        )

    same_numeric_indices = selected == pyranose_restraint["igr2"]
    numeric_index_crosswalk = [
        {
            "atom_index_1_based": index,
            "fructopyranose": fructopyranose["atoms"][index],
            "fructofuranose": fructofuranose["atoms"][index],
            "same_atom_name": (
                fructopyranose["atoms"][index]["atom_name"]
                == fructofuranose["atoms"][index]["atom_name"]
            ),
            "same_element": (
                fructopyranose["atoms"][index]["element"]
                == fructofuranose["atoms"][index]["element"]
            ),
        }
        for index in selected
    ]

    viable_all_heavy = [
        row
        for row in hypotheses
        if row["all_selected_atoms_are_heavy"]
        and row["single_connected_component"]
        and row["retains_five_member_cyclic_core"]
    ]

    report = {
        "schema": "dia-health-one.gluhut-fructofuranose-restraint-candidates.v1",
        "status": "bounded_candidate_replacements_identified_no_correction_selected",
        "scope": "deposited fructofuranose SMD igr2 atom-selection anomaly",
        "not_target_claim": (
            "not NNC2215, PaG, MC-Q1, binding affinity, efficacy, safety, or a drug result"
        ),
        "archive": {
            "filename": archive_path.name,
            "size_bytes": archive_path.stat().st_size,
            "md5": md5,
            "sha256": archive_sha256,
            "crc_verified": True,
        },
        "deposited_fructofuranose_selection": {
            "restraint_path": furanose_restraint["path"],
            "atom_indices_1_based": selected,
            "atoms": selected_atom_table(selected, fructofuranose),
            "cyclic_core_atom_indices_1_based": sorted(ring_core),
            "cyclic_core_atom_names": group_summary["cyclic_core_atom_names"],
            "selected_hydrogen": hydrogens[0],
            "selected_hydrogen_neighbours_1_based": h62_neighbours,
            "selected_hydrogen_heavy_parent": fructofuranose["atoms"][h62_parent],
            "connected_component_count": group_summary["connected_component_count"],
        },
        "adjacent_unselected_guest_heavy_atoms": [
            fructofuranose["atoms"][index]
            for index in adjacent_guest_heavy_atoms
        ],
        "same_numeric_index_comparison": {
            "fructopyranose_and_fructofuranose_igr2_indices_identical": (
                same_numeric_indices
            ),
            "indices": selected,
            "crosswalk": numeric_index_crosswalk,
            "all_indices_preserve_atom_name": all(
                row["same_atom_name"] for row in numeric_index_crosswalk
            ),
            "all_indices_preserve_element": all(
                row["same_element"] for row in numeric_index_crosswalk
            ),
        },
        "bounded_hypotheses": hypotheses,
        "viable_all_heavy_connected_hypotheses": viable_all_heavy,
        "decision": {
            "deposited_H62_has_one_heavy_parent": True,
            "deposited_H62_parent_atom_name": fructofuranose["atoms"][h62_parent][
                "atom_name"
            ],
            "same_numeric_indices_reused_between_fructose_templates": (
                same_numeric_indices
            ),
            "same_numeric_indices_preserve_chemical_identity": all(
                row["same_atom_name"] and row["same_element"]
                for row in numeric_index_crosswalk
            ),
            "at_least_one_all_heavy_connected_alternative_exists": bool(
                viable_all_heavy
            ),
            "automatic_correction_selected": False,
            "author_clarification_required": True,
            "new_smd_or_umbrella_generation_authorized": False,
            "candidate_ranking_authorized": False,
            "next_action": (
                "ask the authors whether H62 was intentional and which of the bounded "
                "five-ring or ring-plus-exocyclic-heavy-atom definitions was intended; "
                "do not emit a corrected production input without that answer"
            ),
        },
        "limitations": [
            "Topology adjacency can identify chemically coherent alternatives but cannot recover author intent.",
            "No coordinates, centre-of-mass values, pulling work or PMF are calculated.",
            "A viable graph replacement is not evidence that the replacement reproduces the published path.",
            "No corrected restraint file is generated by this audit.",
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
        default=Path(
            "research/data/gluhut-fructofuranose-restraint-candidates.json"
        ),
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    report = audit(args.archive, args.output_json)
    print(
        json.dumps(
            {
                "deposited": report["deposited_fructofuranose_selection"],
                "numeric_index_comparison": report[
                    "same_numeric_index_comparison"
                ],
                "viable_hypotheses": report[
                    "viable_all_heavy_connected_hypotheses"
                ],
                "decision": report["decision"],
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
