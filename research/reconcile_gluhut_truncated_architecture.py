#!/usr/bin/env python3
"""Reconcile the compiled GluHUT host with the published truncated architecture.

The public 2026 simulation paper states that GluHUT was truncated at the outer
amide linker, replacing each polycarboxylate dendrimer with a methyl group. The
2019 receptor paper shows two triethyl aromatic caps connected to three
benzamide pillars by six urea bridges.

This audit checks whether the deposited AMBER host encodes that exact graph
architecture and molecular formula. It does not reproduce RESP derivation,
validate force-field parameters, calculate affinity, or authorize transfer to
PaG, MC-Q1, or complete NNC2215.
"""
from __future__ import annotations

import argparse
import collections
import hashlib
import json
import stat
import sys
import zipfile
from pathlib import Path, PurePosixPath
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent))
import audit_gluhut_smd_restraints as amber_map  # noqa: E402

ARCHIVE_MD5 = "9f3c6733dc608d2bd37c488c560a5443"
ARCHIVE_SIZE_BYTES = 37931874
CHUNK_SIZE = 1024 * 1024

EXPECTED_HOST_FORMULA = {"C": 60, "H": 75, "N": 15, "O": 9}
EXPECTED_MODULES = {
    "two_triethyl_tris_methylene_caps_with_three_NH_each": {
        "multiplicity": 2,
        "formula": {"C": 15, "H": 24, "N": 3},
    },
    "three_N_methyl_diaminobenzamide_pillars": {
        "multiplicity": 3,
        "formula": {"C": 8, "H": 9, "N": 3, "O": 1},
    },
    "six_urea_carbonyls": {
        "multiplicity": 6,
        "formula": {"C": 1, "O": 1},
    },
}

ELEMENT_SYMBOLS = {
    1: "H",
    5: "B",
    6: "C",
    7: "N",
    8: "O",
    9: "F",
    15: "P",
    16: "S",
    17: "Cl",
    35: "Br",
    53: "I",
}


class ReconciliationError(RuntimeError):
    pass


def digest_file(path: Path, algorithm: str) -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(CHUNK_SIZE), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_hash(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
        ).encode("utf-8")
    ).hexdigest()


def safe_member(name: str) -> PurePosixPath:
    if "\\" in name:
        raise ReconciliationError(f"backslash in member path: {name!r}")
    path = PurePosixPath(name)
    if path.is_absolute() or any(part in {"", ".", ".."} for part in path.parts):
        raise ReconciliationError(f"unsafe member path: {name!r}")
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


def add_formula(target: dict[str, int], source: dict[str, int], multiplier: int = 1) -> None:
    for element, count in source.items():
        target[element] = target.get(element, 0) + count * multiplier


def module_formula_total(modules: dict[str, dict[str, Any]]) -> dict[str, int]:
    total: dict[str, int] = {}
    for module in modules.values():
        add_formula(total, module["formula"], int(module["multiplicity"]))
    return dict(sorted(total.items()))


def formula_from_atoms(atoms: list[dict[str, Any]]) -> dict[str, int]:
    counts: collections.Counter[str] = collections.Counter()
    for atom in atoms:
        symbol = atom["element"]
        if symbol not in ELEMENT_SYMBOLS.values():
            raise ReconciliationError(f"unsupported element symbol: {symbol}")
        counts[symbol] += 1
    return dict(sorted(counts.items()))


def formula_string(counts: dict[str, int]) -> str:
    order: list[str] = []
    if counts.get("C"):
        order.append("C")
        if counts.get("H"):
            order.append("H")
    order.extend(sorted(element for element in counts if element not in order))
    return "".join(
        element + (str(counts[element]) if counts[element] != 1 else "")
        for element in order
    )


def adjacency(topology: dict[str, Any]) -> dict[int, set[int]]:
    result = {index: set() for index in topology["atoms"]}
    for left, right in topology["bonds"]:
        result[left].add(right)
        result[right].add(left)
    return result


def components(indices: set[int], graph: dict[int, set[int]]) -> list[list[int]]:
    remaining = set(indices)
    result: list[list[int]] = []
    while remaining:
        start = min(remaining)
        remaining.remove(start)
        stack = [start]
        component: list[int] = []
        while stack:
            atom = stack.pop()
            component.append(atom)
            for neighbour in sorted(graph[atom], reverse=True):
                if neighbour in remaining and neighbour in indices:
                    remaining.remove(neighbour)
                    stack.append(neighbour)
        result.append(sorted(component))
    return sorted(result, key=lambda item: (min(item), len(item)))


def induced_edges(indices: set[int], topology: dict[str, Any]) -> list[list[int]]:
    return [
        [left, right]
        for left, right in topology["bonds"]
        if left in indices and right in indices
    ]


def hydrogen_neighbours(atom: int, graph: dict[int, set[int]], topology: dict[str, Any]) -> list[int]:
    return sorted(
        neighbour
        for neighbour in graph[atom]
        if topology["atoms"][neighbour]["element"] == "H"
    )


def heavy_neighbours(atom: int, graph: dict[int, set[int]], topology: dict[str, Any]) -> list[int]:
    return sorted(
        neighbour
        for neighbour in graph[atom]
        if topology["atoms"][neighbour]["element"] != "H"
    )


def aromatic_ring_components(
    host_atoms: set[int], topology: dict[str, Any], graph: dict[int, set[int]]
) -> list[list[int]]:
    aromatic = {
        atom
        for atom in host_atoms
        if topology["atoms"][atom]["element"] == "C"
        and topology["atoms"][atom]["amber_atom_type"].lower() == "ca"
    }
    rings = components(aromatic, graph)
    for ring in rings:
        ring_set = set(ring)
        edges = induced_edges(ring_set, topology)
        degrees = {
            atom: sum(neighbour in ring_set for neighbour in graph[atom])
            for atom in ring
        }
        if len(ring) != 6 or len(edges) != 6 or set(degrees.values()) != {2}:
            raise ReconciliationError(
                f"aromatic component is not a six-member ring: {ring}, edges={edges}, degrees={degrees}"
            )
    return rings


def classify_ring(
    ring_id: str,
    ring: list[int],
    host_atoms: set[int],
    topology: dict[str, Any],
    graph: dict[int, set[int]],
) -> dict[str, Any]:
    ring_set = set(ring)
    external_heavy = sorted(
        {
            neighbour
            for atom in ring
            for neighbour in graph[atom]
            if neighbour in host_atoms
            and neighbour not in ring_set
            and topology["atoms"][neighbour]["element"] != "H"
        }
    )
    ring_hydrogens = sorted(
        {
            neighbour
            for atom in ring
            for neighbour in graph[atom]
            if neighbour in host_atoms
            and topology["atoms"][neighbour]["element"] == "H"
        }
    )
    external_elements = collections.Counter(
        topology["atoms"][atom]["element"] for atom in external_heavy
    )
    if len(external_heavy) == 6 and external_elements == {"C": 6} and not ring_hydrogens:
        classification = "triethyl_tris_methylene_cap"
    elif (
        len(external_heavy) == 3
        and external_elements == {"C": 1, "N": 2}
        and len(ring_hydrogens) == 3
    ):
        classification = "diaminobenzamide_pillar"
    else:
        classification = "unclassified_aromatic_ring"
    return {
        "ring_id": ring_id,
        "classification": classification,
        "atom_indices_1_based": ring,
        "atom_names": [topology["atoms"][atom]["atom_name"] for atom in ring],
        "external_heavy_atom_indices_1_based": external_heavy,
        "external_heavy_elements": dict(sorted(external_elements.items())),
        "ring_hydrogen_atom_indices_1_based": ring_hydrogens,
        "ring_hydrogen_count": len(ring_hydrogens),
    }


def carbonyl_groups(
    host_atoms: set[int], topology: dict[str, Any], graph: dict[int, set[int]]
) -> list[dict[str, Any]]:
    groups: list[dict[str, Any]] = []
    for atom in sorted(host_atoms):
        current = topology["atoms"][atom]
        if current["element"] != "C":
            continue
        heavy = [
            neighbour
            for neighbour in graph[atom]
            if neighbour in host_atoms
            and topology["atoms"][neighbour]["element"] != "H"
        ]
        elements = sorted(topology["atoms"][neighbour]["element"] for neighbour in heavy)
        if "O" not in elements:
            continue
        if elements == ["N", "N", "O"]:
            classification = "urea_carbonyl"
        elif elements == ["C", "N", "O"]:
            classification = "amide_carbonyl"
        else:
            classification = "other_oxygenated_carbon"
        groups.append(
            {
                "atom_1_based": atom,
                "atom_name": current["atom_name"],
                "classification": classification,
                "heavy_neighbour_indices_1_based": sorted(heavy),
                "heavy_neighbour_elements": elements,
            }
        )
    return groups


def ring_membership(rings: list[dict[str, Any]]) -> dict[int, str]:
    result: dict[int, str] = {}
    for ring in rings:
        for atom in ring["atom_indices_1_based"]:
            if atom in result:
                raise ReconciliationError(f"atom {atom} belongs to multiple aromatic rings")
            result[atom] = ring["ring_id"]
    return result


def trace_urea_connections(
    urea_groups: list[dict[str, Any]],
    rings: list[dict[str, Any]],
    topology: dict[str, Any],
    graph: dict[int, set[int]],
) -> list[dict[str, Any]]:
    membership = ring_membership(rings)
    ring_class = {ring["ring_id"]: ring["classification"] for ring in rings}
    connections: list[dict[str, Any]] = []
    for group in urea_groups:
        carbonyl = group["atom_1_based"]
        nitrogen_atoms = [
            atom
            for atom in group["heavy_neighbour_indices_1_based"]
            if topology["atoms"][atom]["element"] == "N"
        ]
        endpoints: list[dict[str, Any]] = []
        for nitrogen in nitrogen_atoms:
            other_heavy = [
                neighbour
                for neighbour in heavy_neighbours(nitrogen, graph, topology)
                if neighbour != carbonyl
            ]
            if len(other_heavy) != 1:
                raise ReconciliationError(
                    f"urea N {nitrogen} has unexpected heavy neighbours {other_heavy}"
                )
            attached = other_heavy[0]
            if attached in membership:
                ring_id = membership[attached]
                endpoint = {
                    "nitrogen_atom_1_based": nitrogen,
                    "attachment_atom_1_based": attached,
                    "attachment_mode": "direct_to_aromatic_ring",
                    "ring_id": ring_id,
                    "ring_classification": ring_class[ring_id],
                }
            else:
                next_rings = sorted(
                    {
                        membership[neighbour]
                        for neighbour in graph[attached]
                        if neighbour in membership
                    }
                )
                if len(next_rings) != 1:
                    raise ReconciliationError(
                        f"urea N {nitrogen} attachment {attached} does not trace to one ring: {next_rings}"
                    )
                ring_id = next_rings[0]
                endpoint = {
                    "nitrogen_atom_1_based": nitrogen,
                    "attachment_atom_1_based": attached,
                    "attachment_mode": "via_methylene_to_aromatic_ring",
                    "ring_id": ring_id,
                    "ring_classification": ring_class[ring_id],
                }
            endpoints.append(endpoint)
        roof = [
            endpoint
            for endpoint in endpoints
            if endpoint["ring_classification"] == "triethyl_tris_methylene_cap"
        ]
        pillar = [
            endpoint
            for endpoint in endpoints
            if endpoint["ring_classification"] == "diaminobenzamide_pillar"
        ]
        if len(roof) != 1 or len(pillar) != 1:
            raise ReconciliationError(
                f"urea {carbonyl} does not connect one cap and one pillar: {endpoints}"
            )
        connections.append(
            {
                "urea_carbonyl_atom_1_based": carbonyl,
                "roof_ring_id": roof[0]["ring_id"],
                "pillar_ring_id": pillar[0]["ring_id"],
                "endpoints": endpoints,
            }
        )
    return sorted(
        connections,
        key=lambda row: (row["roof_ring_id"], row["pillar_ring_id"]),
    )


def audit_cap_branches(
    cap: dict[str, Any], topology: dict[str, Any], graph: dict[int, set[int]]
) -> dict[str, Any]:
    ring = set(cap["atom_indices_1_based"])
    branches: list[dict[str, Any]] = []
    for carbon in cap["external_heavy_atom_indices_1_based"]:
        other_heavy = [
            neighbour
            for neighbour in heavy_neighbours(carbon, graph, topology)
            if neighbour not in ring
        ]
        hydrogen_count = len(hydrogen_neighbours(carbon, graph, topology))
        if len(other_heavy) != 1:
            branch_type = "unexpected"
        else:
            neighbour = other_heavy[0]
            neighbour_element = topology["atoms"][neighbour]["element"]
            if neighbour_element == "N":
                branch_type = "methylene_to_urea_N"
            elif (
                neighbour_element == "C"
                and len(hydrogen_neighbours(neighbour, graph, topology)) == 3
            ):
                branch_type = "ethyl_CH2_to_CH3"
            else:
                branch_type = "unexpected"
        branches.append(
            {
                "branch_carbon_atom_1_based": carbon,
                "branch_carbon_name": topology["atoms"][carbon]["atom_name"],
                "branch_carbon_hydrogen_count": hydrogen_count,
                "other_heavy_neighbours_1_based": other_heavy,
                "classification": branch_type,
            }
        )
    counts = collections.Counter(branch["classification"] for branch in branches)
    return {
        "ring_id": cap["ring_id"],
        "branches": branches,
        "branch_classification_counts": dict(sorted(counts.items())),
        "matches_three_ethyl_three_methylene_pattern": counts
        == {"ethyl_CH2_to_CH3": 3, "methylene_to_urea_N": 3},
    }


def audit_pillar_amide(
    pillar: dict[str, Any],
    amide_groups: list[dict[str, Any]],
    topology: dict[str, Any],
    graph: dict[int, set[int]],
) -> dict[str, Any]:
    ring = set(pillar["atom_indices_1_based"])
    attached_amides = [
        group
        for group in amide_groups
        if any(atom in ring for atom in group["heavy_neighbour_indices_1_based"])
    ]
    if len(attached_amides) != 1:
        return {
            "ring_id": pillar["ring_id"],
            "attached_amide_count": len(attached_amides),
            "matches_N_methyl_benzamide": False,
            "attached_amides": attached_amides,
        }
    group = attached_amides[0]
    carbonyl = group["atom_1_based"]
    nitrogens = [
        atom
        for atom in group["heavy_neighbour_indices_1_based"]
        if topology["atoms"][atom]["element"] == "N"
    ]
    if len(nitrogens) != 1:
        raise ReconciliationError(f"amide {carbonyl} has unexpected N list {nitrogens}")
    nitrogen = nitrogens[0]
    other_heavy = [
        neighbour
        for neighbour in heavy_neighbours(nitrogen, graph, topology)
        if neighbour != carbonyl
    ]
    methyl = None
    if len(other_heavy) == 1:
        candidate = other_heavy[0]
        if (
            topology["atoms"][candidate]["element"] == "C"
            and len(hydrogen_neighbours(candidate, graph, topology)) == 3
        ):
            methyl = candidate
    return {
        "ring_id": pillar["ring_id"],
        "attached_amide_count": 1,
        "amide_carbonyl_atom_1_based": carbonyl,
        "amide_nitrogen_atom_1_based": nitrogen,
        "amide_N_methyl_carbon_atom_1_based": methyl,
        "matches_N_methyl_benzamide": methyl is not None,
    }


def analyse_topology(path: str, payload: bytes) -> dict[str, Any]:
    topology = amber_map.topology_atoms(path, payload)
    graph = adjacency(topology)
    host_atoms = {
        atom
        for atom, record in topology["atoms"].items()
        if record["residue_label"] == "MOL"
    }
    if not host_atoms:
        raise ReconciliationError("no MOL host residue found")
    host_formula = formula_from_atoms([topology["atoms"][atom] for atom in sorted(host_atoms)])

    raw_ring_components = aromatic_ring_components(host_atoms, topology, graph)
    rings = [
        classify_ring(f"ring_{index}", ring, host_atoms, topology, graph)
        for index, ring in enumerate(raw_ring_components, start=1)
    ]
    caps = [ring for ring in rings if ring["classification"] == "triethyl_tris_methylene_cap"]
    pillars = [ring for ring in rings if ring["classification"] == "diaminobenzamide_pillar"]

    carbonyls = carbonyl_groups(host_atoms, topology, graph)
    ureas = [group for group in carbonyls if group["classification"] == "urea_carbonyl"]
    amides = [group for group in carbonyls if group["classification"] == "amide_carbonyl"]
    other_carbonyls = [
        group for group in carbonyls if group["classification"] == "other_oxygenated_carbon"
    ]

    connections = trace_urea_connections(ureas, rings, topology, graph)
    observed_pairs = sorted(
        {(row["roof_ring_id"], row["pillar_ring_id"]) for row in connections}
    )
    expected_pairs = sorted(
        (cap["ring_id"], pillar["ring_id"])
        for cap in caps
        for pillar in pillars
    )

    cap_audits = [audit_cap_branches(cap, topology, graph) for cap in caps]
    pillar_audits = [
        audit_pillar_amide(pillar, amides, topology, graph) for pillar in pillars
    ]

    architecture_signature = {
        "host_formula": host_formula,
        "ring_classifications": [ring["classification"] for ring in rings],
        "cap_branch_counts": [audit["branch_classification_counts"] for audit in cap_audits],
        "urea_connections": observed_pairs,
        "pillar_amides": [audit["matches_N_methyl_benzamide"] for audit in pillar_audits],
        "carbonyl_counts": {
            "urea": len(ureas),
            "amide": len(amides),
            "other": len(other_carbonyls),
        },
    }

    return {
        "path": path,
        "sha256": hashlib.sha256(payload).hexdigest(),
        "system_name": topology["system_name"],
        "host_formula": host_formula,
        "host_formula_string": formula_string(host_formula),
        "rings": rings,
        "cap_audits": cap_audits,
        "pillar_audits": pillar_audits,
        "carbonyl_groups": carbonyls,
        "urea_connections": connections,
        "expected_K2_3_pairs": expected_pairs,
        "observed_urea_pairs": observed_pairs,
        "architecture_signature_sha256": canonical_hash(architecture_signature),
        "checks": {
            "host_formula_matches": host_formula == EXPECTED_HOST_FORMULA,
            "five_aromatic_six_member_rings": len(rings) == 5,
            "two_caps": len(caps) == 2,
            "three_pillars": len(pillars) == 3,
            "six_urea_carbonyls": len(ureas) == 6,
            "three_amide_carbonyls": len(amides) == 3,
            "no_other_oxygenated_carbons": not other_carbonyls,
            "both_caps_have_three_ethyl_three_methylene_pattern": all(
                audit["matches_three_ethyl_three_methylene_pattern"]
                for audit in cap_audits
            ),
            "all_pillars_have_one_N_methyl_benzamide": all(
                audit["matches_N_methyl_benzamide"] for audit in pillar_audits
            ),
            "ureas_form_complete_K2_3_cap_pillar_graph": observed_pairs
            == expected_pairs
            and len(connections) == 6,
        },
    }


def audit(archive_path: Path, output_json: Path) -> dict[str, Any]:
    if archive_path.stat().st_size != ARCHIVE_SIZE_BYTES:
        raise ReconciliationError(
            f"size mismatch: expected {ARCHIVE_SIZE_BYTES}, got {archive_path.stat().st_size}"
        )
    md5 = digest_file(archive_path, "md5")
    if md5 != ARCHIVE_MD5:
        raise ReconciliationError(f"MD5 mismatch: expected {ARCHIVE_MD5}, got {md5}")
    archive_sha256 = digest_file(archive_path, "sha256")

    systems: list[dict[str, Any]] = []
    seen: set[str] = set()
    with zipfile.ZipFile(archive_path) as archive:
        bad = archive.testzip()
        if bad:
            raise ReconciliationError(f"ZIP CRC failure: {bad}")
        for info in archive.infolist():
            safe_member(info.filename)
            if info.filename in seen:
                raise ReconciliationError(f"duplicate member: {info.filename}")
            seen.add(info.filename)
            if info.flag_bits & 0x1:
                raise ReconciliationError(f"encrypted member: {info.filename}")
            kind = member_kind(info)
            if kind in {"symlink", "special"}:
                raise ReconciliationError(f"unsupported {kind}: {info.filename}")
            if kind != "file" or PurePosixPath(info.filename).suffix.lower() not in {
                ".parm7",
                ".prmtop",
            }:
                continue
            systems.append(analyse_topology(info.filename, archive.read(info)))

    systems.sort(key=lambda row: row["system_name"])
    if len(systems) != 4:
        raise ReconciliationError(f"expected four topologies, found {len(systems)}")

    module_formula = module_formula_total(EXPECTED_MODULES)
    signatures = {system["architecture_signature_sha256"] for system in systems}
    all_checks = {
        key: all(system["checks"][key] for system in systems)
        for key in systems[0]["checks"]
    }
    report = {
        "schema": "dia-health-one.gluhut-truncated-architecture-reconciliation.v1",
        "status": "compiled_graph_matches_published_truncated_architecture",
        "scope": "public truncated GluHUT host graph in four deposited AMBER systems",
        "not_target_claim": (
            "not complete receptor 2, not NNC2215, not PaG or MC-Q1 affinity, "
            "not efficacy, safety, or a drug result"
        ),
        "primary_sources": [
            {
                "doi": "10.1038/s41557-018-0155-z",
                "evidence": (
                    "2019 receptor design and synthesis figures show two triethyl aromatic "
                    "caps, three benzamide-containing pillars, six urea bridges, and three "
                    "outer solubilising amide substituents."
                ),
            },
            {
                "doi": "10.1039/D6SC00738D",
                "evidence": (
                    "2026 methods state that the polycarboxylate dendrimer was removed at "
                    "the outer amide linker and replaced with methyl for the simulation model."
                ),
            },
        ],
        "archive": {
            "filename": archive_path.name,
            "size_bytes": archive_path.stat().st_size,
            "md5": md5,
            "sha256": archive_sha256,
            "crc_verified": True,
        },
        "published_module_formula_bookkeeping": {
            "modules": EXPECTED_MODULES,
            "calculated_total": module_formula,
            "calculated_total_string": formula_string(module_formula),
            "matches_expected_host_formula": module_formula == EXPECTED_HOST_FORMULA,
        },
        "systems": systems,
        "cross_system": {
            "architecture_signature_identical_across_four_sugars": len(signatures) == 1,
            "architecture_signature_sha256": systems[0]["architecture_signature_sha256"],
            "all_checks": all_checks,
        },
        "decision": {
            "compiled_formula_reconciled_to_published_truncated_architecture": (
                module_formula == EXPECTED_HOST_FORMULA
                and all_checks["host_formula_matches"]
            ),
            "two_cap_three_pillar_six_urea_architecture_reconciled": all(
                all_checks[key]
                for key in (
                    "five_aromatic_six_member_rings",
                    "two_caps",
                    "three_pillars",
                    "six_urea_carbonyls",
                    "three_amide_carbonyls",
                    "no_other_oxygenated_carbons",
                    "both_caps_have_three_ethyl_three_methylene_pattern",
                    "all_pillars_have_one_N_methyl_benzamide",
                    "ureas_form_complete_K2_3_cap_pillar_graph",
                )
            ),
            "compiled_graph_matches_published_truncated_host": all(all_checks.values()),
            "exact_author_atom_names_and_source_build_files_available": False,
            "resp_gaff_derivation_reconstructable": False,
            "force_field_chemically_validated": False,
            "parent_pag_calibration_authorized": False,
            "mcq1_parameter_transfer_authorized": False,
            "candidate_ranking_authorized": False,
            "next_action": (
                "request source mol2/RESP/frcmod/tleap and starting coordinates; validate "
                "the compiled charges and GAFF terms against those files before parent PaG work"
            ),
        },
        "limitations": [
            "The published architecture is reconciled at molecular-graph and formula level, not by an author-supplied atom-name map.",
            "AMBER atom types and numerical terms remain compiled outputs whose source derivation is unavailable.",
            "Formula and motif agreement do not validate partial charges, protonation, conformational sampling, or transferability.",
            "The simulated object is the methyl-truncated public host, not full receptor 2 or the NNC2215 conjugate."
        ],
    }
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument(
        "--output-json",
        type=Path,
        default=Path(
            "research/data/gluhut-truncated-architecture-reconciliation.json"
        ),
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    report = audit(args.archive, args.output_json)
    print(
        json.dumps(
            {
                "module_bookkeeping": report["published_module_formula_bookkeeping"],
                "cross_system": report["cross_system"],
                "decision": report["decision"],
                "reference_system": report["systems"][0],
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
