#!/usr/bin/env python3
"""Audit molecular graph chemistry in the deposited GluHUT AMBER systems.

This script resolves residue connectivity, formulas, partial-charge partitions,
compiled host motifs and graph invariants from the public SMD topologies. It is
an audit of deposited compiled systems. It does not validate RESP derivation,
predict affinity, or authorize transfer to PaG, MC-Q1 or NNC2215.
"""
from __future__ import annotations

import argparse
import collections
import hashlib
import json
import math
import stat
import sys
import zipfile
from pathlib import Path, PurePosixPath
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent))
import audit_gluhut_compiled_systems as amber  # noqa: E402

ARCHIVE_MD5 = "9f3c6733dc608d2bd37c488c560a5443"
ARCHIVE_SIZE_BYTES = 37931874
CHUNK_SIZE = 1024 * 1024
CHARGE_TOLERANCE_E = 1.0e-5

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


class ChemistryAuditError(RuntimeError):
    pass


def digest_file(path: Path, algorithm: str) -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(CHUNK_SIZE), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


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
        raise ChemistryAuditError(f"backslash in member path: {name!r}")
    path = PurePosixPath(name)
    if path.is_absolute() or any(part in {"", ".", ".."} for part in path.parts):
        raise ChemistryAuditError(f"unsafe member path: {name!r}")
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


def formula_from_atomic_numbers(atomic_numbers: list[int]) -> tuple[str, dict[str, int]]:
    counts: collections.Counter[str] = collections.Counter()
    for atomic_number in atomic_numbers:
        try:
            symbol = ELEMENT_SYMBOLS[int(atomic_number)]
        except KeyError as exc:
            raise ChemistryAuditError(
                f"unsupported atomic number in formula: {atomic_number}"
            ) from exc
        counts[symbol] += 1
    order: list[str] = []
    if counts.get("C"):
        order.append("C")
        if counts.get("H"):
            order.append("H")
    order.extend(sorted(symbol for symbol in counts if symbol not in order))
    formula = "".join(
        symbol + (str(counts[symbol]) if counts[symbol] != 1 else "")
        for symbol in order
    )
    return formula, dict(sorted(counts.items()))


def residue_spans(
    sections: dict[str, list[Any]], pointers: dict[str, int]
) -> list[dict[str, Any]]:
    natom = pointers["NATOM"]
    nres = pointers["NRES"]
    labels = [str(value) for value in amber.required(sections, "RESIDUE_LABEL", nres)[:nres]]
    starts = [int(value) for value in amber.required(sections, "RESIDUE_POINTER", nres)[:nres]]
    if not starts or starts[0] != 1:
        raise ChemistryAuditError("first residue does not begin at atom 1")
    rows: list[dict[str, Any]] = []
    for index, (label, start) in enumerate(zip(labels, starts), start=1):
        end = starts[index] - 1 if index < len(starts) else natom
        rows.append(
            {
                "index_1_based": index,
                "label": label,
                "first_atom_1_based": start,
                "last_atom_1_based": end,
                "atom_count": end - start + 1,
            }
        )
    return rows


def atom_to_residue_map(residues: list[dict[str, Any]], natom: int) -> list[int]:
    result = [0] * (natom + 1)
    for residue in residues:
        for atom in range(
            residue["first_atom_1_based"], residue["last_atom_1_based"] + 1
        ):
            if result[atom]:
                raise ChemistryAuditError(f"atom {atom} belongs to multiple residues")
            result[atom] = residue["index_1_based"]
    if any(value == 0 for value in result[1:]):
        raise ChemistryAuditError("one or more atoms have no residue assignment")
    return result


def all_bonds(
    sections: dict[str, list[Any]], pointers: dict[str, int]
) -> list[dict[str, Any]]:
    natom = pointers["NATOM"]
    force_constants = [float(value) for value in sections.get("BOND_FORCE_CONSTANT", [])]
    equilibrium = [float(value) for value in sections.get("BOND_EQUIL_VALUE", [])]
    rows: list[dict[str, Any]] = []
    seen: set[tuple[int, int]] = set()
    for flag in ("BONDS_INC_HYDROGEN", "BONDS_WITHOUT_HYDROGEN"):
        values = [int(value) for value in sections.get(flag, [])]
        if len(values) % 3:
            raise ChemistryAuditError(f"{flag} is not divisible by three")
        for offset in range(0, len(values), 3):
            left = abs(values[offset]) // 3 + 1
            right = abs(values[offset + 1]) // 3 + 1
            parameter = values[offset + 2]
            if not (1 <= left <= natom and 1 <= right <= natom):
                raise ChemistryAuditError(f"bond atom outside range in {flag}")
            key = tuple(sorted((left, right)))
            if key in seen:
                raise ChemistryAuditError(f"duplicate bond {key}")
            seen.add(key)
            if not 1 <= parameter <= len(force_constants):
                raise ChemistryAuditError(f"bond parameter outside table in {flag}")
            rows.append(
                {
                    "atoms_1_based": list(key),
                    "parameter_index_1_based": parameter,
                    "force_constant": force_constants[parameter - 1],
                    "equilibrium_A": equilibrium[parameter - 1],
                    "source_flag": flag,
                }
            )
    return sorted(rows, key=lambda row: row["atoms_1_based"])


def adjacency_from_bonds(natom: int, bonds: list[dict[str, Any]]) -> list[set[int]]:
    adjacency = [set() for _ in range(natom + 1)]
    for bond in bonds:
        left, right = bond["atoms_1_based"]
        adjacency[left].add(right)
        adjacency[right].add(left)
    return adjacency


def components(atom_indices: set[int], adjacency: list[set[int]]) -> list[list[int]]:
    remaining = set(atom_indices)
    result: list[list[int]] = []
    while remaining:
        start = min(remaining)
        stack = [start]
        component: list[int] = []
        remaining.remove(start)
        while stack:
            atom = stack.pop()
            component.append(atom)
            for neighbour in sorted(adjacency[atom], reverse=True):
                if neighbour in remaining and neighbour in atom_indices:
                    remaining.remove(neighbour)
                    stack.append(neighbour)
        result.append(sorted(component))
    return sorted(result, key=lambda row: (len(row), row), reverse=True)


def residue_summary(
    residue: dict[str, Any],
    atom_names: list[str],
    atom_types: list[str],
    atomic_numbers: list[int],
    charges_e: list[float],
) -> dict[str, Any]:
    start = residue["first_atom_1_based"] - 1
    end = residue["last_atom_1_based"]
    formula, element_counts = formula_from_atomic_numbers(atomic_numbers[start:end])
    return {
        **residue,
        "formula": formula,
        "element_counts": element_counts,
        "charge_e": sum(charges_e[start:end]),
        "atom_names": atom_names[start:end],
        "amber_atom_types": atom_types[start:end],
    }


def classify_carbonyl_motifs(
    host_atoms: set[int],
    adjacency: list[set[int]],
    bonds_by_pair: dict[tuple[int, int], dict[str, Any]],
    atomic_numbers_1: list[int],
    atom_names_1: list[str],
    atom_types_1: list[str],
) -> list[dict[str, Any]]:
    motifs: list[dict[str, Any]] = []
    for atom in sorted(host_atoms):
        if atomic_numbers_1[atom] != 6:
            continue
        oxygen_neighbours = [
            neighbour
            for neighbour in adjacency[atom]
            if neighbour in host_atoms and atomic_numbers_1[neighbour] == 8
        ]
        carbonyl_oxygens = []
        for oxygen in oxygen_neighbours:
            bond = bonds_by_pair[tuple(sorted((atom, oxygen)))]
            if bond["equilibrium_A"] <= 1.30:
                carbonyl_oxygens.append(oxygen)
        if not carbonyl_oxygens:
            continue
        heavy_neighbours = [
            neighbour
            for neighbour in adjacency[atom]
            if neighbour in host_atoms
            and atomic_numbers_1[neighbour] != 1
            and neighbour not in carbonyl_oxygens
        ]
        heavy_elements = sorted(ELEMENT_SYMBOLS[atomic_numbers_1[n]] for n in heavy_neighbours)
        if heavy_elements == ["N", "N"]:
            motif = "urea_carbonyl"
        elif heavy_elements == ["C", "N"]:
            motif = "amide_carbonyl"
        elif heavy_elements == ["C"]:
            motif = "aldehyde_or_acyl_carbonyl"
        elif heavy_elements == ["O", "O"]:
            motif = "carbonate_like_carbonyl"
        else:
            motif = "other_carbonyl"
        motifs.append(
            {
                "carbon_atom_1_based": atom,
                "carbon_atom_name": atom_names_1[atom],
                "carbon_atom_type": atom_types_1[atom],
                "carbonyl_oxygen_atoms_1_based": carbonyl_oxygens,
                "heavy_neighbour_atoms_1_based": heavy_neighbours,
                "heavy_neighbour_elements": heavy_elements,
                "classification": motif,
            }
        )
    return motifs


def methyl_groups(
    host_atoms: set[int],
    adjacency: list[set[int]],
    atomic_numbers_1: list[int],
    atom_names_1: list[str],
    atom_types_1: list[str],
) -> list[dict[str, Any]]:
    rows = []
    for atom in sorted(host_atoms):
        if atomic_numbers_1[atom] != 6:
            continue
        neighbours = [n for n in adjacency[atom] if n in host_atoms]
        hydrogen_count = sum(atomic_numbers_1[n] == 1 for n in neighbours)
        heavy = [n for n in neighbours if atomic_numbers_1[n] != 1]
        if hydrogen_count == 3 and len(heavy) == 1:
            rows.append(
                {
                    "carbon_atom_1_based": atom,
                    "atom_name": atom_names_1[atom],
                    "amber_atom_type": atom_types_1[atom],
                    "attached_heavy_atom_1_based": heavy[0],
                    "attached_heavy_element": ELEMENT_SYMBOLS[
                        atomic_numbers_1[heavy[0]]
                    ],
                }
            )
    return rows


def analyse_topology(path: str, payload: bytes) -> dict[str, Any]:
    sections = amber.parse_parm7_text(payload.decode("utf-8"))
    pointers = amber.pointer_map(sections)
    natom = pointers["NATOM"]
    atom_names = [str(v) for v in amber.required(sections, "ATOM_NAME", natom)[:natom]]
    atom_types = [
        str(v) for v in amber.required(sections, "AMBER_ATOM_TYPE", natom)[:natom]
    ]
    atomic_numbers = [
        int(v) for v in amber.required(sections, "ATOMIC_NUMBER", natom)[:natom]
    ]
    raw_charges = [float(v) for v in amber.required(sections, "CHARGE", natom)[:natom]]
    charges_e = [value / amber.AMBER_CHARGE_SCALE for value in raw_charges]
    residues = residue_spans(sections, pointers)
    residue_by_atom = atom_to_residue_map(residues, natom)
    bonds = all_bonds(sections, pointers)
    adjacency = adjacency_from_bonds(natom, bonds)
    bonds_by_pair = {
        tuple(bond["atoms_1_based"]): bond for bond in bonds
    }

    residue_rows = [
        residue_summary(
            residue,
            atom_names,
            atom_types,
            atomic_numbers,
            charges_e,
        )
        for residue in residues
    ]
    host_residue = max(residue_rows, key=lambda row: row["atom_count"])
    host_index = host_residue["index_1_based"]
    host_atoms = {
        atom for atom in range(1, natom + 1) if residue_by_atom[atom] == host_index
    }
    guest_atoms = set(range(1, natom + 1)) - host_atoms

    cross_residue_bonds = []
    for bond in bonds:
        left, right = bond["atoms_1_based"]
        left_res = residue_by_atom[left]
        right_res = residue_by_atom[right]
        if left_res == right_res:
            continue
        cross_residue_bonds.append(
            {
                **bond,
                "left": {
                    "atom_1_based": left,
                    "atom_name": atom_names[left - 1],
                    "element": ELEMENT_SYMBOLS[atomic_numbers[left - 1]],
                    "amber_atom_type": atom_types[left - 1],
                    "residue_index_1_based": left_res,
                    "residue_label": residues[left_res - 1]["label"],
                },
                "right": {
                    "atom_1_based": right,
                    "atom_name": atom_names[right - 1],
                    "element": ELEMENT_SYMBOLS[atomic_numbers[right - 1]],
                    "amber_atom_type": atom_types[right - 1],
                    "residue_index_1_based": right_res,
                    "residue_label": residues[right_res - 1]["label"],
                },
            }
        )

    combined_guest_atomic_numbers = [
        atomic_numbers[atom - 1] for atom in sorted(guest_atoms)
    ]
    guest_formula, guest_element_counts = formula_from_atomic_numbers(
        combined_guest_atomic_numbers
    )
    guest_charge = sum(charges_e[atom - 1] for atom in guest_atoms)
    guest_components = components(guest_atoms, adjacency)
    host_components = components(host_atoms, adjacency)

    host_formula, host_element_counts = formula_from_atomic_numbers(
        [atomic_numbers[atom - 1] for atom in sorted(host_atoms)]
    )
    host_bonds = [
        bond
        for bond in bonds
        if set(bond["atoms_1_based"]).issubset(host_atoms)
    ]
    host_heavy_atoms = {
        atom for atom in host_atoms if atomic_numbers[atom - 1] != 1
    }
    host_heavy_bond_count = sum(
        left in host_heavy_atoms and right in host_heavy_atoms
        for left, right in (bond["atoms_1_based"] for bond in host_bonds)
    )
    host_heavy_components = components(host_heavy_atoms, adjacency)
    host_cycle_rank = (
        host_heavy_bond_count
        - len(host_heavy_atoms)
        + len(host_heavy_components)
    )

    atomic_numbers_1 = [0] + atomic_numbers
    atom_names_1 = [""] + atom_names
    atom_types_1 = [""] + atom_types
    motifs = classify_carbonyl_motifs(
        host_atoms,
        adjacency,
        bonds_by_pair,
        atomic_numbers_1,
        atom_names_1,
        atom_types_1,
    )
    motif_counts = collections.Counter(row["classification"] for row in motifs)
    methyls = methyl_groups(
        host_atoms,
        adjacency,
        atomic_numbers_1,
        atom_names_1,
        atom_types_1,
    )
    host_type_counts = collections.Counter(
        atom_types[atom - 1] for atom in host_atoms
    )
    host_heavy_degree_counts = collections.Counter(
        sum(
            neighbour in host_heavy_atoms
            for neighbour in adjacency[atom]
        )
        for atom in host_heavy_atoms
    )

    graph_nodes = [
        {
            "atom_1_based": atom,
            "name": atom_names[atom - 1],
            "atomic_number": atomic_numbers[atom - 1],
            "amber_atom_type": atom_types[atom - 1],
            "charge_e": charges_e[atom - 1],
            "residue_index_1_based": residue_by_atom[atom],
            "residue_label": residues[residue_by_atom[atom] - 1]["label"],
        }
        for atom in range(1, natom + 1)
    ]
    graph_edges = [
        {
            "atoms_1_based": bond["atoms_1_based"],
            "equilibrium_A": bond["equilibrium_A"],
            "force_constant": bond["force_constant"],
        }
        for bond in bonds
    ]

    return {
        "path": path,
        "sha256": sha256_bytes(payload),
        "natom": natom,
        "total_charge_e": sum(charges_e),
        "residues": residue_rows,
        "cross_residue_bonds": cross_residue_bonds,
        "host": {
            "residue_index_1_based": host_index,
            "residue_label": host_residue["label"],
            "atom_count": len(host_atoms),
            "heavy_atom_count": len(host_heavy_atoms),
            "formula": host_formula,
            "element_counts": host_element_counts,
            "charge_e": sum(charges_e[atom - 1] for atom in host_atoms),
            "connected_component_count": len(host_components),
            "heavy_graph_bond_count": host_heavy_bond_count,
            "heavy_graph_cycle_rank": host_cycle_rank,
            "amber_atom_type_counts": dict(sorted(host_type_counts.items())),
            "heavy_atom_degree_counts": {
                str(key): value for key, value in sorted(host_heavy_degree_counts.items())
            },
            "carbonyl_motifs": motifs,
            "carbonyl_motif_counts": dict(sorted(motif_counts.items())),
            "methyl_groups": methyls,
            "methyl_group_count": len(methyls),
        },
        "combined_guest": {
            "atom_count": len(guest_atoms),
            "formula": guest_formula,
            "element_counts": guest_element_counts,
            "charge_e": guest_charge,
            "connected_component_count": len(guest_components),
            "component_atom_counts": [len(row) for row in guest_components],
            "residue_indices_1_based": sorted(
                {residue_by_atom[atom] for atom in guest_atoms}
            ),
            "residue_labels": [
                residues[index - 1]["label"]
                for index in sorted({residue_by_atom[atom] for atom in guest_atoms})
            ],
            "is_neutral_within_tolerance": abs(guest_charge) <= CHARGE_TOLERANCE_E,
            "is_single_covalent_component": len(guest_components) == 1,
        },
        "graph": {
            "node_count": len(graph_nodes),
            "edge_count": len(graph_edges),
            "nodes_sha256": canonical_hash(graph_nodes),
            "edges_sha256": canonical_hash(graph_edges),
            "full_graph_sha256": canonical_hash(
                {"nodes": graph_nodes, "edges": graph_edges}
            ),
        },
    }


def infer_system(path: str) -> str:
    lower = path.lower()
    for name in ("fructofuranose", "fructopyranose", "galactose", "glucose"):
        if name in lower:
            return name
    return PurePosixPath(path).parent.name.lower()


def audit(archive_path: Path, output_json: Path) -> dict[str, Any]:
    if archive_path.stat().st_size != ARCHIVE_SIZE_BYTES:
        raise ChemistryAuditError(
            f"size mismatch: expected {ARCHIVE_SIZE_BYTES}, got {archive_path.stat().st_size}"
        )
    md5 = digest_file(archive_path, "md5")
    if md5 != ARCHIVE_MD5:
        raise ChemistryAuditError(f"MD5 mismatch: expected {ARCHIVE_MD5}, got {md5}")
    archive_sha256 = digest_file(archive_path, "sha256")

    systems: list[dict[str, Any]] = []
    seen: set[str] = set()
    with zipfile.ZipFile(archive_path) as archive:
        bad = archive.testzip()
        if bad:
            raise ChemistryAuditError(f"ZIP CRC failure: {bad}")
        for info in archive.infolist():
            safe_member(info.filename)
            if info.filename in seen:
                raise ChemistryAuditError(f"duplicate member: {info.filename}")
            seen.add(info.filename)
            if info.flag_bits & 0x1:
                raise ChemistryAuditError(f"encrypted member: {info.filename}")
            kind = member_kind(info)
            if kind in {"symlink", "special"}:
                raise ChemistryAuditError(f"unsupported {kind} member: {info.filename}")
            if kind != "file" or PurePosixPath(info.filename).suffix.lower() not in {
                ".parm7",
                ".prmtop",
            }:
                continue
            system = analyse_topology(info.filename, archive.read(info))
            system["system_name"] = infer_system(info.filename)
            systems.append(system)

    systems.sort(key=lambda row: row["system_name"])
    if len(systems) != 4:
        raise ChemistryAuditError(f"expected four topologies, found {len(systems)}")

    reference = systems[0]
    host_formulas = {row["host"]["formula"] for row in systems}
    host_motif_counts = {
        canonical_hash(row["host"]["carbonyl_motif_counts"]) for row in systems
    }
    guest_formulas = {row["combined_guest"]["formula"] for row in systems}
    guest_neutral = all(
        row["combined_guest"]["is_neutral_within_tolerance"] for row in systems
    )
    guest_connected = all(
        row["combined_guest"]["is_single_covalent_component"] for row in systems
    )
    roh_covalent_bonds = []
    for system in systems:
        matching = [
            bond
            for bond in system["cross_residue_bonds"]
            if "ROH" in {bond["left"]["residue_label"], bond["right"]["residue_label"]}
        ]
        roh_covalent_bonds.append(
            {
                "system_name": system["system_name"],
                "bonds": matching,
                "bond_count": len(matching),
            }
        )

    report = {
        "schema": "dia-health-one.gluhut-graph-chemistry-audit.v1",
        "status": "compiled_graph_chemistry_audit_not_source_parameter_validation",
        "scope": "four public GluHUT sugar benchmark topologies",
        "not_target_claim": (
            "not NNC2215, PaG, MC-Q1, affinity, efficacy, safety, or a drug result"
        ),
        "archive": {
            "filename": archive_path.name,
            "size_bytes": archive_path.stat().st_size,
            "md5": md5,
            "sha256": archive_sha256,
            "crc_verified": True,
        },
        "systems": systems,
        "cross_system_summary": {
            "host_formula_identical": len(host_formulas) == 1,
            "host_formulas": sorted(host_formulas),
            "host_carbonyl_motif_counts_identical": len(host_motif_counts) == 1,
            "reference_host_carbonyl_motif_counts": reference["host"][
                "carbonyl_motif_counts"
            ],
            "guest_formulas": sorted(guest_formulas),
            "all_combined_guests_neutral": guest_neutral,
            "all_combined_guests_single_covalent_component": guest_connected,
            "roh_cross_residue_covalent_bonds": roh_covalent_bonds,
        },
        "decision": {
            "roh_is_counterion": False if guest_connected else None,
            "roh_is_covalently_part_of_guest": guest_connected
            and all(row["bond_count"] >= 1 for row in roh_covalent_bonds),
            "combined_guest_formula_matches_hexose": guest_formulas == {"C6H12O6"},
            "compiled_host_graph_can_be_fingerprinted": True,
            "source_graph_fully_reconciled": False,
            "candidate_ranking_authorized": False,
            "pag_parameterization_authorized": False,
            "mcq1_parameter_transfer_authorized": False,
            "next_action": (
                "compare the compiled host formula, six-urea/amide motif counts and graph "
                "fingerprint with the published truncated GluHUT drawing and source synthesis; "
                "request source mol2/frcmod/RESP and exact starting coordinates"
            ),
        },
        "limitations": [
            "AMBER topology bonds and atom types define a compiled graph but do not preserve the original quantum-chemistry or RESP derivation.",
            "Bond order is not explicitly reconstructed; carbonyl motifs use equilibrium bond length and neighbour elements.",
            "The host graph is a truncated public GluHUT benchmark, not the complete NNC2215 conjugate.",
            "Internal consistency and formula agreement do not prove chemical correctness or transferability.",
        ],
    }
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return report


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument(
        "--output-json",
        type=Path,
        default=Path("research/data/gluhut-graph-chemistry-audit.json"),
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    report = audit(args.archive, args.output_json)
    print(
        json.dumps(
            {
                "cross_system_summary": report["cross_system_summary"],
                "decision": report["decision"],
                "host": report["systems"][0]["host"],
                "guests": [
                    {
                        "system_name": row["system_name"],
                        "residues": row["residues"],
                        "cross_residue_bonds": row["cross_residue_bonds"],
                        "combined_guest": row["combined_guest"],
                    }
                    for row in report["systems"]
                ],
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
