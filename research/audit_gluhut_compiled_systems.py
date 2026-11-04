#!/usr/bin/env python3
"""Audit compiled AMBER GluHUT host–sugar systems from the public SMD archive.

The audit parses deposited parm7 topologies and restart coordinates, identifies
the largest residue as the synthetic host, compares host atom properties,
nonbonded coefficients, exclusions, and bonded terms across sugars, and records
what is reusable for benchmark reproduction.

It does not reconstruct RESP/GAFF/GLYCAM parameter generation, calculate an
affinity, validate chemical correctness, or authorize PaG/MC-Q1 ranking.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import stat
import zipfile
from pathlib import Path, PurePosixPath
from typing import Any

ARCHIVE_MD5 = "9f3c6733dc608d2bd37c488c560a5443"
ARCHIVE_SIZE_BYTES = 37931874
AMBER_CHARGE_SCALE = 18.2223
CHUNK_SIZE = 1024 * 1024
FLOAT_TOLERANCE = 1.0e-9

POINTER_NAMES = [
    "NATOM", "NTYPES", "NBONH", "MBONA", "NTHETH", "MTHETA", "NPHIH",
    "MPHIA", "NHPARM", "NPARM", "NNB", "NRES", "NBONA", "NTHETA",
    "NPHIA", "NUMBND", "NUMANG", "NPTRA", "NATYP", "NPHB", "IFPERT",
    "NBPER", "NGPER", "NDPER", "MBPER", "MGPER", "MDPER", "IFBOX",
    "NMXRS", "IFCAP", "NUMEXTRA", "NCOPY",
]
FORMAT_RE = re.compile(
    r"\(\s*(?:(?P<scale>[+-]?\d+)P)?(?P<count>\d+)"
    r"(?P<kind>[aAiIeEfF])(?P<width>\d+)"
    r"(?:\.(?P<precision>\d+))?\s*\)"
)


class AuditError(RuntimeError):
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
    payload = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")
    return sha256_bytes(payload)


def safe_member(name: str) -> PurePosixPath:
    if "\\" in name:
        raise AuditError(f"backslash in member path: {name!r}")
    path = PurePosixPath(name)
    if path.is_absolute() or any(part in {"", ".", ".."} for part in path.parts):
        raise AuditError(f"unsafe member path: {name!r}")
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


def parse_fortran_format(line: str) -> tuple[int, str, int]:
    match = FORMAT_RE.fullmatch(line.strip().removeprefix("%FORMAT"))
    if not match:
        raise AuditError(f"unsupported AMBER format line: {line!r}")
    return (
        int(match.group("count")),
        match.group("kind").upper(),
        int(match.group("width")),
    )


def parse_fixed_width_lines(lines: list[str], kind: str, width: int) -> list[Any]:
    values: list[Any] = []
    for line in lines:
        for start in range(0, len(line), width):
            field = line[start : start + width]
            if not field.strip():
                continue
            if kind == "A":
                values.append(field.strip())
            elif kind == "I":
                values.append(int(field))
            elif kind in {"E", "F"}:
                values.append(float(field.replace("D", "E").replace("d", "e")))
            else:
                raise AuditError(f"unsupported field kind: {kind}")
    return values


def parse_parm7_text(text: str) -> dict[str, list[Any]]:
    lines = text.splitlines()
    sections: dict[str, list[Any]] = {}
    index = 0
    while index < len(lines):
        line = lines[index]
        if not line.startswith("%FLAG "):
            index += 1
            continue
        flag = line[6:].strip()
        index += 1
        while index < len(lines) and lines[index].startswith("%COMMENT"):
            index += 1
        if index >= len(lines) or not lines[index].startswith("%FORMAT"):
            raise AuditError(f"missing %FORMAT for {flag}")
        _, kind, width = parse_fortran_format(lines[index])
        index += 1
        data_lines: list[str] = []
        while index < len(lines) and not lines[index].startswith("%FLAG "):
            if not lines[index].startswith("%COMMENT"):
                data_lines.append(lines[index])
            index += 1
        sections[flag] = parse_fixed_width_lines(data_lines, kind, width)
    if "POINTERS" not in sections:
        raise AuditError("parm7 is missing POINTERS")
    return sections


def required(sections: dict[str, list[Any]], flag: str, count: int | None = None) -> list[Any]:
    if flag not in sections:
        raise AuditError(f"parm7 is missing required flag {flag}")
    values = sections[flag]
    if count is not None and len(values) < count:
        raise AuditError(f"{flag}: expected at least {count} entries, found {len(values)}")
    return values


def pointer_map(sections: dict[str, list[Any]]) -> dict[str, int]:
    pointers = [int(value) for value in required(sections, "POINTERS", 31)]
    return {
        name: pointers[index]
        for index, name in enumerate(POINTER_NAMES)
        if index < len(pointers)
    }


def residue_records(
    sections: dict[str, list[Any]],
    pointers: dict[str, int],
    charges_e: list[float],
) -> list[dict[str, Any]]:
    natom = pointers["NATOM"]
    nres = pointers["NRES"]
    labels = [str(value) for value in required(sections, "RESIDUE_LABEL", nres)[:nres]]
    starts = [int(value) for value in required(sections, "RESIDUE_POINTER", nres)[:nres]]
    if not starts or starts[0] != 1:
        raise AuditError(f"unexpected first residue pointer: {starts[:1]}")
    if any(value < 1 or value > natom for value in starts):
        raise AuditError("residue pointer outside atom range")
    if any(right <= left for left, right in zip(starts, starts[1:])):
        raise AuditError("residue pointers are not strictly increasing")
    rows = []
    for index, (label, start) in enumerate(zip(labels, starts)):
        end = starts[index + 1] - 1 if index + 1 < len(starts) else natom
        charge = sum(charges_e[start - 1 : end])
        rows.append(
            {
                "index_1_based": index + 1,
                "label": label,
                "first_atom_1_based": start,
                "last_atom_1_based": end,
                "atom_count": end - start + 1,
                "charge_e": charge,
            }
        )
    return rows


def connectivity_records(
    sections: dict[str, list[Any]],
    host_atoms: set[int],
) -> dict[str, list[dict[str, Any]]]:
    result: dict[str, list[dict[str, Any]]] = {}

    bond_k = [float(v) for v in sections.get("BOND_FORCE_CONSTANT", [])]
    bond_r0 = [float(v) for v in sections.get("BOND_EQUIL_VALUE", [])]
    bonds: list[dict[str, Any]] = []
    for flag in ("BONDS_INC_HYDROGEN", "BONDS_WITHOUT_HYDROGEN"):
        values = [int(v) for v in sections.get(flag, [])]
        if len(values) % 3:
            raise AuditError(f"{flag} length is not divisible by 3")
        for offset in range(0, len(values), 3):
            a = abs(values[offset]) // 3 + 1
            b = abs(values[offset + 1]) // 3 + 1
            parameter = values[offset + 2]
            if a in host_atoms and b in host_atoms:
                if not 1 <= parameter <= len(bond_k):
                    raise AuditError(f"{flag} bond parameter index outside table")
                atoms = sorted((a, b))
                bonds.append(
                    {
                        "atoms": atoms,
                        "force_constant": bond_k[parameter - 1],
                        "equilibrium_A": bond_r0[parameter - 1],
                    }
                )
    result["bonds"] = sorted(
        bonds, key=lambda row: (row["atoms"], row["force_constant"], row["equilibrium_A"])
    )

    angle_k = [float(v) for v in sections.get("ANGLE_FORCE_CONSTANT", [])]
    angle_t0 = [float(v) for v in sections.get("ANGLE_EQUIL_VALUE", [])]
    angles: list[dict[str, Any]] = []
    for flag in ("ANGLES_INC_HYDROGEN", "ANGLES_WITHOUT_HYDROGEN"):
        values = [int(v) for v in sections.get(flag, [])]
        if len(values) % 4:
            raise AuditError(f"{flag} length is not divisible by 4")
        for offset in range(0, len(values), 4):
            atoms = tuple(abs(values[offset + i]) // 3 + 1 for i in range(3))
            parameter = values[offset + 3]
            if all(atom in host_atoms for atom in atoms):
                if not 1 <= parameter <= len(angle_k):
                    raise AuditError(f"{flag} angle parameter index outside table")
                reverse = tuple(reversed(atoms))
                canonical = min(atoms, reverse)
                angles.append(
                    {
                        "atoms": list(canonical),
                        "force_constant": angle_k[parameter - 1],
                        "equilibrium_rad": angle_t0[parameter - 1],
                    }
                )
    result["angles"] = sorted(
        angles, key=lambda row: (row["atoms"], row["force_constant"], row["equilibrium_rad"])
    )

    dih_k = [float(v) for v in sections.get("DIHEDRAL_FORCE_CONSTANT", [])]
    dih_period = [float(v) for v in sections.get("DIHEDRAL_PERIODICITY", [])]
    dih_phase = [float(v) for v in sections.get("DIHEDRAL_PHASE", [])]
    scee = [float(v) for v in sections.get("SCEE_SCALE_FACTOR", [])]
    scnb = [float(v) for v in sections.get("SCNB_SCALE_FACTOR", [])]
    dihedrals: list[dict[str, Any]] = []
    for flag in ("DIHEDRALS_INC_HYDROGEN", "DIHEDRALS_WITHOUT_HYDROGEN"):
        values = [int(v) for v in sections.get(flag, [])]
        if len(values) % 5:
            raise AuditError(f"{flag} length is not divisible by 5")
        for offset in range(0, len(values), 5):
            raw_atoms = values[offset : offset + 4]
            atoms = tuple(abs(value) // 3 + 1 for value in raw_atoms)
            parameter = values[offset + 4]
            if all(atom in host_atoms for atom in atoms):
                if not 1 <= parameter <= len(dih_k):
                    raise AuditError(f"{flag} dihedral parameter index outside table")
                dihedrals.append(
                    {
                        "atoms": list(atoms),
                        "third_atom_negative": raw_atoms[2] < 0,
                        "fourth_atom_negative": raw_atoms[3] < 0,
                        "force_constant": dih_k[parameter - 1],
                        "periodicity": dih_period[parameter - 1],
                        "phase_rad": dih_phase[parameter - 1],
                        "scee": scee[parameter - 1] if parameter <= len(scee) else None,
                        "scnb": scnb[parameter - 1] if parameter <= len(scnb) else None,
                    }
                )
    result["dihedrals"] = sorted(
        dihedrals,
        key=lambda row: (
            row["atoms"],
            row["third_atom_negative"],
            row["fourth_atom_negative"],
            row["force_constant"],
            row["periodicity"],
            row["phase_rad"],
        ),
    )
    return result


def exclusion_records(
    sections: dict[str, list[Any]],
    natom: int,
    host_atoms: set[int],
) -> list[list[int]]:
    counts = [int(v) for v in required(sections, "NUMBER_EXCLUDED_ATOMS", natom)[:natom]]
    flat = [int(v) for v in sections.get("EXCLUDED_ATOMS_LIST", [])]
    cursor = 0
    result: list[list[int]] = []
    for atom_index, count in enumerate(counts, start=1):
        values = flat[cursor : cursor + count]
        cursor += count
        if atom_index in host_atoms:
            result.append(sorted(value for value in values if value in host_atoms and value != 0))
    if cursor != len(flat):
        raise AuditError(
            f"EXCLUDED_ATOMS_LIST consumption mismatch: used {cursor}, found {len(flat)}"
        )
    return result


def lj_records(
    sections: dict[str, list[Any]],
    pointers: dict[str, int],
    host_slice: slice,
) -> dict[str, Any]:
    natom = pointers["NATOM"]
    ntypes = pointers["NTYPES"]
    type_indices = [
        int(v) for v in required(sections, "ATOM_TYPE_INDEX", natom)[:natom]
    ]
    atom_types = [
        str(v) for v in required(sections, "AMBER_ATOM_TYPE", natom)[:natom]
    ]
    matrix = [
        int(v)
        for v in required(
            sections, "NONBONDED_PARM_INDEX", ntypes * ntypes
        )[: ntypes * ntypes]
    ]
    acoef = [float(v) for v in sections.get("LENNARD_JONES_ACOEF", [])]
    bcoef = [float(v) for v in sections.get("LENNARD_JONES_BCOEF", [])]

    name_to_indices: dict[str, set[int]] = {}
    for name, index in zip(atom_types[host_slice], type_indices[host_slice]):
        name_to_indices.setdefault(name, set()).add(index)
    inconsistent = {
        name: sorted(indices)
        for name, indices in name_to_indices.items()
        if len(indices) != 1
    }
    if inconsistent:
        raise AuditError(f"host atom type names map to multiple indices: {inconsistent}")
    name_to_index = {name: next(iter(indices)) for name, indices in name_to_indices.items()}

    names = sorted(name_to_index)
    rows: list[dict[str, Any]] = []
    for left_pos, left_name in enumerate(names):
        for right_name in names[left_pos:]:
            left_index = name_to_index[left_name]
            right_index = name_to_index[right_name]
            parameter_index = matrix[(left_index - 1) * ntypes + (right_index - 1)]
            row: dict[str, Any] = {
                "types": [left_name, right_name],
                "type_indices": [left_index, right_index],
                "nonbonded_parameter_index": parameter_index,
            }
            if parameter_index > 0:
                if parameter_index > len(acoef) or parameter_index > len(bcoef):
                    raise AuditError("Lennard-Jones parameter index outside A/B tables")
                row["acoef"] = acoef[parameter_index - 1]
                row["bcoef"] = bcoef[parameter_index - 1]
            else:
                row["acoef"] = None
                row["bcoef"] = None
            rows.append(row)
    physical_rows = [
        {
            "types": row["types"],
            "acoef": row["acoef"],
            "bcoef": row["bcoef"],
        }
        for row in rows
    ]
    return {
        "indexed": rows,
        "physical": physical_rows,
        "host_type_name_to_global_index": {
            name: name_to_index[name] for name in sorted(name_to_index)
        },
    }


def topology_summary(path: str, payload: bytes) -> dict[str, Any]:
    try:
        text = payload.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise AuditError(f"{path} is not UTF-8 text") from exc
    sections = parse_parm7_text(text)
    pointers = pointer_map(sections)
    natom = pointers["NATOM"]
    atom_names = [str(v) for v in required(sections, "ATOM_NAME", natom)[:natom]]
    amber_types = [
        str(v) for v in required(sections, "AMBER_ATOM_TYPE", natom)[:natom]
    ]
    type_indices = [
        int(v) for v in required(sections, "ATOM_TYPE_INDEX", natom)[:natom]
    ]
    raw_charges = [float(v) for v in required(sections, "CHARGE", natom)[:natom]]
    charges_e = [value / AMBER_CHARGE_SCALE for value in raw_charges]
    masses = [float(v) for v in required(sections, "MASS", natom)[:natom]]
    atomic_numbers = [
        int(v) for v in sections.get("ATOMIC_NUMBER", [0] * natom)[:natom]
    ]
    if len(atomic_numbers) < natom:
        atomic_numbers.extend([0] * (natom - len(atomic_numbers)))
    residues = residue_records(sections, pointers, charges_e)
    host = max(residues, key=lambda row: row["atom_count"])
    guest_candidates = [row for row in residues if row is not host]
    guest = max(guest_candidates, key=lambda row: row["atom_count"]) if guest_candidates else None

    host_start = host["first_atom_1_based"] - 1
    host_end = host["last_atom_1_based"]
    host_slice = slice(host_start, host_end)
    host_atoms = set(range(host["first_atom_1_based"], host["last_atom_1_based"] + 1))

    physical_atom_properties = {
        "atom_names": atom_names[host_slice],
        "charges_e": charges_e[host_slice],
        "masses": masses[host_slice],
        "atomic_numbers": atomic_numbers[host_slice],
        "amber_atom_types": amber_types[host_slice],
    }
    compiled_indexing_properties = {
        "atom_type_indices": type_indices[host_slice],
        "tree_chain_classification": [
            str(v)
            for v in sections.get("TREE_CHAIN_CLASSIFICATION", [""] * natom)[host_slice]
        ],
        "join_array": [
            int(v) for v in sections.get("JOIN_ARRAY", [0] * natom)[host_slice]
        ],
        "irotat": [
            int(v) for v in sections.get("IROTAT", [0] * natom)[host_slice]
        ],
    }
    connectivity = connectivity_records(sections, host_atoms)
    exclusions = exclusion_records(sections, natom, host_atoms)
    lj = lj_records(sections, pointers, host_slice)

    physical_component_hashes = {
        "atom_properties": canonical_hash(physical_atom_properties),
        "bonds": canonical_hash(connectivity["bonds"]),
        "angles": canonical_hash(connectivity["angles"]),
        "dihedrals": canonical_hash(connectivity["dihedrals"]),
        "exclusions": canonical_hash(exclusions),
        "host_lj_pairs": canonical_hash(lj["physical"]),
    }
    indexing_component_hashes = {
        "compiled_indexing_properties": canonical_hash(compiled_indexing_properties),
        "indexed_host_lj_pairs": canonical_hash(lj["indexed"]),
    }
    physical_fingerprint = canonical_hash(physical_component_hashes)
    indexing_fingerprint = canonical_hash(indexing_component_hashes)

    residue_rows = []
    for row in residues:
        first = row["first_atom_1_based"] - 1
        end = row["last_atom_1_based"]
        residue_rows.append(
            {
                **row,
                "atom_names": atom_names[first:end],
                "amber_atom_types": amber_types[first:end],
                "atom_names_sha256": canonical_hash(atom_names[first:end]),
                "charges_sha256": canonical_hash(charges_e[first:end]),
            }
        )

    return {
        "path": path,
        "sha256": sha256_bytes(payload),
        "pointers": pointers,
        "natom": natom,
        "ntypes": pointers["NTYPES"],
        "nres": pointers["NRES"],
        "total_charge_e": sum(charges_e),
        "residues": residue_rows,
        "host": {
            "residue_index_1_based": host["index_1_based"],
            "label": host["label"],
            "first_atom_1_based": host["first_atom_1_based"],
            "last_atom_1_based": host["last_atom_1_based"],
            "atom_count": host["atom_count"],
            "charge_e": host["charge_e"],
            "physical_atom_properties": physical_atom_properties,
            "compiled_indexing_properties": compiled_indexing_properties,
            "host_lj": lj,
            "interaction_counts": {
                "bonds": len(connectivity["bonds"]),
                "angles": len(connectivity["angles"]),
                "dihedrals": len(connectivity["dihedrals"]),
                "excluded_atom_rows": len(exclusions),
                "lj_type_pairs": len(lj["physical"]),
            },
            "physical_component_hashes": physical_component_hashes,
            "indexing_component_hashes": indexing_component_hashes,
            "physical_parameter_fingerprint_sha256": physical_fingerprint,
            "compiled_indexing_fingerprint_sha256": indexing_fingerprint,
        },
        "guest": guest,
        "section_names": sorted(sections),
    }


def parse_restart(path: str, payload: bytes) -> dict[str, Any]:
    try:
        lines = payload.decode("utf-8").splitlines()
    except UnicodeDecodeError as exc:
        raise AuditError(f"{path} restart is not UTF-8 text") from exc
    if len(lines) < 2:
        raise AuditError(f"{path} restart is too short")
    header = lines[1].split()
    if not header:
        raise AuditError(f"{path} restart has no atom count")
    natom = int(header[0])
    values: list[float] = []
    for line in lines[2:]:
        values.extend(float(token.replace("D", "E").replace("d", "e")) for token in line.split())
    coordinate_count = 3 * natom
    if len(values) < coordinate_count:
        raise AuditError(
            f"{path} restart has {len(values)} values for {coordinate_count} coordinates"
        )
    remainder = len(values) - coordinate_count
    velocities_present = remainder >= 3 * natom
    after_velocity = remainder - 3 * natom if velocities_present else remainder
    if after_velocity not in {0, 3, 6}:
        raise AuditError(
            f"{path} restart has unsupported trailing value count {after_velocity}"
        )
    return {
        "path": path,
        "sha256": sha256_bytes(payload),
        "title": lines[0].strip(),
        "natom": natom,
        "time": float(header[1]) if len(header) > 1 else None,
        "coordinate_value_count": coordinate_count,
        "velocities_present": velocities_present,
        "box_value_count": after_velocity,
        "total_numeric_values_after_header": len(values),
    }


def max_abs_difference(left: list[float], right: list[float]) -> float | None:
    if len(left) != len(right):
        return None
    return max((abs(a - b) for a, b in zip(left, right)), default=0.0)


def compare_hosts(topologies: list[dict[str, Any]]) -> dict[str, Any]:
    if not topologies:
        raise AuditError("no topologies to compare")
    reference = topologies[0]
    reference_host = reference["host"]
    pairwise = []
    all_physical_equal = True
    all_indexing_equal = True
    for topology in topologies[1:]:
        host = topology["host"]
        ref_props = reference_host["physical_atom_properties"]
        props = host["physical_atom_properties"]
        ref_indexing = reference_host["compiled_indexing_properties"]
        indexing = host["compiled_indexing_properties"]
        charge_delta = max_abs_difference(ref_props["charges_e"], props["charges_e"])
        mass_delta = max_abs_difference(ref_props["masses"], props["masses"])
        physical_component_equal = {
            key: (
                reference_host["physical_component_hashes"][key]
                == host["physical_component_hashes"][key]
            )
            for key in reference_host["physical_component_hashes"]
        }
        indexing_component_equal = {
            key: (
                reference_host["indexing_component_hashes"][key]
                == host["indexing_component_hashes"][key]
            )
            for key in reference_host["indexing_component_hashes"]
        }
        physical_equal = (
            reference_host["atom_count"] == host["atom_count"]
            and reference_host["label"] == host["label"]
            and ref_props["atom_names"] == props["atom_names"]
            and ref_props["amber_atom_types"] == props["amber_atom_types"]
            and ref_props["atomic_numbers"] == props["atomic_numbers"]
            and charge_delta is not None
            and charge_delta <= FLOAT_TOLERANCE
            and mass_delta is not None
            and mass_delta <= FLOAT_TOLERANCE
            and all(physical_component_equal.values())
            and (
                reference_host["physical_parameter_fingerprint_sha256"]
                == host["physical_parameter_fingerprint_sha256"]
            )
        )
        indexing_equal = (
            ref_indexing == indexing
            and all(indexing_component_equal.values())
            and (
                reference_host["compiled_indexing_fingerprint_sha256"]
                == host["compiled_indexing_fingerprint_sha256"]
            )
        )
        row = {
            "reference_path": reference["path"],
            "comparison_path": topology["path"],
            "host_atom_count_equal": reference_host["atom_count"] == host["atom_count"],
            "host_label_equal": reference_host["label"] == host["label"],
            "atom_names_equal": ref_props["atom_names"] == props["atom_names"],
            "amber_atom_types_equal": (
                ref_props["amber_atom_types"] == props["amber_atom_types"]
            ),
            "atomic_numbers_equal": (
                ref_props["atomic_numbers"] == props["atomic_numbers"]
            ),
            "maximum_absolute_charge_difference_e": charge_delta,
            "maximum_absolute_mass_difference": mass_delta,
            "physical_component_hashes_equal": physical_component_equal,
            "physical_parameter_fingerprint_equal": (
                reference_host["physical_parameter_fingerprint_sha256"]
                == host["physical_parameter_fingerprint_sha256"]
            ),
            "compiled_indexing_properties_equal": ref_indexing == indexing,
            "indexing_component_hashes_equal": indexing_component_equal,
            "compiled_indexing_fingerprint_equal": (
                reference_host["compiled_indexing_fingerprint_sha256"]
                == host["compiled_indexing_fingerprint_sha256"]
            ),
            "physical_host_parameters_equal": physical_equal,
            "compiled_indexing_equal": indexing_equal,
        }
        all_physical_equal = all_physical_equal and physical_equal
        all_indexing_equal = all_indexing_equal and indexing_equal
        pairwise.append(row)
    return {
        "reference_path": reference["path"],
        "pairwise": pairwise,
        "physical_host_parameters_identical_across_all_topologies": (
            all_physical_equal
        ),
        "compiled_indexing_identical_across_all_topologies": all_indexing_equal,
        "host_atom_count": reference_host["atom_count"],
        "host_label": reference_host["label"],
        "host_charge_e": reference_host["charge_e"],
        "reference_physical_component_hashes": (
            reference_host["physical_component_hashes"]
        ),
        "reference_indexing_component_hashes": (
            reference_host["indexing_component_hashes"]
        ),
        "reference_physical_parameter_fingerprint_sha256": (
            reference_host["physical_parameter_fingerprint_sha256"]
        ),
        "reference_compiled_indexing_fingerprint_sha256": (
            reference_host["compiled_indexing_fingerprint_sha256"]
        ),
    }


def infer_system_name(path: str) -> str:
    lower = path.lower()
    for name in ("fructofuranose", "fructopyranose", "galactose", "glucose"):
        if name in lower:
            return name
    return PurePosixPath(path).parent.name.lower()


def audit(archive_path: Path, output_json: Path) -> dict[str, Any]:
    if archive_path.stat().st_size != ARCHIVE_SIZE_BYTES:
        raise AuditError(
            f"size mismatch: expected {ARCHIVE_SIZE_BYTES}, "
            f"got {archive_path.stat().st_size}"
        )
    md5 = digest_file(archive_path, "md5")
    if md5 != ARCHIVE_MD5:
        raise AuditError(f"MD5 mismatch: expected {ARCHIVE_MD5}, got {md5}")
    archive_sha256 = digest_file(archive_path, "sha256")

    topology_summaries: list[dict[str, Any]] = []
    restart_summaries: list[dict[str, Any]] = []
    protocol_members: list[dict[str, Any]] = []
    seen: set[str] = set()

    with zipfile.ZipFile(archive_path) as archive:
        bad = archive.testzip()
        if bad:
            raise AuditError(f"ZIP CRC failure: {bad}")
        for info in archive.infolist():
            safe_member(info.filename)
            if info.filename in seen:
                raise AuditError(f"duplicate member: {info.filename}")
            seen.add(info.filename)
            if info.flag_bits & 0x1:
                raise AuditError(f"encrypted member: {info.filename}")
            kind = member_kind(info)
            if kind in {"symlink", "special"}:
                raise AuditError(f"unsupported {kind} member: {info.filename}")
            if kind != "file":
                continue
            suffix = PurePosixPath(info.filename).suffix.lower()
            payload = archive.read(info)
            if suffix in {".parm7", ".prmtop"}:
                summary = topology_summary(info.filename, payload)
                summary["system_name"] = infer_system_name(info.filename)
                topology_summaries.append(summary)
            elif suffix in {".rst", ".rst7", ".inpcrd"}:
                summary = parse_restart(info.filename, payload)
                summary["system_name"] = infer_system_name(info.filename)
                restart_summaries.append(summary)
            elif suffix in {".in", ".mdin", ".txt"}:
                try:
                    text = payload.decode("utf-8")
                    line_count = len(text.splitlines())
                    head = text.splitlines()[:40]
                except UnicodeDecodeError:
                    line_count = None
                    head = []
                protocol_members.append(
                    {
                        "path": info.filename,
                        "system_name": infer_system_name(info.filename),
                        "sha256": sha256_bytes(payload),
                        "line_count": line_count,
                        "head": head,
                    }
                )

    topology_summaries.sort(key=lambda row: row["system_name"])
    restart_summaries.sort(key=lambda row: row["path"])
    protocol_members.sort(key=lambda row: row["path"])
    if len(topology_summaries) != 4:
        raise AuditError(f"expected 4 topologies, found {len(topology_summaries)}")

    host_comparison = compare_hosts(topology_summaries)
    topology_by_system = {row["system_name"]: row for row in topology_summaries}
    restart_matches = []
    for restart in restart_summaries:
        same_system = topology_by_system.get(restart["system_name"])
        compatible = [
            topology["path"]
            for topology in topology_summaries
            if topology["natom"] == restart["natom"]
        ]
        restart_matches.append(
            {
                **restart,
                "same_system_topology": same_system["path"] if same_system else None,
                "same_system_natom_matches": (
                    same_system is not None and same_system["natom"] == restart["natom"]
                ),
                "all_natom_compatible_topologies": compatible,
                "natom_alone_is_unique_match": len(compatible) == 1,
            }
        )

    guest_rows = []
    for topology in topology_summaries:
        guest = topology["guest"]
        guest_rows.append(
            {
                "system_name": topology["system_name"],
                "topology_path": topology["path"],
                "natom": topology["natom"],
                "total_charge_e": topology["total_charge_e"],
                "guest_residue": guest,
                "residue_labels": [row["label"] for row in topology["residues"]],
                "residue_atom_counts": [row["atom_count"] for row in topology["residues"]],
                "residue_charges_e": [row["charge_e"] for row in topology["residues"]],
            }
        )

    physical_host_identical = host_comparison[
        "physical_host_parameters_identical_across_all_topologies"
    ]
    report = {
        "schema": "dia-health-one.gluhut-compiled-system-audit.v1",
        "status": "compiled_system_audit_complete_source_parameter_provenance_missing",
        "scope": "public GluHUT SMD compiled AMBER systems",
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
        "topology_count": len(topology_summaries),
        "restart_count": len(restart_summaries),
        "protocol_document_count": len(protocol_members),
        "topologies": topology_summaries,
        "host_comparison": host_comparison,
        "guest_and_system_summary": guest_rows,
        "restart_matches": restart_matches,
        "protocol_documents": protocol_members,
        "decision": {
            "physical_host_parameters_identical_across_sugars": (
                physical_host_identical
            ),
            "compiled_indexing_identical_across_sugars": host_comparison[
                "compiled_indexing_identical_across_all_topologies"
            ],
            "compiled_public_host_can_be_used_as_common_benchmark_baseline": (
                physical_host_identical
            ),
            "source_resp_gaff_glycam_generation_reconstructable": False,
            "candidate_ranking_authorized": False,
            "pag_parameterization_authorized": False,
            "mcq1_parameter_transfer_authorized": False,
            "next_action": (
                "use the common compiled host only for public benchmark reproduction; "
                "request source RESP/GAFF/GLYCAM files and independently review the "
                "host graph, charge derivation, and guest parameters before PaG or MC-Q1"
            ),
        },
        "limitations": [
            "A compiled parm7 contains numerical charges and terms but not the full derivation, quantum-chemistry inputs, atom typing history, or original source parameter files.",
            "Exact equality across deposited sugar systems would show internal reuse, not chemical validation.",
            "All four systems can share NATOM, so atom count alone cannot identify a restart; directory and protocol provenance remain necessary.",
            "The public host is a truncated GluHUT benchmark system, not the complete NNC2215 insulin conjugate.",
            "No compiled parameter is assumed transferable to PaG or MC-Q1.",
        ],
    }
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return report


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser()
    result.add_argument("--archive", type=Path, required=True)
    result.add_argument(
        "--output-json",
        type=Path,
        default=Path("research/data/gluhut-compiled-system-audit.json"),
    )
    return result


def main() -> int:
    args = parser().parse_args()
    report = audit(args.archive, args.output_json)
    print(
        json.dumps(
            {
                "topology_count": report["topology_count"],
                "restart_count": report["restart_count"],
                "host_comparison": report["host_comparison"],
                "decision": report["decision"],
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
