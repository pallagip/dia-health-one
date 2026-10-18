#!/usr/bin/env python3
"""Audit the author-supplied ChimeraX session without treating it as chemistry truth.

The public .cxs file is an LZ4-compressed, concatenated MessagePack stream.
Decompress it first with:

    lz4 -d -f input.cxs output.session

Then run this script with the decompressed session and downloaded source file:

    python audit_nnc2215_chimerax_session.py \
      --session output.session --cxs input.cxs

The parser intentionally extracts only structure identity, atom/residue
membership, element counts, bond-edge counts, and the two cached IDATM type
strings. It does not infer bond orders, protonation, charges, or parameters.
Requires msgpack==1.2.3 (see chimerax-session-audit-requirements.txt).
"""

from __future__ import annotations

import argparse
import collections
import datetime as dt
import hashlib
import json
import math
import struct
import sys
from pathlib import Path
from typing import Any

import msgpack


REPO_ROOT = Path(__file__).resolve().parents[1]
FIGSHARE_DOI = "10.6084/m9.figshare.26526460.v1"
FIGSHARE_FILE_ID = 48308389
FIGSHARE_FILENAME = "NNC2215 in open and closed forms docked into insulin receptor.cxs"
FIGSHARE_URL = f"https://ndownloader.figshare.com/files/{FIGSHARE_FILE_ID}"
FIGSHARE_MD5 = "33143f4fad4e80bf4accd09718339aa7"
CHIMERAX_COMMIT = "b18ee9afe9a1536d0044139dcdc4c0ea539d2f7f"
CHIMERAX_SOURCE_BASE = (
    "https://github.com/RBVI/ChimeraX/blob/" + CHIMERAX_COMMIT
    + "/src/bundles/atomic_lib/atomic_cpp/atomstruct_cpp/"
)
ELEMENTS = {
    1: "H",
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
TARGET_RESIDUES = {"PHX", "LYX", "MC"}
RESIDUE_NAME_ALIASES = {"BGLC": "BGL"}
MODEL_NAMES = {
    "2215_inactive_correctaminoacids.pdb": "inactive_closed",
    "2215_active-alt_glucose_rec_correctaminopacids.pdb": "active_open",
}


def _unpack_nested(data: bytes) -> Any:
    return msgpack.unpackb(
        data,
        raw=False,
        strict_map_key=False,
        ext_hook=lambda code, payload: msgpack.ExtType(code, payload),
    )


def _array_values(value: Any) -> list[int | float]:
    if not isinstance(value, msgpack.ExtType) or value.code != 1:
        raise ValueError("Expected ChimeraX numeric-array extension type 1")
    packed = _unpack_nested(value.data)
    dtype = packed.get(b"dtype")
    shape = packed.get(b"shape")
    payload = packed.get(b"data")
    if not isinstance(dtype, str) or not isinstance(shape, list) or not isinstance(payload, bytes):
        raise ValueError("Unrecognized ChimeraX numeric-array payload")
    formats = {"<i4": "i", "<f4": "f", "<f8": "d"}
    if dtype not in formats:
        raise ValueError(f"Unsupported array dtype: {dtype}")
    count = math.prod(shape)
    if len(payload) != count * struct.calcsize(formats[dtype]):
        raise ValueError("Numeric-array byte length does not match its shape")
    return list(struct.unpack("<" + str(count) + formats[dtype], payload))


def _source_file_digest(path: Path) -> dict[str, Any]:
    digest = hashlib.md5()
    sha256 = hashlib.sha256()
    size = 0
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            size += len(chunk)
            digest.update(chunk)
            sha256.update(chunk)
    return {"bytes": size, "md5": digest.hexdigest(), "sha256": sha256.hexdigest()}


def _read_pdb_identities(path: Path) -> list[tuple[str, str, str, int, str]]:
    identities = []
    for line in path.read_text(encoding="ascii", errors="replace").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        try:
            serial = int(line[6:11])
            resseq = int(line[22:26])
        except ValueError:
            continue
        element = line[76:78].strip().title()
        if not element:
            name = line[12:16].strip()
            element = "".join(c for c in name if c.isalpha())[:1].title()
        identities.append(
            (
                line[17:20].strip(),
                line[12:16].strip(),
                line[21:22].strip(),
                resseq,
                element,
            )
        )
    return identities


def _parse_structure(model_state: dict[str, Any], state: dict[str, Any], model_name: str) -> dict[str, Any]:
    version = state.get("version")
    if version != 19:
        raise ValueError(f"Expected ChimeraX AtomicStructure session version 19; found {version}")

    ints = [_array_values(item) for item in state["ints"]]
    atom_misc = state["misc"][1]
    atom_names = atom_misc[0]
    atom_types = atom_misc[1:]
    residue_names, chain_ids, _mmcif_chain_ids = state["misc"][5]
    atom_count = len(atom_names)
    atomic_numbers = ints[1][:atom_count]
    if len(atomic_numbers) != atom_count or len(atom_types) != atom_count:
        raise ValueError("Atom names, element numbers, and type metadata have inconsistent lengths")

    # Pinned ChimeraX session format v19: two external residue integers
    # (number/insertion code), two RGBA records of four integers each,
    # nine Residue session integers, then the residue's atom indices.
    rgba_ints_per_color = 4
    residue_base_ints = 9
    residue_header_size = 2 + (2 * rgba_ints_per_color) + residue_base_ints
    residue_ints = ints[5]
    cursor = 0
    atom_residues: list[tuple[str, str, int, int] | None] = [None] * atom_count
    residue_inventory: dict[str, collections.Counter[str]] = collections.defaultdict(collections.Counter)
    residue_type_counts: dict[str, collections.Counter[str]] = collections.defaultdict(collections.Counter)
    explicit_types: collections.Counter[str] = collections.Counter()
    computed_types: collections.Counter[str] = collections.Counter()
    target_explicit_types: dict[str, collections.Counter[str]] = collections.defaultdict(collections.Counter)
    target_computed_types: dict[str, collections.Counter[str]] = collections.defaultdict(collections.Counter)
    residue_atom_counts = {}

    for residue_index, (residue_name, chain_id) in enumerate(zip(residue_names, chain_ids)):
        number = int(residue_ints[cursor])
        insertion_code = int(residue_ints[cursor + 1])
        num_atoms_offset = cursor + residue_header_size - 1
        num_residue_atoms = int(residue_ints[num_atoms_offset])
        atom_start = cursor + residue_header_size
        atom_indices = residue_ints[atom_start : atom_start + num_residue_atoms]
        if len(atom_indices) != num_residue_atoms:
            raise ValueError(f"Residue {residue_index} atom membership is truncated")
        residue_atom_counts[f"{chain_id}:{residue_name}{number}"] = num_residue_atoms
        residue = (residue_name, chain_id, number, insertion_code)
        for atom_index in atom_indices:
            if atom_index < 0 or atom_index >= atom_count or atom_residues[atom_index] is not None:
                raise ValueError(f"Invalid or repeated atom index {atom_index} in residue table")
            atom_residues[atom_index] = residue
            element = ELEMENTS.get(int(atomic_numbers[atom_index]), f"Z{int(atomic_numbers[atom_index])}")
            residue_inventory[residue_name][element] += 1
            computed, explicit = atom_types[atom_index]
            if computed:
                computed_types[str(computed)] += 1
            if explicit:
                explicit_types[str(explicit)] += 1
            if residue_name in TARGET_RESIDUES:
                if computed:
                    target_computed_types[residue_name][str(computed)] += 1
                if explicit:
                    target_explicit_types[residue_name][str(explicit)] += 1
        cursor = atom_start + num_residue_atoms

    if cursor != len(residue_ints):
        raise ValueError(f"Residue table consumed {cursor} of {len(residue_ints)} integers")
    if any(residue is None for residue in atom_residues):
        raise ValueError("Some atoms do not belong to a serialized residue")

    bond_ints = ints[2]
    bond_count = int(bond_ints[0])
    bond_pairs = bond_ints[1 : 1 + 2 * bond_count]
    if len(bond_pairs) != 2 * bond_count:
        raise ValueError("Bond endpoint table is truncated")
    for endpoint in bond_pairs:
        if endpoint < 0 or endpoint >= atom_count:
            raise ValueError(f"Bond endpoint {endpoint} is outside the atom table")

    identities = []
    for index, (atom_name, residue, atomic_number) in enumerate(zip(atom_names, atom_residues, atomic_numbers)):
        residue_name, chain_id, number, _insertion_code = residue  # type: ignore[misc]
        identities.append(
            (residue_name, str(atom_name).strip(), str(chain_id).strip(), number,
             ELEMENTS.get(int(atomic_number), f"Z{int(atomic_number)}"))
        )

    target_heavy = collections.Counter()
    for residue_name in TARGET_RESIDUES:
        for element, count in residue_inventory[residue_name].items():
            if element != "H":
                target_heavy[element] += count

    return {
        "state": model_name,
        "model_name": model_state.get("name"),
        "source_title_records": state["misc"][0][2].get("TITLE", []),
        "chimerax_structure_session_version": version,
        "atom_count": atom_count,
        "residue_count": len(residue_names),
        "bond_edge_count": bond_count,
        "atom_identity_sequence": identities,
        "element_counts_by_target_residue": {
            residue: dict(sorted(residue_inventory[residue].items()))
            for residue in sorted(TARGET_RESIDUES)
        },
        "target_attachment_heavy_atom_counts": dict(sorted(target_heavy.items())),
        "computed_idatm_type_counts_by_target_residue": {
            residue: dict(sorted(target_computed_types[residue].items()))
            for residue in sorted(TARGET_RESIDUES)
        },
        "explicit_idatm_type_counts_by_target_residue": {
            residue: dict(sorted(target_explicit_types[residue].items()))
            for residue in sorted(TARGET_RESIDUES)
        },
        "explicit_idatm_type_counts_all_atoms": dict(sorted(explicit_types.items())),
        "computed_idatm_type_counts_all_atoms": dict(sorted(computed_types.items())),
        "pdb_link_records": [
            line.strip()
            for line in state["misc"][0][2].get("LINK", [])
        ],
    }


def _extract_models(session_path: Path) -> dict[str, dict[str, Any]]:
    models = {}
    with session_path.open("rb") as handle:
        header = handle.readline()
        if not header.startswith(b"# ChimeraX Session version "):
            raise ValueError("Decompressed file lacks a ChimeraX session header")
        unpacker = msgpack.Unpacker(
            handle,
            raw=False,
            strict_map_key=False,
            ext_hook=lambda code, payload: msgpack.ExtType(code, payload),
            max_buffer_size=200_000_000,
        )
        for item in unpacker:
            if not isinstance(item, dict) or "structure state" not in item:
                continue
            outer = item["structure state"]
            if not isinstance(outer, dict) or "structure state" not in outer:
                continue
            model_state = outer.get("model state", {})
            model_name = model_state.get("name")
            canonical_state = MODEL_NAMES.get(model_name)
            if canonical_state is None:
                continue
            raw_state = outer["structure state"]
            if not isinstance(raw_state, msgpack.ExtType) or raw_state.code != 11:
                raise ValueError(f"Unexpected AtomicStructure state extension for {model_name}")
            state = _unpack_nested(raw_state.data)
            models[canonical_state] = _parse_structure(model_state, state, canonical_state)
    expected = {"active_open", "inactive_closed"}
    if set(models) != expected:
        raise ValueError(f"Expected models {sorted(expected)}; found {sorted(models)}")
    return models


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--session", type=Path, required=True, help="LZ4-decompressed .cxs stream")
    parser.add_argument("--cxs", type=Path, help="Original compressed .cxs file for checksum verification")
    parser.add_argument("--output", type=Path, default=REPO_ROOT / "research/data/nnc2215-chimerax-session-audit.json")
    args = parser.parse_args()

    models = _extract_models(args.session)
    for state, filename in (
        ("active_open", "nnc2215-open-model.pdb"),
        ("inactive_closed", "nnc2215-closed-model.pdb"),
    ):
        pdb_path = REPO_ROOT / "research/structures" / filename
        pdb_identities = _read_pdb_identities(pdb_path)
        cxs_identities = models[state].pop("atom_identity_sequence")
        normalized_cxs_identities = [
            (RESIDUE_NAME_ALIASES.get(residue, residue), atom, chain, number, element)
            for residue, atom, chain, number, element in cxs_identities
        ]
        models[state]["compared_pdb"] = filename
        models[state]["pdb_coordinate_atom_count"] = len(pdb_identities)
        models[state]["pdb_atom_identity_sequence_match_before_alias"] = cxs_identities == pdb_identities
        models[state]["pdb_atom_identity_name_aliases_applied"] = RESIDUE_NAME_ALIASES
        models[state]["pdb_atom_identity_sequence_match"] = normalized_cxs_identities == pdb_identities
        if normalized_cxs_identities != pdb_identities:
            first_differences = []
            for index, (left, right) in enumerate(zip(normalized_cxs_identities, pdb_identities)):
                if left != right:
                    first_differences.append({"index": index, "session": left, "pdb": right})
                    if len(first_differences) == 10:
                        break
            models[state]["first_atom_identity_differences"] = first_differences

    cxs_digest = _source_file_digest(args.cxs) if args.cxs else None
    source_hash_verified = None
    if cxs_digest is not None:
        source_hash_verified = cxs_digest["md5"] == FIGSHARE_MD5

    result = {
        "schema": "dia-health-one.nnc2215-chimerax-session-audit.v1",
        "run": {
            "run_id": "DHO-NNC2215-CXS-" + dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ"),
            "created_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
            "script": "research/audit_nnc2215_chimerax_session.py",
            "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "runtime_python": sys.version.split()[0],
            "msgpack_version": msgpack.__version__,
            "nvidia_gpu_used": False,
        },
        "source": {
            "dataset_doi": FIGSHARE_DOI,
            "dataset_url": "https://springernature.figshare.com/articles/dataset/3D_models_of_NNC2215/26526460",
            "dataset_license": "CC BY 4.0",
            "session_filename": FIGSHARE_FILENAME,
            "session_file_id": FIGSHARE_FILE_ID,
            "session_download_url": FIGSHARE_URL,
            "published_file_bytes": 79145092,
            "published_file_md5": FIGSHARE_MD5,
            "local_cxs_digest": cxs_digest,
            "local_cxs_md5_matches_published": source_hash_verified,
            "decompressed_session_bytes": args.session.stat().st_size,
            "decompressed_session_sha256": hashlib.sha256(args.session.read_bytes()).hexdigest(),
            "chimerax_source_commit": CHIMERAX_COMMIT,
            "chimerax_structure_source": CHIMERAX_SOURCE_BASE + "Structure.cpp",
            "chimerax_residue_source": CHIMERAX_SOURCE_BASE + "Residue.h",
            "chimerax_atom_source": CHIMERAX_SOURCE_BASE + "Atom.cpp",
            "chimerax_bond_source": CHIMERAX_SOURCE_BASE + "Bond.h",
        },
        "serialization_findings": {
            "bond_edges_saved_as_atom_index_pairs": True,
            "bond_order_field_saved": False,
            "basis": "At the pinned ChimeraX source commit, Structure::session_info serializes the number of bonds followed by endpoint atom-index pairs and then Bond session state; Bond has no bond-order property. Atom session state stores computed and explicit IDATM type strings separately.",
            "source_files": [
                CHIMERAX_SOURCE_BASE + "Structure.cpp",
                CHIMERAX_SOURCE_BASE + "Bond.h",
                CHIMERAX_SOURCE_BASE + "Atom.cpp",
            ],
        },
        "models": models,
        "formula_comparison": {
            "published_expected_attachment_heavy_atoms": {"C": 72, "N": 15, "O": 22},
            "model_target_residue_heavy_atoms": {state: model["target_attachment_heavy_atom_counts"] for state, model in models.items()},
            "interpretation": "The CXS snapshot preserves the same model atom inventory. Its computed IDATM types are metadata derived by ChimeraX; they are not an atom-mapped source chemical definition. The two-carbon excess in LYX+MC versus the published formula-derived B29 component remains unresolved.",
        },
        "readiness": {
            "model_coordinates_match_pdb_atom_identity_order": all(m["pdb_atom_identity_sequence_match"] for m in models.values()),
            "explicit_idatm_types_present_for_target_residues": any(
                any(bool(type_counts) for type_counts in m["explicit_idatm_type_counts_by_target_residue"].values())
                for m in models.values()
            ),
            "bond_orders_authoritatively_resolved": False,
            "two_carbon_discrepancy_resolved": False,
            "force_field_parameters_validated": False,
            "production_md_ready": False,
            "target_binding_md_started": False,
            "clinical_use": False,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(args.output), "models": {k: {"atoms": v["atom_count"], "bonds": v["bond_edge_count"], "pdb_identity_match": v["pdb_atom_identity_sequence_match"], "computed_type_counts": v["computed_idatm_type_counts_by_target_residue"], "explicit_target_type_counts": v["explicit_idatm_type_counts_by_target_residue"]} for k, v in models.items()}, "cxs_md5_matches_published": source_hash_verified, "production_md_ready": result["readiness"]["production_md_ready"]}, indent=2))


if __name__ == "__main__":
    main()
