#!/usr/bin/env python3
"""Run a bounded GAFF2/AM1-BCC pilot on the B1-linked Phe peptide context.

The local neutral fragment contains the B1 glucoside, PheB1, ValB2, and an
AsnB3-NH2 truncation cap. It is a toolchain feasibility result, not a
validated insulin force field or simulation-ready system.
"""

from __future__ import annotations

import argparse
from collections import Counter
from decimal import Decimal
import json
import math
import os
from pathlib import Path
import shutil
import tempfile
from typing import Any

from rdkit import Chem, rdBase
from rdkit.Chem import rdMolDescriptors

import audit_nnc2215_full_conjugate_graph as full_graph_audit
from audit_nnc2215_ambertools_parameterization_pilot import (
    copy_normalized_tool_output,
    frcmod_inventory,
    read_lock_package,
    read_mol2,
    run_checked,
    sha256,
    write_json,
)


ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
B1_GRAPH = ROOT / "structures/nnc2215-b1-glucoside-attachment.smi"
B1_MAP = DATA / "nnc2215-b1-glucoside-atom-map.json"
B1_TOPOLOGY = DATA / "nnc2215-b1-topology-audit.json"
FULL_GRAPH = DATA / "nnc2215-full-conjugate-graph.json"
FORMULA_AUDIT = DATA / "nnc2215-model-formula-audit.json"
B29_GRAPH = ROOT / "structures/nnc2215-b29-linker-roof-product.smi"
B29_MAP = DATA / "nnc2215-b29-linker-roof-product-map.json"
RECONCILIATION = DATA / "nnc2215-carbon-reconciliation-candidate.json"
OPEN_PDB = ROOT / "structures/nnc2215-open-source-corrected-candidate.pdb"
LOCK_PATH = ROOT / "pag-chemistry-environment-macos-arm64.lock.txt"
FULL_GRAPH_SCRIPT = ROOT / "audit_nnc2215_full_conjugate_graph.py"
B29_PILOT_HELPERS = ROOT / "audit_nnc2215_ambertools_parameterization_pilot.py"
DEFAULT_OUTPUT_DIR = DATA / "nnc2215-b1-gaff2-boundary-pilot"
DEFAULT_JSON = DATA / "nnc2215-b1-gaff2-boundary-pilot.json"

EXPECTED_HEAVY_COUNT = 37
EXPECTED_FORMULA = "C24H37N3O10"
EXPECTED_FORMAL_CHARGE = 0
B1_SOURCE_MAP_IDS = set(range(2, 20))
EXPECTED_PROTEIN_IDENTITIES = {
    ("B", 1, "PHE", name)
    for name in ("N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ")
} | {
    ("B", 2, "VAL", name)
    for name in ("N", "CA", "C", "O", "CB", "CG1", "CG2")
} | {("B", 3, "ASN", "N")}


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def read_pdb_heavy_atoms(path: Path) -> tuple[dict[tuple[str, int, str, str, str], list[dict[str, Any]]], dict[int, dict[str, Any]]]:
    by_identity: dict[tuple[str, int, str, str, str], list[dict[str, Any]]] = {}
    by_serial: dict[int, dict[str, Any]] = {}
    for line_number, line in enumerate(path.read_text(encoding="ascii").splitlines(), start=1):
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        serial = int(line[6:11])
        element = line[76:78].strip().upper()
        if not element:
            raise ValueError(f"PDB atom at line {line_number} has a blank element field")
        if element in {"H", "D"}:
            continue
        atom = {
            "serial": serial,
            "record_type": line[0:6].strip(),
            "atom_name": line[12:16].strip(),
            "altloc": line[16].strip(),
            "residue_name": line[17:20].strip(),
            "chain": line[21].strip(),
            "residue_number": int(line[22:26]),
            "insertion_code": line[26].strip(),
            "element": element,
            "coordinates_angstrom": [float(line[30:38]), float(line[38:46]), float(line[46:54])],
        }
        if serial in by_serial:
            raise ValueError(f"Duplicate PDB atom serial {serial}")
        by_serial[serial] = atom
        key = (atom["chain"], atom["residue_number"], atom["residue_name"], atom["atom_name"], atom["altloc"])
        by_identity.setdefault(key, []).append(atom)
    return by_identity, by_serial


def one_pdb_atom(
    by_identity: dict[tuple[str, int, str, str, str], list[dict[str, Any]]],
    chain: str,
    residue_number: int,
    residue_name: str,
    atom_name: str,
) -> dict[str, Any]:
    matches = by_identity.get((chain, residue_number, residue_name, atom_name, ""), [])
    if len(matches) != 1:
        raise ValueError(
            f"Expected one PDB identity {chain}:{residue_number}:{residue_name}:{atom_name}; found {len(matches)}"
        )
    return matches[0]


def selected_parent_atom(atom: Chem.Atom) -> tuple[str, Any] | None:
    if atom.HasProp("_NNCOrigin") and atom.GetProp("_NNCOrigin") == "mapped_source_fragment":
        if atom.GetProp("_NNCComponent") == "B1":
            map_id = int(atom.GetProp("_NNCSourceMap"))
            if map_id in B1_SOURCE_MAP_IDS:
                return "source_map", map_id
        return None
    required = ("_NNCChain", "_NNCResidueNumber", "_NNCResidueName", "_NNCAtomName")
    if not all(atom.HasProp(name) for name in required):
        return None
    identity = (
        atom.GetProp("_NNCChain"),
        int(atom.GetProp("_NNCResidueNumber")),
        atom.GetProp("_NNCResidueName"),
        atom.GetProp("_NNCAtomName"),
    )
    if identity in EXPECTED_PROTEIN_IDENTITIES:
        return "protein_identity", identity
    return None


def extract_local_fragment(parent: Chem.Mol) -> Chem.Mol:
    keep_indices: set[int] = set()
    source_maps: set[int] = set()
    protein_identities: set[tuple[str, int, str, str]] = set()
    for atom in parent.GetAtoms():
        selected = selected_parent_atom(atom)
        if selected is None:
            continue
        keep_indices.add(atom.GetIdx())
        kind, value = selected
        if kind == "source_map":
            source_maps.add(int(value))
        else:
            protein_identities.add(value)
    if source_maps != B1_SOURCE_MAP_IDS:
        raise ValueError(f"Selected B1 source maps differ from 2..19: {sorted(source_maps)}")
    if protein_identities != EXPECTED_PROTEIN_IDENTITIES:
        missing = sorted(EXPECTED_PROTEIN_IDENTITIES - protein_identities)
        extra = sorted(protein_identities - EXPECTED_PROTEIN_IDENTITIES)
        raise ValueError(f"Protein context identities differ; missing={missing}, extra={extra}")
    if len(keep_indices) != EXPECTED_HEAVY_COUNT:
        raise ValueError(f"Selected {len(keep_indices)} heavy atoms, expected {EXPECTED_HEAVY_COUNT}")

    editable = Chem.RWMol(parent)
    for atom_index in reversed(range(parent.GetNumAtoms())):
        if atom_index not in keep_indices:
            editable.RemoveAtom(atom_index)
    fragment = editable.GetMol()
    Chem.SanitizeMol(fragment)
    for atom in fragment.GetAtoms():
        atom.UpdatePropertyCache(strict=True)
    Chem.AssignStereochemistry(fragment, cleanIt=True, force=True)

    formula = rdMolDescriptors.CalcMolFormula(fragment)
    if formula != EXPECTED_FORMULA:
        raise ValueError(f"Extracted fragment formula {formula}; expected {EXPECTED_FORMULA}")
    if fragment.GetNumHeavyAtoms() != EXPECTED_HEAVY_COUNT:
        raise ValueError("Extracted fragment heavy-atom count changed after sanitization")
    if Chem.GetFormalCharge(fragment) != EXPECTED_FORMAL_CHARGE:
        raise ValueError(f"Extracted fragment formal charge {Chem.GetFormalCharge(fragment)}; expected 0")
    if len(Chem.GetMolFrags(fragment, asMols=False)) != 1:
        raise ValueError("Extracted B1 local fragment is not connected")

    map_to_atom = {
        int(atom.GetProp("_NNCSourceMap")): atom.GetIdx()
        for atom in fragment.GetAtoms()
        if atom.HasProp("_NNCSourceMap") and atom.GetProp("_NNCComponent") == "B1"
    }
    phe_n = next(
        atom.GetIdx() for atom in fragment.GetAtoms()
        if atom.HasProp("_NNCChain")
        and atom.GetProp("_NNCChain") == "B"
        and atom.GetProp("_NNCResidueNumber") == "1"
        and atom.GetProp("_NNCAtomName") == "N"
    )
    val_c = next(
        atom.GetIdx() for atom in fragment.GetAtoms()
        if atom.HasProp("_NNCChain")
        and atom.GetProp("_NNCChain") == "B"
        and atom.GetProp("_NNCResidueNumber") == "2"
        and atom.GetProp("_NNCAtomName") == "C"
    )
    cap_n = next(
        atom.GetIdx() for atom in fragment.GetAtoms()
        if atom.HasProp("_NNCChain")
        and atom.GetProp("_NNCChain") == "B"
        and atom.GetProp("_NNCResidueNumber") == "3"
        and atom.GetProp("_NNCAtomName") == "N"
    )
    if fragment.GetBondBetweenAtoms(phe_n, map_to_atom[2]) is None:
        raise ValueError("PheB1 N–B1 map 2 source attachment edge is missing")
    if fragment.GetBondBetweenAtoms(val_c, cap_n) is None:
        raise ValueError("ValB2 C–AsnB3 N cap edge is missing")
    if fragment.GetAtomWithIdx(phe_n).GetTotalNumHs() != 1:
        raise ValueError("PheB1 conjugated nitrogen must have one graph-derived hydrogen")
    if fragment.GetAtomWithIdx(cap_n).GetTotalNumHs() != 2:
        raise ValueError("AsnB3 N cap must sanitize as neutral NH2 after truncation")

    center_records = Chem.FindMolChiralCenters(
        fragment, includeUnassigned=True, includeCIP=True, useLegacyImplementation=False
    )
    cip_by_idx = dict(center_records)
    expected_centers: dict[str, str] = {}
    for atom in fragment.GetAtoms():
        if atom.HasProp("_NNCComponent") and atom.GetProp("_NNCComponent") == "B1":
            map_id = int(atom.GetProp("_NNCSourceMap"))
            expected = {9: "R", 11: "R", 13: "S", 15: "S", 17: "R"}.get(map_id)
            if expected:
                expected_centers[f"B1 map {map_id}"] = expected
        if atom.HasProp("_NNCChain") and atom.GetProp("_NNCChain") == "B":
            identity = (atom.GetProp("_NNCResidueNumber"), atom.GetProp("_NNCAtomName"))
            if identity == ("1", "CA"):
                expected_centers["PheB1 CA"] = "S"
            elif identity == ("2", "CA"):
                expected_centers["ValB2 CA"] = "S"
    actual_centers: dict[str, str] = {}
    for label, expected in expected_centers.items():
        if label.startswith("B1 map "):
            map_id = int(label.removeprefix("B1 map "))
            atom_index = map_to_atom[map_id]
        else:
            residue_number, atom_name = ("1", "CA") if label == "PheB1 CA" else ("2", "CA")
            atom_index = next(
                atom.GetIdx() for atom in fragment.GetAtoms()
                if atom.HasProp("_NNCChain") and atom.GetProp("_NNCChain") == "B"
                and atom.GetProp("_NNCResidueNumber") == residue_number
                and atom.GetProp("_NNCAtomName") == atom_name
            )
        observed = cip_by_idx.get(atom_index)
        if observed != expected:
            raise ValueError(f"Stereochemistry at {label} is {observed!r}; expected CIP {expected}")
        actual_centers[label] = observed
    if len(expected_centers) != 7 or len(center_records) != 7:
        raise ValueError(f"Expected exactly seven retained stereocentres; found {len(center_records)}")
    return fragment


def map_coordinates(
    fragment: Chem.Mol,
    mapping_entry: dict[str, Any],
    pdb_by_identity: dict[tuple[str, int, str, str, str], list[dict[str, Any]]],
) -> tuple[Chem.Mol, list[dict[str, Any]]]:
    map_records = mapping_entry["map_id_to_target"]
    if set(map_records) != {str(map_id) for map_id in range(1, 20)}:
        raise ValueError("B1 topology mapping must cover map IDs 1..19 exactly")
    target_serials = [int(record["target_atom"]["serial"]) for record in map_records.values()]
    if len(set(target_serials)) != 19:
        raise ValueError("B1 mapping must use 19 distinct target serials")
    map1_target = map_records["1"]["target_atom"]
    if (map1_target["chain"], int(map1_target["residue_number"]), map1_target["residue"], map1_target["name"], map1_target["element"]) != (
        "B", 1, "PHE", "N", "N"
    ):
        raise ValueError("B1 map 1 must resolve to PheB1 N")

    conformer = Chem.Conformer(fragment.GetNumAtoms())
    conformer.Set3D(True)
    atom_coordinates: list[dict[str, Any]] = []
    used_serials: set[int] = set()
    for atom in fragment.GetAtoms():
        selected = selected_parent_atom(atom)
        if selected is None:
            raise ValueError(f"Unexpected atom in extracted fragment at index {atom.GetIdx()}")
        kind, value = selected
        if kind == "source_map":
            map_id = int(value)
            target = map_records[str(map_id)]["target_atom"]
            key = (
                target["chain"], int(target["residue_number"]), target["residue"], target["name"], target.get("altloc", "")
            )
            matches = pdb_by_identity.get(key, [])
            if len(matches) != 1:
                raise ValueError(f"B1 map {map_id} target identity is not unique in PDB: {key}; found {len(matches)}")
            pdb_atom = matches[0]
            if pdb_atom["serial"] != int(target["serial"]):
                raise ValueError(f"B1 map {map_id} serial differs from the topology audit")
        else:
            chain, residue_number, residue_name, atom_name = value
            pdb_atom = one_pdb_atom(pdb_by_identity, chain, residue_number, residue_name, atom_name)
        if pdb_atom["element"] != atom.GetSymbol().upper():
            raise ValueError(
                f"Element mismatch for selected atom {atom.GetIdx()}: graph={atom.GetSymbol()}, PDB={pdb_atom['element']}"
            )
        if pdb_atom["serial"] in used_serials:
            raise ValueError(f"Two fragment atoms map to PDB serial {pdb_atom['serial']}")
        used_serials.add(pdb_atom["serial"])
        conformer.SetAtomPosition(atom.GetIdx(), pdb_atom["coordinates_angstrom"])
        atom.SetProp("_NNC_PDBSerial", str(pdb_atom["serial"]))
        atom.SetProp("_NNC_PDBIdentity", f"{pdb_atom['chain']}:{pdb_atom['residue_number']}:{pdb_atom['residue_name']}:{pdb_atom['atom_name']}")
        atom_coordinates.append({
            "atom_index_zero_based": atom.GetIdx(),
            "element": atom.GetSymbol(),
            "source_map_id": int(value) if kind == "source_map" else None,
            "pdb_serial": pdb_atom["serial"],
            "pdb_identity": {
                "chain": pdb_atom["chain"],
                "residue_number": pdb_atom["residue_number"],
                "residue_name": pdb_atom["residue_name"],
                "atom_name": pdb_atom["atom_name"],
                "altloc": pdb_atom["altloc"],
            },
            "coordinates_angstrom": pdb_atom["coordinates_angstrom"],
        })
    if len(used_serials) != EXPECTED_HEAVY_COUNT:
        raise ValueError("Heavy-atom coordinate assignment did not cover all selected atoms")
    fragment.AddConformer(conformer, assignId=True)
    hydrogenated = Chem.AddHs(fragment, addCoords=True)
    if hydrogenated.GetNumAtoms() - hydrogenated.GetNumHeavyAtoms() != 37:
        raise ValueError("Explicit hydrogen count differs from the expected C24H37N3O10 formula")
    return hydrogenated, atom_coordinates


def mol2_atom_names(path: Path) -> list[str]:
    names: list[str] = []
    in_atoms = False
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("@<TRIPOS>"):
            in_atoms = line.strip() == "@<TRIPOS>ATOM"
            continue
        if in_atoms and line.strip() and not line.startswith("#"):
            fields = line.split()
            names.append(fields[1])
    return names


def validate_mol2_correspondence(
    molecule: Chem.Mol,
    output_atoms: list[dict[str, Any]],
    output_bonds: list[tuple[int, int, str]],
    coordinate_tolerance_angstrom: float = 0.002,
) -> dict[str, Any]:
    expected_count = molecule.GetNumAtoms()
    observed_ids = [int(atom["id"]) for atom in output_atoms]
    if observed_ids != list(range(1, expected_count + 1)):
        raise ValueError("MOL2 atom IDs must be unique, sequential, and in prepared-SDF atom order")
    if len(output_atoms) != expected_count or len(output_bonds) != molecule.GetNumBonds():
        raise ValueError("MOL2 atom/bond count differs from the prepared SDF")
    if not molecule.GetNumConformers():
        raise ValueError("Prepared molecule has no coordinates for MOL2 correspondence validation")

    conformer = molecule.GetConformer()
    coordinate_deltas = [
        math.dist(
            tuple(conformer.GetAtomPosition(index)),
            tuple(output_atoms[index]["coordinates_angstrom"]),
        )
        for index in range(expected_count)
    ]
    maximum_coordinate_delta = max(coordinate_deltas, default=0.0)
    if maximum_coordinate_delta > coordinate_tolerance_angstrom:
        raise ValueError(
            "MOL2 coordinates do not match prepared-SDF atom order: "
            f"maximum delta {maximum_coordinate_delta:.6f} A exceeds "
            f"{coordinate_tolerance_angstrom:.6f} A"
        )

    bond_type_labels = {
        Chem.BondType.SINGLE: "1",
        Chem.BondType.DOUBLE: "2",
        Chem.BondType.TRIPLE: "3",
        Chem.BondType.AROMATIC: "ar",
    }
    expected_bonds: dict[tuple[int, int], str] = {}
    for bond in molecule.GetBonds():
        bond_type = "ar" if bond.GetIsAromatic() else bond_type_labels.get(bond.GetBondType())
        if bond_type is None:
            raise ValueError(f"Unsupported prepared-SDF bond type: {bond.GetBondType()}")
        edge = tuple(sorted((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())))
        expected_bonds[edge] = bond_type

    observed_bond_map: dict[tuple[int, int], str] = {}
    for begin, end, bond_type in output_bonds:
        if begin < 0 or end < 0 or begin >= expected_count or end >= expected_count:
            raise ValueError(f"MOL2 bond references an invalid atom ID: {begin + 1}-{end + 1}")
        edge = tuple(sorted((begin, end)))
        if edge in observed_bond_map:
            raise ValueError(f"MOL2 contains a duplicate bond for atom IDs {edge[0] + 1}-{edge[1] + 1}")
        normalized_type = "1" if bond_type.lower() == "am" else bond_type.lower()
        observed_bond_map[edge] = normalized_type
    if observed_bond_map != expected_bonds:
        raise ValueError("MOL2 bond graph or bond orders differ from the prepared SDF")

    expected_elements = [atom.GetSymbol().upper() for atom in molecule.GetAtoms()]
    observed_elements = [str(atom["element"]).upper() for atom in output_atoms]
    if observed_elements != expected_elements:
        raise ValueError("MOL2 atom element order differs from the prepared SDF")
    return {
        "atom_ids_sequential_in_prepared_sdf_order": True,
        "atom_coordinates_match_prepared_sdf": True,
        "coordinate_tolerance_angstrom": coordinate_tolerance_angstrom,
        "maximum_coordinate_delta_angstrom": round(maximum_coordinate_delta, 6),
        "bond_graph_and_orders_match_prepared_sdf": True,
        "element_order_matches_prepared_sdf": True,
    }


def atom_origin_record(atom: Chem.Atom) -> dict[str, Any]:
    record: dict[str, Any] = {"atom_index_zero_based": atom.GetIdx(), "element": atom.GetSymbol()}
    if atom.HasProp("_NNCSourceMap") and atom.GetProp("_NNCComponent") == "B1":
        map_id = int(atom.GetProp("_NNCSourceMap"))
        record.update({"component": "B1", "source_map_id": map_id})
        if atom.HasProp("_NNCMapRole"):
            record["source_map_role"] = atom.GetProp("_NNCMapRole")
    if atom.HasProp("_NNC_PDBSerial"):
        record["pdb_serial"] = int(atom.GetProp("_NNC_PDBSerial"))
    if atom.HasProp("_NNC_PDBIdentity"):
        record["pdb_identity"] = atom.GetProp("_NNC_PDBIdentity")
    if atom.GetSymbol() == "H" and atom.GetNeighbors():
        parent = atom.GetNeighbors()[0]
        record["parent_atom_index_zero_based"] = parent.GetIdx()
        if parent.HasProp("_NNCSourceMap") and parent.GetProp("_NNCComponent") == "B1":
            record["parent_source_map_id"] = int(parent.GetProp("_NNCSourceMap"))
        if parent.HasProp("_NNC_PDBSerial"):
            record["parent_pdb_serial"] = int(parent.GetProp("_NNC_PDBSerial"))
    return record


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--output-json", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--timeout-seconds", type=int, default=7200)
    args = parser.parse_args()

    input_paths = (
        B1_GRAPH, B1_MAP, B1_TOPOLOGY, FULL_GRAPH, FORMULA_AUDIT, B29_GRAPH, B29_MAP,
        RECONCILIATION, OPEN_PDB, LOCK_PATH, FULL_GRAPH_SCRIPT, B29_PILOT_HELPERS,
    )
    missing = [str(path) for path in input_paths if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"Missing required pilot inputs: {missing}")
    input_hashes_before = {path.name: sha256(path) for path in input_paths}

    map_legend = read_json(B1_MAP)
    topology = read_json(B1_TOPOLOGY)
    topology_hashes = topology["input_hashes_sha256"]
    for key, path in (("source_graph", B1_GRAPH), ("map_legend", B1_MAP)):
        recorded = topology_hashes[key]
        if recorded["path"] != str(path.relative_to(ROOT)):
            raise ValueError(f"B1 topology audit records a different {key} path")
        if recorded["sha256"] != input_hashes_before[path.name]:
            raise ValueError(f"B1 {key} hash differs from the topology audit")
    open_state = next(state for state in topology["candidate_states"] if state["state"] == "active_open")
    exact_mappings = open_state.get("all_exact_mappings", [])
    if open_state.get("exact_mapping_count") != 1 or len(exact_mappings) != 1:
        raise ValueError("Expected exactly one audited B1 mapping in the active/open candidate")
    mapping_entry = exact_mappings[0]
    reconciliation = read_json(RECONCILIATION)
    reconciliation_model = next(
        model for model in reconciliation["models"] if model["state"] == "active_open"
    )
    expected_pdb_sha = open_state["expected_candidate_pdb_sha256_from_reconciliation"]
    if input_hashes_before[OPEN_PDB.name] != expected_pdb_sha:
        raise ValueError("Active/open PDB hash differs from the B1 topology audit")
    if reconciliation_model["output_sha256"] != input_hashes_before[OPEN_PDB.name]:
        raise ValueError("Active/open PDB hash differs from the coordinate-reconciliation record")
    if reconciliation_model["output_file"] != str(OPEN_PDB.relative_to(ROOT)):
        raise ValueError("Coordinate-reconciliation record points to a different active/open candidate")
    formula_audit = read_json(FORMULA_AUDIT)
    full_graph_report = read_json(FULL_GRAPH)
    if full_graph_report.get("graph", {}).get("heavy_atom_count") != 507:
        raise ValueError("Full conjugate graph has an unexpected heavy-atom count")

    amberhome = Path(os.environ.get("AMBERHOME", "")).resolve()
    antechamber = shutil.which("antechamber")
    parmchk2 = shutil.which("parmchk2")
    if not amberhome.is_dir() or not antechamber or not parmchk2:
        raise RuntimeError("Run this script inside the pinned AmberTools environment")
    amber_lock = read_lock_package(LOCK_PATH, "ambertools")

    parent, _assembly = full_graph_audit.build_full_graph(
        formula_audit,
        map_legend,
        read_json(B29_MAP),
    )
    fragment = extract_local_fragment(parent)
    pdb_by_identity, pdb_by_serial = read_pdb_heavy_atoms(OPEN_PDB)
    hydrogenated, coordinate_rows = map_coordinates(fragment, mapping_entry, pdb_by_identity)
    input_formula = rdMolDescriptors.CalcMolFormula(hydrogenated)
    if input_formula != EXPECTED_FORMULA:
        raise ValueError(f"Hydrogenated input formula {input_formula}; expected {EXPECTED_FORMULA}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    input_sdf = args.output_dir / "b1_boundary_fragment.sdf"
    mol2_path = args.output_dir / "b1_boundary_fragment.mol2"
    frcmod_path = args.output_dir / "b1_boundary_fragment.frcmod"
    antechamber_log = args.output_dir / "antechamber.log"
    parmchk_log = args.output_dir / "parmchk2.log"
    Chem.MolToMolFile(hydrogenated, str(input_sdf))

    mol2_copied_this_run = False
    frcmod_copied_this_run = False
    with tempfile.TemporaryDirectory(prefix="dho-b1-boundary-") as temp_name:
        temp = Path(temp_name)
        local_sdf = temp / "input.sdf"
        local_mol2 = temp / "output.mol2"
        local_frcmod = temp / "output.frcmod"
        shutil.copyfile(input_sdf, local_sdf)
        ante_command = [
            antechamber,
            "-i", local_sdf.name,
            "-fi", "sdf",
            "-o", local_mol2.name,
            "-fo", "mol2",
            "-c", "bcc",
            "-nc", str(EXPECTED_FORMAL_CHARGE),
            "-at", "gaff2",
            "-pl", "20",
            "-ek", "maxcyc=0",
            "-rn", "B1P",
            "-s", "1",
        ]
        ante = run_checked(ante_command, temp, args.timeout_seconds)
        antechamber_log.write_text("STDOUT\n" + ante["stdout"] + "\nSTDERR\n" + ante["stderr"], encoding="utf-8")
        if ante["returncode"] == 0 and local_mol2.is_file():
            parm_command = [
                parmchk2,
                "-i", local_mol2.name,
                "-f", "mol2",
                "-o", local_frcmod.name,
                "-s", "2",
                "-frc", "ff14SB",
            ]
            parm = run_checked(parm_command, temp, args.timeout_seconds)
            parmchk_log.write_text("STDOUT\n" + parm["stdout"] + "\nSTDERR\n" + parm["stderr"], encoding="utf-8")
        else:
            parm = {"returncode": None, "stdout": "", "stderr": "Skipped after Antechamber failure", "timed_out": False}
            parmchk_log.write_text(parm["stderr"] + "\n", encoding="utf-8")
        if local_mol2.is_file():
            copy_normalized_tool_output(local_mol2, mol2_path)
            mol2_copied_this_run = True
        if local_frcmod.is_file():
            copy_normalized_tool_output(local_frcmod, frcmod_path)
            frcmod_copied_this_run = True

    molecule_records: list[dict[str, Any]] = []
    partial_charge_sum: float | None = None
    frcmod_result: dict[str, Any] = {"generated": False}
    mol2_result: dict[str, Any] = {"generated": False}
    output_validation_error: str | None = None
    tool_runs_succeeded = (
        ante["returncode"] == 0
        and not ante["timed_out"]
        and parm["returncode"] == 0
        and not parm["timed_out"]
    )
    mol2_candidate_validated = False
    if tool_runs_succeeded and not (mol2_copied_this_run and frcmod_copied_this_run):
        output_validation_error = "Successful tool runs did not produce both a MOL2 and frcmod output"

    if mol2_copied_this_run:
        mol2_metadata = {
            "path": str(mol2_path.relative_to(ROOT)),
            "sha256": sha256(mol2_path),
        }
        if not tool_runs_succeeded:
            mol2_result = {
                "generated": True,
                "candidate_accepted": False,
                "validation_status": "not_assessed_after_tool_failure",
                **mol2_metadata,
            }
        elif output_validation_error is None:
            try:
                output_atoms, output_bonds = read_mol2(mol2_path)
                names = mol2_atom_names(mol2_path)
                if len(names) != len(output_atoms):
                    raise ValueError("MOL2 atom-name count differs from the atom table")
                mapping_validation = validate_mol2_correspondence(hydrogenated, output_atoms, output_bonds)
                partial_charge_sum = float(sum(Decimal(atom["partial_charge_text"]) for atom in output_atoms))
                if not math.isclose(partial_charge_sum, EXPECTED_FORMAL_CHARGE, abs_tol=0.01):
                    raise ValueError(
                        f"MOL2 charge sum {partial_charge_sum:.6f} is outside the 0.01 e bookkeeping tolerance"
                    )
                atom_type_counts = Counter(atom["atom_type"] for atom in output_atoms)
                for index, (mol2_atom, mol2_name) in enumerate(zip(output_atoms, names, strict=True)):
                    record = atom_origin_record(hydrogenated.GetAtomWithIdx(index))
                    record.update({
                        "mol2_atom_name": mol2_name,
                        "gaff2_type": mol2_atom["atom_type"],
                        "partial_charge_e": mol2_atom["partial_charge_e"],
                    })
                    molecule_records.append(record)
                mol2_result = {
                    "generated": True,
                    "candidate_accepted": True,
                    "validation_status": "passed",
                    **mol2_metadata,
                    "atom_count": len(output_atoms),
                    "bond_count": len(output_bonds),
                    "heavy_atom_count": hydrogenated.GetNumHeavyAtoms(),
                    "mapping_validation": mapping_validation,
                    "partial_charge_sum_e": partial_charge_sum,
                    "net_charge_residual_e": partial_charge_sum - EXPECTED_FORMAL_CHARGE,
                    "charge_sum_within_0_01e_tolerance": True,
                    "charge_correction_applied": False,
                    "gaff2_atom_type_counts": dict(sorted(atom_type_counts.items())),
                    "atom_parameters": molecule_records,
                }
                mol2_candidate_validated = True
            except Exception as error:
                output_validation_error = f"MOL2 validation failed: {type(error).__name__}: {error}"
                mol2_result = {
                    "generated": True,
                    "candidate_accepted": False,
                    "validation_status": "failed",
                    **mol2_metadata,
                    "error": output_validation_error,
                }
                molecule_records = []
                partial_charge_sum = None
        else:
            mol2_result = {
                "generated": True,
                "candidate_accepted": False,
                "validation_status": "not_assessed",
                **mol2_metadata,
            }
    if frcmod_copied_this_run:
        frcmod_metadata = {
            "path": str(frcmod_path.relative_to(ROOT)),
            "sha256": sha256(frcmod_path),
        }
        if tool_runs_succeeded and mol2_candidate_validated:
            frcmod_result = {
                "generated": True,
                "candidate_accepted": True,
                "validation_status": "tool_and_MOL2_correspondence_checks_passed",
                **frcmod_metadata,
                **frcmod_inventory(frcmod_path),
            }
        else:
            frcmod_result = {
                "generated": True,
                "candidate_accepted": False,
                "validation_status": "not_assessed_with_a_complete_validated_MOL2_candidate",
                **frcmod_metadata,
            }

    pdb_hash_after = sha256(OPEN_PDB)
    if pdb_hash_after != input_hashes_before[OPEN_PDB.name]:
        raise ValueError("Active/open model-built candidate changed during the B1 pilot")
    input_hashes_after = {path.name: sha256(path) for path in input_paths}
    if input_hashes_after != input_hashes_before:
        raise ValueError("At least one checked source input changed during the B1 pilot")

    successful = (
        tool_runs_succeeded
        and mol2_candidate_validated
        and frcmod_result.get("candidate_accepted") is True
        and partial_charge_sum is not None
    )
    report = {
        "schema": "dia-health-one.nnc2215-b1-gaff2-boundary-pilot.v1",
        "status": "PILOT_ONLY; COMPLETE_CONJUGATE_FORCE_FIELD_NOT_VALIDATED" if successful else "PARTIAL_OR_FAILED_PILOT",
        "purpose": "Test pinned AmberTools GAFF2/AM1-BCC feasibility for a source-mapped B1 glucoside in a local PheB1/ValB2 peptide context.",
        "method": {
            "input_geometry": "active/open publisher-supplied model-built heavy-atom coordinates; no geometry optimization; explicit hydrogens added by RDKit",
            "fragment_definition": "B1 source maps 2–19 + PheB1 + ValB2 + AsnB3 N as a neutral NH2 cap",
            "declared_formal_charge": EXPECTED_FORMAL_CHARGE,
            "antechamber_charge_method": "AM1-BCC (Antechamber bcc)",
            "atom_type_set": "GAFF2",
            "bonded_parameter_lookup": "ff14SB available to parmchk2 as a comparison set; generated terms are not validated",
            "geometry_optimization": "disabled with maxcyc=0",
            "pdb_hydrogen_records_used": False,
            "nvidia_gpu_used": False,
        },
        "inputs": {
            "audit_script": {"path": str(Path(__file__).relative_to(ROOT)), "sha256": sha256(Path(__file__))},
            "full_graph_builder": {"path": str(FULL_GRAPH_SCRIPT.relative_to(ROOT)), "sha256": input_hashes_before[FULL_GRAPH_SCRIPT.name]},
            "parameterization_helper_script": {"path": str(B29_PILOT_HELPERS.relative_to(ROOT)), "sha256": input_hashes_before[B29_PILOT_HELPERS.name]},
            "source_graph": {"path": str(B1_GRAPH.relative_to(ROOT)), "sha256": input_hashes_before[B1_GRAPH.name]},
            "source_map_legend": {"path": str(B1_MAP.relative_to(ROOT)), "sha256": input_hashes_before[B1_MAP.name]},
            "b1_topology_audit": {"path": str(B1_TOPOLOGY.relative_to(ROOT)), "sha256": input_hashes_before[B1_TOPOLOGY.name]},
            "full_conjugate_graph": {"path": str(FULL_GRAPH.relative_to(ROOT)), "sha256": input_hashes_before[FULL_GRAPH.name]},
            "formula_audit": {"path": str(FORMULA_AUDIT.relative_to(ROOT)), "sha256": input_hashes_before[FORMULA_AUDIT.name]},
            "b29_graph_and_map": {
                "graph": {"path": str(B29_GRAPH.relative_to(ROOT)), "sha256": input_hashes_before[B29_GRAPH.name]},
                "map_legend": {"path": str(B29_MAP.relative_to(ROOT)), "sha256": input_hashes_before[B29_MAP.name]},
                "purpose": "Input required by the existing full-graph builder; the selected fragment contains no B29 atoms.",
            },
            "coordinate_reconciliation": {"path": str(RECONCILIATION.relative_to(ROOT)), "sha256": input_hashes_before[RECONCILIATION.name]},
            "active_open_model_candidate": {
                "path": str(OPEN_PDB.relative_to(ROOT)),
                "sha256": input_hashes_before[OPEN_PDB.name],
                "sha256_expected_by_b1_topology_audit": expected_pdb_sha,
                "hash_unchanged_after_run": pdb_hash_after == input_hashes_before[OPEN_PDB.name],
            },
            "ambertools_lock": {"path": str(LOCK_PATH.relative_to(ROOT)), "sha256": input_hashes_before[LOCK_PATH.name]},
            "ambertools_package": amber_lock,
        },
        "fragment": {
            "formula": input_formula,
            "heavy_atom_count": fragment.GetNumHeavyAtoms(),
            "atom_count_with_hydrogens": hydrogenated.GetNumAtoms(),
            "formal_charge": Chem.GetFormalCharge(fragment),
            "element_counts_heavy": dict(sorted(Counter(atom.GetSymbol() for atom in fragment.GetAtoms()).items())),
            "connected_component_count": len(Chem.GetMolFrags(fragment, asMols=False)),
            "stereochemistry_cip": {
                "PheB1_CA": "S",
                "ValB2_CA": "S",
                "B1_source_maps": {"9": "R", "11": "R", "13": "S", "15": "S", "17": "R"},
            },
            "boundary_hydrogen_counts": {"PheB1_N_graph_hydrogens": 1, "AsnB3_NH2_cap_graph_hydrogens": 2},
            "cap_interpretation": "The AsnB3 N is retained at its model coordinate after downstream atoms are removed; its valence is saturated as neutral NH2. This is an explicit truncation cap, not the full Asn residue environment.",
            "coordinate_mapping": {
                "state": "active_open",
                "mapping_id": mapping_entry["mapping_id"],
                "exact_mapping_count": open_state["exact_mapping_count"],
                "mapping_index_chosen_by_geometry": False,
                "source_map_1_is_protein_phe_n": True,
                "heavy_atom_coordinate_count": len(coordinate_rows),
                "heavy_atoms": coordinate_rows,
            },
        },
        "software": {
            "python_version": os.sys.version.split()[0],
            "rdkit_version": rdBase.rdkitVersion,
            "amberhome": str(amberhome),
            "antechamber_path": antechamber,
            "parmchk2_path": parmchk2,
        },
        "antechamber": {
            "command_template": [
                "antechamber", "-i", "input.sdf", "-fi", "sdf", "-o", "output.mol2", "-fo", "mol2",
                "-c", "bcc", "-nc", "0", "-at", "gaff2", "-pl", "20", "-ek", "maxcyc=0", "-rn", "B1P", "-s", "1",
            ],
            "returncode": ante["returncode"],
            "timed_out": ante["timed_out"],
            "log": {"path": str(antechamber_log.relative_to(ROOT)), "sha256": sha256(antechamber_log)},
        },
        "parmchk2": {
            "command_template": [
                "parmchk2", "-i", "output.mol2", "-f", "mol2", "-o", "output.frcmod", "-s", "2", "-frc", "ff14SB",
            ],
            "returncode": parm["returncode"],
            "timed_out": parm["timed_out"],
            "log": {"path": str(parmchk_log.relative_to(ROOT)), "sha256": sha256(parmchk_log)},
        },
        "sdf": {"path": str(input_sdf.relative_to(ROOT)), "sha256": sha256(input_sdf)},
        "mol2": mol2_result,
        "frcmod": frcmod_result,
        "output_validation_error": output_validation_error,
        "partial_charge_bookkeeping": {
            "formal_charge": EXPECTED_FORMAL_CHARGE,
            "partial_charge_sum_e": partial_charge_sum,
            "residual_e": None if partial_charge_sum is None else partial_charge_sum - EXPECTED_FORMAL_CHARGE,
            "charge_correction_applied": False,
            "interpretation": "AM1-BCC result for this fragment geometry and declared cap; not a validated final charge assignment.",
        },
        "limits": [
            "The active/open coordinates are model-built, and the pilot does not establish that this conformation is experimentally resolved.",
            "The AsnB3 NH2 cap omits the remaining AsnB3 environment and makes the AM1-BCC charges cap-dependent.",
            "The PheB1 hydrogen-record discrepancy is not resolved; PDB hydrogens are not used for proton placement.",
            "Generated GAFF2 types, AM1-BCC charges, and parmchk2 terms are feasibility outputs only and are not validated for transferability or protein-bound cross-boundary use.",
            "B29 mapping/protonation, full insulin microstate and total charge, complete force-field construction, target pose, MD, affinity, and clinical outcome remain unresolved.",
        ],
        "validated_conjugate_force_field_count": 0,
    }
    write_json(args.output_json, report)
    print(json.dumps({
        "status": report["status"],
        "fragment_formula": input_formula,
        "heavy_atom_count": fragment.GetNumHeavyAtoms(),
        "atom_count_with_hydrogens": hydrogenated.GetNumAtoms(),
        "partial_charge_sum_e": partial_charge_sum,
        "generated_parameter_records": frcmod_result.get("generated_parameter_record_count"),
        "validated_conjugate_force_field_count": 0,
        "output_json": str(args.output_json),
    }, indent=2))


if __name__ == "__main__":
    main()
