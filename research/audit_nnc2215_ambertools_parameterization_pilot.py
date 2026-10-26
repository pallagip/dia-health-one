#!/usr/bin/env python3
"""Run a bounded GAFF2/AM1-BCC pilot on the source-mapped B29 component.

The graph is protonated as a declared tri-carboxylate (-3) scenario. Its eight
equivalent Lewis labelings are parameterized separately to reveal sensitivity
to arbitrary carboxylate oxygen labels. This is not a complete insulin
parameterization and does not validate the generated parameters.
"""

from __future__ import annotations

import argparse
from collections import Counter
from decimal import Decimal
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
from typing import Any

from rdkit import Chem, rdBase
from rdkit.Chem import AllChem, rdMolDescriptors


ROOT = Path(__file__).resolve().parent
RESEARCH = ROOT
GRAPH_PATH = RESEARCH / "structures/nnc2215-b29-linker-roof-product.smi"
MAP_PATH = RESEARCH / "data/nnc2215-b29-linker-roof-product-map.json"
TOPOLOGY_PATH = RESEARCH / "data/nnc2215-b29-topology-audit.json"
OPEN_PDB = RESEARCH / "structures/nnc2215-open-source-corrected-candidate.pdb"
LOCK_PATH = RESEARCH / "pag-chemistry-environment-macos-arm64.lock.txt"
DEFAULT_OUTPUT_DIR = RESEARCH / "data/nnc2215-b29-gaff2-pilot"
DEFAULT_JSON = RESEARCH / "data/nnc2215-b29-gaff2-pilot.json"
CARBOXYLATES = ((52, 53, 54), (60, 61, 62), (68, 69, 70))
EXPECTED_FORMULA = "C62H69N16O14-3"
EXPECTED_FORMAL_CHARGE = -3


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def read_pdb_heavy_coordinates(path: Path) -> dict[int, tuple[float, float, float]]:
    result: dict[int, tuple[float, float, float]] = {}
    for line in path.read_text(encoding="ascii").splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            serial = int(line[6:11])
            element = line[76:78].strip().upper()
            if element in {"H", "D"}:
                continue
            result[serial] = (
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54]),
            )
    return result


def read_lock_package(lock_path: Path, name: str) -> dict[str, str]:
    pattern = re.compile(rf"/{re.escape(name)}-([^/#]+)\.conda#([0-9a-f]+)$")
    matches = [pattern.search(line.strip()) for line in lock_path.read_text().splitlines()]
    records = [match for match in matches if match]
    if len(records) != 1:
        raise ValueError(f"Expected one {name} conda package in {lock_path}, found {len(records)}")
    match = records[0]
    package_md5 = match.group(2)
    if len(package_md5) != 32:
        raise ValueError(f"Expected a 32-character explicit-lock MD5 for {name}, found {package_md5}")
    return {"version_build": match.group(1), "package_md5": package_md5}


def mapped_atom_indices(molecule: Chem.Mol) -> dict[int, int]:
    result: dict[int, int] = {}
    for atom in molecule.GetAtoms():
        map_id = atom.GetAtomMapNum()
        if map_id <= 0 or map_id in result:
            raise ValueError(f"Source graph requires unique positive atom maps; observed {map_id}")
        result[map_id] = atom.GetIdx()
    return result


def make_carboxylate_variant(
    neutral_graph: Chem.Mol,
    source_index_by_map: dict[int, int],
    anionic_oxygen_choices: tuple[int, int, int],
) -> Chem.Mol:
    molecule = Chem.RWMol(neutral_graph)
    if len(anionic_oxygen_choices) != len(CARBOXYLATES):
        raise ValueError("Exactly one anionic oxygen choice is required for each carboxylate")

    for (carbon_map, first_oxygen_map, second_oxygen_map), anionic_map in zip(
        CARBOXYLATES, anionic_oxygen_choices, strict=True
    ):
        if anionic_map not in {first_oxygen_map, second_oxygen_map}:
            raise ValueError(f"Map {anionic_map} is not an oxygen in carboxylate {carbon_map}")
        carbon_idx = source_index_by_map[carbon_map]
        oxygen_indices = (source_index_by_map[first_oxygen_map], source_index_by_map[second_oxygen_map])
        for oxygen_idx in oxygen_indices:
            oxygen = molecule.GetAtomWithIdx(oxygen_idx)
            oxygen.SetFormalCharge(-1 if oxygen.GetAtomMapNum() == anionic_map else 0)
            oxygen.SetNoImplicit(True)
            oxygen.SetNumExplicitHs(0)
            molecule.GetBondBetweenAtoms(carbon_idx, oxygen_idx).SetBondType(
                Chem.BondType.SINGLE if oxygen.GetAtomMapNum() == anionic_map else Chem.BondType.DOUBLE
            )

    result = molecule.GetMol()
    Chem.SanitizeMol(result)
    if Chem.GetFormalCharge(result) != EXPECTED_FORMAL_CHARGE:
        raise ValueError(f"Variant formal charge is {Chem.GetFormalCharge(result)}, expected -3")
    return result


def attach_open_candidate_coordinates(
    molecule: Chem.Mol,
    map_to_serial: dict[str, int],
    coordinates_by_serial: dict[int, tuple[float, float, float]],
    source_index_by_map: dict[int, int],
) -> Chem.Mol:
    conformer = Chem.Conformer(molecule.GetNumAtoms())
    conformer.Set3D(True)
    for map_id, atom_index in source_index_by_map.items():
        serial = int(map_to_serial[str(map_id)])
        if serial not in coordinates_by_serial:
            raise ValueError(f"No heavy-atom coordinates for PDB serial {serial}")
        conformer.SetAtomPosition(atom_index, coordinates_by_serial[serial])
    molecule.AddConformer(conformer, assignId=True)
    with_hydrogens = Chem.AddHs(molecule, addCoords=True)
    return with_hydrogens


def read_mol2(path: Path) -> tuple[list[dict[str, Any]], list[tuple[int, int, str]]]:
    section = ""
    atoms: list[dict[str, Any]] = []
    bonds: list[tuple[int, int, str]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("@<TRIPOS>"):
            section = line.removeprefix("@<TRIPOS>")
            continue
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split()
        if section == "ATOM":
            atom_type = fields[5]
            element = atom_type[0].upper()
            atoms.append(
                {
                    "id": int(fields[0]),
                    "name": fields[1],
                    "element": element,
                    "coordinates_angstrom": [float(value) for value in fields[2:5]],
                    "atom_type": atom_type,
                    "partial_charge_e": float(fields[8]),
                    "partial_charge_text": fields[8],
                }
            )
        elif section == "BOND":
            bonds.append((int(fields[1]) - 1, int(fields[2]) - 1, fields[3]))
    if not atoms or not bonds:
        raise ValueError(f"No MOL2 atoms or bonds found in {path}")
    return atoms, bonds


def frcmod_inventory(path: Path) -> dict[str, Any]:
    sections = {"MASS", "BOND", "ANGLE", "DIHE", "IMPROPER", "NONBON"}
    current = ""
    records: list[dict[str, str]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if stripped in sections:
            current = stripped
        elif current and stripped and not stripped.startswith(("!", "#")):
            records.append({"section": current, "line": stripped})
    penalties = [
        float(value)
        for record in records
        for value in re.findall(r"penalty score\s*=\s*([0-9.]+)", record["line"], flags=re.IGNORECASE)
    ]
    return {
        "generated_parameter_record_count": len(records),
        "records_with_attention_marker": sum("ATTN" in item["line"].upper() for item in records),
        "penalty_score_record_count": len(penalties),
        "maximum_penalty_score": max(penalties, default=None),
        "records_with_penalty_score_at_least_20": sum(value >= 20 for value in penalties),
        "generated_records": records,
    }


def run_checked(command: list[str], cwd: Path, timeout_seconds: int) -> dict[str, Any]:
    try:
        completed = subprocess.run(
            command,
            cwd=cwd,
            check=False,
            capture_output=True,
            text=True,
            timeout=timeout_seconds,
            env=os.environ.copy(),
        )
        return {
            "returncode": completed.returncode,
            "stdout": completed.stdout,
            "stderr": completed.stderr,
            "timed_out": False,
        }
    except subprocess.TimeoutExpired as error:
        return {
            "returncode": None,
            "stdout": error.stdout.decode() if isinstance(error.stdout, bytes) else error.stdout or "",
            "stderr": error.stderr.decode() if isinstance(error.stderr, bytes) else error.stderr or "",
            "timed_out": True,
        }


def write_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def copy_normalized_tool_output(source: Path, destination: Path) -> None:
    """Save whitespace-normalized MOL2/frcmod text without changing its records."""
    lines = source.read_text(encoding="utf-8").splitlines()
    normalized = "\n".join(line.rstrip() for line in lines).rstrip() + "\n"
    destination.write_text(normalized, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--output-json", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--variant-limit", type=int, default=8)
    parser.add_argument("--timeout-seconds", type=int, default=7200)
    args = parser.parse_args()
    if not 1 <= args.variant_limit <= 8:
        raise ValueError("--variant-limit must be between 1 and 8")

    for path in (GRAPH_PATH, MAP_PATH, TOPOLOGY_PATH, OPEN_PDB, LOCK_PATH):
        if not path.is_file():
            raise FileNotFoundError(path)

    graph_tokens = GRAPH_PATH.read_text(encoding="utf-8").split()
    if not graph_tokens:
        raise ValueError("The source-mapped B29 graph is empty")
    neutral_graph = Chem.MolFromSmiles(graph_tokens[0])
    if neutral_graph is None:
        raise ValueError("RDKit could not parse the source-mapped B29 graph")
    Chem.SanitizeMol(neutral_graph)
    if Chem.GetFormalCharge(neutral_graph) != 0 or neutral_graph.GetNumHeavyAtoms() != 92:
        raise ValueError("Unexpected neutral B29 product graph charge or heavy-atom count")
    source_index_by_map = mapped_atom_indices(neutral_graph)
    if set(source_index_by_map) != set(range(1, 93)):
        raise ValueError("Expected source maps 1 through 92 exactly")

    map_legend = read_json(MAP_PATH)
    if map_legend["graph_sha256"] != sha256(GRAPH_PATH):
        raise ValueError("Source graph hash does not match its approved map legend")
    declared_carboxylates = map_legend["map_id_legend"]["inherited_roof_map_ids_1_to_83"][
        "three_carboxyl_groups"
    ]
    if tuple(tuple(group) for group in declared_carboxylates) != CARBOXYLATES:
        raise ValueError("Declared carboxylate map groups differ from the map legend")

    topology = read_json(TOPOLOGY_PATH)
    open_state = next(state for state in topology["states"] if state["state"] == "active_open")
    exact_maps = open_state["complete_product_exact_graph_comparison"]["exact_mappings"]
    if len(exact_maps) != 16:
        raise ValueError(f"Expected 16 documented model-graph embeddings, found {len(exact_maps)}")
    expected_open_pdb_sha = open_state["coordinate_candidate_sha256_expected_by_reconciliation"]
    open_graph_sha = sha256(OPEN_PDB)
    if open_graph_sha != expected_open_pdb_sha:
        raise ValueError("The active/open PDB hash differs from the topology audit's coordinate candidate")
    selected_mapping = exact_maps[0]["map_id_to_pdb_serial"]
    coordinate_audit = read_pdb_heavy_coordinates(OPEN_PDB)

    amberhome = Path(os.environ.get("AMBERHOME", "")).resolve()
    antechamber = shutil.which("antechamber")
    parmchk2 = shutil.which("parmchk2")
    if not amberhome.is_dir() or not antechamber or not parmchk2:
        raise RuntimeError("Run this script inside the pinned AmberTools environment")
    amber_lock = read_lock_package(LOCK_PATH, "ambertools")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    variants: list[dict[str, Any]] = []
    all_charge_vectors: list[list[float]] = []
    source_atom_symbols: list[str] | None = None
    source_bond_pairs: set[tuple[int, int]] | None = None
    all_choices = list(itertools.product((53, 54), (61, 62), (69, 70)))[: args.variant_limit]

    for ordinal, choices in enumerate(all_choices, start=1):
        variant_id = f"resonance_{ordinal:02d}"
        mol = make_carboxylate_variant(neutral_graph, source_index_by_map, choices)
        if rdMolDescriptors.CalcMolFormula(mol) != EXPECTED_FORMULA:
            raise ValueError(
                f"{variant_id} formula is {rdMolDescriptors.CalcMolFormula(mol)}, expected {EXPECTED_FORMULA}"
            )
        mol = attach_open_candidate_coordinates(
            mol, selected_mapping, coordinate_audit, source_index_by_map
        )
        atom_symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]
        bond_pairs = {
            tuple(sorted((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())))
            for bond in mol.GetBonds()
        }
        if source_atom_symbols is None:
            source_atom_symbols = atom_symbols
            source_bond_pairs = bond_pairs
        elif atom_symbols != source_atom_symbols or bond_pairs != source_bond_pairs:
            raise ValueError("Resonance variants changed atom ordering or graph connectivity")

        final_sdf = args.output_dir / f"{variant_id}.sdf"
        final_mol2 = args.output_dir / f"{variant_id}.mol2"
        final_frcmod = args.output_dir / f"{variant_id}.frcmod"
        final_antechamber_log = args.output_dir / f"{variant_id}.antechamber.log"
        final_parmchk_log = args.output_dir / f"{variant_id}.parmchk2.log"
        Chem.MolToMolFile(mol, str(final_sdf))

        with tempfile.TemporaryDirectory(prefix=f"dho-{variant_id}-") as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            local_sdf = temp_dir / "input.sdf"
            local_mol2 = temp_dir / "output.mol2"
            local_frcmod = temp_dir / "output.frcmod"
            shutil.copyfile(final_sdf, local_sdf)
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
                "-rn", "B29",
                "-s", "1",
            ]
            ante = run_checked(ante_command, temp_dir, args.timeout_seconds)
            final_antechamber_log.write_text(
                "STDOUT\n" + ante["stdout"] + "\nSTDERR\n" + ante["stderr"],
                encoding="utf-8",
            )
            if ante["returncode"] == 0 and local_mol2.is_file():
                parm_command = [
                    parmchk2,
                    "-i", local_mol2.name,
                    "-f", "mol2",
                    "-o", local_frcmod.name,
                    "-s", "2",
                    "-frc", "ff14SB",
                ]
                parm = run_checked(parm_command, temp_dir, args.timeout_seconds)
                final_parmchk_log.write_text(
                    "STDOUT\n" + parm["stdout"] + "\nSTDERR\n" + parm["stderr"],
                    encoding="utf-8",
                )
            else:
                parm = {"returncode": None, "stdout": "", "stderr": "Skipped after Antechamber failure", "timed_out": False}

            if local_mol2.is_file():
                copy_normalized_tool_output(local_mol2, final_mol2)
            if local_frcmod.is_file():
                copy_normalized_tool_output(local_frcmod, final_frcmod)

        record: dict[str, Any] = {
            "variant_id": variant_id,
            "anionic_oxygen_map_ids": list(choices),
            "formula": rdMolDescriptors.CalcMolFormula(mol),
            "formal_charge": Chem.GetFormalCharge(mol),
            "atom_count_with_hydrogens": mol.GetNumAtoms(),
            "bond_count_with_hydrogens": mol.GetNumBonds(),
            "input_sdf": {"path": str(final_sdf.relative_to(ROOT)), "sha256": sha256(final_sdf)},
            "antechamber": {
                "command_template": [
                    "antechamber", "-i", "input.sdf", "-fi", "sdf", "-o", "output.mol2", "-fo", "mol2",
                    "-c", "bcc", "-nc", "-3", "-at", "gaff2", "-pl", "20", "-ek", "maxcyc=0", "-rn", "B29", "-s", "1",
                ],
                "returncode": ante["returncode"],
                "timed_out": ante["timed_out"],
                "log_sha256": sha256(final_antechamber_log),
            },
            "parmchk2": {
                "command_template": [
                    "parmchk2", "-i", "output.mol2", "-f", "mol2", "-o", "output.frcmod", "-s", "2", "-frc", "ff14SB",
                ],
                "returncode": parm["returncode"],
                "timed_out": parm["timed_out"],
                "log_sha256": sha256(final_parmchk_log),
            },
        }

        if final_mol2.is_file():
            output_atoms, output_bonds = read_mol2(final_mol2)
            if len(output_atoms) != mol.GetNumAtoms() or len(output_bonds) != mol.GetNumBonds():
                raise ValueError(f"{variant_id} MOL2 atom/bond count differs from its SDF")
            output_symbols = [atom["element"] for atom in output_atoms]
            if output_symbols != atom_symbols:
                raise ValueError(f"{variant_id} MOL2 atom element order differs from its SDF")
            charges = [atom["partial_charge_e"] for atom in output_atoms]
            charge_sum = float(sum(Decimal(atom["partial_charge_text"]) for atom in output_atoms))
            atom_type_counts = Counter(atom["atom_type"] for atom in output_atoms)
            if not math.isclose(charge_sum, EXPECTED_FORMAL_CHARGE, abs_tol=0.01):
                raise ValueError(f"{variant_id} partial-charge sum {charge_sum:.6f} differs from -3")
            record["mol2"] = {
                "path": str(final_mol2.relative_to(ROOT)),
                "sha256": sha256(final_mol2),
                "atom_count": len(output_atoms),
                "bond_count": len(output_bonds),
                "partial_charge_sum_e": charge_sum,
                "net_charge_residual_e": charge_sum - EXPECTED_FORMAL_CHARGE,
                "charge_sum_within_0_01e_tolerance": math.isclose(
                    charge_sum, EXPECTED_FORMAL_CHARGE, abs_tol=0.01
                ),
                "atom_type_counts": dict(sorted(atom_type_counts.items())),
                "charged_atoms_by_source_map": [
                    {
                        "source_map_id": neutral_graph.GetAtomWithIdx(index).GetAtomMapNum()
                        if index < neutral_graph.GetNumAtoms() else None,
                        "atom_type": output_atoms[index]["atom_type"],
                        "partial_charge_e": output_atoms[index]["partial_charge_e"],
                    }
                    for index in range(min(len(output_atoms), neutral_graph.GetNumAtoms()))
                ],
            }
            all_charge_vectors.append(charges)
        else:
            record["mol2"] = {"generated": False}

        if final_frcmod.is_file():
            record["frcmod"] = {
                "path": str(final_frcmod.relative_to(ROOT)),
                "sha256": sha256(final_frcmod),
                **frcmod_inventory(final_frcmod),
            }
        else:
            record["frcmod"] = {"generated": False}
        variants.append(record)

    charge_sensitivity: dict[str, Any]
    if len(all_charge_vectors) == len(all_choices) and all_charge_vectors:
        charge_sensitivity = {
            "compared_variant_count": len(all_charge_vectors),
            "largest_per_atom_absolute_charge_range_e": max(
                max(vector[index] for vector in all_charge_vectors)
                - min(vector[index] for vector in all_charge_vectors)
                for index in range(len(all_charge_vectors[0]))
            ),
            "per_source_atom": [
                {
                    "atom_index_zero_based": index,
                    "source_map_id": neutral_graph.GetAtomWithIdx(index).GetAtomMapNum()
                    if index < neutral_graph.GetNumAtoms() else None,
                    "min_partial_charge_e": min(vector[index] for vector in all_charge_vectors),
                    "max_partial_charge_e": max(vector[index] for vector in all_charge_vectors),
                    "range_e": max(vector[index] for vector in all_charge_vectors)
                    - min(vector[index] for vector in all_charge_vectors),
                }
                for index in range(len(all_charge_vectors[0]))
                if max(vector[index] for vector in all_charge_vectors)
                - min(vector[index] for vector in all_charge_vectors) > 1e-5
            ],
        }
    else:
        charge_sensitivity = {"compared_variant_count": len(all_charge_vectors), "complete": False}

    successful = all(
        item["antechamber"]["returncode"] == 0
        and item["parmchk2"]["returncode"] == 0
        and "path" in item["mol2"]
        and "path" in item["frcmod"]
        for item in variants
    ) and len(variants) == 8
    if sha256(OPEN_PDB) != open_graph_sha:
        raise ValueError("The open model-built coordinate candidate changed during the pilot")
    report = {
        "schema": "dia-health-one.nnc2215-b29-gaff2-parameterization-pilot.v1",
        "status": "PILOT_ONLY; COMPLETE_CONJUGATE_FORCE_FIELD_NOT_VALIDATED" if successful else "PARTIAL_OR_FAILED_PILOT",
        "purpose": "Test AmberTools GAFF2/AM1-BCC feasibility on an isolated source-mapped B29 conjugate component under the conditional tri-carboxylate charge state.",
        "inputs": {
            "audit_script": {"path": str(Path(__file__).relative_to(ROOT)), "sha256": sha256(Path(__file__))},
            "product_graph": {"path": str(GRAPH_PATH.relative_to(RESEARCH)), "sha256": sha256(GRAPH_PATH)},
            "product_map_legend": {"path": str(MAP_PATH.relative_to(RESEARCH)), "sha256": sha256(MAP_PATH)},
            "topology_audit": {"path": str(TOPOLOGY_PATH.relative_to(RESEARCH)), "sha256": sha256(TOPOLOGY_PATH)},
            "open_model_coordinate_candidate": {
                "path": str(OPEN_PDB.relative_to(RESEARCH)),
                "sha256": open_graph_sha,
                "sha256_expected_by_topology_audit": expected_open_pdb_sha,
                "mapping_id": exact_maps[0]["mapping_id"],
                "one_of_sixteen_graph_embeddings": True,
                "coordinates_are_model_built": True,
                "hash_unchanged_after_run": True,
            },
            "environment_lock": {"path": str(LOCK_PATH.relative_to(RESEARCH)), "sha256": sha256(LOCK_PATH)},
            "ambertools_package": amber_lock,
        },
        "method": {
            "tools": "AmberTools 26 Antechamber and parmchk2",
            "charge_method": "AM1-BCC (Antechamber bcc)",
            "atom_type_set": "GAFF2",
            "bonded_parameter_lookups": "GAFF2 with ff14SB available as a comparison set; parmchk2-generated terms are not validated",
            "declared_charge_state": "three B29 carboxylates deprotonated; neutral NH2 source endpoint cap; net formal charge -3",
            "input_geometry": "heavy-atom coordinates copied from the public active/open model-built candidate using the first documented exact graph embedding; explicit hydrogens added by RDKit",
            "sqm_geometry_optimization": "disabled with maxcyc=0 to avoid changing the selected model geometry in a vacuum; the resulting charges remain geometry- and model-dependent",
            "resonance_variants": 8,
            "carboxylate_map_groups": [list(group) for group in CARBOXYLATES],
            "rdkit_version": rdBase.rdkitVersion,
            "pdb_or_source_coordinates_modified": False,
        },
        "component_graph": {
            "formula": EXPECTED_FORMULA,
            "formal_charge": EXPECTED_FORMAL_CHARGE,
            "heavy_atom_count": neutral_graph.GetNumHeavyAtoms(),
            "resonance_variant_count_completed": len(variants),
            "atom_count_with_hydrogens": variants[0]["atom_count_with_hydrogens"] if variants else None,
            "heavy_atom_input_mapping_is_authoritative": False,
        },
        "variants": variants,
        "carboxylate_resonance_charge_sensitivity": charge_sensitivity,
        "assessment": {
            "all_eight_variants_completed": successful,
            "gaff2_partial_charges_assigned": successful,
            "bonded_parameter_gaps_inventoried": successful,
            "validated_conjugate_force_fields": 0,
            "full_insulin_topology_or_openmm_system_created": False,
            "md_run": False,
            "clinical_use": False,
            "limits": [
                "The B29 source endpoint NH2 cap is not the attached LysB29 NZ-in-residue environment; atom charges at the protein junction are therefore not final.",
                "The chosen open coordinates are publisher-supplied model-built coordinates, and their source-to-coordinate map is one of 16 exact element-adjacency embeddings.",
                "AM1-BCC/GAFF2 plus parmchk2 output is a candidate fragment assignment, not a reviewed force field.",
                "The eight outputs test arbitrary carboxylate Lewis-label sensitivity; no charge averaging or atom-type normalization is performed.",
                "Custom LysB29 boundary-residue charges and cross-boundary bonded terms remain unassigned.",
                "PheB1/B1 fragment parameters, insulin termini, protein compatibility, full system assembly, and force-field validation remain open.",
                "No solvent, ions, target receptor complex, MD, affinity, pharmacology, or clinical conclusion is included.",
            ],
        },
    }
    write_json(args.output_json, report)
    print(json.dumps({
        "status": report["status"],
        "variants_completed": len(variants),
        "net_formal_charge": EXPECTED_FORMAL_CHARGE,
        "max_charge_sensitivity_e": charge_sensitivity.get("largest_per_atom_absolute_charge_range_e"),
        "validated_conjugate_force_fields": 0,
        "output_json": str(args.output_json),
    }, indent=2))


if __name__ == "__main__":
    main()
