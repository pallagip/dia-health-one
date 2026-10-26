#!/usr/bin/env python3
"""Audit OpenFF Sage parameter-label coverage on source-derived NNC2215 graphs.

This is a chemical-perception coverage check. It does not assign partial
charges, create coordinates or an OpenMM system, validate transferability, or
run molecular dynamics.
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import re
import sys
from typing import Any

from openff.toolkit import ForceField, Molecule, Topology
from openff.toolkit.utils.exceptions import UnassignedValenceParameterException
from openmm import __version__ as openmm_version
from rdkit import Chem, rdBase
from rdkit.Chem import rdMolDescriptors
from openff import toolkit as openff_toolkit
from openff import interchange as openff_interchange
from openff import units as openff_units


ROOT = Path(__file__).resolve().parent
GRAPH_PATH = ROOT / "data" / "nnc2215-full-conjugate-graph.json"
ROOF_PATH = ROOT / "structures" / "nnc2215-roof-propyl-azide.smi"
ROOF_MAP_PATH = ROOT / "data" / "nnc2215-roof-macrocycle-atom-map.json"
LOCK_PATH = ROOT / "pag-chemistry-environment-macos-arm64.lock.txt"
OUTPUT_PATH = ROOT / "data" / "nnc2215-forcefield-label-coverage.json"

EXPECTED_VERSIONS = {
    "openff-toolkit": "0.19.0",
    "openff-interchange": "0.5.5",
    "openff-units": "0.4.0",
    "openff-forcefields": "2026.09.0",
    "openmm": "8.6.1",
    "rdkit": "2026.03.1",
}
EXPECTED_FULL_FORMULA = "C325H461N79O97S6"
EXPECTED_ROOF_FORMULA = "C56H63N15O12"
CARBOXYLATE_OXYGEN_MAP_IDS = (54, 62, 70)
LABEL_BUCKETS = (
    "Constraints",
    "Bonds",
    "Angles",
    "ProperTorsions",
    "ImproperTorsions",
    "vdW",
    "Electrostatics",
)


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def versions_in_lock(lock_raw: bytes) -> dict[str, str]:
    filenames = [
        line.strip().split("#", 1)[0].rsplit("/", 1)[-1].removesuffix(".conda")
        for line in lock_raw.decode("utf-8").splitlines()
        if line.startswith("http")
    ]
    locked: dict[str, str] = {}
    for package in EXPECTED_VERSIONS:
        matches = []
        for filename in filenames:
            if not filename.startswith(f"{package}-"):
                continue
            version = filename[len(package) + 1 :].split("-", 1)[0]
            if re.fullmatch(r"[0-9][A-Za-z0-9.]*", version):
                matches.append(version)
        if len(matches) != 1:
            raise ValueError(f"Expected one exact lock entry for {package}; found {matches}")
        locked[package] = matches[0]

    mismatches = {
        package: {"expected": expected, "locked": locked.get(package)}
        for package, expected in EXPECTED_VERSIONS.items()
        if locked.get(package) != expected
    }
    if mismatches:
        raise ValueError(f"Environment version declaration disagrees with lockfile: {mismatches}")
    return locked


def load_locked_environment() -> dict[str, str]:
    prefix = os.environ.get("CONDA_PREFIX")
    if not prefix:
        raise ValueError("Run inside the locked conda environment (for example, with micromamba run).")

    metadata_dir = Path(prefix) / "conda-meta"
    observed: dict[str, str] = {}
    metadata_records = []
    for path in metadata_dir.glob("*.json"):
        metadata = json.loads(path.read_text(encoding="utf-8"))
        metadata_records.append(metadata)
    for package in EXPECTED_VERSIONS:
        matches = [item for item in metadata_records if item.get("name") == package]
        if len(matches) != 1:
            raise ValueError(f"Expected one conda metadata record for {package}; found {len(matches)}")
        observed[package] = matches[0]["version"]

    mismatches = {
        package: {"expected": expected, "observed": observed.get(package)}
        for package, expected in EXPECTED_VERSIONS.items()
        if observed.get(package) != expected
    }
    if mismatches:
        raise ValueError(f"Environment does not match the audited lock versions: {mismatches}")

    return observed


def build_full_graph(graph: dict[str, Any]) -> Chem.Mol:
    atoms = graph["atoms"]
    bonds = graph["bonds"]
    if len(atoms) != 507 or len(bonds) != 527:
        raise ValueError(f"Unexpected source graph dimensions: {len(atoms)} atoms, {len(bonds)} bonds")

    editable = Chem.RWMol()
    node_to_index: dict[str, int] = {}
    for atom_data in atoms:
        atom = Chem.Atom(atom_data["element"])
        atom.SetFormalCharge(int(atom_data["formal_charge"]))
        atom.SetIsAromatic(bool(atom_data["is_aromatic"]))
        atom.SetNoImplicit(True)
        atom.SetNumExplicitHs(int(atom_data["hydrogen_count_not_as_nodes"]))
        if atom.GetAtomicNum() != int(atom_data["atomic_number"]):
            raise ValueError(f"Element/atomic-number mismatch in {atom_data['node_id']}")
        node_to_index[atom_data["node_id"]] = editable.AddAtom(atom)

    bond_types = {
        "SINGLE": Chem.BondType.SINGLE,
        "DOUBLE": Chem.BondType.DOUBLE,
        "TRIPLE": Chem.BondType.TRIPLE,
        "AROMATIC": Chem.BondType.AROMATIC,
    }
    for bond_data in bonds:
        left, right = bond_data["atom_node_ids"]
        bond_type_name = bond_data["bond_type"]
        if bond_type_name not in bond_types:
            raise ValueError(f"Unsupported bond type {bond_type_name!r} in {bond_data['edge_id']}")
        editable.AddBond(node_to_index[left], node_to_index[right], bond_types[bond_type_name])
        bond = editable.GetBondBetweenAtoms(node_to_index[left], node_to_index[right])
        bond.SetIsAromatic(bool(bond_data["is_aromatic"]))
        bond.SetIsConjugated(bool(bond_data["is_conjugated"]))

    molecule = editable.GetMol()
    Chem.SanitizeMol(molecule)
    formula = rdMolDescriptors.CalcMolFormula(molecule)
    if formula != EXPECTED_FULL_FORMULA or formula != graph["formula"]["rdkit_formula"]:
        raise ValueError(f"Reconstructed full-graph formula {formula} does not match source ledger")
    if molecule.GetNumHeavyAtoms() != 507 or molecule.GetNumAtoms() != 507:
        raise ValueError("Source graph must contain 507 heavy atoms and no hydrogen nodes before H materialization")
    return molecule


def load_roof_graph(raw: bytes) -> Chem.Mol:
    lines = [
        line.strip()
        for line in raw.decode("utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    if len(lines) != 1:
        raise ValueError(f"Expected one mapped roof SMILES line; found {len(lines)}")
    fields = lines[0].split()
    if len(fields) != 2:
        raise ValueError("Expected mapped roof SMILES followed by a molecule name")

    params = Chem.SmilesParserParams()
    params.removeHs = False
    molecule = Chem.MolFromSmiles(fields[0], params)
    if molecule is None:
        raise ValueError("RDKit could not parse the mapped roof graph")
    Chem.SanitizeMol(molecule)

    map_ids = [atom.GetAtomMapNum() for atom in molecule.GetAtoms() if atom.GetAtomicNum() > 1]
    if len(map_ids) != 83 or set(map_ids) != set(range(1, 84)):
        raise ValueError("Mapped roof graph must have unique and complete heavy-atom map IDs 1–83")
    formula = rdMolDescriptors.CalcMolFormula(molecule)
    if formula != EXPECTED_ROOF_FORMULA:
        raise ValueError(f"Roof graph formula {formula} does not match {EXPECTED_ROOF_FORMULA}")
    return molecule


def make_carboxylate_state(neutral: Chem.Mol, deprotonated_map_ids: tuple[int, ...]) -> Chem.Mol:
    molecule = Chem.Mol(neutral)
    for map_id in deprotonated_map_ids:
        matches = [atom for atom in molecule.GetAtoms() if atom.GetAtomMapNum() == map_id]
        if len(matches) != 1 or matches[0].GetSymbol() != "O":
            raise ValueError(f"Expected one mapped acidic oxygen at map ID {map_id}")
        oxygen = matches[0]
        if oxygen.GetTotalNumHs() != 1 or oxygen.GetFormalCharge() != 0:
            raise ValueError(f"Map ID {map_id} is not a neutral COOH oxygen in the source graph")
        oxygen.SetFormalCharge(-1)
        oxygen.SetNumExplicitHs(0)

    Chem.SanitizeMol(molecule)
    expected_charge = -len(deprotonated_map_ids)
    charge = sum(atom.GetFormalCharge() for atom in molecule.GetAtoms())
    if charge != expected_charge:
        raise ValueError(f"Roof state charge {charge} != expected {expected_charge}")
    return molecule


def count_labels(force_field: ForceField, rdkit_molecule: Chem.Mol) -> dict[str, int]:
    with_hydrogens = Chem.AddHs(rdkit_molecule)
    off_molecule = Molecule.from_rdkit(
        with_hydrogens,
        allow_undefined_stereo=True,
        hydrogens_are_explicit=True,
    )
    labeled = force_field.label_molecules(Topology.from_molecules(off_molecule))
    if len(labeled) != 1:
        raise ValueError(f"Expected labels for one molecule; received {len(labeled)} results")
    molecule_labels = labeled[0]
    counts = {
        bucket: len(molecule_labels.get(bucket, {}))
        for bucket in LABEL_BUCKETS
    }
    if counts["Bonds"] == 0 or counts["Angles"] == 0 or counts["vdW"] == 0:
        raise ValueError(f"Sage did not return expected chemical parameter labels: {counts}")
    return counts


def main() -> int:
    try:
        versions = load_locked_environment()
        if openff_toolkit.__version__ != versions["openff-toolkit"]:
            raise ValueError("OpenFF Toolkit module version disagrees with conda metadata")
        if openff_interchange.__version__ != versions["openff-interchange"]:
            raise ValueError("OpenFF Interchange module version disagrees with conda metadata")
        if openff_units.__version__ != versions["openff-units"]:
            raise ValueError("OpenFF Units module version disagrees with conda metadata")
        if openmm_version != versions["openmm"] or rdBase.rdkitVersion != versions["rdkit"]:
            raise ValueError("OpenMM or RDKit module version disagrees with conda metadata")

        graph_raw = GRAPH_PATH.read_bytes()
        graph = json.loads(graph_raw)
        roof_raw = ROOF_PATH.read_bytes()
        roof_map_raw = ROOF_MAP_PATH.read_bytes()
        lock_raw = LOCK_PATH.read_bytes()
        locked_versions = versions_in_lock(lock_raw)
        if locked_versions != versions:
            raise ValueError("Installed package versions do not match the explicit environment lock")
        full_graph = build_full_graph(graph)
        roof = load_roof_graph(roof_raw)
        roof_map = json.loads(roof_map_raw)
        expected_carboxyls = roof_map["map_id_legend"]["three_carboxyl_groups"]
        listed_oxygen_ids = tuple(group[2] for group in expected_carboxyls)
        if listed_oxygen_ids != CARBOXYLATE_OXYGEN_MAP_IDS:
            raise ValueError(f"Carboxylate oxygen map IDs differ from the source map legend: {listed_oxygen_ids}")

        force_field = ForceField("openff-2.3.0.offxml")
        full_counts = count_labels(force_field, full_graph)

        roof_states: list[dict[str, Any]] = []
        for mask in range(1 << len(CARBOXYLATE_OXYGEN_MAP_IDS)):
            deprotonated = tuple(
                map_id
                for bit, map_id in enumerate(CARBOXYLATE_OXYGEN_MAP_IDS)
                if mask & (1 << bit)
            )
            state_graph = make_carboxylate_state(roof, deprotonated)
            state_formula = rdMolDescriptors.CalcMolFormula(state_graph)
            labels = count_labels(force_field, state_graph)
            formal_charge = sum(atom.GetFormalCharge() for atom in state_graph.GetAtoms())
            roof_states.append({
                "deprotonated_acid_oxygen_map_ids": list(deprotonated),
                "formula_before_hydrogen_materialization": state_formula,
                "heavy_atom_count": state_graph.GetNumHeavyAtoms(),
                "atom_count_after_hydrogen_materialization": Chem.AddHs(state_graph).GetNumAtoms(),
                "formal_charge": formal_charge,
                "sage_label_counts": labels,
                "label_coverage": "pass",
            })

        failures = [state for state in roof_states if state["label_coverage"] != "pass"]
        if failures:
            raise ValueError(f"Roof fragment label coverage failed: {failures}")

        output = {
            "schema": "dia-health-one.nnc2215-forcefield-label-coverage.v1",
            "scope": "OpenFF Sage 2.3.0 chemical-perception label coverage only",
            "inputs": {
                "full_conjugate_graph": {
                    "path": str(GRAPH_PATH.relative_to(ROOT)),
                    "sha256": sha256_bytes(graph_raw),
                },
                "roof_mapped_smiles": {
                    "path": str(ROOF_PATH.relative_to(ROOT)),
                    "sha256": sha256_bytes(roof_raw),
                },
                "roof_atom_map_legend": {
                    "path": str(ROOF_MAP_PATH.relative_to(ROOT)),
                    "sha256": sha256_bytes(roof_map_raw),
                },
                "environment_lock": {
                    "path": str(LOCK_PATH.relative_to(ROOT)),
                    "sha256": sha256_bytes(lock_raw),
                },
                "audit_script_sha256": sha256_file(Path(__file__).resolve()),
            },
            "software": {
                "python_version": sys.version.split()[0],
                "locked_package_versions": versions,
                "sage_force_field": "openff-2.3.0.offxml",
            },
            "full_neutral_constitutional_graph": {
                "formula": rdMolDescriptors.CalcMolFormula(full_graph),
                "heavy_atom_count": full_graph.GetNumHeavyAtoms(),
                "atom_count_after_hydrogen_materialization": Chem.AddHs(full_graph).GetNumAtoms(),
                "formal_charge": sum(atom.GetFormalCharge() for atom in full_graph.GetAtoms()),
                "sage_label_counts": full_counts,
                "label_coverage": "pass",
                "coordinates_present": False,
                "stereochemistry_fully_specified": False,
            },
            "roof_macrocycle_fragment_carboxyl_states": {
                "oxygen_map_ids_tested": list(CARBOXYLATE_OXYGEN_MAP_IDS),
                "state_count": len(roof_states),
                "states": roof_states,
                "interpretation": "All eight combinations label successfully; no pKa weighting or dominant-state inference was performed.",
            },
            "not_established": [
                "partial atomic charges or electrostatic validation",
                "validated parameters or transferability for the insulin–conjugate hybrid",
                "compatibility of Sage conjugate parameters with the protein force-field files",
                "complete stereochemistry, protonation, or experimental atom mapping",
                "3D coordinates, geometry quality, minimization, or an OpenMM System",
                "molecular dynamics, binding pose, affinity, or free energy",
                "pharmacological efficacy or clinical safety",
            ],
            "assessment": {
                "status": "SAGE_LABEL_COVERAGE_ONLY",
                "partial_charges_assigned": False,
                "coordinates_created": False,
                "openmm_system_created": False,
                "molecular_dynamics_run": False,
                "validated_hybrid_force_field": False,
                "clinical_use": False,
            },
        }
        OUTPUT_PATH.write_text(json.dumps(output, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        print(json.dumps({
            "audit": str(OUTPUT_PATH),
            "status": output["assessment"]["status"],
            "full_graph_formula": output["full_neutral_constitutional_graph"]["formula"],
            "full_graph_atom_count": output["full_neutral_constitutional_graph"]["atom_count_after_hydrogen_materialization"],
            "full_graph_label_counts": full_counts,
            "roof_state_count": len(roof_states),
            "roof_states_all_labeled": all(state["label_coverage"] == "pass" for state in roof_states),
            "toolkit_version": versions["openff-toolkit"],
            "sage_force_field": output["software"]["sage_force_field"],
        }, indent=2))
        return 0
    except (OSError, ValueError, json.JSONDecodeError, UnassignedValenceParameterException) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
