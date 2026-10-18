#!/usr/bin/env python3
"""Audit a GAFF2/AM1-BCC PaG carboxylate parameter-file export.

This checks atom-order/connectivity correspondence, stereochemistry from the
exported coordinates, expected carboxylate resonance representation, total
charge, and OpenMM Amber-file loading. It is a ligand-only preparation audit;
it does not validate protein compatibility, binding energetics, or NNC2215.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path

from rdkit import Chem, rdBase
from rdkit.Chem import rdMolDescriptors


ROOT = Path(__file__).parent
DEFAULT_SDF = ROOT / "ligands/pag-pubchem-cid83944-carboxylate-3d.sdf"
DEFAULT_MOL2 = ROOT / "ligands/pag-gaff2/pag-carboxylate-gaff2.mol2"
DEFAULT_FRCMOD = ROOT / "ligands/pag-gaff2/pag-carboxylate-gaff2.frcmod"
DEFAULT_PRMTOP = ROOT / "ligands/pag-gaff2/pag-carboxylate-gaff2.prmtop"
DEFAULT_INPCRD = ROOT / "ligands/pag-gaff2/pag-carboxylate-gaff2.inpcrd"
DEFAULT_OUTPUT = ROOT / "data/pag-gaff2-validation.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_mol2(path: Path) -> tuple[list[dict], list[tuple[int, int, str]]]:
    section = ""
    atoms: list[dict] = []
    bonds: list[tuple[int, int, str]] = []
    for line in path.read_text().splitlines():
        if line.startswith("@<TRIPOS>"):
            section = line.removeprefix("@<TRIPOS>")
            continue
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split()
        if section == "ATOM":
            atoms.append(
                {
                    "id": int(fields[0]),
                    "name": fields[1],
                    "element": fields[1][0].upper() + fields[1][1:].lower()
                    if fields[1][1:2].islower()
                    else fields[1][0].upper(),
                    "xyz": tuple(float(value) for value in fields[2:5]),
                    "atom_type": fields[5],
                    "partial_charge": float(fields[8]),
                }
            )
        elif section == "BOND":
            bonds.append((int(fields[1]) - 1, int(fields[2]) - 1, fields[3]))
    if not atoms or not bonds:
        raise ValueError(f"No atom or bond records found in {path}.")
    return atoms, bonds


def cip_labels_from_tags(molecule: Chem.Mol) -> list[dict]:
    Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)
    return [
        {"atom_index_zero_based": atom.GetIdx(), "cip": atom.GetProp("_CIPCode")}
        for atom in molecule.GetAtoms()
        if atom.HasProp("_CIPCode")
    ]


def cip_labels_from_coordinates(molecule: Chem.Mol, coordinates: list[tuple]) -> list[dict]:
    result = Chem.Mol(molecule)
    if result.GetNumAtoms() != len(coordinates):
        raise ValueError("Coordinate count differs from the SDF atom count.")
    conformer = result.GetConformer()
    for index, xyz in enumerate(coordinates):
        conformer.SetAtomPosition(index, xyz)
    for atom in result.GetAtoms():
        atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    Chem.AssignStereochemistryFrom3D(result, replaceExistingTags=True)
    Chem.AssignStereochemistry(result, cleanIt=True, force=True)
    return [
        {"atom_index_zero_based": atom.GetIdx(), "cip": atom.GetProp("_CIPCode")}
        for atom in result.GetAtoms()
        if atom.HasProp("_CIPCode")
    ]


def normalized_source_bond_type(bond: Chem.Bond) -> str:
    return "AROMATIC" if bond.GetIsAromatic() else str(bond.GetBondType()).upper()


def normalized_mol2_bond_type(bond_type: str) -> str:
    return {
        "1": "SINGLE",
        "2": "DOUBLE",
        "3": "TRIPLE",
        "ar": "AROMATIC",
    }.get(bond_type.lower(), bond_type.upper())


def audit_amber_openmm(prmtop: Path, inpcrd: Path) -> tuple[float, float, str]:
    import openmm
    from openmm import unit
    from openmm.app import AmberInpcrdFile, AmberPrmtopFile, NoCutoff

    topology = AmberPrmtopFile(str(prmtop))
    coordinates = AmberInpcrdFile(str(inpcrd))
    system = topology.createSystem(nonbondedMethod=NoCutoff, constraints=None)
    nonbonded = next(
        force
        for force in system.getForces()
        if isinstance(force, openmm.NonbondedForce)
    )
    total_charge = sum(
        nonbonded.getParticleParameters(index)[0].value_in_unit(unit.elementary_charge)
        for index in range(system.getNumParticles())
    )
    integrator = openmm.VerletIntegrator(1.0 * unit.femtoseconds)
    platform = openmm.Platform.getPlatformByName("CPU")
    context = openmm.Context(system, integrator, platform, {"Threads": "1"})
    context.setPositions(coordinates.positions)
    energy = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
        unit.kilocalories_per_mole
    )
    del context, integrator
    return float(total_charge), float(energy), openmm.__version__


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sdf", type=Path, default=DEFAULT_SDF)
    parser.add_argument("--mol2", type=Path, default=DEFAULT_MOL2)
    parser.add_argument("--frcmod", type=Path, default=DEFAULT_FRCMOD)
    parser.add_argument("--prmtop", type=Path, default=DEFAULT_PRMTOP)
    parser.add_argument("--inpcrd", type=Path, default=DEFAULT_INPCRD)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    paths = [args.sdf, args.mol2, args.frcmod, args.prmtop, args.inpcrd]
    for path in paths:
        if not path.is_file():
            raise FileNotFoundError(path)

    source = Chem.MolFromMolFile(str(args.sdf), removeHs=False, sanitize=True)
    if source is None:
        raise ValueError(f"Could not read source PaG carboxylate SDF: {args.sdf}")
    source_cip = cip_labels_from_tags(Chem.Mol(source))
    source_formula = rdMolDescriptors.CalcMolFormula(source)
    source_charge = Chem.GetFormalCharge(source)
    if source.GetNumAtoms() != 39 or source_charge != -1:
        raise ValueError("Expected a 39-atom, singly deprotonated PaG conformer.")

    atoms, mol2_bonds = read_mol2(args.mol2)
    if len(atoms) != source.GetNumAtoms() or len(mol2_bonds) != source.GetNumBonds():
        raise ValueError("GAFF2 Mol2 atom/bond counts differ from the input SDF.")
    if [atom["element"] for atom in atoms] != [
        atom.GetSymbol() for atom in source.GetAtoms()
    ]:
        raise ValueError("GAFF2 Mol2 atom element order differs from the source SDF.")

    source_bonds = {
        tuple(sorted((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))): normalized_source_bond_type(bond)
        for bond in source.GetBonds()
    }
    mol2_bond_map = {
        tuple(sorted((left, right))): normalized_mol2_bond_type(kind)
        for left, right, kind in mol2_bonds
    }
    bond_differences = [
        {
            "atom_names": [atoms[left]["name"], atoms[right]["name"]],
            "source_sdf_type": source_bonds.get((left, right)),
            "gaff2_mol2_type": mol2_bond_map.get((left, right)),
        }
        for left, right in sorted(set(source_bonds) | set(mol2_bond_map))
        if source_bonds.get((left, right)) != mol2_bond_map.get((left, right))
    ]
    if set(source_bonds) != set(mol2_bond_map):
        raise ValueError("GAFF2 Mol2 connectivity differs from the source SDF.")

    # The source SDF localizes the glucuronate charge and C=O bond on one
    # oxygen. GAFF2 represents this carboxylate with two single C-O bonds.
    carboxylates = []
    for carbon in source.GetAtoms():
        if carbon.GetAtomicNum() != 6:
            continue
        oxygens = [neighbor for neighbor in carbon.GetNeighbors() if neighbor.GetAtomicNum() == 8]
        if len(oxygens) != 2:
            continue
        oxygen_charges = sum(oxygen.GetFormalCharge() for oxygen in oxygens)
        bond_types = [
            source.GetBondBetweenAtoms(carbon.GetIdx(), oxygen.GetIdx()).GetBondType()
            for oxygen in oxygens
        ]
        if oxygen_charges == -1 and sorted(str(kind) for kind in bond_types) == ["DOUBLE", "SINGLE"]:
            carboxylates.append((carbon, oxygens))
    if len(carboxylates) != 1:
        raise ValueError(f"Expected one localized glucuronate carboxylate; found {len(carboxylates)}.")
    carboxyl_carbon, carboxylate_oxygens = carboxylates[0]
    double_bond_oxygen = next(
        oxygen
        for oxygen in carboxylate_oxygens
        if source.GetBondBetweenAtoms(carboxyl_carbon.GetIdx(), oxygen.GetIdx()).GetBondType()
        == Chem.BondType.DOUBLE
    )
    expected_resonance_difference = [
        {
            "atom_names": [
                atoms[double_bond_oxygen.GetIdx()]["name"],
                atoms[carboxyl_carbon.GetIdx()]["name"],
            ],
            "source_sdf_type": "DOUBLE",
            "gaff2_mol2_type": "SINGLE",
        }
    ]
    output_carboxylate_types = [
        mol2_bond_map[
            tuple(sorted((carboxyl_carbon.GetIdx(), oxygen.GetIdx())))
        ]
        for oxygen in carboxylate_oxygens
    ]
    if output_carboxylate_types != ["SINGLE", "SINGLE"]:
        raise ValueError(
            f"GAFF2 carboxylate bonds are not both single bonds: {output_carboxylate_types}"
        )
    if bond_differences != expected_resonance_difference:
        raise ValueError(f"Unexpected bond-order differences: {bond_differences}")

    output_cip = cip_labels_from_coordinates(
        source, [atom["xyz"] for atom in atoms]
    )
    if output_cip != source_cip:
        raise ValueError(
            f"Stereochemistry differs between source and GAFF2 geometry: "
            f"{source_cip} versus {output_cip}"
        )

    mol2_charge = sum(atom["partial_charge"] for atom in atoms)
    topology_charge, initial_energy, openmm_version = audit_amber_openmm(
        args.prmtop, args.inpcrd
    )
    if not math.isclose(mol2_charge, -1.0, abs_tol=2e-4):
        raise ValueError(f"GAFF2 Mol2 partial charges sum to {mol2_charge}, expected -1.")
    if not math.isclose(topology_charge, -1.0, abs_tol=2e-4):
        raise ValueError(f"Amber topology charges sum to {topology_charge}, expected -1.")
    if not math.isfinite(initial_energy):
        raise ValueError("OpenMM initial energy is not finite.")

    report = {
        "record_type": "pag_gaff2_ligand_parameter_file_audit",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "source": {
            "pubchem_cid": 83944,
            "prepared_conformer": args.sdf.name,
            "formula": source_formula,
            "formal_charge": source_charge,
            "stereocenters_cip": source_cip,
        },
        "parameterization": {
            "method": "GAFF2 / AM1-BCC",
            "toolchain": "AmberTools 26.0",
            "mol2_atom_count": len(atoms),
            "mol2_bond_count": len(mol2_bonds),
            "atom_element_order_matches_source": True,
            "connectivity_matches_source": True,
            "output_geometry_cip_matches_source": True,
            "output_geometry_cip": output_cip,
            "mol2_partial_charge_sum_e": mol2_charge,
            "expected_carboxylate_resonance_bond_difference": bond_differences,
        },
        "openmm_file_check": {
            "openmm_version": openmm_version,
            "prmtop_particle_count": len(atoms),
            "topology_charge_sum_e": topology_charge,
            "initial_potential_energy_kcal_mol": initial_energy,
            "energy_is_finite": math.isfinite(initial_energy),
            "platform": "CPU",
            "scope": "ligand-only Amber parameter and coordinate loading check",
        },
        "files": {
            path.name: {"path": str(path), "sha256": sha256(path)}
            for path in paths
        },
        "limitations": [
            "No NNC2215 conjugate or protein is included.",
            "No host–guest pose, binding energy, or affinity was calculated.",
            "GAFF2/AM1-BCC is a candidate parameterization, not a validated transfer to the NNC2215 conjugate.",
            "The carboxylate uses a delocalized GAFF2 single-bond representation where the source SDF uses one localized C=O bond.",
        ],
        "clinical_use": False,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
