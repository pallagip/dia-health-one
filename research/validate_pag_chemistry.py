#!/usr/bin/env python3
"""Validate a bounded PaG small-molecule parameterization candidate.

This checks a PubChem chemical-identity record and a separately constructed
carboxylate state. It does not parameterize or simulate the NNC2215 conjugate.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors


EXPECTED_FORMULA = "C14H17NO8"
EXPECTED_SOURCE_CHARGE = 0
MODELED_CARBOXYLATE_CHARGE = -1


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def cip_labels(molecule: Chem.Mol) -> list[tuple[int, str]]:
    Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)
    return [
        (atom.GetIdx(), atom.GetProp("_CIPCode"))
        for atom in molecule.GetAtoms()
        if atom.HasProp("_CIPCode")
    ]


def deprotonate_glucuronate(source: Chem.Mol) -> Chem.Mol:
    """Remove the explicit acidic H and assign -1 to the glucuronate oxygen."""
    candidates: list[tuple[int, int]] = []
    for atom in source.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        carbonyl_o = [
            neighbor
            for neighbor in atom.GetNeighbors()
            if neighbor.GetAtomicNum() == 8
            and source.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx()).GetBondType()
            == Chem.BondType.DOUBLE
        ]
        hydroxyl_o = [
            neighbor
            for neighbor in atom.GetNeighbors()
            if neighbor.GetAtomicNum() == 8
            and source.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx()).GetBondType()
            == Chem.BondType.SINGLE
            and any(h.GetAtomicNum() == 1 for h in neighbor.GetNeighbors())
        ]
        if carbonyl_o and hydroxyl_o:
            for oxygen in hydroxyl_o:
                hydrogen = next(
                    h for h in oxygen.GetNeighbors() if h.GetAtomicNum() == 1
                )
                candidates.append((oxygen.GetIdx(), hydrogen.GetIdx()))

    if len(candidates) != 1:
        raise ValueError(
            f"Expected one explicit carboxylic-acid OH; found {len(candidates)}."
        )

    oxygen_idx, hydrogen_idx = candidates[0]
    editable = Chem.RWMol(source)
    editable.GetAtomWithIdx(oxygen_idx).SetFormalCharge(-1)
    editable.RemoveAtom(hydrogen_idx)
    result = editable.GetMol()
    Chem.SanitizeMol(result)
    if Chem.GetFormalCharge(result) != MODELED_CARBOXYLATE_CHARGE:
        raise ValueError("Constructed glucuronate did not have net charge -1.")
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        type=Path,
        default=Path(__file__).parent / "ligands/pag-pubchem-cid83944-2d.sdf",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=Path(__file__).parent / "ligands"
    )
    parser.add_argument("--force-field", default="openff-2.3.0.offxml")
    args = parser.parse_args()

    source = Chem.MolFromMolFile(str(args.input), removeHs=False, sanitize=True)
    if source is None:
        raise ValueError(f"Could not parse SDF: {args.input}")
    source_formula = rdMolDescriptors.CalcMolFormula(source)
    source_charge = Chem.GetFormalCharge(source)
    if source_formula != EXPECTED_FORMULA:
        raise ValueError(f"Unexpected PubChem source formula: {source_formula}")
    if source_charge != EXPECTED_SOURCE_CHARGE:
        raise ValueError(f"Expected neutral PubChem source form, got {source_charge}.")

    source_stereo = cip_labels(Chem.Mol(source))
    modeled = deprotonate_glucuronate(source)
    modeled_formula = rdMolDescriptors.CalcMolFormula(modeled)
    modeled_stereo = cip_labels(Chem.Mol(modeled))
    if source_stereo != modeled_stereo:
        raise ValueError("CIP stereochemistry changed during carboxylate construction.")

    modeled = Chem.RemoveAllHs(modeled)
    modeled = Chem.AddHs(modeled)
    params = AllChem.ETKDGv3()
    params.randomSeed = 2215
    params.useSmallRingTorsions = True
    if AllChem.EmbedMolecule(modeled, params) != 0:
        raise RuntimeError("RDKit ETKDGv3 failed to generate a 3D conformer.")
    mmff = AllChem.MMFFGetMoleculeProperties(modeled, mmffVariant="MMFF94s")
    if mmff is None:
        raise RuntimeError("RDKit MMFF94s parameters unavailable for the PaG model.")
    status = AllChem.MMFFOptimizeMolecule(
        modeled, mmffVariant="MMFF94s", maxIters=2000
    )
    if status != 0:
        raise RuntimeError(f"RDKit MMFF94s geometry optimization failed: status={status}")

    # Imports are intentionally delayed so the identity and stereochemistry
    # checks can be run independently from the force-field toolchain.
    from openff.toolkit import ForceField, Molecule

    off_molecule = Molecule.from_rdkit(
        modeled, allow_undefined_stereo=False, hydrogens_are_explicit=True
    )
    off_molecule.assign_partial_charges(partial_charge_method="am1bcc")
    charge_sum = float(off_molecule.partial_charges.m_as("elementary_charge").sum())
    if abs(charge_sum - MODELED_CARBOXYLATE_CHARGE) > 1e-5:
        raise ValueError(f"AM1-BCC charge sum is {charge_sum}, expected -1.")

    force_field = ForceField(args.force_field)
    openmm_system = force_field.create_openmm_system(off_molecule.to_topology())
    if openmm_system.getNumParticles() != off_molecule.n_atoms:
        raise ValueError("OpenMM system particle count differs from PaG atom count.")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    output_sdf = args.output_dir / "pag-pubchem-cid83944-carboxylate-3d.sdf"
    writer = Chem.SDWriter(str(output_sdf))
    modeled.SetProp("_Name", "PaG carboxylate; modeled from PubChem CID 83944")
    modeled.SetProp("protonation_assumption", "glucuronate carboxylate; net charge -1")
    modeled.SetProp("coordinate_method", "RDKit ETKDGv3 seed 2215; MMFF94s minimization")
    writer.write(modeled)
    writer.close()

    import openff.toolkit
    import openmm
    from rdkit import rdBase

    report = {
        "record_type": "pag_small_molecule_parameterization_preflight",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "source": {
            "pubchem_cid": 83944,
            "url": "https://pubchem.ncbi.nlm.nih.gov/compound/83944",
            "file": args.input.name,
            "sha256": sha256(args.input),
            "formula": source_formula,
            "formal_charge": source_charge,
            "coordinates": "2D source coordinates; not used as a simulation conformer",
        },
        "modeled_ligand": {
            "name": "paracetamol glucuronide carboxylate",
            "formula": modeled_formula,
            "formal_charge": MODELED_CARBOXYLATE_CHARGE,
            "protonation_assumption": "deprotonated glucuronic-acid carboxyl group near physiological pH",
            "stereochemistry_cip_before_after": source_stereo,
            "coordinates": "RDKit ETKDGv3 seed 2215, MMFF94s minimized",
            "output_sdf": output_sdf.name,
            "output_sha256": sha256(output_sdf),
            "partial_charge_method": "AM1-BCC",
            "partial_charge_sum_e": charge_sum,
            "force_field": args.force_field,
            "force_field_parameterization": "OpenFF toolkit small-molecule system creation completed",
            "openmm_particles": openmm_system.getNumParticles(),
            "not_checked": [
                "NNC2215 macrocycle/linker/glucoside parameters",
                "protein–ligand compatibility",
                "PaG binding affinity",
                "published 500 ns receptor-surrogate simulation",
                "clinical exposure or response",
            ],
        },
        "software": {
            "rdkit": rdBase.rdkitVersion,
            "openff_toolkit": openff.toolkit.__version__,
            "openmm": openmm.__version__,
        },
        "clinical_use": False,
    }
    report_path = args.output_dir / "pag-chemistry-preflight.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
