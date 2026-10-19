#!/usr/bin/env python3
"""Run a short explicit-water PaG compatibility pilot, not a binding study."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import platform as host_platform
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import openmm
from openff.toolkit import Molecule
from openmm import app, unit
from openmmforcefields.generators import SMIRNOFFTemplateGenerator


ROOT = Path(__file__).resolve().parent
INPUT = ROOT / "ligands/pag-pubchem-cid83944-carboxylate-3d.sdf"
OUTPUT = ROOT / "data/pag-solvation-preflight"
LOCK = ROOT / "pag-chemistry-environment-macos-arm64.lock.txt"
FORCE_FIELD = "openff-2.3.0.offxml"
PLATFORM_NAME = "OpenCL"
PLATFORM_PROPERTIES = {"Precision": "single", "DeviceIndex": "0"}
TEMPERATURE_K = 298.15
PRESSURE_ATM = 1.0
TIMESTEP_FS = 2.0
EQUILIBRATION_STEPS = 5_000  # 10 ps
PRODUCTION_STEPS = 50_000  # 100 ps
REPORT_INTERVAL = 1_000
REPLICATE_SEEDS = (2215001, 2215002)
PADDING_NM = 1.0
SALT_CONCENTRATION_M = 0.15


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def system_is_finite(simulation: app.Simulation) -> dict[str, float]:
    state = simulation.context.getState(getEnergy=True)
    potential = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    kinetic = state.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)
    degrees_of_freedom = (
        3 * simulation.system.getNumParticles()
        - simulation.system.getNumConstraints()
        - 3  # CMMotionRemover removes center-of-mass translation.
    )
    temperature = state.getKineticEnergy() / (
        0.5 * degrees_of_freedom * unit.MOLAR_GAS_CONSTANT_R
    )
    temperature_k = temperature.value_in_unit(unit.kelvin)
    values = {
        "potential_energy_kj_mol": float(potential),
        "kinetic_energy_kj_mol": float(kinetic),
        "temperature_K": float(temperature_k),
    }
    if not all(math.isfinite(value) for value in values.values()):
        raise RuntimeError(f"Non-finite state diagnostic: {values}")
    return values


def inspect_report(path: Path) -> dict[str, object]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise RuntimeError(f"No thermodynamic samples in {path}")
    numeric = []
    for row in rows:
        converted = {}
        for key, value in row.items():
            if key in {"Step", "Time (ps)"}:
                continue
            if value is None or value == "":
                continue
            try:
                number = float(value)
            except ValueError:
                continue
            if not math.isfinite(number):
                raise RuntimeError(f"Non-finite {key} in {path}: {value}")
            converted[key] = number
        numeric.append(converted)
    temps = [row.get("Temperature (K)") for row in numeric if "Temperature (K)" in row]
    energies = [row.get("Potential Energy (kJ/mole)") for row in numeric if "Potential Energy (kJ/mole)" in row]
    volumes = [
        row.get("Box Volume (nm^3)")
        for row in numeric
        if "Box Volume (nm^3)" in row
    ]
    return {
        "samples": len(rows),
        "temperature_mean_K": float(np.mean(temps)) if temps else None,
        "temperature_min_K": float(np.min(temps)) if temps else None,
        "temperature_max_K": float(np.max(temps)) if temps else None,
        "potential_energy_min_kj_mol": float(np.min(energies)) if energies else None,
        "potential_energy_max_kj_mol": float(np.max(energies)) if energies else None,
        "volume_min_nm3": float(np.min(volumes)) if volumes else None,
        "volume_max_nm3": float(np.max(volumes)) if volumes else None,
    }


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    run_id = "DHO-PAG-SOLVATION-" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")

    molecule = Molecule.from_file(INPUT)
    if molecule.n_atoms != 39 or molecule.total_charge.m_as("elementary_charge") != -1:
        raise ValueError("Expected the checked 39-atom, net -1 PaG carboxylate model.")
    molecule.assign_partial_charges(partial_charge_method="am1bcc")

    topology = molecule.to_topology().to_openmm()
    positions = molecule.conformers[0].to_openmm()
    modeller = app.Modeller(topology, positions)
    smirnoff = SMIRNOFFTemplateGenerator(
        molecules=[molecule], forcefield=FORCE_FIELD
    )
    forcefield = app.ForceField("amber14-all.xml", "amber14/tip3p.xml")
    forcefield.registerTemplateGenerator(smirnoff.generator)
    modeller.addSolvent(
        forcefield,
        model="tip3p",
        padding=PADDING_NM * unit.nanometer,
        ionicStrength=SALT_CONCENTRATION_M * unit.molar,
        neutralize=True,
    )
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=app.HBonds,
        hydrogenMass=4.0 * unit.amu,
    )

    box_dimensions_nm = modeller.topology.getUnitCellDimensions().value_in_unit(
        unit.nanometer
    )
    box_volume_nm3 = float(np.prod(box_dimensions_nm))
    box_volume_l = box_volume_nm3 * 1e-24
    box_implied_pag_concentration_m = 1.0 / (
        unit.AVOGADRO_CONSTANT_NA.value_in_unit(unit.mole**-1) * box_volume_l
    )
    residues: dict[str, int] = {}
    for residue in modeller.topology.residues():
        residues[residue.name] = residues.get(residue.name, 0) + 1

    initial_pdb = OUTPUT / f"{run_id}-initial.pdb"
    with initial_pdb.open("w") as handle:
        app.PDBFile.writeFile(modeller.topology, modeller.positions, handle)
    system_xml = OUTPUT / f"{run_id}-system.xml"
    system_xml.write_text(openmm.XmlSerializer.serialize(system))

    replicate_results = []
    platform = openmm.Platform.getPlatformByName(PLATFORM_NAME)
    for replicate_index, seed in enumerate(REPLICATE_SEEDS, start=1):
        rep_name = f"replicate-{replicate_index:02d}"
        log_path = OUTPUT / f"{run_id}-{rep_name}.csv"
        final_path = OUTPUT / f"{run_id}-{rep_name}-final.pdb"
        replicate_system = openmm.XmlSerializer.deserialize(
            openmm.XmlSerializer.serialize(system)
        )
        barostat = openmm.MonteCarloBarostat(
            PRESSURE_ATM * unit.atmosphere,
            TEMPERATURE_K * unit.kelvin,
            25,
        )
        barostat.setRandomNumberSeed(seed + 100)
        replicate_system.addForce(barostat)
        integrator = openmm.LangevinMiddleIntegrator(
            TEMPERATURE_K * unit.kelvin,
            1.0 / unit.picosecond,
            TIMESTEP_FS * unit.femtosecond,
        )
        integrator.setRandomNumberSeed(seed)
        simulation = app.Simulation(
            modeller.topology,
            replicate_system,
            integrator,
            platform,
            PLATFORM_PROPERTIES,
        )
        simulation.context.setPositions(modeller.positions)
        simulation.minimizeEnergy(maxIterations=1_000)
        simulation.context.setVelocitiesToTemperature(TEMPERATURE_K * unit.kelvin, seed)
        post_minimization_initial_state = system_is_finite(simulation)
        simulation.reporters.append(
            app.StateDataReporter(
                str(log_path),
                REPORT_INTERVAL,
                step=True,
                time=True,
                potentialEnergy=True,
                kineticEnergy=True,
                temperature=True,
                volume=True,
                separator=",",
            )
        )
        simulation.step(EQUILIBRATION_STEPS)
        simulation.step(PRODUCTION_STEPS)
        final_diagnostics = system_is_finite(simulation)
        report = inspect_report(log_path)
        with final_path.open("w") as handle:
            app.PDBFile.writeFile(
                modeller.topology,
                simulation.context.getState(getPositions=True).getPositions(),
                handle,
            )
        replicate_results.append(
            {
                "replicate": rep_name,
                "integrator_seed": seed,
                "barostat_seed": seed + 100,
                "equilibration_steps": EQUILIBRATION_STEPS,
                "production_steps": PRODUCTION_STEPS,
                "production_duration_ps": PRODUCTION_STEPS * TIMESTEP_FS / 1000.0,
                "post_minimization_initial_state": post_minimization_initial_state,
                "final_state": final_diagnostics,
                "thermodynamic_samples": report,
                "log_file": log_path.name,
                "log_sha256": sha256(log_path),
                "final_structure_file": final_path.name,
                "final_structure_sha256": sha256(final_path),
            }
        )
        del simulation
        del integrator

    report = {
        "schema": "dia-health-one.pag-solvation-preflight.v1",
        "run": {
            "id": run_id,
            "started_at_utc": run_id.removeprefix("DHO-PAG-SOLVATION-"),
            "completed_at_utc": datetime.now(timezone.utc).isoformat(),
            "status": "completed",
            "purpose": "PaG small-molecule explicit-water parameter compatibility smoke test; not a binding or affinity study",
            "clinical_use": False,
            "script_sha256": sha256(Path(__file__)),
            "input_sha256": sha256(INPUT),
            "environment_lock_sha256": sha256(LOCK),
        },
        "chemical_model": {
            "source": "PubChem CID 83944; neutral 2D identity record",
            "modeled_state": "PaG glucuronate carboxylate, net charge -1",
            "protonation": "explicit modeling assumption near physiological pH; not an experimental pKa determination",
            "stereochemistry": "five CIP centers retained from the PubChem source record",
            "conformer": "RDKit ETKDGv3 seed 2215 and MMFF94s minimization",
            "partial_charge_method": "AM1-BCC assigned by AmberTools 26.0 through OpenFF Toolkit",
            "ligand_force_field": FORCE_FIELD,
            "water_force_field": "Amber14/TIP3P",
            "mixed_force_field_note": "OpenFF ligand parameters combined with Amber14/TIP3P solvent through the SMIRNOFF template generator; this is a small-molecule compatibility setup, not a validated NNC2215 parameterization.",
        },
        "system": {
            "atom_count": modeller.topology.getNumAtoms(),
            "residue_counts": residues,
            "box_dimensions_nm_initial": [float(x) for x in box_dimensions_nm],
            "box_volume_nm3_initial": box_volume_nm3,
            "one_ligand_box_implied_concentration_M": box_implied_pag_concentration_m,
            "salt_added_M": SALT_CONCENTRATION_M,
            "neutralize": True,
            "nonbonded_method": "PME",
            "nonbonded_cutoff_nm": 1.0,
            "constraints": "HBonds",
            "hydrogen_mass_amu": 4.0,
            "clinical_concentration_interpretation": "The one-ligand periodic-box concentration is a finite-size consequence and is not equivalent to the 70 micromol/L experimental assay or plasma exposure.",
        },
        "protocol": {
            "ensemble": "NPT",
            "temperature_K": TEMPERATURE_K,
            "pressure_atm": PRESSURE_ATM,
            "timestep_fs": TIMESTEP_FS,
            "equilibration_ps": EQUILIBRATION_STEPS * TIMESTEP_FS / 1000.0,
            "production_ps_per_replicate": PRODUCTION_STEPS * TIMESTEP_FS / 1000.0,
            "replicate_count": len(REPLICATE_SEEDS),
            "independent_velocity_seeds": list(REPLICATE_SEEDS),
            "platform": PLATFORM_NAME,
            "platform_properties": PLATFORM_PROPERTIES,
            "host_os": host_platform.platform(),
            "host_architecture": host_platform.machine(),
            "host_processor": "Apple M4; OpenCL device 0",
            "limitations": [
                "Only PaG in water and salt was simulated; no NNC2215 or macrocycle is present.",
                "Two 100 ps replicas assess setup and gross stability only; they do not establish conformational equilibrium.",
                "No ligand binding, receptor affinity, or published PaG MD result is reproduced.",
                "Water-mixed ligand parameters do not validate transferability to the macrocycle/protein system.",
                "This Apple OpenCL run is not an NVIDIA benchmark.",
            ],
        },
        "replicates": replicate_results,
        "artifacts": {
            "initial_system_pdb": initial_pdb.name,
            "initial_system_pdb_sha256": sha256(initial_pdb),
            "openmm_system_xml": system_xml.name,
            "openmm_system_xml_sha256": sha256(system_xml),
        },
    }
    output_json = OUTPUT / f"{run_id}.json"
    output_json.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
