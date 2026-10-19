#!/usr/bin/env python3
"""Replay the PaG-only solvent setup on an NVIDIA GPU as a portability preflight.

This reuses the serialized OpenMM system and initial coordinates from the local
PaG-in-water pilot. It does not contain NNC2215, a receptor, or a binding
observable, so it is not a reproduction of the published PaG-binding study.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform as host_platform
import statistics
import subprocess
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import openmm
from openmm import app, unit


TEMPERATURE_K = 298.15
PRESSURE_ATM = 1.0
TIMESTEP_FS = 2.0
FRICTION_PER_PS = 1.0
EQUILIBRATION_STEPS = 5_000  # 10 ps
PRODUCTION_STEPS = 50_000  # 100 ps
REPORT_INTERVAL = 1_000
REPLICATE_SEEDS = (2215001, 2215002)
PLATFORM_NAME = "OpenCL"
REQUESTED_PLATFORM_PROPERTIES = {"DeviceIndex": "0", "Precision": "single"}
EXPECTED_OPENMM_VERSION = "8.6.1"
EXPECTED_PDB_SHA256 = "1b5bf3cff8f6c433069d621636e89150000e5fc3dd7a3c5d76e15d9d8436fd4b"
EXPECTED_SYSTEM_SHA256 = "1c5d1de602db30728bb342cc2607c5b3720c93421fe7bace6a103eb6a632082f"


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def stage_input(source: str, destination: Path) -> dict[str, str]:
    """Copy a local file or download a public URL into the run directory."""
    if source.startswith(("https://", "http://")):
        request = urllib.request.Request(
            source, headers={"User-Agent": "DiaHealthOne-research/1.0"}
        )
        with urllib.request.urlopen(request, timeout=60) as response:
            payload = response.read()
        destination.write_bytes(payload)
        return {"source": source, "sha256": sha256_bytes(payload)}
    source_path = Path(source).resolve()
    if not source_path.is_file():
        raise FileNotFoundError(source_path)
    payload = source_path.read_bytes()
    destination.write_bytes(payload)
    return {"source": str(source_path), "sha256": sha256_bytes(payload)}


def gpu_inventory() -> dict[str, str]:
    result = subprocess.run(
        [
            "nvidia-smi",
            "--query-gpu=name,uuid,driver_version,memory.total",
            "--format=csv,noheader",
        ],
        capture_output=True,
        text=True,
        check=True,
    )
    first = result.stdout.strip().splitlines()[0]
    model, uuid, driver, memory = (part.strip() for part in first.split(",", 3))
    return {
        "model": model,
        "uuid": uuid,
        "driver_version": driver,
        "memory_total": memory,
    }


def finite_diagnostics(simulation: app.Simulation) -> dict[str, float]:
    state = simulation.context.getState(getEnergy=True)
    potential = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    kinetic = state.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)
    dof = (
        3 * simulation.system.getNumParticles()
        - simulation.system.getNumConstraints()
        - 3  # CMMotionRemover removes center-of-mass translation.
    )
    temperature = state.getKineticEnergy() / (
        0.5 * dof * unit.MOLAR_GAS_CONSTANT_R
    )
    values = {
        "potential_energy_kj_mol": float(potential),
        "kinetic_energy_kj_mol": float(kinetic),
        "temperature_K": float(temperature.value_in_unit(unit.kelvin)),
    }
    if not all(math.isfinite(value) for value in values.values()):
        raise RuntimeError(f"Non-finite state diagnostics: {values}")
    return values


def report_summary(path: Path) -> dict[str, float | int]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise RuntimeError(f"No thermodynamic samples in {path}")
    columns = (
        "Temperature (K)",
        "Potential Energy (kJ/mole)",
        "Box Volume (nm^3)",
    )
    values: dict[str, list[float]] = {column: [] for column in columns}
    for row in rows:
        for column in columns:
            number = float(row[column])
            if not math.isfinite(number):
                raise RuntimeError(f"Non-finite {column} in {path}: {number}")
            values[column].append(number)
    return {
        "samples": len(rows),
        "temperature_mean_K": statistics.mean(values["Temperature (K)"]),
        "temperature_min_K": min(values["Temperature (K)"]),
        "temperature_max_K": max(values["Temperature (K)"]),
        "potential_energy_min_kj_mol": min(
            values["Potential Energy (kJ/mole)"]
        ),
        "potential_energy_max_kj_mol": max(
            values["Potential Energy (kJ/mole)"]
        ),
        "volume_min_nm3": min(values["Box Volume (nm^3)"]),
        "volume_max_nm3": max(values["Box Volume (nm^3)"]),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-pdb", required=True, help="Local path or public URL")
    parser.add_argument("--system-xml", required=True, help="Local path or public URL")
    parser.add_argument("--output-dir", type=Path, default=Path("pag-t4-preflight"))
    args = parser.parse_args()

    started_at = datetime.now(timezone.utc)
    run_id = "DHO-PAG-GPU-" + started_at.strftime("%Y%m%dT%H%M%SZ")
    output_dir = args.output_dir.resolve()
    run_dir = output_dir / run_id
    run_dir.mkdir(parents=True, exist_ok=False)

    pdb_path = run_dir / "initial.pdb"
    xml_path = run_dir / "system.xml"
    pdb_source = stage_input(args.input_pdb, pdb_path)
    xml_source = stage_input(args.system_xml, xml_path)
    if pdb_source["sha256"] != EXPECTED_PDB_SHA256:
        raise ValueError(
            "Initial PDB checksum mismatch: "
            f"{pdb_source['sha256']} != {EXPECTED_PDB_SHA256}"
        )
    if xml_source["sha256"] != EXPECTED_SYSTEM_SHA256:
        raise ValueError(
            "System XML checksum mismatch: "
            f"{xml_source['sha256']} != {EXPECTED_SYSTEM_SHA256}"
        )
    if openmm.__version__ != EXPECTED_OPENMM_VERSION:
        raise RuntimeError(
            f"Expected OpenMM {EXPECTED_OPENMM_VERSION}, found {openmm.__version__}"
        )

    pdb = app.PDBFile(str(pdb_path))
    system = openmm.XmlSerializer.deserialize(xml_path.read_text(encoding="utf-8"))
    if system.getNumParticles() != pdb.topology.getNumAtoms():
        raise ValueError(
            "System/PDB particle mismatch: "
            f"{system.getNumParticles()} vs {pdb.topology.getNumAtoms()}"
        )
    if pdb.topology.getPeriodicBoxVectors() is None:
        raise ValueError("Input PDB does not contain periodic box vectors")

    gpu = gpu_inventory()
    platform = openmm.Platform.getPlatformByName(PLATFORM_NAME)
    supported = set(platform.getPropertyNames())
    properties = {
        name: value
        for name, value in REQUESTED_PLATFORM_PROPERTIES.items()
        if name in supported
    }
    measurements = []
    for index, seed in enumerate(REPLICATE_SEEDS, start=1):
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
            FRICTION_PER_PS / unit.picosecond,
            TIMESTEP_FS * unit.femtosecond,
        )
        integrator.setRandomNumberSeed(seed)
        simulation = app.Simulation(
            pdb.topology, replicate_system, integrator, platform, properties
        )
        device_name = platform.getPropertyValue(simulation.context, "DeviceName")
        if gpu["model"].lower() not in device_name.lower():
            raise RuntimeError(
                f"GPU mismatch: nvidia-smi={gpu['model']!r}, OpenMM={device_name!r}"
            )
        simulation.context.setPositions(pdb.positions)
        simulation.minimizeEnergy(maxIterations=1_000)
        simulation.context.setVelocitiesToTemperature(TEMPERATURE_K * unit.kelvin, seed)
        initial_state = finite_diagnostics(simulation)

        replicate_name = f"replicate-{index:02d}"
        csv_path = run_dir / f"{replicate_name}.csv"
        pdb_out = run_dir / f"{replicate_name}-final.pdb"
        simulation.reporters.append(
            app.StateDataReporter(
                str(csv_path),
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
        final_state = finite_diagnostics(simulation)
        summary = report_summary(csv_path)
        with pdb_out.open("w", encoding="utf-8") as stream:
            app.PDBFile.writeFile(
                pdb.topology,
                simulation.context.getState(getPositions=True).getPositions(),
                stream,
                keepIds=True,
            )
        measurements.append(
            {
                "replicate": replicate_name,
                "integrator_seed": seed,
                "barostat_seed": seed + 100,
                "equilibration_steps": EQUILIBRATION_STEPS,
                "production_steps": PRODUCTION_STEPS,
                "production_ps": PRODUCTION_STEPS * TIMESTEP_FS / 1000.0,
                "device_name": device_name,
                "initial_state_after_minimization": initial_state,
                "final_state": final_state,
                "thermodynamic_summary": summary,
                "log_file": csv_path.name,
                "log_sha256": sha256_file(csv_path),
                "final_structure_file": pdb_out.name,
                "final_structure_sha256": sha256_file(pdb_out),
            }
        )
        del simulation, integrator

    completed_at = datetime.now(timezone.utc)
    box = pdb.topology.getPeriodicBoxVectors()
    box_nm = [
        [float(x) for x in vector.value_in_unit(unit.nanometer)] for vector in box
    ]
    a, b, c = box_nm
    box_volume_nm3 = abs(
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )
    box_volume_l = box_volume_nm3 * 1e-24
    implied_pag_concentration_m = 1.0 / (
        unit.AVOGADRO_CONSTANT_NA.value_in_unit(unit.mole**-1) * box_volume_l
    )
    residues: dict[str, int] = {}
    for residue in pdb.topology.residues():
        residues[residue.name] = residues.get(residue.name, 0) + 1
    result = {
        "schema": "dia-health-one.pag-gpu-portability-preflight.v1",
        "run": {
            "id": run_id,
            "started_at_utc": started_at.isoformat(),
            "completed_at_utc": completed_at.isoformat(),
            "status": "completed",
            "purpose": "PaG-only water-system portability preflight on an NVIDIA GPU; not binding MD",
            "nvidia_gpu_used": True,
            "clinical_use": False,
        },
        "inputs": {
            "initial_pdb": pdb_source,
            "serialized_openmm_system": xml_source,
            "coordinate_precision_note": "The reused initial PDB is a decimal-coordinate serialization of the earlier Apple M4 system builder output.",
        },
        "environment": {
            "openmm": openmm.__version__,
            "python": host_platform.python_version(),
            "os": host_platform.platform(),
            "architecture": host_platform.machine(),
            "gpu": gpu,
        },
        "system": {
            "atom_count": system.getNumParticles(),
            "residue_counts": residues,
            "box_vectors_nm": box_nm,
            "box_volume_nm3": box_volume_nm3,
            "one_ligand_box_implied_concentration_M": implied_pag_concentration_m,
            "platform": PLATFORM_NAME,
            "platform_properties": properties,
            "force_field_source": "Serialized OpenFF 2.3.0 PaG + Amber14/TIP3P pilot system from research/data/pag-solvation-preflight.",
        },
        "protocol": {
            "ensemble": "NPT",
            "temperature_K": TEMPERATURE_K,
            "pressure_atm": PRESSURE_ATM,
            "timestep_fs": TIMESTEP_FS,
            "friction_per_ps": FRICTION_PER_PS,
            "equilibration_ps": EQUILIBRATION_STEPS * TIMESTEP_FS / 1000.0,
            "production_ps_per_replicate": PRODUCTION_STEPS * TIMESTEP_FS / 1000.0,
            "independent_seeds": list(REPLICATE_SEEDS),
            "barostat_interval_steps": 25,
            "constraint_and_force_parameters": "Loaded without change from the serialized OpenMM System XML.",
        },
        "measurements": measurements,
        "limitations": [
            "This system contains one PaG carboxylate in water and salt; it contains no insulin, NNC2215 macrocycle, or receptor.",
            "The run checks portability and short-term numerical stability only. It does not reproduce the published binding assay or 500 ns Desmond simulation.",
            "The implied one-ligand periodic-box concentration is not the 70 micromolar assay condition.",
            "A short solvent run does not validate binding, affinity, force-field transfer to NNC2215, or clinical effects.",
            "This free temporary GPU runtime is not an NVIDIA Inception or partner allocation.",
        ],
    }
    manifest_path = run_dir / "run.json"
    manifest_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    print(f"\nMANIFEST_PATH={manifest_path}")


if __name__ == "__main__":
    main()
