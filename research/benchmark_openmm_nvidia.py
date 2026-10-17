#!/usr/bin/env python3
"""Benchmark generic 1UBQ dynamics with OpenMM OpenCL on an NVIDIA GPU.

This is a hardware/engine preflight. It is not an NNC2215 or PaG simulation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform as host_platform
import statistics
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path

import openmm
from openmm import LocalEnergyMinimizer, Platform, app, unit


ROOT = Path(__file__).resolve().parent
INPUT = ROOT / "benchmarks" / "1ubq.pdb"
DEFAULT_OUTPUT = ROOT / "data" / "nvidia-engine-preflight.json"
SOURCE_URL = "https://www.rcsb.org/structure/1UBQ"
TIMESTEP_FS = 2.0
TEMPERATURE_K = 300.0
FRICTION_PER_PS = 1.0
WARMUP_STEPS = 1_000
TIMED_STEPS = 50_000
REPLICATE_SEEDS = (2215001, 2215002, 2215003)
PADDING_NM = 1.0
CUTOFF_NM = 1.0


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def gpu_inventory() -> dict[str, str]:
    result = subprocess.run(
        [
            "nvidia-smi",
            "--query-gpu=name,driver_version,memory.total",
            "--format=csv,noheader",
        ],
        capture_output=True,
        text=True,
        check=True,
    )
    first_gpu = result.stdout.strip().splitlines()[0]
    name, driver, memory = (part.strip() for part in first_gpu.split(",", 2))
    return {
        "model": name,
        "driver_version": driver,
        "memory_total": memory,
    }


def available_platforms() -> list[dict[str, object]]:
    return [
        {
            "name": Platform.getPlatform(index).getName(),
            "properties": list(Platform.getPlatform(index).getPropertyNames()),
        }
        for index in range(Platform.getNumPlatforms())
    ]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    input_path = args.input.resolve()
    output_path = args.output.resolve()
    if not input_path.is_file():
        raise FileNotFoundError(input_path)

    started_at = datetime.now(timezone.utc)
    run_id = "DHO-NVIDIA-OPENMM-" + started_at.strftime("%Y%m%dT%H%M%SZ")
    gpu = gpu_inventory()
    platform = Platform.getPlatformByName("OpenCL")
    available_properties = set(platform.getPropertyNames())
    requested_properties = {"DeviceIndex": "0", "Precision": "single"}
    properties = {
        name: value
        for name, value in requested_properties.items()
        if name in available_properties
    }
    if "DeviceName" not in available_properties:
        raise RuntimeError("OpenMM OpenCL does not expose its selected device name")

    pdb = app.PDBFile(str(input_path))
    forcefield = app.ForceField("amber14-all.xml", "amber14/tip3p.xml")
    modeller = app.Modeller(pdb.topology, pdb.positions)
    modeller.addHydrogens(forcefield)
    modeller.addSolvent(
        forcefield,
        model="tip3p",
        padding=PADDING_NM * unit.nanometer,
    )
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=CUTOFF_NM * unit.nanometer,
        constraints=app.HBonds,
        rigidWater=True,
        ewaldErrorTolerance=0.0005,
    )
    integrator = openmm.LangevinMiddleIntegrator(
        TEMPERATURE_K * unit.kelvin,
        FRICTION_PER_PS / unit.picosecond,
        TIMESTEP_FS * unit.femtoseconds,
    )
    simulation = app.Simulation(
        modeller.topology,
        system,
        integrator,
        platform,
        properties,
    )
    simulation.context.setPositions(modeller.positions)
    LocalEnergyMinimizer.minimize(
        simulation.context,
        10.0 * unit.kilojoule_per_mole / unit.nanometer,
        1_000,
    )
    minimized_positions = simulation.context.getState(getPositions=True).getPositions()
    device_name = platform.getPropertyValue(simulation.context, "DeviceName")
    if gpu["model"].lower() not in device_name.lower():
        raise RuntimeError(
            f"OpenMM selected a different device, found GPU={gpu['model']!r}, "
            f"OpenMM device={device_name!r}"
        )

    measurements = []
    for seed in REPLICATE_SEEDS:
        integrator.setRandomNumberSeed(seed)
        simulation.context.setPositions(minimized_positions)
        simulation.context.setVelocitiesToTemperature(
            TEMPERATURE_K * unit.kelvin,
            seed,
        )
        simulation.step(WARMUP_STEPS)
        warmup_state = simulation.context.getState(getEnergy=True)
        initial_energy = warmup_state.getPotentialEnergy().value_in_unit(
            unit.kilojoule_per_mole
        )
        if not math.isfinite(initial_energy):
            raise RuntimeError(f"Non-finite energy after warmup for seed {seed}")

        started = time.perf_counter()
        simulation.step(TIMED_STEPS)
        # Synchronize device work before stopping the timer.
        final_state = simulation.context.getState(getEnergy=True)
        elapsed = time.perf_counter() - started
        final_energy = final_state.getPotentialEnergy().value_in_unit(
            unit.kilojoule_per_mole
        )
        if not math.isfinite(final_energy) or elapsed <= 0:
            raise RuntimeError(
                f"Invalid timed result for seed {seed}: elapsed={elapsed}, energy={final_energy}"
            )
        simulated_ns = TIMED_STEPS * TIMESTEP_FS / 1_000_000.0
        measurements.append(
            {
                "seed": seed,
                "timed_steps": TIMED_STEPS,
                "warmup_steps": WARMUP_STEPS,
                "simulated_ns": simulated_ns,
                "elapsed_seconds": round(elapsed, 4),
                "ns_per_day": round(simulated_ns / elapsed * 86_400.0, 3),
                "potential_energy_after_warmup_kj_mol": round(initial_energy, 3),
                "potential_energy_after_timing_kj_mol": round(final_energy, 3),
            }
        )

    rates = [item["ns_per_day"] for item in measurements]
    completed_at = datetime.now(timezone.utc)
    result = {
        "schema": "dia-health-one.openmm-gpu-preflight.v1",
        "run": {
            "id": run_id,
            "started_at_utc": started_at.isoformat(),
            "completed_at_utc": completed_at.isoformat(),
            "status": "completed",
            "script": "research/benchmark_openmm_nvidia.py",
            "script_sha256": sha256(Path(__file__).resolve()),
            "input": str(input_path),
            "input_sha256": sha256(input_path),
            "environment": {
                "openmm": openmm.__version__,
                "python": host_platform.python_version(),
                "os": host_platform.platform(),
                "architecture": host_platform.machine(),
            },
            "gpu": gpu,
            "purpose": "Generic 1UBQ OpenMM OpenCL throughput preflight; not NNC2215/PaG target MD.",
            "clinical_use": False,
        },
        "system": {
            "source_id": "1UBQ",
            "source_url": SOURCE_URL,
            "structure_class": "public experimental X-ray ubiquitin structure",
            "force_field": "Amber14-all",
            "water_model": "TIP3P",
            "solvent_padding_nm": PADDING_NM,
            "salt_added": False,
            "atom_count": system.getNumParticles(),
            "nonbonded_method": "PME",
            "nonbonded_cutoff_nm": CUTOFF_NM,
            "constraints": "HBonds",
        },
        "protocol": {
            "engine": "OpenMM",
            "version": openmm.__version__,
            "platform": "OpenCL",
            "platform_properties": properties,
            "device_name": device_name,
            "available_platforms": available_platforms(),
            "integrator": "LangevinMiddleIntegrator",
            "ensemble": "NVT",
            "temperature_K": TEMPERATURE_K,
            "friction_per_ps": FRICTION_PER_PS,
            "timestep_fs": TIMESTEP_FS,
            "minimizer_tolerance_kj_mol_nm": 10.0,
            "warmup_steps": WARMUP_STEPS,
            "timed_steps_per_replicate": TIMED_STEPS,
            "independent_velocity_seeds": list(REPLICATE_SEEDS),
            "synchronization_in_timed_interval": True,
        },
        "measurements": measurements,
        "summary": {
            "replicates": len(measurements),
            "median_ns_per_day": round(statistics.median(rates), 3),
            "mean_ns_per_day": round(statistics.mean(rates), 3),
            "sample_sd_ns_per_day": round(statistics.stdev(rates), 3),
            "min_ns_per_day": round(min(rates), 3),
            "max_ns_per_day": round(max(rates), 3),
        },
        "limitations": [
            "This generic ubiquitin system is not NNC2215, PaG, or an insulin-receptor complex.",
            "OpenMM used the OpenCL platform on the Tesla T4; a CUDA platform was not present in this runtime.",
            "The Amber14-all/TIP3P setup differs from the earlier Apple M4 benchmark and is not a controlled cross-machine comparison.",
            "Measured throughput is not target-system performance and does not validate conjugate parameters or predict binding affinity.",
            "The free Colab runtime is temporary; it is not an NVIDIA Inception or partner GPU allocation.",
        ],
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
