#!/usr/bin/env python3
"""Run a small, generic protein-engine benchmark on the local machine.

This is an engine and environment preflight, not an NNC2215/PaG simulation.
It runs the same solvated 1UBQ coordinates through GROMACS and OpenMM, then
records one short throughput measurement per available engine/platform.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform as py_platform
import re
import shutil
import subprocess
import sys
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path

import openmm
from openmm import LocalEnergyMinimizer, Platform, unit
from openmm.app import GromacsGroFile, GromacsTopFile, HBonds, PME


ROOT = Path(__file__).resolve().parent
PDB = ROOT / "benchmarks" / "1ubq.pdb"
RESULT = ROOT / "data" / "engine-preflight.json"
CHART = ROOT / "assets" / "engine-preflight.svg"
RCSB_URL = "https://www.rcsb.org/structure/1UBQ"
FORCE_FIELD = "amber99sb-ildn"
WATER_MODEL = "tip3p"
BOX_PADDING_NM = 1.0
CUTOFF_NM = 1.0
TIMESTEP_PS = 0.002
TEMPERATURE_K = 300.0
FRICTION_PER_PS = 1.0
SEED = 1729
THREADS = 4


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(command: list[str], cwd: Path, stdin: str | None = None) -> str:
    completed = subprocess.run(
        command,
        cwd=cwd,
        input=stdin,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    if completed.returncode:
        raise RuntimeError(
            f"Command failed ({completed.returncode}): {' '.join(command)}\n"
            f"{completed.stdout[-8000:]}"
        )
    return completed.stdout


def gromacs_preflight(work: Path, gmx: str, production_steps: int) -> dict:
    gmx_root = Path(gmx).resolve().parents[1]
    top_dir = gmx_root / "share" / "gromacs" / "top"
    if not top_dir.is_dir():
        raise RuntimeError(f"Could not find GROMACS topology directory: {top_dir}")

    run(
        [gmx, "pdb2gmx", "-f", str(PDB), "-o", "processed.gro", "-p", "topol.top",
         "-ff", FORCE_FIELD, "-water", WATER_MODEL, "-ignh"],
        work,
        stdin="1\n1\n",
    )
    run([gmx, "editconf", "-f", "processed.gro", "-o", "boxed.gro", "-c", "-d",
         str(BOX_PADDING_NM), "-bt", "cubic"], work)
    run([gmx, "solvate", "-cp", "boxed.gro", "-cs", "spc216.gro", "-o",
         "solvated.gro", "-p", "topol.top"], work)

    minimization_mdp = """integrator = steep
emtol = 1000.0
emstep = 0.01
nsteps = 5000
cutoff-scheme = Verlet
nstlist = 20
rlist = 1.0
coulombtype = PME
rcoulomb = 1.0
fourierspacing = 0.12
pme-order = 4
vdwtype = cutoff
rvdw = 1.0
pbc = xyz
constraints = h-bonds
constraint-algorithm = lincs
lincs-iter = 1
lincs-order = 4
"""
    (work / "em.mdp").write_text(minimization_mdp, encoding="utf-8")
    run([gmx, "grompp", "-f", "em.mdp", "-c", "solvated.gro", "-p", "topol.top",
         "-o", "em.tpr", "-maxwarn", "1"], work)
    em_log = run([gmx, "mdrun", "-deffnm", "em", "-ntmpi", "1", "-ntomp",
                  str(THREADS)], work)
    em_steps_match = re.search(r"Steepest Descents converged to Fmax\s*=.*?in (\d+) steps", em_log, re.S)
    if not em_steps_match:
        em_steps_match = re.search(r"converged to Fmax.*?in (\d+) steps", (work / "em.log").read_text(), re.S)
    em_steps = int(em_steps_match.group(1)) if em_steps_match else None

    benchmark_mdp = f"""integrator = md
dt = {TIMESTEP_PS}
nsteps = {production_steps}
continuation = no
constraint-algorithm = lincs
constraints = h-bonds
lincs-iter = 1
lincs-order = 4
cutoff-scheme = Verlet
nstlist = 20
rlist = 1.0
coulombtype = PME
rcoulomb = 1.0
fourierspacing = 0.12
pme-order = 4
vdwtype = cutoff
rvdw = 1.0
pbc = xyz
tcoupl = v-rescale
tc-grps = System
tau-t = 1.0
ref-t = {TEMPERATURE_K}
pcoupl = no
gen-vel = yes
gen-temp = {TEMPERATURE_K}
gen-seed = {SEED}
comm-mode = Linear
nstxout = 0
nstvout = 0
nstfout = 0
nstenergy = 100
nstlog = 100
nstxout-compressed = 0
"""
    (work / "bench.mdp").write_text(benchmark_mdp, encoding="utf-8")
    run([gmx, "grompp", "-f", "bench.mdp", "-c", "em.gro", "-p", "topol.top",
         "-o", "bench.tpr", "-maxwarn", "1"], work)
    log = run([gmx, "mdrun", "-deffnm", "bench", "-ntmpi", "1", "-ntomp",
               str(THREADS)], work)
    performance = re.search(
        r"Performance:\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)", log
    )
    if not performance:
        raise RuntimeError("Could not read the GROMACS Performance line from mdrun output")
    gro = GromacsGroFile(str(work / "em.gro"))
    topology = GromacsTopFile(
        str(work / "topol.top"),
        periodicBoxVectors=gro.getPeriodicBoxVectors(),
        includeDir=str(top_dir),
    )
    version, build_report = gromacs_build_report(gmx)
    return {
        "engine": "GROMACS",
        "version": version,
        "platform": "CPU",
        "threads": THREADS,
        "gpu_support_in_build": "GPU support: enabled" in build_report,
        "steps": production_steps,
        "timestep_ps": TIMESTEP_PS,
        "ns_per_day": float(performance.group(1)),
        "hours_per_ns": float(performance.group(2)),
        "milliseconds_per_step": float(performance.group(3)),
        "million_atom_steps_per_second": float(performance.group(4)),
        "minimization_steps": em_steps,
        "solvated_atom_count": len(list(topology.topology.atoms())),
        "timing_note": "GROMACS mdrun-reported wall time, including engine run setup.",
    }


def gromacs_build_report(gmx: str) -> tuple[str, str]:
    output = run([gmx, "--version"], Path.cwd())
    match = re.search(r"GROMACS version:\s*(.+)", output)
    version = match.group(1).strip() if match else output.splitlines()[0].strip()
    return version, output


def openmm_preflight(work: Path, platform_name: str, steps: int, warmup_steps: int) -> dict:
    gmx = shutil.which("gmx")
    if not gmx:
        raise RuntimeError("GROMACS executable is not available on PATH")
    gmx_root = Path(gmx).resolve().parents[1]
    gro = GromacsGroFile(str(work / "em.gro"))
    top = GromacsTopFile(
        str(work / "topol.top"),
        periodicBoxVectors=gro.getPeriodicBoxVectors(),
        includeDir=str(gmx_root / "share" / "gromacs" / "top"),
    )
    system = top.createSystem(
        nonbondedMethod=PME,
        nonbondedCutoff=CUTOFF_NM * unit.nanometer,
        constraints=HBonds,
        rigidWater=True,
        ewaldErrorTolerance=0.0005,
    )
    integrator = openmm.LangevinMiddleIntegrator(
        TEMPERATURE_K * unit.kelvin,
        FRICTION_PER_PS / unit.picosecond,
        TIMESTEP_PS * unit.picoseconds,
    )
    integrator.setRandomNumberSeed(SEED)

    platform = Platform.getPlatformByName(platform_name)
    available = set(platform.getPropertyNames())
    properties: dict[str, str] = {}
    if platform_name == "CPU":
        if "Threads" in available:
            properties["Threads"] = str(THREADS)
        if "DeterministicForces" in available:
            properties["DeterministicForces"] = "true"
    if platform_name == "OpenCL":
        if "Precision" in available:
            properties["Precision"] = "single"
        if "DeviceIndex" in available:
            properties["DeviceIndex"] = "0"

    simulation = openmm.app.Simulation(top.topology, system, integrator, platform, properties)
    simulation.context.setPositions(gro.positions)
    simulation.context.setVelocitiesToTemperature(TEMPERATURE_K * unit.kelvin, SEED)
    LocalEnergyMinimizer.minimize(simulation.context, 10.0 * unit.kilojoule_per_mole / unit.nanometer, 1000)
    simulation.step(warmup_steps)
    initial_energy = simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
        unit.kilojoule_per_mole
    )
    if not math.isfinite(initial_energy):
        raise RuntimeError(f"{platform_name} produced a non-finite potential energy")

    started = time.perf_counter()
    simulation.step(steps)
    elapsed = time.perf_counter() - started
    simulated_ns = steps * TIMESTEP_PS / 1000.0
    throughput = simulated_ns / elapsed * 86400.0
    atoms = system.getNumParticles()
    simulation.context.getState(getEnergy=True)
    del simulation
    return {
        "engine": "OpenMM",
        "version": openmm.__version__,
        "platform": platform_name,
        "platform_properties": properties,
        "threads": THREADS if platform_name == "CPU" else None,
        "steps": steps,
        "warmup_steps_excluded_from_timing": warmup_steps,
        "timestep_ps": TIMESTEP_PS,
        "timed_seconds": round(elapsed, 6),
        "ns_per_day": round(throughput, 3),
        "million_atom_steps_per_second": round(atoms * steps / elapsed / 1_000_000, 3),
        "solvated_atom_count": atoms,
        "initial_potential_energy_kj_mol_after_warmup": round(initial_energy, 3),
        "timing_note": "OpenMM wall time for production steps only, after 500 warmup steps; setup and minimization excluded.",
    }


def render_chart(records: list[dict], atom_count: int, run_date: str) -> None:
    width, height = 1000, 390
    left, right = 290, 920
    chart_top = 98
    row_gap = 70
    max_value = max(item["ns_per_day"] for item in records) * 1.12
    labels = []
    for index, item in enumerate(records):
        y = chart_top + index * row_gap
        if item["engine"] == "GROMACS":
            version = item["version"].split("-", 1)[0]
            label = f"GROMACS {version} · CPU (4 threads)"
            color = "#8573b8"
        elif item["platform"] == "CPU":
            label = f"OpenMM {item['version']} · CPU, 4 threads"
            color = "#b4a6d3"
        else:
            label = f"OpenMM {item['version']} · {item['platform']} / Apple M4"
            color = "#9ebd77"
        bar_width = (item["ns_per_day"] / max_value) * (right - left)
        labels.extend([
            f'<text x="32" y="{y + 5}" font-family="Arial, sans-serif" font-size="16" fill="#40394b">{label}</text>',
            f'<rect x="{left}" y="{y - 13}" width="{right - left}" height="18" rx="9" fill="#eeeaf2"/>',
            f'<rect x="{left}" y="{y - 13}" width="{bar_width:.1f}" height="18" rx="9" fill="{color}"/>',
            f'<text x="{min(left + bar_width + 12, 940):.1f}" y="{y + 5}" font-family="Arial, sans-serif" font-size="15" font-weight="700" fill="#40394b">{item["ns_per_day"]:.2f} ns/day</text>',
        ])
    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" role="img" aria-labelledby="title desc">',
        '<title id="title">Local molecular-dynamics engine preflight on ubiquitin</title>',
        '<desc id="desc">Single short throughput runs for GROMACS CPU, OpenMM CPU, and OpenMM OpenCL on an Apple M4 using the same solvated ubiquitin test system. This is not an NNC2215 or NVIDIA run.</desc>',
        '<rect width="100%" height="100%" rx="14" fill="#ffffff"/>',
        f'<text x="32" y="35" font-family="Arial, sans-serif" font-size="11" letter-spacing="2" fill="#75678f">LOCAL ENGINE PREFLIGHT · {run_date}</text>',
        '<text x="32" y="64" font-family="Arial, sans-serif" font-size="20" font-weight="600" fill="#24212b">Short-run throughput on one standard protein system</text>',
        *labels,
        '<line x1="32" y1="319" x2="968" y2="319" stroke="#ece9e4"/>',
        f'<text x="32" y="345" font-family="Arial, sans-serif" font-size="11" fill="#77737e">1UBQ ubiquitin · {atom_count:,} solvated atoms · 5,000 timed steps · single run per engine</text>',
        '<text x="32" y="367" font-family="Arial, sans-serif" font-size="10" fill="#918c96">Different integrators and timing boundaries; indicative only. No target MD or NVIDIA GPU allocation.</text>',
        '</svg>',
    ]
    CHART.parent.mkdir(parents=True, exist_ok=True)
    CHART.write_text("\n".join(svg), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--steps", type=int, default=5000, help="Timed production steps per engine")
    parser.add_argument("--warmup-steps", type=int, default=500, help="OpenMM warmup steps excluded from timing")
    args = parser.parse_args()
    if args.steps < 1 or args.warmup_steps < 0:
        parser.error("--steps must be positive and --warmup-steps cannot be negative")
    if not PDB.is_file():
        raise SystemExit(f"Missing benchmark input: {PDB}")

    gmx = shutil.which("gmx")
    if not gmx:
        raise SystemExit("GROMACS executable 'gmx' is not available on PATH")
    started_at = datetime.now(timezone.utc)
    with tempfile.TemporaryDirectory(prefix="dho-engine-preflight-") as temp_dir:
        work = Path(temp_dir)
        gromacs = gromacs_preflight(work, gmx, args.steps)
        openmm_records = []
        platform_status = []
        for name in ("CPU", "OpenCL"):
            try:
                Platform.getPlatformByName(name)
            except Exception as error:
                platform_status.append({"platform": name, "available": False, "reason": str(error)})
                continue
            openmm_records.append(openmm_preflight(work, name, args.steps, args.warmup_steps))
            platform_status.append({"platform": name, "available": True})

        records = [gromacs, *openmm_records]
        atom_counts = {record["solvated_atom_count"] for record in records}
        if len(atom_counts) != 1:
            raise RuntimeError(f"Engine topology atom counts differ: {sorted(atom_counts)}")
        atom_count = atom_counts.pop()
        completed_at = datetime.now(timezone.utc)
        cpu_model = ""
        if sys.platform == "darwin":
            cpu_model = subprocess.run(
                ["sysctl", "-n", "machdep.cpu.brand_string"],
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.DEVNULL,
                check=False,
            ).stdout.strip()
        result = {
            "schema": "dia-health-one.engine-preflight.v1",
            "run": {
                "id": "DHO-ENGINE-PREFLIGHT-" + completed_at.strftime("%Y%m%dT%H%M%SZ"),
                "started_at_utc": started_at.isoformat(),
                "completed_at_utc": completed_at.isoformat(),
                "script_sha256": sha256(Path(__file__).resolve()),
                "input_sha256": sha256(PDB),
                "software_environment": "research/engine-benchmark-environment.yml",
                "python_version": sys.version.split()[0],
                "host_os": py_platform.platform(),
                "host_architecture": py_platform.machine(),
                "host_processor": cpu_model or py_platform.processor() or "not reported",
                "nvidia_gpu_available": False,
                "nvidia_gpu_used": False,
                "clinical_use": False,
            },
            "source_system": {
                "structure_id": "1UBQ",
                "source_url": RCSB_URL,
                "structure_class": "public experimental X-ray structure of ubiquitin",
                "protein_residues": 76,
                "force_field": "Amber99SB-ILDN",
                "water_model": "TIP3P",
                "water_box_coordinates": "GROMACS spc216.gro three-site water coordinates",
                "box_padding_nm": BOX_PADDING_NM,
                "added_salt": False,
                "solvated_atom_count": atom_count,
            },
            "protocol": {
                "purpose": "engine throughput and installation preflight only",
                "timed_steps": args.steps,
                "openmm_warmup_steps_excluded_from_timing": args.warmup_steps,
                "timestep_ps": TIMESTEP_PS,
                "temperature_K": TEMPERATURE_K,
                "random_seed": SEED,
                "threads": THREADS,
                "gromacs_integrator": "md with v-rescale thermostat",
                "openmm_integrator": "LangevinMiddleIntegrator",
                "limitations": [
                    "One short run per engine/platform; timing is not a statistically sampled benchmark.",
                    "GROMACS and OpenMM use different thermostat implementations.",
                    "GROMACS mdrun timing includes run setup; OpenMM timing excludes setup, minimization, and warmup.",
                    "The test protein is not NNC2215, PaG, or an insulin-receptor complex.",
                    "The OpenCL result is from Apple M4 hardware, not an NVIDIA GPU.",
                    "Throughput does not validate force-field quality or predict binding affinity.",
                ],
            },
            "available_openmm_platforms": platform_status,
            "measurements": records,
            "interpretation": {
                "status": "local generic engine preflight complete",
                "target_nnc2215_pag_md_run": False,
                "nvidia_acceleration_demonstrated": False,
                "scale_up_recommendation": "Repeat a controlled benchmark on the allocated NVIDIA system; do not use this preflight to select production MD parameters or infer NNC2215/PaG energetics.",
            },
        }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    render_chart(records, atom_count, completed_at.strftime("%Y-%m-%d UTC"))
    print(json.dumps({
        "run_id": result["run"]["id"],
        "system_atoms": atom_count,
        "measurements": [
            {"engine": item["engine"], "platform": item["platform"], "ns_per_day": item["ns_per_day"]}
            for item in records
        ],
        "output": str(RESULT.relative_to(ROOT)),
        "chart": str(CHART.relative_to(ROOT)),
    }, indent=2))


if __name__ == "__main__":
    main()
