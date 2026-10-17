#!/usr/bin/env python3
"""Plot thermodynamic traces from the accepted PaG-only water pilot."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from io import StringIO
from pathlib import Path
from statistics import mean, stdev
from datetime import datetime, timezone

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parent
RUN_ID = "DHO-PAG-SOLVATION-20261017T191929Z"
RUN_DIR = ROOT / "data/pag-solvation-preflight"
OUTPUT = ROOT / "assets/pag-solvation-preflight.svg"
MANIFEST = RUN_DIR / f"{RUN_ID}.json"
ANALYSIS = RUN_DIR / f"{RUN_ID}-analysis.json"
REPLICATES = (
    ("Replicate 1", "#6757c5"),
    ("Replicate 2", "#66894f"),
)


def read_samples(path: Path) -> list[dict[str, float]]:
    text = path.read_text()
    if text.startswith("#"):
        text = text[1:]
    rows = list(csv.DictReader(StringIO(text)))
    return [
        {
            "time_ps": float(row["Time (ps)"]),
            "temperature_K": float(row["Temperature (K)"]),
            "volume_nm3": float(row["Box Volume (nm^3)"]),
        }
        for row in rows
    ]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def describe(values: list[float]) -> dict[str, float]:
    return {
        "mean": mean(values),
        "min": min(values),
        "max": max(values),
        "sample_sd": stdev(values) if len(values) > 1 else 0.0,
    }


def main() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "svg.fonttype": "none",
            "axes.edgecolor": "#d8d4ce",
            "axes.labelcolor": "#55515a",
            "xtick.color": "#77727c",
            "ytick.color": "#77727c",
            "text.color": "#25232a",
        }
    )
    fig, (temperature_ax, volume_ax) = plt.subplots(
        2,
        1,
        figsize=(10, 4.25),
        sharex=True,
        gridspec_kw={"height_ratios": [1.15, 1]},
        layout="constrained",
    )
    fig.patch.set_facecolor("#fbfaf8")

    manifest = json.loads(MANIFEST.read_text())
    equilibration_ps = manifest["protocol"]["equilibration_ps"]
    replicate_summaries = []
    for index, (label, color) in enumerate(REPLICATES, start=1):
        samples = read_samples(RUN_DIR / f"{RUN_ID}-replicate-{index:02d}.csv")
        times = [sample["time_ps"] for sample in samples]
        temperatures = [sample["temperature_K"] for sample in samples]
        volumes = [sample["volume_nm3"] for sample in samples]
        production = [sample for sample in samples if sample["time_ps"] > equilibration_ps]
        if not production:
            raise RuntimeError(f"No production samples after {equilibration_ps} ps.")
        if not all(
            math.isfinite(value)
            for sample in samples
            for value in sample.values()
        ):
            raise RuntimeError(f"Non-finite report sample in replicate {index}.")
        temperature_ax.plot(times, temperatures, color=color, lw=1.45, label=label)
        volume_ax.plot(times, volumes, color=color, lw=1.35, label=label)
        replicate_summaries.append(
            {
                "replicate": index,
                "log_file": f"{RUN_ID}-replicate-{index:02d}.csv",
                "log_sha256": sha256(RUN_DIR / f"{RUN_ID}-replicate-{index:02d}.csv"),
                "all_samples": len(samples),
                "equilibration_samples": sum(
                    sample["time_ps"] <= equilibration_ps for sample in samples
                ),
                "production_samples": len(production),
                "production_temperature_K": describe(
                    [sample["temperature_K"] for sample in production]
                ),
                "production_box_volume_nm3": describe(
                    [sample["volume_nm3"] for sample in production]
                ),
                "temperature_target_K": manifest["protocol"]["temperature_K"],
                "all_samples_finite": True,
            }
        )

    for axis in (temperature_ax, volume_ax):
        axis.axvspan(0, 10, color="#ece9e2", alpha=0.72, zorder=0)
        axis.set_xlim(0, 112)
        axis.grid(axis="y", color="#e8e5e0", linewidth=0.7)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
        axis.spines["left"].set_linewidth(0.7)
        axis.spines["bottom"].set_linewidth(0.7)
        axis.tick_params(labelsize=8, length=0, pad=6)

    temperature_ax.axhline(
        298.15, color="#827c87", lw=0.9, linestyle=(0, (3, 3)), zorder=1
    )
    temperature_ax.set_ylim(265, 335)
    temperature_ax.set_ylabel("Temperature (K)", fontsize=9, labelpad=9)
    temperature_ax.text(
        0.012,
        0.87,
        "TEMPERATURE",
        transform=temperature_ax.transAxes,
        fontsize=7,
        fontweight="bold",
        color="#79727e",
    )
    temperature_ax.text(
        0.988,
        0.87,
        "Target 298.15 K",
        transform=temperature_ax.transAxes,
        fontsize=7,
        ha="right",
        color="#77727c",
    )
    temperature_ax.legend(
        loc="upper right",
        bbox_to_anchor=(1.0, 1.24),
        ncol=2,
        frameon=False,
        fontsize=7,
        handlelength=1.8,
        columnspacing=1.4,
    )

    volume_ax.set_ylim(11.9, 13.8)
    volume_ax.set_ylabel("Box volume (nm³)", fontsize=9, labelpad=9)
    volume_ax.set_xlabel("Simulation time (ps)", fontsize=9, labelpad=8)
    volume_ax.text(
        0.012,
        0.82,
        "BOX VOLUME",
        transform=volume_ax.transAxes,
        fontsize=7,
        fontweight="bold",
        color="#79727e",
    )
    volume_ax.text(
        5,
        13.63,
        "10 ps equilibration",
        fontsize=6.5,
        color="#77727c",
        ha="left",
    )
    volume_ax.set_xticks([0, 20, 40, 60, 80, 100])
    fig.suptitle(
        "PaG-only explicit-water setup · two short replicas",
        x=0.055,
        y=1.035,
        ha="left",
        fontsize=12,
        fontweight="bold",
    )
    fig.text(
        0.055,
        -0.03,
        "298 K · 0.15 M NaCl · Apple M4 OpenCL  |  No NNC2215 or binding measurement",
        fontsize=7,
        color="#77727c",
    )
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT, format="svg", bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    analysis = {
        "record_type": "pag_solvation_preflight_postprocessing",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "run_id": RUN_ID,
        "run_manifest": MANIFEST.name,
        "run_manifest_sha256": sha256(MANIFEST),
        "analysis_script_sha256": sha256(Path(__file__)),
        "chart": str(OUTPUT.relative_to(ROOT)),
        "chart_sha256": sha256(OUTPUT),
        "equilibration_ps": equilibration_ps,
        "replicates": replicate_summaries,
        "interpretation": "Descriptive checks for finite thermodynamic traces and sampled temperatures after the equilibration interval. This postprocessing does not test PaG binding, conformational equilibrium, or transferability to the NNC2215 conjugate.",
        "clinical_use": False,
    }
    ANALYSIS.write_text(json.dumps(analysis, indent=2) + "\n")
    pointer_path = RUN_DIR / "latest.json"
    pointer = json.loads(pointer_path.read_text())
    pointer.update(
        {
            "updated_at_utc": datetime.now(timezone.utc).isoformat(),
            "analysis": ANALYSIS.name,
            "analysis_sha256": sha256(ANALYSIS),
            "figure": str(OUTPUT.relative_to(ROOT)),
            "figure_sha256": sha256(OUTPUT),
        }
    )
    pointer_path.write_text(json.dumps(pointer, indent=2) + "\n")
    print(json.dumps(analysis, indent=2))
    print(OUTPUT)


if __name__ == "__main__":
    main()
