#!/usr/bin/env python3
"""Reconstruct the published GluHUT/glucose PMF from public CV time series.

This reproduces an analogue-system analysis from the Zenodo source archive.
It does not simulate NNC2215 or PaG and must not be presented as a target run.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import shutil
import statistics
import subprocess
import tempfile
import urllib.request
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).parent
ARCHIVE_URL = "https://zenodo.org/api/records/19664426/files/WHAM_PMF.zip/content"
WHAM_REPOSITORY = "https://github.com/agrossfield/wham"
WHAM_SOURCE_COMMIT = "43bdd6fcd1916ae2b169cb52ac326b7eea567abb"
WHAM_VERSION = "2.0.11"
ARTICLE_URL = "https://pubs.rsc.org/en/content/articlehtml/2026/sc/d6sc00738d"
DATASET_URL = "https://zenodo.org/records/19664426"
RUNS = (1, 2, 3)
WINDOW_COUNT = 51
KCAL_MOL_A2_IN_WHAM_METADATA = 100.0
PUBLISHED_SPRING_CONSTANT_KCAL_MOL_A2 = 50.0


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def obtain_archive(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        request = urllib.request.Request(
            ARCHIVE_URL,
            headers={"User-Agent": "DiaHealthOne-Reproducibility/1.0"},
        )
        with urllib.request.urlopen(request, timeout=120) as response, path.open("wb") as out:
            shutil.copyfileobj(response, out)
    with zipfile.ZipFile(path) as archive:
        bad_member = archive.testzip()
        if bad_member:
            raise ValueError(f"Zenodo ZIP member failed CRC check: {bad_member}")
    return path


def parse_window(path: str) -> tuple[float, str]:
    match = re.search(r"/(\d+)-(\d)-test-dist_vs_t\.dat$", path)
    if not match:
        raise ValueError(f"Unexpected GluHUT window filename: {path}")
    center = float(f"{match.group(1)}.{match.group(2)}")
    return center, path.rsplit("/", 1)[-1]


def pmf_from_file(path: Path) -> list[tuple[float, float]]:
    values = []
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split()
        coordinate, free_energy = float(fields[0]), float(fields[1])
        if math.isfinite(coordinate) and math.isfinite(free_energy):
            values.append((coordinate, free_energy))
    if len(values) != 400:
        raise ValueError(f"Expected 400 PMF bins in {path}; found {len(values)}.")
    minimum = min(energy for _, energy in values)
    return [(coordinate, energy - minimum) for coordinate, energy in values]


def build_plot(
    coordinates: list[float],
    profiles: list[list[float]],
    mean_profile: list[float],
    sd_profile: list[float],
    binding_energies: list[float],
    output: Path,
) -> None:
    colors = ["#b9c5d5", "#8aa1b9", "#6d8199"]
    navy = "#18344b"
    teal = "#16877f"
    fig, ax = plt.subplots(figsize=(12.8, 5.8), dpi=160)
    fig.patch.set_facecolor("#f6f8f8")
    ax.set_facecolor("#f6f8f8")

    for profile, color, number in zip(profiles, colors, RUNS):
        ax.plot(coordinates, profile, color=color, linewidth=1.1, alpha=0.82,
                label=f"Independent PMF {number}")
    ax.fill_between(
        coordinates,
        [value - error for value, error in zip(mean_profile, sd_profile)],
        [value + error for value, error in zip(mean_profile, sd_profile)],
        color=teal,
        alpha=0.13,
        linewidth=0,
        label="Between-run ±1 SD",
    )
    ax.plot(coordinates, mean_profile, color=teal, linewidth=3.0,
            label="Three-run mean")
    ax.axvspan(13, 15, color="#d8e9e5", alpha=0.62, zorder=0)
    ax.text(14, ax.get_ylim()[1] * 0.93, "BULK REFERENCE\n13–15 Å", ha="center",
            va="top", fontsize=8, color="#47796f", fontweight="bold")

    ax.set_xlim(0, 15)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Host–guest centre-of-mass distance (Å)", fontsize=10, color=navy, labelpad=9)
    ax.set_ylabel("Relative free energy (kcal mol⁻¹)", fontsize=10, color=navy, labelpad=9)
    ax.set_title("GluHUT · glucose unbinding profile", loc="left", fontsize=19,
                 fontweight="bold", color=navy, pad=20)
    ax.text(0, 1.015,
            "Published umbrella windows re-analysed with Grossfield WHAM 2.0.11 · 3 independent PMFs",
            transform=ax.transAxes, fontsize=9, color="#637587", va="bottom")
    ax.grid(axis="y", color="#dce3e6", linewidth=0.8)
    ax.grid(axis="x", visible=False)
    ax.tick_params(axis="both", colors="#637587", labelsize=8, length=0, pad=7)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_color("#bac8cd")
    ax.legend(loc="upper left", frameon=False, ncol=2, fontsize=8,
              bbox_to_anchor=(0, 1.0), borderaxespad=0.15)

    avg = statistics.mean(binding_energies)
    sd = statistics.stdev(binding_energies)
    ax.text(0.98, 0.34,
            f"ΔG = {avg:.2f} ± {sd:.2f} kcal mol⁻¹\n"
            "source reports −12.5 ± 0.9",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=9,
            color=navy,
            bbox={"boxstyle": "round,pad=0.55", "facecolor": "white",
                  "edgecolor": "#dce3e6", "linewidth": 0.8})
    fig.text(0.01, 0.015,
             "Public GluHUT analogue data · not NNC2215/PaG · mean ± SD across independent PMFs",
             color="#637587", fontsize=8)
    fig.subplots_adjust(left=0.09, right=0.98, top=0.81, bottom=0.15)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, format="svg", facecolor=fig.get_facecolor())
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--archive",
        type=Path,
        default=Path("/tmp/dho-gluhut-wham-pmf.zip"),
        help="Zenodo WHAM_PMF.zip; downloaded automatically if missing.",
    )
    parser.add_argument(
        "--wham-binary",
        type=Path,
        required=True,
        help="Grossfield WHAM 2.0.11 binary built from the recorded source commit.",
    )
    parser.add_argument(
        "--output-json", type=Path,
        default=ROOT / "data/gluhut-glucose-wham-reproduction.json",
    )
    parser.add_argument(
        "--output-svg", type=Path,
        default=ROOT / "assets/gluhut-glucose-pmf.svg",
    )
    args = parser.parse_args()
    archive_path = obtain_archive(args.archive)
    if not args.wham_binary.is_file():
        raise FileNotFoundError(args.wham_binary)

    source_prefix = "WHAM_PMF/"
    input_profiles: list[list[tuple[float, float]]] = []
    binding_energies: list[float] = []
    frame_counts: list[int] = []
    convergence_iterations: list[int | None] = []
    with zipfile.ZipFile(archive_path) as archive, tempfile.TemporaryDirectory(
        prefix="dho-gluhut-wham-"
    ) as temp_name:
        temp_root = Path(temp_name)
        example_meta_name = (
            "WHAM_PMF/example-meta-file-histogram-collection-glucose-US-1.dat"
        )
        example_meta = archive.read(example_meta_name).decode()
        for run in RUNS:
            source_dir = f"{source_prefix}histogram-collection-glucose-US-{run}/"
            names = sorted(
                [
                    name for name in archive.namelist()
                    if name.startswith(source_dir) and name.endswith("-test-dist_vs_t.dat")
                ],
                key=lambda name: parse_window(name)[0],
            )
            if len(names) != WINDOW_COUNT:
                raise ValueError(f"Expected 51 glucose umbrella windows for run {run}; found {len(names)}.")
            window_metadata = []
            run_dir = temp_root / source_dir
            run_dir.mkdir(parents=True, exist_ok=True)
            for name in names:
                center, basename = parse_window(name)
                if not math.isclose(center / 0.3, round(center / 0.3), abs_tol=1e-7):
                    raise ValueError(f"Unexpected umbrella center spacing: {center} Å")
                data = archive.read(name)
                local_file = run_dir / basename
                local_file.write_bytes(data)
                frame_counts.append(sum(1 for line in data.splitlines() if line.strip()))
                window_metadata.append(
                    f"../histogram-collection-glucose-US-{run}/{basename} "
                    f"{center:.1f} {KCAL_MOL_A2_IN_WHAM_METADATA:.1f}"
                )

            if run == 1:
                # Use the published example metadata verbatim where possible.
                parsed_example = [line.split() for line in example_meta.splitlines() if line.strip()]
                parsed_generated = [line.split() for line in window_metadata]
                if len(parsed_example) != WINDOW_COUNT:
                    raise ValueError("The published example metadata does not contain 51 windows.")
                for example, generated in zip(parsed_example, parsed_generated):
                    if not (
                        math.isclose(float(example[1]), float(generated[1]), abs_tol=1e-9)
                        and math.isclose(float(example[2]), float(generated[2]), abs_tol=1e-9)
                    ):
                        raise ValueError("Generated run-1 centers/force constants differ from source metadata.")
                    if Path(example[0]).name != Path(generated[0]).name:
                        raise ValueError("Generated run-1 filenames differ from source metadata.")
                metadata_text = example_meta
                metadata_source = "Zenodo example metadata, validated against extracted run-1 windows"
            else:
                metadata_text = "\n".join(window_metadata) + "\n"
                metadata_source = "generated using the Zenodo run-1 metadata convention"
            metadata_path = run_dir / f"metadata-glucose-US-{run}.dat"
            metadata_path.write_text(metadata_text)
            output_name = f"glucose-pmf-US-{run}.dat"
            process = subprocess.run(
                [
                    str(args.wham_binary), "0.0", "15.0", "400", "0.000001",
                    "300", "0", metadata_path.name, output_name,
                ],
                cwd=run_dir,
                check=False,
                capture_output=True,
                text=True,
                timeout=3600,
            )
            if process.returncode != 0:
                raise RuntimeError(
                    f"Grossfield WHAM failed for glucose run {run}: "
                    f"{process.stderr[-2000:]}"
                )
            iteration_numbers = [
                int(value)
                for value in re.findall(r"#Iteration\s+(\d+):", process.stdout)
            ]
            convergence_iterations.append(max(iteration_numbers) if iteration_numbers else None)
            profile = pmf_from_file(run_dir / output_name)
            input_profiles.append(profile)
            bulk = [
                energy for coordinate, energy in profile
                if 13.0 <= coordinate <= 15.0
            ]
            minimum = min(energy for _, energy in profile)
            binding_energies.append(minimum - statistics.mean(bulk))

    coordinates = [coordinate for coordinate, _ in input_profiles[0]]
    if any([coordinate for coordinate, _ in profile] != coordinates for profile in input_profiles[1:]):
        raise ValueError("WHAM output grids differ across independent runs.")
    profiles = [[energy for _, energy in profile] for profile in input_profiles]
    mean_profile = [statistics.mean(values) for values in zip(*profiles)]
    sd_profile = [statistics.stdev(values) for values in zip(*profiles)]
    reproduced_mean = statistics.mean(binding_energies)
    reproduced_sd = statistics.stdev(binding_energies)
    published_mean = -12.5
    published_sd = 0.9

    build_plot(
        coordinates, profiles, mean_profile, sd_profile,
        binding_energies, args.output_svg,
    )
    report = {
        "record_type": "public_analogue_pmf_reproduction",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "system": "GluHUT host with D-glucose; public methodological analogue only",
        "not_target_claim": "This is not an NNC2215–PaG simulation or a target reproduction.",
        "sources": {
            "article": {
                "citation": "Eades et al., Chem. Sci. 2026, 17, 10896–10906",
                "doi": "10.1039/D6SC00738D",
                "url": ARTICLE_URL,
                "license": "CC BY 3.0",
                "published_glucose_binding_free_energy_kcal_mol": {
                    "mean": published_mean,
                    "sd": published_sd,
                },
            },
            "dataset": {
                "doi": "10.5281/zenodo.19664426",
                "url": DATASET_URL,
                "license": "CC BY 4.0",
                "file": "WHAM_PMF.zip",
                "download_url": ARCHIVE_URL,
                "archive_sha256": sha256_file(archive_path),
            },
        },
        "analysis": {
            "method": "Weighted Histogram Analysis Method on deposited reaction-coordinate time series",
            "software": {
                "name": "Grossfield WHAM",
                "version": WHAM_VERSION,
                "repository": WHAM_REPOSITORY,
                "source_commit": WHAM_SOURCE_COMMIT,
                "binary_sha256": sha256_file(args.wham_binary),
            },
            "parameters": {
                "reaction_coordinate": "distance between sugar-ring COM and GluHUT cavity COM, Å",
                "range_angstrom": [0.0, 15.0],
                "bins": 400,
                "temperature_kelvin": 300,
                "tolerance": 1e-6,
                "periodic": False,
                "wham_metadata_spring_kcal_mol_angstrom2": KCAL_MOL_A2_IN_WHAM_METADATA,
                "article_reported_amber_spring_kcal_mol_angstrom2": PUBLISHED_SPRING_CONSTANT_KCAL_MOL_A2,
                "spring_convention_note": "Zenodo example metadata uses 100; the paper reports Amber restraint constant 50. Grossfield WHAM documentation explains that Amber omits the 1/2 factor, so WHAM metadata uses twice the Amber value.",
                "bulk_region_angstrom": [13.0, 15.0],
                "binding_free_energy_definition": "global PMF minimum minus mean PMF over the 13–15 Å bulk region",
            },
            "replicates": [
                {
                    "run": run,
                    "windows": WINDOW_COUNT,
                    "frames_across_windows": sum(frame_counts[(run - 1) * WINDOW_COUNT:run * WINDOW_COUNT]),
                    "wham_convergence_iteration": convergence_iterations[run - 1],
                    "binding_free_energy_kcal_mol": binding_energies[run - 1],
                    "metadata_source": (
                        "Zenodo example metadata, validated against extracted run-1 windows"
                        if run == 1 else "generated using the Zenodo run-1 metadata convention"
                    ),
                }
                for run in RUNS
            ],
            "result": {
                "mean_binding_free_energy_kcal_mol": reproduced_mean,
                "sd_across_independent_pmfs_kcal_mol": reproduced_sd,
                "published_mean_kcal_mol": published_mean,
                "published_sd_kcal_mol": published_sd,
                "difference_from_published_mean_kcal_mol": reproduced_mean - published_mean,
                "mean_within_published_one_sd": abs(reproduced_mean - published_mean) <= published_sd,
                "mean_profile": [
                    {"coordinate_angstrom": x, "mean_kcal_mol": mean, "sd_kcal_mol": sd}
                    for x, mean, sd in zip(coordinates, mean_profile, sd_profile)
                ],
                "independent_profiles_kcal_mol": [
                    {"run": run, "values": profile}
                    for run, profile in zip(RUNS, profiles)
                ],
            },
        },
        "compute": {
            "execution": "CPU WHAM post-processing; no molecular dynamics was run for this reproduction",
            "nvidia_gpu_used": False,
        },
        "limitations": [
            "GluHUT is a glucose-binding receptor analogue and is not the NNC2215 insulin conjugate.",
            "No PaG was simulated or analysed in this PMF reproduction.",
            "The deposited trajectories are original authors' simulation outputs; this workflow redoes their WHAM analysis.",
            "This source-data reproduction does not establish NNC2215–PaG binding, clinical effects, or force-field validity for the insulin conjugate.",
        ],
        "clinical_use": False,
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({
        "output_json": str(args.output_json),
        "output_svg": str(args.output_svg),
        "replicate_binding_free_energies_kcal_mol": binding_energies,
        "reproduced_mean_kcal_mol": reproduced_mean,
        "reproduced_sd_kcal_mol": reproduced_sd,
        "published_mean_kcal_mol": published_mean,
        "published_sd_kcal_mol": published_sd,
        "difference_kcal_mol": reproduced_mean - published_mean,
    }, indent=2))


if __name__ == "__main__":
    main()
