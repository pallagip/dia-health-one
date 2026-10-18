#!/usr/bin/env python3
"""Fit an exploratory Hill curve to the published Nature 2024 source values.

This is a descriptive reanalysis of the five published glucose concentrations
and three independent Figure 3d experiments. It is not a binding model or a
test of a clinical glucose threshold.
"""

from __future__ import annotations

import csv
from datetime import datetime, timezone
import hashlib
import json
import platform
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import scipy
from scipy.optimize import least_squares


ROOT = Path(__file__).resolve().parent
INPUT = ROOT / "data" / "source-data.csv"
OUTPUT = ROOT / "data" / "hill-fit-summary.json"
PLOT = ROOT / "assets" / "hill-fit.svg"
BOOTSTRAP_REPLICATES = 2_000
BOOTSTRAP_SEED = 2215
EC50_UPPER_BOUND_MM = 1_000.0


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_source() -> tuple[list[dict[str, object]], np.ndarray, np.ndarray, list[int]]:
    with INPUT.open(newline="", encoding="utf-8") as handle:
        raw_rows = list(csv.DictReader(handle))
    rows = [row for row in raw_rows if row["analyte"] == "NNC2215"]
    experiments = sorted({int(row["experiment"]) for row in rows})
    glucose = np.array(sorted({float(row["glucose_mmol_L"]) for row in rows}))
    values = np.array(
        [
            [
                float(
                    next(
                        row["relative_affinity_percent"]
                        for row in rows
                        if int(row["experiment"]) == experiment
                        and float(row["glucose_mmol_L"]) == level
                    )
                )
                for level in glucose
            ]
            for experiment in experiments
        ],
        dtype=float,
    )
    if values.shape != (3, 5) or glucose[0] != 0:
        raise ValueError(
            f"Expected 3 experiments × 5 concentrations including 0 mM; got {values.shape} and {glucose.tolist()}"
        )
    return rows, glucose, values, experiments


def model_curve(glucose: np.ndarray, baselines: np.ndarray, tops: np.ndarray, ec50: float, hill: float) -> np.ndarray:
    positive = np.maximum(glucose, 0.0)
    powered_glucose = positive**hill
    powered_midpoint = ec50**hill
    fraction = powered_glucose / (powered_midpoint + powered_glucose)
    return baselines[:, None] + (tops[:, None] - baselines[:, None]) * fraction[None, :]


def fit_shared_curve(glucose: np.ndarray, values: np.ndarray):
    """Use each observed 0 mM point as baseline; share EC50 and slope."""
    baselines = values[:, 0].copy()
    experiment_count = values.shape[0]
    lower = np.r_[baselines + 1e-8, np.log(0.01), np.log(0.1)]
    upper = np.r_[np.full(experiment_count, 10_000.0), np.log(EC50_UPPER_BOUND_MM), np.log(10.0)]
    initial_tops = np.maximum(values[:, -1] * 1.5, baselines + 1.0)
    initial = np.r_[initial_tops, np.log(10.0), np.log(1.0)]

    def predict(parameters: np.ndarray) -> np.ndarray:
        tops = parameters[:experiment_count]
        ec50, hill = np.exp(parameters[-2:])
        return model_curve(glucose, baselines, tops, ec50, hill)

    result = least_squares(
        lambda parameters: (predict(parameters) - values).ravel(),
        initial,
        bounds=(lower, upper),
        max_nfev=50_000,
    )
    if not result.success:
        raise RuntimeError(f"Hill fit did not converge: {result.message}")
    tops = result.x[:experiment_count]
    ec50, hill = np.exp(result.x[-2:])
    fitted = predict(result.x)
    return result, tops, float(ec50), float(hill), fitted


def percentile_range(values: np.ndarray) -> dict[str, float]:
    lower, median, upper = np.percentile(values, [2.5, 50, 97.5])
    return {"lower_2_5_percent": float(lower), "median": float(median), "upper_97_5_percent": float(upper)}


def main() -> None:
    rows, glucose, values, experiments = load_source()
    result, tops, ec50, hill, fitted = fit_shared_curve(glucose, values)
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    plot_glucose = np.linspace(0.0, 20.0, 201)
    bootstrap_ec50: list[float] = []
    bootstrap_hill: list[float] = []
    bootstrap_curve: list[np.ndarray] = []
    failed_bootstraps = 0

    for _ in range(BOOTSTRAP_REPLICATES):
        sampled_indices = rng.integers(0, len(experiments), size=len(experiments))
        sampled_values = values[sampled_indices]
        try:
            _, sampled_tops, sampled_ec50, sampled_hill, _ = fit_shared_curve(glucose, sampled_values)
        except (RuntimeError, ValueError, FloatingPointError):
            failed_bootstraps += 1
            continue
        bootstrap_ec50.append(sampled_ec50)
        bootstrap_hill.append(sampled_hill)
        curve_by_experiment = model_curve(
            plot_glucose,
            sampled_values[:, 0],
            sampled_tops,
            sampled_ec50,
            sampled_hill,
        )
        bootstrap_curve.append(np.mean(curve_by_experiment, axis=0))

    if len(bootstrap_ec50) < BOOTSTRAP_REPLICATES * 0.95:
        raise RuntimeError(
            f"Only {len(bootstrap_ec50)} of {BOOTSTRAP_REPLICATES} bootstrap fits converged."
        )

    bootstrap_curve_array = np.asarray(bootstrap_curve)
    curve_lower, curve_median, curve_upper = np.percentile(
        bootstrap_curve_array, [2.5, 50, 97.5], axis=0
    )
    mean_fit = np.mean(fitted, axis=0)
    predicted_mean = np.mean(
        model_curve(plot_glucose, values[:, 0], tops, ec50, hill), axis=0
    )
    observed_mean = np.mean(values, axis=0)
    observed_sd = np.std(values, axis=0, ddof=1)
    residuals = fitted - values
    residual_rmse = float(np.sqrt(np.mean(residuals[:, 1:] ** 2)))
    ec50_range = percentile_range(np.asarray(bootstrap_ec50))
    hill_range = percentile_range(np.asarray(bootstrap_hill))
    ec50_upper_hits = int(
        np.count_nonzero(np.asarray(bootstrap_ec50) >= EC50_UPPER_BOUND_MM * (1.0 - 1e-6))
    )
    ec50_upper_hit_fraction = ec50_upper_hits / len(bootstrap_ec50)

    curve_rows = []
    for index, level in enumerate(plot_glucose):
        curve_rows.append(
            {
                "glucose_mmol_L": float(level),
                "fitted_mean": float(predicted_mean[index]),
                "bootstrap_median_mean": float(curve_median[index]),
                "bootstrap_lower_2_5_percent": float(curve_lower[index]),
                "bootstrap_upper_97_5_percent": float(curve_upper[index]),
            }
        )

    created_at = datetime.now(timezone.utc)
    point_responses = {
        str(experiment): float(top)
        for experiment, top in zip(experiments, tops, strict=True)
    }
    source_rows = [
        {
            "experiment": int(row["experiment"]),
            "glucose_mmol_L": float(row["glucose_mmol_L"]),
            "relative_affinity_percent": float(row["relative_affinity_percent"]),
            "source_sheet": row["source_sheet"],
            "source_cell": row["source_cell"],
        }
        for row in rows
    ]
    record = {
        "schema": "dia-health-one.nnc2215-hill-fit.v1",
        "run": {
            "run_id": "DHO-NNC2215-HILL-"
            + created_at.strftime("%Y%m%dT%H%M%SZ"),
            "created_at_utc": created_at.isoformat(),
            "script": "research/fit_nnc2215_hill_response.py",
            "script_sha256": sha256(Path(__file__).resolve()),
            "input": "research/data/source-data.csv",
            "input_sha256": sha256(INPUT),
            "engine": "SciPy least_squares and deterministic cluster bootstrap",
            "software": {
                "python": platform.python_version(),
                "numpy": np.__version__,
                "scipy": scipy.__version__,
                "matplotlib": matplotlib.__version__,
            },
            "bootstrap_seed": BOOTSTRAP_SEED,
            "bootstrap_requested": BOOTSTRAP_REPLICATES,
            "bootstrap_converged": len(bootstrap_ec50),
            "bootstrap_failed": failed_bootstraps,
        },
        "source": {
            "doi": "10.1038/s41586-024-08042-3",
            "article": "https://www.nature.com/articles/s41586-024-08042-3",
            "source_data": "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-024-08042-3/MediaObjects/41586_2024_8042_MOESM5_ESM.xlsx",
            "sheet": "Figure 3d",
            "assay": "hIR-A affinity relative to human insulin; 1.5% HSA",
            "analyte": "NNC2215",
            "experiments": experiments,
            "glucose_mmol_L": glucose.tolist(),
            "rows": source_rows,
            "values_are_source_workbook_values": True,
            "digitized_from_figure": False,
        },
        "model": {
            "equation": "R_i(G) = b_i + (T_i - b_i) * G^n / (EC50^n + G^n)",
            "baseline": "b_i is the measured 0 mmol/L value for experiment i; it is fixed, not fitted.",
            "upper_response": "T_i is fitted separately for each independent experiment.",
            "shared_parameters": ["EC50", "apparent Hill coefficient n"],
            "fit_method": "Unweighted nonlinear least squares on the 15 experiment-level source values; the zero-glucose values define the fixed baselines.",
            "parameter_bounds": {
                "upper_response_percent": [0.0, 10_000.0],
                "EC50_mmol_L": [0.01, 1_000.0],
                "apparent_hill_coefficient": [0.1, 10.0],
            },
        },
        "derived": {
            "point_estimate": {
                "ec50_mmol_L": ec50,
                "apparent_hill_coefficient": hill,
                "experiment_baselines_percent": {
                    str(experiment): float(baseline)
                    for experiment, baseline in zip(experiments, values[:, 0], strict=True)
                },
                "experiment_upper_responses_percent": point_responses,
                "rmse_nonzero_glucose_points_percent": residual_rmse,
                "residuals_percent": {
                    str(experiment): [float(value) for value in experiment_residuals]
                    for experiment, experiment_residuals in zip(experiments, residuals, strict=True)
                },
            },
            "cluster_bootstrap": {
                "method": "Resample the three independent experiments with replacement; refit experiment-specific upper responses and shared EC50/Hill slope for each sample.",
                "resamples": len(bootstrap_ec50),
                "percentile_ranges_not_formal_confidence_intervals": {
                    "ec50_mmol_L": ec50_range,
                    "apparent_hill_coefficient": hill_range,
                },
                "parameter_bound_diagnostics": {
                    "ec50_upper_bound_mmol_L": EC50_UPPER_BOUND_MM,
                    "ec50_upper_bound_hit_count": ec50_upper_hits,
                    "ec50_upper_bound_hit_fraction": ec50_upper_hit_fraction,
                    "upper_percentile_is_censored": ec50_upper_hit_fraction > 0,
                },
            },
            "observed_summary": [
                {
                    "glucose_mmol_L": float(level),
                    "mean_relative_affinity_percent": float(mean_value),
                    "sample_sd_percent": float(sd),
                    "n_independent_experiments": len(experiments),
                }
                for level, mean_value, sd in zip(glucose, observed_mean, observed_sd, strict=True)
            ],
            "curve": curve_rows,
        },
        "interpretation": {
            "status": "exploratory; EC50 weakly identified",
            "summary": "The fitted midpoint varies widely when the three independent experiments are resampled. Some fits hit the 1000 mmol/L EC50 optimization ceiling, so the upper bootstrap tail is censored and no reliable midpoint interval is identified.",
            "limits": [
                "Only three independent experiment clusters and five glucose concentrations are available.",
                "The observed data do not establish a high-glucose response plateau; EC50 and upper response are therefore strongly correlated.",
                "Bootstrap percentile ranges are fragile with three clusters and are not formal inferential confidence intervals.",
                "The apparent Hill coefficient is a descriptive steepness parameter, not a binding-site count or proof of cooperativity.",
                "No 6 mmol/L activation threshold is inferred; the secondary product-design target remains separate.",
                "This fit describes published receptor-affinity data; it is not a molecular simulation, clinical model, or dosing recommendation.",
            ],
        },
    }
    OUTPUT.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.labelcolor": "#5f5a66",
            "text.color": "#302541",
            "xtick.color": "#77737e",
            "ytick.color": "#77737e",
            "axes.edgecolor": "#d9d5de",
            "svg.fonttype": "none",
        }
    )
    figure, axis = plt.subplots(figsize=(9, 4.8), constrained_layout=True)
    figure.patch.set_facecolor("#ffffff")
    axis.set_facecolor("#ffffff")
    axis.fill_between(
        plot_glucose,
        curve_lower,
        curve_upper,
        color="#cde6aa",
        alpha=0.62,
        linewidth=0,
        label="95% experiment-bootstrap range",
    )
    axis.plot(plot_glucose, predicted_mean, color="#503a82", linewidth=2.2, label="Shared-slope Hill fit")
    axis.errorbar(
        glucose,
        observed_mean,
        yerr=observed_sd,
        fmt="o",
        markersize=5.2,
        color="#302541",
        ecolor="#6c53a2",
        elinewidth=1.1,
        capsize=3,
        zorder=4,
        label="Source mean ± sample SD",
    )
    replicate_colors = ["#a58fc7", "#8d78b4", "#c1b4d4"]
    for index, (experiment, color) in enumerate(zip(experiments, replicate_colors, strict=True)):
        axis.scatter(
            glucose + (index - 1) * 0.14,
            values[index],
            s=20,
            facecolors="white",
            edgecolors=color,
            linewidths=1.05,
            zorder=3,
        )
    axis.set_xlim(-0.6, 20.6)
    axis.set_ylim(bottom=0)
    axis.set_xticks(glucose)
    axis.set_xlabel("Glucose concentration (mmol/L)", labelpad=9)
    axis.set_ylabel("hIR-A affinity (% of human insulin)", labelpad=9)
    axis.grid(axis="y", color="#ece9ef", linewidth=0.8)
    axis.set_axisbelow(True)
    for spine in ("top", "right"):
        axis.spines[spine].set_visible(False)
    axis.spines["bottom"].set_color("#d9d5de")
    axis.spines["left"].set_color("#d9d5de")
    axis.legend(
        handles=[
            Line2D([0], [0], marker="o", linestyle="none", markerfacecolor="white", markeredgecolor="#a58fc7", label="3 independent experiments"),
            Line2D([0], [0], marker="o", color="#302541", markerfacecolor="#302541", label="Mean ± sample SD"),
            Line2D([0], [0], color="#503a82", linewidth=2.2, label="Exploratory Hill fit"),
            Line2D([0], [0], color="#cde6aa", linewidth=7, label="95% experiment-bootstrap range"),
        ],
        loc="upper left",
        frameon=False,
        ncol=2,
        fontsize=8,
        handlelength=1.7,
        columnspacing=1.2,
    )
    figure.savefig(PLOT, format="svg", facecolor=figure.get_facecolor())
    plt.close(figure)

    print(f"Wrote {OUTPUT.relative_to(ROOT.parent)}")
    print(f"Wrote {PLOT.relative_to(ROOT.parent)}")
    print(
        f"EC50={ec50:.3f} mmol/L (experiment-bootstrap 2.5th percentile "
        f"{ec50_range['lower_2_5_percent']:.3f}; 97.5th percentile "
        f"{ec50_range['upper_97_5_percent']:.3f}{' censored at fit ceiling' if ec50_upper_hit_fraction > 0 else ''}); "
        f"apparent Hill slope={hill:.3f} "
        f"({hill_range['lower_2_5_percent']:.3f}–{hill_range['upper_97_5_percent']:.3f})"
    )
    print(
        f"Input SHA-256: {sha256(INPUT)}; converged bootstrap fits: {len(bootstrap_ec50)}; "
        f"EC50 upper-bound hits: {ec50_upper_hits} ({ec50_upper_hit_fraction:.1%})"
    )


if __name__ == "__main__":
    main()
