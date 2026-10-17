#!/usr/bin/env python3
"""Recalculate the published NNC2215 Figure 3d benchmark from source values."""

import csv
import hashlib
import json
import platform
import statistics
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data" / "source-data.csv"
RESULTS = ROOT / "data" / "reproduction.json"
CHART = ROOT / "assets" / "glucose-response.svg"
SOURCE = "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-024-08042-3/MediaObjects/41586_2024_8042_MOESM5_ESM.xlsx"
WORKBOOK_SHA256 = "09a9b6aa4103f9d61af7df03886d59f16d6b09f69f2fe768260548061c7dc81d"


def fmean(values):
    return statistics.mean(values)


def sample_sd(values):
    return statistics.stdev(values) if len(values) > 1 else 0.0


def load_data():
    with DATA.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        row["experiment"] = int(row["experiment"])
        row["glucose_mmol_L"] = float(row["glucose_mmol_L"])
        row["relative_affinity_percent"] = float(row["relative_affinity_percent"])
    return rows


def summarize(rows, analyte):
    series = []
    for glucose in (0, 3, 5, 10, 20):
        values = [r["relative_affinity_percent"] for r in rows
                  if r["analyte"] == analyte and r["glucose_mmol_L"] == glucose]
        series.append({"glucose_mmol_L": glucose, "values": values,
                       "n": len(values), "mean": fmean(values), "sd": sample_sd(values)})
    return series


def replicate_fold(rows, low, high):
    paired = {}
    for r in rows:
        if r["analyte"] == "NNC2215" and r["glucose_mmol_L"] in (low, high):
            paired.setdefault(r["experiment"], {})[r["glucose_mmol_L"]] = r["relative_affinity_percent"]
    ratios = [v[high] / v[low] for _, v in sorted(paired.items())]
    mean = fmean(ratios)
    group_means = summarize(rows, "NNC2215")
    means = {x["glucose_mmol_L"]: x["mean"] for x in group_means}
    return {"from_glucose_mmol_L": low, "to_glucose_mmol_L": high,
            "replicate_ratios": ratios, "mean_within_experiment_ratio": mean,
            "sd_within_experiment_ratio": sample_sd(ratios),
            "ratio_of_group_means": means[high] / means[low]}


def svg_chart(nnc):
    width, height = 900, 480
    left, right, top, bottom = 84, 856, 54, 374
    y_max = 11.5
    x = lambda g: left + (g / 20) * (right - left)
    y = lambda v: bottom - (v / y_max) * (bottom - top)
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" role="img" aria-labelledby="title desc">',
           '<title id="title">NNC2215 insulin receptor affinity across glucose concentrations</title>',
           '<desc id="desc">Three independent experiment values and mean with sample standard deviation at five glucose concentrations. Insulin degludec values are listed in the adjacent control summary.</desc>',
           '<rect width="100%" height="100%" fill="#fbfaf8"/>',
           '<g font-family="Arial, sans-serif" fill="#292633">',
           '<text x="84" y="27" font-size="12" letter-spacing="2" fill="#706b7c">PUBLISHED SOURCE DATA · FIGURE 3D</text>']
    for tick in (0, 2, 4, 6, 8, 10):
        yy = y(tick)
        out.append(f'<line x1="{left}" y1="{yy:.1f}" x2="{right}" y2="{yy:.1f}" stroke="#e9e6e1" stroke-width="1"/>')
        out.append(f'<text x="{left-17}" y="{yy+4:.1f}" text-anchor="end" font-size="12" fill="#817d88">{tick}</text>')
    for glucose in (0, 3, 5, 10, 20):
        xx = x(glucose)
        out.append(f'<line x1="{xx:.1f}" y1="{top}" x2="{xx:.1f}" y2="{bottom}" stroke="#eeece8" stroke-width="1"/>')
        out.append(f'<text x="{xx:.1f}" y="{bottom+25}" text-anchor="middle" font-size="12" fill="#817d88">{glucose}</text>')
    out.append(f'<text x="{(left+right)/2:.1f}" y="{height-20}" text-anchor="middle" font-size="12" fill="#6c6875">D-glucose concentration (mmol/L)</text>')
    out.append(f'<text transform="translate(22 {(top+bottom)/2:.1f}) rotate(-90)" text-anchor="middle" font-size="12" fill="#6c6875">Relative hIR-A affinity (% of human insulin)</text>')
    points = " ".join(f'{x(d["glucose_mmol_L"]):.1f},{y(d["mean"]):.1f}' for d in nnc)
    out.append(f'<polyline points="{points}" fill="none" stroke="#65509b" stroke-width="3" stroke-linecap="round" stroke-linejoin="round"/>')
    for d in nnc:
        xx, yy = x(d["glucose_mmol_L"]), y(d["mean"])
        half = 3.8
        for i, val in enumerate(d["values"]):
            jitter = (i - 1) * 7
            out.append(f'<circle cx="{xx+jitter:.1f}" cy="{y(val):.1f}" r="4.2" fill="#c7bce0" stroke="#fff" stroke-width="1.4"><title>Experiment {i+1}: {val:.4f}%</title></circle>')
        out.append(f'<line x1="{xx:.1f}" y1="{y(d["mean"]+d["sd"]):.1f}" x2="{xx:.1f}" y2="{y(d["mean"]-d["sd"]):.1f}" stroke="#49367a" stroke-width="1.5"/>')
        out.append(f'<line x1="{xx-half:.1f}" y1="{y(d["mean"]+d["sd"]):.1f}" x2="{xx+half:.1f}" y2="{y(d["mean"]+d["sd"]):.1f}" stroke="#49367a" stroke-width="1.5"/>')
        out.append(f'<line x1="{xx-half:.1f}" y1="{y(d["mean"]-d["sd"]):.1f}" x2="{xx+half:.1f}" y2="{y(d["mean"]-d["sd"]):.1f}" stroke="#49367a" stroke-width="1.5"/>')
        out.append(f'<circle cx="{xx:.1f}" cy="{yy:.1f}" r="6.3" fill="#65509b" stroke="#fff" stroke-width="2"/>')
    out.extend(['<circle cx="84" cy="431" r="4" fill="#c7bce0"/>',
                '<text x="96" y="435" font-size="11" fill="#6c6875">Independent experiments (n=3)</text>',
                '<line x1="326" y1="431" x2="354" y2="431" stroke="#49367a" stroke-width="1.5"/>',
                '<circle cx="340" cy="431" r="5.5" fill="#65509b" stroke="#fff" stroke-width="1.5"/>',
                '<text x="365" y="435" font-size="11" fill="#6c6875">Mean ± sample SD</text>',
                '</g></svg>'])
    CHART.parent.mkdir(parents=True, exist_ok=True)
    CHART.write_text("\n".join(out), encoding="utf-8")


def main():
    rows = load_data()
    sha = hashlib.sha256(DATA.read_bytes()).hexdigest()
    nnc = summarize(rows, "NNC2215")
    control = summarize(rows, "Insulin degludec")
    folds = [replicate_fold(rows, 3, 20), replicate_fold(rows, 0, 20)]
    folds[0]["published_reported_fold"] = 3.2
    folds[1]["published_reported_fold"] = 12.5
    for item in folds:
        item["matches_published_rounding"] = round(item["mean_within_experiment_ratio"], 1) == item["published_reported_fold"]
    result = {
        "schema": "dia-health-one.reproduction.v2",
        "run": {"id": "DHO-NNC2215-FIG3D-" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ"),
                "created_at_utc": datetime.now(timezone.utc).isoformat(),
                "engine": "Python standard library", "engine_version": sys.version.split()[0],
                "analysis": "source-data reanalysis; means, sample SD, paired experiment fold changes",
                "input_sha256": sha,
                "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                "environment": {"platform": platform.platform(), "architecture": platform.machine(),
                                "dependencies": "Python standard library only", "accelerator": "CPU", "gpu_simulation": False},
                "clinical_use": False},
        "source": {"doi": "10.1038/s41586-024-08042-3", "article": "https://www.nature.com/articles/s41586-024-08042-3",
                   "source_data": SOURCE, "source_data_sha256": WORKBOOK_SHA256, "accessed_utc_date": "2026-10-17",
                   "sheet": "Figure 3d", "assay": "hIR-A affinity relative to human insulin; 1.5% HSA",
                   "replicates": "Three independent experiments; source workbook values, not digitised chart estimates"},
        "inputs": rows,
        "series": {"NNC2215": nnc, "Insulin_degludec": control},
        "fold_changes": folds,
        "interpretation": {"hill_fit": "A separate exploratory shared-slope Hill fit is recorded in data/hill-fit-summary.json. Its experiment-bootstrap EC50 upper tail is censored at the 1000 mmol/L fit ceiling; no reliable midpoint or 6 mmol/L threshold is inferred.",
                           "scope": "Reanalysis of published in-vitro data only. No molecular dynamics, binding free-energy, or patient-response model was run."}
    }
    RESULTS.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    svg_chart(nnc)
    print(json.dumps({"run_id": result["run"]["id"], "input_sha256": sha,
                      "fold_changes": folds, "means_percent": [round(x["mean"], 4) for x in nnc],
                      "sample_sd_percent": [round(x["sd"], 4) for x in nnc],
                      "output": str(RESULTS.relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    main()
