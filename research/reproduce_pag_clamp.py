#!/usr/bin/env python3
"""Recalculate descriptive ratios from published pig-clamp group summaries.

This does not reanalyse individual animals or reproduce the paper's statistical
tests. The input contains the published group-level AUC means and SEM values.
"""

from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parent
INPUT = ROOT / "data" / "pag-clamp-summary.csv"
OUTPUT = ROOT / "data" / "pag-clamp-reproduction.json"


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def main() -> None:
    input_bytes = INPUT.read_bytes()
    with INPUT.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    controls = [row for row in rows if row["treatment"] == "Saline"]
    if len(controls) != 1:
        raise ValueError("Expected exactly one saline control row.")

    control = controls[0]
    control_mean = float(control["auc_gir_300_420_mg_kg"])
    comparisons = []
    for row in rows:
        if row["treatment"] == "Saline":
            continue
        mean_auc = float(row["auc_gir_300_420_mg_kg"])
        comparisons.append(
            {
                "treatment": row["treatment"],
                "dose": row["dose"],
                "auc_mean_mg_kg": mean_auc,
                "sem_mg_kg": float(row["sem_mg_kg"]),
                "n": int(row["n"]),
                "published_p_value_vs_saline": float(row["p_value_vs_saline"]),
                "derived_ratio_of_group_means_vs_saline": mean_auc / control_mean,
            }
        )

    script_bytes = Path(__file__).read_bytes()
    result = {
        "schema": "dia-health-one.pag-clamp-summary.v1",
        "run": {
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
            "analysis": "descriptive ratios of published group-level AUC means",
            "input_sha256": sha256(input_bytes),
            "script_sha256": sha256(script_bytes),
            "individual_animal_data_used": False,
            "statistical_tests_recomputed": False,
        },
        "source": {
            "doi": "10.1021/acsptsci.6c00069",
            "article": "https://pmc.ncbi.nlm.nih.gov/articles/PMC13270365/",
            "source_location": "Main article Table 1; clamp at 3 mmol/L glucose",
            "endpoint": "GIR AUC from 300 to 420 min, mg/kg",
            "reported_uncertainty": "SEM",
        },
        "published_group_summaries": rows,
        "saline_control": {
            "auc_mean_mg_kg": control_mean,
            "sem_mg_kg": float(control["sem_mg_kg"]),
            "n": int(control["n"]),
        },
        "derived_comparisons": comparisons,
        "interpretation": (
            "Ratios divide published group means and are descriptive. They are not "
            "paired animal-level effects and do not reproduce the published p-values."
        ),
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {OUTPUT.relative_to(ROOT)}")
    print(f"Input SHA-256: {result['run']['input_sha256']}")
    for comparison in comparisons:
        print(
            f"{comparison['treatment']}: "
            f"{comparison['derived_ratio_of_group_means_vs_saline']:.3f}x vs saline "
            f"(published p={comparison['published_p_value_vs_saline']:.2f})"
        )


if __name__ == "__main__":
    main()
