#!/usr/bin/env python3
"""Apply the published linear PaG coefficient to reported concentration anchors."""

from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parent
SCENARIOS = ROOT / "data" / "pag-cmax-scenarios.csv"
PARAMETERS = ROOT / "data" / "pag-model-parameters.csv"
OUTPUT = ROOT / "data" / "pag-sensitivity-reproduction.json"


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def main() -> None:
    scenario_bytes = SCENARIOS.read_bytes()
    parameter_bytes = PARAMETERS.read_bytes()
    with SCENARIOS.open(newline="", encoding="utf-8") as handle:
        scenarios = list(csv.DictReader(handle))
    with PARAMETERS.open(newline="", encoding="utf-8") as handle:
        parameters = list(csv.DictReader(handle))

    pag_coefficients = [row for row in parameters if row["parameter"] == "CoeffPaG"]
    if len(pag_coefficients) != 1:
        raise ValueError("Expected exactly one CoeffPaG parameter row.")
    coefficient = float(pag_coefficients[0]["estimate"])

    results = []
    for scenario in scenarios:
        concentration = float(scenario["pa_g_cmax_umol_l"])
        relative_increase_percent = concentration * coefficient
        results.append(
            {
                "scenario": scenario["scenario"],
                "concentration_umol_l": concentration,
                "source_class": scenario["source_class"],
                "source_context": scenario["source_context"],
                "source_location": scenario["source_location"],
                "model_implied_si_increase_percent": relative_increase_percent,
                "model_implied_si_multiplier": 1 + relative_increase_percent / 100,
            }
        )

    script_bytes = Path(__file__).read_bytes()
    result = {
        "schema": "dia-health-one.pag-sensitivity.v1",
        "run": {
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
            "analysis": "published Equation 4b applied to published PaG concentration anchors",
            "scenarios_sha256": sha256(scenario_bytes),
            "parameters_sha256": sha256(parameter_bytes),
            "script_sha256": sha256(script_bytes),
            "compute": "CPU arithmetic",
            "clinical_use": False,
        },
        "source": {
            "doi": "10.1021/acsptsci.6c00069",
            "article": "https://pmc.ncbi.nlm.nih.gov/articles/PMC13270365/",
            "model_equation": "SI = SIbase * (1 + ConcPaG * CoeffPaG / 100), Appendix S4 Equation 4b",
            "coefficient": coefficient,
            "coefficient_unit": "% increase from SIbase per umol/L PaG",
            "coefficient_estimate_cv_percent": float(pag_coefficients[0]["estimate_cv_percent"]),
            "coefficient_source": "Supporting Information Table S2",
        },
        "results": results,
        "interpretation": (
            "These values are the fitted model's linear SI component evaluated at published "
            "PaG concentration anchors. They are not receptor-affinity measurements, GIR "
            "predictions, or individual-patient outcomes."
        ),
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {OUTPUT.relative_to(ROOT)}")
    print(f"CoeffPaG: {coefficient:.4f}% per umol/L PaG")
    for row in results:
        print(
            f"{row['scenario']}: {row['model_implied_si_increase_percent']:.1f}% "
            f"SI increase ({row['model_implied_si_multiplier']:.3f}x baseline)"
        )


if __name__ == "__main__":
    main()
