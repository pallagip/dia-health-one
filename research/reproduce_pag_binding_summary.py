#!/usr/bin/env python3
"""Recalculate ratios from rounded, published PaG binding summaries."""

from __future__ import annotations

import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
INPUT = ROOT / "data" / "pag-binding-summary.csv"
OUTPUT = ROOT / "data" / "pag-binding-summary-reproduction.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    with INPUT.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    kd = {
        row["condition"]: float(row["reported_value"])
        for row in rows
        if row["measure"] == "Kd"
    }
    affinity = {
        row["condition"]: float(row["reported_value"])
        for row in rows
        if row["measure"] == "Relative insulin-receptor affinity"
    }

    glucose_kd = kd["D-glucose binding to NNC2215"]
    pag_kd = kd["PaG binding to NNC2215"]
    baseline_affinity = affinity["zero d-glucose without PaG"]
    pag_affinity = affinity["zero d-glucose with 70 umol/L PaG"]

    result = {
        "record_type": "published_pag_binding_summary_reanalysis",
        "schema": "dia-health-one.pag-binding-summary-reproduction.v1",
        "run": {
            "run_id": "DHO-PAG-BINDING-SUMMARY-"
            + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ"),
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
            "script": "research/reproduce_pag_binding_summary.py",
            "script_sha256": sha256(Path(__file__).resolve()),
            "input": "research/data/pag-binding-summary.csv",
            "input_sha256": sha256(INPUT),
        },
        "source": {
            "doi": "10.1021/acsptsci.6c00069",
            "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC13270365/",
            "basis": "Published Kd estimates and rounded zero-glucose receptor-affinity means transcribed from the article.",
        },
        "reported_values": rows,
        "derived": {
            "pag_lower_kd_fold_vs_glucose": {
                "value": glucose_kd / pag_kd,
                "equation": "Kd_glucose / Kd_PaG",
                "unit": "fold",
                "interpretation": "Ratio of reported Kd means; a lower Kd indicates stronger binding.",
            },
            "pag_receptor_affinity_fold_at_zero_glucose": {
                "value": pag_affinity / baseline_affinity,
                "equation": "relative_affinity_70uM_PaG / relative_affinity_no_PaG",
                "unit": "fold",
                "interpretation": "Ratio of rounded, published hIR-A relative-affinity means.",
            },
            "pag_receptor_affinity_change_at_zero_glucose": {
                "value": pag_affinity - baseline_affinity,
                "equation": "relative_affinity_70uM_PaG - relative_affinity_no_PaG",
                "unit": "percentage points of human-insulin relative affinity",
            },
        },
        "limits": [
            "These calculations redivide published summary values; they do not fit the native-MS or receptor-binding data.",
            "The receptor-affinity inputs are rounded means. No uncertainty is calculated because the underlying replicate values and covariance are not transcribed here.",
            "This is an in-vitro summary reanalysis, not a molecular simulation or patient-specific prediction.",
        ],
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {OUTPUT.relative_to(ROOT.parent)}")
    print(
        "PaG lower Kd vs glucose: "
        f"{glucose_kd / pag_kd:.3f}x; zero-glucose receptor-affinity ratio: "
        f"{pag_affinity / baseline_affinity:.3f}x"
    )


if __name__ == "__main__":
    main()
