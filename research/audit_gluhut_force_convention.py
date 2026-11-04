#!/usr/bin/env python3
"""Resolve AMBER/Grossfield force-constant conventions in the GluHUT deposit.

AMBER's distance-restraint parabolic region is written as

    E = rk2 * (R - R0)^2

whereas the conventional harmonic form is

    E = 1/2 * K * (R - R0)^2.

Therefore K = 2 * rk2.  Grossfield WHAM metadata uses the conventional K.
This audit compares the deposited SMD ``rk2`` values and the independently
reproduced WHAM metadata with the values reported in the GluHUT paper.

It resolves a scalar convention only.  It does not repair the inconsistent
sugar atom selections, reconstruct starting coordinates, validate source force
fields, or authorize PaG/MC-Q1 calculations.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

PAPER_DOI = "10.1039/D6SC00738D"
PAPER_SMD_REPORTED_KCAL_MOL_A2 = 200.0
PAPER_UMBRELLA_REPORTED_KCAL_MOL_A2 = 50.0
AMBER_MANUAL_ENERGY_FORM = "E = rk2 * (R - R0)^2"
STANDARD_HARMONIC_ENERGY_FORM = "E = 0.5 * K * (R - R0)^2"


class ForceConventionError(RuntimeError):
    pass


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def close(left: float, right: float, tolerance: float = 1.0e-12) -> bool:
    return math.isclose(left, right, rel_tol=0.0, abs_tol=tolerance)


def collect_smd_rk2(mapping: dict[str, Any]) -> list[float]:
    values: list[float] = []
    for row in mapping.get("mappings", []):
        value = row.get("parameters", {}).get("rk2")
        if value is None:
            raise ForceConventionError(
                f"missing rk2 for mapping {row.get('system_name', '<unknown>')}"
            )
        values.append(float(value))
    if not values:
        raise ForceConventionError("no SMD mapping rows found")
    return values


def energy_penalties(coefficient_rk2: float) -> list[dict[str, float]]:
    rows = []
    standard_k = 2.0 * coefficient_rk2
    for displacement in (0.1, 0.5, 1.0):
        amber_energy = coefficient_rk2 * displacement**2
        standard_energy = 0.5 * standard_k * displacement**2
        if not close(amber_energy, standard_energy):
            raise ForceConventionError("harmonic convention conversion is inconsistent")
        rows.append(
            {
                "displacement_A": displacement,
                "amber_energy_kcal_mol": amber_energy,
                "standard_harmonic_energy_kcal_mol": standard_energy,
            }
        )
    return rows


def audit(
    mapping_path: Path,
    wham_crosscheck_path: Path,
    output_path: Path,
) -> dict[str, Any]:
    mapping = json.loads(mapping_path.read_text(encoding="utf-8"))
    wham = json.loads(wham_crosscheck_path.read_text(encoding="utf-8"))

    smd_rk2_values = collect_smd_rk2(mapping)
    unique_smd_rk2 = sorted(set(smd_rk2_values))
    if len(unique_smd_rk2) != 1:
        raise ForceConventionError(
            f"deposited SMD mappings use multiple rk2 values: {unique_smd_rk2}"
        )
    deposited_smd_rk2 = unique_smd_rk2[0]
    deposited_wham_k = float(
        wham["implementation"]["force_constant_kcal_mol_A2"]
    )

    smd_standard_k = 2.0 * deposited_smd_rk2
    umbrella_expected_wham_k = 2.0 * PAPER_UMBRELLA_REPORTED_KCAL_MOL_A2

    report = {
        "schema": "dia-health-one.gluhut-force-convention.v1",
        "status": "scalar_force_convention_resolved_other_provenance_gates_open",
        "scope": "published and deposited GluHUT SMD/umbrella restraint coefficients",
        "not_target_claim": (
            "not NNC2215, PaG, MC-Q1, binding affinity, efficacy, safety, or a drug result"
        ),
        "inputs": {
            "smd_mapping_path": str(mapping_path),
            "smd_mapping_sha256": sha256(mapping_path),
            "wham_crosscheck_path": str(wham_crosscheck_path),
            "wham_crosscheck_sha256": sha256(wham_crosscheck_path),
            "paper_doi": PAPER_DOI,
        },
        "definitions": {
            "amber_parabolic_energy": AMBER_MANUAL_ENERGY_FORM,
            "standard_harmonic_energy": STANDARD_HARMONIC_ENERGY_FORM,
            "conversion": "K_standard = 2 * rk2_amber",
            "amber_distance_units": "rk2 in kcal mol^-1 A^-2",
            "grossfield_wham_metadata_convention": (
                "conventional harmonic K in 0.5*K*(R-R0)^2"
            ),
        },
        "smd": {
            "paper_reported_numeric_value_kcal_mol_A2": (
                PAPER_SMD_REPORTED_KCAL_MOL_A2
            ),
            "deposited_rk2_values_kcal_mol_A2": smd_rk2_values,
            "unique_deposited_rk2_kcal_mol_A2": deposited_smd_rk2,
            "standard_harmonic_K_kcal_mol_A2": smd_standard_k,
            "paper_numeric_value_matches_deposited_rk2": close(
                PAPER_SMD_REPORTED_KCAL_MOL_A2, deposited_smd_rk2
            ),
            "paper_value_if_interpreted_as_standard_K_matches": close(
                PAPER_SMD_REPORTED_KCAL_MOL_A2, smd_standard_k
            ),
            "energy_examples": energy_penalties(deposited_smd_rk2),
            "interpretation": (
                "The paper's SMD value 200 matches the deposited AMBER rk2 input. "
                "In the conventional 0.5*K*x^2 notation, the same potential has K=400."
            ),
        },
        "umbrella": {
            "paper_reported_numeric_value_kcal_mol_A2": (
                PAPER_UMBRELLA_REPORTED_KCAL_MOL_A2
            ),
            "expected_amber_rk2_kcal_mol_A2": (
                PAPER_UMBRELLA_REPORTED_KCAL_MOL_A2
            ),
            "expected_grossfield_wham_K_kcal_mol_A2": umbrella_expected_wham_k,
            "deposited_reproduction_wham_K_kcal_mol_A2": deposited_wham_k,
            "wham_K_matches_twice_paper_value": close(
                deposited_wham_k, umbrella_expected_wham_k
            ),
            "energy_examples": energy_penalties(
                PAPER_UMBRELLA_REPORTED_KCAL_MOL_A2
            ),
            "interpretation": (
                "The paper's umbrella value 50 is consistent with an AMBER rk2 "
                "coefficient of 50, while Grossfield WHAM requires K=100."
            ),
        },
        "decision": {
            "amber_to_standard_harmonic_conversion_resolved": True,
            "smd_paper_and_deposit_numerically_consistent_under_amber_convention": (
                close(PAPER_SMD_REPORTED_KCAL_MOL_A2, deposited_smd_rk2)
            ),
            "umbrella_paper_and_wham_metadata_consistent": close(
                deposited_wham_k, umbrella_expected_wham_k
            ),
            "force_constant_scalar_convention_validated": (
                close(PAPER_SMD_REPORTED_KCAL_MOL_A2, deposited_smd_rk2)
                and close(deposited_wham_k, umbrella_expected_wham_k)
            ),
            "paper_uses_unambiguous_standard_half_k_notation": False,
            "new_smd_or_umbrella_generation_authorized": False,
            "parent_pag_calibration_authorized": False,
            "candidate_ranking_authorized": False,
            "next_action": (
                "retain the resolved coefficient conversion, but obtain author "
                "clarification/corrected guest groups and source coordinates before "
                "prospective trajectory generation"
            ),
        },
        "limitations": [
            "The paper uses the phrase force constant without explicitly distinguishing the AMBER rk2 coefficient from conventional 0.5*K*x^2 notation.",
            "This audit resolves scalar equivalence only; it does not validate atom groups, mass weighting, target interpolation, trajectory generation, or starting coordinates.",
            "The deposited fructofuranose group inconsistency remains unresolved.",
            "Agreement of coefficients does not validate the host or guest force-field derivation.",
        ],
        "references": [
            {
                "source": "Eades et al., Chemical Science 2026",
                "doi": PAPER_DOI,
                "relevant_report": (
                    "SMD: 200 kcal mol^-1 A^-2; umbrella sampling: 50 kcal mol^-1 A^-2"
                ),
            },
            {
                "source": "AMBER users manual, NMR restraint functional form",
                "relevant_report": "parabolic energy is rk2*(R-r2)^2",
            },
            {
                "source": "AMBER users manual, umbrella sampling convention",
                "relevant_report": (
                    "the conventionally defined force constant is twice rk2"
                ),
            },
        ],
    }

    if not report["decision"]["force_constant_scalar_convention_validated"]:
        raise ForceConventionError(
            "paper/deposit force coefficients do not satisfy the AMBER/Grossfield conversion"
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return report


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--mapping-json",
        type=Path,
        default=Path("research/data/gluhut-smd-restraint-mapping.json"),
    )
    parser.add_argument(
        "--wham-json",
        type=Path,
        default=Path("research/data/gluhut-glucose-numpy-wham-crosscheck.json"),
    )
    parser.add_argument(
        "--output-json",
        type=Path,
        default=Path("research/data/gluhut-force-convention.json"),
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    report = audit(args.mapping_json, args.wham_json, args.output_json)
    print(
        json.dumps(
            {
                "smd": report["smd"],
                "umbrella": report["umbrella"],
                "decision": report["decision"],
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
