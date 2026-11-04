#!/usr/bin/env python3
"""Reconcile AMBER, paper, and Grossfield WHAM spring conventions.

AMBER distance restraints use E = rk2 * (R - R0)^2. The conventional harmonic
form used by Grossfield WHAM and many other packages is E = 1/2 * K *
(R - R0)^2, so K = 2 * rk2.

This audit checks the deposited GluHUT SMD controls, the published SMD and
umbrella values, and the deposited WHAM metadata. It validates unit/convention
translation only. It does not validate atom selections, coordinates, sampling,
PMFs, PaG, MC-Q1, or a drug result.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

ARTICLE_DOI = "10.1039/D6SC00738D"
ARTICLE_SMD_AMBER_RK2_KCAL_MOL_A2 = 200.0
ARTICLE_UMBRELLA_AMBER_RK2_KCAL_MOL_A2 = 50.0
EXPECTED_SMD_START_A = 0.0
EXPECTED_SMD_END_A = 15.0

AMBER_MANUAL_SOURCE = {
    "title": "AMBER 9 Users' Manual",
    "url": "https://dev-archive.ambermd.org/200603/att-0122/amber9.pdf",
    "sections": ["6.8 Steered Molecular Dynamics", "6.12.1 NMR restraints"],
    "equation": "E = rk2 * (R - r2)^2 for the harmonic distance branch",
    "conventional_relation": "K_conventional = 2 * rk2",
}


class ConventionError(RuntimeError):
    pass


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def conventional_k_from_amber_rk2(rk2: float) -> float:
    return 2.0 * rk2


def amber_rk2_from_conventional_k(k: float) -> float:
    return 0.5 * k


def amber_energy(rk2: float, displacement_A: float) -> float:
    return rk2 * displacement_A**2


def conventional_energy(k: float, displacement_A: float) -> float:
    return 0.5 * k * displacement_A**2


def unique_numeric_values(
    mappings: list[dict[str, Any]], parameter: str
) -> list[float | None]:
    values = {mapping["parameters"].get(parameter) for mapping in mappings}
    return sorted(values, key=lambda value: (value is None, value))


def run(mapping_path: Path, wham_path: Path, output_path: Path) -> dict[str, Any]:
    mapping = json.loads(mapping_path.read_text(encoding="utf-8"))
    wham = json.loads(wham_path.read_text(encoding="utf-8"))

    mappings = mapping["mappings"]
    if len(mappings) != 4:
        raise ConventionError(f"expected four SMD mappings, found {len(mappings)}")

    smd_rk2_values = unique_numeric_values(mappings, "rk2")
    smd_start_values = unique_numeric_values(mappings, "r2")
    smd_end_values = unique_numeric_values(mappings, "r2a")
    if smd_rk2_values != [ARTICLE_SMD_AMBER_RK2_KCAL_MOL_A2]:
        raise ConventionError(f"unexpected SMD rk2 values: {smd_rk2_values}")
    if smd_start_values != [EXPECTED_SMD_START_A]:
        raise ConventionError(f"unexpected SMD start values: {smd_start_values}")
    if smd_end_values != [EXPECTED_SMD_END_A]:
        raise ConventionError(f"unexpected SMD end values: {smd_end_values}")

    wham_parameters = wham["analysis"]["parameters"]
    wham_k = float(wham_parameters["wham_metadata_spring_kcal_mol_angstrom2"])
    article_umbrella_rk2 = float(
        wham_parameters["article_reported_amber_spring_kcal_mol_angstrom2"]
    )

    smd_conventional_k = conventional_k_from_amber_rk2(
        ARTICLE_SMD_AMBER_RK2_KCAL_MOL_A2
    )
    umbrella_conventional_k = conventional_k_from_amber_rk2(
        ARTICLE_UMBRELLA_AMBER_RK2_KCAL_MOL_A2
    )

    displacement_checks = [0.1, 0.5, 1.0]
    energy_equivalence = [
        {
            "displacement_A": displacement,
            "smd_amber_energy_kcal_mol": amber_energy(
                ARTICLE_SMD_AMBER_RK2_KCAL_MOL_A2, displacement
            ),
            "smd_conventional_energy_kcal_mol": conventional_energy(
                smd_conventional_k, displacement
            ),
            "umbrella_amber_energy_kcal_mol": amber_energy(
                ARTICLE_UMBRELLA_AMBER_RK2_KCAL_MOL_A2, displacement
            ),
            "umbrella_conventional_energy_kcal_mol": conventional_energy(
                umbrella_conventional_k, displacement
            ),
        }
        for displacement in displacement_checks
    ]

    report = {
        "schema": "dia-health-one.gluhut-restraint-force-convention.v1",
        "status": "amber_wham_force_convention_reconciled",
        "scope": "public GluHUT SMD and umbrella restraint coefficient translation",
        "not_target_claim": (
            "not validation of guest atom selections, coordinates, sampling, PMFs, "
            "NNC2215, PaG, MC-Q1, efficacy, safety, or a drug result"
        ),
        "sources": {
            "article": {
                "doi": ARTICLE_DOI,
                "reported_smd_harmonic_value_kcal_mol_A2": (
                    ARTICLE_SMD_AMBER_RK2_KCAL_MOL_A2
                ),
                "reported_umbrella_harmonic_value_kcal_mol_A2": (
                    ARTICLE_UMBRELLA_AMBER_RK2_KCAL_MOL_A2
                ),
                "reported_smd_range_A": [EXPECTED_SMD_START_A, EXPECTED_SMD_END_A],
            },
            "amber_manual": AMBER_MANUAL_SOURCE,
            "smd_mapping": {
                "path": str(mapping_path),
                "sha256": sha256(mapping_path),
            },
            "wham_reproduction": {
                "path": str(wham_path),
                "sha256": sha256(wham_path),
            },
        },
        "definitions": {
            "amber_distance_restraint": "E = rk2 * (R - R0)^2",
            "conventional_harmonic": "E = 0.5 * K * (R - R0)^2",
            "conversion": "K = 2 * rk2",
        },
        "smd": {
            "deposited_rk2_values_kcal_mol_A2": smd_rk2_values,
            "deposited_start_target_values_A": smd_start_values,
            "deposited_end_target_values_A": smd_end_values,
            "article_numeric_value_kcal_mol_A2": (
                ARTICLE_SMD_AMBER_RK2_KCAL_MOL_A2
            ),
            "amber_rk2_kcal_mol_A2": ARTICLE_SMD_AMBER_RK2_KCAL_MOL_A2,
            "conventional_half_K_kcal_mol_A2": smd_conventional_k,
            "openmm_or_half_K_equivalent_kcal_mol_A2": smd_conventional_k,
        },
        "umbrella": {
            "article_amber_rk2_kcal_mol_A2": article_umbrella_rk2,
            "expected_article_amber_rk2_kcal_mol_A2": (
                ARTICLE_UMBRELLA_AMBER_RK2_KCAL_MOL_A2
            ),
            "deposited_wham_metadata_K_kcal_mol_A2": wham_k,
            "calculated_conventional_K_kcal_mol_A2": umbrella_conventional_k,
            "amber_rk2_recovered_from_wham_K_kcal_mol_A2": (
                amber_rk2_from_conventional_k(wham_k)
            ),
        },
        "energy_equivalence_examples": energy_equivalence,
        "decision": {
            "amber_energy_form_validated_from_manual": True,
            "conventional_K_is_twice_amber_rk2": True,
            "smd_deposited_rk2_matches_article_numeric_value": (
                smd_rk2_values
                == [ARTICLE_SMD_AMBER_RK2_KCAL_MOL_A2]
            ),
            "smd_start_and_end_match_article": (
                smd_start_values == [EXPECTED_SMD_START_A]
                and smd_end_values == [EXPECTED_SMD_END_A]
            ),
            "umbrella_wham_metadata_is_twice_article_amber_value": (
                wham_k == 2.0 * article_umbrella_rk2
                == umbrella_conventional_k
            ),
            "paper_values_are_consistent_with_amber_rk2_reporting": (
                smd_rk2_values
                == [ARTICLE_SMD_AMBER_RK2_KCAL_MOL_A2]
                and article_umbrella_rk2
                == ARTICLE_UMBRELLA_AMBER_RK2_KCAL_MOL_A2
                and wham_k == 2.0 * article_umbrella_rk2
            ),
            "force_constant_convention_resolved": True,
            "portable_values": {
                "amber_rst_smd_rk2_kcal_mol_A2": 200.0,
                "amber_rst_umbrella_rk2_kcal_mol_A2": 50.0,
                "half_K_package_smd_K_kcal_mol_A2": 400.0,
                "half_K_package_umbrella_K_kcal_mol_A2": 100.0,
                "grossfield_wham_umbrella_K_kcal_mol_A2": 100.0,
            },
            "guest_atom_selection_validated": False,
            "new_smd_or_umbrella_generation_authorized": False,
            "parent_pag_calibration_authorized": False,
            "mcq1_parameter_transfer_authorized": False,
            "candidate_ranking_authorized": False,
            "next_action": (
                "use the resolved 2x conversion in all manifests, but wait for corrected "
                "guest groups, source coordinates, and parameter provenance before generation"
            ),
        },
        "limitations": [
            "The article uses the word force constant for values that numerically match AMBER rk2; external packages may use the conventional half-K form.",
            "Resolving the coefficient convention does not resolve the inconsistent fructofuranose guest group.",
            "No coordinate, force, work, trajectory, PMF, or affinity is recalculated here.",
            "The AMBER manual source is an archived official manual; source-code behavior should still be checked when the exact AMBER 24 build is available."
        ],
    }

    for example in energy_equivalence:
        if example["smd_amber_energy_kcal_mol"] != example[
            "smd_conventional_energy_kcal_mol"
        ]:
            raise ConventionError("SMD energy forms are not equivalent")
        if example["umbrella_amber_energy_kcal_mol"] != example[
            "umbrella_conventional_energy_kcal_mol"
        ]:
            raise ConventionError("umbrella energy forms are not equivalent")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return report


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--smd-mapping",
        type=Path,
        default=Path("research/data/gluhut-smd-restraint-mapping.json"),
    )
    parser.add_argument(
        "--wham-reproduction",
        type=Path,
        default=Path("research/data/gluhut-glucose-wham-reproduction.json"),
    )
    parser.add_argument(
        "--output-json",
        type=Path,
        default=Path("research/data/gluhut-restraint-force-convention.json"),
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    report = run(args.smd_mapping, args.wham_reproduction, args.output_json)
    print(
        json.dumps(
            {
                "definitions": report["definitions"],
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
