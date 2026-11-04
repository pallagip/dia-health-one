#!/usr/bin/env python3
"""Merge the compiled-charge symmetry result into public status.json.

The updater preserves all existing status fields and only adds bounded claims
that are present in the deterministic charge-symmetry output. It intentionally
does not change MC-Q1's evidence level or enable PaG/candidate computation.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

STATUS_SCHEMA = "dia-health-one.public-confidence-status.v10"


class StatusUpdateError(RuntimeError):
    pass


def update_status(status: dict[str, Any], audit: dict[str, Any]) -> dict[str, Any]:
    decision = audit["decision"]
    graph = audit["graph_symmetry"]
    summary = audit["summary"]
    cross_system = audit["cross_system"]

    if decision["source_resp_derivation_reconstructable"] is not False:
        raise StatusUpdateError("source RESP derivation must remain unavailable")
    if decision["direct_transfer_to_pag_or_mcq1_authorized"] is not False:
        raise StatusUpdateError("direct parameter transfer must remain disabled")
    if decision["candidate_ranking_authorized"] is not False:
        raise StatusUpdateError("candidate ranking must remain disabled")
    if cross_system["host_charge_vectors_identical"] is not True:
        raise StatusUpdateError("host charge vectors differ across systems")

    status["schema"] = STATUS_SCHEMA
    status["updated"] = "2026-11-04"
    benchmark = status["public_benchmark_calibration"]
    benchmark.update(
        {
            "compiled_host_charge_vectors_identical_across_four_sugars": True,
            "compiled_host_charge_vector_sha256": cross_system[
                "host_charge_vector_sha256"
            ],
            "parameter_labeled_heavy_graph_automorphism_count": graph[
                "automorphism_count"
            ],
            "parameter_labeled_heavy_graph_orbit_count": graph[
                "heavy_orbit_count"
            ],
            "parameter_labeled_heavy_graph_nontrivial_orbit_count": graph[
                "nontrivial_heavy_orbit_count"
            ],
            "parameter_labeled_heavy_graph_orbit_size_distribution": graph[
                "heavy_orbit_size_distribution"
            ],
            "maximum_graph_equivalent_heavy_charge_range_e": summary[
                "maximum_nontrivial_heavy_orbit_charge_range_e"
            ],
            "maximum_graph_equivalent_hydrogen_charge_range_e": summary[
                "maximum_nontrivial_hydrogen_group_charge_range_e"
            ],
            "compiled_charges_exactly_equal_within_preserved_graph_symmetries": decision[
                "compiled_charges_are_exactly_graph_symmetry_averaged"
            ],
            "large_graph_equivalent_charge_asymmetry_detected": decision[
                "large_graph_equivalent_charge_asymmetry_detected"
            ],
            "symmetry_averaging_rms_per_atom_adjustment_e": summary[
                "symmetry_averaging_adjustment"
            ]["rms_per_atom_adjustment_e"],
            "automatic_charge_symmetrisation_authorized": decision[
                "automatic_charge_symmetrisation_authorized"
            ],
            "direct_compiled_charge_transfer_to_pag_or_mcq1_authorized": decision[
                "direct_transfer_to_pag_or_mcq1_authorized"
            ],
            "source_resp_equivalencing_constraints_available": False,
            "charge_symmetry_record": "../GLUHUT_CHARGE_SYMMETRY.md",
        }
    )

    status["aws_execution"]["current_gate"] = (
        "AWS-1 public analysis reproduced; host architecture, force convention, "
        "all deposited coordinates, endpoint/frame provenance, and compiled-charge "
        "symmetry audited; source RESP derivation, complete trajectory provenance, "
        "start coordinates and parent PaG calibration remain pending"
    )
    status["aws_execution"]["next_action"] = (
        "obtain source RESP conformers, ESP grids, equivalencing constraints and final "
        "mol2, plus full trajectory/start-coordinate provenance; then calibrate the "
        "parent glucose-versus-PaG system"
    )

    next_changes = status["next_confidence_changes"]
    resp_item = (
        "obtain source RESP conformers, ESP grids, equivalencing constraints, final "
        "mol2 and GAFF/GLYCAM assembly provenance"
    )
    if resp_item not in next_changes:
        next_changes.insert(0, resp_item)

    # Preserve the central evidence boundary even if an earlier status file drifts.
    status["primary_candidate"]["evidence_level"] = "E2"
    status["primary_candidate"]["candidate_ranking_authorized"] = False
    status["working_drug_discovered"] = False
    status["measured_dia_health_one_hit"] = False
    status["aws_execution"]["candidate_screen_authorized"] = False
    benchmark["candidate_ranking_authorized"] = False
    benchmark["pag_calibration_completed"] = False
    return status


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--status-json",
        type=Path,
        default=Path("research/status/status.json"),
    )
    parser.add_argument(
        "--audit-json",
        type=Path,
        default=Path("research/data/gluhut-charge-symmetry.json"),
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    status = json.loads(args.status_json.read_text(encoding="utf-8"))
    audit = json.loads(args.audit_json.read_text(encoding="utf-8"))
    updated = update_status(status, audit)
    args.status_json.write_text(
        json.dumps(updated, indent=2, sort_keys=False) + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "schema": updated["schema"],
                "primary_candidate": updated["primary_candidate"],
                "charge_symmetry": {
                    key: updated["public_benchmark_calibration"][key]
                    for key in (
                        "compiled_host_charge_vectors_identical_across_four_sugars",
                        "parameter_labeled_heavy_graph_automorphism_count",
                        "parameter_labeled_heavy_graph_orbit_count",
                        "maximum_graph_equivalent_heavy_charge_range_e",
                        "maximum_graph_equivalent_hydrogen_charge_range_e",
                        "compiled_charges_exactly_equal_within_preserved_graph_symmetries",
                        "large_graph_equivalent_charge_asymmetry_detected",
                        "direct_compiled_charge_transfer_to_pag_or_mcq1_authorized",
                    )
                },
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
