#!/usr/bin/env python3
"""Inventory Amber residue-template candidates for the NNC2215 hybrid graph.

This is a source-provenance and template-name inventory only. It does not
assign charges, parameterize the conjugate, build an OpenMM System, create
coordinates, or run molecular dynamics.
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import re
import sys
import xml.etree.ElementTree as ET
from collections import defaultdict
from typing import Any

import openmm
from openmm import app
from openff import interchange as openff_interchange
from openff import toolkit as openff_toolkit
from openff import units as openff_units
from rdkit import rdBase


ROOT = Path(__file__).resolve().parent
GRAPH_PATH = ROOT / "data" / "nnc2215-full-conjugate-graph.json"
SAGE_COVERAGE_PATH = ROOT / "data" / "nnc2215-forcefield-label-coverage.json"
LOCK_PATH = ROOT / "pag-chemistry-environment-macos-arm64.lock.txt"
OUTPUT_PATH = ROOT / "data" / "nnc2215-hybrid-feasibility.json"

EXPECTED_LOCK_VERSIONS = {
    "python": "3.12.15",
    "openff-toolkit": "0.19.0",
    "openff-interchange": "0.5.5",
    "openff-units": "0.4.0",
    "openff-forcefields": "2026.09.0",
    "openmm": "8.6.1",
    "rdkit": "2026.03.1",
    "ambertools": "26.0",
    "openff-nagl": "0.6.1",
    "openff-nagl-models": "2026.09.0",
}

EXPECTED_FORMULA = "C325H461N79O97S6"
EXPECTED_INSULIN_RESIDUES = 50
EXPECTED_RESIDUE_NAME_COUNT = 17
EXPECTED_LITERAL_TEMPLATE_NAME_COUNT = 16
EXPECTED_MAPPED_FRAGMENT_HEAVY_ATOMS = 109
EXPECTED_SG_SG_EDGES = 3
EXPECTED_JOINS = {
    ("B", 1, "PHE", "N", "B1", 2),
    ("B", 29, "LYS", "NZ", "B29", 84),
}
STATUS = "AMBER_TEMPLATE_CANDIDATES_IDENTIFIED; CUSTOM_CONJUGATE_PARAMETERS_REQUIRED"


def sha256_bytes(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def locked_versions_from_file(raw: bytes) -> dict[str, str]:
    filenames = []
    for line in raw.decode("utf-8").splitlines():
        line = line.strip()
        if not line.startswith("http"):
            continue
        filename = line.split("#", 1)[0].rsplit("/", 1)[-1]
        if filename.endswith(".tar.bz2"):
            filename = filename[:-8]
        elif filename.endswith(".conda"):
            filename = filename[:-6]
        filenames.append(filename)

    versions: dict[str, str] = {}
    for package, expected in EXPECTED_LOCK_VERSIONS.items():
        matches = []
        for filename in filenames:
            prefix = f"{package}-"
            if not filename.startswith(prefix):
                continue
            candidate = filename[len(prefix) :].split("-", 1)[0]
            if re.fullmatch(r"[0-9][A-Za-z0-9.]*", candidate):
                matches.append(candidate)
        if len(matches) != 1:
            raise ValueError(f"Expected one lock entry for {package}; found {matches}")
        versions[package] = matches[0]
        if matches[0] != expected:
            raise ValueError(
                f"Pinned lock version mismatch for {package}: {matches[0]} != {expected}"
            )
    return versions


def installed_locked_versions() -> dict[str, str]:
    prefix = os.environ.get("CONDA_PREFIX")
    if not prefix:
        raise ValueError("Run this audit inside the pinned micromamba environment")

    metadata_dir = Path(prefix) / "conda-meta"
    records = []
    for path in metadata_dir.glob("*.json"):
        records.append(json.loads(path.read_text(encoding="utf-8")))

    installed: dict[str, str] = {}
    for package, expected in EXPECTED_LOCK_VERSIONS.items():
        matches = [record for record in records if record.get("name") == package]
        if len(matches) != 1:
            raise ValueError(f"Expected one installed record for {package}; found {len(matches)}")
        version = str(matches[0]["version"])
        if version != expected:
            raise ValueError(f"Installed {package} version {version} != pinned {expected}")
        installed[package] = version

    module_versions = {
        "openff-toolkit": openff_toolkit.__version__,
        "openff-interchange": openff_interchange.__version__,
        "openff-units": openff_units.__version__,
        "openmm": openmm.__version__,
        "rdkit": rdBase.rdkitVersion,
    }
    for package, observed in module_versions.items():
        if observed != installed[package]:
            raise ValueError(
                f"Imported {package} module version {observed} != installed {installed[package]}"
            )
    python_version = sys.version.split()[0]
    if python_version != installed["python"]:
        raise ValueError(f"Python runtime {python_version} != pinned {installed['python']}")
    return installed


def load_amber_templates(root_xml: Path) -> tuple[set[str], list[dict[str, str]]]:
    data_dir = root_xml.parent.resolve()
    root_xml = root_xml.resolve()
    if not root_xml.is_relative_to(data_dir):
        raise ValueError(f"Amber root XML escapes OpenMM data directory: {root_xml}")

    pending = [root_xml]
    visited: dict[str, str] = {}
    template_names: set[str] = set()
    while pending:
        path = pending.pop().resolve()
        if not path.is_relative_to(data_dir):
            raise ValueError(f"Amber include escapes OpenMM data directory: {path}")
        relative = path.relative_to(data_dir).as_posix()
        if relative in visited:
            continue
        if not path.is_file():
            raise ValueError(f"Amber XML include does not exist: {relative}")

        raw = path.read_bytes()
        visited[relative] = sha256_bytes(raw)
        try:
            xml_root = ET.fromstring(raw)
        except ET.ParseError as exc:
            raise ValueError(f"Could not parse Amber XML {relative}: {exc}") from exc
        template_names.update(
            node.attrib["name"]
            for node in xml_root.iter("Residue")
            if "name" in node.attrib
        )
        for include in xml_root.iter("Include"):
            include_name = include.attrib.get("file")
            if not include_name:
                raise ValueError(f"Amber XML include in {relative} has no file attribute")
            candidate = (data_dir / include_name).resolve()
            if not candidate.is_relative_to(data_dir):
                raise ValueError(
                    f"Amber include escapes OpenMM data directory: {include_name}"
                )
            if not candidate.is_file():
                raise ValueError(f"Missing Amber XML include {include_name} from {relative}")
            pending.append(candidate)

    file_hashes = [
        {"path": path, "sha256": digest}
        for path, digest in sorted(visited.items())
    ]
    return template_names, file_hashes


def residue_key(provenance: dict[str, Any]) -> tuple[str, int, str]:
    return (
        str(provenance["chain"]),
        int(provenance["residue_number"]),
        str(provenance["residue_name"]),
    )


def residue_record(key: tuple[str, int, str]) -> dict[str, Any]:
    chain, number, name = key
    return {"chain": chain, "residue_number": number, "residue_name": name}


def classify_graph(
    graph: dict[str, Any], amber_template_names: set[str]
) -> dict[str, Any]:
    atoms = graph["atoms"]
    bonds = graph["bonds"]
    node_atoms = {atom["node_id"]: atom for atom in atoms}
    insulin_residues: dict[tuple[str, int, str], list[dict[str, Any]]] = defaultdict(list)
    mapped_atoms_by_component: dict[str, list[dict[str, Any]]] = defaultdict(list)

    for atom in atoms:
        provenance = atom.get("provenance", {})
        origin = provenance.get("origin")
        if origin == "insulin_chain_template":
            insulin_residues[residue_key(provenance)].append(atom)
        elif origin == "mapped_source_fragment":
            component = str(provenance["component"])
            mapped_atoms_by_component[component].append(atom)
        else:
            raise ValueError(f"Unclassified graph atom origin at {atom.get('node_id')}: {origin}")

    if len(insulin_residues) != EXPECTED_INSULIN_RESIDUES:
        raise ValueError(f"Expected 50 insulin residues; found {len(insulin_residues)}")
    observed_names = sorted({key[2] for key in insulin_residues})
    if len(observed_names) != EXPECTED_RESIDUE_NAME_COUNT:
        raise ValueError(f"Expected 17 residue-name labels; found {len(observed_names)}")

    mapped_count = sum(len(group) for group in mapped_atoms_by_component.values())
    if mapped_count != EXPECTED_MAPPED_FRAGMENT_HEAVY_ATOMS:
        raise ValueError(f"Expected 109 mapped-fragment heavy atoms; found {mapped_count}")

    atom_by_residue_name: dict[tuple[str, int, str], dict[str, dict[str, Any]]] = {}
    for key, group in insulin_residues.items():
        atom_names = [atom["provenance"]["atom_name"] for atom in group]
        if len(atom_names) != len(set(atom_names)):
            raise ValueError(f"Duplicate atom name in insulin residue {key}")
        atom_by_residue_name[key] = {
            atom["provenance"]["atom_name"]: atom for atom in group
        }

    literal_names = sorted(name for name in observed_names if name in amber_template_names)
    if len(literal_names) != EXPECTED_LITERAL_TEMPLATE_NAME_COUNT:
        raise ValueError(
            f"Expected 16 observed literal template names; found {len(literal_names)}: {literal_names}"
        )
    missing_literal_names = sorted(
        name for name in observed_names if name not in amber_template_names
    )
    if missing_literal_names != ["HIS"]:
        raise ValueError(f"Expected HIS as the only missing literal name; found {missing_literal_names}")
    if not {"HID", "HIE", "HIP", "CYX"} <= amber_template_names:
        absent = sorted({"HID", "HIE", "HIP", "CYX"} - amber_template_names)
        raise ValueError(f"Missing required Amber variant templates: {absent}")

    sg_sg_edges: list[tuple[tuple[str, int, str], tuple[str, int, str], str]] = []
    joins: list[dict[str, Any]] = []
    for bond in bonds:
        left_id, right_id = bond["atom_node_ids"]
        left = node_atoms[left_id]
        right = node_atoms[right_id]
        left_prov = left["provenance"]
        right_prov = right["provenance"]

        if (
            left_prov.get("origin") == "insulin_chain_template"
            and right_prov.get("origin") == "insulin_chain_template"
            and left["element"] == "S"
            and right["element"] == "S"
            and left_prov.get("atom_name") == "SG"
            and right_prov.get("atom_name") == "SG"
        ):
            sg_sg_edges.append((residue_key(left_prov), residue_key(right_prov), bond["edge_id"]))

        if {left_prov.get("origin"), right_prov.get("origin")} == {
            "insulin_chain_template",
            "mapped_source_fragment",
        }:
            protein_atom, mapped_atom = (
                (left, right)
                if left_prov["origin"] == "insulin_chain_template"
                else (right, left)
            )
            protein_prov = protein_atom["provenance"]
            mapped_prov = mapped_atom["provenance"]
            joins.append({
                "edge_id": bond["edge_id"],
                "protein": {
                    "chain": protein_prov["chain"],
                    "residue_number": protein_prov["residue_number"],
                    "residue_name": protein_prov["residue_name"],
                    "atom_name": protein_prov["atom_name"],
                    "node_id": protein_atom["node_id"],
                },
                "conjugate": {
                    "component": mapped_prov["component"],
                    "source_map_id": mapped_prov["source_map_id"],
                    "map_role": mapped_prov.get("map_role", ""),
                    "element": mapped_atom["element"],
                    "node_id": mapped_atom["node_id"],
                },
                "bond_type": bond["bond_type"],
            })

    observed_joins = {
        (
            join["protein"]["chain"],
            int(join["protein"]["residue_number"]),
            join["protein"]["residue_name"],
            join["protein"]["atom_name"],
            join["conjugate"]["component"],
            int(join["conjugate"]["source_map_id"]),
        )
        for join in joins
    }
    if observed_joins != EXPECTED_JOINS or len(joins) != 2:
        raise ValueError(f"Unexpected insulin-to-conjugate joins: {sorted(observed_joins)}")

    if len(sg_sg_edges) != EXPECTED_SG_SG_EDGES:
        raise ValueError(f"Expected exactly three graph SG-SG disulfide edges; found {len(sg_sg_edges)}")
    cys_residues = sorted(key for key in insulin_residues if key[2] == "CYS")
    cysteine_partners: dict[tuple[str, int, str], tuple[str, int, str]] = {}
    for left, right, _edge_id in sg_sg_edges:
        if left in cysteine_partners or right in cysteine_partners:
            raise ValueError("A cysteine participates in multiple graph SG-SG edges")
        cysteine_partners[left] = right
        cysteine_partners[right] = left
    if len(cys_residues) != 6 or set(cysteine_partners) != set(cys_residues):
        raise ValueError("Expected all six CYS residues to participate in three disulfide edges")

    histidine_residues = sorted(key for key in insulin_residues if key[2] == "HIS")
    if len(histidine_residues) != 2:
        raise ValueError(f"Expected two neutral-template HIS residues; found {len(histidine_residues)}")
    his_ne2_hydrogen_counts: dict[tuple[str, int, str], int] = {}
    for key in histidine_residues:
        ne2 = atom_by_residue_name[key].get("NE2")
        if ne2 is None:
            raise ValueError(f"HIS residue {key} has no NE2 atom in the graph ledger")
        his_ne2_hydrogen_counts[key] = int(ne2["hydrogen_count_not_as_nodes"])
    if set(his_ne2_hydrogen_counts.values()) != {1}:
        raise ValueError(f"Expected one NE2-bound hydrogen for each HIS; found {his_ne2_hydrogen_counts}")

    modified_keys = {("B", 1, "PHE"), ("B", 29, "LYS")}
    if not modified_keys <= set(insulin_residues):
        raise ValueError("Expected modified PHEB1 and LYSB29 residues in the insulin graph")

    sorted_joins = sorted(
        joins,
        key=lambda join: (
            join["protein"]["chain"],
            int(join["protein"]["residue_number"]),
        ),
    )
    sorted_disulfides = sorted(
        [
            {
                "edge_id": edge_id,
                "residue_1": residue_record(left),
                "residue_2": residue_record(right),
                "amber_variant_candidate": "CYX",
            }
            for left, right, edge_id in sg_sg_edges
        ],
        key=lambda edge: (edge["residue_1"]["chain"], edge["residue_1"]["residue_number"]),
    )

    residue_rows = []
    for key in sorted(insulin_residues):
        chain, number, name = key
        is_his = name == "HIS"
        is_disulfide_cys = key in cysteine_partners
        is_modified = key in modified_keys
        variant_candidates: list[str] = []
        evidence: dict[str, Any] = {}
        if is_his:
            variant_candidates = ["HID", "HIE", "HIP"]
            evidence = {
                "NE2_hydrogen_count_in_neutral_graph": his_ne2_hydrogen_counts[key],
                "candidate_supported_by_hydrogen_ledger": "HIE",
                "interpretation": "HIE is a candidate; the ledger does not establish an exact Amber residue-template match.",
            }
        elif is_disulfide_cys:
            variant_candidates = ["CYX"]
            partner = cysteine_partners[key]
            evidence = {
                "graph_disulfide_partner": residue_record(partner),
                "interpretation": "CYX is a disulfide-cysteine candidate; exact residue-template matching remains unresolved.",
            }

        row: dict[str, Any] = {
            **residue_record(key),
            "graph_atom_count": len(insulin_residues[key]),
            "literal_residue_name_template_available": name in amber_template_names,
            "literal_template_name": name if name in amber_template_names else None,
            "variant_template_candidates_available": variant_candidates,
            "variant_evidence": evidence or None,
            "modified_site": is_modified,
            "custom_conjugate_parameters_required": is_modified,
            "ordinary_residue_template_covers_attached_conjugate": False if is_modified else None,
            "exact_template_match_established": False,
        }
        residue_rows.append(row)

    component_rows = []
    for component in sorted(mapped_atoms_by_component):
        group = mapped_atoms_by_component[component]
        atom_records = sorted(
            [
                {
                    "source_map_id": int(atom["provenance"]["source_map_id"]),
                    "node_id": atom["node_id"],
                    "element": atom["element"],
                    "map_role": atom["provenance"].get("map_role", ""),
                }
                for atom in group
            ],
            key=lambda atom: (atom["source_map_id"], atom["node_id"]),
        )
        component_rows.append({
            "component": component,
            "heavy_atom_count": len(group),
            "mapped_atoms": atom_records,
        })

    return {
        "insulin_residues": residue_rows,
        "insulin_residue_name_labels": observed_names,
        "literal_template_residue_name_labels": literal_names,
        "literal_template_residue_name_label_count": len(literal_names),
        "residues_with_literal_template_names": sum(
            name in amber_template_names for _, _, name in insulin_residues
        ),
        "histidine_variant_candidates": {
            "observed_residues": [residue_record(key) for key in histidine_residues],
            "available_templates": ["HID", "HIE", "HIP"],
            "NE2_hydrogen_counts": [
                {
                    **residue_record(key),
                    "hydrogen_count": his_ne2_hydrogen_counts[key],
                }
                for key in histidine_residues
            ],
            "candidate_supported_by_graph_hydrogen_ledger": "HIE",
            "exact_template_match_established": False,
        },
        "disulfide_cysteines": {
            "residues": [
                {
                    **residue_record(key),
                    "amber_variant_candidate": "CYX",
                    "partner": residue_record(cysteine_partners[key]),
                }
                for key in cys_residues
            ],
            "sg_sg_edge_count": len(sg_sg_edges),
            "edges": sorted_disulfides,
        },
        "modified_sites": [
            {
                **residue_record(key),
                "literal_template_name_available": key[2] in amber_template_names,
                "custom_conjugate_parameters_required": True,
                "ordinary_residue_template_covers_attached_conjugate": False,
            }
            for key in sorted(modified_keys)
        ],
        "insulin_to_conjugate_joins": sorted_joins,
        "mapped_source_fragments": {
            "heavy_atom_count": mapped_count,
            "components": component_rows,
        },
    }


def main() -> int:
    try:
        graph_raw = GRAPH_PATH.read_bytes()
        sage_raw = SAGE_COVERAGE_PATH.read_bytes()
        lock_raw = LOCK_PATH.read_bytes()
        graph = json.loads(graph_raw)
        sage_coverage = json.loads(sage_raw)
        lock_versions = locked_versions_from_file(lock_raw)
        installed_versions = installed_locked_versions()
        if lock_versions != installed_versions:
            raise ValueError("Installed packages do not match the pinned environment lock")

        graph_digest = sha256_bytes(graph_raw)
        sage_graph_input = sage_coverage["inputs"]["full_conjugate_graph"]
        if sage_graph_input["sha256"] != graph_digest:
            raise ValueError("Sage coverage JSON does not refer to the current full-conjugate graph")
        sage_graph_summary = sage_coverage["full_neutral_constitutional_graph"]
        if (
            int(sage_graph_summary["heavy_atom_count"]) != 507
            or int(sage_graph_summary["atom_count_after_hydrogen_materialization"]) != 968
        ):
            raise ValueError("Sage coverage dimensions disagree with the expected NNC2215 graph")

        openmm_data_dir = Path(app.__file__).resolve().parent / "data"
        amber_root = openmm_data_dir / "amber14-all.xml"
        if not amber_root.is_file():
            raise ValueError(f"Installed OpenMM amber14-all.xml not found at {amber_root}")
        amber_names, amber_xml_hashes = load_amber_templates(amber_root)
        classifications = classify_graph(graph, amber_names)

        graph_summary = graph["graph"]
        formula = graph["formula"]["rdkit_formula"]
        if formula != EXPECTED_FORMULA:
            raise ValueError(f"Unexpected full graph formula: {formula}")
        if int(graph_summary["heavy_atom_count"]) != 507 or int(graph_summary["edge_count"]) != 527:
            raise ValueError("Unexpected source graph atom or bond count")

        output = {
            "schema": "dia-health-one.nnc2215-hybrid-feasibility.v1",
            "scope": "Amber template-name and source-provenance feasibility inventory only",
            "inputs": {
                "full_conjugate_graph": {
                    "path": str(GRAPH_PATH.relative_to(ROOT)),
                    "sha256": graph_digest,
                },
                "sage_label_coverage": {
                    "path": str(SAGE_COVERAGE_PATH.relative_to(ROOT)),
                    "sha256": sha256_bytes(sage_raw),
                },
                "pinned_environment_lock": {
                    "path": str(LOCK_PATH.relative_to(ROOT)),
                    "sha256": sha256_bytes(lock_raw),
                },
            },
            "software": {
                "python_version": sys.version.split()[0],
                "locked_package_versions": installed_versions,
                "sage_force_field": sage_coverage["software"]["sage_force_field"],
                "nagl_model": "openff-gnn-am1bcc-1.0.0.pt",
                "openmm_amber_xml_root": "amber14-all.xml",
            },
            "amber_xml_inventory": {
                "data_directory_basis": "openmm.app.__file__ parent / data",
                "included_xml_files": amber_xml_hashes,
                "residue_template_names": sorted(amber_names),
                "residue_template_name_count": len(amber_names),
            },
            "graph": {
                "representation": graph["representation"]["type"],
                "formula": formula,
                "heavy_atom_count": int(graph_summary["heavy_atom_count"]),
                "bond_count": int(graph_summary["edge_count"]),
                "atom_count_after_hydrogen_expansion": int(
                    sage_graph_summary["atom_count_after_hydrogen_materialization"]
                ),
                "formal_charge": int(sage_graph_summary["formal_charge"]),
                "coordinates_present": False,
                "sage_label_coverage_context": {
                    "scope": sage_coverage["scope"],
                    "label_coverage": sage_graph_summary["label_coverage"],
                    "electrostatics_label_count": int(
                        sage_graph_summary["sage_label_counts"].get("Electrostatics", 0)
                    ),
                    "partial_charges_assigned_by_this_audit": False,
                },
            },
            "insulin_template_coverage": {
                "residue_count": len(classifications["insulin_residues"]),
                "observed_residue_name_label_count": len(
                    classifications["insulin_residue_name_labels"]
                ),
                "insulin_residue_name_labels": classifications["insulin_residue_name_labels"],
                "literal_template_residue_name_label_count": classifications[
                    "literal_template_residue_name_label_count"
                ],
                "literal_template_residue_name_labels": classifications[
                    "literal_template_residue_name_labels"
                ],
                "residues_with_literal_template_names": classifications[
                    "residues_with_literal_template_names"
                ],
                "residues": classifications["insulin_residues"],
                "histidine_variant_candidates": classifications["histidine_variant_candidates"],
                "disulfide_cysteines": classifications["disulfide_cysteines"],
                "modified_sites": classifications["modified_sites"],
                "exact_template_matching": "unresolved; name and variant candidate availability is not an exact template match",
                "termini_and_charge_compatibility": "unresolved for the source-derived neutral graph and intended physiological charge states",
            },
            "conjugate_scope": {
                "mapped_fragment_heavy_atom_count": classifications[
                    "mapped_source_fragments"
                ]["heavy_atom_count"],
                "mapped_source_fragments_by_component": classifications[
                    "mapped_source_fragments"
                ]["components"],
                "insulin_to_conjugate_joins": classifications[
                    "insulin_to_conjugate_joins"
                ],
                "custom_conjugate_parameters_required": True,
                "ordinary_phe_or_lys_templates_parameterize_the_attached_conjugate": False,
            },
            "nagl_model_size_domain": {
                "source": "https://docs.openforcefield.org/projects/nagl-models/en/stable/models/openff-gnn-am1bcc-1.0.0/",
                "published_test_heavy_atom_distribution": {
                    "80_to_89_heavy_atoms": 2,
                    "90_or_more_heavy_atoms": 0,
                    "highest_test_heavy_atom_bin": "80-89",
                },
                "published_test_total_atom_count_maximum": 169,
                "published_peptide_training_context": {
                    "training_peptide_fragments": "1-4 amino-acid residues",
                    "training_peptide_fragment_count": 22690,
                    "interpretation": "This describes the peptide training subset; the model also uses non-peptide training datasets.",
                },
                "current_graph_heavy_atoms": 507,
                "current_graph_atoms_after_hydrogen_expansion": 968,
                "whole_graph_nagl_use": "extrapolation beyond the published test heavy-atom size coverage",
            },
            "custom_parameterization_precedents": {
                "ambertools_custom_polymer_residue": {
                    "source": "https://ambermd.org/tutorials/basic/tutorial5/index.php",
                    "workflow": ["Antechamber", "prepgen", "parmchk2", "LEaP"],
                    "classification": "official modified-polymer-residue workflow precedent; not NNC2215 parameterization evidence",
                    "generated_charges_and_types_require_review": True,
                },
                "openff_sage_ff14sb_interchange_example": {
                    "source": "https://docs.openforcefield.org/en/latest/examples/openforcefield/openff-interchange/protein_ligand/protein_ligand.html",
                    "classification": "experimental combination example for separate protein and ligand components",
                    "demonstrates_covalent_conjugate_parameterization": False,
                },
                "isolated_pag_gaff2_result": {
                    "source": "data/pag-gaff2-validation.json",
                    "classification": "ligand-only toolchain precedent; does not parameterize the NNC2215 conjugate or hybrid",
                },
            },
            "not_established": [
                "exact Amber residue-template matching or terminal-template compatibility",
                "physiological protonation microstate or charge-state compatibility",
                "partial charges or validated cross-boundary force-field terms",
                "custom parameters for PHEB1, LYSB29, or the mapped conjugate fragments",
                "unique B29 coordinate-based atom mapping; the source graph preserves its documented mapping ambiguity",
                "complete stereochemistry or coordinate-based geometry validation",
                "an OpenMM System, simulation-ready topology, or molecular dynamics result",
                "binding pose, affinity, pharmacological efficacy, or clinical inference",
            ],
            "assessment": {
                "status": STATUS,
                "insulin_residues_classified": len(classifications["insulin_residues"]),
                "insulin_residue_name_labels": len(classifications["insulin_residue_name_labels"]),
                "literal_template_name_labels": classifications[
                    "literal_template_residue_name_label_count"
                ],
                "mapped_fragment_heavy_atoms": classifications[
                    "mapped_source_fragments"
                ]["heavy_atom_count"],
                "insulin_to_conjugate_join_count": len(
                    classifications["insulin_to_conjugate_joins"]
                ),
                "validated_conjugate_force_field_count": 0,
                "partial_charges_generated": False,
                "force_field_parameters_generated": False,
                "openmm_system_generated": False,
                "coordinates_generated_or_used": False,
                "molecular_dynamics_run": False,
            },
        }

        OUTPUT_PATH.write_text(
            json.dumps(output, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        print(json.dumps({
            "audit": str(OUTPUT_PATH),
            "status": STATUS,
            "insulin_residue_count": len(classifications["insulin_residues"]),
            "residue_name_label_count": len(classifications["insulin_residue_name_labels"]),
            "literal_template_name_count": classifications[
                "literal_template_residue_name_label_count"
            ],
            "mapped_fragment_heavy_atom_count": classifications[
                "mapped_source_fragments"
            ]["heavy_atom_count"],
            "join_count": len(classifications["insulin_to_conjugate_joins"]),
            "validated_conjugate_force_field_count": 0,
        }, indent=2, sort_keys=True))
        return 0
    except (OSError, ValueError, KeyError, json.JSONDecodeError, ET.ParseError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
