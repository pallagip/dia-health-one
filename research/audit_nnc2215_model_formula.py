#!/usr/bin/env python3
"""Compare publisher-model atom inventories with the published NNC2215 formula.

This audit counts coordinate-record elements only. It does not reconstruct a
chemical graph, infer bond orders, or assign protonation/charge states.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
MANIFEST_PATH = ROOT / "data" / "structure-manifest.json"
OUTPUT_PATH = ROOT / "data" / "nnc2215-model-formula-audit.json"
PUBLISHED_FORMULA = {"C": 325, "H": 461, "N": 79, "O": 97, "S": 6}
INSULIN_SEQUENCES = {
    "A": "GIVEQCCTSICSLYQLENYCN",
    "B": "FVNQHLCGSHLVEALYLVCGERGFFYTPK",
}
AMINO_ACID_FORMULAE = {
    "A": {"C": 3, "H": 7, "N": 1, "O": 2},
    "R": {"C": 6, "H": 14, "N": 4, "O": 2},
    "N": {"C": 4, "H": 8, "N": 2, "O": 3},
    "D": {"C": 4, "H": 7, "N": 1, "O": 4},
    "C": {"C": 3, "H": 7, "N": 1, "O": 2, "S": 1},
    "Q": {"C": 5, "H": 10, "N": 2, "O": 3},
    "E": {"C": 5, "H": 9, "N": 1, "O": 4},
    "G": {"C": 2, "H": 5, "N": 1, "O": 2},
    "H": {"C": 6, "H": 9, "N": 3, "O": 2},
    "I": {"C": 6, "H": 13, "N": 1, "O": 2},
    "L": {"C": 6, "H": 13, "N": 1, "O": 2},
    "K": {"C": 6, "H": 14, "N": 2, "O": 2},
    "M": {"C": 5, "H": 11, "N": 1, "O": 2, "S": 1},
    "F": {"C": 9, "H": 11, "N": 1, "O": 2},
    "P": {"C": 5, "H": 9, "N": 1, "O": 2},
    "S": {"C": 3, "H": 7, "N": 1, "O": 3},
    "T": {"C": 4, "H": 9, "N": 1, "O": 3},
    "W": {"C": 11, "H": 12, "N": 2, "O": 2},
    "Y": {"C": 9, "H": 11, "N": 1, "O": 3},
    "V": {"C": 5, "H": 11, "N": 1, "O": 2},
}
AMINO_ACID_NAMES = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
SOURCE_URL = (
    "https://media.springernature.com/original/springer-static/esm/"
    "art%3A10.1038%2Fs41586-024-08042-3/MediaObjects/"
    "41586_2024_8042_MOESM1_ESM.pdf"
)
PATENT_URL = "https://patents.google.com/patent/WO2023144240A1/en"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def formula_text(counts: dict[str, int]) -> str:
    order = ("C", "H", "N", "O", "S", "P", "F", "CL", "BR", "I")
    elements = [element for element in order if counts.get(element, 0)]
    elements.extend(sorted(
        element for element in set(counts) - set(elements)
        if counts.get(element, 0)
    ))
    return "".join(
        element + (str(counts[element]) if counts[element] != 1 else "")
        for element in elements
    )


def difference(observed: Counter[str], expected: dict[str, int]) -> dict[str, int]:
    elements = sorted(set(observed) | set(expected))
    return {
        element: observed[element] - expected.get(element, 0)
        for element in elements
        if observed[element] - expected.get(element, 0) != 0
    }


def derive_desb30_formula(disulfide_bond_count: int) -> Counter[str]:
    """Derive the two-chain desB30 insulin formula from residue formulas."""
    result: Counter[str] = Counter()
    for sequence in INSULIN_SEQUENCES.values():
        chain_formula: Counter[str] = Counter()
        for residue in sequence:
            chain_formula.update(AMINO_ACID_FORMULAE[residue])
        # Condense each separate peptide chain, preserving its terminal groups.
        chain_formula.subtract(Counter(H=2 * (len(sequence) - 1), O=len(sequence) - 1))
        result.update(chain_formula)
    result.subtract(Counter(H=2 * disulfide_bond_count))
    return result


def source_component_decomposition() -> dict[str, object]:
    """Derive the two attachment formulas from documented NNC2215 building blocks."""
    b1_peracetyl_acid = Counter(C=18, H=26, O=13)
    b1_deacetylated_acid = Counter(b1_peracetyl_acid)
    b1_deacetylated_acid.subtract(Counter(C=8, H=8, O=4))
    b1_attachment = Counter(b1_deacetylated_acid)
    b1_attachment.subtract(Counter(H=2, O=1))  # amide condensation removes H2O

    b29_active_ester = Counter(C=10, H=11, N=1, O=5)
    nhs_leaving_group = Counter(C=4, H=5, N=1, O=3)
    b29_linker_net = Counter(b29_active_ester)
    b29_linker_net.subtract(nhs_leaving_group)
    macrocycle_roof_propyl_azide = Counter(C=56, H=63, N=15, O=12)
    b29_attachment = Counter(b29_linker_net)
    b29_attachment.update(macrocycle_roof_propyl_azide)

    total_attachment = Counter(b1_attachment)
    total_attachment.update(b29_attachment)
    expected_from_published_formula = Counter(PUBLISHED_FORMULA)
    expected_from_published_formula.subtract(
        derive_desb30_formula(disulfide_bond_count=3)
    )
    # Counter subtraction retains zero-valued keys (here, S: 0), while the
    # independently assembled attachment formula has no such keys. Normalize
    # before comparing formula dictionaries.
    expected_from_published_formula += Counter()
    return {
        "sources": [
            {
                "role": "NNC2215 synthesis and component formulae",
                "title": "Glucose-sensitive insulin with attenuation of hypoglycaemia — Supplementary Information",
                "url": SOURCE_URL,
                "locations": [
                    "Methods §§1–4, PDF pp. 2–11",
                    "B1 peracetyl glycoside acid formula C18H26O13, p. 3",
                    "B29 active ester formula C10H11NO5, p. 2",
                    "Macrocycle roof propyl azide formula C56H63N15O12, p. 10",
                    "NNC2215 formula C325H461N79O97S6, p. 11",
                ],
            },
            {
                "role": "Independent synthesis-description cross-check",
                "publication": "WO2023144240A1, Example 2",
                "url": PATENT_URL,
                "reported_product": "B1-glucoside-EG-Ac B29-MC-propyl-triazole-propanoxycarbonyl desB30 human insulin (INS2)",
                "reported_calculated_ions": {
                    "m_plus_4h_mz": 1805.8,
                    "m_plus_5h_mz": 1444.8,
                },
                "interpretation": "The patent describes the same named B1/B29 building blocks and reports the same calculated mass-to-charge values as the Nature SI. This supports product-identity cross-checking; it does not supply an atom-mapped topology for the PDB residues.",
            },
        ],
        "derivation": {
            "b1_peracetyl_acid_formula": formula_text(dict(b1_peracetyl_acid)),
            "b1_deacetylation_count": 4,
            "formula_removed_per_o_deacetylation": "C2H2O",
            "b1_deacetylated_acid_formula": formula_text(dict(b1_deacetylated_acid)),
            "b1_net_attachment_formula_after_condensation": formula_text(dict(b1_attachment)),
            "b1_net_attachment_element_counts": dict(sorted(b1_attachment.items())),
            "b29_active_ester_formula": formula_text(dict(b29_active_ester)),
            "n_hydroxysuccinimide_leaving_group_formula": formula_text(dict(nhs_leaving_group)),
            "b29_linker_formula_after_nhs_departure": formula_text(dict(b29_linker_net)),
            "macrocycle_roof_propyl_azide_formula": formula_text(dict(macrocycle_roof_propyl_azide)),
            "b29_net_attachment_formula_after_click": formula_text(dict(b29_attachment)),
            "b29_net_attachment_element_counts": dict(sorted(b29_attachment.items())),
            "expected_total_attachment_formula": formula_text(dict(total_attachment)),
            "expected_total_attachment_element_counts": dict(sorted(total_attachment.items())),
            "published_formula_minus_derived_parent_formula": formula_text(
                dict(expected_from_published_formula)
            ),
            "published_formula_minus_derived_parent_element_counts": dict(
                sorted(expected_from_published_formula.items())
            ),
            "source_decomposition_matches_published_formula_difference": (
                dict(total_attachment) == dict(expected_from_published_formula)
            ),
            "assumptions": [
                "Four O-acetyl groups are removed during the reported saponification.",
                "B1 amide formation removes H2O from the deacetylated acid and insulin amino group.",
                "The B29 NHS active ester transfers the C6H6O2 linker composition after NHS leaves.",
                "The azide–alkyne cycloaddition incorporates both source molecules without atom loss.",
            ],
        },
    }


def sequence_and_terminus_inventory(path: Path) -> dict[str, object]:
    residues: dict[str, dict[tuple[str, str], str]] = defaultdict(dict)
    oxt_atoms: list[dict[str, str]] = []
    disulfide_records: list[str] = []
    with path.open(encoding="ascii") as stream:
        for line in stream:
            if line[:6].strip() == "SSBOND":
                disulfide_records.append(line.rstrip())
            if line[:6].strip() != "ATOM":
                continue
            chain = line[21].strip()
            residue_name = line[17:20].strip()
            residue_number = line[22:26].strip()
            atom_name = line[12:16].strip()
            if chain in INSULIN_SEQUENCES:
                residues[chain][(int(residue_number), line[26].strip())] = residue_name
            if atom_name == "OXT":
                oxt_atoms.append({
                    "chain": chain,
                    "residue": residue_name,
                    "residue_number": residue_number,
                    "element": line[76:78].strip().upper(),
                })
    sequences = {
        chain: "".join(AMINO_ACID_NAMES[name] for _, name in sorted(items.items()))
        for chain, items in residues.items()
    }
    if sequences != INSULIN_SEQUENCES:
        raise ValueError(f"Insulin chain sequence mismatch in {path}: {sequences}")
    return {
        "chain_sequences": sequences,
        "disulfide_bond_count": len(disulfide_records),
        "disulfide_records": disulfide_records,
        "chain_a_asn21_oxt_present": any(
            item["chain"] == "A" and item["residue"] == "ASN"
            and item["residue_number"].strip() == "21" and item["element"] == "O"
            for item in oxt_atoms
        ),
        "oxt_atoms": oxt_atoms,
    }


def audit_model(model: dict[str, object]) -> dict[str, object]:
    path = ROOT / str(model["path"])
    observed_sha = sha256(path)
    if observed_sha != model["sha256"]:
        raise ValueError(f"Source PDB checksum mismatch for {path}")

    all_atoms: Counter[str] = Counter()
    by_record: dict[str, Counter[str]] = defaultdict(Counter)
    by_residue: dict[str, Counter[str]] = defaultdict(Counter)
    excluded_free_glucose: Counter[str] = Counter()
    atom_records = 0
    altloc_records = 0

    with path.open(encoding="ascii") as stream:
        for line in stream:
            record = line[:6].strip()
            if record not in {"ATOM", "HETATM"}:
                continue
            atom_records += 1
            altloc = line[16].strip()
            if altloc:
                altloc_records += 1
            element = line[76:78].strip().upper()
            residue = line[17:20].strip()
            if not element:
                raise ValueError(f"Missing element field in {path}: {line.rstrip()}")
            if residue == "BGL":
                excluded_free_glucose[element] += 1
                continue
            all_atoms[element] += 1
            by_record[record][element] += 1
            by_residue[residue][element] += 1

    if altloc_records:
        raise ValueError(f"Alternate-location atom records need explicit handling: {path}")

    expected_atom_records = int(model["coordinate_atom_records"])
    if atom_records != expected_atom_records:
        raise ValueError(
            f"Coordinate-record count changed for {path}: "
            f"{atom_records} != manifest {expected_atom_records}"
        )

    observed = dict(sorted(all_atoms.items()))
    expected = dict(PUBLISHED_FORMULA)
    heavy_observed = {element: count for element, count in observed.items() if element != "H"}
    heavy_expected = {element: count for element, count in expected.items() if element != "H"}

    protein_sequence = sequence_and_terminus_inventory(path)
    protein_heavy = {
        element: count for element, count in by_record["ATOM"].items() if element != "H"
    }
    hetero_heavy = {
        element: count for element, count in by_record["HETATM"].items() if element != "H"
    }
    return {
        "state": model["state"],
        "file": model["path"],
        "sha256": observed_sha,
        "coordinate_atom_record_count_including_free_glucose": atom_records,
        "excluded_residue": {
            "label": "BGL",
            "reason": "Free D-glucose is modelled in the open state but is not part of the NNC2215 molecule.",
            "element_counts": dict(sorted(excluded_free_glucose.items())),
        },
        "model_molecule_atom_inventory": {
            "formula_from_records": formula_text(observed),
            "element_counts": observed,
            "heavy_atom_formula_from_records": formula_text(heavy_observed),
            "heavy_element_counts": heavy_observed,
            "by_record_type": {
                record: dict(sorted(counts.items()))
                for record, counts in sorted(by_record.items())
            },
            "by_residue": {
                residue: dict(sorted(counts.items()))
                for residue, counts in sorted(by_residue.items())
            },
            "protein_atom_records": {
                "element_counts": dict(sorted(by_record["ATOM"].items())),
                "heavy_element_counts": protein_heavy,
                **protein_sequence,
            },
            "covalent_nonstandard_hetero_atom_records_excluding_free_glucose": {
                "element_counts": dict(sorted(by_record["HETATM"].items())),
                "heavy_element_counts": hetero_heavy,
                "by_residue": {
                    residue: {
                        element: count for element, count in sorted(counts.items())
                        if element != "H"
                    }
                    for residue, counts in sorted(by_residue.items())
                    if residue in {"PHX", "LYX", "MC"}
                },
            },
        },
        "difference_from_published_formula": {
            "element_count_delta_model_minus_source": difference(all_atoms, expected),
            "heavy_element_count_delta_model_minus_source": difference(
                Counter(heavy_observed), heavy_expected
            ),
        },
    }


def main() -> None:
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    models = [audit_model(model) for model in manifest["files"]]
    disulfide_counts = {
        item["model_molecule_atom_inventory"]["protein_atom_records"]["disulfide_bond_count"]
        for item in models
    }
    if disulfide_counts != {3}:
        raise ValueError(f"Expected three disulfide records in both models, got {disulfide_counts}")
    disulfide_bond_count = disulfide_counts.pop()
    parent_formula = derive_desb30_formula(disulfide_bond_count)
    expected_attachment_formula = Counter(PUBLISHED_FORMULA)
    expected_attachment_formula.subtract(parent_formula)
    source_components = source_component_decomposition()
    source_derivation = source_components["derivation"]
    if not source_derivation["source_decomposition_matches_published_formula_difference"]:
        raise ValueError("Source building-block formulae do not sum to the published formula difference")
    expected_attachment_heavy = {
        element: count for element, count in expected_attachment_formula.items()
        if element != "H" and count
    }
    expected_b1_heavy = {
        element: count for element, count in source_derivation[
            "b1_net_attachment_element_counts"
        ].items() if element != "H" and count
    }
    expected_b29_heavy = {
        element: count for element, count in source_derivation[
            "b29_net_attachment_element_counts"
        ].items() if element != "H" and count
    }
    observed_attachment_heavy = [
        item["model_molecule_atom_inventory"][
            "covalent_nonstandard_hetero_atom_records_excluding_free_glucose"
        ]["heavy_element_counts"]
        for item in models
    ]
    attachment_deltas = [
        difference(Counter(observed), expected_attachment_heavy)
        for observed in observed_attachment_heavy
    ]
    carbon_deltas = [
        item["difference_from_published_formula"]["heavy_element_count_delta_model_minus_source"].get("C", 0)
        for item in models
    ]
    component_count_comparison = []
    for item in models:
        by_residue = item["model_molecule_atom_inventory"][
            "covalent_nonstandard_hetero_atom_records_excluding_free_glucose"
        ]["by_residue"]
        b1_observed = Counter({
            element: count for element, count in by_residue.get("PHX", {}).items()
            if element != "H"
        })
        b29_observed: Counter[str] = Counter()
        for residue in ("LYX", "MC"):
            b29_observed.update({
                element: count for element, count in by_residue.get(residue, {}).items()
                if element != "H"
            })
        component_count_comparison.append({
            "state": item["state"],
            "b1_source_component_vs_pdb_phx": {
                "source_expected_heavy_element_counts": expected_b1_heavy,
                "pdb_heavy_element_counts": dict(sorted(b1_observed.items())),
                "pdb_minus_source_delta": difference(b1_observed, expected_b1_heavy),
            },
            "b29_source_component_vs_pdb_lyx_mc_pool": {
                "source_expected_heavy_element_counts": expected_b29_heavy,
                "pdb_heavy_element_counts": dict(sorted(b29_observed.items())),
                "pdb_minus_source_delta": difference(b29_observed, expected_b29_heavy),
                "scope_note": "Residue labels group the B29-linked pool for element-count comparison only; this is not an atom-by-atom graph mapping.",
            },
        })
    report = {
        "schema": "dia-health-one.nnc2215-model-formula-audit.v1",
        "run": {
            "run_id": "DHO-NNC2215-FORMULA-" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ"),
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
            "script": "research/audit_nnc2215_model_formula.py",
            "script_sha256": sha256(Path(__file__).resolve()),
            "nvidia_gpu_used": False,
        },
        "published_reference": {
            "paper_doi": "10.1038/s41586-024-08042-3",
            "supporting_information_url": SOURCE_URL,
            "location": "Methods §4, Preparation of insulin conjugate NNC2215, PDF pp. 10–11",
            "reported_hrms_formula": "C325H461N79O97S6",
            "reported_formula_element_counts": PUBLISHED_FORMULA,
        },
        "source_component_decomposition": source_components,
        "source_structure_manifest": "research/data/structure-manifest.json",
        "derived_parent_reference": {
            "label": "desB30 human insulin before the two conjugations",
            "chain_sequences": INSULIN_SEQUENCES,
            "formula": formula_text(dict(parent_formula)),
            "formula_element_counts": dict(sorted(parent_formula.items())),
            "disulfide_bond_count_from_each_pdb": disulfide_bond_count,
            "method": "Sum standard neutral amino-acid formulas for the 21-residue A chain and 29-residue desB30 B chain; subtract H2O for each peptide bond within each chain; subtract two H atoms per SSBOND record. This is a sequence-derived reference, not a separately measured value in the cited SI.",
        },
        "attachment_group_heavy_atom_comparison": {
            "expected_from_published_formula_minus_derived_parent": expected_attachment_heavy,
            "model_non_bgl_hetero_heavy_element_counts_by_state": {
                item["state"]: observed_attachment_heavy[index]
                for index, item in enumerate(models)
            },
            "model_minus_expected_heavy_element_deltas_by_state": {
                item["state"]: attachment_deltas[index]
                for index, item in enumerate(models)
            },
            "interpretation": "The source-reported B1 and B29 building-block formulae sum to the published attachment difference C72H85N15O22. The PDB PHX heavy-element inventory matches the source-derived B1 attachment C10O8. The B29-linked PDB LYX/MC pool contains C64N15O14 versus source-derived C62N15O14: nitrogen and oxygen match, with two excess carbon coordinate records. This narrows the inventory discrepancy to the B29-linked component pool, but does not establish which atom records are unsupported or missing from the source atom mapping.",
        },
        "model_source_component_count_comparison": component_count_comparison,
        "method": {
            "description": "Verify source-file SHA-256 values and coordinate-record counts; count ATOM/HETATM element fields; derive the desB30 parent formula from validated chain sequences and three disulfides; calculate attachment formulae from source-reported building-block compositions and reaction stoichiometry; and exclude the explicitly modelled free-glucose BGL residue before comparing model component inventories.",
            "chemistry_inference": False,
            "source_formula_arithmetic_performed": True,
            "bond_order_or_charge_assignment": False,
            "hydrogen_counts_interpreted_as_formula_validation": False,
            "limitation": "This combines source-reported building-block formula arithmetic with PDB coordinate-record counts; it is not a validated molecular graph or atom map. The count comparison narrows the two-carbon residual to the B29-linked LYX/MC component pool but does not identify the specific atom records, establish which model atoms are chemically intended, or validate bond orders, stereochemistry, protonation, charges, or parameters. Hydrogen coordinate counts are not used as formula-validation evidence.",
        },
        "model_results": models,
        "assessment": {
            "same_extra_heavy_carbon_count_in_both_states": len(set(carbon_deltas)) == 1 and carbon_deltas[0] == 2,
            "unresolved_coordinate_inventory_mismatch": any(carbon_deltas),
            "source_building_block_decomposition_matches_published_formula": source_derivation["source_decomposition_matches_published_formula_difference"],
            "b1_heavy_inventory_matches_source_count_in_both_states": all(
                not entry["b1_source_component_vs_pdb_phx"]["pdb_minus_source_delta"]
                for entry in component_count_comparison
            ),
            "b29_pool_has_two_excess_carbon_records_in_both_states": all(
                entry["b29_source_component_vs_pdb_lyx_mc_pool"]["pdb_minus_source_delta"] == {"C": 2}
                for entry in component_count_comparison
            ),
            "interpretation": "The source-reported B1/B29 building-block formula arithmetic reproduces the published NNC2215 attachment difference. The PDB PHX B1 component matches its expected heavy-element count; the PDB LYX/MC B29-linked pool has two additional carbon coordinate records, with N/O counts matching source composition. This is an element-count localization, not an atom mapping. The open model's one-oxygen shortfall is the missing ASN A21 OXT terminal record; the closed model contains it. Hydrogen counts differ between states and are not treated as formula evidence.",
            "production_md_ready": False,
            "target_binding_md_started": False,
            "clinical_use": False,
        },
        "next_gate": "Refine the ambiguous B29 source-to-coordinate mappings, resolve the two-carbon excess in the B29-linked PDB pool, and reconcile the remaining hydrogen/link-record differences; then join the audited B1 and B29 components to the insulin backbone and validate the complete graph, protonation, charges, and force-field parameters before target binding MD.",
    }
    OUTPUT_PATH.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT_PATH.relative_to(ROOT)),
        "published_formula": report["published_reference"]["reported_hrms_formula"],
        "states": [
            {
                "state": item["state"],
                "model_inventory": item["model_molecule_atom_inventory"]["formula_from_records"],
                "heavy_delta": item["difference_from_published_formula"]["heavy_element_count_delta_model_minus_source"],
            }
            for item in models
        ],
        "assessment": report["assessment"]["interpretation"],
        "attachment_group_heavy_atom_comparison": report["attachment_group_heavy_atom_comparison"],
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
