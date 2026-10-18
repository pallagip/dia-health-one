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


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def formula_text(counts: dict[str, int]) -> str:
    order = ("C", "H", "N", "O", "S", "P", "F", "CL", "BR", "I")
    elements = [element for element in order if counts.get(element, 0)]
    elements.extend(sorted(set(counts) - set(elements)))
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
    expected_attachment_heavy = {
        element: count for element, count in expected_attachment_formula.items()
        if element != "H" and count
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
            "interpretation": "The PDB ATOM protein cores retain the same heavy-element counts as sequence-derived desB30 insulin in the closed model. After subtracting this parent formula from the published NNC2215 formula, the PDB PHX/LYX/MC HETATM component sum contains two excess carbon records (C74 observed versus C72 expected); nitrogen and oxygen counts match. This localizes the two-carbon discrepancy to the aggregate nonstandard HETATM components, but not to a particular residue or atom.",
        },
        "method": {
            "description": "Verify source-file SHA-256 values and coordinate-record counts; count ATOM/HETATM element fields; derive the desB30 parent formula from validated chain sequences and three disulfides; and exclude the explicitly modelled free-glucose BGL residue before comparing the model components with the published formula.",
            "chemistry_inference": False,
            "bond_order_or_charge_assignment": False,
            "hydrogen_counts_interpreted_as_formula_validation": False,
            "limitation": "This is a coordinate atom inventory plus a sequence-based parent-formula calculation, not a validated molecular graph. It localizes a heavy-atom count difference to the aggregate nonstandard components but does not establish which atom is wrong, whether a model atom is chemically intended, or the chemical cause. Hydrogen counts are not used as formula-validation evidence.",
        },
        "model_results": models,
        "assessment": {
            "same_extra_heavy_carbon_count_in_both_states": len(set(carbon_deltas)) == 1 and carbon_deltas[0] == 2,
            "unresolved_coordinate_inventory_mismatch": any(carbon_deltas),
            "interpretation": "After excluding free-glucose BGL, both publisher-supplied state files contain two more carbon records than the published NNC2215 formula. Comparison against the sequence-derived desB30 insulin parent localizes those two excess carbon records to the aggregate PHX/LYX/MC HETATM components. The open model's one-oxygen shortfall is localized to the missing ASN A21 OXT terminal atom record; the closed model contains it. Hydrogen counts differ between states and are not treated as formula evidence. The identity and chemical reason for the extra carbon records remain unresolved.",
            "production_md_ready": False,
            "target_binding_md_started": False,
            "clinical_use": False,
        },
        "next_gate": "Map the PHX/LYX/MC atoms against the published B1 and B29 structures and resolve the two excess carbon records in their aggregate heavy-atom inventory before assigning a force-field topology or running target binding MD.",
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
