#!/usr/bin/env python3
"""Audit source sequence and PDB connection records in public NNC2215 models.

This validates coordinate-file provenance, insulin sequence coverage, backbone
atom presence, native disulfide records, and documented modification-site
links. It does not assign chemical bond orders or validate a force field.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "research" / "data" / "structure-manifest.json"
OUTPUT = ROOT / "research" / "data" / "nnc2215-model-validation.json"

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
EXPECTED_SEQUENCES = {
    "A": "GIVEQCCTSICSLYQLENYCN",
    "B": "FVNQHLCGSHLVEALYLVCGERGFFYTPK",
}
BACKBONE = {"N", "CA", "C", "O"}
EXPECTED_DISULFIDES = {
    frozenset({("A", 6), ("A", 11)}),
    frozenset({("A", 7), ("B", 7)}),
    frozenset({("A", 20), ("B", 19)}),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_pdb(path: Path) -> dict:
    seqres: dict[str, list[str]] = defaultdict(list)
    residues: dict[tuple[str, int, str, str], dict] = {}
    atom_counts_by_chain: Counter[str] = Counter()
    record_counts: Counter[str] = Counter()
    links: list[dict] = []
    disulfides: list[dict] = []
    duplicate_atoms: list[dict] = []
    invalid_coordinates: list[str] = []

    for line_number, line in enumerate(path.read_text(encoding="ascii", errors="replace").splitlines(), 1):
        record = line[:6].strip()
        record_counts[record] += 1

        if record == "SEQRES":
            chain = line[11:12].strip()
            seqres[chain].extend(line[19:70].split())
            continue

        if record in {"ATOM", "HETATM"}:
            atom_name = line[12:16].strip()
            altloc = line[16:17].strip()
            residue_name = line[17:20].strip()
            chain = line[21:22].strip()
            try:
                residue_number = int(line[22:26])
                x, y, z = (float(line[start:end]) for start, end in ((30, 38), (38, 46), (46, 54)))
            except ValueError as error:
                invalid_coordinates.append(f"line {line_number}: {error}")
                continue
            insertion_code = line[26:27].strip()
            key = (chain, residue_number, insertion_code, residue_name)
            residue = residues.setdefault(key, {"atoms": set(), "elements": Counter(), "record": record})
            atom_key = (atom_name, altloc)
            if atom_key in residue["atoms"]:
                duplicate_atoms.append({"residue": list(key), "atom": atom_name, "altloc": altloc})
            residue["atoms"].add(atom_key)
            element = line[76:78].strip().upper()
            if not element:
                element = "".join(character for character in atom_name if character.isalpha())[:1].upper()
            residue["elements"][element] += 1
            atom_counts_by_chain[chain] += 1
            if not all(map(lambda value: value == value and abs(value) != float("inf"), (x, y, z))):
                invalid_coordinates.append(f"line {line_number}: non-finite coordinates")
            continue

        if record == "SSBOND":
            fields = line.split()
            if len(fields) >= 8:
                disulfides.append({
                    "residue_1": (fields[3], int(fields[4])),
                    "residue_2": (fields[6], int(fields[7])),
                    "raw_record": line.rstrip(),
                })
            continue

        if record == "LINK":
            fields = line.split()
            if len(fields) >= 9:
                links.append({
                    "endpoint_1": (fields[3], int(fields[4]), fields[2], fields[1]),
                    "endpoint_2": (fields[7], int(fields[8]), fields[6], fields[5]),
                    "raw_record": line.rstrip(),
                })

    return {
        "seqres": {chain: names for chain, names in sorted(seqres.items())},
        "residues": residues,
        "atom_counts_by_chain": dict(sorted(atom_counts_by_chain.items())),
        "record_counts": dict(sorted(record_counts.items())),
        "links": links,
        "disulfides": disulfides,
        "duplicate_atoms": duplicate_atoms,
        "invalid_coordinates": invalid_coordinates,
    }


def sequence_from_names(names: list[str]) -> str:
    try:
        return "".join(THREE_TO_ONE[name] for name in names)
    except KeyError as error:
        raise ValueError(f"Unknown standard amino-acid name in SEQRES: {error.args[0]}") from error


def model_checks(source: dict, parsed: dict) -> dict:
    path = ROOT / "research" / source["path"]
    seqres = parsed["seqres"]
    residues = parsed["residues"]
    findings: list[str] = []
    sequence_results = {}
    missing_backbone = []
    coordinate_sequences = {}

    for chain, expected in EXPECTED_SEQUENCES.items():
        reported = sequence_from_names(seqres.get(chain, []))
        sequence_results[chain] = {
            "expected": expected,
            "pdb_seqres": reported,
            "matches": reported == expected,
        }
        chain_residues = sorted(
            ((number, name, value) for (res_chain, number, _icode, name), value in residues.items()
             if res_chain == chain and name in THREE_TO_ONE),
            key=lambda item: item[0],
        )
        coordinate_sequence = "".join(THREE_TO_ONE[name] for _, name, _ in chain_residues)
        coordinate_sequences[chain] = coordinate_sequence
        if coordinate_sequence != expected:
            sequence_results[chain]["coordinate_sequence"] = coordinate_sequence
            sequence_results[chain]["coordinate_matches"] = False
        else:
            sequence_results[chain]["coordinate_matches"] = True

        for index, (number, name, residue) in enumerate(chain_residues, 1):
            if number != index:
                findings.append(f"Chain {chain} sequence position {index} is numbered {number} in coordinates.")
            observed = {atom for atom, _altloc in residue["atoms"]}
            missing = sorted(BACKBONE - observed)
            if missing:
                missing_backbone.append({"chain": chain, "residue_number": number,
                                         "residue_name": name, "missing_atoms": missing})

    expected_link_groups = {
        "B29 Lys NZ to linker LYX CD1": frozenset({
            ("B", 29, "LYS", "NZ"), ("B", 30, "LYX", "CD1")
        }),
        "linker to macrocycle roof": frozenset({
            ("B", 30, "LYX", "CD6"), ("C", 1, "MC", "C31")
        }),
        "B1 Phe N to glucoside linker": frozenset({
            ("B", 1, "PHE", "N"), ("B", 0, "PHX", "C2")
        }),
    }
    observed_link_groups = {
        frozenset({tuple(item["endpoint_1"]), tuple(item["endpoint_2"])})
        for item in parsed["links"]
    }
    b29_link_present = expected_link_groups["B29 Lys NZ to linker LYX CD1"] in observed_link_groups
    b1_link_present = expected_link_groups["B1 Phe N to glucoside linker"] in observed_link_groups
    observed_roof_links = [
        item for item in parsed["links"]
        if {item["endpoint_1"][0:3], item["endpoint_2"][0:3]}
        == {("B", 30, "LYX"), ("C", 1, "MC")}
    ]
    roof_atom_names = []
    for item in observed_roof_links:
        roof_atom_names.append({
            "lyx_attachment_atom": next(endpoint[3] for endpoint in (item["endpoint_1"], item["endpoint_2"])
                                          if endpoint[0:3] == ("B", 30, "LYX")),
            "macrocycle_attachment_atom": next(endpoint[3] for endpoint in (item["endpoint_1"], item["endpoint_2"])
                                                 if endpoint[0:3] == ("C", 1, "MC")),
            "record": item["raw_record"],
        })

    observed_disulfides = {
        frozenset({tuple(item["residue_1"]), tuple(item["residue_2"])})
        for item in parsed["disulfides"]
    }
    standard_chain_ids = {chain for chain, *_ in residues if chain in EXPECTED_SEQUENCES}
    nonstandard_residues = []
    for (chain, number, insertion, name), residue in sorted(residues.items()):
        if name not in THREE_TO_ONE:
            nonstandard_residues.append({
                "chain": chain, "residue_number": number, "insertion_code": insertion,
                "residue_name": name,
                "atom_count": len(residue["atoms"]),
                "heavy_atom_count": sum(count for element, count in residue["elements"].items()
                                         if element not in {"H", "D"}),
                "element_counts": dict(sorted(residue["elements"].items())),
            })

    glucose_present = any(name == "BGL" for _chain, _number, _icode, name in residues)
    is_open = source["state"] == "active_open"
    state_ligand_matches = glucose_present if is_open else not glucose_present
    if not state_ligand_matches:
        findings.append("The free-glucose residue label does not match the source state's open/closed designation.")
    if not all(item["matches"] and item["coordinate_matches"] for item in sequence_results.values()):
        findings.append("One or more insulin chains differ from the documented desB30 human insulin sequences.")
    if missing_backbone:
        findings.append("One or more insulin residues are missing backbone atoms.")
    if not b29_link_present:
        findings.append("The PDB LINK records do not contain the expected B29 Lys-NZ to LYX linker connection.")
    if not b1_link_present:
        findings.append("The PDB LINK records do not contain the expected B1 Phe-N to PHX glucoside-linker connection.")
    if observed_disulfides != EXPECTED_DISULFIDES:
        findings.append("The PDB disulfide records differ from the three native insulin disulfide pairs.")
    if parsed["duplicate_atoms"]:
        findings.append("Duplicate atom names were found within a residue.")
    if parsed["invalid_coordinates"]:
        findings.append("Invalid coordinate values were found.")

    return {
        "state": source["state"],
        "file": source["path"],
        "sha256": sha256(path),
        "coordinate_atom_records": sum(parsed["record_counts"].get(k, 0) for k in ("ATOM", "HETATM")),
        "seqres_sequences": sequence_results,
        "coordinate_sequences": coordinate_sequences,
        "backbone_atoms_complete": not missing_backbone,
        "missing_backbone_atoms": missing_backbone,
        "native_disulfide_pairs_match": observed_disulfides == EXPECTED_DISULFIDES,
        "disulfide_records": [item["raw_record"] for item in parsed["disulfides"]],
        "b29_lys_link_present": b29_link_present,
        "b1_phe_glucoside_link_present": b1_link_present,
        "macrocycle_roof_link_atom_names": roof_atom_names,
        "free_glucose_residue_present": glucose_present,
        "open_closed_ligand_state_label_matches": state_ligand_matches,
        "record_counts": parsed["record_counts"],
        "atom_counts_by_chain": parsed["atom_counts_by_chain"],
        "nonstandard_residues": nonstandard_residues,
        "duplicate_atoms": parsed["duplicate_atoms"],
        "invalid_coordinates": parsed["invalid_coordinates"],
        "findings": findings,
        "source_model_interpretation": source["coordinate_status"],
    }


def main() -> None:
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    model_results = []
    for source in manifest["files"]:
        path = ROOT / "research" / source["path"]
        local_digest = sha256(path)
        if local_digest != source["sha256"]:
            raise SystemExit(f"Source manifest SHA-256 mismatch for {source['path']}: {local_digest}")
        model_results.append(model_checks(source, parse_pdb(path)))

    open_model, closed_model = model_results
    open_roof = open_model["macrocycle_roof_link_atom_names"]
    closed_roof = closed_model["macrocycle_roof_link_atom_names"]
    same_roof_atom_name = bool(open_roof and closed_roof) and (
        open_roof[0]["lyx_attachment_atom"] == closed_roof[0]["lyx_attachment_atom"]
        and open_roof[0]["macrocycle_attachment_atom"] == closed_roof[0]["macrocycle_attachment_atom"]
    )
    result = {
        "record_type": "nnc2215_source_model_validation",
        "schema": "dia-health-one.nnc2215-model-validation.v1",
        "run": {
            "run_id": "DHO-NNC2215-MODEL-VALIDATION-" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ"),
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
            "script": "research/validate_nnc2215_models.py",
            "script_sha256": sha256(Path(__file__).resolve()),
            "python_version": platform.python_version(),
            "host_platform": platform.platform(),
            "nvidia_gpu_used": False,
        },
        "source": {
            "dataset_doi": manifest["dataset_doi"],
            "dataset_url": manifest["dataset_url"],
            "license": manifest["license"],
            "article_doi": manifest["related_article"]["doi"],
            "article_url": manifest["related_article"]["url"],
            "chemistry_reference": "The 2024 Nature paper describes the B29 Lys Nε macrocycle linkage and the B1 Phe N-terminal O1-glucoside modification.",
        },
        "expected_insulin_sequence": EXPECTED_SEQUENCES,
        "models": model_results,
        "cross_model_comparison": {
            "macrocycle_roof_attachment_atom_names_match": same_roof_atom_name,
            "open_model_roof_link": open_roof,
            "closed_model_roof_link": closed_roof,
            "interpretation": (
                "The open and closed PDBs use different LYX atom names for the linker-to-macrocycle connection. "
                "This may be representational, but it must be reconciled to the chemical source before assigning a shared topology."
            ),
        },
        "validation_summary": {
            "source_hashes_match_manifest": True,
            "insulin_sequences_match_desB30_human_insulin": all(
                all(item["matches"] and item["coordinate_matches"] for item in model["seqres_sequences"].values())
                for model in model_results
            ),
            "backbone_atoms_complete": all(model["backbone_atoms_complete"] for model in model_results),
            "three_native_disulfides_match": all(model["native_disulfide_pairs_match"] for model in model_results),
            "b29_macrocycle_and_b1_glucoside_link_records_present": all(
                model["b29_lys_link_present"] and model["b1_phe_glucoside_link_present"]
                for model in model_results
            ),
            "chemical_graph_and_bond_orders_validated": False,
            "protonation_and_partial_charges_validated": False,
            "force_field_parameters_validated": False,
            "insulin_receptor_coordinates_present": False,
            "production_simulation_ready": False,
        },
        "readiness_boundary": (
            "Source, sequence, backbone, disulfide, and modification-site records were checked. "
            "These are author-supplied model-built NNC2215 conformations, not experimentally resolved NNC2215 coordinates. "
            "No receptor coordinates are included. The nonstandard LYX/PHX/MC chemical graphs, bond orders, stereochemistry, "
            "protonation, charges, and force-field parameters remain unvalidated; no target simulation is authorized by this record."
        ),
        "models_have_findings": any(model["findings"] for model in model_results),
        "cross_model_topology_requires_review": not same_roof_atom_name,
        "clinical_use": False,
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT.relative_to(ROOT)),
        "run_id": result["run"]["run_id"],
        "validation_summary": result["validation_summary"],
        "open_roof_link": open_roof,
        "closed_roof_link": closed_roof,
        "cross_model_topology_requires_review": result["cross_model_topology_requires_review"],
        "findings": {model["state"]: model["findings"] for model in model_results},
    }, indent=2))


if __name__ == "__main__":
    main()
