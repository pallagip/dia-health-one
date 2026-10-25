#!/usr/bin/env python3
"""Validate the source-derived NNC2215 roof macrocycle graph and write its audit."""

from __future__ import annotations

import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

from rdkit import Chem, rdBase
from rdkit.Chem import rdMolDescriptors


ROOT = Path(__file__).resolve().parent
GRAPH_PATH = ROOT / "structures" / "nnc2215-roof-propyl-azide.smi"
MAP_PATH = ROOT / "data" / "nnc2215-roof-macrocycle-atom-map.json"
AUDIT_PATH = ROOT / "data" / "nnc2215-roof-macrocycle-graph-audit.json"
EXPECTED_FORMULA = "C56H63N15O12"
EXPECTED_HEAVY_ATOMS = 83
EXPECTED_MAP_IDS = set(range(1, EXPECTED_HEAVY_ATOMS + 1))


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def load_graph() -> tuple[str, str, bytes]:
    raw = GRAPH_PATH.read_bytes()
    lines = [line.strip() for line in raw.decode("utf-8").splitlines() if line.strip() and not line.lstrip().startswith("#")]
    if len(lines) != 1:
        raise ValueError(f"Expected exactly one non-comment SMILES line, found {len(lines)}")
    fields = lines[0].split()
    if len(fields) != 2:
        raise ValueError("Expected a SMILES string followed by one molecule name")
    return fields[0], fields[1], raw


def main() -> int:
    try:
        smiles, molecule_name, graph_bytes = load_graph()
        params = Chem.SmilesParserParams()
        params.sanitize = False
        mol = Chem.MolFromSmiles(smiles, params)
        if mol is None:
            raise ValueError("RDKit could not parse the mapped SMILES")
        Chem.SanitizeMol(mol)

        formula = rdMolDescriptors.CalcMolFormula(mol)
        heavy_atoms = mol.GetNumHeavyAtoms()
        formal_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
        mapped_heavy_atoms = [atom.GetAtomMapNum() for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
        mapped_id_set = set(mapped_heavy_atoms)
        duplicate_ids = sorted({map_id for map_id in mapped_heavy_atoms if mapped_heavy_atoms.count(map_id) > 1})
        missing_ids = sorted(EXPECTED_MAP_IDS - mapped_id_set)
        unexpected_ids = sorted(mapped_id_set - EXPECTED_MAP_IDS)

        failures: list[str] = []
        if formula != EXPECTED_FORMULA:
            failures.append(f"formula {formula} != {EXPECTED_FORMULA}")
        if heavy_atoms != EXPECTED_HEAVY_ATOMS:
            failures.append(f"heavy-atom count {heavy_atoms} != {EXPECTED_HEAVY_ATOMS}")
        if len(mapped_heavy_atoms) != EXPECTED_HEAVY_ATOMS:
            failures.append(f"mapped heavy-atom count {len(mapped_heavy_atoms)} != {EXPECTED_HEAVY_ATOMS}")
        if duplicate_ids:
            failures.append(f"duplicate atom-map IDs: {duplicate_ids}")
        if missing_ids or unexpected_ids:
            failures.append(f"map IDs differ from 1..83; missing={missing_ids}, unexpected={unexpected_ids}")
        if formal_charge != 0:
            failures.append(f"formal charge {formal_charge} != 0 for the source-drawn neutral-acid representation")
        if failures:
            raise ValueError("; ".join(failures))

        map_data = json.loads(MAP_PATH.read_text(encoding="utf-8"))
        canonical_mapped_smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
        chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True, includeCIP=True)
        audit = {
            "schema": "dia-health-one.nnc2215-roof-macrocycle-graph-audit.v1",
            "generated_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z"),
            "molecule_name": molecule_name,
            "source": map_data["source"],
            "input": {
                "graph_file": str(GRAPH_PATH.relative_to(ROOT)),
                "graph_sha256": sha256(graph_bytes),
                "atom_map_file": str(MAP_PATH.relative_to(ROOT)),
                "atom_map_sha256": sha256(MAP_PATH.read_bytes()),
                "canonical_mapped_smiles": canonical_mapped_smiles,
            },
            "software": {"rdkit_version": rdBase.rdkitVersion, "python_version": sys.version.split()[0]},
            "chemistry": {
                "formula": formula,
                "heavy_atom_count": heavy_atoms,
                "formal_charge": formal_charge,
                "ring_count_sssr": mol.GetRingInfo().NumRings(),
                "atom_map_count": len(mapped_heavy_atoms),
                "atom_map_ids_unique_and_complete_1_to_83": True,
                "specified_chiral_centers": [
                    {"atom_index": index, "cip_code": code} for index, code in chiral_centers if code != "?"
                ],
                "unassigned_potential_chiral_centers": [
                    index for index, code in chiral_centers if code == "?"
                ],
            },
            "representation_assumptions": map_data["representation_assumptions"],
            "limitations": map_data["limitations"],
            "checks": {
                "rdkit_parse_and_sanitize": "pass",
                "formula_matches_source": "pass",
                "heavy_atom_count_matches_source_graph": "pass",
                "positive_unique_complete_atom_map_ids": "pass",
                "neutral_formal_charge_in_depicted_form": "pass",
            },
            "readiness": {
                "isolated_macrocycle_graph_candidate": True,
                "complete_nnc2215_conjugate_graph": False,
                "pdb_atom_mapping_validated": False,
                "simulation_protonation_and_charges_validated": False,
                "force_field_parameters_validated": False,
                "production_md_ready": False,
            },
            "clinical_use": False,
            "nvidia_gpu_used": False,
        }
        AUDIT_PATH.write_text(json.dumps(audit, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        print(json.dumps({"audit": str(AUDIT_PATH), "formula": formula, "heavy_atoms": heavy_atoms, "map_ids": len(mapped_heavy_atoms), "formal_charge": formal_charge, "rdkit": rdBase.rdkitVersion}, indent=2))
        return 0
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
