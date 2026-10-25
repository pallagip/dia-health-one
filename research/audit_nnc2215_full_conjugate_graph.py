#!/usr/bin/env python3
"""Assemble and audit the source-derived neutral NNC2215 conjugate graph.

The constitutional graph is built from pinned RDKit FASTA residue templates,
the mapped B1/B29 source SMILES, and the three declared insulin disulfides.
Corrected candidate PDBs are read only for sequence, SSBOND, LINK, and input
integrity cross-checks. Coordinates are never used to build or map this graph.
"""

from __future__ import annotations

from collections import Counter
import hashlib
import json
from pathlib import Path
import platform
import re
from typing import Any

from rdkit import Chem, rdBase
from rdkit.Chem import rdMolDescriptors


ROOT = Path(__file__).resolve().parents[1]
RESEARCH = ROOT / "research"
DATA = RESEARCH / "data"
STRUCTURES = RESEARCH / "structures"
OUTPUT_PATH = DATA / "nnc2215-full-conjugate-graph.json"

PATHS = {
    "formula_audit": DATA / "nnc2215-model-formula-audit.json",
    "b1_graph": STRUCTURES / "nnc2215-b1-glucoside-attachment.smi",
    "b1_map_legend": DATA / "nnc2215-b1-glucoside-atom-map.json",
    "b1_topology_audit": DATA / "nnc2215-b1-topology-audit.json",
    "b29_graph": STRUCTURES / "nnc2215-b29-linker-roof-product.smi",
    "b29_map_legend": DATA / "nnc2215-b29-linker-roof-product-map.json",
    "b29_topology_audit": DATA / "nnc2215-b29-topology-audit.json",
    "reconciliation": DATA / "nnc2215-carbon-reconciliation-candidate.json",
    "boundary_audit": DATA / "nnc2215-attachment-boundary-audit.json",
    "conditional_carboxylate_audit": DATA / "nnc2215-b29-topology-audit.json",
    "active_open_source_pdb": STRUCTURES / "nnc2215-open-model.pdb",
    "inactive_closed_source_pdb": STRUCTURES / "nnc2215-closed-model.pdb",
    "active_open_candidate_pdb": STRUCTURES / "nnc2215-open-source-corrected-candidate.pdb",
    "inactive_closed_candidate_pdb": STRUCTURES / "nnc2215-closed-source-corrected-candidate.pdb",
}

STATE_ORDER = ("active_open", "inactive_closed")
STATE_DATA = {
    "active_open": {
        "source_path_key": "active_open_source_pdb",
        "candidate_path_key": "active_open_candidate_pdb",
        "formula_state": "active_open",
    },
    "inactive_closed": {
        "source_path_key": "inactive_closed_source_pdb",
        "candidate_path_key": "inactive_closed_candidate_pdb",
        "formula_state": "inactive_closed",
    },
}

EXPECTED_SEQUENCES = {
    "A": "GIVEQCCTSICSLYQLENYCN",
    "B": "FVNQHLCGSHLVEALYLVCGERGFFYTPK",
}
ONE_TO_THREE = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS",
    "Q": "GLN", "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE",
    "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
    "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
}
THREE_TO_ONE = {three: one for one, three in ONE_TO_THREE.items()}
PUBLISHED_FORMULA = "C325H461N79O97S6"
PUBLISHED_COUNTS = {"C": 325, "H": 461, "N": 79, "O": 97, "S": 6}
DISULFIDES = (
    (("A", 6, "SG"), ("A", 11, "SG")),
    (("A", 7, "SG"), ("B", 7, "SG")),
    (("A", 20, "SG"), ("B", 19, "SG")),
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(ROOT).as_posix()


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object in {relative(path)}")
    return value


def formula_counts_text(counts: dict[str, int]) -> str:
    order = ("C", "H", "N", "O", "S", "P", "F", "CL", "BR", "I")
    elements = [element for element in order if counts.get(element, 0)]
    elements.extend(sorted(set(counts) - set(elements)))
    return "".join(
        element + (str(counts[element]) if counts[element] != 1 else "")
        for element in elements
        if counts.get(element, 0)
    )


def smiles_from_file(path: Path) -> str:
    tokens = path.read_text(encoding="utf-8").split()
    if not tokens:
        raise ValueError(f"Empty mapped SMILES in {relative(path)}")
    return tokens[0]


def atom_map_index(molecule: Chem.Mol, label: str) -> dict[int, int]:
    result: dict[int, int] = {}
    for atom in molecule.GetAtoms():
        map_id = atom.GetAtomMapNum()
        if map_id <= 0 or map_id in result:
            raise ValueError(f"Invalid or duplicated atom map {map_id} in {label}")
        result[map_id] = atom.GetIdx()
    return result


def atom_map_roles(legend: dict[str, Any], component: str) -> dict[int, str]:
    roles: dict[int, str] = {}
    if component == "B1":
        for raw_id, record in legend.get("map_id_legend", {}).items():
            if isinstance(record, dict):
                roles[int(raw_id)] = str(record.get("role", "mapped B1 source atom"))
            else:
                roles[int(raw_id)] = str(record)
        return roles

    map_legend = legend.get("map_id_legend", {})

    def collect(value: Any, path: str) -> None:
        if isinstance(value, dict):
            for key, child in value.items():
                if str(key).isdigit():
                    roles[int(key)] = str(child)
                else:
                    collect(child, f"{path}/{key}")
        elif isinstance(value, list):
            for child in value:
                collect(child, path)
        elif isinstance(value, int):
            roles[value] = path

    collect(map_legend, "B29 product source map")
    for map_id in range(1, 93):
        roles.setdefault(map_id, "B29 product source map")
    return roles


def validate_hashes(
    formula_audit: dict[str, Any],
    b1_audit: dict[str, Any],
    b29_audit: dict[str, Any],
    reconciliation: dict[str, Any],
    boundary_audit: dict[str, Any],
) -> tuple[dict[str, Any], dict[str, str]]:
    """Verify source/candidate inputs against all recorded upstream digests."""
    expected: dict[str, list[dict[str, str]]] = {key: [] for key in PATHS}

    def add(key: str, value: Any, source: str) -> None:
        if isinstance(value, str) and re.fullmatch(r"[0-9a-f]{64}", value):
            record = {"sha256": value, "recorded_by": source}
            if record not in expected[key]:
                expected[key].append(record)

    b1_hashes = b1_audit.get("input_hashes_sha256", {})
    add("b1_graph", b1_hashes.get("source_graph", {}).get("sha256"), "B1 topology audit")
    add("b1_map_legend", b1_hashes.get("map_legend", {}).get("sha256"), "B1 topology audit")
    add("reconciliation", b1_hashes.get("coordinate_reconciliation", {}).get("sha256"), "B1 topology audit")

    b29_hash_records = b29_audit.get("inputs", {}).get("sha256_file_integrity_records", {})
    add("b29_graph", b29_audit.get("source_product_graph", {}).get("graph_sha256"), "B29 topology audit")
    add("b29_map_legend", b29_audit.get("source_product_graph", {}).get("map_legend_sha256"), "B29 topology audit")
    add("b29_graph", b29_hash_records.get("product_smiles", {}).get("sha256_expected"), "B29 topology audit input integrity")
    add("reconciliation", b29_hash_records.get("coordinate_reconciliation_record", {}).get("sha256_expected"), "B29 topology audit input integrity")

    boundary_records = boundary_audit.get("inputs", {}).get("sha256_file_integrity_records", {})
    boundary_keys = {
        "b1_graph": "b1_source_graph",
        "b1_map_legend": "b1_map_legend",
        "b29_graph": "b29_source_graph",
        "b29_map_legend": "b29_map_legend",
        "reconciliation": "coordinate_reconciliation",
    }
    for path_key, audit_key in boundary_keys.items():
        record = boundary_records.get(audit_key, {})
        add(path_key, record.get("sha256_expected"), "attachment-boundary audit")
        if record.get("sha256_expected") is None:
            add(path_key, record.get("sha256_observed"), "attachment-boundary audit observed digest")

    reconciliation_models = {
        item.get("state"): item for item in reconciliation.get("models", [])
    }
    formula_models = {item.get("state"): item for item in formula_audit.get("model_results", [])}
    b1_states = {item.get("state"): item for item in b1_audit.get("candidate_states", [])}
    b29_states = {item.get("state"): item for item in b29_audit.get("states", [])}
    boundary_candidate_hashes = boundary_audit.get("inputs", {}).get("candidate_pdb_before_after_hashes", {})

    for state in STATE_ORDER:
        details = STATE_DATA[state]
        reconciliation_model = reconciliation_models.get(state, {})
        formula_model = formula_models.get(state, {})
        add(details["source_path_key"], reconciliation_model.get("source_sha256"), "coordinate reconciliation")
        add(details["source_path_key"], formula_model.get("sha256"), "model formula audit")
        add(details["candidate_path_key"], reconciliation_model.get("output_sha256"), "coordinate reconciliation")
        add(details["candidate_path_key"], b1_states.get(state, {}).get("candidate_pdb_sha256"), "B1 topology audit")
        b29_state = b29_states.get(state, {})
        add(details["candidate_path_key"], b29_state.get("coordinate_candidate_sha256_expected_by_reconciliation"), "B29 topology audit")
        boundary_state = boundary_candidate_hashes.get(state, {})
        add(details["candidate_path_key"], boundary_state.get("sha256_expected_by_reconciliation"), "attachment-boundary audit")
        add(details["candidate_path_key"], boundary_state.get("sha256_after_parsing"), "attachment-boundary audit post-parse")
        b1_candidate = next((item for item in b1_audit.get("input_hashes_sha256", {}).get("corrected_candidate_pdbs", []) if item.get("state") == state), {})
        add(details["candidate_path_key"], b1_candidate.get("observed_sha256"), "B1 topology audit input hash")
        b29_integrity_key = f"{state}_corrected_candidate_pdb"
        add(details["candidate_path_key"], b29_hash_records.get(b29_integrity_key, {}).get("sha256_expected"), "B29 topology audit input integrity")

    hash_records: dict[str, Any] = {}
    observed: dict[str, str] = {}
    for key, path in PATHS.items():
        actual = sha256(path)
        observed[key] = actual
        expected_values = expected[key]
        comparisons = [
            {**record, "matches": record["sha256"] == actual}
            for record in expected_values
        ]
        if expected_values and any(not item["matches"] for item in comparisons):
            sources = ", ".join(item["recorded_by"] for item in comparisons if not item["matches"])
            raise ValueError(f"SHA-256 mismatch for {relative(path)} against {sources}")
        hash_records[key] = {
            "path": relative(path),
            "observed_sha256": actual,
            "upstream_expected_sha256": comparisons,
            "all_recorded_expected_hashes_match": all(item["matches"] for item in comparisons),
            "integrity_basis": "upstream expected digests" if comparisons else "observed digest only",
        }

    audit_paths = {
        "formula_audit": PATHS["formula_audit"],
        "b1_topology_audit": PATHS["b1_topology_audit"],
        "b29_topology_audit": PATHS["b29_topology_audit"],
        "reconciliation": PATHS["reconciliation"],
        "boundary_audit": PATHS["boundary_audit"],
    }
    hash_records["audit_record_files"] = {
        key: {"path": relative(path), "observed_sha256": sha256(path)}
        for key, path in audit_paths.items()
    }
    return hash_records, observed


def make_chain(sequence: str, chain: str) -> Chem.Mol:
    molecule = Chem.MolFromFASTA(sequence)
    if molecule is None:
        raise ValueError(f"RDKit MolFromFASTA failed for insulin chain {chain}")
    for atom in molecule.GetAtoms():
        info = atom.GetPDBResidueInfo()
        if info is None:
            raise ValueError(f"RDKit FASTA atom {atom.GetIdx()} has no residue metadata")
        atom.SetProp("_NNCOrigin", "protein_template")
        atom.SetProp("_NNCChain", chain)
        atom.SetProp("_NNCResidueNumber", str(info.GetResidueNumber()))
        atom.SetProp("_NNCResidueName", info.GetResidueName().strip())
        atom.SetProp("_NNCAtomName", info.GetName().strip())
        atom.SetProp("_NNCSource", "RDKit MolFromFASTA neutral amino-acid template")
    for bond in molecule.GetBonds():
        bond.SetProp("_NNCProvenance", f"RDKit MolFromFASTA peptide/residue template, chain {chain}")
    return molecule


def find_template_atom(molecule: Chem.Mol, chain: str, residue: int, name: str) -> int:
    matches = [
        atom.GetIdx() for atom in molecule.GetAtoms()
        if atom.HasProp("_NNCChain")
        and atom.GetProp("_NNCChain") == chain
        and atom.HasProp("_NNCResidueNumber")
        and int(atom.GetProp("_NNCResidueNumber")) == residue
        and atom.HasProp("_NNCAtomName")
        and atom.GetProp("_NNCAtomName") == name
    ]
    if len(matches) != 1:
        raise ValueError(f"Expected unique FASTA atom {chain}:{residue}:{name}; found {len(matches)}")
    return matches[0]


def parse_mapped_source(path: Path, component: str, legend: dict[str, Any]) -> tuple[Chem.Mol, dict[int, int], dict[int, str]]:
    molecule = Chem.MolFromSmiles(smiles_from_file(path))
    if molecule is None:
        raise ValueError(f"RDKit could not parse {component} mapped SMILES")
    Chem.SanitizeMol(molecule)
    index = atom_map_index(molecule, component)
    roles = atom_map_roles(legend, component)
    for atom in molecule.GetAtoms():
        map_id = atom.GetAtomMapNum()
        atom.SetProp("_NNCOrigin", "mapped_source_fragment")
        atom.SetProp("_NNCComponent", component)
        atom.SetProp("_NNCSourceMap", str(map_id))
        atom.SetProp("_NNCSource", relative(path))
        atom.SetProp("_NNCMapRole", roles.get(map_id, f"{component} mapped source atom"))
    for bond in molecule.GetBonds():
        left = molecule.GetAtomWithIdx(bond.GetBeginAtomIdx()).GetAtomMapNum()
        right = molecule.GetAtomWithIdx(bond.GetEndAtomIdx()).GetAtomMapNum()
        bond.SetProp("_NNCProvenance", f"Mapped {component} source SMILES edge map {left}–{right}")
    return molecule, index, roles


def remove_source_cap(molecule: Chem.Mol, cap_map: int, endpoint_map: int, component: str) -> Chem.Mol:
    maps = atom_map_index(molecule, component)
    if cap_map not in maps or endpoint_map not in maps:
        raise ValueError(f"{component} cap map {cap_map} or endpoint map {endpoint_map} is missing")
    cap_atom = molecule.GetAtomWithIdx(maps[cap_map])
    neighbor_maps = sorted(
        atom.GetAtomMapNum() for atom in cap_atom.GetNeighbors()
    )
    if neighbor_maps != [endpoint_map]:
        raise ValueError(f"{component} cap map {cap_map} must be adjacent only to map {endpoint_map}; got {neighbor_maps}")
    rw_mol = Chem.RWMol(molecule)
    rw_mol.RemoveAtom(maps[cap_map])
    result = rw_mol.GetMol()
    remaining = atom_map_index(result, component)
    if endpoint_map not in remaining or cap_map in remaining:
        raise ValueError(f"Failed to remove only {component} cap map {cap_map}")
    return result


def attach_fragment(
    combined: Chem.Mol,
    fragment: Chem.Mol,
    endpoint_map: int,
    protein_atom_index: int,
    component: str,
) -> Chem.Mol:
    offset = combined.GetNumAtoms()
    merged = Chem.CombineMols(combined, fragment)
    endpoint_indices = [
        atom.GetIdx() for atom in merged.GetAtoms()
        if atom.HasProp("_NNCSourceMap") and int(atom.GetProp("_NNCSourceMap")) == endpoint_map
    ]
    if len(endpoint_indices) != 1:
        raise ValueError(f"Expected one appended {component} map {endpoint_map}; found {len(endpoint_indices)}")
    endpoint_index = endpoint_indices[0]
    if endpoint_index < offset:
        raise ValueError(f"{component} source map {endpoint_map} was not appended")
    rw_mol = Chem.RWMol(merged)
    rw_mol.AddBond(protein_atom_index, endpoint_index, Chem.BondType.SINGLE)
    bond = rw_mol.GetBondBetweenAtoms(protein_atom_index, endpoint_index)
    bond.SetProp(
        "_NNCProvenance",
        f"{component} protein attachment boundary; source cap replaced at map {endpoint_map}",
    )
    return rw_mol.GetMol()


def formula_from_atom_ledger(molecule: Chem.Mol) -> tuple[dict[str, int], dict[str, int]]:
    full: Counter[str] = Counter()
    heavy: Counter[str] = Counter()
    for atom in molecule.GetAtoms():
        if atom.GetAtomicNum() == 1:
            raise ValueError("Hydrogens must remain implicit in the graph ledger")
        element = atom.GetSymbol().upper()
        heavy[element] += 1
        full[element] += 1
        full["H"] += atom.GetTotalNumHs(includeNeighbors=False)
    return dict(sorted(full.items())), dict(sorted(heavy.items()))


def atom_record(atom: Chem.Atom) -> dict[str, Any]:
    explicit_valence = atom.GetValence(Chem.ValenceType.EXPLICIT)
    implicit_valence = atom.GetValence(Chem.ValenceType.IMPLICIT)
    record: dict[str, Any] = {
        "node_id": f"atom_{atom.GetIdx() + 1:04d}",
        "atom_index_zero_based": atom.GetIdx(),
        "element": atom.GetSymbol(),
        "atomic_number": atom.GetAtomicNum(),
        "formal_charge": atom.GetFormalCharge(),
        "is_aromatic": atom.GetIsAromatic(),
        "implicit_hydrogen_count": atom.GetNumImplicitHs(),
        "bracket_hydrogen_count": atom.GetNumExplicitHs(),
        "hydrogen_count_not_as_nodes": atom.GetTotalNumHs(includeNeighbors=False),
        "explicit_valence": explicit_valence,
        "total_valence": explicit_valence + implicit_valence,
    }
    if atom.HasProp("_NNCOrigin") and atom.GetProp("_NNCOrigin") == "protein_template":
        record["provenance"] = {
            "origin": "insulin_chain_template",
            "template": "RDKit MolFromFASTA neutral amino-acid template",
            "chain": atom.GetProp("_NNCChain"),
            "residue_number": int(atom.GetProp("_NNCResidueNumber")),
            "residue_name": atom.GetProp("_NNCResidueName"),
            "atom_name": atom.GetProp("_NNCAtomName"),
            "source": "sequence from research/data/nnc2215-model-formula-audit.json",
        }
    elif atom.HasProp("_NNCOrigin") and atom.GetProp("_NNCOrigin") == "mapped_source_fragment":
        record["provenance"] = {
            "origin": "mapped_source_fragment",
            "component": atom.GetProp("_NNCComponent"),
            "source_map_id": int(atom.GetProp("_NNCSourceMap")),
            "map_role": atom.GetProp("_NNCMapRole"),
            "source": atom.GetProp("_NNCSource"),
        }
    else:
        raise ValueError(f"Atom {atom.GetIdx()} has no recognized provenance")
    return record


def bond_record(bond: Chem.Bond) -> dict[str, Any]:
    if bond.HasProp("_NNCProvenance"):
        provenance = bond.GetProp("_NNCProvenance")
    else:
        left = bond.GetBeginAtom()
        right = bond.GetEndAtom()
        if left.HasProp("_NNCChain") and right.HasProp("_NNCChain"):
            provenance = f"RDKit MolFromFASTA neutral amino-acid template, chain {left.GetProp('_NNCChain')}"
        elif left.HasProp("_NNCComponent") and right.HasProp("_NNCComponent"):
            provenance = f"Mapped {left.GetProp('_NNCComponent')} source SMILES"
        else:
            raise ValueError(f"Bond {bond.GetIdx()} is missing provenance")
    return {
        "edge_id": f"bond_{bond.GetIdx() + 1:04d}",
        "atom_node_ids": [
            f"atom_{bond.GetBeginAtomIdx() + 1:04d}",
            f"atom_{bond.GetEndAtomIdx() + 1:04d}",
        ],
        "atom_indices_zero_based": [bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()],
        "bond_type": str(bond.GetBondType()),
        "bond_order": bond.GetBondTypeAsDouble(),
        "is_aromatic": bond.GetIsAromatic(),
        "is_conjugated": bond.GetIsConjugated(),
        "provenance": provenance,
    }


def parse_endpoint(line: str, second: bool = False) -> dict[str, Any]:
    if second:
        start = 42
    else:
        start = 12
    return {
        "atom_name": line[start : start + 4].strip(),
        "altloc": line[start + 4].strip(),
        "residue_name": line[start + 5 : start + 8].strip(),
        "chain": line[start + 9].strip(),
        "residue_number": int(line[start + 10 : start + 14]),
        "insertion_code": line[start + 14].strip(),
    }


def endpoint_key(record: dict[str, Any]) -> tuple[Any, ...]:
    return (
        record.get("chain", ""),
        int(record.get("residue_number", 0)),
        record.get("insertion_code", ""),
        record.get("residue_name", ""),
        record.get("atom_name", record.get("name", "")),
        record.get("altloc", ""),
    )


def parse_pdb(path: Path) -> dict[str, Any]:
    atoms: list[dict[str, Any]] = []
    residue_names: dict[str, dict[tuple[int, str], str]] = {"A": {}, "B": {}}
    links: list[dict[str, Any]] = []
    ssbonds: list[dict[str, Any]] = []
    oxt_records: list[dict[str, Any]] = []
    for line_number, line in enumerate(path.read_text(encoding="ascii").splitlines(), start=1):
        record_type = line[:6].strip()
        if record_type in {"ATOM", "HETATM"}:
            atom_name = line[12:16].strip()
            residue_name = line[17:20].strip()
            chain = line[21].strip()
            residue_number = int(line[22:26])
            insertion_code = line[26].strip()
            atom = {
                "record_type": record_type,
                "serial": int(line[6:11]),
                "atom_name": atom_name,
                "altloc": line[16].strip(),
                "residue_name": residue_name,
                "chain": chain,
                "residue_number": residue_number,
                "insertion_code": insertion_code,
                "element": line[76:78].strip().upper(),
            }
            atoms.append(atom)
            if record_type == "ATOM" and chain in residue_names:
                key = (residue_number, insertion_code)
                prior = residue_names[chain].get(key)
                if prior is not None and prior != residue_name:
                    raise ValueError(f"Conflicting residue names at {chain}:{residue_number} in {relative(path)}")
                residue_names[chain][key] = residue_name
            if atom_name == "OXT":
                oxt_records.append({
                    "chain": chain,
                    "residue_name": residue_name,
                    "residue_number": residue_number,
                    "element": atom["element"],
                })
        elif record_type == "LINK":
            links.append({
                "line_number": line_number,
                "first_endpoint": parse_endpoint(line),
                "second_endpoint": parse_endpoint(line, second=True),
                "raw_record": line.rstrip(),
            })
        elif record_type == "SSBOND":
            ssbonds.append({
                "line_number": line_number,
                "first": {
                    "residue_name": line[11:14].strip(),
                    "chain": line[15].strip(),
                    "residue_number": int(line[17:21]),
                    "insertion_code": line[21].strip(),
                },
                "second": {
                    "residue_name": line[25:28].strip(),
                    "chain": line[29].strip(),
                    "residue_number": int(line[31:35]),
                    "insertion_code": line[35].strip(),
                },
                "raw_record": line.rstrip(),
            })
    sequences = {
        chain: "".join(THREE_TO_ONE[residue] for _, residue in sorted(residues.items()))
        for chain, residues in residue_names.items()
    }
    return {
        "atoms": atoms,
        "sequences": sequences,
        "links": links,
        "ssbonds": ssbonds,
        "oxt_records": oxt_records,
    }


def pdb_disulfide_pairs(records: list[dict[str, Any]]) -> set[tuple[tuple[str, int], tuple[str, int]]]:
    pairs = set()
    for record in records:
        first = record["first"]
        second = record["second"]
        left = (first["chain"], first["residue_number"])
        right = (second["chain"], second["residue_number"])
        pairs.add(tuple(sorted((left, right))))
    return pairs


def expected_disulfide_pairs() -> set[tuple[tuple[str, int], tuple[str, int]]]:
    return {
        tuple(sorted(((left[0], left[1]), (right[0], right[1]))))
        for left, right in DISULFIDES
    }


def boundary_link_cross_checks(
    parsed_pdb: dict[str, Any],
    boundary_state: dict[str, Any],
    state: str,
) -> list[dict[str, Any]]:
    joins_by_id = {
        item.get("join_id"): item
        for item in boundary_state.get("attachment_boundary_joins", [])
    }
    required = {
        "b1_glucoside_attachment": [1, 2],
        "b29_lys_attachment": [84, 92],
        "roof_tail_attachment": [78, 79],
    }
    results: list[dict[str, Any]] = []
    for join_id, expected_maps in required.items():
        join = joins_by_id.get(join_id)
        if join is None:
            raise ValueError(f"Boundary audit lacks {join_id} for {state}")
        source_edge = join.get("source_graph_edge", {})
        if source_edge.get("map_ids") != expected_maps:
            raise ValueError(f"Unexpected source map edge for {join_id} in {state}")
        identity_records = join.get("atom_identities", [])
        if len(identity_records) != 2:
            raise ValueError(f"Expected two linked atom identities for {join_id} in {state}")
        expected_pair = tuple(sorted(endpoint_key(item) for item in identity_records))
        matching_links = [
            link for link in parsed_pdb["links"]
            if tuple(sorted((
                endpoint_key(link["first_endpoint"]),
                endpoint_key(link["second_endpoint"]),
            ))) == expected_pair
        ]
        if len(matching_links) != 1:
            raise ValueError(f"Expected one PDB LINK for {join_id} in {state}; found {len(matching_links)}")
        if len(join.get("link_records", [])) != 1:
            raise ValueError(f"Boundary audit should retain one LINK record for {join_id} in {state}")
        results.append({
            "join_id": join_id,
            "source_map_ids": expected_maps,
            "pdb_link_record": matching_links[0]["raw_record"],
            "pdb_endpoints": [matching_links[0]["first_endpoint"], matching_links[0]["second_endpoint"]],
            "boundary_audit_identity_match": True,
            "coordinate_values_used": False,
        })
    return results


def build_full_graph(
    formula_audit: dict[str, Any],
    b1_legend: dict[str, Any],
    b29_legend: dict[str, Any],
) -> tuple[Chem.Mol, dict[str, Any]]:
    observed_sequences = []
    for model in formula_audit.get("model_results", []):
        sequences = model.get("model_molecule_atom_inventory", {}).get("protein_atom_records", {}).get("chain_sequences")
        if sequences:
            observed_sequences.append(sequences)
    if not observed_sequences or any(sequences != EXPECTED_SEQUENCES for sequences in observed_sequences):
        raise ValueError(f"Formula audit insulin sequences differ from required A/B sequences: {observed_sequences}")
    chain_sequences = observed_sequences[0]

    chain_a = make_chain(chain_sequences["A"], "A")
    chain_b = make_chain(chain_sequences["B"], "B")
    combined = Chem.CombineMols(chain_a, chain_b)

    b1, b1_maps, b1_roles = parse_mapped_source(PATHS["b1_graph"], "B1", b1_legend)
    b29, b29_maps, b29_roles = parse_mapped_source(PATHS["b29_graph"], "B29", b29_legend)
    if set(b1_maps) != set(range(1, 20)):
        raise ValueError("B1 source graph must retain exact atom map IDs 1–19")
    if set(b29_maps) != set(range(1, 93)):
        raise ValueError("B29 source graph must retain exact atom map IDs 1–92")
    b1_map1 = b1.GetAtomWithIdx(b1_maps[1])
    if b1_map1.GetSymbol() != "*" or sorted(atom.GetAtomMapNum() for atom in b1_map1.GetNeighbors()) != [2]:
        raise ValueError("B1 map 1 must be the wildcard protein-boundary placeholder adjacent only to map 2")
    if "PheB1" not in b1_roles[1] or "N" not in b1_roles[1]:
        raise ValueError("B1 map legend does not identify map 1 as the PheB1 N placeholder")
    if b1.GetAtomWithIdx(b1_maps[2]).GetSymbol() != "C":
        raise ValueError("B1 map 2 must be the source carbonyl carbon")
    b29_map92 = b29.GetAtomWithIdx(b29_maps[92])
    if b29_map92.GetSymbol() != "N" or sorted(atom.GetAtomMapNum() for atom in b29_map92.GetNeighbors()) != [84]:
        raise ValueError("B29 map 92 must be the Lys NZ cap adjacent only to map 84")
    if "LysB29 NZ" not in b29_roles[92]:
        raise ValueError("B29 map legend does not identify map 92 as the LysB29 NZ cap")
    if b29.GetAtomWithIdx(b29_maps[84]).GetSymbol() != "C":
        raise ValueError("B29 map 84 must be the source carbamate carbonyl carbon")

    b1_no_cap = remove_source_cap(b1, 1, 2, "B1")
    b29_no_cap = remove_source_cap(b29, 92, 84, "B29")
    b1_protein_n = find_template_atom(combined, "B", 1, "N")
    b29_protein_nz = find_template_atom(combined, "B", 29, "NZ")
    combined = attach_fragment(combined, b1_no_cap, 2, b1_protein_n, "B1")
    combined = attach_fragment(combined, b29_no_cap, 84, b29_protein_nz, "B29")

    rw_mol = Chem.RWMol(combined)
    disulfide_records: list[dict[str, Any]] = []
    for left, right in DISULFIDES:
        left_index = find_template_atom(rw_mol, *left)
        right_index = find_template_atom(rw_mol, *right)
        left_atom = rw_mol.GetAtomWithIdx(left_index)
        right_atom = rw_mol.GetAtomWithIdx(right_index)
        if left_atom.GetSymbol() != "S" or right_atom.GetSymbol() != "S":
            raise ValueError(f"Disulfide endpoints are not sulfur atoms: {left}, {right}")
        if rw_mol.GetBondBetweenAtoms(left_index, right_index) is not None:
            raise ValueError(f"Duplicate disulfide edge: {left}, {right}")
        rw_mol.AddBond(left_index, right_index, Chem.BondType.SINGLE)
        bond = rw_mol.GetBondBetweenAtoms(left_index, right_index)
        bond.SetProp("_NNCProvenance", "Declared insulin disulfide source: "
                     f"{left[0]}{left[1]} SG–{right[0]}{right[1]} SG")
        disulfide_records.append({
            "atom_identities": [
                {"chain": left[0], "residue_number": left[1], "atom_name": left[2]},
                {"chain": right[0], "residue_number": right[1], "atom_name": right[2]},
            ],
            "bond_order": "single",
            "provenance": "declared insulin disulfide connectivity",
        })
    molecule = rw_mol.GetMol()
    Chem.SanitizeMol(molecule)
    for atom in molecule.GetAtoms():
        atom.UpdatePropertyCache(strict=True)

    join_operations = [
        {
            "component": "B1",
            "removed_cap_map_id": 1,
            "removed_cap_role": b1_roles[1],
            "source_endpoint_map_id": 2,
            "protein_endpoint": {"chain": "B", "residue_number": 1, "residue_name": "PHE", "atom_name": "N"},
            "bond_order": "single",
            "interpretation": "map 1 placeholder was removed and its single edge to map 2 was replaced by the PheB1 N–map 2 amide bond",
        },
        {
            "component": "B29",
            "removed_cap_map_id": 92,
            "removed_cap_role": b29_roles[92],
            "source_endpoint_map_id": 84,
            "protein_endpoint": {"chain": "B", "residue_number": 29, "residue_name": "LYS", "atom_name": "NZ"},
            "bond_order": "single",
            "interpretation": "map 92 NZ cap was removed and its single edge to map 84 was replaced by the LysB29 NZ–map 84 carbamate bond",
        },
    ]
    return molecule, {"join_operations": join_operations, "disulfides": disulfide_records}


def main() -> None:
    formula_audit = load_json(PATHS["formula_audit"])
    b1_audit = load_json(PATHS["b1_topology_audit"])
    b29_audit = load_json(PATHS["b29_topology_audit"])
    reconciliation = load_json(PATHS["reconciliation"])
    boundary_audit = load_json(PATHS["boundary_audit"])
    b1_legend = load_json(PATHS["b1_map_legend"])
    b29_legend = load_json(PATHS["b29_map_legend"])

    formula_hash = validate_hashes(
        formula_audit, b1_audit, b29_audit, reconciliation, boundary_audit
    )
    hash_records, input_hashes = formula_hash
    molecule, assembly = build_full_graph(formula_audit, b1_legend, b29_legend)

    formula_counts, heavy_counts = formula_from_atom_ledger(molecule)
    formula = formula_counts_text(formula_counts)
    rdkit_formula = rdMolDescriptors.CalcMolFormula(molecule)
    if formula_counts != PUBLISHED_COUNTS or formula != PUBLISHED_FORMULA or rdkit_formula != PUBLISHED_FORMULA:
        raise ValueError(
            "Assembled graph formula mismatch: atom ledger="
            f"{formula} {formula_counts}, RDKit={rdkit_formula}, expected={PUBLISHED_FORMULA}"
        )
    if formula_audit.get("published_reference", {}).get("reported_hrms_formula") != PUBLISHED_FORMULA:
        raise ValueError("Formula audit's published reference differs from the required NNC2215 formula")

    atom_count = molecule.GetNumAtoms()
    edge_count = molecule.GetNumBonds()
    component_count = len(Chem.GetMolFrags(molecule, asMols=False))
    cycle_rank = edge_count - atom_count + component_count
    if atom_count != 507 or component_count != 1:
        raise ValueError(f"Unexpected graph size/connectivity: atoms={atom_count}, components={component_count}")
    if cycle_rank != 21:
        raise ValueError(f"Unexpected full-graph cycle rank {cycle_rank}; expected 21")
    total_charge = sum(atom.GetFormalCharge() for atom in molecule.GetAtoms())
    if total_charge != 0:
        raise ValueError(f"Neutral graph has nonzero formal charge {total_charge}")

    mapping_states = []
    b1_states = {item["state"]: item for item in b1_audit.get("candidate_states", [])}
    b29_states = {item["state"]: item for item in b29_audit.get("states", [])}
    boundary_states = {item["state"]: item for item in boundary_audit.get("states", [])}
    reconciliation_models = {item["state"]: item for item in reconciliation.get("models", [])}
    formula_models = {item["state"]: item for item in formula_audit.get("model_results", [])}
    pdb_checks = []
    candidate_hashes_after_read: dict[str, str] = {}

    for state in STATE_ORDER:
        state_paths = STATE_DATA[state]
        candidate_path = PATHS[state_paths["candidate_path_key"]]
        candidate_hash_before = sha256(candidate_path)
        parsed = parse_pdb(candidate_path)
        candidate_hash_after = sha256(candidate_path)
        if candidate_hash_before != candidate_hash_after:
            raise ValueError(f"Candidate PDB changed while being read: {relative(candidate_path)}")
        candidate_hashes_after_read[relative(candidate_path)] = candidate_hash_after

        if parsed["sequences"] != EXPECTED_SEQUENCES:
            raise ValueError(f"Candidate PDB sequence mismatch in {state}: {parsed['sequences']}")
        disulfide_pairs = pdb_disulfide_pairs(parsed["ssbonds"])
        if len(parsed["ssbonds"]) != 3 or disulfide_pairs != expected_disulfide_pairs():
            raise ValueError(f"Candidate PDB disulfide records differ from declared three bridges in {state}")
        link_checks = boundary_link_cross_checks(parsed, boundary_states[state], state)
        if len(link_checks) != 3:
            raise ValueError(f"Candidate PDB did not resolve all three required attachment LINKs in {state}")

        b1_count = b1_states[state].get("exact_mapping_count")
        b29_comparison = b29_states[state].get("complete_product_exact_graph_comparison", {})
        b29_count = b29_comparison.get("exact_mapping_count")
        b29_exact_mappings = b29_comparison.get("exact_mappings", [])
        if b1_count != 1 or b29_count != 16 or len(b29_exact_mappings) != 16:
            raise ValueError(f"Mapping-count records changed in {state}: B1={b1_count}, B29={b29_count}")
        mapping_states.append({
            "state": state,
            "b1_exact_mapping_count": b1_count,
            "b29_exact_mapping_count": b29_count,
            "b29_full_ambiguity_preserved": len(b29_exact_mappings) == 16,
            "coordinate_mappings_used_to_assemble_graph": False,
        })

        reconciliation_model = reconciliation_models[state]
        formula_model = formula_models[state]
        a21_oxt_present = any(
            record["chain"] == "A"
            and record["residue_name"] == "ASN"
            and record["residue_number"] == 21
            and record["element"] == "O"
            for record in parsed["oxt_records"]
        )
        formula_audit_a21_oxt = formula_model.get("model_molecule_atom_inventory", {}).get("protein_atom_records", {}).get("chain_a_asn21_oxt_present")

        pdb_checks.append({
            "state": state,
            "candidate_pdb": relative(candidate_path),
            "candidate_pdb_sha256_before_parse": candidate_hash_before,
            "candidate_pdb_sha256_after_parse": candidate_hash_after,
            "candidate_pdb_unchanged": True,
            "sequence_by_chain": parsed["sequences"],
            "sequence_matches_formula_audit_and_required_sequences": True,
            "ssbond_record_count": len(parsed["ssbonds"]),
            "ssbond_pairs": [
                {"first": record["first"], "second": record["second"]}
                for record in parsed["ssbonds"]
            ],
            "ssbond_context_matches_declared_three_disulfides": True,
            "attachment_link_checks": link_checks,
            "a21_asn_oxt_coordinate_record_present": a21_oxt_present,
            "source_model_formula_audit_a21_oxt_record_present": formula_audit_a21_oxt,
            "a21_oxt_context": (
                "corrected open candidate includes ASN A21 OXT although the original model-formula audit records it absent"
                if state == "active_open" and a21_oxt_present and formula_audit_a21_oxt is False
                else "open candidate retains the source-model omission"
                if state == "active_open" and not a21_oxt_present
                else "closed candidate has the terminal OXT record"
                if state == "inactive_closed" and a21_oxt_present
                else "observed terminal-coordinate record; not part of graph assembly"
            ),
            "coordinates_used_for_graph_or_atom_mapping": False,
        })

    # Verify that the source PDB inputs used by the existing formula audit and
    # coordinate reconciliation are intact, without parsing them for graph data.
    source_file_hash_checks = []
    for state in STATE_ORDER:
        source_path = PATHS[STATE_DATA[state]["source_path_key"]]
        actual = sha256(source_path)
        source_file_hash_checks.append({
            "state": state,
            "path": relative(source_path),
            "observed_sha256": actual,
            "formula_audit_sha256": formula_models[state].get("sha256"),
            "reconciliation_sha256": reconciliation_models[state].get("source_sha256"),
            "all_expected_hashes_match": actual == formula_models[state].get("sha256") == reconciliation_models[state].get("source_sha256"),
        })
        if not source_file_hash_checks[-1]["all_expected_hashes_match"]:
            raise ValueError(f"Source model PDB hash mismatch in {state}")

    nodes = [atom_record(atom) for atom in molecule.GetAtoms()]
    bonds = [bond_record(bond) for bond in molecule.GetBonds()]
    if len(nodes) != atom_count or len(bonds) != edge_count:
        raise ValueError("Serialized atom/bond ledger does not match the assembled RDKit graph")
    if any("provenance" not in node for node in nodes):
        raise ValueError("At least one graph atom lacks provenance")
    if sum(node["hydrogen_count_not_as_nodes"] for node in nodes) != PUBLISHED_COUNTS["H"]:
        raise ValueError("Per-atom hydrogen ledger does not reconcile with the published formula")

    report = {
        "schema": "dia-health-one.nnc2215-full-conjugate-graph.v1",
        "run": {
            "script": relative(Path(__file__)),
            "script_sha256": sha256(Path(__file__)),
            "python_version": platform.python_version(),
            "rdkit_version": rdBase.rdkitVersion,
            "coordinate_files_modified": False,
            "coordinates_used_for_graph_assembly_or_mapping": False,
        },
        "representation": {
            "type": "source-derived neutral constitutional graph",
            "fasta_template_convention": "Pinned RDKit MolFromFASTA neutral amino-acid residue templates, with N-terminal NH2 and C-terminal carboxylic-acid OXT/OH termini before conjugate attachment.",
            "hydrogen_representation": "Hydrogens are not graph nodes; per-atom RDKit implicit and bracket-hydrogen counts are recorded.",
            "conditional_tricarboxylate_audit": {
                "path": relative(PATHS["conditional_carboxylate_audit"]),
                "representation_is_separate": True,
                "merged_into_this_neutral_graph": False,
                "interpretation": "A conditional pH 7.4 local B29 tri-carboxylate representation; it does not assign whole-system protonation or total charge.",
            },
        },
        "sequences": {
            "source_audit": relative(PATHS["formula_audit"]),
            "chain_A": EXPECTED_SEQUENCES["A"],
            "chain_B_desB30": EXPECTED_SEQUENCES["B"],
            "sequence_audit_match": True,
        },
        "input_hashes": hash_records,
        "source_model_pdb_hash_checks": source_file_hash_checks,
        "assembly": {
            "source_components": [
                {"component": "insulin chain A", "method": "RDKit MolFromFASTA", "sequence": EXPECTED_SEQUENCES["A"]},
                {"component": "insulin chain B desB30", "method": "RDKit MolFromFASTA", "sequence": EXPECTED_SEQUENCES["B"]},
                {"component": "B1 glucoside", "source_graph": relative(PATHS["b1_graph"]), "map_legend": relative(PATHS["b1_map_legend"]), "source_maps_preserved": True},
                {"component": "B29 linker-roof", "source_graph": relative(PATHS["b29_graph"]), "map_legend": relative(PATHS["b29_map_legend"]), "source_maps_preserved": True},
            ],
            "join_operations": assembly["join_operations"],
            "disulfides": assembly["disulfides"],
            "only_declared_disulfides_added": True,
        },
        "formula": {
            "rdkit_formula": rdkit_formula,
            "atom_ledger_formula": formula,
            "element_counts": formula_counts,
            "heavy_element_counts": heavy_counts,
            "reported_formula": PUBLISHED_FORMULA,
            "reported_formula_element_counts": PUBLISHED_COUNTS,
            "formula_matches_published_reference": True,
            "interpretation": "Internal source-graph consistency with the formula reported in the Nature supporting information; not independent experimental confirmation.",
            "total_formal_charge": total_charge,
        },
        "graph": {
            "heavy_atom_count": atom_count,
            "explicit_hydrogen_nodes": 0,
            "implicit_and_bracket_hydrogens_ledger_total": formula_counts["H"],
            "edge_count": edge_count,
            "connected_component_count": component_count,
            "cycle_rank_E_minus_V_plus_components": cycle_rank,
            "rdkit_sanitization": "pass",
            "valence_status": "pass; RDKit sanitization and strict atom property-cache update completed",
        },
        "mapping_ambiguity": {
            "states": mapping_states,
            "b29_coordinate_map_selected": False,
            "coordinate_values_used_to_select_or_rank_mappings": False,
        },
        "pdb_cross_checks": pdb_checks,
        "atoms": nodes,
        "bonds": bonds,
        "checks": {
            "source_sequences_match_expected": True,
            "one_connected_full_graph": component_count == 1,
            "neutral_formula_matches_reported_formula": formula == PUBLISHED_FORMULA,
            "per_atom_provenance_complete": all("provenance" in node for node in nodes),
            "implicit_hydrogen_ledger_matches_formula": sum(node["hydrogen_count_not_as_nodes"] for node in nodes) == formula_counts["H"],
            "three_declared_disulfides_present": len(assembly["disulfides"]) == 3,
            "b1_and_b29_boundary_caps_removed_and_replaced": len(assembly["join_operations"]) == 2,
            "candidate_pdb_sequences_links_and_disulfides_match": all(
                item["sequence_matches_formula_audit_and_required_sequences"]
                and item["ssbond_context_matches_declared_three_disulfides"]
                and len(item["attachment_link_checks"]) == 3
                for item in pdb_checks
            ),
            "b1_one_and_b29_sixteen_exact_mappings_per_state_retained": all(
                item["b1_exact_mapping_count"] == 1 and item["b29_exact_mapping_count"] == 16
                for item in mapping_states
            ),
            "source_and_corrected_candidate_hashes_match_upstream_records": all(
                item["all_recorded_expected_hashes_match"] for key, item in hash_records.items()
                if key != "audit_record_files"
            ) and all(item["all_expected_hashes_match"] for item in source_file_hash_checks),
            "conditional_tri_carboxylate_remains_separate": True,
            "full_system_protonation_force_field_or_md_resolved": False,
        },
        "assessment": {
            "status": "source-derived neutral constitutional graph assembled and internally audited",
            "formula_and_connectivity": "The graph is connected, valence-sanitized, and formula-consistent with the reported C325H461N79O97S6 composition.",
            "scope": "The corrected model-built PDBs were used only for named sequence, SSBOND, LINK, and hash cross-checks.",
            "not_established": [
                "a unique B29 coordinate atom map; sixteen exact mappings per state remain",
                "PheB1 hydrogen ownership; the H2/H2A discrepancy remains open",
                "neutral macrocycle oxygen-label assignment or a physiological-pH full-system microstate",
                "histidine tautomer choice beyond the pinned RDKit neutral template",
                "roof stereochemistry beyond the source graph",
                "whole-system charge, force-field parameters, structural relaxation, or molecular dynamics",
                "drug efficacy or clinical readiness",
            ],
        },
    }
    OUTPUT_PATH.write_text(json.dumps(report, indent=2, sort_keys=False) + "\n", encoding="utf-8")
    print(f"Wrote {relative(OUTPUT_PATH)}")
    print(f"Formula: {rdkit_formula} (matches reported formula)")
    print(f"Graph: {atom_count} heavy atoms, {edge_count} bonds, {component_count} component, cycle rank {cycle_rank}")
    print("Mapping ambiguity retained: B1 1/state; B29 16/state")
    print("PDB sequence, 3 disulfides, and 3 attachment LINK checks: pass in both states")
    print("Coordinates used for graph assembly/mapping: no")


if __name__ == "__main__":
    main()
