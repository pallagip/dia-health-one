#!/usr/bin/env python3
"""Audit the source-mapped B1 glucoside graph against corrected NNC2215 models.

PDB LINK and CONECT records are treated as adjacency only. They do not supply
bond orders, formal charges, or stereochemistry. Coordinates are measured only
after graph mapping, as a geometry screen; they never select a mapping.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import platform
from typing import Any

import networkx as nx
from rdkit import Chem, rdBase
from rdkit.Chem import rdMolDescriptors


ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
GRAPH_PATH = ROOT / "structures" / "nnc2215-b1-glucoside-attachment.smi"
LEGEND_PATH = DATA / "nnc2215-b1-glucoside-atom-map.json"
RECONCILIATION_PATH = DATA / "nnc2215-carbon-reconciliation-candidate.json"
MANIFEST_PATH = DATA / "structure-manifest.json"
OUTPUT_PATH = DATA / "nnc2215-b1-topology-audit.json"

EXPECTED_MAP_IDS = set(range(1, 20))
EXPECTED_GLUCOSE_CIP = {"C1": "R", "C2": "R", "C3": "S", "C4": "S", "C5": "R"}
EXPECTED_STEREO_MAPS = {"C1": 9, "C2": 11, "C3": 13, "C4": 15, "C5": 17}
EXPECTED_SUBSTITUENT_HEAVY = Counter({"C": 10, "O": 8})
EXPECTED_TARGET_HEAVY = Counter({"C": 10, "O": 8, "N": 1})
EXPECTED_SELECTED_NODE_COUNT = 19
TARGET_HEAVY_ELEMENTS = {"C", "O", "N", "P", "S", "F"}
STATE_ORDER = ("active_open", "inactive_closed")

# Broad, heuristic screening intervals, adapted from the B29 audit. These are
# not experimental bond-length ranges and do not affect graph mapping.
COVALENT_RADII_ANGSTROM = {"C": 0.76, "N": 0.71, "O": 0.66, "P": 1.07, "S": 1.05, "F": 0.57}
BOND_CLASS_OFFSETS_ANGSTROM = {
    "single": [-0.18, 0.25],
    "double": [-0.35, -0.06],
    "triple": [-0.45, -0.18],
    "aromatic": [-0.20, 0.05],
}

AA_FORMULAE: dict[str, Counter[str]] = {
    "A": Counter(C=3, H=7, N=1, O=2), "R": Counter(C=6, H=14, N=4, O=2),
    "N": Counter(C=4, H=8, N=2, O=3), "D": Counter(C=4, H=7, N=1, O=4),
    "C": Counter(C=3, H=7, N=1, O=2, S=1), "Q": Counter(C=5, H=10, N=2, O=3),
    "E": Counter(C=5, H=9, N=1, O=4), "G": Counter(C=2, H=5, N=1, O=2),
    "H": Counter(C=6, H=9, N=3, O=2), "I": Counter(C=6, H=13, N=1, O=2),
    "L": Counter(C=6, H=13, N=1, O=2), "K": Counter(C=6, H=14, N=2, O=2),
    "M": Counter(C=5, H=11, N=1, O=2, S=1), "F": Counter(C=9, H=11, N=1, O=2),
    "P": Counter(C=5, H=9, N=1, O=2), "S": Counter(C=3, H=7, N=1, O=3),
    "T": Counter(C=4, H=9, N=1, O=3), "W": Counter(C=11, H=12, N=2, O=2),
    "Y": Counter(C=9, H=11, N=1, O=3), "V": Counter(C=5, H=11, N=1, O=2),
}
B1_B22_SEQUENCE = "FVNQHLCGSHLVEALYLVCGER"
MONOISOTOPIC_MASS_DA = {
    "C": 12.0,
    "H": 1.00782503223,
    "N": 14.00307400443,
    "O": 15.99491461957,
    "S": 31.9720711744,
}
PROTON_MASS_DA = 1.007276466621


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def formula_text(counts: Counter[str] | dict[str, int]) -> str:
    order = ("C", "H", "N", "O", "S", "P", "F", "CL", "BR", "I")
    keys = [element for element in order if counts.get(element, 0)]
    keys.extend(sorted(element for element in set(counts) - set(keys) if counts.get(element, 0)))
    return "".join(element + (str(counts[element]) if counts[element] != 1 else "") for element in keys)


def normalized_counts(counts: Counter[str] | dict[str, int]) -> dict[str, int]:
    return {key: int(value) for key, value in sorted(counts.items()) if value}


def bond_label(bond: Chem.Bond) -> tuple[str, float]:
    if bond.GetIsAromatic():
        return "aromatic", 1.5
    order = bond.GetBondTypeAsDouble()
    for value, label in ((1.0, "single"), (2.0, "double"), (3.0, "triple")):
        if math.isclose(order, value):
            return label, value
    raise ValueError(f"Unsupported source bond order: {bond.GetBondType()}")


def graph_summary(graph: nx.Graph) -> dict[str, Any]:
    components = nx.number_connected_components(graph) if graph.number_of_nodes() else 0
    return {
        "node_count": graph.number_of_nodes(),
        "element_counts": dict(sorted(Counter(d["element"] for _, d in graph.nodes(data=True)).items())),
        "edge_count": graph.number_of_edges(),
        "connected_component_count": components,
        "cycle_rank_E_minus_V_plus_components": graph.number_of_edges() - graph.number_of_nodes() + components,
    }


def parse_source_graph() -> tuple[Chem.Mol, nx.Graph, dict[str, Any]]:
    legend = load_json(LEGEND_PATH)
    smiles_tokens = GRAPH_PATH.read_text(encoding="utf-8").split()
    if not smiles_tokens:
        raise ValueError("Mapped B1 SMILES file is empty")
    molecule = Chem.MolFromSmiles(smiles_tokens[0])
    if molecule is None:
        raise ValueError("RDKit failed to parse mapped B1 SMILES")
    Chem.SanitizeMol(molecule)

    atoms_by_map: dict[int, Chem.Atom] = {}
    for atom in molecule.GetAtoms():
        map_id = atom.GetAtomMapNum()
        if map_id <= 0 or map_id in atoms_by_map:
            raise ValueError(f"Atom map IDs must be unique positive integers; duplicate/invalid {map_id}")
        atoms_by_map[map_id] = atom
    if set(atoms_by_map) != EXPECTED_MAP_IDS:
        raise ValueError(f"Mapped graph map IDs differ from expected 1..19: {sorted(atoms_by_map)}")

    map_legend = legend.get("map_id_legend")
    expected_keys = {str(map_id) for map_id in EXPECTED_MAP_IDS}
    if not isinstance(map_legend, dict) or set(map_legend) != expected_keys:
        raise ValueError("Map legend must contain exactly string keys 1 through 19")
    for map_id, atom in atoms_by_map.items():
        entry = map_legend[str(map_id)]
        if map_id == 1:
            if atom.GetAtomicNum() != 0 or entry.get("element_in_smiles") != "*":
                raise ValueError("Map 1 must be the mapped dummy attachment placeholder")
            identity = entry.get("model_match_identity", "")
            if entry.get("model_match_element") != "N" or "PheB1" not in identity or "backbone N" not in identity:
                raise ValueError("Map 1 legend must map the placeholder to the PheB1 backbone N")
        elif entry.get("element") != atom.GetSymbol():
            raise ValueError(f"Map {map_id} element {atom.GetSymbol()} conflicts with legend {entry.get('element')}")
    if map_legend["1"].get("role", "").count("PheB1 N") != 1:
        raise ValueError("Map 1 must be the single PheB1-N boundary placeholder")

    # Require the legend's listed connectivity to be a complete exact account
    # of the mapped graph (bond order is source metadata, never PDB metadata).
    listed_bonds: dict[tuple[int, int], str] = {}
    for item in legend.get("key_connectivity", []):
        pair = tuple(sorted(int(value) for value in item["atoms"]))
        if pair in listed_bonds:
            raise ValueError(f"Duplicate edge in map legend: {pair}")
        listed_bonds[pair] = item["bond"]
    observed_bonds: dict[tuple[int, int], str] = {}
    for bond in molecule.GetBonds():
        left = molecule.GetAtomWithIdx(bond.GetBeginAtomIdx()).GetAtomMapNum()
        right = molecule.GetAtomWithIdx(bond.GetEndAtomIdx()).GetAtomMapNum()
        observed_bonds[tuple(sorted((left, right)))] = bond_label(bond)[0]
    if listed_bonds != observed_bonds:
        raise ValueError("The map legend connectivity/bond classes do not exactly match the mapped SMILES")

    graph = nx.Graph()
    for map_id, atom in atoms_by_map.items():
        element = map_legend[str(map_id)].get("model_match_element") if map_id == 1 else atom.GetSymbol()
        graph.add_node(map_id, element=element, source_smiles_element=atom.GetSymbol())
    for bond in molecule.GetBonds():
        left = molecule.GetAtomWithIdx(bond.GetBeginAtomIdx()).GetAtomMapNum()
        right = molecule.GetAtomWithIdx(bond.GetEndAtomIdx()).GetAtomMapNum()
        label, order = bond_label(bond)
        graph.add_edge(left, right, bond_class=label, bond_order=order)

    transcription = legend["transcription"]
    stereo = transcription["stereochemistry"]
    # The source map's recorded CIP assignments were generated by RDKit's
    # legacy CIP labeler. Keep that convention for the source-legend check,
    # and separately record the newer CIP engine's labels for transparency.
    centers = Chem.FindMolChiralCenters(
        molecule, includeUnassigned=True, includeCIP=True, useLegacyImplementation=True
    )
    observed_cip: dict[str, str] = {}
    for atom_index, cip in centers:
        atom = molecule.GetAtomWithIdx(atom_index)
        if cip == "?":
            raise ValueError(f"Unassigned glucose stereochemistry at map {atom.GetAtomMapNum()}")
        observed_cip[str(atom.GetAtomMapNum())] = cip
    modern_centers = Chem.FindMolChiralCenters(
        molecule, includeUnassigned=True, includeCIP=True, useLegacyImplementation=False
    )
    modern_cip = {
        str(molecule.GetAtomWithIdx(atom_index).GetAtomMapNum()): cip
        for atom_index, cip in modern_centers
    }
    map_to_glucose = {str(map_id): source for source, map_id in stereo["stereocenter_atom_maps"].items()}
    source_cip_legend = stereo["rdkit_cip_by_source_glucose_carbon"]
    if {key: source_cip_legend.get(key) for key in EXPECTED_GLUCOSE_CIP} != EXPECTED_GLUCOSE_CIP or source_cip_legend.get("C6") != "not stereogenic (CH2OH)":
        raise ValueError("Source map's glucose CIP legend differs from expected C1-C5 R/R/S/S/R")
    if stereo["stereocenter_atom_maps"] != EXPECTED_STEREO_MAPS:
        raise ValueError("Source map's glucose stereocenter map legend differs from expected maps")
    expected_by_map = {str(EXPECTED_STEREO_MAPS[label]): cip for label, cip in EXPECTED_GLUCOSE_CIP.items()}
    if observed_cip != expected_by_map or map_to_glucose.keys() != expected_by_map.keys():
        raise ValueError(f"RDKit glucose CIP assignments {observed_cip} conflict with source legend {expected_by_map}")

    element_inventory = Counter(atom.GetSymbol() for atom in molecule.GetAtoms() if atom.GetAtomicNum() > 1)
    element_inventory.subtract(Counter({"*": 1}))
    if +element_inventory != EXPECTED_SUBSTITUENT_HEAVY:
        raise ValueError(f"Substituent heavy inventory excluding map-1 placeholder is {+element_inventory}")

    acid_mol = acid_from_mapped_amide(molecule)
    acid_formula = rdMolDescriptors.CalcMolFormula(acid_mol)
    if acid_formula != "C10H18O9":
        raise ValueError(f"Derived deprotected acid formula is {acid_formula}, expected C10H18O9")
    increment_counts = Counter(C=10, H=18, O=9)
    increment_counts.subtract(Counter(H=2, O=1))
    increment_formula = formula_text(increment_counts)
    if increment_formula != "C10H16O8":
        raise ValueError(f"Amide condensation increment is {increment_formula}, expected C10H16O8")

    return molecule, graph, {
        "smiles": smiles_tokens[0],
        "mapped_atom_count_including_placeholder": len(atoms_by_map),
        "map_ids_unique_positive_exact_1_to_19": True,
        "map_legend_exact_and_element_consistent": True,
        "rdkit_sanitization": "pass",
        "map_1_placeholder": {
            "map_id": 1,
            "smiles_element": "*",
            "target_match_element": map_legend["1"]["model_match_element"],
            "counts_toward_substituent_formula_or_inventory": False,
        },
        "source_substituent_heavy_inventory_excluding_placeholder": normalized_counts(EXPECTED_SUBSTITUENT_HEAVY),
        "rdkit_mapped_graph_formula_including_placeholder": rdMolDescriptors.CalcMolFormula(molecule),
        "deprotected_acid_formula": acid_formula,
        "net_amide_condensation_formula_increment": increment_formula,
        "formula_derivation": {
            "acid_formula": "C10H18O9",
            "condensation_loss": "H2O",
            "net_increment": "C10H16O8",
            "acid_graph_construction": "Remove mapped dummy map 1, replace its bond to carbonyl C map 2 by a terminal hydroxyl oxygen, sanitize with RDKit, and count implicit hydrogens.",
            "placeholder_formula_contribution": "zero; the PheB1 nitrogen belongs to insulin, not the substituent increment",
        },
        "glucose_stereochemistry": {
            "rdkit_cip_by_source_carbon": {
                label: observed_cip[str(map_id)] for label, map_id in EXPECTED_STEREO_MAPS.items()
            },
            "rdkit_cip_algorithm": "RDKit legacy CIP labeler (`useLegacyImplementation=True`), matching the source map artifact's recorded convention",
            "rdkit_new_cip_engine_observation_by_source_carbon": {
                label: modern_cip[str(map_id)] for label, map_id in EXPECTED_STEREO_MAPS.items()
            },
            "cip_method_difference_note": "The newer RDKit CIP engine reports different R/S labels for some centers on this same sanitized stereochemical graph. Both outputs are recorded; the source-map legend comparison uses its legacy RDKit convention. No stereochemical atom tags are changed.",
            "source_legend_cip_by_source_carbon": dict(stereo["rdkit_cip_by_source_glucose_carbon"]),
            "stereocenter_atom_maps": dict(stereo["stereocenter_atom_maps"]),
            "center_count": len(centers),
            "all_five_expected_centers_assigned": True,
            "C6": "not stereogenic (CH2OH)",
        },
        "graph_summary": graph_summary(graph),
        "map_legend_atom_roles": {
            str(map_id): {
                "source_label": map_legend[str(map_id)].get("source_label"),
                "role": map_legend[str(map_id)].get("role"),
                "element": graph.nodes[map_id]["element"],
            }
            for map_id in sorted(graph.nodes)
        },
    }


def acid_from_mapped_amide(molecule: Chem.Mol) -> Chem.Mol:
    rw_mol = Chem.RWMol(molecule)
    placeholder = next(atom for atom in rw_mol.GetAtoms() if atom.GetAtomMapNum() == 1)
    carbonyl = next(atom for atom in rw_mol.GetAtoms() if atom.GetAtomMapNum() == 2)
    carbonyl_index = carbonyl.GetIdx()
    rw_mol.RemoveAtom(placeholder.GetIdx())
    # Removing index 0 shifts later indices, so find mapped carbonyl again.
    carbonyl_index = next(atom.GetIdx() for atom in rw_mol.GetAtoms() if atom.GetAtomMapNum() == 2)
    hydroxyl = Chem.Atom("O")
    hydroxyl_index = rw_mol.AddAtom(hydroxyl)
    rw_mol.AddBond(carbonyl_index, hydroxyl_index, Chem.BondType.SINGLE)
    acid = rw_mol.GetMol()
    Chem.SanitizeMol(acid)
    return acid


def derive_fragment_mass() -> dict[str, Any]:
    residue_sum: Counter[str] = Counter()
    for residue in B1_B22_SEQUENCE:
        residue_sum.update(AA_FORMULAE[residue])
    peptide = Counter(residue_sum)
    peptide.subtract(Counter(H=2 * (len(B1_B22_SEQUENCE) - 1), O=len(B1_B22_SEQUENCE) - 1))
    modified = Counter(peptide)
    modified.update(Counter(C=10, H=16, O=8))
    modified += Counter()
    neutral_mass = sum(MONOISOTOPIC_MASS_DA[element] * count for element, count in modified.items())
    mz = (neutral_mass + 2 * PROTON_MASS_DA) / 2
    calculated_source = 1376.2
    return {
        "sequence": B1_B22_SEQUENCE,
        "residue_count": len(B1_B22_SEQUENCE),
        "digestion_and_redox_assumptions": [
            "Trypsin cleavage after ArgB22 gives a free B1 N terminus and a free B22 carboxyl terminus.",
            "Peptide formula is the sum of neutral amino-acid formulae minus H2O for each of 21 peptide bonds.",
            "TCEP-reduced cysteines are modeled as separate thiols; no disulfide subtraction is applied to this isolated fragment.",
            "The B1 glycoside adds net C10H16O8 after amide condensation.",
        ],
        "neutral_unmodified_fragment_formula": formula_text(peptide),
        "net_attachment_increment": "C10H16O8",
        "neutral_modified_fragment_formula": formula_text(modified),
        "neutral_monoisotopic_mass_da": neutral_mass,
        "charge_state": "[M+2H]2+",
        "monoisotopic_mass_convention": {
            "atomic_masses_da": MONOISOTOPIC_MASS_DA,
            "proton_mass_da": PROTON_MASS_DA,
            "calculation": "(neutral monoisotopic mass + 2 × proton mass) / 2; proton mass is used directly, without an electron-mass correction",
        },
        "calculated_m_z": mz,
        "calculated_m_z_rounded_one_decimal": round(mz, 1),
        "source_calculated_m_z": calculated_source,
        "source_observed_m_z": 1376.6,
        "source_calculated_rounding_match": round(mz, 1) == calculated_source,
        "difference_from_source_calculated_m_z": mz - calculated_source,
        "benchmark_source": "Nature Supporting Information §4, p. 11",
    }


def parse_link_endpoint(line: str, first: bool) -> dict[str, str]:
    fields = (
        (line[12:16], line[16:17], line[17:20], line[21:22], line[22:26], line[26:27])
        if first else
        (line[42:46], line[46:47], line[47:50], line[51:52], line[52:56], line[56:57])
    )
    atom_name, altloc, residue, chain, residue_number, insertion_code = (value.strip() for value in fields)
    return {
        "name": atom_name,
        "altloc": altloc,
        "residue": residue,
        "chain": chain,
        "residue_number": residue_number,
        "insertion_code": insertion_code,
    }


def atom_identity(atom: dict[str, Any]) -> tuple[str, ...]:
    return tuple(atom[key] for key in ("chain", "residue_number", "insertion_code", "residue", "altloc", "name"))


def endpoint_identity(endpoint: dict[str, str]) -> tuple[str, ...]:
    return tuple(endpoint[key] for key in ("chain", "residue_number", "insertion_code", "residue", "altloc", "name"))


def atom_record(atom: dict[str, Any], include_coordinates: bool = True) -> dict[str, Any]:
    result = {key: atom[key] for key in (
        "serial", "record_type", "chain", "residue_number", "insertion_code", "residue", "altloc", "name", "element"
    )}
    if include_coordinates:
        result["coordinates_angstrom"] = atom["coordinates_angstrom"]
    return result


def distance_angstrom(left: list[float], right: list[float]) -> float:
    return math.sqrt(sum((left[index] - right[index]) ** 2 for index in range(3)))


def parse_pdb(path: Path) -> dict[str, Any]:
    atoms: dict[int, dict[str, Any]] = {}
    raw_conect: list[tuple[int, int]] = []
    raw_links: list[tuple[str, dict[str, str], dict[str, str]]] = []
    for line_number, line in enumerate(path.read_text(encoding="ascii").splitlines(), start=1):
        if line.startswith(("ATOM  ", "HETATM")):
            try:
                serial = int(line[6:11])
                coordinates = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
            except ValueError as error:
                raise ValueError(f"Malformed atom record at {path}:{line_number}: {error}") from error
            if serial in atoms:
                raise ValueError(f"Duplicate PDB serial {serial} at {path}:{line_number}")
            element = line[76:78].strip().upper()
            if not element:
                raise ValueError(f"Missing element at {path}:{line_number}")
            atoms[serial] = {
                "serial": serial,
                "record_type": line[:6].strip(),
                "name": line[12:16].strip(),
                "altloc": line[16:17].strip(),
                "residue": line[17:20].strip(),
                "chain": line[21:22].strip(),
                "residue_number": line[22:26].strip(),
                "insertion_code": line[26:27].strip(),
                "element": element,
                "coordinates_angstrom": coordinates,
            }
        elif line.startswith("CONECT"):
            serials = [
                int(line[index:index + 5])
                for index in range(6, len(line), 5)
                if line[index:index + 5].strip()
            ]
            if serials:
                raw_conect.extend((serials[0], neighbor) for neighbor in serials[1:])
        elif line.startswith("LINK  "):
            raw_links.append((line, parse_link_endpoint(line, True), parse_link_endpoint(line, False)))

    for left, right in raw_conect:
        if left not in atoms or right not in atoms:
            raise ValueError(f"Dangling CONECT reference {left}-{right} in {path}")
    conect_pairs = {tuple(sorted((left, right))) for left, right in raw_conect if left != right}
    edge_sources: dict[tuple[int, int], set[str]] = defaultdict(set)
    for pair in conect_pairs:
        edge_sources[pair].add("CONECT")

    identity_to_serials: dict[tuple[str, ...], list[int]] = defaultdict(list)
    for serial, atom in atoms.items():
        identity_to_serials[atom_identity(atom)].append(serial)
    resolved_links = []
    for raw_line, first, second in raw_links:
        resolved = []
        for endpoint in (first, second):
            matches = identity_to_serials.get(endpoint_identity(endpoint), [])
            if len(matches) != 1:
                raise ValueError(f"Unresolved LINK endpoint {endpoint}; matches={matches} in {path}")
            resolved.append(matches[0])
        if resolved[0] == resolved[1]:
            raise ValueError(f"LINK resolves to the same atom twice in {path}: {raw_line}")
        pair = tuple(sorted(resolved))
        edge_sources[pair].add("LINK")
        resolved_links.append({
            "record": raw_line,
            "first_endpoint": first,
            "second_endpoint": second,
            "resolved_serials_in_record_order": resolved,
            "edge_pair": list(pair),
        })
    return {
        "atoms": atoms,
        "raw_conect_reference_count": len(raw_conect),
        "unique_undirected_conect_edge_count": len(conect_pairs),
        "edge_sources": edge_sources,
        "resolved_links": resolved_links,
    }


def selected_target(parsed: dict[str, Any]) -> tuple[nx.Graph, dict[int, dict[str, Any]], dict[str, Any]]:
    atoms = parsed["atoms"]
    phx = {
        serial: atom for serial, atom in atoms.items()
        if atom["chain"] == "B" and atom["residue_number"].strip() == "0"
        and atom["insertion_code"] == "" and atom["residue"] == "PHX"
        and atom["element"] in TARGET_HEAVY_ELEMENTS
    }
    phe_ns = [
        (serial, atom) for serial, atom in atoms.items()
        if atom["chain"] == "B" and atom["residue_number"].strip() == "1"
        and atom["insertion_code"] == "" and atom["residue"] == "PHE"
        and atom["name"] == "N" and atom["element"] == "N"
    ]
    if len(phe_ns) != 1:
        raise ValueError(f"Expected one PHE B1 backbone N, found {len(phe_ns)}")
    selected = dict(phx)
    n_serial, n_atom = phe_ns[0]
    selected[n_serial] = n_atom
    identities = [atom_identity(atom) for atom in selected.values()]
    if len(set(identities)) != len(identities):
        raise ValueError("Duplicate selected atom identities in PHX B0 + PHE B1 N target")
    if len(selected) != EXPECTED_SELECTED_NODE_COUNT:
        raise ValueError(f"Unexpected selected-node count {len(selected)}; expected {EXPECTED_SELECTED_NODE_COUNT}")
    inventory = Counter(atom["element"] for atom in selected.values())
    if inventory != EXPECTED_TARGET_HEAVY:
        raise ValueError(f"Unexpected selected element inventory {dict(inventory)}")

    graph = nx.Graph()
    for serial, atom in selected.items():
        graph.add_node(serial, element=atom["element"], atom=atom)
    for pair, sources in parsed["edge_sources"].items():
        left, right = pair
        if left in selected and right in selected:
            graph.add_edge(left, right, adjacency_sources=sorted(sources))

    n_identity = atom_identity(n_atom)
    c2_matches = [
        (serial, atom) for serial, atom in phx.items() if atom["name"] == "C2" and atom["element"] == "C"
    ]
    if len(c2_matches) != 1:
        raise ValueError(f"Expected one PHX B0 C2 attachment carbon, found {len(c2_matches)}")
    c2_serial, c2_atom = c2_matches[0]
    matching_links = []
    for link in parsed["resolved_links"]:
        if set(link["resolved_serials_in_record_order"]) == {n_serial, c2_serial}:
            matching_links.append(link)
    target_links = [
        link for link in parsed["resolved_links"]
        if any(serial in selected for serial in link["resolved_serials_in_record_order"])
        and any(serial in phx for serial in link["resolved_serials_in_record_order"])
    ]
    if len(matching_links) != 1 or len(target_links) != 1:
        raise ValueError(f"Expected a unique PHE B1 N–PHX B0 C2 LINK; matching={len(matching_links)}, target_links={len(target_links)}")
    if not graph.has_edge(n_serial, c2_serial) or "LINK" not in graph.edges[n_serial, c2_serial]["adjacency_sources"]:
        raise ValueError("The selected PHE B1 N–PHX B0 C2 LINK is not present as target adjacency")

    components = nx.number_connected_components(graph)
    summary = graph_summary(graph)
    summary["element_counts_by_residue"] = {
        "PHX B 0": dict(sorted(Counter(atom["element"] for atom in phx.values()).items())),
        "PHE B 1 N": {"N": 1},
    }
    summary["selected_atoms"] = [atom_record(atom) for _, atom in sorted(selected.items())]
    summary["selected_adjacency_edges"] = [
        {
            "serials": list(pair),
            "atom_names": [selected[pair[0]]["name"], selected[pair[1]]["name"]],
            "adjacency_sources": sorted(sources),
        }
        for pair, sources in sorted(parsed["edge_sources"].items())
        if pair[0] in selected and pair[1] in selected
    ]
    summary["unique_n_to_c2_link"] = {
        "count": 1,
        "record": matching_links[0]["record"],
        "resolved_serials_in_record_order": matching_links[0]["resolved_serials_in_record_order"],
        "phe_b1_n": atom_record(n_atom),
        "phx_b0_c2": atom_record(c2_atom),
        "interpretation": "LINK is adjacency-only; no bond order or charge is inferred.",
    }
    summary["all_selected_resolved_links"] = target_links
    return graph, selected, {
        "summary": summary,
        "phe_b1_n_serial": n_serial,
        "phx_serials": set(phx),
        "selected_serials": set(selected),
        "selected_c2_serial": c2_serial,
        "phe_b1_n_identity": n_identity,
    }


def hydrogen_observations(parsed: dict[str, Any], target_info: dict[str, Any]) -> dict[str, Any]:
    atoms = parsed["atoms"]
    phx_h = [
        atom for atom in atoms.values()
        if atom["chain"] == "B" and atom["residue_number"].strip() == "0"
        and atom["insertion_code"] == "" and atom["residue"] == "PHX"
        and atom["element"] in {"H", "D"}
    ]
    phe_b1_h = [
        atom for atom in atoms.values()
        if atom["chain"] == "B" and atom["residue_number"].strip() == "1"
        and atom["insertion_code"] == "" and atom["residue"] == "PHE"
        and atom["element"] in {"H", "D"}
    ]
    n_serial = target_info["phe_b1_n_serial"]
    n_atom = atoms[n_serial]
    conect_hydrogen_serials = sorted({
        right if left == n_serial else left
        for pair, sources in parsed["edge_sources"].items()
        if "CONECT" in sources
        for left, right in [pair]
        if left == n_serial or right == n_serial
        if atoms[right if left == n_serial else left]["element"] in {"H", "D"}
    })
    n_h = [atom_record(atoms[serial]) for serial in conect_hydrogen_serials]
    nearest = sorted(
        [
            {
                "record_owner": "PHX B 0",
                "atom": atom_record(atom),
                "distance_to_phe_b1_n_angstrom": distance_angstrom(atom["coordinates_angstrom"], n_atom["coordinates_angstrom"]),
                "explicit_conect_neighbor_of_phe_b1_n": atom["serial"] in conect_hydrogen_serials,
            }
            for atom in phx_h
        ],
        key=lambda item: (item["distance_to_phe_b1_n_angstrom"], item["atom"]["serial"]),
    )
    return {
        "phx_b0_hydrogen_record_count": len(phx_h),
        "phx_b0_hydrogen_records": [atom_record(atom) for atom in sorted(phx_h, key=lambda atom: atom["serial"])],
        "phe_b1_hydrogen_record_count": len(phe_b1_h),
        "phe_b1_hydrogen_records": [atom_record(atom) for atom in sorted(phe_b1_h, key=lambda atom: atom["serial"])],
        "phe_b1_n_conect_hydrogen_neighbors": n_h,
        "nearest_phx_owned_hydrogens_to_phe_b1_n": nearest[:3],
        "ownership_and_connectivity_interpretation": (
            "Hydrogen ownership is the ATOM/HETATM residue record owner; explicit adjacency is reported only from CONECT. "
            "Coordinate proximity does not establish an N-H bond."
        ),
    }


def source_symmetry_groups(graph: nx.Graph) -> tuple[list[dict[int, int]], dict[tuple[int, ...], list[int]]]:
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        graph, graph, node_match=lambda left, right: left["element"] == right["element"]
    )
    automorphisms = [dict(mapping) for mapping in matcher.isomorphisms_iter()]
    orbits: dict[tuple[int, ...], list[int]] = defaultdict(list)
    remaining = set(graph.nodes)
    while remaining:
        seed = min(remaining)
        orbit = {auto[seed] for auto in automorphisms}
        orbits[tuple(sorted(orbit))] = sorted(orbit)
        remaining -= orbit
    return automorphisms, orbits


def target_mappings(source: nx.Graph, target: nx.Graph) -> list[dict[int, int]]:
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        source,
        target,
        node_match=lambda source_node, target_node: source_node["element"] == target_node["element"],
    )
    mappings = [dict(mapping) for mapping in matcher.isomorphisms_iter()]
    mappings.sort(key=lambda mapping: tuple(mapping[map_id] for map_id in sorted(mapping)))
    return mappings


def interval_for_edge(left_element: str, right_element: str, bond_class: str) -> list[float]:
    if left_element not in COVALENT_RADII_ANGSTROM or right_element not in COVALENT_RADII_ANGSTROM:
        raise ValueError(f"No screening radius for source edge {left_element}-{right_element}")
    low, high = BOND_CLASS_OFFSETS_ANGSTROM[bond_class]
    radii = COVALENT_RADII_ANGSTROM[left_element] + COVALENT_RADII_ANGSTROM[right_element]
    return [round(radii + low, 2), round(radii + high, 2)]


def mapping_result(
    mapping: dict[int, int],
    index: int,
    source: nx.Graph,
    target: nx.Graph,
    target_atoms: dict[int, dict[str, Any]],
    automorphisms: list[dict[int, int]],
) -> tuple[dict[str, Any], tuple[int, ...]]:
    map_ids = sorted(source.nodes)
    symmetry_forms = [tuple(mapping[auto[map_id]] for map_id in map_ids) for auto in automorphisms]
    symmetry_key = min(symmetry_forms)
    map_records = {
        str(map_id): {
            "target_atom": atom_record(target_atoms[mapping[map_id]]),
            "source_element": source.nodes[map_id]["element"],
        }
        for map_id in map_ids
    }
    edge_observations = []
    for source_left, source_right, edge in sorted(source.edges(data=True), key=lambda item: tuple(sorted(item[:2]))):
        target_left, target_right = mapping[source_left], mapping[source_right]
        pair = tuple(sorted((target_left, target_right)))
        if not target.has_edge(*pair):
            raise ValueError(f"Mapped source edge {source_left}-{source_right} absent from target adjacency {pair}")
        left_atom, right_atom = target_atoms[target_left], target_atoms[target_right]
        dist = distance_angstrom(left_atom["coordinates_angstrom"], right_atom["coordinates_angstrom"])
        bond_class = edge["bond_class"]
        screen = interval_for_edge(source.nodes[source_left]["element"], source.nodes[source_right]["element"], bond_class)
        edge_observations.append({
            "source_map_ids": [source_left, source_right],
            "target_serials": [target_left, target_right],
            "target_atom_names": [left_atom["name"], right_atom["name"]],
            "element_pair": [source.nodes[source_left]["element"], source.nodes[source_right]["element"]],
            "source_bond_class_annotation": bond_class,
            "source_bond_order_numeric_annotation": edge["bond_order"],
            "target_pdb_adjacency_sources": target.edges[pair]["adjacency_sources"],
            "raw_distance_angstrom": dist,
            "broad_heuristic_screen_interval_angstrom": screen,
            "geometry_screen_only": True,
            "screen_result": "in_range" if screen[0] <= dist <= screen[1] else "out_of_range",
        })
    return {
        "mapping_id": f"mapping_{index:03d}",
        "map_id_to_target": map_records,
        "mapped_source_edge_distances": edge_observations,
        "geometry_screen_summary": {
            "edge_count": len(edge_observations),
            "in_range_count": sum(edge["screen_result"] == "in_range" for edge in edge_observations),
            "out_of_range_count": sum(edge["screen_result"] == "out_of_range" for edge in edge_observations),
            "minimum_distance_angstrom": min(edge["raw_distance_angstrom"] for edge in edge_observations),
            "maximum_distance_angstrom": max(edge["raw_distance_angstrom"] for edge in edge_observations),
            "interpretation": "Distances are observed after graph mapping; no mapping was selected or rejected by geometry.",
        },
    }, symmetry_key


def audit_state(
    state: str,
    model: dict[str, Any],
    source_graph: nx.Graph,
    source_summary: dict[str, Any],
    source_hashes: dict[str, str],
    automorphisms: list[dict[int, int]],
    symmetry_orbits: dict[tuple[int, ...], list[int]],
) -> dict[str, Any]:
    pdb_path = ROOT / model["output_file"]
    record: dict[str, Any] = {
        "state": state,
        "candidate_pdb": model["output_file"],
        "candidate_pdb_sha256": sha256(pdb_path),
        "expected_candidate_pdb_sha256_from_reconciliation": model["output_sha256"],
        "candidate_pdb_hash_matches_reconciliation": sha256(pdb_path) == model["output_sha256"],
        "source_graph_sha256": source_hashes["source_graph"],
        "map_legend_sha256": source_hashes["map_legend"],
        "failure_details": [],
    }
    try:
        if not record["candidate_pdb_hash_matches_reconciliation"]:
            raise ValueError(
                "Corrected candidate PDB hash differs from reconciliation: "
                f"expected {model['output_sha256']}, observed {record['candidate_pdb_sha256']}"
            )
        parsed = parse_pdb(pdb_path)
        target, selected_atoms, target_info = selected_target(parsed)
        record["pdb_record_checks"] = {
            "duplicate_serials_rejected": True,
            "all_atom_elements_present": True,
            "all_conect_references_resolved": True,
            "all_link_endpoints_resolved": True,
            "selected_atom_identities_unique": True,
            "raw_conect_reference_count": parsed["raw_conect_reference_count"],
            "unique_undirected_conect_edge_count_filewide": parsed["unique_undirected_conect_edge_count"],
        }
        record["target_graph"] = target_info["summary"]
        record["hydrogen_records"] = hydrogen_observations(parsed, target_info)
        mappings = target_mappings(source_graph, target)
        mapping_records = []
        symmetry_groups: dict[tuple[int, ...], list[str]] = defaultdict(list)
        for index, mapping in enumerate(mappings, start=1):
            result, symmetry_key = mapping_result(
                mapping, index, source_graph, target, selected_atoms, automorphisms
            )
            mapping_records.append(result)
            symmetry_groups[symmetry_key].append(result["mapping_id"])
        if not mappings:
            record["correspondence_status"] = "unresolved_no_exact_element_adjacency_mapping"
            record["failure_details"].append(
                "No complete element-labeled graph isomorphism covers every source and target node and edge."
            )
        else:
            record["correspondence_status"] = "exact_element_labeled_adjacency_mapping"
        record["source_graph"] = source_summary["graph_summary"]
        record["target_graph_summary"] = graph_summary(target)
        record["exact_mapping_count"] = len(mappings)
        record["all_exact_mappings"] = mapping_records
        record["source_graph_symmetry"] = {
            "automorphism_count_element_and_adjacency_preserving": len(automorphisms),
            "source_node_orbits": [list(orbit) for orbit in sorted(symmetry_orbits)],
            "mapping_groups_by_source_graph_symmetry": [
                {
                    "symmetry_group_id": f"symmetry_{index:02d}",
                    "mapping_ids": sorted(mapping_ids),
                    "representative_target_serial_vector_by_sorted_source_map": list(key),
                }
                for index, (key, mapping_ids) in enumerate(sorted(symmetry_groups.items()), start=1)
            ],
            "bond_orders_and_coordinates_not_used_for_symmetry_or_mapping": True,
        }
    except Exception as error:  # Preserve a state-specific audit failure in the JSON.
        record["correspondence_status"] = "unresolved_audit_failure"
        record["failure_details"].append(f"{type(error).__name__}: {error}")
    return record


def main() -> None:
    legend = load_json(LEGEND_PATH)
    reconciliation = load_json(RECONCILIATION_PATH)
    manifest = load_json(MANIFEST_PATH)
    source_hashes = {
        "source_graph": sha256(GRAPH_PATH),
        "map_legend": sha256(LEGEND_PATH),
        "coordinate_reconciliation": sha256(RECONCILIATION_PATH),
        "structure_manifest": sha256(MANIFEST_PATH),
    }
    molecule, source_graph, source_summary = parse_source_graph()
    source_automorphisms, source_symmetry_orbits = source_symmetry_groups(source_graph)
    fragment_mass = derive_fragment_mass()
    if fragment_mass["neutral_modified_fragment_formula"] != "C120H187N31O39S2":
        raise ValueError(f"Unexpected B1-B22 modified-fragment formula: {fragment_mass['neutral_modified_fragment_formula']}")
    if not fragment_mass["source_calculated_rounding_match"]:
        raise ValueError(
            f"Derived [M+2H]2+ {fragment_mass['calculated_m_z']} does not round to source 1376.2"
        )

    model_by_state = {model["state"]: model for model in reconciliation["models"]}
    if set(model_by_state) != set(STATE_ORDER):
        raise ValueError(f"Coordinate reconciliation states differ from expected: {sorted(model_by_state)}")
    manifest_by_state = {entry["state"]: entry for entry in manifest["files"]}
    states = [
        audit_state(
            state,
            model_by_state[state],
            source_graph,
            source_summary,
            source_hashes,
            source_automorphisms,
            source_symmetry_orbits,
        )
        for state in STATE_ORDER
    ]
    for state_record in states:
        original_model = manifest_by_state.get(state_record["state"])
        if original_model is None:
            state_record["failure_details"].append("No source structure manifest entry for state")
            state_record["correspondence_status"] = "unresolved_audit_failure"
        else:
            state_record["source_model_context"] = {
                "source_model_path": original_model["path"],
                "source_model_sha256": original_model["sha256"],
                "coordinate_status": original_model["coordinate_status"],
            }

    exact_count = sum(state["correspondence_status"] == "exact_element_labeled_adjacency_mapping" for state in states)
    output = {
        "schema": "dia-health-one.nnc2215-b1-topology-audit.v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "software": {
            "python_version": platform.python_version(),
            "rdkit_version": rdBase.rdkitVersion,
            "networkx_version": nx.__version__,
            "operating_platform": platform.platform(),
            "nvidia_gpu_used": False,
        },
        "input_hashes_sha256": {
            "source_graph": {"path": str(GRAPH_PATH.relative_to(ROOT)), "sha256": source_hashes["source_graph"]},
            "map_legend": {"path": str(LEGEND_PATH.relative_to(ROOT)), "sha256": source_hashes["map_legend"]},
            "coordinate_reconciliation": {"path": str(RECONCILIATION_PATH.relative_to(ROOT)), "sha256": source_hashes["coordinate_reconciliation"]},
            "structure_manifest": {"path": str(MANIFEST_PATH.relative_to(ROOT)), "sha256": source_hashes["structure_manifest"]},
            "corrected_candidate_pdbs": [
                {
                    "state": model["state"],
                    "path": model["output_file"],
                    "observed_sha256": sha256(ROOT / model["output_file"]),
                    "reconciliation_expected_sha256": model["output_sha256"],
                    "hash_matches": sha256(ROOT / model["output_file"]) == model["output_sha256"],
                }
                for model in (model_by_state[state] for state in STATE_ORDER)
            ],
        },
        "source_graph_validation": source_summary,
        "fragment_mass_check": fragment_mass,
        "source_graph_symmetry": {
            "automorphism_count_element_and_adjacency_preserving": len(source_automorphisms),
            "node_orbits": [list(orbit) for orbit in sorted(source_symmetry_orbits)],
            "bond_order_and_stereochemistry_not_used_to_select_coordinate_mapping": True,
        },
        "candidate_states": states,
        "overall_correspondence_status": (
            "both_candidates_exact" if exact_count == len(STATE_ORDER)
            else "partial_exact_correspondence" if exact_count
            else "unresolved"
        ),
        "candidate_pdbs_modified_by_audit": False,
        "limits": [
            "The corrected candidate coordinates are derivatives of publisher-supplied author models, not experimentally resolved NNC2215 structures.",
            "PDB LINK and CONECT records establish adjacency only; this audit does not infer their bond order, formal charge, or stereochemistry.",
            "Source bond classes are annotations transferred after exact graph mapping. Heuristic distance intervals are a geometry screen only and do not choose mappings.",
            "PheB1 backbone atoms beyond the N boundary endpoint are outside the B1 attachment component graph.",
            "The B1–B22 mass check is a fragment-composition benchmark under stated cleavage/reduction assumptions; it does not prove atom-by-atom coordinates.",
            "Physiological protonation, full-conjugate graph, formal charges, custom force-field parameters, and simulation stability remain unvalidated.",
            "No target-binding molecular dynamics readiness is established by this B1 topology audit.",
        ],
    }
    OUTPUT_PATH.write_text(json.dumps(output, indent=2, sort_keys=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT_PATH.relative_to(ROOT)),
        "overall_correspondence_status": output["overall_correspondence_status"],
        "fragment_formula": fragment_mass["neutral_modified_fragment_formula"],
        "fragment_mz": fragment_mass["calculated_m_z"],
        "candidate_results": [
            {
                "state": item["state"],
                "status": item["correspondence_status"],
                "mapping_count": item.get("exact_mapping_count", 0),
                "failure_details": item["failure_details"],
            }
            for item in states
        ],
    }, indent=2))


if __name__ == "__main__":
    main()
