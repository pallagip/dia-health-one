#!/usr/bin/env python3
"""Check candidate mappings from numbered Nature schemes to model PDB atoms.

This audits a bounded set of CONECT edges only. It does not infer bond orders,
charges, stereochemistry, protonation, or a production-ready force field.
"""

from __future__ import annotations

from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
MANIFEST = ROOT / "data" / "structure-manifest.json"
OUTPUT = ROOT / "data" / "nnc2215-linker-mapping-audit.json"
NATURE_FIG7 = (
    "https://media.springernature.com/original/springer-static/esm/"
    "art%3A10.1038%2Fs41586-024-08042-3/MediaObjects/"
    "41586_2024_8042_Fig12_ESM.jpg"
)
NATURE_FIG8 = (
    "https://media.springernature.com/original/springer-static/esm/"
    "art%3A10.1038%2Fs41586-024-08042-3/MediaObjects/"
    "41586_2024_8042_Fig13_ESM.jpg"
)
SOURCE_COMPONENTS = {
    "b29_linker": {"C": 6, "H": 6, "O": 2},
    "macrocycle_roof_propyl_azide": {"C": 56, "H": 63, "N": 15, "O": 12},
}
EXPECTED_B29 = Counter(SOURCE_COMPONENTS["b29_linker"])
EXPECTED_B29.update(SOURCE_COMPONENTS["macrocycle_roof_propyl_azide"])

STATES = {
    "active_open": {
        "path": "structures/nnc2215-open-model.pdb",
        "oxygen_carbonyl": "OD",
        "propyl_c22": "CD6",
    },
    "inactive_closed": {
        "path": "structures/nnc2215-closed-model.pdb",
        "oxygen_carbonyl": "O",
        "propyl_c22": "C",
    },
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_pdb(path: Path) -> tuple[dict[int, dict[str, str]], set[tuple[int, int]]]:
    atoms: dict[int, dict[str, str]] = {}
    edges: set[tuple[int, int]] = set()
    for line in path.read_text(encoding="ascii").splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            serial = int(line[6:11])
            if serial in atoms:
                raise ValueError(f"Duplicate atom serial {serial} in {path}")
            atoms[serial] = {
                "name": line[12:16].strip(),
                "residue": line[17:20].strip(),
                "chain": line[21].strip(),
                "residue_number": line[22:26].strip(),
                "element": line[76:78].strip().upper(),
            }
        elif line.startswith("CONECT"):
            serials = [
                int(line[index : index + 5])
                for index in range(6, len(line), 5)
                if line[index : index + 5].strip()
            ]
            if serials:
                edges.update(
                    tuple(sorted((serials[0], neighbor)))
                    for neighbor in serials[1:]
                )
    missing = sorted({serial for edge in edges for serial in edge if serial not in atoms})
    if missing:
        raise ValueError(f"CONECT references absent serials {missing} in {path}")
    return atoms, edges


def residue_atom(
    atoms: dict[int, dict[str, str]],
    chain: str,
    residue_number: str,
    residue: str,
    name: str,
) -> int:
    matches = [
        serial
        for serial, atom in atoms.items()
        if atom["chain"] == chain
        and atom["residue_number"] == residue_number
        and atom["residue"] == residue
        and atom["name"] == name
    ]
    if len(matches) != 1:
        raise ValueError(
            f"Expected one {chain}:{residue_number}:{residue}:{name}; got {matches}"
        )
    return matches[0]


def audit_state(label: str, state: dict[str, str], expected_hash: str) -> dict[str, object]:
    path = ROOT / state["path"]
    observed_hash = sha256(path)
    if observed_hash != expected_hash:
        raise ValueError(f"Source checksum mismatch for {path}")
    atoms, edges = parse_pdb(path)

    def atom(chain: str, residue_number: str, residue: str, name: str) -> int:
        return residue_atom(atoms, chain, residue_number, residue, name)

    lys_nz = atom("B", "29", "LYS", "NZ")
    l = lambda name: atom("B", "30", "LYX", name)
    mc = lambda name: atom("C", "1", "MC", name)
    map_ids = {
        "insulin_lys_b29_nz": lys_nz,
        "source_b29_carbon_3_carbonyl": l("CD1"),
        "source_b29_carbon_4": l("CD2"),
        "source_b29_carbon_5": l("C1"),
        "source_b29_carbon_6": l("C2"),
        "source_b29_carbon_7_alkyne": l("CD3"),
        "source_b29_carbon_8_alkyne": l("CD4"),
        "source_macrocycle_c21_candidate": mc("C31"),
        "source_macrocycle_c22_candidate": l(state["propyl_c22"]),
        "source_macrocycle_c23_candidate": l("CD5"),
        "source_azide_n_attached_to_c23_candidate": l("ND3"),
        "source_macrocycle_ring_anchor_candidate": mc("C71"),
    }

    bonds = [
        ("insulin_lys_b29_nz", "source_b29_carbon_3_carbonyl"),
        ("source_b29_carbon_3_carbonyl", "source_b29_carbonyl_oxygen"),
        ("source_b29_carbon_3_carbonyl", "source_b29_oxygen_link"),
        ("source_b29_oxygen_link", "source_b29_carbon_4"),
        ("source_b29_carbon_4", "source_b29_carbon_5"),
        ("source_b29_carbon_5", "source_b29_carbon_6"),
        ("source_b29_carbon_6", "source_b29_carbon_7_alkyne"),
        ("source_b29_carbon_7_alkyne", "source_b29_carbon_8_alkyne"),
        ("source_b29_carbon_7_alkyne", "source_triazole_n1_candidate"),
        ("source_triazole_n1_candidate", "source_triazole_n2_candidate"),
        ("source_triazole_n2_candidate", "source_azide_n_attached_to_c23_candidate"),
        ("source_azide_n_attached_to_c23_candidate", "source_b29_carbon_8_alkyne"),
        ("source_azide_n_attached_to_c23_candidate", "source_macrocycle_c23_candidate"),
        ("source_macrocycle_c23_candidate", "source_macrocycle_c22_candidate"),
        ("source_macrocycle_c22_candidate", "source_macrocycle_c21_candidate"),
        ("source_macrocycle_c21_candidate", "source_macrocycle_ring_anchor_candidate"),
    ]
    map_ids.update({
        "source_b29_carbonyl_oxygen": l(state["oxygen_carbonyl"]),
        "source_b29_oxygen_link": l("O1"),
        "source_triazole_n1_candidate": l("ND1"),
        "source_triazole_n2_candidate": l("ND2"),
    })
    missing_edges = [
        {"left": left, "right": right}
        for left, right in bonds
        if tuple(sorted((map_ids[left], map_ids[right]))) not in edges
    ]
    if missing_edges:
        raise ValueError(f"Candidate linker/tail adjacency failed in {label}: {missing_edges}")

    component_counts: dict[str, dict[str, int]] = {}
    for chain, residue_number, residue in (("B", "30", "LYX"), ("C", "1", "MC")):
        counts = Counter(
            atoms[serial]["element"]
            for serial in atoms
            if atoms[serial]["chain"] == chain
            and atoms[serial]["residue_number"] == residue_number
            and atoms[serial]["residue"] == residue
        )
        component_counts[residue] = dict(sorted(counts.items()))
    observed_heavy = Counter()
    for counts in component_counts.values():
        observed_heavy.update({key: value for key, value in counts.items() if key != "H"})
    expected_heavy = Counter({k: v for k, v in EXPECTED_B29.items() if k != "H"})
    delta = {
        element: observed_heavy[element] - expected_heavy[element]
        for element in sorted(set(observed_heavy) | set(expected_heavy))
        if observed_heavy[element] != expected_heavy[element]
    }

    source_to_model = {
        key: {
            "pdb_atom": atoms[serial],
            "serial": serial,
        }
        for key, serial in map_ids.items()
    }
    return {
        "state": label,
        "model_file": state["path"],
        "sha256": observed_hash,
        "candidate_source_to_model_atom_map": source_to_model,
        "candidate_adjacency_edges_verified": len(bonds),
        "candidate_adjacency_edges_missing": missing_edges,
        "observed_element_counts_by_model_residue": component_counts,
        "observed_b29_pool_heavy_elements": dict(sorted(observed_heavy.items())),
        "source_b29_pool_heavy_elements": dict(sorted(expected_heavy.items())),
        "observed_minus_source_heavy_element_delta": delta,
        "candidate_macrocycle_partition": {
            "source_macrocycle_c21_c23_tail_atoms": 3,
            "candidate_tail_atoms_in_lyx": 2,
            "candidate_tail_atom_in_mc": 1,
            "source_macrocycle_carbons_expected_in_mc_after_split": 54,
            "observed_mc_carbon_records": component_counts["MC"].get("C", 0),
            "candidate_mc_carbon_excess": component_counts["MC"].get("C", 0) - 54,
        },
    }


def main() -> None:
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    expected_by_path = {entry["path"]: entry["sha256"] for entry in manifest["files"]}
    states = [
        audit_state(name, state, expected_by_path[state["path"]])
        for name, state in STATES.items()
    ]
    result = {
        "schema": "dia-health-one.nnc2215-linker-mapping-audit.v1",
        "run": {
            "run_id": "DHO-NNC2215-LINK-MAP-" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ"),
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
            "script": "research/audit_nnc2215_linker_mapping.py",
            "script_sha256": sha256(Path(__file__)),
            "nvidia_gpu_used": False,
        },
        "sources": {
            "paper_doi": "10.1038/s41586-024-08042-3",
            "figure_7a_b29_alkyne_linker_url": NATURE_FIG7,
            "figure_8i_macrocycle_roof_propyl_azide_url": NATURE_FIG8,
            "figure_caption_evidence": "Figure 7a numbers carbons 3–8 in O-succinimidyl pentyn-1-oxycarbonyl; Figure 8i numbers the macrocycle roof propyl-azide chain carbons 21–23.",
            "source_b29_linker_formula_after_nhs_departure": "C6H6O2",
            "source_macrocycle_roof_propyl_azide_formula": "C56H63N15O12",
            "expected_b29_attachment_heavy_elements_after_click": dict(
                sorted({k: v for k, v in EXPECTED_B29.items() if k != "H"}.items())
            ),
        },
        "states": states,
        "assessment": {
            "candidate_b29_linker_mapping_supported_by_pdb_adjacency": all(
                item["candidate_adjacency_edges_verified"] == 16
                for item in states
            ),
            "candidate_roof_propyl_tail_mapping_supported_by_pdb_adjacency": all(
                item["candidate_macrocycle_partition"]["candidate_mc_carbon_excess"] == 2
                for item in states
            ),
            "interpretation": "The numbered source schemes and serial-level PDB adjacency support a candidate map for the six-carbon B29 alkyne linker and the three-carbon macrocycle propyl-azide tail. Under this map, two propyl-tail carbons are in LYX and one is in MC; the source macrocycle formula therefore leaves 54 carbons expected in the MC residue, while both model files contain 56. The two-carbon inventory surplus is localized to the MC residue under this candidate partition, but the two unsupported/duplicated MC atom records are not identified.",
            "not_established": [
                "a complete atom-mapped graph for either conjugate state",
                "which individual MC carbon records account for the two-carbon excess",
                "bond orders, formal charges, stereochemistry, or protonation",
                "parameter quality or force-field compatibility",
                "production-MD readiness",
            ],
            "target_binding_md_started": False,
            "clinical_use": False,
        },
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(OUTPUT.relative_to(ROOT)), "run": result["run"], "assessment": result["assessment"]}, indent=2))


if __name__ == "__main__":
    main()
