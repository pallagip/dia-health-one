#!/usr/bin/env python3
"""Create a checksum-checked PDB record inventory for the published NNC2215 models."""

from __future__ import annotations

from collections import Counter, defaultdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform


ROOT = Path(__file__).resolve().parents[1]
MANIFEST_PATH = ROOT / "research" / "data" / "structure-manifest.json"
OUTPUT_PATH = ROOT / "research" / "data" / "structure-preflight.json"
STANDARD_AMINO_ACIDS = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}
WATER_NAMES = {"HOH", "WAT"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def inspect_pdb(path: Path) -> dict:
    record_counts: Counter[str] = Counter()
    atom_record_counts: Counter[str] = Counter()
    chain_atom_counts: Counter[str] = Counter()
    residues: dict[str, set[tuple[str, str]]] = defaultdict(set)
    residue_names: dict[str, set[str]] = defaultdict(set)
    nonstandard_hetero_names: set[str] = set()
    disulfides: list[str] = []
    covalent_links: list[str] = []

    with path.open(encoding="ascii", errors="replace") as handle:
        for line in handle:
            record = line[:6].strip()
            record_counts[record] += 1
            if record in {"ATOM", "HETATM"}:
                resname = line[17:20].strip()
                chain = line[21:22].strip() or "_"
                residue_id = line[22:26].strip() + line[26:27].strip()
                atom_record_counts[record] += 1
                chain_atom_counts[chain] += 1
                residues[chain].add((residue_id, resname))
                residue_names[chain].add(resname)
                if record == "HETATM" and resname not in STANDARD_AMINO_ACIDS | WATER_NAMES:
                    nonstandard_hetero_names.add(resname)
            elif record == "SSBOND":
                disulfides.append(line.rstrip())
            elif record == "LINK":
                covalent_links.append(line.rstrip())

    atom_total = sum(atom_record_counts.values())
    if atom_total == 0:
        raise ValueError(f"No ATOM/HETATM coordinates found in {path}")

    return {
        "coordinate_atom_records": atom_total,
        "atom_records_by_type": dict(sorted(atom_record_counts.items())),
        "chain_atom_record_counts": dict(sorted(chain_atom_counts.items())),
        "residue_count_by_chain": {
            chain: len(items) for chain, items in sorted(residues.items())
        },
        "residue_names_by_chain": {
            chain: sorted(names) for chain, names in sorted(residue_names.items())
        },
        "nonstandard_hetero_residue_names": sorted(nonstandard_hetero_names),
        "connectivity_record_counts": {
            name: record_counts[name] for name in ("SSBOND", "LINK", "CONECT")
        },
        "ssbond_records": disulfides,
        "link_records": covalent_links,
    }


def main() -> None:
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    structures = []

    for source in manifest["files"]:
        path = ROOT / "research" / source["path"]
        local_sha256 = sha256(path)
        if local_sha256 != source["sha256"]:
            raise SystemExit(f"SHA-256 mismatch for {path}: {local_sha256}")

        inventory = inspect_pdb(path)
        if inventory["coordinate_atom_records"] != source["coordinate_atom_records"]:
            raise SystemExit(f"Coordinate count mismatch for {path}")
        if sorted(inventory["chain_atom_record_counts"]) != source["chain_ids"]:
            raise SystemExit(f"Chain ID mismatch for {path}")
        if inventory["nonstandard_hetero_residue_names"] != source["nonstandard_hetero_residue_labels"]:
            raise SystemExit(f"Nonstandard residue-label mismatch for {path}")

        structures.append({
            "state": source["state"],
            "file": source["path"],
            "sha256": local_sha256,
            **inventory,
        })

    result = {
        "record_type": "structure_coordinate_preflight",
        "schema": "dia-health-one.structure-preflight.v1",
        "run": {
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
            "script": "research/inspect_structure_models.py",
            "script_sha256": sha256(Path(__file__).resolve()),
            "runtime": platform.python_version(),
        },
        "source_manifest": "research/data/structure-manifest.json",
        "source_dataset_doi": manifest["dataset_doi"],
        "structures": structures,
        "interpretation": (
            "This is a coordinate and PDB-record inventory cross-checked against the source manifest. "
            "It does not validate chemical graphs, missing atoms, protonation states, force-field "
            "parameters, or simulation readiness."
        ),
    }
    OUTPUT_PATH.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {OUTPUT_PATH.relative_to(ROOT)}")
    for item in structures:
        print(
            f"{item['state']}: {item['coordinate_atom_records']} coordinate records; "
            f"chains {', '.join(item['chain_atom_record_counts'])}; "
            f"{item['connectivity_record_counts']['SSBOND']} SSBOND, "
            f"{item['connectivity_record_counts']['LINK']} LINK, "
            f"{item['connectivity_record_counts']['CONECT']} CONECT records"
        )


if __name__ == "__main__":
    main()
