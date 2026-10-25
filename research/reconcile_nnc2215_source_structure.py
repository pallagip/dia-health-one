#!/usr/bin/env python3
"""Create source-grounded coordinate correction candidates for NNC2215.

The public open/closed PDBs encode two extra methylene atoms between two
macrocycle aryl carbons and the two C18 methyl groups shown in the published
chemical structure. This script removes those atoms, moves the methyls into
their direct aryl-substituent positions, repairs the open-model terminal
oxygen record, makes duplicate atom names unique, and records the evidence.
It does not assign bond orders, formal charges, protonation states, or
force-field terms.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent
MANIFEST = ROOT / "data" / "structure-manifest.json"
OUTPUT_JSON = ROOT / "data" / "nnc2215-carbon-reconciliation-candidate.json"

STATES = {
    "active_open": {
        "source": ROOT / "structures" / "nnc2215-open-model.pdb",
        "output": ROOT / "structures" / "nnc2215-open-source-corrected-candidate.pdb",
    },
    "inactive_closed": {
        "source": ROOT / "structures" / "nnc2215-closed-model.pdb",
        "output": ROOT / "structures" / "nnc2215-closed-source-corrected-candidate.pdb",
    },
}

EXPECTED_NNC_HEAVY = {"C": 325, "N": 79, "O": 97, "S": 6}

CORRECTIONS = (
    # ring aryl carbon, source-identified methyl carbon, intervening excess CH2
    ("C73", "C28", "C27", ("H23", "H24"), ("H25", "H26", "H27")),
    ("C75", "C30", "C29", ("H28", "H29"), ("H30", "H31", "H32")),
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def unit(vector: tuple[float, float, float]) -> tuple[float, float, float]:
    norm = math.sqrt(sum(value * value for value in vector))
    if norm < 1e-8:
        raise ValueError("Cannot normalize a zero-length vector.")
    return tuple(value / norm for value in vector)


def add(a: tuple[float, float, float], b: tuple[float, float, float]):
    return tuple(x + y for x, y in zip(a, b))


def sub(a: tuple[float, float, float], b: tuple[float, float, float]):
    return tuple(x - y for x, y in zip(a, b))


def scale(vector: tuple[float, float, float], factor: float):
    return tuple(value * factor for value in vector)


def dot(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return sum(x * y for x, y in zip(a, b))


def cross(a: tuple[float, float, float], b: tuple[float, float, float]):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def read_pdb(path: Path):
    lines = path.read_text(encoding="ascii").splitlines()
    atom_lines: dict[int, str] = {}
    atoms: dict[int, dict] = {}
    edges: set[tuple[int, int]] = set()
    for line in lines:
        record = line[:6].strip()
        if record in {"ATOM", "HETATM"}:
            serial = int(line[6:11])
            if serial in atoms:
                raise ValueError(f"Duplicate atom serial {serial} in {path}")
            atom_lines[serial] = line
            atoms[serial] = {
                "name": line[12:16].strip(),
                "residue": line[17:20].strip(),
                "chain": line[21].strip(),
                "resno": line[22:26].strip(),
                "element": line[76:78].strip().upper(),
                "xyz": tuple(float(line[i : i + 8]) for i in (30, 38, 46)),
            }
        elif record == "CONECT":
            serials = [
                int(line[i : i + 5])
                for i in range(6, len(line), 5)
                if line[i : i + 5].strip()
            ]
            if serials:
                edges.update(tuple(sorted((serials[0], other))) for other in serials[1:])
    for line in lines:
        if line.startswith("LINK  "):
            left = find_atom(
                atoms,
                line[12:16].strip(),
                line[17:20].strip(),
                line[21].strip(),
                line[22:26].strip(),
            )
            right = find_atom(
                atoms,
                line[42:46].strip(),
                line[47:50].strip(),
                line[51].strip(),
                line[52:56].strip(),
            )
            edges.add(tuple(sorted((left, right))))
        elif line.startswith("SSBOND"):
            left = find_atom(atoms, "SG", "CYS", line[15].strip(), line[17:21].strip())
            right = find_atom(atoms, "SG", "CYS", line[29].strip(), line[31:35].strip())
            edges.add(tuple(sorted((left, right))))
    invalid = sorted({serial for edge in edges for serial in edge if serial not in atoms})
    if invalid:
        raise ValueError(f"Connectivity references absent serials {invalid}")
    return lines, atom_lines, atoms, edges


def find_atom(atoms, name, residue, chain, resno) -> int:
    matches = [
        serial
        for serial, atom in atoms.items()
        if atom["name"] == name
        and atom["residue"] == residue
        and atom["chain"] == chain
        and atom["resno"].strip() == str(resno).strip()
    ]
    if len(matches) != 1:
        raise ValueError(
            f"Expected one atom {chain}:{resno}:{residue}:{name}; found {matches}"
        )
    return matches[0]


def neighbors(serial: int, edges: set[tuple[int, int]]) -> set[int]:
    return {
        right if left == serial else left
        for left, right in edges
        if serial in (left, right)
    }


def format_xyz(line: str, xyz: tuple[float, float, float]) -> str:
    return line[:30] + "".join(f"{value:8.3f}" for value in xyz) + line[54:]


def set_methyl_hydrogens(
    xyz_carbon: tuple[float, float, float],
    xyz_ring: tuple[float, float, float],
    phase_degrees: float = 0.0,
) -> list[tuple[float, float, float]]:
    axis = unit(sub(xyz_carbon, xyz_ring))
    reference = (1.0, 0.0, 0.0) if abs(axis[0]) < 0.8 else (0.0, 1.0, 0.0)
    first = unit(cross(axis, reference))
    second = cross(axis, first)
    # The axis points from the ring carbon to the methyl carbon. For a
    # tetrahedral methyl, each C-H vector therefore has a +1/3 projection on
    # this axis (the C-ring bond points in the opposite direction).
    axial = 1.0 / 3.0
    radial = math.sqrt(8.0 / 9.0)
    result = []
    for degrees in (0.0, 120.0, 240.0):
        angle = math.radians(degrees + phase_degrees)
        direction = add(
            scale(axis, axial),
            add(scale(first, radial * math.cos(angle)), scale(second, radial * math.sin(angle))),
        )
        result.append(add(xyz_carbon, scale(direction, 1.09)))
    return result


def optimize_methyl_hydrogen_phase(
    xyz_carbon: tuple[float, float, float],
    xyz_ring: tuple[float, float, float],
    atoms: dict[int, dict],
    excluded_serials: set[int],
) -> tuple[list[tuple[float, float, float]], float, list[dict]]:
    """Rotate the equivalent methyl Hs to reduce severe local overlaps."""
    cutoffs = {"H": 1.85, "C": 2.20, "N": 2.05, "O": 2.00, "S": 2.35}
    best = None
    for phase in range(120):
        hydrogen_xyzs = set_methyl_hydrogens(xyz_carbon, xyz_ring, float(phase))
        score = 0.0
        contacts = []
        for hydrogen_xyz in hydrogen_xyzs:
            for serial, atom in atoms.items():
                element = atom["element"]
                if serial in excluded_serials or element not in cutoffs:
                    continue
                distance = math.dist(hydrogen_xyz, atom["xyz"])
                cutoff = cutoffs[element]
                overlap = max(0.0, cutoff - distance)
                score += overlap * overlap * (2.0 if element == "H" else 1.0)
                if overlap > 0.0:
                    contacts.append(
                        {
                            "other_serial": serial,
                            "other_name": atom["name"],
                            "other_residue": f"{atom['chain']}:{atom['resno'].strip()}:{atom['residue']}",
                            "other_element": element,
                            "distance_angstrom": round(distance, 3),
                            "cutoff_angstrom": cutoff,
                        }
                    )
        candidate = (score, phase, hydrogen_xyzs, contacts)
        if best is None or candidate[:2] < best[:2]:
            best = candidate
    assert best is not None
    score, phase, hydrogen_xyzs, contacts = best
    contacts.sort(key=lambda item: item["distance_angstrom"])
    return hydrogen_xyzs, float(phase), contacts[:8]


def write_conect(edges: set[tuple[int, int]]) -> list[str]:
    adjacency: dict[int, set[int]] = defaultdict(set)
    for left, right in edges:
        adjacency[left].add(right)
        adjacency[right].add(left)
    records = []
    for source in sorted(adjacency):
        targets = sorted(adjacency[source])
        for offset in range(0, len(targets), 4):
            records.append(
                f"CONECT{source:5d}"
                + "".join(f"{target:5d}" for target in targets[offset : offset + 4])
            )
    return records


def inventory(atoms: dict[int, dict], excluded_residue: str | None = None) -> dict[str, int]:
    counts = Counter(
        atom["element"]
        for atom in atoms.values()
        if atom["element"] and atom["residue"] != excluded_residue
    )
    return dict(sorted(counts.items()))


def apply_state(label: str, state: dict, expected_sha: str) -> dict:
    source = state["source"]
    output = state["output"]
    observed_sha = sha256(source)
    if observed_sha != expected_sha:
        raise ValueError(f"Source checksum mismatch for {source}: {observed_sha}")
    lines, atom_lines, atoms, edges = read_pdb(source)
    original_inventory = inventory(atoms, excluded_residue="BGL")
    mc = {
        atom["name"]: serial
        for serial, atom in atoms.items()
        if atom["residue"] == "MC" and atom["chain"] == "C" and atom["resno"].strip() == "1"
    }
    removed_serials: set[int] = set()
    moved: list[dict] = []
    changed_atom_lines = dict(atom_lines)
    new_edges = set(edges)
    inserted_atom_lines: list[str] = []
    terminal_oxygen_reconstruction = None
    terminal_atom_insertion_after_serial = None
    planned_removed_serials: set[int] = set()
    for _, _, extra_name, _, _ in CORRECTIONS:
        extra_serial = mc[extra_name]
        planned_removed_serials.add(extra_serial)
        planned_removed_serials.update(
            serial
            for serial in neighbors(extra_serial, edges)
            if atoms[serial]["element"] == "H"
        )

    for ring_name, methyl_name, extra_name, extra_h_names, methyl_h_names in CORRECTIONS:
        ring = mc[ring_name]
        methyl = mc[methyl_name]
        extra = mc[extra_name]
        extra_neighbors = neighbors(extra, edges)
        methyl_neighbors = neighbors(methyl, edges)
        expected_extra = {ring, methyl}
        extra_h = {
            serial for serial in extra_neighbors if atoms[serial]["element"] == "H"
        }
        methyl_h = {
            serial for serial in methyl_neighbors if atoms[serial]["element"] == "H"
        }
        if not expected_extra.issubset(extra_neighbors):
            raise ValueError(f"{extra_name} is not between {ring_name} and {methyl_name}")
        if len(extra_h) != 2 or len(methyl_h) != 3:
            raise ValueError(f"Unexpected hydrogen ownership around {extra_name}/{methyl_name}")
        if {atoms[serial]["name"] for serial in extra_h} != set(extra_h_names):
            raise ValueError(f"Unexpected hydrogen names on {extra_name}")
        if {atoms[serial]["name"] for serial in methyl_h} != set(methyl_h_names):
            raise ValueError(f"Unexpected hydrogen names on {methyl_name}")
        if tuple(sorted((ring, extra))) not in edges or tuple(sorted((extra, methyl))) not in edges:
            raise ValueError(f"Missing source-model bonds around {extra_name}")

        ring_xyz = atoms[ring]["xyz"]
        extra_xyz = atoms[extra]["xyz"]
        new_methyl_xyz = add(ring_xyz, scale(unit(sub(extra_xyz, ring_xyz)), 1.50))
        hydrogen_xyzs, hydrogen_phase, close_contacts = optimize_methyl_hydrogen_phase(
            new_methyl_xyz,
            ring_xyz,
            atoms,
            planned_removed_serials | {ring, methyl} | methyl_h,
        )
        ordered_methyl_h = sorted(methyl_h, key=lambda serial: atoms[serial]["name"])
        changed_atom_lines[methyl] = format_xyz(atom_lines[methyl], new_methyl_xyz)
        atoms[methyl]["xyz"] = new_methyl_xyz
        for serial, xyz in zip(ordered_methyl_h, hydrogen_xyzs):
            changed_atom_lines[serial] = format_xyz(atom_lines[serial], xyz)
            atoms[serial]["xyz"] = xyz

        deleted = {extra, *extra_h}
        removed_serials.update(deleted)
        new_edges = {
            edge for edge in new_edges if edge[0] not in deleted and edge[1] not in deleted
        }
        new_edges.add(tuple(sorted((ring, methyl))))
        moved.append(
            {
                "removed_extra_atom": {
                    "name": extra_name,
                    "serial": extra,
                    "element": atoms[extra]["element"],
                    "removed_hydrogens": sorted(atoms[serial]["name"] for serial in extra_h),
                },
                "source_methyl_atom": {"name": methyl_name, "serial": methyl},
                "methyl_hydrogen_serials": sorted(methyl_h),
                "aryl_atom": {"name": ring_name, "serial": ring},
                "new_bond": [ring_name, methyl_name],
                "methyl_hydrogens_retained": sorted(atoms[serial]["name"] for serial in methyl_h),
                "new_aryl_methyl_distance_angstrom": round(math.dist(ring_xyz, new_methyl_xyz), 3),
                "methyl_c_h_distances_angstrom": [
                    round(math.dist(new_methyl_xyz, xyz), 3) for xyz in hydrogen_xyzs
                ],
                "methyl_hydrogen_phase_degrees": hydrogen_phase,
                "closest_local_contacts_below_screening_cutoffs": close_contacts,
            }
        )

    if label == "active_open":
        a21 = {
            atom["name"]: serial
            for serial, atom in atoms.items()
            if atom["residue"] == "ASN"
            and atom["chain"] == "A"
            and atom["resno"].strip() == "21"
        }
        if "OXT" in a21:
            raise ValueError("The open model unexpectedly already contains ASN A21 OXT")
        for required_name in ("CA", "C", "O", "HXT"):
            if required_name not in a21:
                raise ValueError(f"The open ASN A21 terminal group lacks {required_name}")

        ca_serial, c_serial, o_serial, hxt_serial = (
            a21["CA"], a21["C"], a21["O"], a21["HXT"]
        )
        terminal_atom_insertion_after_serial = max(a21.values())
        ca_xyz = atoms[ca_serial]["xyz"]
        c_xyz = atoms[c_serial]["xyz"]
        o_xyz = atoms[o_serial]["xyz"]
        hxt_xyz = atoms[hxt_serial]["xyz"]
        hxt_c_distance = math.dist(hxt_xyz, c_xyz)
        if atoms[hxt_serial]["element"] != "H" or not 0.95 <= hxt_c_distance <= 1.20:
            raise ValueError(
                "Open ASN A21 HXT does not match the short terminal C-H record "
                f"seen in the source model ({hxt_c_distance:.3f} Å)"
            )

        # Rebuild a planar carboxylate from the observed CA-C-O plane. The
        # source HXT lies at the missing OXT site but is only 1.08 Å from C;
        # retaining it beside an OXT would create a severe overlap.
        ca_direction = unit(sub(ca_xyz, c_xyz))
        o_direction = unit(sub(o_xyz, c_xyz))
        oxt_direction = unit(scale(add(ca_direction, o_direction), -1.0))
        oxt_xyz = add(c_xyz, scale(oxt_direction, 1.25))
        oxt_serial = terminal_atom_insertion_after_serial + 1
        o_template = atom_lines[o_serial]
        oxt_line = (
            o_template[:6]
            + f"{oxt_serial:5d}"
            + o_template[11:12]
            + " OXT"
            + o_template[16:30]
            + "".join(f"{value:8.3f}" for value in oxt_xyz)
            + o_template[54:]
        )
        atoms[oxt_serial] = {
            "name": "OXT",
            "residue": "ASN",
            "chain": "A",
            "resno": "  21",
            "element": "O",
            "xyz": oxt_xyz,
        }
        inserted_atom_lines.append(oxt_line)
        removed_serials.add(hxt_serial)
        new_edges = {
            edge for edge in new_edges if hxt_serial not in edge
        }
        new_edges.add(tuple(sorted((c_serial, oxt_serial))))
        terminal_oxygen_reconstruction = {
            "removed_source_atom": "ASN A21 HXT",
            "removed_source_element": "H",
            "source_C_HXT_distance_angstrom": round(hxt_c_distance, 3),
            "added_atom": "ASN A21 OXT",
            "added_serial": oxt_serial,
            "added_element": "O",
            "reconstruction": "placed in the CA-C-O plane opposite the CA and O vectors, at a 1.25 Å C-OXT distance",
            "new_C_OXT_distance_angstrom": round(math.dist(c_xyz, oxt_xyz), 3),
            "O_C_OXT_angle_degrees": round(
                math.degrees(
                    math.acos(
                        max(-1.0, min(1.0, dot(unit(sub(o_xyz, c_xyz)), unit(sub(oxt_xyz, c_xyz)))))
                    )
                ),
                2,
            ),
            "protonation_assumption": "terminal carboxylate; no terminal OXT hydrogen",
        }

    atom_name_mapping = []
    by_residue: dict[tuple[str, str, str], dict[str, list[int]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for serial, atom in atoms.items():
        if serial in removed_serials:
            continue
        residue_key = (
            atom["chain"],
            atom["residue"],
            atom["resno"].strip(),
        )
        by_residue[residue_key][atom["name"]].append(serial)

    for residue_key, names_to_serials in by_residue.items():
        reserved_names = set(names_to_serials)
        used_names = set(reserved_names)
        for old_name, serials in sorted(names_to_serials.items()):
            if len(serials) < 2:
                continue
            suffixes = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"
            stem = old_name[:3]
            for index, serial in enumerate(sorted(serials)):
                candidate_name = None
                for suffix in suffixes:
                    trial = f"{stem}{suffix}"
                    if len(trial) <= 4 and trial not in used_names:
                        candidate_name = trial
                        break
                if candidate_name is None:
                    raise ValueError(
                        f"Cannot disambiguate duplicate atom name {old_name} in {residue_key}"
                    )
                used_names.add(candidate_name)
                atoms[serial]["name"] = candidate_name
                if serial not in changed_atom_lines:
                    changed_atom_lines[serial] = atom_lines[serial]
                old_line = changed_atom_lines[serial]
                changed_atom_lines[serial] = (
                    old_line[:12] + f"{candidate_name:>4s}" + old_line[16:]
                )
                atom_name_mapping.append(
                    {
                        "serial": serial,
                        "residue": list(residue_key),
                        "source_atom_name": old_name,
                        "candidate_atom_name": candidate_name,
                    }
                )

    final_names = [
        (atom["chain"], atom["residue"], atom["resno"].strip(), atom["name"])
        for serial, atom in atoms.items()
        if serial not in removed_serials
    ]
    if len(final_names) != len(set(final_names)):
        raise ValueError("Candidate still contains duplicate atom names within a residue")

    output_atoms = {
        serial: atom for serial, atom in atoms.items() if serial not in removed_serials
    }
    output_serial_map = {
        serial: (
            serial + 1
            if terminal_atom_insertion_after_serial is not None
            and serial != terminal_oxygen_reconstruction["added_serial"]
            and serial > terminal_atom_insertion_after_serial
            else serial
        )
        for serial in output_atoms
    }
    corrected_inventory = inventory(output_atoms, excluded_residue="BGL")
    corrected_heavy_inventory = {
        key: value for key, value in corrected_inventory.items() if key != "H"
    }
    if corrected_heavy_inventory != EXPECTED_NNC_HEAVY:
        raise ValueError(
            "Corrected full-conjugate heavy-atom inventory did not match the "
            f"reported formula: {corrected_heavy_inventory}"
        )
    if corrected_inventory.get("C") != original_inventory.get("C", 0) - 2:
        raise ValueError("The correction did not remove exactly two carbon records.")
    pool = Counter(
        atom["element"]
        for atom in output_atoms.values()
        if (atom["residue"], atom["chain"], atom["resno"].strip())
        in {("LYX", "B", "30"), ("MC", "C", "1")}
        and atom["element"] != "H"
    )
    if dict(sorted(pool.items())) != {"C": 62, "N": 15, "O": 14}:
        raise ValueError(f"Corrected B29 pool did not match source formula: {dict(pool)}")

    retained = [
        line
        for line in lines
        if not line.startswith("CONECT")
        and not (
            line.startswith(("ATOM  ", "HETATM"))
            and int(line[6:11]) in removed_serials
        )
    ]
    remarks = [
        "REMARK 950 SOURCE-BASED CANDIDATE: removed two excess MC methylenes.",
        "REMARK 950 Reconnected MC C28/C30 directly to MC C73/C75 as methyls.",
        "REMARK 950 Rebuilt methyl H geometry; unique names assigned within residues.",
        "REMARK 950 See candidate audit JSON for source evidence and corrections.",
        "REMARK 950 Bond orders, charges and force-field terms remain unvalidated.",
    ]
    if terminal_oxygen_reconstruction is not None:
        remarks.insert(3, "REMARK 950 Open ASN A21 OXT rebuilt from the CA-C-O plane.")
    insert_at = next((i for i, line in enumerate(retained) if line.startswith("ATOM  ") or line.startswith("HETATM")), 0)
    retained[insert_at:insert_at] = remarks
    for index, line in enumerate(retained):
        if line.startswith(("ATOM  ", "HETATM")):
            serial = int(line[6:11])
            retained[index] = changed_atom_lines[serial]
    for added_line in inserted_atom_lines:
        insertion_index = max(
            index
            for index, line in enumerate(retained)
            if line.startswith("ATOM  ")
            and line[21].strip() == "A"
            and line[17:20].strip() == "ASN"
            and line[22:26].strip() == "21"
        ) + 1
        retained.insert(insertion_index, added_line)

    # Inserting OXT before chain A's TER requires shifting subsequent atom and
    # TER serials so PDB records remain ordered and unique.
    for index, line in enumerate(retained):
        record = line[:6].strip()
        if record in {"ATOM", "HETATM"}:
            old_serial = int(line[6:11])
            new_serial = output_serial_map[old_serial]
            retained[index] = line[:6] + f"{new_serial:5d}" + line[11:]
        elif record == "TER" and terminal_atom_insertion_after_serial is not None:
            ter_serial = int(line[6:11])
            if ter_serial > terminal_atom_insertion_after_serial:
                retained[index] = line[:6] + f"{ter_serial + 1:5d}" + line[11:]

    remapped_edges = {
        tuple(sorted((output_serial_map[left], output_serial_map[right])))
        for left, right in new_edges
    }
    conect = write_conect(remapped_edges)
    end_index = next((i for i, line in enumerate(retained) if line.startswith("END")), len(retained))
    retained[end_index:end_index] = conect
    output.write_text("\n".join(retained) + "\n", encoding="ascii")

    # Re-open the exact file written above and verify the serialized candidate,
    # including PDB connectivity references and coordinate geometry.
    written_lines, _, written_atoms, written_edges = read_pdb(output)
    if len(written_atoms) != len(output_atoms):
        raise ValueError("Written candidate atom count changed during PDB serialization")
    written_heavy_inventory = {
        key: value
        for key, value in inventory(written_atoms, excluded_residue="BGL").items()
        if key != "H"
    }
    if written_heavy_inventory != EXPECTED_NNC_HEAVY:
        raise ValueError(f"Written candidate heavy inventory is incorrect: {written_heavy_inventory}")
    written_name_counts = Counter(
        (atom["chain"], atom["residue"], atom["resno"].strip(), atom["name"])
        for atom in written_atoms.values()
    )
    written_duplicate_names = sum(count - 1 for count in written_name_counts.values() if count > 1)
    if written_duplicate_names:
        raise ValueError(f"Written candidate has {written_duplicate_names} duplicate atom names")
    written_atom_serials = [
        int(line[6:11]) for line in output.read_text(encoding="ascii").splitlines()
        if line.startswith(("ATOM  ", "HETATM"))
    ]
    serials_monotonic = all(
        left < right for left, right in zip(written_atom_serials, written_atom_serials[1:])
    )
    if not serials_monotonic:
        raise ValueError("Written candidate atom serials are not strictly increasing")
    max_line_length = max(map(len, written_lines), default=0)
    if max_line_length > 80:
        raise ValueError(f"Written candidate has a {max_line_length}-column PDB record")

    methyl_bond_checks = []
    for correction in moved:
        ring_serial = output_serial_map[correction["aryl_atom"]["serial"]]
        methyl_serial = output_serial_map[correction["source_methyl_atom"]["serial"]]
        bond = tuple(sorted((ring_serial, methyl_serial)))
        if bond not in written_edges:
            raise ValueError(f"Written PDB is missing aryl-methyl edge {bond}")
        ring_methyl_distance = math.dist(
            written_atoms[ring_serial]["xyz"], written_atoms[methyl_serial]["xyz"]
        )
        hydrogen_distances = [
            math.dist(written_atoms[methyl_serial]["xyz"], written_atoms[serial]["xyz"])
            for serial in (
                output_serial_map[source_serial]
                for source_serial in correction["methyl_hydrogen_serials"]
            )
        ]
        if not math.isclose(ring_methyl_distance, 1.50, abs_tol=0.01):
            raise ValueError(f"Written aryl-methyl distance is {ring_methyl_distance:.3f} Å")
        if any(not math.isclose(distance, 1.09, abs_tol=0.01) for distance in hydrogen_distances):
            raise ValueError(f"Written methyl C-H distances are incorrect: {hydrogen_distances}")
        methyl_bond_checks.append(
            {
                "aryl_methyl_distance_angstrom": round(ring_methyl_distance, 3),
                "methyl_c_h_distances_angstrom": [round(distance, 3) for distance in hydrogen_distances],
                "bond_record_present": True,
                "methyl_hydrogen_overlap_screen_passed": not correction["closest_local_contacts_below_screening_cutoffs"],
            }
        )

    terminal_oxygen_check = None
    if terminal_oxygen_reconstruction is not None:
        oxygen_serial = terminal_oxygen_reconstruction["added_serial"]
        carbon_serial = find_atom(written_atoms, "C", "ASN", "A", "21")
        oxygen_bond = tuple(sorted((carbon_serial, oxygen_serial)))
        if oxygen_bond not in written_edges:
            raise ValueError("Written candidate is missing ASN A21 C-OXT connectivity")
        terminal_oxygen_check = {
            "C_OXT_distance_angstrom": round(
                math.dist(written_atoms[carbon_serial]["xyz"], written_atoms[oxygen_serial]["xyz"]), 3
            ),
            "O_C_OXT_angle_degrees": terminal_oxygen_reconstruction["O_C_OXT_angle_degrees"],
            "OXT_bond_record_present": True,
            "source_HXT_removed": not any(
                atom["name"] == "HXT" and atom["residue"] == "ASN" and atom["chain"] == "A" and atom["resno"].strip() == "21"
                for atom in written_atoms.values()
            ),
        }
        if terminal_oxygen_check["C_OXT_distance_angstrom"] != 1.25:
            raise ValueError("Written ASN A21 C-OXT distance is incorrect")

    coordinate_candidate_checks = {
        "serialized_pdb_reopened": True,
        "atom_count": len(written_atoms),
        "atom_serials_strictly_increasing": serials_monotonic,
        "pdb_records_at_most_80_columns": max_line_length <= 80,
        "conect_link_ssbond_references_valid": True,
        "heavy_inventory_matches_reported_formula": True,
        "duplicate_atom_names_within_residue": written_duplicate_names,
        "aryl_methyl_geometry": methyl_bond_checks,
        "asn_a21_terminal_oxygen": terminal_oxygen_check,
    }

    return {
        "state": label,
        "source_file": str(source.relative_to(ROOT)),
        "source_sha256": observed_sha,
        "output_file": str(output.relative_to(ROOT)),
        "output_sha256": sha256(output),
        "source_heavy_inventory_excluding_free_glucose": {
            key: value for key, value in original_inventory.items() if key != "H"
        },
        "corrected_heavy_inventory_excluding_free_glucose": corrected_heavy_inventory,
        "corrected_b29_linked_pool": dict(sorted(pool.items())),
        "corrections": moved,
        "duplicate_atom_name_mapping": atom_name_mapping,
        "terminal_oxygen_reconstruction": terminal_oxygen_reconstruction,
        "coordinate_candidate_checks": coordinate_candidate_checks,
        "full_conjugate_heavy_inventory_difference": {
            "source_formula_expected": EXPECTED_NNC_HEAVY,
            "coordinate_model_excluding_free_glucose": corrected_heavy_inventory,
            "difference": {
                element: corrected_heavy_inventory.get(element, 0) - expected
                for element, expected in EXPECTED_NNC_HEAVY.items()
            },
        },
    }


def main() -> None:
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    expected_hashes = {item["path"]: item["sha256"] for item in manifest["files"]}
    results = []
    for label, state in STATES.items():
        relative = str(state["source"].relative_to(ROOT))
        results.append(apply_state(label, state, expected_hashes[relative]))
    result = {
        "schema": "dia-health-one.nnc2215-coordinate-reconciliation-candidate.v2",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "source_chemistry": {
            "paper_doi": "10.1038/s41586-024-08042-3",
            "figure": "Extended Data Fig. 8i",
            "supplementary_formula": "C56H63N15O12",
            "nmr_assignment": "6H, H-18, corresponding to the two aryl methyl substituents",
            "source_formula_and_synthesis": "C325H461N79O97S6 for NNC2215; B29 linker C6H6O2; macrocycle roof propyl azide C56H63N15O12",
            "coordinate_source": {
                "authors": "Eva Johansson; Daniele Granata",
                "title": "3D models of NNC2215",
                "doi": "10.6084/m9.figshare.26526460",
                "license": "CC BY 4.0",
                "modification_note": "Candidate PDBs are modified derivatives of the author-supplied model coordinates, not original or experimentally solved structures.",
            },
            "sources": [
                "https://www.nature.com/articles/s41586-024-08042-3",
                "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-024-08042-3/MediaObjects/41586_2024_8042_MOESM1_ESM.pdf",
                "https://springernature.figshare.com/articles/dataset/3D_models_of_NNC2215/26526460",
            ],
        },
        "mapping_inference": {
            "pdb_extra_atoms": ["MC C27", "MC C29"],
            "source_methyl_atoms": ["MC C28", "MC C30"],
            "source_aryl_atoms": ["MC C73", "MC C75"],
            "candidate_correction": "Delete C27/C29 and their four hydrogens; connect C28/C30 directly to C73/C75; retain C28/C30 as methyl carbons and rebuild their hydrogen coordinates.",
            "basis": "The source 2D chemical structure depicts two direct ring methyls, the SI assigns a six-proton H-18 signal to the two methyl groups, and this correction brings the B29 linked heavy-atom pool to the source-derived C62N15O14 in both models.",
            "status": "source-supported mapping inference; requires independent chemistry review before production parameterization",
        },
        "models": results,
        "readiness": {
            "coordinate_candidate_checks_passed": True,
            "two_carbon_discrepancy_resolved_as_candidate_mapping": True,
            "reported_full_conjugate_heavy_inventory_matched": True,
            "full_molecular_graph_validated": False,
            "bond_orders_assigned": False,
            "formal_charges_assigned": False,
            "protonation_inventory_validated": False,
            "force_field_parameters_generated": False,
            "target_md_ready": False,
        },
        "limitations": [
            "PDB CONECT records provide adjacency only and do not encode bond order or formal charge.",
            "Coordinates for the methyl hydrogens are rebuilt geometrically; heavy atom ring coordinates are retained.",
            "The open-model ASN A21 terminal OXT was reconstructed from the CA-C-O plane after removing a misplaced HXT record at the OXT position; this geometry and the terminal carboxylate microstate are model assumptions.",
            "Hydrogen coordinates and the full protonation inventory remain unvalidated and should be rebuilt from the resolved molecular graph.",
            "This is a source-reconciled coordinate candidate, not a parameterized simulation system.",
        ],
        "clinical_use": False,
    }
    OUTPUT_JSON.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(OUTPUT_JSON.relative_to(ROOT)), "readiness": result["readiness"], "models": results}, indent=2))


if __name__ == "__main__":
    main()
