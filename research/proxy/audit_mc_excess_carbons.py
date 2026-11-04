#!/usr/bin/env python3
"""Locate the two excess MC carbon records in the public NNC2215 model PDBs.

Compares the PDB heavy-atom graph of residue MC (+ LYX tail) with the
macrocycle roof propyl-azide drawn and numbered in Nature 2024 Extended Data
Fig. 8i (formula C56H63N15O12). For every hexasubstituted benzene ring the
script lists the substituent types and the carbon count of each alkyl
substituent. Read-only; writes a JSON evidence record.
"""
from __future__ import annotations
import hashlib, json, sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
import networkx as nx

ROOT = Path(__file__).resolve().parents[1]
FILES = {"active_open": ROOT / "structures/nnc2215-open-model.pdb",
         "inactive_closed": ROOT / "structures/nnc2215-closed-model.pdb"}

def parse(path):
    atoms, g = {}, nx.Graph()
    for line in path.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            serial = int(line[6:11]); el = line[76:78].strip() or line[12:14].strip()[0]
            atoms[serial] = dict(name=line[12:16].strip(), res=line[17:20].strip(), chain=line[21],
                                 resseq=line[22:26].strip(), el=el)
            g.add_node(serial)
        elif line.startswith("CONECT"):
            fields = [line[i:i+5] for i in range(6, len(line), 5)]
            ids = [int(f) for f in fields if f.strip()]
            for j in ids[1:]:
                g.add_edge(ids[0], j)
    return atoms, g

def alkyl_carbons(g, atoms, start, ring):
    """Carbon count of a pure-hydrocarbon substituent (None if it has heteroatoms)."""
    seen, stack, n = set(ring) | {start}, [start], 0
    hetero = []
    while stack:
        a = stack.pop()
        if atoms[a]["el"] == "C": n += 1
        elif atoms[a]["el"] != "H": hetero.append(atoms[a]["el"])
        for b in g[a]:
            if b not in seen:
                seen.add(b)
                if atoms[b]["el"] != "H": stack.append(b)
                if len(seen) > 400: break
    return n, hetero

def first_hetero_path(g, atoms, start, ring, maxdepth=6):
    path, prev, cur = [atoms[start]["el"]], None, start
    return path

def audit(state, path):
    atoms, g = parse(path)
    heavy = [a for a in atoms if atoms[a]["el"] != "H"]
    pool = [a for a in heavy if atoms[a]["res"] in ("MC", "LYX")]
    sub = g.subgraph(heavy).copy()
    pool_g = sub.subgraph(pool).copy()
    rings = []
    for cyc in nx.cycle_basis(pool_g):
        if len(cyc) == 6 and all(atoms[a]["el"] == "C" for a in cyc):
            subs = []
            for a in cyc:
                for b in sub[a]:
                    if b in cyc: continue
                    el = atoms[b]["el"]
                    nC, het = alkyl_carbons(sub, atoms, b, cyc) if el == "C" else (None, None)
                    nb = sorted(atoms[x]["el"] for x in sub[b] if x not in cyc)
                    nH = sum(1 for x in g[b] if atoms[x]["el"] == "H")
                    kind = ("alkyl" if het == [] else
                            "CH2-N" if el == "C" and "N" in nb and nH == 2 else
                            "C(=O)O" if el == "C" and nb.count("O") == 2 else
                            "N" if el == "N" else f"C:{''.join(nb)}")
                    if het and "N" in nb and nH == 2: kind = "CH2-N"
                    # propyl-azide/triazole tail: chain of CH2 ending in N, not an arm CH2-NH-C(=O)
                    subs.append(dict(ring_atom=atoms[a]["name"], atom=atoms[b]["name"], serial=b, el=el,
                                     kind=kind, alkyl_carbons=nC if het == [] else None))
            hcount = sum(1 for a in cyc for x in g[a] if atoms[x]["el"] == "H")
            rings.append(dict(ring_atoms=[atoms[a]["name"] for a in cyc], ring_H=hcount,
                              residues=sorted({atoms[a]["res"] for a in cyc}), substituents=subs))
    hexasub = [r for r in rings if r["ring_H"] == 0]
    return dict(state=state, file=str(path.relative_to(ROOT)),
                sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                pool_heavy_counts=dict(Counter(atoms[a]["el"] for a in pool)),
                aromatic_carbocycles=len(rings), hexasubstituted_rings=hexasub,
                trisubstituted_aryl_rings=len([r for r in rings if r["ring_H"] == 3]))

def main():
    out = {"schema": "dia-health-one.nnc2215-mc-excess-carbon-audit.v1",
           "created_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
           "source_reference": {
               "figure": "Nature 2024 Extended Data Fig. 8i (41586_2024_8042_Fig13_ESM.jpg)",
               "url": "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-024-08042-3/MediaObjects/41586_2024_8042_Fig13_ESM.jpg",
               "expected_roof": "benzene with 3 x CH2NH (C15), 2 x CH3 (C18), 1 x (CH2)3N3 (C21-C23)",
               "expected_floor": "benzene with 3 x CH2NH (C5), 3 x CH2CH3 (C2-C1)"},
           "states": [audit(s, p) for s, p in FILES.items()]}
    summary = []
    for st in out["states"]:
        for r in st["hexasubstituted_rings"]:
            alk = sorted(s["alkyl_carbons"] for s in r["substituents"] if s["alkyl_carbons"] is not None)
            summary.append(dict(state=st["state"], ring=r["ring_atoms"], alkyl_substituent_carbons=alk,
                                non_alkyl=[s["kind"] for s in r["substituents"] if s["alkyl_carbons"] is None]))
    out["ring_summary"] = summary
    dest = ROOT / "data/nnc2215-mc-excess-carbon-audit.json"
    dest.write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(summary, indent=1))

if __name__ == "__main__":
    main()
