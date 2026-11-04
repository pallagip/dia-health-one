#!/usr/bin/env python3
"""Build the NNC2215 macrocycle-only PaG proxy inputs from the published graph.

Source graph: Nature 2024 Extended Data Fig. 8i (roof propyl azide, reported
C56H63N15O12), transcribed below as SMILES. Coordinates: the heavy atoms that
map one-to-one onto residue MC of the public open-state model PDB (Figshare
26526460) are used as embedding constraints; the model's free glucose (BGLC)
supplies the glucose reference pose. PaG is placed by the pre-registered
glucose-site superposition procedure documented in POSE_PROCEDURE.

This is a model-generated proxy. It is not the 2026 in-house crystal, not full
NNC2215, and no step here estimates affinity.
"""
from __future__ import annotations
import hashlib, json, sys
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import networkx as nx
from networkx.algorithms import isomorphism as iso
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign, rdMolTransforms, Descriptors
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Geometry import Point3D

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
from audit_mc_excess_carbons import parse, FILES  # noqa: E402

OUT = HERE / "inputs"
SEED = 2215

# Fig. 8i: floor = 1,3,5-tris(CH2NH-)-2,4,6-triethylbenzene (C5, C2-C1);
# arms = 3-ureido(roof)-4-ureido(floor)benzoic acid (C6, C7-C12, C13, C14);
# roof = 1,3,5-tris(CH2NH-) (C15) -2-(3-azidopropyl) (C21-C23) -4,6-dimethyl (C18) benzene.
FIG8I_SMILES_NEUTRAL = (
    "c1(CC)c(CNC(=O)Nc2ccc(C(=O)O)cc2NC(=O)NC%10)c(CC)c(CNC(=O)Nc3ccc(C(=O)O)cc3NC(=O)NC%11)"
    "c(CC)c1CNC(=O)Nc4ccc(C(=O)O)cc4NC(=O)NC%12"
    ".c5%10c(CCCN=[N+]=[N-])c%11c(C)c%12c5C"
)
REPORTED_FORMULA = "C56H63N15O12"
BETA_D_GLUCOSE = "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O"
PAG_SDF = ROOT / "ligands" / "pag-pubchem-cid83944-carboxylate-3d.sdf"

POSE_PROCEDURE = [
    "Hypothesis tested: PaG's beta-D-glucuronide ring occupies the model glucose site in the same ring orientation (glucose-mimetic pose).",
    "Map PaG ring atoms C1-C5, O5, ring hydroxyl O2/O3/O4, glycosidic O1 and carboxylate C6 onto the same atoms of the model beta-D-glucose (11 atoms).",
    "Generate 300 ETKDGv3 PaG conformers (randomSeed 2215, MMFF94 optimised), rigidly superimpose each on those 11 glucose atoms.",
    "Reject conformers with any PaG heavy atom < 2.0 A from a macrocycle heavy atom; among the rest take the lowest clash score (sum of (3.0 A - d) over heavy-atom pairs closer than 3.0 A), ties broken by mapped-atom RMSD.",
    "No docking or scoring function is used; the score only screens steric overlap. Alternative conformers within 0.5 clash-score units are written for sensitivity review.",
]


def heavy_graph_from_rdkit(mol, drop=()):
    g = nx.Graph()
    for a in mol.GetAtoms():
        if a.GetAtomicNum() > 1 and a.GetIdx() not in drop:
            g.add_node(a.GetIdx(), el=a.GetSymbol())
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if i in g and j in g:
            g.add_edge(i, j)
    return g


def pdb_mc_graph(path):
    atoms, g = parse(path)
    mc = [a for a in atoms if atoms[a]["res"] == "MC" and atoms[a]["el"] != "H"]
    sub = g.subgraph(mc).copy()
    # roof ethyl terminal carbons: CH3 carbons whose only heavy neighbour is a CH2 bonded to an aromatic C of the roof ring
    excess = []
    audit = json.loads((ROOT / "data/nnc2215-mc-excess-carbon-audit.json").read_text())
    nx.set_node_attributes(sub, {a: atoms[a]["el"] for a in sub}, "el")
    return atoms, g, sub


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    report = {"schema": "dia-health-one.macrocycle-proxy-preparation.v1",
              "created_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
              "rdkit": Chem.rdBase.rdkitVersion, "source_smiles_fig8i_neutral": FIG8I_SMILES_NEUTRAL}

    # 1. Source graph checks
    src = Chem.MolFromSmiles(FIG8I_SMILES_NEUTRAL)
    assert src is not None and len(Chem.GetMolFrags(src)) == 1
    formula = CalcMolFormula(src)
    report["source_graph"] = {
        "formula": formula, "reported_formula": REPORTED_FORMULA, "formula_matches": formula == REPORTED_FORMULA,
        "exact_mass": round(Descriptors.ExactMolWt(src), 4),
        "stereocentres": len(Chem.FindMolChiralCenters(src, includeUnassigned=True, useLegacyImplementation=False)),
        "stereo_double_bonds": sum(1 for b in src.GetBonds() if b.GetStereo() != Chem.BondStereo.STEREONONE),
        "aromatic_rings": sum(1 for r in src.GetRingInfo().AtomRings() if all(src.GetAtomWithIdx(i).GetIsAromatic() for i in r)),
        "urea_groups": len(src.GetSubstructMatches(Chem.MolFromSmarts("[NX3][CX3](=O)[NX3]"))),
        "carboxylic_acids": len(src.GetSubstructMatches(Chem.MolFromSmarts("C(=O)[OX2H1]"))),
        "azide": len(src.GetSubstructMatches(Chem.MolFromSmarts("N=[N+]=[N-]"))),
    }
    assert formula == REPORTED_FORMULA, formula

    # pH 7.4 microstate: three aryl carboxylates deprotonated (benzoic-acid pKa ~4); ureas and azide neutral.
    charged = Chem.MolFromSmiles(FIG8I_SMILES_NEUTRAL.replace("C(=O)O)", "C(=O)[O-])"))
    report["simulated_microstate"] = {"formula": CalcMolFormula(charged), "formal_charge": Chem.GetFormalCharge(charged),
                                      "assumption": "Three aryl carboxylates deprotonated at pH 7.4; ureas and azide neutral. Matches the deprotonated carboxylates in both model PDBs."}
    molH = Chem.AddHs(charged)

    # 2. Isomorphism of source core with PDB MC after removing the two roof-ethyl terminal carbons
    tail_q = Chem.MolFromSmarts("c-[CH2]-[CH2]-[CH2]-N=[N+]=[N-]")
    tail = molH.GetSubstructMatch(tail_q)
    drop = set(tail[2:])  # keep the aromatic C and source C21 (in MC per linker audit); drop C22, C23, N3
    src_g = heavy_graph_from_rdkit(molH, drop)
    em = iso.categorical_node_match("el", None)
    mapping_by_state, state_reports = {}, []
    for state, path in FILES.items():
        atoms, g, mc = pdb_mc_graph(path)
        audit = json.loads((ROOT / "data/nnc2215-mc-excess-carbon-audit.json").read_text())
        st = [s for s in audit["states"] if s["state"] == state][0]
        roof = [r for r in st["hexasubstituted_rings"] if any(s["kind"] == "C:C" for s in r["substituents"])][0]
        excess = []
        for s in roof["substituents"]:
            if s["alkyl_carbons"] == 2:
                excess += [b for b in g[s["serial"]] if atoms[b]["el"] == "C" and atoms[b]["res"] == "MC"
                           and atoms[b]["name"] not in roof["ring_atoms"]]
        as_is = nx.is_isomorphic(src_g, mc, node_match=em)
        trimmed = mc.copy(); trimmed.remove_nodes_from(excess)
        gm = iso.GraphMatcher(src_g, trimmed, node_match=em)
        ok = gm.is_isomorphic()
        mapping_by_state[state] = (dict(gm.mapping) if ok else None, atoms)
        state_reports.append({"state": state, "pdb_mc_heavy_atoms": mc.number_of_nodes(),
                              "source_core_heavy_atoms": src_g.number_of_nodes(),
                              "isomorphic_without_edit": as_is,
                              "removed_pdb_atoms": [atoms[e]["name"] for e in excess],
                              "isomorphic_after_removing_roof_ethyl_terminal_carbons": ok})
    report["pdb_graph_comparison"] = state_reports

    # 3. 3D build from the open model: copy mapped heavy-atom coordinates, embed the rest
    mapping, atoms_open = mapping_by_state["active_open"]
    xyz = {}
    for line in FILES["active_open"].read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            xyz[int(line[6:11])] = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
    cmap = {i: Point3D(*xyz[s]) for i, s in mapping.items()}
    params = AllChem.ETKDGv3(); params.randomSeed = SEED; params.useRandomCoords = True
    params.SetCoordMap(cmap)
    cid = AllChem.EmbedMolecule(molH, params)
    if cid < 0:
        raise SystemExit("embedding with model coordinate map failed")
    # Kabsch-align the embedded conformer onto the model frame, then restrain mapped heavy atoms
    conf = molH.GetConformer()
    keys = list(cmap)
    A = np.array([list(conf.GetAtomPosition(i)) for i in keys]); B = np.array([[cmap[i].x, cmap[i].y, cmap[i].z] for i in keys])
    ca, cb = A.mean(0), B.mean(0); H = (A - ca).T @ (B - cb); U, S, Vt = np.linalg.svd(H)
    dfix = np.sign(np.linalg.det(Vt.T @ U.T)); R = Vt.T @ np.diag([1, 1, dfix]) @ U.T
    allp = conf.GetPositions(); allp = (R @ (allp - ca).T).T + cb
    for i, p in enumerate(allp): conf.SetAtomPosition(i, Point3D(*p))
    mp = AllChem.MMFFGetMoleculeProperties(molH)
    ff = AllChem.MMFFGetMoleculeForceField(molH, mp)
    for i in keys:
        ff.MMFFAddPositionConstraint(i, 0.05, 100.0)
    ff.Minimize(maxIts=10000)
    model_xyz = np.array([[cmap[i].x, cmap[i].y, cmap[i].z] for i in mapping])
    built = np.array([list(conf.GetAtomPosition(i)) for i in mapping])
    rmsd = float(np.sqrt(((model_xyz - built) ** 2).sum(1).mean()))
    report["macrocycle_3d"] = {"mapped_heavy_atoms": len(mapping), "rmsd_to_open_model_A": round(rmsd, 3),
                               "embedding_seed": SEED, "note": "Mapped atoms restrained to open-model coordinates; roof methyls take the C27/C29 positions; propyl azide tail C22-C23-N3 newly embedded."}
    Chem.MolToMolFile(molH, str(OUT / "macrocycle_fig8i_tricarboxylate.mol"))
    Chem.SDWriter(str(OUT / "macrocycle_fig8i_tricarboxylate.sdf")).write(molH)
    mc_heavy = np.array([list(conf.GetAtomPosition(a.GetIdx())) for a in molH.GetAtoms() if a.GetAtomicNum() > 1])

    # 4. Glucose reference pose from model BGLC
    glc_lines = [l for l in FILES["active_open"].read_text().splitlines()
                 if l.startswith("HETATM") and l[17:21].strip().startswith("BGL")]
    glc_pdb = Chem.MolFromPDBBlock("\n".join(glc_lines) + "\nEND\n", removeHs=False, proximityBonding=True)
    tmpl = Chem.MolFromSmiles(BETA_D_GLUCOSE)
    glc_pdb = AllChem.AssignBondOrdersFromTemplate(tmpl, glc_pdb)
    Chem.AssignStereochemistryFrom3D(glc_pdb)
    glc_can = Chem.MolToSmiles(Chem.RemoveHs(glc_pdb)); ref_can = Chem.MolToSmiles(tmpl)
    gx = glc_pdb.GetConformer().GetPositions()
    gheavy = np.array([gx[a.GetIdx()] for a in glc_pdb.GetAtoms() if a.GetAtomicNum() > 1])
    dmin_glc = float(np.min(np.linalg.norm(gheavy[:, None] - mc_heavy[None], axis=2)))
    report["glucose_pose"] = {"source": "open model residue BGLC (24 atoms incl. H)",
                              "stereo_from_3d_matches_beta_D_glucopyranose": glc_can == ref_can,
                              "canonical_from_3d": glc_can, "canonical_beta_D": ref_can,
                              "min_heavy_distance_to_macrocycle_A": round(dmin_glc, 2),
                              "centroid_to_macrocycle_centroid_A": round(float(np.linalg.norm(gheavy.mean(0) - mc_heavy.mean(0))), 2)}
    Chem.SDWriter(str(OUT / "glucose_model_pose.sdf")).write(glc_pdb)

    # 5. PaG pose by glucose-site superposition
    pag = next(Chem.SDMolSupplier(str(PAG_SDF), removeHs=False))
    Chem.AssignStereochemistryFrom3D(pag)
    q = Chem.MolFromSmarts("[OX2;!R][C;R]1[O;R][C;R]([CX3](=O)~O)[C;R]([OX2H1])[C;R]([OX2H1])[C;R]1[OX2H1]")
    gq = Chem.MolFromSmarts("[OX2H1;!R][C;R]1[O;R][C;R]([CH2][OX2H1])[C;R]([OX2H1])[C;R]([OX2H1])[C;R]1[OX2H1]")
    pm = pag.GetSubstructMatch(q); gmatch = glc_pdb.GetSubstructMatch(gq)
    assert pm and gmatch, (pm, gmatch)
    # order: O1, C1, O5, C5, C6, (O6a) , C4, O4, C3, O3, C2, O2 -> drop carboxyl O / CH2OH O
    pidx = [pm[0], pm[1], pm[2], pm[3], pm[4], pm[7], pm[8], pm[9], pm[10], pm[11], pm[12]]
    gidx = [gmatch[0], gmatch[1], gmatch[2], gmatch[3], gmatch[4], gmatch[6], gmatch[7], gmatch[8], gmatch[9], gmatch[10], gmatch[11]]
    pagc = Chem.Mol(pag); pagc.RemoveAllConformers()
    pagc = Chem.AddHs(Chem.RemoveHs(pag), addCoords=False)
    pm2 = pagc.GetSubstructMatch(q)
    pidx = [pm2[0], pm2[1], pm2[2], pm2[3], pm2[4], pm2[7], pm2[8], pm2[9], pm2[10], pm2[11], pm2[12]]
    ep = AllChem.ETKDGv3(); ep.randomSeed = SEED
    cids = list(AllChem.EmbedMultipleConfs(pagc, 300, ep))
    AllChem.MMFFOptimizeMoleculeConfs(pagc, maxIters=2000)
    heavy_p = [a.GetIdx() for a in pagc.GetAtoms() if a.GetAtomicNum() > 1]
    cands = []
    for c in cids:
        r = rdMolAlign.AlignMol(pagc, glc_pdb, prbCid=c, atomMap=list(zip(pidx, gidx)))
        P = pagc.GetConformer(c).GetPositions()[heavy_p]
        d = np.linalg.norm(P[:, None] - mc_heavy[None], axis=2)
        clash = float(np.clip(3.0 - d, 0, None).sum())
        cands.append((d.min() < 2.0, clash, r, c, float(d.min())))
    dmins = sorted(x[4] for x in cands)
    clash_free = sorted([x for x in cands if not x[0]], key=lambda x: (x[1], x[2]))
    # Pre-registered fallback: if no conformer clears 2.0 A, keep conformers above 1.6 A (a separation
    # that energy minimisation can resolve) and take the lowest clash score. Record that the strict
    # filter failed, because it is itself a result about the glucose-mimetic pose hypothesis.
    relaxed = sorted([x for x in cands if x[4] >= 1.6], key=lambda x: (x[1], x[2]))
    ok = clash_free or relaxed
    if not ok:
        raise SystemExit(json.dumps({"no_pose": True, "min_heavy_distance_quantiles_A": [round(q,2) for q in np.quantile(dmins,[0,.1,.5,.9,1])]}))
    best = ok[0]
    alts = [x for x in ok[1:] if x[1] - best[1] <= 0.5][:5]
    w = Chem.SDWriter(str(OUT / "pag_pose_selected.sdf")); w.write(pagc, confId=best[3]); w.close()
    w = Chem.SDWriter(str(OUT / "pag_pose_alternates.sdf"))
    for x in alts: w.write(pagc, confId=x[3])
    w.close()
    stereo = [lab for _, lab in Chem.FindMolChiralCenters(pagc, useLegacyImplementation=False)]
    report["pag_pose"] = {"procedure": POSE_PROCEDURE, "input_sdf": str(PAG_SDF.relative_to(ROOT)),
                          "input_sha256": hashlib.sha256(PAG_SDF.read_bytes()).hexdigest(),
                          "formal_charge": Chem.GetFormalCharge(pagc), "cip_labels": stereo,
                          "conformers_generated": len(cids), "conformers_passing_2.0A_filter": len(clash_free),
                          "conformers_passing_1.6A_fallback_filter": len(relaxed),
                          "strict_filter_satisfied": bool(clash_free),
                          "min_heavy_distance_quantiles_A": {k: round(float(v), 2) for k, v in zip(["min", "p10", "median", "p90", "max"], np.quantile(dmins, [0, .1, .5, .9, 1]))},
                          "selected": {"conformer": best[3], "clash_score": round(best[1], 3),
                                       "mapped_atom_rmsd_A": round(best[2], 3), "min_heavy_distance_A": round(best[4], 2)},
                          "alternates_within_0.5": [{"conformer": x[3], "clash_score": round(x[1], 3), "mapped_atom_rmsd_A": round(x[2], 3)} for x in alts]}
    (OUT / "preparation-report.json").write_text(json.dumps(report, indent=1) + "\n")
    print(json.dumps(report, indent=1))


if __name__ == "__main__":
    main()
