# NNC2215 chemistry reconciliation

**Updated:** 25 October 2026
**Decision:** A source-supported coordinate correction candidate now matches the reported heavy-atom formula and passes serialized-coordinate checks. The full molecular graph and force-field chemistry remain unresolved; do not generate production parameters or start target MD.

## Evidence reviewed

The 2024 synthesis route says the B29 Lys ε-amino group is first modified with O-succinimidyl pentyn-1-oxycarbonyl at pH 10.5. A B1 O-peracetyl-D-glucoside active ester is coupled near pH 7.5 and then deacetylated. The B29 terminal alkyne is joined to the macrocycle roof propyl azide by copper-catalysed azide–alkyne cycloaddition. The resulting product is reported as C325H461N79O97S6. The article identifies the B29 Lys attachment and B1 Phe N-terminal glycoside attachment; the patent's Example 2 describes the same conjugate and reports the matching calculated mass-to-charge values. See the [2024 article](https://www.nature.com/articles/s41586-024-08042-3), its [Supporting Information](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-024-08042-3/MediaObjects/41586_2024_8042_MOESM1_ESM.pdf), and [WO2023144240A1 Example 2](https://patents.google.com/patent/WO2023144240A1/en).

The coordinate inputs are the author-supplied, model-built open and closed structures by Eva Johansson and Daniele Granata, [“3D models of NNC2215”](https://springernature.figshare.com/articles/dataset/3D_models_of_NNC2215/26526460), Figshare DOI 10.6084/m9.figshare.26526460, licensed CC BY 4.0. The PDBs linked in this report are modified coordinate candidates derived from those files; they are not the original author files or experimentally solved structures.

The numbered schemes identify six B29 linker carbons (source C3–C8) and the three-carbon roof propyl-azide tail (source C21–C23). Formula arithmetic in [the existing audit](data/nnc2215-model-formula-audit.json) gives a desB30 insulin parent of C253H376N64O75S6 and an attachment difference of C72H85N15O22. The B1 attachment accounts for C10O8 heavy atoms. The B29 component consists of the reported C56H63N15O12 macrocycle roof propyl azide plus the C6H6O2 B29 linker after the NHS leaving group departs, giving C62N15O14 heavy atoms after conjugation.

The publisher's model PDBs have a 113-node/120-edge heavy-atom adjacency match for the selected PHX/LYX/MC and attachment scope, after documented open/closed LYX name aliases. [The connectivity audit](data/nnc2215-connectivity-audit.json) confirms that these records give adjacency only; PDB CONECT/LINK data do not encode bond order or formal charge.

## Two-carbon discrepancy

Both source model states contain 64 carbon, 15 nitrogen and 14 oxygen records in the LYX/MC B29-linked pool, against 62 carbon, 15 nitrogen and 14 oxygen expected from the synthesis-derived formula. The two-carbon surplus is present in both open and closed models.

The [candidate linker map](data/nnc2215-linker-mapping-audit.json) follows the six-carbon B29 linker and the three-carbon roof tail through serial-level PDB connectivity. Under that candidate partition, two tail carbons fall in LYX and one falls in MC, leaving 54 source macrocycle carbons expected in MC. The PDB MC residue has 56 carbon records.

The [coordinate-reconciliation candidate](data/nnc2215-carbon-reconciliation-candidate.json) removes MC C27/C29 and their four hydrogens, then connects the retained methyl carbons C28/C30 directly to aryl carbons C73/C75. This matches the published direct aryl methyl substituents and the SI's six-proton H-18 assignment. The B29-linked pool becomes C62N15O14 in both states, and the full heavy-atom inventory becomes C325N79O97S6. This is a source-supported mapping inference, not an authoritative atom map.

The candidate PDBs reopen successfully with all connectivity references resolved, no duplicate atom names within residues, strictly increasing atom serials, and records no longer than 80 columns. Both new aryl-methyl bonds are 1.50 Å and their C-H bonds are 1.09 Å; the local methyl-hydrogen overlap screen passes. In the open model, ASN A21's misplaced HXT record was removed and an OXT was reconstructed in the CA-C-O plane (C-OXT 1.25 Å; O-C-OXT 119.66°). These checks validate file consistency and local geometry only. They do not validate the full graph, stereochemistry, protonation, or force-field energy.

The synthesis schemes and patent establish the intended product, linker chemistry, reaction order and product composition. They do not provide a complete source-atom-to-PDB mapping for these model residues. Therefore the discrepancy is localized but not reconciled.

## Provisional chemical assignments

These are source-informed expectations for a future graph review, not final force-field assignments:

- The B29 Lys ε-amino group is converted to a neutral carbamate in the pentynyl linker. The alkyne is consumed in CuAAC; the expected product is the 1,4-disubstituted 1,2,3-triazole.
- The B1 Phe α-amino group is acylated by the glucoside linker and is neutral in that conjugated form. The glucoside is the deprotected D-glucoside specified by the synthesis scheme; its anomer and every mapped stereocentre must be explicitly verified from the source drawing, rather than inferred from PDB coordinates.
- The macrocycle has three reported carboxylic-acid groups. Their likely dominant state near pH 7.4 is three carboxylates, but the complete charge assignment must be computed from the resolved graph and checked against pKa/microstate assumptions.
- The existing isolated PaG audit supports the singly deprotonated glucuronate, net charge −1, and five CIP labels S, S, R, S, S. It validates the isolated ligand files only; it does not validate the NNC2215 conjugate.
- The full insulin system requires an explicit pH 7.4 protonation inventory, including terminal groups, ionizable side chains, disulfides, and the conjugation sites. A pH label alone does not define that inventory.

The current PDB records also disagree around hydrogen ownership: the closed model assigns a nearby PheB1 N hydrogen to PHX without a CONECT edge, and the LYX CD2 hydrogen-neighbour count differs between states. These are additional unresolved topology records.

## Evidence required to close chemistry

1. A complete atom-mapped graph for B1, the NNC2215 B29 linker, macrocycle, triazole and insulin attachment sites, independently checked against the synthesis schemes and patent.
2. Independent chemistry review of the candidate MC C27/C29 removal and C28/C30 aryl-methyl mapping, plus a one-to-one graph-to-coordinate mapping for both states that explains the differing hydrogen/link records.
3. Explicit bond orders, aromaticity, CuAAC regioisomer, stereochemistry, protonation microstates, formal-charge ledger and integer total charge for each simulated condition.
4. A versioned structure file with the validated graph, atom names and source links. Compare the mapped source graph against the reported formula and HRMS.
5. Only after items 1–4 pass: select one force-field family and validate custom macrocycle/linker parameters, carbohydrate compatibility, PaG, torsions/conformers, charge sums and prepared-system energies.

No production topology or parameter files were generated in this review. Rebuild the candidate coordinates with `python research/reconcile_nnc2215_source_structure.py`. The readiness evidence remains in [the candidate coordinate audit](data/nnc2215-carbon-reconciliation-candidate.json), [model validation record](data/nnc2215-model-validation.json), [formula audit](data/nnc2215-model-formula-audit.json), [connectivity audit](data/nnc2215-connectivity-audit.json), and [linker-map audit](data/nnc2215-linker-mapping-audit.json).
