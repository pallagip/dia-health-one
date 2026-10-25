# NNC2215 B1 Glucoside Topology Audit Design

**Date:** 2026-10-25

**Status:** User approved this design; review this spec before implementation.

## Goal

Build a source-mapped graph for the deprotected B1 glucoside attachment in NNC2215 and determine whether it matches the B1 modification in both corrected, model-built coordinate candidates. Resolve source-supported stereochemistry and distinguish a complete graph match from the PDB's adjacency-only records. This closes the B1 chemistry gate before integrating the B1 and B29 modifications with the full insulin graph.

## Current evidence

- The [2024 Nature article](https://www.nature.com/articles/s41586-024-08042-3) describes an O1-glucoside with a short linker attached to the B-chain N-terminal amino group, PheB1. LC–MS analysis places the glucoside on the B1–B22 fragment.
- The [Supporting Information](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-024-08042-3/MediaObjects/41586_2024_8042_MOESM1_ESM.pdf) documents a beta-carboxymethyl-oxyethyl-O-peracetyl-D-glucoside active ester, coupling to insulin, and deprotection. Its anomeric H1 signal is reported as a doublet with a 7.9 Hz coupling constant, consistent with a beta-D-glucopyranoside.
- The source-supported product path, from the protein outward, is PheB1 N–C(=O)–CH2–O–CH2–CH2–O–C1 of beta-D-glucopyranose. The sugar's C2, C3, C4, and C6 hydroxyls are deprotected in the final conjugate. The net B1 attachment increment is C10H16O8 after amide formation; this formula increment excludes the protein's PheB1 nitrogen.
- The formula audit derives a C10O8 heavy-atom inventory for the B1 substituent, matching the PHX residue inventory in both models. Each corrected candidate contains a PDB LINK from PHX C2 to PheB1 N. These counts and LINK records establish the comparison target, not a source-to-coordinate atom map or bond orders.
- The open and closed coordinates are author-supplied model-built candidates, not experimentally resolved NNC2215 structures. The closed candidate assigns a nearby PheB1 N hydrogen to PHX without a CONECT edge; hydrogen ownership must remain an explicit audit item.
- No source-mapped B1 graph or B1-specific coordinate audit is currently present. The [patent WO2023144240A1](https://patents.google.com/patent/WO2023144240A1/en) independently describes glycoside attachment through a linker to the B1 alpha-amino group.

## Proposed design

1. Transcribe the deprotected B1 product graph from the published synthesis scheme and supporting information. Use stable graph map IDs and retain source numbering where the drawing supports it. Encode the beta-D anomer and all glucose stereocentres explicitly. Record the transformation from the peracetylated active ester to the deprotected acid and then to the PheB1 amide.
2. Represent the graph boundary explicitly. Include the PheB1 nitrogen as the protein attachment atom in the product graph, while treating the N–Cα peptide continuation as out-of-scope context. Record the net substituent increment formula separately from the graph formula that includes the protein atom.
3. Validate graph atom counts, formula arithmetic, connectedness, cycle rank, valence, and stereochemistry with RDKit. Cross-check the source-derived B1 increment and the reported B1–B22 fragment mass without treating the PDB hydrogen inventory as formula evidence.
4. Parse both corrected PDB candidates into heavy-atom adjacency graphs using their CONECT/LINK records. Select the PHX residue plus PheB1 N. Require exact element-labeled graph isomorphism for the complete B1 product component; report all symmetry-equivalent mappings and do not choose among them by distance alone. Transfer source bond labels only as a mapping-derived annotation because PDB CONECT/LINK records do not encode bond order.
5. Record per-bond coordinate distances as a geometry screen, hydrogen-neighbor records around PHX and PheB1 N, open/closed differences, hashes, graph counts, mapping ambiguity, source provenance, and unresolved assumptions. Preserve the PDB inputs byte-for-byte. If either model fails a complete graph match or the source leaves a stereochemical assignment unresolved, report the specific issue and stop short of claiming a validated topology.
6. After the audit has a defensible result, update the chemistry reconciliation report, README, and research dashboard. Keep the result scoped to the B1 component and state that full-conjugate topology, protonation, charges, force-field parameters, and simulation readiness remain later gates.

## Deliverables

- A mapped product graph artifact at `research/structures/nnc2215-b1-glucoside-attachment.smi` and a map/provenance legend at `research/data/nnc2215-b1-glucoside-atom-map.json`.
- A reproducible audit script at `research/audit_nnc2215_b1_topology.py` and its JSON record at `research/data/nnc2215-b1-topology-audit.json`.
- A concise reconciliation, README, and dashboard update that reports the result and limits.

## Acceptance criteria

- The source graph explicitly represents the beta-D glucopyranoside stereochemistry, deprotected linker, and PheB1 N-amide connection, with every atom-map ID unique and source provenance recorded.
- The net substituent increment is C10H16O8, comprising C10O8 heavy atoms; including the PheB1 N anchor adds that protein atom to the audit graph without changing the stated increment.
- RDKit sanitization, valence, formula, and stereochemistry checks pass for the documented representation. The reported B1–B22 LC–MS fragment is used as an independent composition check with digestion/protonation assumptions disclosed.
- Each candidate is compared on all source graph nodes and edges. Exact match status, every mapping, symmetry classes, mismatches, and hydrogen-neighbor differences are reported separately for open and closed models.
- Bond distances are not used to assign bond orders, stereochemistry, or a unique mapping. Input model files remain byte-identical.
- The dashboard and report distinguish a source-supported graph-to-candidate mapping from a validated experimental structure or simulation-ready model.

## Alternatives considered

- **Build and audit the B1 graph against PHX and PheB1 N (selected):** closes the missing, bounded B1 chemistry gate and checks it against both existing candidates before full-conjugate integration.
- **Build the complete insulin-conjugate graph immediately:** would include the B1 and B29 chemistry, the entire protein, disulfides, and current coordinate ambiguities in one large audit. It obscures the B1-specific result and increases the number of coupled assumptions.
- **Prepare force-field inputs now:** would require assigning chemistry to PHX before its atom mapping and stereochemistry are explicit, so resulting energies or dynamics would not be interpretable.

## Out of scope

- Editing or minimizing the published model coordinates.
- Assigning physiological protonation microstates, partial charges, or force-field atom types.
- Claiming the open or closed candidate is an experimentally resolved NNC2215 structure.
- Completing full insulin/B1/B29 conjugate topology or declaring production molecular-dynamics readiness in this B1 audit.
- Inferring efficacy, safety, or clinical benefit from the topology or coordinate comparison.

## Limits

A successful result supports the source-derived B1 attachment graph and its correspondence to the two supplied model-built candidates. It does not validate a unique full-conjugate coordinate map, the B29 scaffold mapping, the protein's simulated protonation, custom parameters, receptor binding, pharmacology, or clinical performance. The next chemistry gate is to merge the resolved B1 component with the B29 component and full insulin graph, then audit the resulting complete conjugate before force-field preparation.
