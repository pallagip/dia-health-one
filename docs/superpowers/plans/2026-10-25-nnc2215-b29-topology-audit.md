# NNC2215 B29 Linker–Roof Product Topology Audit Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a versioned atom-mapped graph for the NNC2215 B29 linker–roof click product and determine whether it exactly matches both corrected open and closed coordinate candidates by heavy-atom element and adjacency.

**Architecture:** Preserve the existing 83-atom roof precursor map, add the six numbered B29 linker carbons, its two oxygens, and the LysB29 NZ endpoint, then validate the source-derived product graph with RDKit. Parse the corrected candidate PDBs into heavy-atom graphs from CONECT/LINK records and use exact graph isomorphism in NetworkX. Transfer source bond labels to coordinates only through complete mappings; keep symmetry-equivalent mappings explicit. Record provenance, distances, hydrogen-neighbor differences, and unresolved assumptions in JSON.

**Tech Stack:** Python, RDKit, NetworkX, PDB fixed-column records, JSON, SHA-256.

---

### Task 1: Transcribe and validate the B29 linker–roof product graph

**Files:**
- Create: `research/structures/nnc2215-b29-linker-roof-product.smi`
- Create: `research/data/nnc2215-b29-linker-roof-product-map.json`
- Read: `research/structures/nnc2215-roof-propyl-azide.smi`
- Read: `research/data/nnc2215-roof-macrocycle-atom-map.json`
- Read: `research/data/nnc2215-linker-mapping-audit.json`
- Read: `research/data/nnc2215-chemistry-source-audit.json`

- [ ] Transcribe the B29 linker and click-product connectivity from Nature Extended Data Fig. 7a/Fig. 8i, their captions and the cited synthesis scheme; preserve roof graph map IDs 1–83.
- [ ] Assign new map IDs 84–89 to source linker carbons C3–C8, 90–91 to the carbonyl and linker oxygens, and 92 to the LysB29 NZ endpoint. Record source roles and any candidate PDB atom names in the legend without calling those names authoritative.
- [ ] Represent the carbamate, terminal alkyne conversion, and 1,4-disubstituted 1,2,3-triazole with a sanitized, explicit bond-order graph. Record the neutral-acid roof representation and the omitted LysB29 CE continuation as representation boundaries; do not infer simulation-pH states.
- [ ] If the primary product drawing does not uniquely establish the click regiochemistry, preserve both plausible mapped product graphs and mark unique topology assignment unresolved. Do not select a regioisomer from PDB adjacency or distances alone.
- [ ] Validate the source graph inventory as C62N16O14, 92 heavy atoms, 99 adjacency edges, one connected component, and cycle rank 8. Check unique complete map IDs, RDKit sanitization, valence, carbamate connectivity, and membership of the three N and two C atoms in the new triazole ring. Investigate and record any deviation instead of adjusting the expected values silently.

### Task 2: Implement the exact topology and coordinate audit

**Files:**
- Create: `research/audit_nnc2215_b29_topology.py`
- Create: `research/data/nnc2215-b29-topology-audit.json`
- Read: `research/data/nnc2215-carbon-reconciliation-candidate.json`
- Read: `research/data/nnc2215-roof-graph-coordinate-mapping-audit.json`
- Read: `research/data/nnc2215-linker-mapping-audit.json`
- Read: `research/structures/nnc2215-open-source-corrected-candidate.pdb`
- Read: `research/structures/nnc2215-closed-source-corrected-candidate.pdb`

- [ ] Load and hash-check the mapped product graph, map legend, reconciliation record, prior roof-map audit, and linker-map audit. Resolve the two corrected PDB paths and expected hashes from the reconciliation record.
- [ ] Parse PDB atoms, CONECT, and LINK records. Build the selected heavy-atom target from LYX B30, MC C1, and Lys B29 NZ; verify that each expected atom resolves once, the NZ endpoint is connected to the linker carbonyl, and no out-of-scope atoms enter the selected graph.
- [ ] Reproduce the roof-only preliminary result: 16 element/adjacency embeddings per state under the documented candidate tail anchors, and two roof automorphisms when source bond labels are preserved. Label these as precursor checks, not full-product assignments.
- [ ] Compare the complete 92-node source graph and each 92-node target graph by exact element-labeled graph isomorphism, requiring every source edge and every target edge to match. Retain every mapping and group mappings by source bond-label-preserving symmetry; do not use a subgraph match as a complete product match.
- [ ] Confirm the candidate graph inventory (C62N16O14), 99 edges, one component, and cycle rank 8 in each state. On every complete mapping, carry source bond orders to target atom pairs as a derived annotation and explicitly state that CONECT/LINK records carry no bond-order evidence.
- [ ] Calculate each proposed heavy-atom bond distance for every mapping in both states. Store the per-bond element/bond-order class, distance in Å, and broad plausibility interval with its cited basis; treat this as a geometry screen and retain outliers without using them alone to decide bond order or mapping.
- [ ] Record hydrogen-neighbor names/connectivity around mapped heavy atoms, open/closed differences, candidate-PDB hashes before and after processing, input graph hashes, software versions, ambiguity classes, and every unresolved assumption. Mark input coordinate files unchanged.
- [ ] If either state has no exact full-graph mapping, has unmatched atoms/edges, or supports multiple chemically distinct regioisomers, report the specific mismatch and stop short of claiming a unique validated topology.
- [ ] Run: `micromamba run -p /Users/andrii/scratch/dia-target-tools/env python research/audit_nnc2215_b29_topology.py`
- [ ] Expected: a reproducible JSON record covering both candidate states, including complete mappings if supported, symmetry classes, bond-distance observations, hydrogen-neighbor differences, hashes, and explicit limits; no PDB or simulation input is modified.

### Task 3: Update the chemistry report and research dashboard

**Files:**
- Modify: `research/NNC2215_CHEMISTRY_RECONCILIATION.md`
- Modify: `research/README.md`
- Modify: `research/index.html`

- [ ] Add the source-mapped B29 product graph and map legend links, the audit JSON and rerun command, and a concise summary of exact mapping, symmetry, graph checks, geometry observations, and open/closed hydrogen differences.
- [ ] State clearly whether the source drawing and both candidate graphs support a unique topology. If ambiguity or any mismatch remains, publish that unresolved result and do not describe the coordinates as validated chemistry.
- [ ] Preserve the existing provenance and limits: publisher-supplied model-built candidates, candidate methyl correction, adjacency-only PDB records, no validated protonation/charges/force field, and no target-binding MD or clinical conclusion.
- [ ] Keep the report and dashboard consistent with the JSON outcome; only add a positive topology-readiness claim if every acceptance condition is met.

### Task 4: Verify the audit artifacts and coordinate immutability

**Files:**
- Verify: `research/data/nnc2215-b29-topology-audit.json`
- Verify: `research/data/nnc2215-b29-linker-roof-product-map.json`
- Verify: both corrected candidate PDBs and all report/dashboard links

- [ ] Run: `python3 -m json.tool research/data/nnc2215-b29-topology-audit.json`
- [ ] Run: `python3 -m json.tool research/data/nnc2215-b29-linker-roof-product-map.json`
- [ ] Run the audit command from Task 2 again and confirm the PDB hashes recorded by the audit match the checked-in coordinate reconciliation record and the input PDB files remain byte-identical.
- [ ] Run: `git diff --check`
- [ ] Confirm every graph, legend, script, and JSON link target exists and that the final written status matches the audit evidence.

---

**Completion boundary:** A successful result supports a source-derived candidate graph for the B29 linker–roof component and its mapping to these two model-built coordinate candidates. It does not establish a complete insulin topology, validated force-field parameters, a production simulation model, pharmacology, efficacy, safety, or a drug product.
