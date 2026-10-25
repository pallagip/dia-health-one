# NNC2215 B1 Glucoside Topology Audit Implementation Plan

> **For agentic workers:** Use the already-approved parallel preference for independent tasks. Steps use checkbox (`- [ ]`) syntax. Do not add a software test suite; verify the chemistry with the reproducible audit and artifact checks below.

**Goal:** Build a source-mapped B1 glucoside attachment graph and audit its exact element/adjacency correspondence to both corrected NNC2215 coordinate candidates.

**Architecture:** One task transcribes the source graph and map legend. In parallel, a second task independently inventories the PDB attachment target and hydrogen records without editing shared files. The primary agent integrates both results in a reproducible RDKit/NetworkX audit, then updates the report and research dashboard from the JSON result.

**Tech Stack:** Python 3, RDKit, NetworkX, fixed-column PDB records, JSON, SHA-256, the existing `dia-target-tools` micromamba environment.

---

## File structure

- Create `research/structures/nnc2215-b1-glucoside-attachment.smi`: source-numbered, atom-mapped deprotected B1 product graph, including a mapped attachment placeholder for PheB1 N.
- Create `research/data/nnc2215-b1-glucoside-atom-map.json`: graph schema, atom roles/source labels, attachment-placeholder mapping, source URLs, transformations, and representation limits.
- Create `research/audit_nnc2215_b1_topology.py`: reproducible source-graph, candidate-PDB, exact-mapping, geometry, and hydrogen-record audit.
- Create `research/data/nnc2215-b1-topology-audit.json`: run output with input hashes, validation findings, all complete mappings, per-state differences, and limits.
- Modify `research/NNC2215_CHEMISTRY_RECONCILIATION.md`, `research/README.md`, and `research/index.html`: concise result summary and links to the graph, legend, audit, and rerun command.
- Read `research/data/nnc2215-model-formula-audit.json`, `research/data/nnc2215-connectivity-audit.json`, `research/data/nnc2215-carbon-reconciliation-candidate.json`, and both corrected candidate PDBs. Preserve all source inputs.

## Parallel work

Run Tasks 1 and 2 at the same time in separate subagents.

- Task 1 owns the two graph artifacts only. It must use the published SI drawing/numbering and report its atom numbering and stereochemical interpretation.
- Task 2 is read-only. It independently reports the selected PHX/PheB1-N target atoms and edges in both PDB candidates, the Phe N–PHX LINK details, and the hydrogen records/ownership difference. It must not edit files or decide bond orders from geometry.
- Task 3 starts after both tasks return. It integrates their findings in the audit script and resolves any discrepancy explicitly.

### Task 1: Transcribe the source-supported B1 product graph

**Files:**
- Create: `research/structures/nnc2215-b1-glucoside-attachment.smi`
- Create: `research/data/nnc2215-b1-glucoside-atom-map.json`
- Read: `research/data/nnc2215-chemistry-source-audit.json`
- Read: `research/data/nnc2215-model-formula-audit.json`
- Read: `research/NNC2215_CHEMISTRY_RECONCILIATION.md`
- Source: Nature Fig. 2a, Extended Data Fig. 7b, Supporting Information §§2 and 4, and patent WO2023144240A1 Example 2.

- [ ] Confirm source numbering and atom roles from Extended Data Fig. 7b and the B1 building-block scheme. Record each numbered carbon and each unnumbered oxygen/nitrogen endpoint in the JSON legend.
- [ ] Encode the final deprotected product path as PheB1 N–C(=O)–CH2–O–CH2–CH2–O–C1 of beta-D-glucopyranose. Keep the PheB1 N attachment endpoint as a mapped placeholder that the audit legend maps to the PDB nitrogen; do not count this endpoint as part of the C10H16O8 substituent increment.
- [ ] Encode all five glucose stereocentres explicitly. Preserve the four free hydroxyls (C2, C3, C4, C6) and document that the source's peracetylated intermediate is deprotected before the final structure.
- [ ] Store the mapped SMILES and a JSON legend containing source label, element, role, stereochemistry, model-match role, source URL/location, and transformation for each atom. Include graph and input SHA-256 values.
- [ ] Confirm the net substituent increment C10H16O8 and heavy-atom inventory C10O8 independently from the formula audit. Record that the audit graph adds one PheB1 nitrogen endpoint for coordinate matching.
- [ ] Record the published B1–B22 fragment ion benchmark (observed m/z 1376.6, calculated m/z 1376.2 for [M+2H]2+). Derive the neutral fragment from sequence `FVNQHLCGSHLVEALYLVCGER`, free N and C termini after trypsin cleavage, reduced cysteines after TCEP treatment, and the C10H16O8 B1 increment; disclose the protonation convention. Do not use PDB hydrogen-record counts as a formula check.

### Task 2: Independently inventory the candidate B1 attachment target

**Files:**
- Read only: `research/structures/nnc2215-open-source-corrected-candidate.pdb`
- Read only: `research/structures/nnc2215-closed-source-corrected-candidate.pdb`
- Read only: `research/data/nnc2215-connectivity-audit.json`
- Read only: `research/data/nnc2215-model-formula-audit.json`

- [ ] For each state, select every heavy atom in PHX chain B residue 0 and the single PheB1 N atom in chain B residue 1. Report atom serial, name, element, and coordinate for each selected node.
- [ ] Resolve the PheB1 N–PHX C2 LINK and every within-target CONECT edge. Report target node and edge counts and any unresolved endpoint; treat these records as adjacency only.
- [ ] Inventory hydrogen records on PHX and PheB1 N, distinguishing record ownership from explicit CONECT adjacency. Report the closed-state nearby PheB1 N hydrogen that is assigned to PHX without a CONECT edge.
- [ ] Return a compact mapping-preflight note to the primary agent. Do not modify PDB files, assign stereochemistry, infer bond order, or create an atom map from distances.

### Task 3: Implement the reproducible exact B1 topology audit

**Files:**
- Create: `research/audit_nnc2215_b1_topology.py`
- Create: `research/data/nnc2215-b1-topology-audit.json`
- Read: both Task 1 graph artifacts and Task 2 preflight findings
- Read: `research/data/nnc2215-carbon-reconciliation-candidate.json`
- Read: `research/data/structure-manifest.json`
- Read: both corrected PDB candidates

- [ ] Follow `research/audit_nnc2215_b29_topology.py` for repository conventions: `ROOT`-relative paths, explicit SHA-256 checks, RDKit sanitization, NetworkX graphs, per-state JSON results, and a `main()` entry point.
- [ ] Parse the mapped SMILES and legend. Require unique positive atom-map IDs, one PheB1-N attachment placeholder, explicit glucose stereochemistry, the expected C10O8 heavy inventory excluding the attachment placeholder, and the C10H16O8 net substituent formula from an explicit derivation.
- [ ] Calculate stereochemical assignments with RDKit and check them against the source map legend. Reject missing or contradictory stereochemistry with a specific error. Include a separate formula calculation for the B1–B22 peptide fragment and compare its [M+2H]2+ m/z to the published calculated value, keeping digestion and protonation assumptions in the JSON.
- [ ] Parse ATOM/HETATM, CONECT, LINK, and hydrogen records using PDB fixed columns. Reject duplicate serials, missing elements, unresolved LINK endpoints, duplicate selected atom identities, and unexpected selected-node counts. Select PHX B 0 plus PHE B 1 N; verify the single N–C2 LINK.
- [ ] Build a source graph keyed by atom-map ID and a target graph keyed by PDB atom identity. Map the source attachment placeholder to element N using the legend. Match nodes by element and require full graph isomorphism, so every source and target node and edge is covered. Preserve all mappings and group them by source-graph automorphism; never choose one using coordinate distances.

```python
matcher = nx.algorithms.isomorphism.GraphMatcher(
    source_graph,
    target_graph,
    node_match=lambda source, target: source["element"] == target["element"],
)
mappings = [dict(mapping) for mapping in matcher.isomorphisms_iter()]
if not mappings:
    raise ValueError("No complete B1 product graph mapping for this candidate")
```

- [ ] After each complete mapping, carry source bond classes as annotations and calculate raw mapped heavy-atom distances. Label distance ranges as a geometry screen only; do not use them to set bond order, stereochemistry, or mapping choice.
- [ ] Record for each state: input PDB hash, source graph/legend hashes, source and target graph summaries (node count, element inventory, edge count, connected components, and cycle rank), exact map count and all map records, symmetry classes, bond-distance observations, PHX/PheB1-N hydrogen-record observations, mismatch details, software versions, timestamp, and `nvidia_gpu_used: false`.
- [ ] If either candidate does not match exactly, preserve the mismatch in JSON and mark the topology correspondence unresolved. Do not alter candidate PDBs.

### Task 4: Run and inspect the chemistry audit artifacts

**Files:**
- Verify: `research/data/nnc2215-b1-topology-audit.json`
- Verify: `research/data/nnc2215-b1-glucoside-atom-map.json`
- Verify: `research/structures/nnc2215-b1-glucoside-attachment.smi`
- Verify: both corrected candidate PDBs

- [ ] Run the audit from repository root:

```bash
micromamba run -p /Users/andrii/scratch/dia-target-tools/env python research/audit_nnc2215_b1_topology.py
```

Expected: JSON output is regenerated; both state records include graph counts, mapping results, hydrogen observations, input hashes, and explicit limits. A failed exact match is a reportable result, not grounds to modify coordinates.

- [ ] Parse both JSON artifacts with `python3 -m json.tool` and confirm the audit's graph/legend hashes match the checked-in files.
- [ ] Re-run the audit and confirm both corrected PDB files remain byte-identical to the hashes in `research/data/nnc2215-carbon-reconciliation-candidate.json`.
- [ ] Run `git diff --check` and resolve all whitespace errors.

### Task 5: Publish the result in the research report and dashboard

**Files:**
- Modify: `research/NNC2215_CHEMISTRY_RECONCILIATION.md`
- Modify: `research/README.md`
- Modify: `research/index.html`

- [ ] Add the graph, map legend, audit JSON, and rerun command to the chemistry reconciliation report and README.
- [ ] Summarize the exact per-state mapping outcome, stereochemistry/formula checks, any symmetry ambiguity, bond-distance screen, and hydrogen ownership difference from the generated JSON only.
- [ ] Add a concise B1 audit card or section to the existing research dashboard and link all new artifacts. State that the candidate coordinates are model-built and that PDB records do not validate bond order or charge.
- [ ] Preserve current limits: B29 scaffold mappings remain ambiguous; full-conjugate stereochemistry/protonation/charges/parameters and target-binding simulation are not validated by the B1 audit.
- [ ] Confirm every local link target exists and run `git diff --check` again.

### Task 6: Commit, merge, and verify publication

**Files:**
- Commit all B1 graph, audit, and report changes on `codex/nnc2215-b1-topology`.

- [ ] Review the complete diff and confirm that source inputs and both corrected PDBs are unchanged.
- [ ] Commit the implementation and documentation with a message such as `feat: audit NNC2215 B1 topology`.
- [ ] Push the branch and open a pull request against `main`; include the source graph result, exact mapping status, and unresolved limits in the PR description.
- [ ] Merge after repository checks succeed, as previously authorized.
- [ ] Verify the main-branch deployment workflows and fetch the live research page and B1 graph/audit links. Confirm the published files match the merged artifacts and that the dashboard describes the result consistently.

## Completion boundary

This plan closes only the source-mapped B1 graph and its relationship to the two candidate structures. After publication, continue toward the full working model by integrating B1, B29, and the desB30 insulin graph; audit the whole conjugate; then resolve pH-specific protonation and charge, custom force-field parameters, and simulation stability. Do not claim target-binding MD readiness until those later gates pass.
