# NNC2215 Roof Graph Coordinate Mapping Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Reproducibly audit candidate subgraph embeddings of the isolated roof macrocycle in the corrected open and closed PDB candidates.

**Architecture:** A focused Python audit script builds an RDKit graph from the mapped SMILES, builds each LYX/MC heavy-atom graph from PDB CONECT/LINK adjacency, and uses NetworkX subgraph monomorphism with documented source-tail anchors. It writes all matches and limits to JSON; the dashboard links to the audit.

**Tech Stack:** Python, RDKit 2026.03.1, NetworkX, PDB fixed-column records, JSON.

---

### Task 1: Build the source-to-coordinate candidate audit

**Files:**
- Create: `research/audit_nnc2215_roof_graph_coordinate_mapping.py`
- Create: `research/data/nnc2215-roof-graph-coordinate-mapping-audit.json`

- [x] Read the mapped SMILES and corrected PDB paths/hashes from the checked-in graph and coordinate reconciliation records.
- [x] Parse PDB atoms, CONECT records, and LINK records; select heavy atoms from LYX B30 and MC C1 and verify their inventory is C62N15O14.
- [x] Build an element-labeled graph for the target, use NetworkX `subgraph_monomorphisms_iter` with the documented candidate map IDs 78–83 fixed to their linker-audit PDB atoms, and retain every valid embedding.
- [x] Record the map ID to PDB atom candidates for every embedding, unmatched fragment atoms, graph/match counts, hashes, software versions, anchor provenance, and explicit adjacency-only limitations.
- [x] Run: `micromamba run -p /Users/andrii/scratch/dia-target-tools/env python research/audit_nnc2215_roof_graph_coordinate_mapping.py`
- [x] Expected: both states have 91 selected heavy atoms, C62N15O14, 16 anchored embeddings, 83 mapped graph atoms, and eight unmatched B29-linker atoms.

### Task 2: Publish the mapping result with limits

**Files:**
- Modify: `research/NNC2215_CHEMISTRY_RECONCILIATION.md`
- Modify: `research/README.md`
- Modify: `research/index.html`
- Modify: `research/MILESTONE_REPORT.md`

- [x] State that both corrected candidates contain the source graph as an adjacency-only subgraph under the previously documented candidate tail anchors.
- [x] State that 16 embeddings remain and that no single authoritative graph-to-PDB atom assignment or bond-order validation follows from this result.
- [x] Link the audit JSON and script from the dashboard and document the eight unmatched B29-linker atoms.
- [x] Confirm that the input PDB hashes still match the coordinate reconciliation record and no PDB file changed during this audit.

### Task 3: Verify artifacts

**Files:**
- Verify: `research/data/nnc2215-roof-graph-coordinate-mapping-audit.json`
- Verify: `research/index.html`
- Verify: all JSON records touched by this plan

- [x] Run: `python3 -m json.tool research/data/nnc2215-roof-graph-coordinate-mapping-audit.json`
- [x] Run: `git diff --check`
- [x] Confirm the graph, audit, map legend, and script link targets exist under the published `research/` tree.
