# NNC2215 Roof Macrocycle Graph Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Publish a source-traceable, RDKit-validated graph of the NNC2215 roof macrocycle propyl azide as the next auditable component of the working model.

**Architecture:** Keep the mapped molecular structure in a single `.smi` source file and put graph checks/provenance generation in a focused Python script. The script writes a machine-readable audit JSON; the dashboard links to both artifacts and reports their bounded scope.

**Tech Stack:** Python 3, RDKit 2026.03.1 from the pinned local chemistry environment, JSON, static HTML.

---

### Task 1: Add the source-derived mapped graph

**Files:**
- Create: `research/structures/nnc2215-roof-propyl-azide.smi`
- Create: `research/data/nnc2215-roof-macrocycle-atom-map.json`

- [x] Transcribe the numbered Extended Data Fig. 8i roof-azide structure into mapped isomeric SMILES, using the source-drawn neutral-acid form and preserving source numbering for C21–C23.
- [x] Record stable map-ID ranges and their structural groups, identify the three numbered tail carbons, and state that the remaining IDs are graph IDs rather than authoritative PDB atom names.
- [x] Record the chosen azide resonance form and state that atom-to-PDB mapping, protonation at simulation pH, and force-field assignments remain unresolved.

### Task 2: Validate the graph and write its audit

**Files:**
- Create: `research/validate_nnc2215_roof_macrocycle_graph.py`
- Create: `research/data/nnc2215-roof-macrocycle-graph-audit.json`

- [x] Parse the `.smi` with RDKit, sanitize it, calculate formula/heavy atoms/formal charge/ring count, and inspect heavy-atom map IDs.
- [x] Fail with a clear message unless the formula is C56H63N15O12, the heavy-atom count is 83, every heavy atom has a unique positive map ID, and the total formal charge is zero.
- [x] Write an audit JSON containing the source URLs, input SHA-256, RDKit version, canonical mapped SMILES, formula, counts, assumptions, and readiness limits.
- [x] Run: `micromamba run -p /Users/andrii/scratch/dia-target-tools/env python research/validate_nnc2215_roof_macrocycle_graph.py`
- [x] Expected: exit code 0 and an audit reporting formula `C56H63N15O12`, 83 heavy atoms, 83 unique maps, and formal charge 0.

### Task 3: Update research documentation and dashboard

**Files:**
- Modify: `research/NNC2215_CHEMISTRY_RECONCILIATION.md`
- Modify: `research/README.md`
- Modify: `research/index.html`

- [x] Add the graph as a completed component-level milestone while retaining full-conjugate, PDB mapping, protonation, parameterization, and target-MD blockers.
- [x] Add dashboard links to the `.smi` and audit JSON and label the graph as a source-derived candidate.
- [x] State that the graph is neutral-acid form and that the azide resonance notation is representational; do not imply force-field or target readiness.
- [x] Inspect the rendered section or source markup and verify both relative links point to files in the published `research/` tree.

### Task 4: Verify the public status correction and record consistency

**Files:**
- Verify: `research/COMPUTE_STATUS.md`
- Verify: `research/NVIDIA_INNOVATION_LAB_BRIEF.md`
- Verify: `research/data/nvidia-access-status.json`
- Verify: `research/index.html`

- [x] Confirm NVIDIA Inception is stated as active based on the 2026-10-19 support email, while the last portal observation remains dated 2026-10-18 and no GPU/credit award is claimed.
- [x] Confirm the obsolete Innovation Lab application brief is marked archived and the public card includes the official Inception program link.
- [x] Run: `python3 -m json.tool research/data/nvidia-access-status.json`
- [x] Run: `git diff --check`
