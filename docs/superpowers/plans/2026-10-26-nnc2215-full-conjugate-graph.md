# NNC2215 Full Conjugate Graph Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Assemble and audit a provenance-tagged neutral constitutional graph for full desB30 insulin–NNC2215 from the existing insulin sequences and mapped B1/B29 source fragments.

**Architecture:** Use pinned RDKit FASTA templates for the two insulin chains, add the three source disulfides, and merge the source-mapped B1 and B29 graphs at their declared protein endpoints. Emit a deterministic JSON ledger and cross-check sequence, disulfides, attachment LINK records, and preserved mapping counts against the existing corrected candidate files and audits.

**Tech Stack:** Python 3.11, RDKit 2026.03.1 in `/Users/andrii/scratch/dia-target-tools/env`, JSON, standard-library PDB parsing, static HTML.

---

## Files

- Create `research/audit_nnc2215_full_conjugate_graph.py`: source-hash validation, RDKit graph assembly, named PDB record parsing, checks, and JSON output.
- Create `research/data/nnc2215-full-conjugate-graph.json`: graph nodes and bonds, provenance, formula, topology, checks, source hashes, limitations, and assessment.
- Modify `research/README.md`: explain the graph, limits, result path, download, and exact rerun command.
- Modify `research/NNC2215_CHEMISTRY_RECONCILIATION.md`: add source-derived full graph status and retain open chemistry gates.
- Modify `research/index.html`: add a result card and JSON download without marking coordinates, force fields, MD, or clinical readiness as complete.
- Create `docs/superpowers/specs/2026-10-26-nnc2215-full-conjugate-graph-design.md` and this plan as tracked design records.

## Task 1: Assemble and audit the full constitutional graph

**Files:** `research/audit_nnc2215_full_conjugate_graph.py`, `research/data/nnc2215-full-conjugate-graph.json`

- [ ] Read chain sequences from `nnc2215-model-formula-audit.json`; assert exact A/B sequences and create separate `Chem.MolFromFASTA` graphs. Preserve RDKit residue number, residue name, and atom name on each template atom and attach explicit chain provenance.
- [ ] Read mapped source SMILES and map legends for B1 and B29. Require map 1 to be the B1 protein-boundary placeholder adjacent to map 2; require map 92 to be the B29 Lys NZ cap adjacent to map 84. Remove each cap, merge remaining fragment atoms, and connect map 2 to PheB1 N and map 84 to LysB29 NZ. Preserve source map IDs on all source-derived nodes.
- [ ] Add only the three declared disulfides: A6 SG–A11 SG, A7 SG–B7 SG, and A20 SG–B19 SG. Sanitize the RDKit graph and record any valence exception as a fatal audit error.
- [ ] Compute formula, heavy-atom element counts, heavy node and edge counts, connected component count, cycle rank, per-atom implicit hydrogen counts, formal charges, and RDKit valence status from the assembled graph. Require the neutral formula to equal `C325H461N79O97S6`.
- [ ] Parse both corrected candidate PDBs without editing them. Verify their A/B residue sequences, the three SSBOND contexts, PheB1 N–PHX C2 and LysB29 NZ–LYX CD1 attachment identities, and the roof-tail LINK identity already recorded by the boundary audit.
- [ ] Cross-check source and corrected-candidate hashes against existing formula/topology/reconciliation/boundary audit records. Require one B1 and sixteen B29 exact mappings per state; serialize counts and ambiguity retained without importing coordinates into graph assembly.
- [ ] Emit schema `dia-health-one.nnc2215-full-conjugate-graph.v1` with per-atom provenance, per-bond endpoints/order/provenance, hashes, formula/topology checks, PDB cross-checks, assessment, and explicit unresolved gates. Label the graph as neutral constitutional and distinguish the conditional tri-carboxylate audit.

## Task 2: Publish reproducibility and scoped conclusions

**Files:** `research/README.md`, `research/NNC2215_CHEMISTRY_RECONCILIATION.md`, `research/index.html`

- [ ] Add the audit command `micromamba run -p /Users/andrii/scratch/dia-target-tools/env python research/audit_nnc2215_full_conjugate_graph.py`, the result JSON path, and a working dashboard download link.
- [ ] State the graph formula and connectivity findings only if the emitted JSON supports them. Describe formula agreement as internal source-graph consistency with the published formula, not independent experimental confirmation.
- [ ] Keep the B29 16-mapping count per state, PheB1 H2/H2A discrepancy, conditional B29 tri-carboxylate representation, protonation/charge, histidine tautomer, roof stereochemistry, force field, MD, and clinical readiness as unresolved or explicitly representation-dependent.
- [ ] Label corrected coordinate candidates as model-built and state that they were used for named sequence/connectivity cross-checks only; do not change coordinates, source SMILES, source maps, or prior audit outputs.

## Task 3: Run the scoped audit and inspect generated evidence

**Files:** all files above.

- [ ] Run the exact audit command in the configured chemistry environment and inspect the full summary and exit code.
- [ ] Inspect the JSON to confirm graph node/bond counts reconcile with its formula and cycle rank, every atom has provenance and implicit-hydrogen data, all required links/disulfides resolve in each state, source hashes match their recorded expected values, and the B29 mapping ambiguity remains 16 per state.
- [ ] Confirm the generated report does not assert full-model, simulation, binding, or clinical readiness and does not derive atom identity from coordinates.
- [ ] Run `git diff --check` and inspect changed prose and dashboard links for accurate scope and relative paths.

## Task 4: Review and integrate

**Files:** all files above.

- [ ] Complete an independent design-compliance review followed by a separate code-quality review; resolve any blocking finding and rerun the scoped audit after code changes.
- [ ] Commit the design, plan, implementation, result data, and documentation with clear messages. Push the feature branch, create a PR, and merge it to `main` after the required review and audit evidence are clean.
- [ ] Verify deployment workflows and fetch the live dashboard and result JSON; compare the live assessment and result schema to the merged output.
- [ ] Send Patrick a milestone update in the existing thread with the PR, production page, JSON, formula/connectivity result, and remaining model limitations.

## Scope guard

Do not add coordinates, optimize conformers, assign a full physiological-pH microstate or total charge, resolve the B29 coordinate map, repair PDB hydrogen records, generate force-field inputs, or run molecular dynamics in this milestone.
