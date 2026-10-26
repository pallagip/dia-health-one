# NNC2215 Attachment-Boundary Audit Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Reconcile the source-derived B1 and B29 component graphs with their local insulin attachment contexts in both corrected candidate coordinate files, while surfacing the PheB1 hydrogen-record discrepancy without changing molecular records.

**Architecture:** Add a read-only PDB and graph audit that cross-checks source maps and existing topology audits, resolves three named heavy-atom joins in each state, and measures their coordinate distances only after identities and edges are fixed. Emit a checksum-bound JSON ledger. Update the research README, chemistry reconciliation, and dashboard with bounded conclusions and rerun/download links.

**Tech Stack:** Python 3, RDKit, standard-library PDB/CONECT/LINK parsing, JSON, static HTML.

## Files

- Add `research/audit_nnc2215_attachment_boundaries.py`: parse and validate the source B1/B29 graphs, mapping ledgers, candidate PDB atoms, named `LINK` records, serial `CONECT` records, local neighbors, and distances.
- Add `research/data/nnc2215-attachment-boundary-audit.json`: store inputs and hashes, per-state boundary rows, observed neighborhoods/H records, discrepancies, and scoped assessment.
- Modify `research/README.md`: document the audit, rerun command, candidate-coordinate scope, and unresolved chemistry.
- Modify `research/NNC2215_CHEMISTRY_RECONCILIATION.md`: add a bounded attachment-boundary closure section and retain the unresolved PheB1-H and full-model gates.
- Modify `research/index.html`: add a dashboard result card and downloadable JSON link without presenting full-model or simulation readiness as complete.

## Task 1: Implement the checksum-bound boundary audit

**Files:** `research/audit_nnc2215_attachment_boundaries.py`, `research/data/nnc2215-attachment-boundary-audit.json`

- [x] Define explicit repository-relative inputs for the corrected open/closed candidates, coordinate reconciliation record, B1/B29 graph/map files, and existing B1/B29 topology audit JSONs. Record SHA-256 for every input actually used.
- [x] Parse candidate PDB atom identity using model/state, record type, chain, residue number/insertion code, residue name, atom name, element, and serial. Reject missing or duplicate identities for required endpoints. Parse `LINK` identities and `CONECT` serial adjacency separately; do not infer bonds from distances.
- [x] Parse the source B1 and B29 SMILES/map legends and recover the three source boundary edges: B1 maps 1–2, B29 maps 84–92, and roof-tail maps 78–79. Check these against the source graph edge sets and legends.
- [x] Cross-check the existing topology audit correspondence and mapping counts before reporting: one exact B1 mapping per state and 16 exact B29 mappings per state. Preserve ambiguity; do not select among B29 embeddings by coordinates.
- [x] Resolve the six expected candidate heavy-atom pairs by full identity in open and closed states: PheB1 N–PHX C2, LYS B29 NZ–LYX B30 CD1, and MC C31–LYX B30 CD6 (open) / C (closed). Require each named `LINK` once and each atom identity once.
- [x] Resolve native context pairs PheB1 N–CA and LysB29 CE–NZ from candidate `CONECT`/`LINK` adjacency and assert a unique observation. Report each attachment and native-context distance as a coordinate observation after resolving identities; do not rank mappings or infer bond order from geometry.
- [x] For PheB1 N and LysB29 NZ, record adjacent heavy-atom identities and CONECT hydrogen identities. Separately measure the PheB1 N to nearby PHX H2A distance when present. Require the observed open H2 and closed absent-neighbor/nearby H2A pattern to be explicitly surfaced as a discrepancy; do not assign the hydrogen to another atom or edit any file.
- [x] Follow-up: compare the open PheB1 N–H2 coordinates with the closed PHX H2A coordinates using N–H distance, local N-centered angles, and nearest heavy atoms. Record the comparison as geometry only; it must not add connectivity or assign hydrogen ownership.
- [x] Check expected neutral local boundary patterns against observed records as observations: PheB1 N context CA + B1 carbonyl + one hydrogen; LysB29 NZ context CE + carbamate carbonyl + one hydrogen. Keep incomplete or inconsistent hydrogen records separate from the heavy-atom closure result.
- [x] Hash the candidate PDBs before and after parsing and fail if either changes. Check existing source graph/map/audit hashes or graph summaries and exact-match counts against the recorded audit inputs; fail on unexpected changes.
- [x] Emit `dia-health-one.nnc2215-attachment-boundary-audit.v1` with source paths/hashes, state-specific boundary rows and distances, local neighbor ledgers, discrepancy list, limits, and an assessment that can pass only when all six unique heavy-atom joins and both native context joins resolve. Keep hydrogen chemistry and full conjugate readiness open.

## Task 2: Update reproducibility and public-facing research text

**Files:** `research/README.md`, `research/NNC2215_CHEMISTRY_RECONCILIATION.md`, `research/index.html`

- [x] Add the exact rerun command for the new script, result JSON path, and dashboard download link.
- [x] State that the three heavy-atom attachment boundaries and two native residue context edges close in the corrected model-built candidates if the generated audit supports that result.
- [x] Report the PheB1 open/closed hydrogen-record discrepancy and nearby closed PHX H2A observation without assigning or correcting hydrogen ownership. Describe the LysB29 NZ H record as coordinate-file connectivity evidence only.
- [x] Record the open H2 / closed PHX H2A distance and local-angle comparison, while keeping the absent closed CONECT edge and unresolved ownership explicit in the README, reconciliation notes, compute status, milestone report, and dashboard.
- [x] Keep the 16 B29 mappings per state and existing carboxylate-label ambiguity intact. Do not imply the local boundary audit makes a unique full source-to-PDB atom map.
- [x] Explicitly leave full conjugate graph/bond orders/stereochemistry, whole-protein protonation and charge, force-field validation, target-binding MD, and clinical readiness unresolved.
- [x] Label the coordinate files as corrected model-built candidates, not experimentally solved structures, and avoid changing any PDB, SMILES, source map, or existing topology-audit result.

## Task 3: Run the scoped audit and review its output

**Files:** all files listed above.

- [x] Run `micromamba run -p /Users/andrii/scratch/dia-target-tools/env python research/audit_nnc2215_attachment_boundaries.py` in the configured chemistry environment.
- [x] Confirm the audit reports six unique attachment `LINK` edges, four unique native-context joins, and the expected open/closed PheB1 hydrogen observation; verify all recorded candidate PDB hashes match the reconciliation record and are unchanged after execution.
- [x] Confirm the closed PHX H2A is nearest to PheB1 N and that the report's geometry values remain descriptive only; verify no PDB identity, coordinates, or CONECT edges were changed.
- [x] Inspect the JSON for resolved atom identities, coordinates/distance units, provenance, ambiguity, failure behavior, and scoped limitations; confirm no coordinates or hydrogen ownership were modified or inferred.
- [x] Run `git diff --check` and inspect the changed README, reconciliation text, and dashboard section for accurate scope and working relative links.

## Task 4: Review and commit

**Files:** all files listed above.

- [x] Review the implementation against `docs/superpowers/specs/2026-10-26-nnc2215-attachment-boundary-design.md`; correct any schema, identity, or scope mismatch and rerun only the scoped audit.
- [ ] Commit the plan and audit work with clear, separate commit messages if the plan has already been committed; otherwise include the complete change in the feature commit. Preserve the existing merged main history and do not amend prior commits.

## Self-review checklist

- The source edge assertions are separated from candidate PDB edge observations.
- Distances are descriptive fields only and cannot affect graph matching, endpoint identity, or hydrogen assignment.
- B29’s 16 embeddings and carboxyl oxygen-label ambiguity are preserved.
- A heavy-atom boundary result cannot mark the full model or simulation as ready.
- Candidate inputs are hashed before and after, and the audit writes only its JSON output.
