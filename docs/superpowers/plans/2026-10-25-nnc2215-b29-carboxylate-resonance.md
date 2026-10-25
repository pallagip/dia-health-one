# NNC2215 B29 Carboxylate Resonance Audit Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a conditional pH 7.4 resonance-normalized B29 mapping analysis without changing the neutral graph audit or coordinate files.

**Architecture:** Extend `audit_nnc2215_b29_topology.py` with a graph copy that normalizes only the three macrocycle carboxylate C–O pairs, then group the existing exact embeddings under that graph's automorphisms. Store the conditional results in the generated JSON and describe them in the research page and chemistry reconciliation note.

**Tech Stack:** Python 3, NetworkX, RDKit, JSON, static HTML.

---

## Files

- Modify `research/audit_nnc2215_b29_topology.py`: declare the three carboxylate map groups; build a resonance-normalized source graph; write conditional analysis fields without changing neutral-state calculations.
- Regenerate `research/data/nnc2215-b29-topology-audit.json`: preserve existing neutral mapping data and add the pH 7.4 scenario, formula/charge delta, normalized equivalence groups, and interpretation.
- Modify `research/NNC2215_CHEMISTRY_RECONCILIATION.md`: distinguish the neutral eight classes from the one normalized class under the declared tri-carboxylate assumption.
- Modify `research/index.html`: explain the conditional equivalence result while leaving the complete chemistry gate open.

## Task 1: Add resonance-normalized mapping classes

**Files:** `research/audit_nnc2215_b29_topology.py`

- [ ] Define `CARBOXYLATE_MAP_GROUPS = ((52, 53, 54), (60, 61, 62), (68, 69, 70))`.
- [ ] Add a helper that copies the source graph, verifies each group is C bonded to exactly the two declared O atoms, and changes only those six edges to a shared `carboxylate_resonance` bond class with numeric order `1.5`.
- [ ] For each state, call `mapping_symmetry_classes` with automorphisms of the normalized graph and the already enumerated exact embeddings. Store this as a separate conditional result; do not replace `exact_mapping_symmetry_classes`.
- [ ] Before emitting JSON, require three groups, neutral formula `C62H72N16O14`, three assumed proton removals, local charge delta `-3`, 16 mappings per state, eight neutral classes per state, and one normalized class per state. Raise a descriptive error if any value differs.
- [ ] Record that the normalized graph is an equivalence representation only: it does not assign PDB oxygen names, partial charges, or a force field.
- [ ] Regenerate the JSON with `micromamba run -p /Users/andrii/scratch/dia-target-tools/env python research/audit_nnc2215_b29_topology.py` and confirm both PDB hashes remain unchanged.

## Task 2: Update chemistry reporting

**Files:** `research/NNC2215_CHEMISTRY_RECONCILIATION.md`, `research/index.html`

- [ ] Update the top decision paragraph and B29 section to state that the neutral source graph retains eight label-transfer classes, all differing only in the three carboxyl oxygen assignments.
- [ ] Add the declared pH 7.4 tri-carboxylate result: the normalized graph yields one chemistry-equivalent class per state, with no coordinate-based atom assignment.
- [ ] Keep the full-conjugate topology, site-specific protonation, charge ledger, and force-field readiness open.
- [ ] Update the dashboard card and local B29 summary to present the neutral result and conditional charge-state result together without marking the overall gate complete.

## Task 3: Review and commit

**Files:** all files listed above.

- [ ] Run `micromamba run -p /Users/andrii/scratch/dia-target-tools/env python research/audit_nnc2215_b29_topology.py` and confirm `assessment.both_states_have_exact_full_element_labeled_adjacency_isomorphism` remains true, neutral classes remain eight per state, normalized classes equal one per state, and both coordinate hashes are unchanged.
- [ ] Run `git diff --check` and inspect the rendered dashboard text in `research/index.html`.
- [ ] Commit with `git add docs/superpowers/specs/2026-10-25-nnc2215-b29-carboxylate-resonance-design.md docs/superpowers/plans/2026-10-25-nnc2215-b29-carboxylate-resonance.md research/audit_nnc2215_b29_topology.py research/data/nnc2215-b29-topology-audit.json research/NNC2215_CHEMISTRY_RECONCILIATION.md research/index.html && git commit -m "feat: audit B29 carboxylate resonance mappings"`.
