# NNC2215 Conditional pH 7.4 Charge-State Ledger Plan

> **For agentic workers:** Use the approved scope in [the design](../specs/2026-10-26-nnc2215-physiological-charge-ledger-design.md); preserve all uncertainty as explicit data.

**Goal:** Record a pH 7.4 formal-charge candidate and enumerate the unresolved ionization states needed before atom typing or force-field preparation.

**Architecture:** Read the merged source-derived neutral graph JSON, validate titratable group identities from provenance, enumerate group-level protonation/tautomer combinations, and emit deterministic JSON. The audit does not modify coordinates or rewrite the neutral graph.

---

## Task 1: Build the charge-state audit

**Files:** `research/audit_nnc2215_physiological_charge_states.py`, `research/data/nnc2215-physiological-charge-states.json`

- [x] Locate four Glu sidechain carboxyls, both chain C termini, all three mapped B29 macrocycle carboxyls, A-chain N terminus, ArgB22, and HisB5/HisB10 from graph provenance; fail closed on identity, neutral-group topology, sequence, and conjugation-boundary mismatches.
- [x] Record fixed conventional pH 7.4 assumptions, uncertain sites, excluded high-pKa sensitivity, and unresolved PheB1 hydrogen ownership.
- [x] Enumerate three binary macrocycle acid states, two A1 amino states, and three states for each histidine, with no assigned probabilities.
- [x] Calculate state formula hydrogen counts and group-summed net formal charges relative to the neutral graph; record the all-carboxylate, protonated-A1, neutral-histidine candidate as −7 and C325H454N79O97S6.
- [x] Write deterministic JSON with input SHA-256, site identities, exact assumptions, state enumeration, corrected full-ring histidine hydrogen counts, checks, and force-field/MD limitations.

## Task 2: Publish the evidence and limitations

**Files:** `research/README.md`, `research/NNC2215_CHEMISTRY_RECONCILIATION.md`, `research/index.html`

- [x] Link the audit JSON and exact rerun command.
- [x] Describe the −7 charge as a conventional candidate under stated assumptions, not a measured microstate or computed population.
- [x] State that formal group charges are not partial atomic charges or force-field parameters.
- [x] Preserve the separate neutral constitutional graph and all open mapping, hydrogen, stereochemistry, and force-field gates.

## Task 3: Verify

- [x] Run the exact scoped audit command in the configured chemistry environment if needed; inspect output, group identities, formula/charge arithmetic, and JSON determinism.
- [x] Run `git diff --check` and inspect the site prose and links.
- [x] Obtain independent scientific/source and code review; resolve blockers before integration.

## Task 4: Integrate and communicate

- [ ] Commit and push the feature branch; open and merge a reviewed PR.
- [ ] Verify the main production deployments and compare the live dashboard and JSON with the merged files.
- [ ] Send Patrick a milestone reply in the existing email thread with concrete links and remaining model gates.
