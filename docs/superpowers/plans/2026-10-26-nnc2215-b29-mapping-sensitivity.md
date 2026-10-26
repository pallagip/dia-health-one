# NNC2215 B29 Mapping Sensitivity Audit Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Determine whether the isolated B29 GAFF2/AM1-BCC candidate assignment is stable across every exact active/open coordinate embedding recorded in the topology audit.

**Architecture:** Extend the reproducible pilot to select one of the 16 verified embeddings by index and record its full map-to-PDB serial mapping. Add an orchestrator that runs the pinned pilot for all 16 embeddings with four concurrent workers and all eight carboxylate Lewis labelings per embedding, then compares candidate charges and atom types by physical PDB serial. Preserve each run's inputs, outputs, logs, and hashes; the aggregate remains a mapping-sensitivity audit, not a validated force field.

**Tech Stack:** Python 3.12, RDKit, AmberTools 26 Antechamber/parmchk2, GAFF2, AM1-BCC, `concurrent.futures`.

---

### Task 1: Make one exact mapping an explicit pilot input

**Files:**
- Modify: `research/audit_nnc2215_ambertools_parameterization_pilot.py`
- Refresh: `research/data/nnc2215-b29-gaff2-pilot.json`

- [x] Add `--mapping-index` with a default of `0` and reject values outside `0..15`.
- [x] Select the requested entry from the topology audit's 16 exact mappings after checking the active/open PDB hash against the topology audit's expected SHA-256.
- [x] Require the selected map to contain exactly source map IDs `1..92` and 92 distinct PDB serials.
- [x] Record the zero-based index, mapping ID, and full map-ID-to-PDB-serial mapping in the report.
- [x] Run the default pilot once and confirm its status and prior eight-variant result are preserved with the updated script hash.

### Task 2: Run the complete 16-by-8 sensitivity matrix

**Files:**
- Create: `research/audit_nnc2215_b29_mapping_sensitivity.py`
- Create: `research/data/nnc2215-b29-gaff2-mapping-sensitivity.json`
- Create: `research/data/nnc2215-b29-gaff2-mapping-sensitivity/`

- [x] Launch one pilot process for each mapping index `0..15` with `ThreadPoolExecutor(max_workers=4)` and a separate output directory for each map. Reuse the default pilot directory for index `0`; write indices `1..15` under `mapping_02` through `mapping_16`.
- [x] Require all 16 reports to show the same graph, map legend, topology, PDB, lockfile, and AmberTools package inputs; require eight successful variants per mapping and 128 total variant records.
- [x] Aggregate MOL2/frcmod hash diversity and compare heavy-atom partial charges and GAFF2 types by PDB serial across all mappings and resonance labelings.
- [x] Record per-serial charge range and type set, each mapping's output hashes, and whether the mapped PDB was unchanged.
- [x] Keep the aggregate status explicitly at `MAPPING_SENSITIVITY_ONLY; COMPLETE_CONJUGATE_FORCE_FIELD_NOT_VALIDATED` and report zero validated conjugate force fields.

### Task 3: Reproduce and inspect the sensitivity result

**Files:** all audit and result files above.

- [x] Run the complete 16-by-8 audit twice in `/private/tmp/dia-target-tools-pinned`; require byte-identical aggregate JSON and all per-map artifact hashes.
- [x] Run `git diff --check` and inspect the aggregate for missing mappings, failed processes, charge residuals, and high-penalty parameter records.
- [x] Confirm that no source PDB or upstream topology file was modified.

### Task 4: Publish the next evidence gate

**Files:**
- Modify: `research/COMPUTE_STATUS.md`
- Modify: `research/README.md`
- Modify: `research/NNC2215_CHEMISTRY_RECONCILIATION.md`
- Modify: `research/MILESTONE_REPORT.md`
- Modify: `research/index.html`

- [x] Add the matrix result using only values recorded in its JSON, including whether outputs are stable or mapping-sensitive.
- [x] Link the aggregate JSON and reproducible script; preserve the unresolved LysB29 cap, PheB1 junction, and complete-force-field gates.
- [x] State explicitly that no full insulin topology, OpenMM system, MD, affinity, or clinical conclusion results from this audit.

### Task 5: Merge and verify publication

- [ ] Commit and push the change, open a pull request to `main`, obtain an independent review, and merge under the user's standing authorization.
- [ ] Confirm the Cloudflare Pages and GitHub Pages production workflows succeed.
- [ ] Fetch the live dashboard and aggregate JSON; require the dashboard card to be present and the live JSON bytes to match the merged file.
- [ ] Send Patrick the milestone summary and live links in the existing email thread.

### Scope limits

This audit tests mapping and resonance sensitivity for an isolated, capped B29 fragment. It does not assign final LysB29 boundary charges, resolve stereochemistry from experimental evidence, parameterize the B1 junction or whole insulin, validate GAFF2 transferability, build a complete topology or OpenMM system, run MD, estimate binding, or support clinical use.
