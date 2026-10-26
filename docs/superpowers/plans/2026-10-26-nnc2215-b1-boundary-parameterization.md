# NNC2215 B1 Local Boundary Parameterization Pilot Plan

> **Scope:** Implement the approved design in `docs/superpowers/specs/2026-10-26-nnc2215-b1-boundary-parameterization-design.md`. This is a reproducible feasibility audit, not force-field validation. Do not add or run a software test suite; the research audit and its recorded checks are the validation artifact.

## Files

- Add `research/audit_nnc2215_b1_ambertools_parameterization_pilot.py`: assemble, coordinate-map, parameterize, and summarize the local fragment.
- Add `research/data/nnc2215-b1-gaff2-boundary-pilot.json`: checksum-bound report for the fragment, mapping, charge bookkeeping, and generated terms.
- Add `research/data/nnc2215-b1-gaff2-boundary-pilot/`: SDF, MOL2, frcmod, and command logs.
- Update `research/README.md`, `research/COMPUTE_STATUS.md`, `research/NNC2215_CHEMISTRY_RECONCILIATION.md`, `research/MILESTONE_REPORT.md`, and `research/index.html` with the scoped result and exact rerun command.

## Task 1: Assemble and validate the source fragment

- [ ] Load the checked full-conjugate graph through its existing builder and retain the B1 source maps, PheB1 and ValB2 template atoms, and AsnB3 N cap atom. Keep only maps 2–19 for the B1 substituent; map 1 is the PheB1 N protein endpoint, not an extra atom.
- [ ] Require one connected fragment with 37 heavy atoms, formula C24H37N3O10, neutral formal charge, the expected element inventory, and complete valid stereochemistry: CIP S at PheB1 CA and ValB2 CA; CIP R, R, S, S, R at B1 source maps 9, 11, 13, 15, 17.
- [ ] Cross-check the open candidate hash against both the B1 topology audit and the coordinate-reconciliation record. Require exactly one B1 embedding and use its map-to-serial correspondence; map remaining PheB1/ValB2/AsnB3 heavy atoms by unique PDB identity.
- [ ] Reject duplicate or missing PDB identities, wrong elements, missing coordinates, graph/map disagreement, unexpected hydrogen counts, or source hashes that do not match the audits. Snapshot the candidate PDB hash before execution.
- [ ] Copy all 37 heavy-atom coordinates without minimization. Recompute implicit hydrogens after trimming downstream AsnB3 atoms, then add explicit hydrogens with RDKit. Record that the resulting AsnB3 N is a neutral NH2 cap and PheB1 N carries one graph-derived H; do not infer PDB hydrogen ownership.

## Task 2: Run the pinned AmberTools feasibility pilot

- [ ] Require the recorded lockfile's AmberTools package and the active `AMBERHOME`, `antechamber`, and `parmchk2` executables.
- [ ] Write the exact prepared SDF, then run Antechamber with AM1-BCC, GAFF2, formal charge 0, and `maxcyc=0`; run parmchk2 with ff14SB available as a comparison lookup set. Use isolated temporary working directories and retain normalized outputs and logs.
- [ ] Record command templates, environment/package versions, input and output SHA-256 values, atom/bond counts, atom-type counts, partial charge per source map/PDB serial, total charge residual, and generated parameter records/penalties. Never apply a charge correction.
- [ ] If Antechamber or parmchk2 fails, preserve logs and mark the pilot partial/failed; do not report an assigned parameter candidate. Check the coordinate PDB hash is unchanged after execution.

## Task 3: Regenerate and inspect the scientific artifact

- [ ] Run the audit from the pinned environment and confirm the neutral fragment and source-map checks pass.
- [ ] Repeat the full audit; compare the aggregate JSON and all retained artifacts byte-for-byte for reproducibility.
- [ ] Confirm the charge sum is within the existing 0.01 e bookkeeping tolerance without correction, and inventory every high-penalty or missing parameter record without treating generated terms as validated.
- [ ] Run `git diff --check`; inspect all reported formula, mapping, cap, coordinate, and scope fields. No separate tests are added or run.

## Task 4: Publish the bounded milestone

- [ ] Update the README, compute status, chemistry reconciliation, milestone report, and dashboard from the generated JSON only. Link the input SDF, output artifacts, JSON, and rerun script.
- [ ] State that the pilot uses active/open model-built coordinates, a single audited B1 embedding, and a neutral AsnB3 NH2 cap; preserve PheB1 hydrogen-record, B29, whole-protein microstate, cross-boundary, and complete force-field limitations.
- [ ] Inspect all new links and check that no copy claims a simulation-ready system, validated charges/terms, MD, affinity, or clinical result.

## Task 5: Review, merge, and deliver

- [ ] Inspect the complete diff and ensure no source graph, candidate coordinate, or existing audit was changed except the intended new output artifacts.
- [ ] Obtain an independent read-only premerge review; correct any material findings and reproduce the affected audit.
- [ ] Commit and push the branch, open a PR to `main`, and merge under the standing user authorization after repository checks succeed.
- [ ] Confirm both production workflows succeed. Fetch the live research page and JSON; verify the published JSON byte-matches the merged artifact and the card links resolve.
- [ ] Reply to Patrick in the existing project thread with the merged PR, live result, key measurements, and explicit feasibility-only scope; verify the sent reply.

## Completion boundary

This closes only the B1 local-fragment toolchain feasibility gate. It does not select a B29 map, fix PheB1 hydrogen records, validate parameter transfer or cross-boundary terms, choose a physiological whole-insulin microstate, build a complete OpenMM system, or establish target-binding MD readiness.
