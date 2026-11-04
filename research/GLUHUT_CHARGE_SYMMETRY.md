# GluHUT compiled-charge symmetry audit

**Date:** 4 November 2026  
**Status:** compiled charge model is exactly equivalent across every symmetry preserved by the parameter-labelled host graph; source RESP derivation remains unavailable  
**Scope:** public methyl-truncated GluHUT host only

## Question

The deposited host graph is chemically repetitive: two caps, three pillars and six urea links. Before treating its compiled charges as a stable public benchmark, graph-equivalent atoms were tested for charge consistency.

This matters because a RESP calculation performed on one asymmetric conformation can assign different charges to chemically equivalent atoms unless equivalencing constraints are used. Large unexplained asymmetry would make direct parameter transfer to a new analogue particularly unsafe.

## Method

The audit used the checksum-verified compiled AMBER topology and constructed the 84-heavy-atom host graph.

Graph automorphisms were enumerated using these node labels:

- element;
- AMBER atom type;
- full degree;
- heavy-atom degree;
- attached-hydrogen count.

Edges were labelled with:

- compiled bond force constant;
- compiled equilibrium bond length.

These labels deliberately test symmetries preserved by the deposited parameterised model, not every possible unlabeled chemical-graph symmetry.

Attached hydrogens were grouped by:

- the automorphism orbit of their parent heavy atom;
- AMBER hydrogen atom type.

## Result

The compiled parameter-labelled heavy graph has:

- 84 heavy atoms;
- 90 heavy-heavy bonds;
- **6 automorphisms**;
- **28 heavy-atom symmetry orbits**;
- every heavy-atom orbit has size **3**.

Across all 28 heavy-atom orbits:

- maximum charge range: **0.0 e**;
- every charge is exactly equal within `1 × 10⁻⁸ e`.

Across every corresponding hydrogen group:

- maximum charge range: **0.0 e**;
- every charge is exactly equal within the same tolerance.

A hypothetical symmetry averaging operation produces:

- total-charge change: **0.0 e**;
- RMS per-atom adjustment: **0.0 e**;
- maximum per-atom adjustment: **0.0 e**.

The complete host charge vector is also byte-equivalent across all four deposited sugar systems.

## Interpretation

The deposited model contains exact threefold charge equivalencing for every graph symmetry retained by its own atom types and bonded parameters. This is strong internal-consistency evidence for the compiled benchmark.

It also means the public host is not carrying large arbitrary RESP differences among the three repeated pillars or the three corresponding positions identified by the parameter-labelled graph.

This result does **not** prove that:

- every chemically possible cap-exchange or larger symmetry is preserved by the compiled atom types and bonded parameters;
- the original RESP calculation used explicit equivalencing constraints;
- the source electrostatic potential or conformer was appropriate;
- the charge model is transferable to PaG, MC-Q1 or complete NNC2215.

The source RESP conformers, ESP grids, equivalencing controls and final source `mol2` remain absent. Exact compiled symmetry cannot substitute for that provenance.

## Decision

1. Accept the compiled charge vector as internally symmetry-consistent for exact reproduction of the public GluHUT sugar benchmark.
2. Do not apply any automatic charge symmetrisation; it would make no numerical change under the tested orbits.
3. Do not transfer the host charges directly to MC-Q1 or other analogues.
4. Request the source RESP conformers, ESP calculations, equivalencing constraints and final `mol2` file.
5. Keep source-matched parent glucose-versus-PaG calibration blocked until the charge derivation and starting-coordinate provenance are sufficiently resolved.
6. Keep MC-Q1 ranking blocked.

## Reproducibility

- Script: `research/audit_gluhut_charge_symmetry.py`.
- Machine-readable result: `research/data/gluhut-charge-symmetry.json`.
- Tests: `tests/test_audit_gluhut_charge_symmetry.py`.
- Workflow: `.github/workflows/gluhut-charge-symmetry.yml`.
- NetworkX version: `3.4.2`.
- Host charge-vector SHA-256: `eedb21f8da70b1803d7ee02d7d37819f11c414cd0aa0b44ff64744f6858d0a3f`.

## Evidence boundary

This is an internal-consistency result for compiled public benchmark charges. It is not source RESP reconstruction, chemical validation, an NNC2215 result, a PaG result, an MC-Q1 result, a binding-affinity result, efficacy, safety or a working-drug result.
