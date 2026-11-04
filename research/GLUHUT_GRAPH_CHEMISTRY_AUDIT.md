# GluHUT compiled graph chemistry audit

**Date:** 4 November 2026  
**Status:** compiled-system chemistry audit complete; source RESP/GAFF/GLYCAM derivation remains unavailable  
**Scope:** public truncated GluHUT benchmark systems only

## Question

Two points in the deposited AMBER systems needed resolution before any parent glucose-versus-PaG calibration:

1. What exact graph and molecular composition are encoded for the 159-atom host?
2. Is the two-atom residue labelled `ROH` a separate counterion, or is it covalently part of the carbohydrate guest?

## Verified archive

The audit used the checksum-verified public archive:

- `SMD_trajectories.zip`;
- 37,931,874 bytes;
- MD5 `9f3c6733dc608d2bd37c488c560a5443`;
- SHA-256 `04c58666328764aeb47bac0d9f8607546f9763892cb9a333b108630df70e3943`;
- ZIP CRCs verified.

All four deposited `parm7` files were parsed directly.

## Host graph

The common synthetic host is one connected 159-atom component:

- formula: **`C60H75N15O9`**;
- heavy atoms: 84;
- compiled partial charge: approximately zero (`−2.98 × 10⁻⁶ e` from rounded AMBER values);
- heavy-atom bonds: 90;
- heavy-graph cycle rank: 7;
- methyl groups: 9;
- carbonyl motifs: **6 urea carbonyls and 3 amide carbonyls**.

The formula, motif inventory and graph-level properties are identical across glucose, galactose, fructopyranose and fructofuranose.

This supports the published description of a truncated GluHUT-like host containing two tris-urea receptor halves and three retained amide linker carbonyls. The exact source drawing and synthetic intermediate still need to be mapped atom by atom before the compiled graph is declared fully reconciled.

## `ROH` is not a counterion

The previous residue-level charge table could be misread as a positively charged sugar balanced by a negatively charged `ROH` residue. The bond graph resolves this.

In every system, atom 160 of the carbohydrate residue is covalently bonded to atom 183, oxygen `O1` of `ROH`:

- glucose and galactose: guest `C1–O1(ROH)`;
- fructopyranose and fructofuranose: guest `C2–O1(ROH)`;
- equilibrium C–O bond length: 1.43 Å;
- force constant: 320 kcal mol⁻¹ Å⁻² in the compiled AMBER table.

The two residue records together form one connected, neutral 24-atom guest with formula:

**`C6H12O6`**

Therefore:

- `ROH` is a GLYCAM-style terminal hydroxyl fragment used to complete the monosaccharide;
- it is not a separate counterion;
- the apparent `+0.194 e` sugar-residue charge and `−0.194 e` `ROH` charge are an internal residue-template partition of one neutral molecule.

This correction removes a false chemistry blocker. It does not validate the derivation of the guest charges.

## What this increases confidence in

- The deposited guest in each topology is a complete neutral hexose, not an ion pair.
- All four systems use one common host graph and one consistent guest-residue construction pattern.
- The public sugar benchmark can be interpreted as a host-versus-neutral-hexose comparison at the compiled-topology level.
- The exact compiled host can now be fingerprinted and compared with a published chemical graph.

## What remains unresolved

The deposit still lacks:

- source host `mol2`, RESP inputs and RESP outputs;
- quantum-chemistry coordinates and electrostatic-potential calculation;
- original GAFF and GLYCAM assembly inputs;
- tleap or equivalent construction scripts;
- standalone simulation-start coordinates;
- a source-matched PaG model;
- validated MC-Q1 parameters.

Compiled consistency is not source-level validation. No charge or bonded term is assumed transferable to PaG, MC-Q1 or complete NNC2215.

## Decision

1. Accept the four public systems as internally consistent neutral-hexose benchmark systems.
2. Correct every record that describes `ROH` as a charge-balancing counterion.
3. Compare the host formula, six-urea/three-amide motif inventory and graph fingerprint against the published truncated GluHUT structure and synthesis.
4. Map the deposited SMD restraint atom groups to exact host and guest atoms.
5. Continue requesting the source RESP/GAFF/GLYCAM package and starting coordinates.
6. Keep PaG parameterization and MC-Q1 ranking disabled until source-graph and parent-calibration gates pass.

## Reproducibility

- Script: `research/audit_gluhut_graph_chemistry.py`.
- Machine-readable output: `research/data/gluhut-graph-chemistry-audit.json`.
- Tests: `tests/test_audit_gluhut_graph_chemistry.py`.
- Workflow: `.github/workflows/gluhut-graph-chemistry-audit.yml`.

## Evidence boundary

This is a graph and composition audit of deposited public benchmark systems. It is not an affinity calculation, NNC2215 result, PaG result, MC-Q1 result, efficacy, safety or working-drug result.
