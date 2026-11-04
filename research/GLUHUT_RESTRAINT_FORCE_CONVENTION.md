# GluHUT restraint-force convention

**Date:** 4 November 2026  
**Status:** AMBER, article and Grossfield WHAM coefficient conventions reconciled  
**Scope:** public GluHUT SMD and umbrella restraint notation

## The ambiguity

The deposited AMBER restraint inputs and the paper report numerical force values that differ from the spring constants used in Grossfield WHAM and many other simulation packages.

The difference is a notation convention, not a physical disagreement.

AMBER distance restraints use:

`E = rk2 × (R − R0)²`

The conventional harmonic form is:

`E = ½ × K × (R − R0)²`

Therefore:

`K = 2 × rk2`

## Resolved values

| Use | AMBER `rk2` | Conventional/half-K package `K` | Grossfield WHAM metadata |
|---|---:|---:|---:|
| steered MD | 200 kcal mol⁻¹ Å⁻² | 400 kcal mol⁻¹ Å⁻² | not applicable |
| umbrella sampling | 50 kcal mol⁻¹ Å⁻² | 100 kcal mol⁻¹ Å⁻² | 100 kcal mol⁻¹ Å⁻² |

The deposited SMD controls contain:

- `rk2 = 200`;
- initial target `r2 = 0 Å`;
- final target `r2a = 15 Å`.

These match the numeric SMD values and range reported in the article under the AMBER convention.

The paper reports an umbrella value of 50. The deposited Grossfield WHAM metadata uses 100. These are the same physical restraint because Grossfield WHAM uses the conventional `½K` form.

## Energy equivalence

For a displacement of 1 Å:

- AMBER SMD: `200 × 1² = 200 kcal/mol`;
- conventional SMD: `½ × 400 × 1² = 200 kcal/mol`;
- AMBER umbrella: `50 × 1² = 50 kcal/mol`;
- conventional/WHAM umbrella: `½ × 100 × 1² = 50 kcal/mol`.

The audit verifies exact equivalence at 0.1, 0.5 and 1.0 Å.

## Portable manifest rule

Every future run manifest must record both values and the energy form explicitly.

Use:

- AMBER `&rst` SMD: `rk2 = 200`;
- AMBER `&rst` umbrella: `rk2 = 50`;
- OpenMM or another package using `½Kx²`: SMD `K = 400`;
- OpenMM or another package using `½Kx²`: umbrella `K = 100`;
- Grossfield WHAM umbrella metadata: `K = 100`.

A bare number such as “force constant 200” is insufficient without naming the energy convention.

## Confidence effect

The factor-of-two discrepancy between the article’s umbrella value and deposited WHAM metadata is resolved. It is not evidence of a mismatched simulation.

This removes a unit/convention blocker for future reproduction. It does not remove the remaining input blockers:

- the fructofuranose guest group includes unexplained hydrogen `H62`;
- source starting coordinates are missing;
- source RESP/GAFF/GLYCAM generation files are missing;
- complete parent glucose-versus-PaG calibration has not been run.

## Decision

1. Mark the AMBER-to-conventional restraint conversion as resolved.
2. Apply the 2× conversion in every future OpenMM, WHAM and cross-engine manifest.
3. Do not generate new SMD or umbrella paths until the guest atom groups and source coordinates are resolved.
4. Keep PaG calibration and MC-Q1 ranking disabled.

## Reproducibility

- Script: `research/audit_gluhut_restraint_force_convention.py`.
- Machine result: `research/data/gluhut-restraint-force-convention.json`.
- Tests: `tests/test_audit_gluhut_restraint_force_convention.py`.
- Workflow: `.github/workflows/gluhut-restraint-force-convention.yml`.
- Public evidence PR: #41.

## Sources

- Eades et al., DOI `10.1039/D6SC00738D`.
- AMBER Users' Manual, steered-MD and NMR restraint sections.
- Grossfield WHAM documentation and the deposited GluHUT metadata.

## Evidence boundary

This resolves a coefficient convention. It does not validate atom selections, force-field chemistry, sampling, PMFs, NNC2215, PaG, MC-Q1, efficacy, safety or a working drug.
