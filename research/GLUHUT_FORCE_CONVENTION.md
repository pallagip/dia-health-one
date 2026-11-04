# GluHUT restraint-force convention audit

**Date:** 4 November 2026  
**Status:** scalar AMBER/Grossfield convention resolved; trajectory-generation gates remain blocked  
**Scope:** published and deposited GluHUT SMD and umbrella restraint coefficients

## Question

The GluHUT paper reports:

- steered MD with a harmonic restraint of **200 kcal mol⁻¹ Å⁻²**;
- umbrella sampling with a force constant of **50 kcal mol⁻¹ Å⁻²**.

The deposited SMD files contain `rk2 = 200`, while the deposited/reproduced Grossfield WHAM metadata use a force constant of `100.0`. These numbers appear inconsistent unless the AMBER and conventional harmonic definitions are separated.

## AMBER functional form

The AMBER distance-restraint manual defines the parabolic region as:

`E = rk2 × (R − R0)²`

The usual harmonic notation is:

`E = ½ K × (R − R0)²`

Therefore:

`K = 2 × rk2`

The AMBER umbrella-sampling documentation explicitly notes that the conventionally defined force constant is twice `rk2`. Grossfield WHAM uses the conventional `K` coefficient.

## Steered MD

The deposited GluHUT SMD inputs use:

`rk2 = 200 kcal mol⁻¹ Å⁻²`

This matches the paper’s reported numeric SMD value of 200 when the paper value is interpreted as the AMBER input coefficient.

The equivalent coefficient in standard `½Kx²` notation is:

`K = 400 kcal mol⁻¹ Å⁻²`

Example penalties for the deposited potential:

| Displacement | AMBER `rk2 × ΔR²` | Standard `½K × ΔR²` |
|---:|---:|---:|
| 0.1 Å | 2 kcal mol⁻¹ | 2 kcal mol⁻¹ |
| 0.5 Å | 50 kcal mol⁻¹ | 50 kcal mol⁻¹ |
| 1.0 Å | 200 kcal mol⁻¹ | 200 kcal mol⁻¹ |

Thus there is no scalar mismatch between the paper’s SMD number and the deposited `rk2`. The nomenclature is ambiguous because “force constant” can refer either to the AMBER coefficient or to conventional `K`.

## Umbrella sampling and WHAM

The paper reports an umbrella value of:

`50 kcal mol⁻¹ Å⁻²`

The independently reproduced WHAM implementation uses:

`K = 100 kcal mol⁻¹ Å⁻²`

This is exactly the expected conversion:

`K_WHAM = 2 × rk2_AMBER = 2 × 50 = 100`

The deposited WHAM metadata are therefore consistent with the paper and AMBER’s restraint-energy convention.

## Decision

The scalar force-constant convention is resolved:

- paper SMD 200 ↔ deposited AMBER `rk2 = 200`;
- equivalent conventional SMD `K = 400`;
- paper umbrella 50 ↔ AMBER `rk2 = 50`;
- Grossfield WHAM `K = 100`.

This closes the force-coefficient conversion question. It does **not** authorize prospective trajectory generation because independent blockers remain:

- guest reaction-coordinate definitions differ across sugars;
- fructofuranose contains unexplained selected hydrogen `H62`;
- source starting coordinates are absent;
- source RESP/GAFF/GLYCAM generation files are absent;
- the exact compiled host graph is not yet fully reconciled to the source synthesis;
- parent glucose-versus-PaG calibration has not passed.

## Program consequence

The authors no longer need to answer whether the numeric factor-of-two between paper umbrella value 50 and WHAM metadata 100 is an error. It is explained by AMBER versus conventional harmonic notation.

Author clarification is still required for:

- the fructofuranose atom selection;
- the nonuniform guest COM definitions;
- source coordinates and assembly scripts;
- parameter-generation provenance.

## Reproducibility

- Script: `research/audit_gluhut_force_convention.py`.
- Machine-readable output: `research/data/gluhut-force-convention.json`.
- Tests: `tests/test_audit_gluhut_force_convention.py`.
- Workflow: `.github/workflows/gluhut-force-convention.yml`.
- SMD source record: `research/data/gluhut-smd-restraint-mapping.json`.
- WHAM source record: `research/data/gluhut-glucose-numpy-wham-crosscheck.json`.

## Primary references

- Eades et al., *Free energy simulations on a biomimetic glucose receptor: understanding the selectivity of GluHUT*, Chemical Science 2026, DOI `10.1039/D6SC00738D`.
- AMBER Users’ Manual, distance-restraint parabolic functional form and umbrella-sampling force-constant convention.
- AMBER-hub `rst` documentation for the meaning of `rk2` and `rk3`.

## Evidence boundary

This resolves a scalar software convention. It is not a binding free-energy result, NNC2215 result, PaG result, MC-Q1 result, efficacy, safety or working-drug result.
