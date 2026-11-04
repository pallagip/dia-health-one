# GluHUT fructofuranose short-window sensitivity

**Date:** 4 November 2026  
**Status:** prespecified sensitivity threshold `FAIL`; calibration remains useful with a route-specific limitation  
**Scope:** deposited public GluHUT fructofuranose US-2 data only

## Question

The deposited fructofuranose US-2 umbrella set contains one severe length outlier:

`WHAM_PMF/histogram-collection-fructofuranose-US-2/12-6-test-dist_vs_t.dat`

- umbrella centre: 12.6 Å;
- frames: 205;
- neighbouring and median window frames: approximately 15,001;
- file size: 4,100 bytes.

The preregistered test asked whether reasonable handling of this window changes the route-level binding free-energy estimate by more than 0.50 kcal/mol.

## Result

| Analysis | Binding free energy, kcal/mol | Change from deposited baseline |
|---|---:|---:|
| Deposited baseline, 400 bins | −7.521 | 0.000 |
| Exclude 12.6 Å window | −6.815 | **+0.707** |
| Discard first 10% of all windows | −7.436 | +0.085 |
| Discard first 20% of all windows | −7.775 | −0.254 |
| 200 bins | −7.246 | +0.275 |
| 800 bins | −7.770 | −0.248 |

All runs completed, but exclusion of the short window exceeded the preregistered 0.50 kcal/mol threshold. The sensitivity decision is therefore **FAIL**, not pass.

## Local overlap

Bhattacharyya histogram overlap was:

- 12.3 ↔ 12.6 Å: 0.0944;
- 12.6 ↔ 12.9 Å: 0.1371;
- direct 12.3 ↔ 12.9 Å after excluding the short window: 0.000133.

The 12.6 Å window is poorly sampled, but it bridges an otherwise nearly disconnected 12.3-to-12.9 Å interval. Its removal therefore produces a material route-level change.

## Interpretation

This does not invalidate the complete public GluHUT calibration:

- the deposited Grossfield workflow reproduced all published sugar-group means within one published standard deviation;
- the qualitative strongest-to-weakest sugar order was reproduced without sugar-specific retuning;
- the independent NumPy WHAM implementation agreed with the three glucose routes within 0.000454 kcal/mol.

It does qualify the fructofuranose US-2 route. That route must not be treated as an equally robust independent path without either:

1. additional sampling around 12.3–12.9 Å; or
2. a prespecified uncertainty treatment that includes the 0.707 kcal/mol omission sensitivity.

The overall fructofuranose group and sugar ranking should be reported with this anomaly attached, not silently corrected after inspection.

## Confidence decision

- Public WHAM analysis implementation: supported.
- Glucose route reproduction: supported independently.
- Lower-affinity qualitative ordering: supported.
- Fructofuranose US-2 route robustness: **not supported under the preregistered 0.50 kcal/mol criterion**.
- MC-Q1 candidate ranking: not authorized.
- PaG calibration: not performed.

## Reproducibility

- Archive MD5: `dbc6723b8adcb85deb8230f63b6a51da`.
- Archive SHA-256: `e0e7b7c04326b8fe86cf71e2193315f34fb865168cf11117b3caf9253119c6eb`.
- Grossfield WHAM source commit: `43bdd6fcd1916ae2b169cb52ac326b7eea567abb`.
- Script: `research/audit_gluhut_short_window_sensitivity.py`.
- Machine-readable result: `research/data/gluhut-fructofuranose-us2-sensitivity.json`.

## Evidence boundary

This is a methodological sensitivity result for a public synthetic host–sugar dataset. It is not NNC2215, PaG, MC-Q1, insulin-receptor activity, efficacy, safety, or a drug result.
