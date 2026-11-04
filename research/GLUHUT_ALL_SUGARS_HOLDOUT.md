# GluHUT all-sugar WHAM hold-out reproduction

**Date:** 4 November 2026  
**Status:** Grossfield hold-out reproduction `PASS WITH ANOMALY REVIEW`; independent NumPy and sensitivity checks remain open  
**Scope:** public GluHUT host–sugar methodological calibration only

## Result

The checksum-verified Zenodo `WHAM_PMF.zip` archive was reprocessed without sugar-specific retuning using Grossfield WHAM 2.0.11 from pinned source commit `43bdd6fcd1916ae2b169cb52ac326b7eea567abb`.

| Sugar/path group | Reproduced kcal/mol | Published kcal/mol | Difference | Within published 1 SD |
|---|---:|---:|---:|---:|
| D-glucose | −12.782 ± 1.017 | −12.5 ± 0.9 | −0.282 | yes |
| D-galactose path B | −8.441 ± 0.968 | −8.1 ± 1.9 | −0.341 | yes |
| D-galactose path A | −6.611, one route | −6.9 ± 0.3 | +0.289 | yes |
| D-fructopyranose | −4.970 ± 1.518 | −5.0 ± 1.3 | +0.030 | yes |
| D-fructofuranose | −7.594 ± 0.488 | −7.7 ± 0.5 | +0.106 | yes |

The reproduced strongest-to-weakest order is identical to the published order:

1. glucose;
2. galactose path B;
3. fructofuranose;
4. galactose path A;
5. fructopyranose.

The output is in:

- `research/data/gluhut-all-sugars-wham-reproduction.json`;
- `research/data/gluhut-all-sugars-wham-summary.csv`.

## What passed

- deposited MD5 and previously recorded SHA-256 verified;
- archive CRC validation passed;
- pinned Grossfield WHAM built and executed;
- all expected glucose, galactose, fructopyranose and fructofuranose route groups were found;
- all group means fall within one published standard deviation;
- the qualitative affinity ranking reproduced without protocol retuning;
- the evidence boundary remained explicit: no NNC2215, PaG or MC-Q1 inference was made.

This closes the Grossfield lower-affinity hold-out component of AWS-1b.

## Corrected archive anomaly

The archive audit identified one severe length outlier:

`WHAM_PMF/histogram-collection-fructofuranose-US-2/12-6-test-dist_vs_t.dat`

- centre: 12.6 Å;
- frames: 205;
- size: 4,100 bytes;
- size fraction of route median: 0.0137.

The earlier browser-preview note described this as fructopyranose. The checksum-verified archive establishes that it belongs to **fructofuranose US-2**.

The baseline WHAM result naturally gives the short window much lower sample weight, but that does not by itself establish adequate overlap around 12.6 Å. A targeted sensitivity analysis must compare:

1. the deposited baseline;
2. the route with the 12.6 Å window excluded;
3. overlap around 12.3–12.9 Å;
4. reasonable trajectory-discard and bin-count alternatives.

## 51-versus-52 preview discrepancy

The actual analysis discovered exactly 51 matched time-series windows per standard route. The browser preview's count of 52 members therefore included a non-window member or directory artifact. The scientific analysis did not silently use 52 windows.

## Confidence decision

The result materially increases confidence in the **public GluHUT analysis pipeline**. It demonstrates that the same deposited workflow reproduces both the glucose result and lower-affinity hold-outs.

It does not increase MC-Q1 above E2 because:

- PaG has not been calibrated;
- the free receptor is not the complete NNC2215 conjugate;
- the independent NumPy solver is still running;
- the fructofuranose short-window sensitivity is unresolved;
- no MC-Q1 material or direct binding exists.

## Remaining gates before MC-Q1 ranking

1. independent NumPy WHAM agreement on the three glucose routes;
2. targeted short-window and bin/discard sensitivity;
3. exact parent free-receptor and PaG chemistry/force-field gate;
4. source-matched parent glucose-versus-PaG calibration;
5. prospective, preregistered MC-P0/MC-Q1/MC-Q1D/MC-P3/MC-P1 comparison.

Candidate ranking remains unauthorized.
