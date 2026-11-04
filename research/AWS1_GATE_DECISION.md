# AWS-1 final gate decision

**Decision date:** 4 November 2026  
**Gate:** reproduce the deposited GluHUT WHAM/PMF analysis before any prospective Dia Health One candidate calculation  
**AWS spend:** `$0`  
**Overall verdict:** **PASS for the deposited WHAM-analysis scope, with mandatory qualifications**  
**Raw trajectory-to-histogram lineage subgate:** **FAIL**  
**MC-Q1 ranking:** not authorized

## Why the overall gate passes

AWS-1 asked whether the public GluHUT analysis could be reproduced from the deposited WHAM inputs and whether the qualitative sugar ranking survived a protocol applied without sugar-specific retuning.

That bounded question has been answered positively:

- the deposited archives and selected inputs match the published checksums and byte sizes;
- archive extraction, input parsing and generated records are deterministic;
- the deposited Grossfield WHAM workflow reproduces the public glucose PMFs;
- all three glucose routes are interpretable;
- an independent NumPy WHAM implementation agrees with the three glucose routes within a maximum absolute difference of `0.000454 kcal/mol`;
- the published qualitative strongest-to-weakest sugar order is reproduced without sugar-specific retuning;
- reasonable discard and bin-count variations do not reverse the qualitative ranking;
- raw deposited WHAM inputs, generated outputs, code, hashes, uncertainty and negative findings are retained.

This is sufficient to close the original deposited-analysis gate as `PASS`.

## Pass-criteria assessment

| AWS-1 criterion | Decision | Evidence boundary |
|---|---|---|
| Deposited checksums and byte sizes verify | PASS | Zenodo manifest, MD5, SHA-256 and ZIP-CRC records |
| Extraction and parsing are deterministic | PASS | Safe-path audit, content-addressed extraction and reproducibility tests |
| One glucose PMF is regenerated | PASS | Deposited Grossfield workflow |
| Three glucose routes provide interpretable uncertainty | PASS | Three route outputs plus independent NumPy cross-check |
| Hold-out sugar ranking agrees with publication | PASS | Glucose, galactose, fructopyranose and fructofuranose ordering reproduced |
| Reasonable analysis choices do not reverse ranking | PASS WITH QUALIFICATION | Qualitative ranking survives; fructofuranose US-2 route robustness fails |
| Inputs, outputs and limitations are retained | PASS | Public machine records, scripts, tests and scientific notes |

## Mandatory qualification: fructofuranose US-2

Fructofuranose US-2 contains a 12.6 Å window with only 205 samples, versus approximately 15,001 in neighbouring deposited WHAM windows.

- deposited route estimate: `−7.521 kcal/mol`;
- omit the 12.6 Å window: `−6.815 kcal/mol`;
- change: `+0.707 kcal/mol`;
- preregistered route-robustness threshold: `0.50 kcal/mol`;
- direct 12.3-to-12.9 Å overlap after omission: `0.000133`.

Therefore the **fructofuranose US-2 route robustness subgate fails**. The complete sugar ranking remains usable only when this route-specific limitation is attached. The anomaly must not be hidden, retrospectively repaired or treated as an equally robust independent path.

## Separate failure: raw trajectory lineage

The public trajectory archives do not provide a complete, authoritative route from raw trajectory files to the deposited 15,001-sample WHAM time series.

The audit found, among other issues:

- 714 umbrella trajectory files across 14 raw route directories but only 13 unique complete route signatures;
- exact duplicate complete galactose path-A routes;
- duplicate galactose path-B 0.0 and 0.3 Å files;
- a galactose 2.25 Å WHAM window with no matching named public trajectory;
- fructofuranose route files containing one or ten frames where the deposited WHAM inputs contain hundreds or approximately 15,001 samples;
- no authoritative frame-selection, time map or route-to-WHAM lineage.

Consequently:

- deposited histogram and WHAM reanalysis is reproducible;
- regeneration of the deposited WHAM inputs from the public trajectory archive is **not** reproducible;
- new SMD or umbrella generation from the deposited controls is not authorized;
- the raw-lineage failure does not convert the narrower deposited-WHAM reproduction into a failure, but it blocks prospective model generation.

## Chemistry and coordinate limitations discovered during AWS-1

The calibration work also established:

- one identical compiled 159-atom host is used across all four sugars;
- compiled host formula: `C60H75N15O9`;
- six urea and three amide carbonyls;
- the compiled graph matches the published truncated receptor architecture;
- the two-record carbohydrate plus `ROH` construction is one neutral covalent `C6H12O6` guest;
- the source RESP/GAFF/GLYCAM derivation and source starting coordinates are absent;
- the host SMD reference group is consistent;
- the sugar SMD groups are not chemically uniform;
- fructopyranose and fructofuranose reuse the same numeric indices despite different atom identities, including carbon `C5` becoming hydrogen `H62` at index 165;
- no corrected fructofuranose restraint can be selected without author clarification or a preregistered prospective sensitivity study.

These findings increase confidence in the evidence boundary. They do not authorize parameter transfer.

## What AWS-1 does and does not establish

AWS-1 establishes that:

- the deposited WHAM analysis is reproducible;
- the public sugar ranking is a valid calibration target at the deposited-data-product level;
- the common compiled host and guest graph can be audited consistently;
- known route and provenance weaknesses are explicit.

AWS-1 does **not** establish that:

- the source RESP/GAFF/GLYCAM derivation is reproducible;
- the raw public trajectories regenerate the deposited WHAM inputs;
- a source-matched PaG model exists;
- the complete NNC2215 conjugate is represented;
- MC-Q1 preserves glucose affinity or rejects PaG;
- any Dia Health One compound is a hit, switch, drug or therapeutic candidate.

## Transition to the next gate

The next gate is not a large candidate screen. It is a source-controlled parent calibration readiness gate.

Before any prospective glucose-versus-PaG free-energy calculation, require:

1. an exact parent receptor chemical graph and atom mapping;
2. source or independently regenerated host charges and nonstandard parameters;
3. reviewed glucose and PaG parameterization, protonation and ion assumptions;
4. source-controlled starting coordinates or a preregistered pose-generation protocol;
5. one chemically uniform reaction coordinate;
6. multiple independent paths with overlap and convergence diagnostics;
7. a predeclared uncertainty and stop rule;
8. explicit separation of the free-receptor/model observable from the complete NNC2215 experimental observable.

Until those conditions pass, AWS candidate spend remains paused and MC-Q1 remains at evidence level E2.

## Final statement

**AWS-1 is closed as PASS for reproducibility of the deposited WHAM analysis, with a failed fructofuranose route-robustness subgate and a failed raw trajectory-lineage subgate. No prospective affinity or drug claim follows from this verdict.**
