# GluHUT SMD coordinate and endpoint provenance

**Date:** 4 November 2026  
**Status:** all deposited coordinates recovered; galactose endpoint heterogeneity and an exact duplicate detected  
**Scope:** twenty public GluHUT SMD BINPOS files

## Result

All twenty deposited Scripps BINPOS files were decoded with a small audited native reader implementing the official CPPTRAJ layout:

- four-byte `fxyz` magic;
- one integer atom count before every frame;
- `3 × natom` single-precision coordinates per frame.

Every file contains:

- 183 atoms per frame, matching its AMBER topology;
- 1,000 recorded coordinate frames;
- finite coordinates;
- plausible first-frame bond lengths;
- a forward host-to-guest centre-of-mass separation path.

Total recovered coordinate frames: **20,000**.

The first recorded frame is not assumed to be the original source restart. It may follow unrecorded minimisation, equilibration or integration steps.

## Fructofuranose guest-group alternatives measured on the deposited paths

The deposited fructofuranose coordinate uses the five-member ring plus isolated hydrogen `H62`. Six graph-coherent all-heavy alternatives were measured post hoc on the same deposited trajectories.

The two highest-information alternatives are:

| Guest-group definition | Mean absolute COM-distance change across runs | Mean of each run's maximum absolute change |
|---|---:|---:|
| five ring heavy atoms only | **0.0260 Å** | **0.0501 Å** |
| replace H62 with its heavy parent C6 | 0.2225 Å | 0.3843 Å |

For the five-ring-only definition, the mean last-frame difference is 0.0082 Å, with individual-run last-frame differences from −0.0364 to +0.0427 Å.

This makes the five-ring-only interpretation geometrically closest to the deposited H62 coordinate on the already generated trajectories. It does **not** establish that five ring atoms were the intended restraint. A different group would change the applied centre of mass during trajectory generation, and post hoc remeasurement cannot reproduce the corresponding work distribution or path.

No corrected restraint or trajectory is selected.

## Endpoint protocols recovered from the files

Endpoints were clustered using a preregistered maximum within-cluster diameter of 0.25 Å.

### Consistent deposited systems

| System | Files | Distinct files | Endpoint cluster |
|---|---:|---:|---:|
| glucose | 5 | 5 | approximately 14.974 Å |
| fructofuranose | 5 | 5 | mean 14.977 Å; range 14.910–15.068 Å |
| fructopyranose | 5 | 5 | mean 20.002 Å; range 19.957–20.067 Å |

The 20 Å fructopyranose endpoint differs from the 15 Å target encoded in the deposited common restraint input. That discrepancy requires protocol provenance; an endpoint alone does not reveal the applied target history or production duration.

### Galactose is heterogeneous

Galactose does not contain five independent trajectories under one endpoint protocol.

Runs 1–3 are distinct `*-B.binpos` files and form one approximately 15 Å cluster:

- mean endpoint: 14.999 Å;
- range: 14.943–15.048 Å;
- distinct coordinate series: 3.

Runs 4 and 5 are `*-A.binpos` files with endpoint 22.126 Å. They are exact duplicates:

- byte-for-byte file SHA-256: `7a689d6eadd31bfaa3b5689d579f698066544494b254a2858911fa800bcc1903`;
- coordinate-series SHA-256: `2190f8472e49e407c95818f89599d08835562b6e12dfaafc73b3c9f13c348f57`;
- identical 1,000-frame coordinate series;
- unique coordinate series represented by runs 4 and 5: **1**, not 2.

The archive therefore contains:

- three distinct approximately 15 Å galactose trajectories;
- two copies of one approximately 22.1 Å trajectory;
- at least two different endpoint protocols or file-generation histories.

The five files must not be pooled as five independent replicates.

## What this does and does not affect

This provenance result applies to the deposited SMD files. It does not automatically invalidate the separately deposited umbrella histograms or their deterministic WHAM reanalysis. A direct file-to-window lineage has not been established.

It does establish that:

- the deposited SMD archive is not a uniform five-run-per-sugar dataset;
- galactose has only three distinct files in the approximately 15 Å cluster;
- the approximately 22.1 Å galactose record has only one unique trajectory, duplicated under two run directories;
- fructopyranose was pulled to approximately 20 Å rather than the approximately 15 Å endpoint seen for glucose and fructofuranose;
- the single deposited restraint file is insufficient to explain all recorded endpoints.

## Decision

1. Preserve file names, A/B labels, endpoint clusters and hashes in every analysis.
2. Do not treat galactose runs 4 and 5 as independent observations.
3. Do not pool the approximately 15 Å and approximately 22.1 Å galactose files under one protocol.
4. Request the intended target history, production length, seeds and file-generation lineage for every run.
5. Keep prospective SMD and umbrella generation blocked until the source protocols and corrected atom groups are established.
6. Keep parent glucose-versus-PaG calibration and MC-Q1 ranking blocked.
7. Do not retrospectively alter deposited WHAM outputs based on this SMD-only provenance result.

## External clarification required

The author request now includes:

- why fructopyranose reaches approximately 20 Å while its deposited restraint file contains `r2a = 15 Å`;
- why galactose contains distinct A and B files with approximately 22.1 Å and 15 Å endpoints;
- whether galactose runs 4 and 5 were intentionally duplicated;
- which trajectories seeded the deposited umbrella paths;
- exact run scripts, restart files, target schedules, random seeds and trajectory-generation lineage.

## Reproducibility

- Native BINPOS reader: `research/binpos_native_reader.py`.
- Coordinate audit: `research/audit_gluhut_binpos_coordinates.py`.
- Coordinate result: `research/data/gluhut-binpos-coordinate-audit.json`.
- Endpoint summarizer: `research/summarize_gluhut_smd_endpoint_provenance.py`.
- Endpoint result: `research/data/gluhut-smd-endpoint-provenance.json`.
- Public evidence PRs: #46 and #47.
- Archive MD5: `9f3c6733dc608d2bd37c488c560a5443`.
- Archive SHA-256: `04c58666328764aeb47bac0d9f8607546f9763892cb9a333b108630df70e3943`.

## Evidence boundary

This is coordinate and file provenance for public SMD data. It is not SMD work, an umbrella PMF, a binding affinity, an NNC2215 result, a PaG result, an MC-Q1 result, efficacy, safety or a working drug.
