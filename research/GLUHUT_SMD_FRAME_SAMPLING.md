# GluHUT SMD frame and time-provenance audit

**Date:** 4 November 2026  
**Status:** public coordinate files are downsampled or otherwise incomplete; exact frame-time mapping is absent  
**Scope:** deposited AMBER SMD control and twenty public BINPOS files

## Question

The deposited coordinate files were described and named as visualisation trajectories. Before using their frame sequence to reconstruct pulling speed, target history or nonequilibrium work, the requested AMBER output cadence was compared with the actual deposited frame count.

## Deposited AMBER control

`SMD.in` specifies:

- `nstlim = 15,000,000` integration steps;
- `dt = 0.002 ps`;
- total production time: **30 ns**;
- `ntwx = 1,000` steps;
- `ntwr = 1,000` steps;
- `ntpr = 1,000` steps.

At the requested coordinate-write cadence, a complete production trajectory would contain:

`15,000,000 / 1,000 = 15,000 coordinate outputs`

For a reported 15 Å target change, the nominal target increment between requested coordinate outputs is:

`15 Å / 15,000 = 0.001 Å per requested output`

## Deposited coordinate files

Every public BINPOS file contains exactly:

- 1,000 frames;
- 183 atoms per frame.

Across twenty files:

- total deposited frames: 20,000;
- requested-output to deposited-frame ratio: **15:1**;
- fifteen filenames explicitly contain `visualisation`;
- no frame-selection script is deposited;
- no per-frame time mapping is deposited;
- no applied target-versus-time series is deposited.

The 15-fold ratio is consistent with retaining approximately every fifteenth requested output, but the archive does not prove that this is how the files were generated. The first or last frames may also have been selected or transformed separately.

## Coordinate progress

A trajectory was classified as consistent with the paper’s approximately 15 Å target change when its reconstructed net COM-coordinate change was within 0.5 Å of 15 Å.

### Consistent files: 13

- glucose: 5 of 5;
- fructofuranose: 5 of 5;
- galactose path B: 3 of 3.

Mean reconstructed net changes:

- glucose: `14.884 Å`;
- fructofuranose: `14.904 Å`;
- galactose path B: approximately 15 Å.

The mean coordinate increment per deposited frame is close to `0.015 Å`, which is numerically consistent with a 15-fold reduction from a nominal `0.001 Å` target increment. This remains an inference, not deposited provenance.

### Inconsistent files: 7

- fructopyranose: 5 of 5, mean net change `19.908 Å`;
- duplicate galactose path A: 2 files, net change `21.396 Å`.

These files cannot result from the same simple 0-to-15 Å target history unless additional undocumented processing, longer target distances or different production controls were used.

## Consequences

The public BINPOS files are sufficient for:

- coordinate visualisation;
- static and frame-order continuity checks;
- reconstructing the deposited COM-coordinate sequence;
- post hoc comparison of alternative atom-group definitions.

They are not sufficient for:

- assigning exact simulation time to each deposited frame;
- reconstructing the applied target position at each deposited frame;
- calculating pulling velocity from file cadence alone;
- reconstructing nonequilibrium work;
- independently regenerating Jarzynski or Crooks estimates;
- establishing which source frames were retained;
- treating endpoint differences as solely molecular response.

## Decision

1. Treat every public BINPOS file as a visualisation or reduced coordinate record unless full frame-selection provenance is supplied.
2. Do not infer a 30 ps deposited frame spacing merely from the 15-fold frame-count ratio.
3. Request the original full trajectories or exact downsampling script and frame-time mapping.
4. Request target histories and production controls specifically for fructopyranose and galactose path A.
5. Do not reconstruct work or corrected PMFs from the public visualisation files.
6. Keep parent glucose-versus-PaG calibration and MC-Q1 ranking blocked.

## Reproducibility

- Audit script: `research/audit_gluhut_smd_frame_sampling.py`.
- Machine-readable result: `research/data/gluhut-smd-frame-sampling.json`.
- Tests: `tests/test_audit_gluhut_smd_frame_sampling.py`.
- Workflow: `.github/workflows/gluhut-smd-frame-sampling.yml`.
- Coordinate input: `research/data/gluhut-binpos-coordinate-audit.json`.
- Archive inventory input: `research/data/gluhut-smd-archive-inventory.json`.

## Evidence boundary

This is frame-count and protocol-provenance evidence. It is not SMD work, PMF, NNC2215 affinity, PaG affinity, MC-Q1 affinity, efficacy, safety or a working-drug result.
