# GluHUT SMD coordinate and file-provenance audit

**Date:** 4 November 2026  
**Status:** all deposited coordinate files parsed; fructofuranose and galactose provenance limitations retained  
**Scope:** public GluHUT SMD BINPOS trajectories only

## Purpose

The deposited restraint files define atom selections but do not show how the resulting centre-of-mass coordinate behaved in the actual SMD trajectories. The twenty public BINPOS files were therefore read directly, mapped to the compiled AMBER topologies, and analysed frame by frame.

This audit answers three questions:

1. Can all deposited binary coordinate files be read reproducibly?
2. How much would the bounded fructofuranose atom-group alternatives change the measured coordinate on the already deposited trajectories?
3. Do the deposited files form one consistent five-run SMD protocol per sugar?

## Native BINPOS reader

MDTraj 1.11.2 no longer exposes a BINPOS reader. A minimal local reader was implemented from the historical Scripps/VMD BINPOS format:

- four-byte ASCII magic `fxyz`;
- repeated int32 atom count;
- `3 × natom` float32 coordinates in Å per frame;
- little- and big-endian support;
- frame-level atom-count checks;
- exact-length and truncation checks;
- non-finite-coordinate rejection.

The reader was validated with round-trip little- and big-endian fixtures and malformed-file tests.

## Complete coordinate acquisition

All twenty public trajectories were read successfully:

- four sugar systems;
- five deposited files per sugar;
- 1,000 frames per file;
- 20,000 total frames;
- 183 atoms per frame;
- topology atom counts matched every file;
- median first-frame compiled bond lengths passed the Å-unit sanity check.

The deposited mass-weighted host-to-guest centre-of-mass coordinate was reconstructed for every frame using the compiled AMBER masses and the deposited `igr1`/`igr2` groups.

## Fructofuranose H62 alternatives on deposited coordinates

The deposited fructofuranose group contains the five furanose-ring heavy atoms plus hydrogen `H62`. Six graph-coherent all-heavy alternatives were evaluated on the same deposited coordinate frames.

These are **post hoc coordinate measurements** on trajectories generated with the H62 restraint. They are not corrected trajectories and do not recover corrected SMD work or a corrected PMF.

### Highest-information alternatives

| Alternative | Mean absolute coordinate difference across five runs | Mean maximum absolute difference |
|---|---:|---:|
| Five ring heavy atoms only | **0.0260 Å** | **0.0501 Å** |
| Replace H62 with parent C6 atom 164 | 0.2225 Å | 0.3843 Å |

Other bounded substitutions changed the post hoc coordinate more:

- C1 atom 177: approximately 0.322 Å mean absolute difference;
- O4 atom 171: approximately 0.307 Å;
- O3 atom 175: approximately 0.306 Å;
- terminal O1 atom 183: approximately 0.429 Å.

### Interpretation

Deleting H62 and using the five ring heavy atoms gives a coordinate numerically closest to the deposited H62 definition on the already generated trajectories. This makes the five-ring-only hypothesis the least disruptive post hoc coordinate interpretation.

It does **not** establish author intent or show that a five-ring restraint would generate the same trajectory. The applied centre of mass, mass distribution, pulling force and nonequilibrium work would differ during a newly generated simulation.

No corrected restraint is selected automatically.

## Reconstructed endpoints

### Systems with one deposited endpoint cluster

| System | Files | Mean final coordinate |
|---|---:|---:|
| Glucose | 5 distinct | 14.974 Å |
| Fructofuranose | 5 distinct | 14.977 Å |
| Fructopyranose | 5 distinct | 20.002 Å |

The fructopyranose endpoint is near 20 Å even though the article describes SMD pulls from 0 to 15 Å. This remains a source-provenance question.

### Galactose file heterogeneity

The five galactose paths are not one set of five independent trajectories under one endpoint protocol.

#### Path B

Runs 1–3 are distinct files:

- endpoint range: 14.943–15.048 Å;
- mean endpoint: 14.999 Å;
- file and coordinate hashes are unique.

#### Path A

Runs 4 and 5 are byte-for-byte identical copies:

- file SHA-256: `7a689d6eadd31bfaa3b5689d579f698066544494b254a2858911fa800bcc1903`;
- coordinate-series SHA-256: `2190f8472e49e407c95818f89599d08835562b6e12dfaafc73b3c9f13c348f57`;
- first coordinate: 0.730 Å;
- final coordinate: 22.126 Å;
- 1,000 identical frames in each deposited copy.

The archive therefore contains:

- three distinct galactose path-B trajectories near 15 Å;
- one unique path-A trajectory near 22.1 Å, deposited twice.

It does not contain five independent galactose SMD trajectories.

## Periodic-imaging test

The anomalous 22.126 Å galactose path-A endpoint was tested for a single box-scale coordinate jump.

Result:

- no frame-to-frame step exceeded the conservative 5 Å discontinuity threshold;
- maximum absolute step: 2.672 Å;
- the final-quarter mean coordinate was 18.257 Å;
- the coordinate progressed to the long endpoint rather than showing one periodic-imaging jump.

Periodic imaging is therefore not supported as the explanation for the endpoint. The data support a different-target, wrong, or mislabeled deposited path.

No automatic unwrapping or correction is justified because BINPOS contains no box vectors and the trajectory does not show a box-scale single-frame jump.

## Consequences for confidence

Supported:

- deterministic reading of all twenty deposited BINPOS coordinate files;
- reconstruction of the deposited mass-weighted reaction coordinate;
- quantitative post hoc comparison of the fructofuranose alternatives;
- exact duplicate detection and endpoint clustering;
- rejection of a simple periodic-imaging explanation for the anomalous galactose path-A file.

Not supported:

- treating five galactose files as five independent runs;
- treating all deposited sugars as one uniform 0-to-15 Å trajectory protocol;
- retrospectively correcting H62, galactose A or fructopyranose trajectories;
- reconstructing nonequilibrium work from coordinate-only BINPOS files;
- generating a corrected PMF;
- starting parent glucose-versus-PaG calibration;
- ranking MC-Q1.

## External clarification required

The authors have been asked to clarify:

- whether the galactose A trajectory was intentionally pulled beyond 15 Å;
- why runs 4 and 5 are identical;
- whether a second independent A trajectory was intended but not deposited;
- whether the fructopyranose 20 Å endpoint is intentional;
- which fructofuranose guest group was intended;
- whether source coordinates, target histories and restraint-generation scripts are available.

## Reproducibility records

- Native reader: `research/binpos_native_reader.py`.
- Coordinate audit: `research/audit_gluhut_binpos_coordinates.py`.
- Execution wrapper: `research/run_gluhut_binpos_coordinate_audit.py`.
- Endpoint audit: `research/summarize_gluhut_smd_endpoint_provenance.py`.
- Step audit: `research/audit_gluhut_binpos_step_discontinuity.py`.
- Coordinate output: `research/data/gluhut-binpos-coordinate-audit.json`.
- Endpoint output: `research/data/gluhut-smd-endpoint-provenance.json`.
- Step output: `research/data/gluhut-binpos-step-discontinuity.json`.
- Workflows and tests are stored with each audit.

## Evidence boundary

This is trajectory-coordinate and deposit-provenance evidence. It is not an SMD-work calculation, PMF, binding-affinity result, NNC2215 result, PaG result, MC-Q1 result, efficacy, safety or working-drug result.
