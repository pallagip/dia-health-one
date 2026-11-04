# AWS-1 GluHUT dataset acquisition manifest

**Updated:** 4 November 2026  
**Gate:** AWS-1 — reproduce deposited GluHUT WHAM/PMF analysis  
**Compute authorization:** no new candidate simulation; no Boltz host–guest screen  
**Initial AWS spend cap:** $25

## Objective

Reproduce one published GluHUT glucose potential of mean force from the authors' deposited WHAM inputs, then test whether the same analysis reproduces the qualitative sugar-affinity ranking without retrospective retuning.

Passing this gate would validate an analysis pipeline. It would not validate MC-Q1, establish PaG selectivity, or discover a drug.

## Public source record

- Record title: *Free Energy Simulations on a Biomimetic Glucose Receptor: Understanding the Selectivity of GluHUT*
- DOI: `10.5281/zenodo.19664426`
- Version: `v1`
- Publication date: `20 April 2026`
- Licence: Creative Commons Attribution 4.0 International
- Deposited size: approximately `1.6 GB`
- Simulation package: AMBER
- Trajectory format: Scripps binpos with AMBER topology files
- Individual archive downloads: available

## Deposited file manifest

| File | Displayed size | Deposited MD5 | AWS-1 decision |
|---|---:|---|---|
| `README.txt` | 1.8 kB | `934eddb24acf743a0c8d7113f6111ce8` | Download first |
| `WHAM_PMF.zip` | 41.6 MB | `dbc6723b8adcb85deb8230f63b6a51da` | Download first |
| `SMD_trajectories.zip` | 37.9 MB | `9f3c6733dc608d2bd37c488c560a5443` | Hold until needed for path provenance |
| `US_Trajectories.zip` | 123.8 MB | `c460e6a1ee97ac531fc284d8f07e381e` | Hold until histogram generation must be audited |
| `Figures.zip` | 29.0 MB | `c3d75833cefc48e907bb21d1f320ba14` | Optional comparison source |
| `Unrestrained_MD_trajectories.zip` | 701.8 MB | `e28eb41df0683b7f5c1cf653916fcece` | Excluded from initial acquisition |
| `GIST.zip` | 678.8 MB | `80c55a03d3c17ba7ba34732acbccf276` | Excluded from initial acquisition |

The initial acquisition is therefore approximately **41.6 MB**, rather than the full 1.6 GB deposit.

## Dataset README summary

The deposited README states that:

- unrestrained MD includes empty GluHUT and complexes with glucose, galactose, fructopyranose and fructofuranose;
- SMD trajectories dissociate each sugar using time-dependent harmonic restraints;
- US trajectories restrain each sugar at defined reaction-coordinate values;
- `WHAM_PMF.zip` contains the inputs and outputs used to construct PMFs;
- `GIST.zip` contains the empty-receptor solvation analysis;
- trajectories were produced with AMBER and are supplied in binpos format with corresponding topology files.

## Confirmed `WHAM_PMF.zip` structure

The public archive preview exposes:

- `WHAM_PMF/run-wham.txt` — 129 bytes;
- `WHAM_PMF/example-meta-file-histogram-collection-glucose-US-1.dat` — 3.5 kB;
- `histogram-collection-glucose-US-1/`;
- `histogram-collection-glucose-US-2/`;
- `histogram-collection-glucose-US-3/`;
- `histogram-collection-galactose-US-1/`;
- `histogram-collection-galactose-US-2/`;
- `histogram-collection-galactose-US-3/`;
- `histogram-collection-galactose-US-routeA/`;
- three fructopyranose collections;
- three fructofuranose collections.

An example exposed window file is:

`WHAM_PMF/histogram-collection-glucose-US-1/0-0-test-dist_vs_t.dat`

The public preview did not expose the exact contents of `run-wham.txt` or the metadata file. Those are the first post-download inspection targets.

## First reproducibility package

The smallest confirmed subset for one glucose PMF is:

1. `WHAM_PMF/run-wham.txt`;
2. `WHAM_PMF/example-meta-file-histogram-collection-glucose-US-1.dat`;
3. the complete `WHAM_PMF/histogram-collection-glucose-US-1/` directory.

The first analysis must not edit those source inputs. Any generated metadata or normalized paths must be stored separately with explicit provenance.

## Immutable acquisition requirements

1. Download `README.txt` and `WHAM_PMF.zip` from the Zenodo record.
2. Verify the deposited MD5 values before extraction.
3. Generate SHA-256 values and record exact byte sizes, final URLs, UTC time and HTTP validators.
4. Store the downloaded files read-only.
5. Audit the ZIP for unsafe paths, duplicate entries, symlinks and path traversal.
6. Record every archive path, uncompressed size and CRC before extraction.
7. Extract into a content-addressed working directory.

## Analysis sequence

### A. One-path glucose reproduction

- read `run-wham.txt` verbatim;
- record WHAM executable/version and all command-line parameters;
- rerun glucose `US-1` from the deposited metadata and histograms;
- compare a regenerated PMF with a deposited output or a documented digitisation of the published figure;
- report the numerical tolerance and reference-alignment rule.

### B. Three-path uncertainty

- run glucose `US-1`, `US-2` and `US-3` independently;
- quantify between-path variation;
- test equilibration-discard and histogram-binning sensitivity;
- require a stable bulk plateau and interpretable overlap.

### C. Hold-out ranking

Apply the unchanged protocol to at least one lower-affinity sugar. The preferred complete gate compares glucose, galactose, fructopyranose and fructofuranose.

## Pass criteria

AWS-1 passes only if:

1. deposited hashes and sizes verify;
2. extraction and parsing are deterministic;
3. one glucose PMF is regenerated from deposited inputs;
4. three glucose routes give interpretable uncertainty;
5. a hold-out sugar ranking agrees qualitatively with the publication;
6. reasonable analysis choices do not reverse the ranking;
7. raw inputs, outputs, commands, deviations and negative results are retained.

## Stop criteria

Stop before generating new trajectories if:

- checksums fail;
- required histogram or metadata files are absent;
- the published PMF cannot be regenerated;
- hidden manual choices materially alter the result;
- between-path uncertainty is as large as the claimed sugar difference;
- the ranking requires sugar-specific retuning.

## Current state

- Public source record identified and inventoried.
- Initial download reduced from 1.6 GB to 41.6 MB plus README.
- Required glucose metadata and histogram collection confirmed in the public preview.
- Exact WHAM command and metadata contents still require checksum-verified acquisition.
- No AWS credits have been spent.
- No candidate ranking is authorized until AWS-1, AWS-2, AWS-3 and parent calibration pass.
