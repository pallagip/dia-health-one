# Dia Health One research reproducibility demo

This directory contains a small, inspectable reanalysis of source values from the NNC2215 glucose-response study, plus a dashboard that presents the published PaG findings separately. It is a computational research demonstration, not a clinical model.

## Run the analysis

From the repository root:

```sh
python3 research/run_reproduction.py
python3 research/reproduce_pag_clamp.py
python3 research/reproduce_pag_sensitivity.py
python3 research/reproduce_pag_binding_summary.py
```

The script uses only the Python standard library. It reads `data/source-data.csv`, calculates the mean and sample standard deviation at each glucose concentration, calculates paired fold changes within each independent experiment, and writes:

- `data/reproduction.json` — machine-readable result and input SHA-256.
- `assets/glucose-response.svg` — original vector chart from the same values.

The JSON result includes the exact source rows used, with experiment IDs and workbook cell references, so the result can be reviewed without reconstructing the input from the chart.

The second script processes the published group-level summaries for the 2026 pig clamp study. It writes `data/pag-clamp-reproduction.json` with the source CSV checksum and descriptive treatment-to-saline ratios. These ratios divide published group means; they are not paired animal-level effects. The study's individual animal measurements and hypothesis tests are not reanalysed by this script.

The third script applies the paper's Appendix S4 Equation 4b and Table S2 `CoeffPaG` estimate to the published human and pig PaG concentration anchors. It writes `data/pag-sensitivity-reproduction.json` with checksums for the parameter and concentration inputs. The result is a reproduction of the model's linear SI component at those concentrations, not the full PK/PD model, receptor-affinity assay, GIR response, or an individual response prediction.

The fourth script recalculates two ratios from the 2026 paper's reported binding summaries: the ratio of glucose and PaG Kd means, and the zero-glucose hIR-A relative-affinity ratio with versus without 70 μmol/L PaG. It writes `data/pag-binding-summary-reproduction.json`. These are arithmetic checks on rounded summary values; the source curves and individual replicate values are not refitted.

The structural audit is reproducible with `python3 research/validate_nnc2215_models.py`. It checks the publisher file hashes, insulin sequences, backbone coverage, native disulfides, and PDB modification-site/link records. It deliberately does not infer bond orders, stereochemistry, protonation, partial charges, or force-field parameters. The present audit flags duplicate names in PHX/LYX/MC residues and different open/closed LYX-to-macrocycle link atom names; these must be reconciled before system preparation. See `data/nnc2215-model-validation.json` and `MILESTONE_REPORT.md`.

`data/target-md-pilot-protocol.json` records a proposed three-condition pilot (apo, glucose reference, PaG) with three 50 ns replicas per condition. It is explicitly not run: the full conjugate topology, PaG starting pose, force-field compatibility, and GPU allocation remain unresolved. It defines a public-input proxy, not an exact rerun of the 2026 macrocycle-analogue simulation.

The dashboard also re-runs the paired calculations in the browser from the CSV and offers a timestamped JSON run manifest. Serve the repository over HTTP to use the browser run button; browsers do not allow `fetch()` of a sibling CSV from a `file://` page.

## Local engine preflight

The engine preflight runs a generic, solvated ubiquitin benchmark through GROMACS CPU and OpenMM CPU/OpenCL. It uses the public [RCSB 1UBQ structure](https://www.rcsb.org/structure/1UBQ), Amber99SB-ILDN, TIP3P, a 1.0 nm solvent padding, a 2 fs timestep, and 5,000 timed steps. The input PDB is checked into `benchmarks/1ubq.pdb`; its SHA-256 is included in `data/engine-preflight.json`.

Recreate the isolated environment and run the workflow from the repository root:

```sh
mamba env create -f research/engine-benchmark-environment.yml
mamba run -n dho-engine-preflight python research/benchmark_engine_preflight.py
```

The run rebuilds and minimizes the test system, records one short throughput measurement per available engine/platform, and regenerates `data/engine-preflight.json` and `assets/engine-preflight.svg`. The latest recorded Apple M4 run measured 37.420 ns/day for GROMACS CPU, 11.141 ns/day for OpenMM CPU, and 126.274 ns/day for OpenMM OpenCL. The pinned environment uses Python 3.13.16, OpenMM 8.6.1, and GROMACS 2026.3. OpenMM OpenCL was available; the GROMACS build had CPU support only. These timings are single-run environment checks, not a controlled comparison: GROMACS and OpenMM use different thermostat implementations and timing boundaries. This is not an NNC2215/PaG run, does not use NVIDIA hardware, and does not validate the conjugate force field or predict binding energetics.

## PaG small-molecule solvation pilot

The PaG identity reference is PubChem CID 83944. Its downloaded neutral 2D SDF is retained unchanged; a separate simulation conformer was built with the glucuronic-acid carboxylate deprotonated (net charge −1), with the five source stereocentres retained. RDKit ETKDGv3/MMFF94s generated and minimized the conformer, AmberTools AM1-BCC assigned charges, and OpenFF Sage 2.3.0 parameterized the 39-atom ligand. These checks cover the isolated PaG molecule only.

The completed pilot ran two independent OpenMM 8.6.1 OpenCL simulations on Apple M4. Each system contained one PaG carboxylate in Amber14/TIP3P water with 0.15 M NaCl; NPT was run at 298.15 K and 1 atm with a 2 fs timestep, 10 ps equilibration, and 100 ps production. The 1,272-atom periodic system implies approximately 0.118 M PaG from its one-ligand box composition. This is not the 70 μM experimental assay concentration. Production temperature means were 298.85 K and 297.54 K; the traces remained finite. Two short replicas establish only that this small-molecule setup ran, not equilibrium sampling, force-field transferability to the insulin conjugate, or binding affinity.

Recreate the chemistry and pilot with the arm64 environment lock:

```sh
micromamba create -p /tmp/dho-pag-chem -f research/pag-chemistry-environment-macos-arm64.lock.txt
micromamba run -p /tmp/dho-pag-chem python research/validate_pag_chemistry.py
micromamba run -p /tmp/dho-pag-chem python research/run_pag_solvation_preflight.py
micromamba run -p /tmp/dho-pag-chem python research/plot_pag_solvation_preflight.py
```

The scripts write a chemistry preflight, timestamped run manifest, raw per-replicate logs and final coordinates, production-only analysis JSON, and SVG figure. The latest run pointer is `data/pag-solvation-preflight/latest.json`. Raw logs and hashes are retained; earlier pilot attempts are marked superseded in their manifests.

## Source and extraction

- Høeg-Jensen et al., “Glucose-sensitive insulin with attenuation of hypoglycaemia,” *Nature* (2024), [article and methods](https://www.nature.com/articles/s41586-024-08042-3), DOI `10.1038/s41586-024-08042-3`.
- The publisher's [Source Data workbook](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-024-08042-3/MediaObjects/41586_2024_8042_MOESM5_ESM.xlsx), worksheet `Figure 3d`.
- Workbook downloaded 2026-10-17 UTC; SHA-256 `09a9b6aa4103f9d61af7df03886d59f16d6b09f69f2fe768260548061c7dc81d`.
- Source assay label: hIR-A affinity relative to human insulin, with 1.5% HSA. Three independent experiments, five D-glucose concentrations (0, 3, 5, 10, 20 mmol/L), and two series (NNC2215 and insulin degludec).
- `source-data.csv` contains the workbook's numeric values, not values digitised from a plot. Each row records the source worksheet and source cell. The workbook itself is linked above and is not copied into this repository.

The analysis reports the mean and sample SD across the three experiment-level values at each concentration. To compare with the paper's reported fold changes, it first divides the high-glucose value by the corresponding low-glucose value within each experiment, then reports the mean and sample SD across those three paired ratios. This gives 3.1576× for 3→20 mmol/L and 12.4720× for 0→20 mmol/L, which round to the reported 3.2× and 12.5×. Dividing the group means instead gives different values (3.1379× and 12.2725×); both calculations are preserved in the output to make the analysis choice visible.

## Exploratory Hill-fit identifiability check

`fit_nnc2215_hill_response.py` fits the 15 NNC2215 source-workbook values with one shared apparent Hill slope and EC50, measured 0 mmol/L baselines fixed per experiment, and an experiment-specific upper response. It uses unweighted nonlinear least squares and 2,000 deterministic bootstrap resamples of the three independent experiments. Run it in the pinned arm64 chemistry environment:

```sh
micromamba run -p /tmp/dho-pag-chem python research/fit_nnc2215_hill_response.py
```

The recorded point estimates are EC50 = 13.43 mmol/L and apparent Hill slope = 1.115. The 2.5th bootstrap percentile for EC50 is 8.06 mmol/L, while 64 of 2,000 resamples hit the 1,000 mmol/L optimization ceiling; the 97.5th percentile is therefore censored by that fit bound. The apparent-slope bootstrap percentile range is 0.675–1.491. These percentiles are exploratory and are not formal confidence intervals: there are only three independent experiment clusters, and the observed series does not show a high-glucose plateau. The result demonstrates that this dataset does not identify a reliable EC50. It does not establish a 6 mmol/L threshold, binding-site count, cooperativity, receptor affinity in a new condition, or clinical response. The script writes the machine-readable fit and provenance record to `data/hill-fit-summary.json` and the figure to `assets/hill-fit.svg`.

To verify that the checked-in CSV still matches the publisher workbook, download the workbook and run:

```sh
python3 research/verify_source_extract.py /path/to/41586_2024_8042_MOESM5_ESM.xlsx
```

The verifier checks the workbook checksum, worksheet name, source-cell references, and all 30 numeric values.

## Published NNC2215 model coordinates

The 2024 study's [Springer Nature Figshare dataset](https://springernature.figshare.com/articles/dataset/3D_models_of_NNC2215/26526460) publishes model-built PDB coordinates for open and closed NNC2215 states. The two PDB files are included under `structures/` under the dataset's CC BY 4.0 licence. Their original filenames, download URLs, source checksum values, local SHA-256 values, and attribution are recorded in `data/structure-manifest.json`.

The open model places free D-glucose in the B29 macrocycle; the closed model places the tethered glucoside there. These are useful author-supplied starting models, not experimentally resolved NNC2215 structures. They do not supply the separate in-house glucose-bound macrocycle crystal coordinates used to initialize the 2026 PaG molecular-dynamics work, and their presence does not establish force-field readiness.

Run `python3 research/inspect_structure_models.py` to verify the local PDB checksums and produce `data/structure-preflight.json`. This inventory records coordinate counts, chain/residue labels, disulfide and covalent-link records, and nonstandard residue names. It is a file/record preflight, not chemical-graph or force-field validation.

## PaG evidence boundary

Sturis et al., [“Interaction between Paracetamol Glucuronide and a Glucose-Sensitive Insulin Molecule (NNC2215): Increased Risk of Hypoglycemia”](https://pubs.acs.org/doi/10.1021/acsptsci.6c00069), *ACS Pharmacology & Translational Science* (2026), DOI `10.1021/acsptsci.6c00069`, is represented on the dashboard through explicitly labelled published summary values. The reported Kd means are 1.91 mM for D-glucose and 0.17 mM for PaG; their ratio (about 11.2) is calculated from those means. The paper reports 0.9% and 5.4% hIR-A relative affinity at zero glucose without PaG and with 70 μM PaG. These rounded summaries are not independently reanalysed raw values.

The PaG study's [Supporting Information](https://acs.figshare.com/articles/journal_contribution/Interaction_between_Paracetamol_Glucuronide_and_a_Glucose-Sensitive_Insulin_Molecule_NNC2215_Increased_Risk_of_Hypoglycemia/32292700) specifies Desmond v6.6 (Schrödinger Release 2021-2), OPLS4, SPC explicit water, 0.15 mol/L NaCl, and a 500 ns production phase. A representative snapshot from the glucose-bound macrocycle run was used to build the PaG starting configuration in Maestro, followed by the same protocol. The input is described as an in-house crystal structure of a related macrocycle; its coordinates are not supplied. The reported method is transcribed in `data/pag-md-protocol.json`; the SI PDF is not redistributed here. Exact MD therefore cannot be rerun from the public input set currently available to this project.

The 2024 NNC2215 paper cites Tromans et al. (2019), [“A biomimetic receptor for glucose”](https://www.nature.com/articles/s41557-018-0155-z), as the source of its macrocycle model. That paper's [Supporting Information](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41557-018-0155-z/MediaObjects/41557_2018_155_MOESM46_ESM.pdf) describes Monte Carlo molecular-mechanics models of a receptor surrogate and glucose complex, but the reviewed PDF supplies figures and methods rather than coordinates. It is useful model provenance, not a replacement for the unpublished 2026 seed coordinates. The relationship and boundary are recorded in `data/pag-md-protocol.json`.

PDB 6PXV and 2VJZ are reference structures used in the 2024 paper's modelling. The open/closed model files in the Figshare dataset are distinct model-built NNC2215 coordinate files, not experimental structures.

The pig-clamp panel reports the article's GIR AUC from 300 to 420 minutes at a 3 mmol/L glucose clamp. Means, SEM, sample sizes, and p-values are transcribed from the main article's Table 1; ratios to saline are derived from the reported means. `data/pag-model-parameters.csv` transcribes the minimal-model values in Supporting Information Table S2. The equations and these selected parameter estimates are not enough to recreate the authors' full PK/PD fit: individual-animal time series and some model inputs are not published with the article.

The model-coefficient panel evaluates `SI = SIbase * (1 + ConcPaG * CoeffPaG / 100)` using `CoeffPaG = 0.7784% per μmol/L`. The 70 μmol/L human anchor is the article's estimated mean maximum after 1000 mg oral paracetamol; 102 and 426 μmol/L are the reported mean maxima after low- and high-dose IV PaG in pigs. Values are equation-derived and remain separate from observed GIR and receptor-binding outcomes.

## Interpretation limits

- The exploratory Hill fit is weakly identified: its EC50 bootstrap upper tail reaches the numerical fit ceiling. No 6 mmol/L threshold is fitted or inferred.
- The only molecular dynamics in this demo are two short simulations of PaG alone in water. No NNC2215–PaG binding simulation, docking, binding free-energy calculation, or patient-response model is run. The published NNC2215 model PDBs are catalogued inputs, not evidence that target binding was simulated.
- GPU acceleration is not needed for source-data calculations. A generic OpenMM 8.6.1 OpenCL preflight has now run on a temporary free Colab Tesla T4: three 100 ps 1UBQ runs measured a median 270.175 ns/day. This is not target-system throughput or an NVIDIA Inception/partner allocation. Inception is approved; its `$1,000` Lambda, `$10,000` AWS, and `$5,000` Azure catalog requests remain pending, while the separate partner applications and Innovation Lab tile remain unresolved. See `data/nvidia-engine-preflight.json`, `data/nvidia-access-status.json`, and `COMPUTE_STATUS.md`.
- The original PaG MD starting coordinates are not public in the supporting information. Any later study using the published NNC2215 models must be labelled as a public-input proxy, not an exact rerun of the PaG simulation.
- The Apple M4 and Colab T4 1UBQ preflights verify that a generic solvated protein system can be constructed and stepped on OpenMM CPU/OpenCL platforms. The T4 timing is not a controlled comparison with Apple M4, does not estimate target-system performance, and does not validate the PaG/NNC2215 setup.
- The new NNC2215 model audit confirms sequence/backbone/disulfide and modification-link records, but flags duplicate atom names and an open/closed linker-attachment naming difference in the nonstandard residues. Chemical graphs and parameters remain unvalidated; see `MILESTONE_REPORT.md` and the audit JSON.
- All source data, derived calculations, model-built structures, and future simulation outputs must remain separately labelled.
