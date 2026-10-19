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

The dashboard's **Run reproduction** button independently recalculates the same paired fold changes in browser JavaScript and provides a run-specific JSON manifest with input and script hashes. Archived browser run `data/browser-runs/DHO-NNC2215-FIG3D-20261018020919.json` matches the Python result's source-data SHA-256 and 3.1576× / 12.4720× fold changes. The latest CPU rerun (`DHO-NNC2215-FIG3D-20261018T111522Z`) independently produced the same checksums and paired fold changes; its full source-cell manifest is archived at `data/analysis-runs/DHO-NNC2215-FIG3D-20261018T111522Z.json`. This confirms the published Figure 3d arithmetic in two implementations; it is CPU analysis, not molecular dynamics or a GPU result.

The second script processes the published group-level summaries for the 2026 pig clamp study. It writes `data/pag-clamp-reproduction.json` with the source CSV checksum and descriptive treatment-to-saline ratios. These ratios divide published group means; they are not paired animal-level effects. The study's individual animal measurements and hypothesis tests are not reanalysed by this script.

The third script applies the paper's Appendix S4 Equation 4b and Table S2 `CoeffPaG` estimate to the published human and pig PaG concentration anchors. It writes `data/pag-sensitivity-reproduction.json` with checksums for the parameter and concentration inputs. The result is a reproduction of the model's linear SI component at those concentrations, not the full PK/PD model, receptor-affinity assay, GIR response, or an individual response prediction.

The fourth script recalculates two ratios from the 2026 paper's reported binding summaries: the ratio of glucose and PaG Kd means, and the zero-glucose hIR-A relative-affinity ratio with versus without 70 μmol/L PaG. It writes `data/pag-binding-summary-reproduction.json`. These are arithmetic checks on rounded summary values; the source curves and individual replicate values are not refitted.

`python3 research/audit_pag_assay_source.py` retrieves the official Europe PMC full-text XML, verifies the DOI and license, extracts the two exact zero-glucose affinity means quoted in the article text, and joins the reviewed ACS Supporting Information inventory. It writes `data/pag-assay-source-audit.json` with source checksums, assay context, and the data-access boundary. The article reports Figure 3 means ± SD from three experiments, but the reviewed article and single SI PDF do not provide numeric curve points, exact SD values, or individual replicates in machine-readable form. The source is licensed CC BY-NC-ND 4.0; the public dashboard links to the article rather than reproducing or digitizing its figure. The six-fold calculation uses only the two rounded text values and is not a reanalysis of the full curve.

`NVIDIA_INNOVATION_LAB_BRIEF.md` records the prepared GPU request, the proposed 450 ns pilot, evidence limits, and the reason a target-specific GPU-hour estimate is not available before a valid target-system benchmark. It is not an application confirmation: the authenticated catalog search currently returns no Innovation Lab tile.

The structural audit is reproducible with `python3 research/validate_nnc2215_models.py`. It checks the publisher file hashes, insulin sequences, backbone coverage, native disulfides, and PDB modification-site/link records. A second audit, `python3 research/audit_nnc2215_connectivity.py`, walks the atom-serial `CONECT` and named `LINK` records. The selected PHX/LYX/MC heavy atoms plus PheB1 N and LysB29 NZ form matching 113-node, 120-edge adjacency graphs after explicitly reconciling two LYX labels: CD6↔C for the MC C31-linked roof carbon, and OD↔O for the oxygen whose only heavy-atom neighbor is LYX CD1. This is an adjacency match, not a validated chemical graph. Hydrogen-level records still differ: the closed model contains an unconnected PHX H2 near PheB1 N and records one fewer LYX CD2 hydrogen neighbor. The wwPDB PDB specification does not provide bond orders or formal-charge assignments; stereochemistry, protonation, charges, and force-field parameters remain unvalidated. See `data/nnc2215-model-validation.json`, `data/nnc2215-connectivity-audit.json`, and `MILESTONE_REPORT.md`.

The formula decomposition is reproducible with `python3 research/audit_nnc2215_model_formula.py`. It verifies publisher PDB hashes and coordinate records, derives the desB30 insulin parent formula from the A/B sequences and three disulfide bonds, and calculates the net B1/B29 attachment formulas from the Nature Supporting Information building blocks and reaction steps. The source-derived attachments sum to C72H85N15O22, matching the published NNC2215 formula minus the sequence-derived parent. In both model states, the PHX heavy-element counts match the B1 component (C10O8); the combined LYX/MC B29-linked pool contains C64N15O14 versus source-derived C62N15O14, leaving two excess carbon coordinate records. This localizes the element-count residual to the B29-linked pool but does not map atoms or validate a molecular graph, bond orders, stereochemistry, charges, or force-field parameters. The patent WO2023144240A1 Example 2 reports the same named conjugate and calculated ions as the Nature supplement; it cross-checks product identity, not the PDB topology. The open model's one-oxygen count difference is the missing ASN A21 OXT record present in the closed model. See `data/nnc2215-model-formula-audit.json` for the source formulas, assumptions, per-state comparisons, hashes, and run ID.

`python3 research/audit_nnc2215_linker_mapping.py` compares bounded CONECT paths in each PDB with the numbered B29 alkyne linker in Nature Extended Data Fig. 7a and macrocycle roof propyl-azide tail in Fig. 8i. The graph supports a candidate map for the six linker carbons (source C3–C8) and the three-carbon tail (source C21–C23) in both states. Under that candidate partition, two tail carbons are in LYX and one is in MC, leaving 54 source macrocycle carbons expected in MC; the PDB MC residue contains 56. The two-carbon surplus is therefore localized to the MC residue under this candidate map, but the two unsupported MC atom records are not identified. The atom numbering is inferred from a unique adjacency path; no complete source-to-PDB atom map is supplied. No atoms are removed or chemically retyped. The result remains insufficient for force-field parameterization or production MD. See `data/nnc2215-linker-mapping-audit.json` for per-edge checks, source/model atom candidates, component counts, and checksums.

The connectivity script only parses the original PDB records and does not alter coordinates or infer chemistry. A preliminary RDKit bond-order perception attempt on the 164-atom HET component did not complete within roughly one minute and was interrupted; no inferred bond-order result was retained or used. A production topology still requires a tractable method checked against source chemistry and reported composition.

`data/nnc2215-macrocycle-reference-audit.json` traces the macrocycle to the 2019 glucose-receptor paper cited by the 2024 NNC2215 article and records the accessible synthesis/modeling evidence. The 2019 SI models a nonacarboxylate surrogate because the full receptor caused software errors. This improves provenance but does not resolve the exact atom mapping, full conjugate graph, or force-field parameters for the NNC2215 PDB residues.

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

### NVIDIA portability preflight

`run_pag_gpu_portability.py` loads the checksum-tracked initial PDB and serialized OpenMM system from the Apple M4 PaG-only pilot, then repeats its short NPT setup on an NVIDIA OpenCL device. It records the GPU/driver, OpenMM and Python versions, input hashes, two independent seeds, thermodynamic summaries, and output hashes. The committed inputs have SHA-256 values `1b5bf3cff8f6c433069d621636e89150000e5fc3dd7a3c5d76e15d9d8436fd4b` (initial PDB) and `1c5d1de602db30728bb342cc2607c5b3720c93421fe7bace6a103eb6a632082f` (System XML).

With an NVIDIA runtime that provides `nvidia-smi` and OpenMM 8.6.1, run:

```sh
python research/run_pag_gpu_portability.py \
  --input-pdb https://feat-reproducibility-dashboa.dia-health-one.pages.dev/research/data/pag-solvation-preflight/DHO-PAG-SOLVATION-20261017T191929Z-initial.pdb \
  --system-xml https://feat-reproducibility-dashboa.dia-health-one.pages.dev/research/data/pag-solvation-preflight/DHO-PAG-SOLVATION-20261017T191929Z-system.xml \
  --output-dir /content/pag-t4-preflight
```

This is a planned portability run until its manifest is present under `data/pag-t4-preflight/`. It contains one PaG molecule in solvent and no insulin or receptor. It checks GPU execution and short-term system stability only; it cannot reproduce PaG binding, the reported 500 ns Desmond study, or an affinity value. A free temporary GPU runtime is not an NVIDIA Inception allocation.

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

The dashboard includes an interactive side-by-side view of the publisher's open and closed PDB coordinates. It uses the vendored 3Dmol.js 2.5.5 browser bundle under BSD-3-Clause; `data/structure-viewer-library.json` records its registry integrity, SHA-256, version, and license, and the matching license text is in `assets/vendor/3Dmol-2.5.5-LICENSE.txt`. The viewer renders source coordinates for inspection and provides no chemical or simulation validation. It does not change the target MD readiness gate.

### ChimeraX receptor-session audit

The same Figshare dataset also provides the 79 MB ChimeraX session [“NNC2215 in open and closed forms docked into insulin receptor.cxs”](https://ndownloader.figshare.com/files/48308389). The audit verifies its MD5 against the publisher's file metadata and parses the embedded active/open and inactive/closed NNC2215 snapshots. Their atom-identity sequences match the corresponding public PDB coordinate records after applying the documented `BGLC` → `BGL` alias to the open model. The session contains 990 atoms / 1,010 bond edges in the open snapshot and 968 atoms / 986 bond edges in the closed snapshot.

This is a session-content audit, not chemical validation. At the pinned ChimeraX source revision, session serialization stores bond endpoint atom-index pairs but no authoritative bond-order field. The target PHX/LYX/MC residues have no explicit IDATM type strings; cached computed types are derived labels. The two-carbon excess in the B29-linked LYX/MC pool remains unresolved, and no force-field parameters are validated. Do not use this session as a production-MD topology. Machine-readable findings are in `data/nnc2215-chimerax-session-audit.json`; the source file is distributed under the dataset's CC BY 4.0 licence.

To reproduce the audit, install `research/chimerax-session-audit-requirements.txt`, download the CXS file, and decompress it with LZ4:

```sh
lz4 -d -f "NNC2215 in open and closed forms docked into insulin receptor.cxs" nnc2215-session.raw
python research/audit_nnc2215_chimerax_session.py \
  --session nnc2215-session.raw \
  --cxs "NNC2215 in open and closed forms docked into insulin receptor.cxs"
```

The script checks the downloaded session's published MD5 and writes the complete audit record. Interpretation of session fields is grounded in ChimeraX source at commit [`b18ee9afe9a1536d0044139dcdc4c0ea539d2f7f`](https://github.com/RBVI/ChimeraX/tree/b18ee9afe9a1536d0044139dcdc4c0ea539d2f7f/src/bundles/atomic_lib/atomic_cpp/atomstruct_cpp). This work ran on the local CPU and did not use an NVIDIA GPU.

## PaG evidence boundary

Sturis et al., [“Interaction between Paracetamol Glucuronide and a Glucose-Sensitive Insulin Molecule (NNC2215): Increased Risk of Hypoglycemia”](https://pubs.acs.org/doi/10.1021/acsptsci.6c00069), *ACS Pharmacology & Translational Science* (2026), DOI `10.1021/acsptsci.6c00069`, is represented on the dashboard through explicitly labelled published summary values. The reported Kd means are 1.91 mM for D-glucose and 0.17 mM for PaG; their ratio (about 11.2) is calculated from those means. The paper reports 0.9% and 5.4% hIR-A relative affinity at zero glucose without PaG and with 70 μM PaG. These rounded summaries are not independently reanalysed raw values.

The PaG study's [Supporting Information](https://acs.figshare.com/articles/journal_contribution/Interaction_between_Paracetamol_Glucuronide_and_a_Glucose-Sensitive_Insulin_Molecule_NNC2215_Increased_Risk_of_Hypoglycemia/32292700) specifies Desmond v6.6 (Schrödinger Release 2021-2), OPLS4, SPC explicit water, 0.15 mol/L NaCl, and a 500 ns production phase. A representative snapshot from the glucose-bound macrocycle run was used to build the PaG starting configuration in Maestro, followed by the same protocol. The input is described as an in-house crystal structure of a related macrocycle; its coordinates are not supplied. The reported method is transcribed in `data/pag-md-protocol.json`; the SI PDF is not redistributed here. Exact MD therefore cannot be rerun from the public input set currently available to this project.

The 2024 NNC2215 paper cites Tromans et al. (2019), [“A biomimetic receptor for glucose”](https://www.nature.com/articles/s41557-018-0155-z), as the source of its macrocycle model. That paper's [Supporting Information](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41557-018-0155-z/MediaObjects/41557_2018_155_MOESM46_ESM.pdf) describes Monte Carlo molecular-mechanics models of a receptor surrogate and glucose complex, but the reviewed PDF supplies figures and methods rather than coordinates. It is useful model provenance, not a replacement for the unpublished 2026 seed coordinates. The relationship and boundary are recorded in `data/pag-md-protocol.json`.

PDB 6PXV and 2VJZ are reference structures used in the 2024 paper's modelling. The open/closed model files in the Figshare dataset are distinct model-built NNC2215 coordinate files, not experimental structures.

The pig-clamp panel reports the article's GIR AUC from 300 to 420 minutes at a 3 mmol/L glucose clamp. Means, SEM, sample sizes, and p-values are transcribed from the main article's Table 1; ratios to saline are derived from the reported means. `data/pag-model-parameters.csv` transcribes the minimal-model values in Supporting Information Table S2. The equations and these selected parameter estimates are not enough to recreate the authors' full PK/PD fit: individual-animal time series and some model inputs are not published with the article.

The model-coefficient panel evaluates `SI = SIbase * (1 + ConcPaG * CoeffPaG / 100)` using `CoeffPaG = 0.7784% per μmol/L`. The 70 μmol/L human anchor is the article's estimated mean maximum after 1000 mg oral paracetamol; 102 and 426 μmol/L are the reported mean maxima after low- and high-dose IV PaG in pigs. Values are equation-derived and remain separate from observed GIR and receptor-binding outcomes.

## Public analogue PMF reanalysis: GluHUT and D-glucose

Eades et al., [“Free energy simulations on a biomimetic glucose receptor: understanding the selectivity of GluHUT”](https://pubs.rsc.org/en/content/articlehtml/2026/sc/d6sc00738d), report a glucose binding free energy of −12.5 ± 0.9 kcal/mol. Their public [Zenodo dataset](https://zenodo.org/records/19664426) includes the umbrella-window time series used for a source-data reproduction of the WHAM analysis. The downloaded `WHAM_PMF.zip` SHA-256 is `e0e7b7c04326b8fe86cf71e2193315f34fb865168cf11117b3caf9253119c6eb`.

`reproduce_gluhut_glucose_pmf.py` extracts the three glucose runs, checks the 51-window metadata, and invokes [Grossfield WHAM](https://github.com/agrossfield/wham) 2.0.11 from source commit `43bdd6fcd1916ae2b169cb52ac326b7eea567abb`. It uses 400 bins over 0–15 Å, 300 K, tolerance `1e-6`, and a 13–15 Å bulk reference. Each of the three independent PMFs contains 51 windows and 765,051 frames across windows. The Zenodo example metadata uses a spring parameter of 100; the paper reports an Amber restraint constant of 50, consistent with WHAM's documented convention that Amber omits the 1/2 factor.

The reanalysed run values are −13.794, −12.790, and −11.760 kcal/mol. Their mean ± sample SD is −12.782 ± 1.017 kcal/mol, 0.282 kcal/mol below the published mean and within its reported one-SD range. These are WHAM calculations on the authors' trajectories, not new MD. The JSON result is `data/gluhut-glucose-wham-reproduction.json`, and the generated figure is `assets/gluhut-glucose-pmf.svg`. Reproduce with Python dependencies for the analysis/plot and a Grossfield WHAM binary built from the pinned commit:

```sh
python research/reproduce_gluhut_glucose_pmf.py \
  --archive /path/to/WHAM_PMF.zip \
  --wham-binary /path/to/wham
```

GluHUT is a biomimetic glucose receptor, not the NNC2215 insulin conjugate. The dataset contains glucose, not PaG. This analogue demonstrates that the published PMF post-processing can be reconstructed; it does not reproduce NNC2215–PaG binding, validate the NNC2215 force field, or predict receptor/clinical outcomes. The WHAM post-processing ran on CPU and did not use an NVIDIA GPU.

## PaG GAFF2 ligand-only file audit

`validate_pag_gaff2.py` compares the prepared, singly deprotonated PaG carboxylate SDF with the AmberTools 26 GAFF2/AM1-BCC Mol2 and Amber files. The source structure has formula C14H16NO8−, 39 atoms, and five CIP centres (`S, S, R, S, S`). The audit confirms matching atom element order and 40-bond connectivity, and geometry-based CIP labels match at all five centres. The one localized SDF C=O bond is represented as a single bond in the GAFF2 carboxylate resonance form; this is the expected carboxylate bond representation, not a connectivity mismatch.

The Mol2 partial-charge sum is −1.000001 e. OpenMM 8.6.1 loads the 39-atom Amber topology with total charge −1.00000099 e and finite initial potential energy 106.799 kcal/mol. Machine-readable details and SHA-256 checksums are in `data/pag-gaff2-validation.json`; the parameter and coordinate files are under `ligands/pag-gaff2/`. This is an isolated-ligand preparation/file check. It does not validate the force field against experiment, transfer parameters to the NNC2215 conjugate, or calculate a pose or affinity. A finite initial energy is a load/geometry sanity check, not a stability or accuracy result.

## Interpretation limits

- The exploratory Hill fit is weakly identified: its EC50 bootstrap upper tail reaches the numerical fit ceiling. No 6 mmol/L threshold is fitted or inferred.
- The only PaG molecular dynamics in this demo are two short simulations of PaG alone in water. Separate generic 1UBQ engine preflights include no project ligands. No NNC2215–PaG binding simulation, docking, binding free-energy calculation, or patient-response model is run. The GluHUT PMF panel reanalyses the source authors’ trajectories; the published NNC2215 model PDBs are catalogued inputs, not evidence that target binding was simulated.
- GPU acceleration is not needed for source-data calculations. A generic OpenMM 8.6.1 OpenCL preflight ran on a temporary free Colab Tesla T4: three 100 ps 1UBQ runs measured a median 270.175 ns/day. This is not target-system throughput or an NVIDIA Inception/partner allocation. Partner programs have separate eligibility and application steps; no GPU allocation is confirmed. HunAIFA’s cross-border guidance intake is submitted. See `data/nvidia-engine-preflight.json` for the preflight record.
- The original PaG MD starting coordinates are not public in the supporting information. Any later study using the published NNC2215 models must be labelled as a public-input proxy, not an exact rerun of the PaG simulation.
- The Apple M4 and Colab T4 1UBQ preflights verify that a generic solvated protein system can be constructed and stepped on OpenMM CPU/OpenCL platforms. The T4 timing is not a controlled comparison with Apple M4, does not estimate target-system performance, and does not validate the PaG/NNC2215 setup.
- The new NNC2215 model audit confirms sequence/backbone/disulfide and modification-link records, but flags duplicate atom names and an open/closed linker-attachment naming difference in the nonstandard residues. Chemical graphs and parameters remain unvalidated; see `MILESTONE_REPORT.md` and the audit JSON.
- All source data, derived calculations, model-built structures, and future simulation outputs must remain separately labelled.
