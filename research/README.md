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

The structural audit is reproducible with `python3 research/validate_nnc2215_models.py`. It checks the publisher file hashes, insulin sequences, backbone coverage, native disulfides, and PDB modification-site/link records. A second audit, `python3 research/audit_nnc2215_connectivity.py`, walks the atom-serial `CONECT` and named `LINK` records. The selected PHX/LYX/MC heavy atoms plus PheB1 N and LysB29 NZ form matching 113-node, 120-edge adjacency graphs after explicitly reconciling two LYX labels: CD6↔C for the MC C31-linked roof carbon, and OD↔O for the oxygen whose only heavy-atom neighbor is LYX CD1. This is an adjacency match, not a validated chemical graph. Hydrogen-level records still differ: the closed model contains an unconnected PHX H2 near PheB1 N and records one fewer LYX CD2 hydrogen neighbor. The wwPDB PDB specification does not provide bond orders or formal-charge assignments; stereochemistry, protonation, charges, and force-field parameters remain unvalidated. See `data/nnc2215-model-validation.json` and `data/nnc2215-connectivity-audit.json`.

The source-mapped B1 glucoside graph and coordinate audit are in `structures/nnc2215-b1-glucoside-attachment.smi`, `data/nnc2215-b1-glucoside-atom-map.json`, and `data/nnc2215-b1-topology-audit.json`. The audit matches all 19 selected atoms and 19 adjacency edges exactly in each model, with one mapping per state. It independently checks the glucose stereochemistry using RDKit's accurate CIP labels and a detached beta-D-glucose reference structure; the derived B1–B22 fragment formula is C120H187N31O39S2 and [M+2H]2+ is 1376.15947, rounding to the Nature SI value 1376.2. The models provide adjacency only: this does not validate bond orders, protonation, charge, force-field parameters, or simulation readiness. Re-run with `micromamba run -p /Users/andrii/scratch/dia-target-tools/env python research/audit_nnc2215_b1_topology.py`.

The formula decomposition is reproducible with `python3 research/audit_nnc2215_model_formula.py`. It verifies original publisher PDB hashes and coordinate records, derives the desB30 insulin parent formula from the A/B sequences and three disulfide bonds, and calculates the net B1/B29 attachment formulas from the Nature Supporting Information building blocks and reaction steps. The source-derived attachments sum to C72H85N15O22, matching the published NNC2215 formula minus the sequence-derived parent. In both original publisher models, the PHX heavy-element counts match the B1 component (C10O8); the combined LYX/MC B29-linked pool contains C64N15O14 versus source-derived C62N15O14, leaving two excess carbon coordinate records. A separate candidate correction removes MC C27/C29 and reconnects retained methyl carbons C28/C30 directly to aryl carbons C73/C75; this source-supported inference matches the product drawing and formula, but is not an author-supplied atom map. The open original model also lacks ASN A21 OXT, which is reconstructed in the corrected candidate. See `data/nnc2215-model-formula-audit.json` for the original source formulas, assumptions, per-state comparisons, hashes, and run ID.

`python3 research/audit_nnc2215_linker_mapping.py` compares bounded CONECT paths in each original PDB with the numbered B29 alkyne linker in Nature Extended Data Fig. 7a and macrocycle roof propyl-azide tail in Fig. 8i. The graph supports a candidate map for the six linker carbons (source C3–C8) and the three-carbon tail (source C21–C23) in both states. Under that candidate partition, two tail carbons are in LYX and one is in MC, leaving 54 source macrocycle carbons expected in MC; the original PDB MC residue contains 56. A separate pairwise-contraction enumeration over 18 degree-two carbon vertices found MC C27/C29 to be the only contracted pair yielding the source-derived graph in both states, within that defined edit class. The corrected coordinate candidates remove C27/C29 and reconnect C28/C30 directly to C73/C75. This is source-supported but still an inferred coordinate repair, not an author-provided atom map. See `data/nnc2215-linker-mapping-audit.json` and `data/nnc2215-carbon-reconciliation-candidate.json` for the original records and correction evidence. Neither establishes force-field parameterization or production MD.

The connectivity script only parses the original PDB records and does not alter coordinates or infer chemistry. A preliminary RDKit bond-order perception attempt on the 164-atom HET component did not complete within roughly one minute and was interrupted; no inferred bond-order result was retained or used. A production topology still requires a tractable method checked against source chemistry and reported composition.

`data/nnc2215-macrocycle-reference-audit.json` traces the macrocycle to the 2019 glucose-receptor paper cited by the 2024 NNC2215 article and records the accessible synthesis/modeling evidence. The 2019 SI models a nonacarboxylate surrogate because the full receptor caused software errors. This improves provenance but does not resolve the unique B29 coordinate mapping or force-field parameters for the NNC2215 PDB residues; the neutral constitutional graph assembled from the available sources is documented below.

`data/target-md-pilot-protocol.json` records a proposed three-condition pilot (apo, glucose reference, PaG) with three 50 ns replicas per condition. It is explicitly not run: the full conjugate topology, PaG starting pose, and force-field compatibility remain unresolved. It defines a public-input proxy, not an exact rerun of the 2026 macrocycle-analogue simulation.

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

This is a planned portability run until its manifest is present under `data/pag-t4-preflight/`. It contains one PaG molecule in solvent and no insulin or receptor. It checks GPU execution and short-term system stability only; it cannot reproduce PaG binding, the reported 500 ns Desmond study, or an affinity value.

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

### Source-derived roof macrocycle graph

The numbered roof propyl-azide drawing in Nature Extended Data Fig. 8i has been transcribed as [`structures/nnc2215-roof-propyl-azide.smi`](structures/nnc2215-roof-propyl-azide.smi). RDKit 2026.03.1 parses and sanitizes the graph as C56H63N15O12 with 83 heavy atoms, 83 unique atom-map IDs, and formal charge 0 in the neutral-acid drawing form. [`data/nnc2215-roof-macrocycle-atom-map.json`](data/nnc2215-roof-macrocycle-atom-map.json) explains the map IDs and source numbering; [`data/nnc2215-roof-macrocycle-graph-audit.json`](data/nnc2215-roof-macrocycle-graph-audit.json) records the source, checks, and hashes. Recreate the audit with `micromamba run -p /Users/andrii/scratch/dia-target-tools/env python research/validate_nnc2215_roof_macrocycle_graph.py`.

This is an isolated component graph, not a completed insulin-conjugate topology. The azide is represented in one conventional resonance form, the source drawing does not assign stereocentres, and graph map IDs do not establish PDB atom correspondence. Simulation-pH protonation, charges, force-field parameters, and target-MD readiness remain unresolved.

An element-labeled adjacency search maps this graph into both corrected candidate PDBs under the six candidate C21–C23/azide-tail anchors. It finds 16 embeddings per state, with 54 of 83 graph map IDs remaining ambiguous across the embeddings. Each selected LYX B30 + MC C1 fragment has 91 heavy atoms (C62N15O14); eight unmatched LYX atoms are assigned to the B29 linker only under the existing candidate partition. See [`data/nnc2215-roof-graph-coordinate-mapping-audit.json`](data/nnc2215-roof-graph-coordinate-mapping-audit.json) and recreate it with `micromamba run -p /Users/andrii/scratch/dia-target-tools/env python research/audit_nnc2215_roof_graph_coordinate_mapping.py`. PDB connectivity records do not provide authoritative bond orders. The result does not establish a unique source-to-PDB atom map, geometry, protonation, force-field parameters, or target-MD readiness. The audit verified the candidate PDB hashes and did not edit either PDB.

### B29 linker–roof product graph and exact candidate audit

The source-derived mapped B29 linker–roof product component is available as [`structures/nnc2215-b29-linker-roof-product.smi`](structures/nnc2215-b29-linker-roof-product.smi), with source roles and map provenance in [`data/nnc2215-b29-linker-roof-product-map.json`](data/nnc2215-b29-linker-roof-product-map.json). The selected component includes the LysB29 NZ endpoint, so its C62N16O14 heavy-element inventory (92 heavy atoms) includes one N beyond the C62N15O14 linker–roof appendage count, which excludes NZ. It has 99 edges, one connected component, and cycle rank 8; Lys CE and the rest of insulin are outside the selected graph.

The B29 product-component topology/coordinate audit is recorded in [`data/nnc2215-b29-topology-audit.json`](data/nnc2215-b29-topology-audit.json). Rerun from the repository root with [`audit_nnc2215_b29_topology.py`](audit_nnc2215_b29-topology.py):

```sh
micromamba run -p /Users/andrii/scratch/dia-target-tools/env python research/audit_nnc2215_b29_topology.py
```

The source-derived B29 product component exactly matches the element-labeled adjacency of each corrected, model-built open and closed candidate PDB over LYX B30, MC C1, and LysB29 NZ. Each state has 16 complete mappings; map IDs 78–92 are fixed across them. The neutral graph has eight bond-label transfer classes per state, differing only by independent carbonyl-O/OH swaps at the three macrocycle carboxyls. Under the separately declared pH 7.4 tri-carboxylate scenario, resonance-normalizing the two C–O edges in each carboxylate reduces those eight classes to one chemistry-equivalent class per state. This does not choose PDB oxygen names or reduce the 16 graph mappings. The source drawing in [Nature Fig. 2a](https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fs41586-024-08042-3/MediaObjects/41586_2024_8042_Fig2_HTML.png) supports the 1,4-disubstituted triazole regiochemistry; the remaining mappings do not represent alternative click regioisomers. The corrected PDBs are model-built, not experimentally resolved, and PDB CONECT/LINK records encode adjacency only. A conditional group-level pH 7.4 ledger is recorded below; a reviewed simulation-ready protonation model, atom-localized charge model, and force field remain open.

The pH 7.4 carboxylate result is a declared graph scenario, not a pKa calculation or measured microstate. The isolated B29 component changes from C62H72N16O14 to C62H69N16O14 with three removed protons and local charge change −3. No coordinates or PDB oxygen labels are changed; the rest of insulin is outside this audit. Details and machine-readable results are in `data/nnc2215-b29-topology-audit.json`.

### B1/B29 attachment-boundary audit — 26 October 2026

[`data/nnc2215-attachment-boundary-audit.json`](data/nnc2215-attachment-boundary-audit.json) reconciles the source component boundary edges with the local insulin contexts in both corrected, model-built candidate PDBs. It checks B1 map 1–2 at PheB1 N–PHX C2, B29 map 84–92 at LysB29 NZ–LYX B30 CD1, and roof maps 78–79 at MC C31–LYX B30 CD6 (open) / C (closed), plus native PheB1 N–CA and LysB29 CE–NZ adjacency. Rerun from the repository root with `micromamba run -p /Users/andrii/scratch/dia-target-tools/env python research/audit_nnc2215_attachment_boundaries.py`.

The coordinate records show the PheB1 N hydrogen discrepancy: the open candidate has a CONECT neighbor H2; the closed candidate has no hydrogen neighbor at N, while PHX H2A is 1.006 Å away without a CONECT edge. LysB29 NZ has the recorded HZ2 neighbor in both states. These are connectivity and proximity observations, not validated hydrogen ownership. This audit closes only the local heavy-atom attachment boundaries. The separate neutral constitutional graph and conditional group-level pH 7.4 charge ledger are reported below; the 16 B29 mappings per state, oxygen-label ambiguity, roof stereochemistry, simulation-ready protonation/charge assignment, force-field parameters, and target MD remain open. It does not edit coordinates or assign the disputed hydrogen.

The neutral-graph geometry screen covers 1,584 mapped edges per state. Distances span 1.2116–1.5611 Å (open) and 1.2116–1.5548 Å (closed); 14 of 16 mappings in each state have one to three heuristic short carboxyl C–O outliers, while two mappings per state fall fully within the broad screen. These distances cannot establish bond order or choose a mapping. Model hydrogen records differ at roof C22/C23 and linker C4–C6; the closed candidate has one fewer recorded hydrogen at linker C4. Those PDB hydrogen records are not validated chemistry.

This work updates source-derived component graphs, a neutral full-conjugate graph, and candidate-coordinate audits. It does not establish a production simulation model, target-binding MD, or a clinical conclusion. The audits ran on CPU.

### Source-derived full insulin–NNC2215 neutral graph

[`data/nnc2215-full-conjugate-graph.json`](data/nnc2215-full-conjugate-graph.json) records the assembled neutral constitutional graph of the full desB30 insulin–NNC2215 conjugate. It joins RDKit 2026.03.1 `MolFromFASTA` neutral residue templates for chains A and B with the mapped B1 glucoside and B29 linker–roof fragments and the three disulfides A6–A11, A7–B7, and A20–B19. The source caps are replaced at PheB1 N and LysB29 NZ as specified by the map legends. The ledger records provenance for every atom and bond: residue-template or source-map origins, explicit attachment-join operations, or declared insulin-disulfide connectivity. Hydrogens are implicit in the graph and have a per-atom ledger.

The graph has formula C325H461N79O97S6, 507 heavy atoms, 527 heavy-atom edges, one connected component, cycle rank 21, and passes RDKit sanitization. This is internal consistency with the formula reported in the Nature SI, not independent experimental confirmation. The corrected open and closed model-built PDBs match the two chain sequences, three disulfide records, and three named attachment LINKs. The candidate checks retain one exact B1 mapping and sixteen exact B29 mappings per state; no coordinates are used to assemble the graph or choose an atom mapping.

The graph uses pinned neutral residue templates and neutral terminal conventions (A-chain N terminus NH2 and C-terminal carboxylic acid; B-chain N terminus modified at B1 and C terminus LysB29 carboxylic acid). The separate conditional ledger below records a conventional group-level pH 7.4 candidate; it does not establish a measured full-system microstate or atom-localized charges. The conditional B29 tri-carboxylate result remains a separate local graph representation. The PheB1 H2/H2A ownership discrepancy, roof stereochemistry, neutral B29 oxygen-label ambiguity, unique B29 coordinate map, force-field parameters, and target MD remain open. The open corrected candidate includes ASN A21 OXT although the original open-model formula inventory reports it absent; the report records both observations. Candidate coordinates were not changed by this audit.

Re-run from the repository root with:

```sh
micromamba run -p /Users/andrii/scratch/dia-target-tools/env python research/audit_nnc2215_full_conjugate_graph.py
```

### OpenFF Sage label-coverage audit

[`data/nnc2215-forcefield-label-coverage.json`](data/nnc2215-forcefield-label-coverage.json) records label coverage from OpenFF Sage 2.3.0 on the neutral 507-heavy-atom full constitutional graph and on all eight protonation combinations of the three carboxyl groups in the isolated roof macrocycle fragment. The full graph received labels for 461 constraints, 988 bonds, 1,753 angles, 2,569 proper torsions, 225 improper torsions, and 968 vdW terms; the label operation returned no electrostatics terms. Constraint, bond, angle, proper-torsion, and vdW counts match independently derived topology expectations for the full graph and each roof state. Improper-torsion counts are recorded, but their completeness is not independently derived. This is label matching, not partial-charge assignment or validation of a hybrid insulin–conjugate force field. Roof-fragment states are not pKa-weighted, and the full graph remains stereochemically incomplete and coordinate-free.

The audit records input, script, and lockfile SHA-256 values and checks exact package versions from the macOS ARM64 lock. To recreate its environment and result from the repository root:

```sh
micromamba create -p /tmp/dho-pag-chem --file research/pag-chemistry-environment-macos-arm64.lock.txt
micromamba run -p /tmp/dho-pag-chem python research/audit_nnc2215_forcefield_label_coverage.py
```

The result is limited to matching Sage labels against topology-derived expectations for constraints, bonds, angles, proper torsions, and vdW interactions; it reports but does not independently validate improper-torsion completeness. It does not validate atom mapping, stereochemistry, protonation populations, atom-localized charges, protein-force-field compatibility, geometry, an OpenMM system, target MD, binding affinity, or clinical activity.

### Hybrid Amber template feasibility inventory — 26 October 2026

[`data/nnc2215-hybrid-feasibility.json`](data/nnc2215-hybrid-feasibility.json) inventories the full source-derived graph against the residue-template names in the pinned OpenMM Amber14 XML include tree. It classifies 50 insulin residues: 16 of 17 observed residue-name labels have literal Amber template names. The graph's NE2-bound histidine hydrogens support HIE as a candidate by inference; the six disulfide cysteines are CYX candidates. The source-mapped appendages contribute 109 heavy atoms, and the audit records the two cross-boundary joins: PheB1 N to B1 map 2 and LysB29 NZ to B29 map 84.

This is a template-name and provenance inventory, not exact residue-template matching. Plain PHE/LYS templates do not parameterize the attached glucoside, linker, triazole, or roof, and the PheB1/LysB29 junction terms still require custom treatment. Terminal and physiological charge compatibility, atom-level charges, custom parameters, and an OpenMM system remain unresolved. Status: `AMBER_TEMPLATE_CANDIDATES_IDENTIFIED; CUSTOM_CONJUGATE_PARAMETERS_REQUIRED`; the validated conjugate-force-field count remains zero.

Whole-graph NAGL charge assignment would extrapolate from its published test-size coverage: the largest reported heavy-atom test bin is 80–89 atoms, with two molecules in that bin and none above it, while this graph has 507 heavy atoms. The highest reported total-atom test bin is 160–169; the source table reports bins rather than an exact maximum. The official AmberTools modified-polymer-residue workflow is a relevant precedent, but its generated charges and types require review; the documented OpenFF Sage/ff14SB Interchange example combines separate components experimentally and does not validate this covalent conjugate. See the audit's source links for the [NAGL model](https://docs.openforcefield.org/projects/nagl-models/en/stable/models/openff-gnn-am1bcc-1.0.0/), [AmberTools tutorial](https://ambermd.org/tutorials/basic/tutorial5/index.php), and [OpenFF Interchange example](https://docs.openforcefield.org/en/latest/examples/openforcefield/openff-interchange/protein_ligand/protein_ligand.html).

Recreate the audit from the repository root in an environment built from the recorded lock:

```sh
micromamba create -p /tmp/dho-hybrid-feasibility --file research/pag-chemistry-environment-macos-arm64.lock.txt
micromamba run -p /tmp/dho-hybrid-feasibility python research/audit_nnc2215_hybrid_feasibility.py
```

The run records input, lock, and every included Amber XML checksum. It inventories template-name candidates only; it does not validate exact residue matching, assign charges, create a complete parameter set or OpenMM system, or run dynamics.

### Isolated B29 GAFF2 parameterization pilot — 26 October 2026

[`data/nnc2215-b29-gaff2-pilot.json`](data/nnc2215-b29-gaff2-pilot.json) records a bounded AmberTools 26 GAFF2/AM1-BCC feasibility pilot on the source-mapped B29 linker–roof component (92 heavy atoms) in the declared tri-carboxylate state, formal charge −3. It uses heavy-atom coordinates from the publisher-supplied active/open model-built candidate, with one of the 16 exact element-adjacency graph embeddings, and explicit hydrogens added by RDKit. Antechamber's SQM geometry optimization was disabled with `maxcyc=0`; the charges therefore depend on the selected model geometry.

All eight carboxylate Lewis labelings produced identical MOL2 and frcmod outputs; all 161 atoms received GAFF2 types. The raw MOL2 partial-charge sum is −2.996005 e, a +0.003995 e residual from the formal −3; no correction was applied. `parmchk2` emitted 26 additional parameter records. None carries an `ATTN` marker, but two have penalty scores of 41.9 and need expert review. The raw input SDFs, MOL2 files, frcmod files, and run logs are in [`data/nnc2215-b29-gaff2-pilot/`](data/nnc2215-b29-gaff2-pilot/).

This supports toolchain feasibility for the isolated capped B29 fragment only. The source NH2 cap is not the attached LysB29 NZ environment; charges at that boundary are not final. The exact map to model coordinates, B1 junction, insulin terminals, custom cross-boundary terms, full topology, and complete charge reconciliation remain open. No validated conjugate force field, OpenMM system, MD, affinity calculation, or clinical inference results from this pilot. The representative [input SDF](data/nnc2215-b29-gaff2-pilot/resonance_01.sdf), [MOL2](data/nnc2215-b29-gaff2-pilot/resonance_01.mol2), and [frcmod](data/nnc2215-b29-gaff2-pilot/resonance_01.frcmod) are available; the JSON indexes all eight variants. The official [AmberTools modified-residue tutorial](https://ambermd.org/tutorials/basic/tutorial5/index.php) describes the Antechamber/prepgen/parmchk2 workflow and cautions that generated charges and terms require review; see also the [Amber24 reference manual](https://ambermd.org/doc12/Amber24.pdf). Re-run from a lock-derived environment with `micromamba run -p <environment-created-from-lock> python research/audit_nnc2215_ambertools_parameterization_pilot.py`.

### Conditional pH 7.4 charge-state ledger

[`data/nnc2215-physiological-charge-states.json`](data/nnc2215-physiological-charge-states.json) records a group-level formal-charge envelope from the neutral graph. Its conventional starting candidate treats GluA4/GluA17/GluB13/GluB21, both chain C termini, and all three B29 macrocycle carboxyls as deprotonated; ArgB22 and the free GlyA1 N terminus are protonated; HisB5 and HisB10 retain the neutral tautomers in the pinned RDKit templates. Under those stated assumptions, the candidate is C325H454N79O97S6 with net formal charge −7.

The ledger enumerates 144 combinations by varying the three macrocycle acids, the GlyA1 N terminus, and each histidine among two neutral tautomers or a protonated state. It assigns no pKa values or state probabilities; the charge range is −8 to −2 under the fixed-group assumptions. Carboxylate charge is recorded at the functional-group level without choosing an oxygen label. This is not an atom-localized partial-charge model, a force field, or an MD-ready state. Tyrosine deprotonation, local pKa shifts, counterions, and solvent remain outside this ledger. Re-run with `python3 research/audit_nnc2215_physiological_charge_states.py`.

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
- GPU acceleration is not needed for source-data calculations. A generic OpenMM 8.6.1 OpenCL preflight ran on a temporary free Colab Tesla T4: three 100 ps 1UBQ runs measured a median 270.175 ns/day. This is not target-system throughput. HunAIFA’s cross-border guidance intake is submitted. See `data/nvidia-engine-preflight.json` for the preflight record.
- The original PaG MD starting coordinates are not public in the supporting information. Any later study using the published NNC2215 models must be labelled as a public-input proxy, not an exact rerun of the PaG simulation.
- The Apple M4 and Colab T4 1UBQ preflights verify that a generic solvated protein system can be constructed and stepped on OpenMM CPU/OpenCL platforms. The T4 timing is not a controlled comparison with Apple M4, does not estimate target-system performance, and does not validate the PaG/NNC2215 setup.
- The new NNC2215 model audit confirms sequence/backbone/disulfide and modification-link records, but flags duplicate atom names and an open/closed linker-attachment naming difference in the nonstandard residues. Chemical graphs and parameters remain unvalidated; see the audit JSON.
- All source data, derived calculations, model-built structures, and future simulation outputs must remain separately labelled.
