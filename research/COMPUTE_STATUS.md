# Compute and molecular-simulation readiness

**Updated:** 2026-10-26

**Current NVIDIA status:** Dia Health One is an active NVIDIA Inception member. NVIDIA Inception support confirmed membership on 2026-10-19 and stated that DGX Cloud Innovation Lab is no longer available. The last authenticated portal check was 2026-10-18 at approximately 09:59 UTC; it showed six partner-benefit requests pending. This is a dated portal observation, not a fresh status check. Membership and pending requests do not constitute a credit award or GPU allocation. No NVIDIA or partner GPU allocation is confirmed.

**Current simulation status:** OpenFF Sage 2.3.0 labels the source-derived neutral 507-heavy-atom full graph and all eight carboxylate combinations in the isolated roof fragment under the pinned environment. A pinned Amber14 inventory identified literal template-name candidates for 16 of 17 insulin residue labels, six CYX disulfide candidates, and both conjugate joins. A bounded AmberTools 26 GAFF2/AM1-BCC pilot assigned candidate types and charges to the isolated 92-heavy-atom B29 fragment in the conditional −3 state; all eight carboxylate Lewis labelings produced identical MOL2 and frcmod files. The output charge sum is −2.996005 e (0.003995 e from the formal −3), and parmchk2 emitted 26 additional parameter records, including two penalty scores of 41.9. These isolated-fragment outputs need charge reconciliation and expert review; the LysB29 boundary, PheB1 junction, exact mapping, and full hybrid force field remain open. No validated conjugate force field or OpenMM target system exists. Published-data analyses run reproducibly on the local CPU; generic engine preflights ran on Apple M4 OpenCL and a temporary free Colab Tesla T4; and a PaG-only explicit-water setup pilot completed as two short OpenMM OpenCL replicas. None is target-binding MD. The published PaG MD method is documented, but its in-house starting coordinates are not public. No target-binding MD has been run.

The dated access and reproduction notes below are historical observations. Where they describe the Innovation Lab as an available program or say NVIDIA had not replied, the 2026-10-19 NVIDIA support response supersedes them. Google rejected the separate Cloud Program application on 19 October because the company website was inaccessible or inactive at review time; the public site is now reachable and a reprocessing reply was sent that day. The separate Lambda, AWS, Azure, and AI Factory France partner applications had not been submitted at the last recorded review. The Azure handoff stopped at a multi-factor prompt for a different identity; no code was entered. The six NVIDIA benefit requests are partner requests, not awards or allocations.

## OpenFF Sage label-coverage audit — 26 October 2026

OpenFF Sage 2.3.0 label coverage was run in the macOS ARM64 environment recorded by `research/pag-chemistry-environment-macos-arm64.lock.txt` (OpenFF Toolkit 0.19.0, Interchange 0.5.5, Units 0.4.0, OpenFF force fields 2026.09.0, OpenMM 8.6.1, RDKit 2026.03.1). The reconstructed neutral constitutional graph retained formula C325H461N79O97S6 and 507 heavy atoms; after explicit hydrogen materialization it contained 968 atoms. Sage assigned labels to 461 constraints, 988 bonds, 1,753 angles, 2,569 proper torsions, 225 improper torsions, and 968 vdW entries; this method returned zero electrostatics labels. Constraint, bond, angle, proper-torsion, and vdW counts match graph-derived expectations for all nine graphs.

The isolated roof macrocycle fragment also matched those topology-derived expected counts in all eight combinations of its carboxyl protonation states (formal charges 0 through −3). Improper-torsion counts are recorded but not independently checked for completeness. These fragment checks do not estimate pKa values or state populations and do not replace the separate full-system charge-state ledger. The machine-readable result records hashes for the graph, mapped roof SMILES, map legend, lockfile, and audit script at `research/data/nnc2215-forcefield-label-coverage.json`; reproduce it with `micromamba run -p <environment-created-from-lock> python research/audit_nnc2215_forcefield_label_coverage.py`.

This result is parameter-label coverage only. No partial charges, coordinates, OpenMM system, or dynamics were produced. Improper-torsion completeness, atom mapping, roof stereochemistry, a physiological full-system microstate, transferability, and compatibility with the protein force field remain unvalidated. The validated-conjugate-force-field count remains zero; target-binding MD, affinity, and clinical effects remain unestablished.

## Hybrid Amber template feasibility inventory — 26 October 2026

`research/audit_nnc2215_hybrid_feasibility.py` parsed `amber14-all.xml` and its included Amber XML files from the exact pinned OpenMM environment and matched their template names against the source-provenance ledger for all 50 insulin residues. Sixteen of 17 observed labels have literal template names; the two histidines are HIE candidates inferred from their NE2-bound graph hydrogens, and all six cysteines are CYX candidates because they participate in three SG–SG edges. The full graph contains 109 mapped-conjugate heavy atoms and exactly two insulin-to-conjugate joins: PHE B1 N–B1 map 2 and LYS B29 NZ–B29 map 84.

This identifies candidates by name and graph context only. It does not establish exact template matching or terminal/charge compatibility. Modified PheB1 and LysB29 and their conjugate fragments require custom parameters, including cross-boundary terms. The status is `AMBER_TEMPLATE_CANDIDATES_IDENTIFIED; CUSTOM_CONJUGATE_PARAMETERS_REQUIRED`; no charges or parameters were generated, and the validated conjugate-force-field count remains zero. Whole-graph NAGL use would be a large size extrapolation (507 heavy atoms versus two published test molecules in the 80–89 heavy-atom bin and none above it). The AmberTools custom-residue workflow is a precedent requiring manual review, while the Sage/ff14SB Interchange example covers separate components experimentally, not this covalent conjugate.

The checksum-recorded result and XML inventory are in `research/data/nnc2215-hybrid-feasibility.json`. Reproduce from a lock-derived environment with `micromamba run -p <environment-created-from-lock> python research/audit_nnc2215_hybrid_feasibility.py`. This milestone produced no coordinates or OpenMM system and ran no molecular dynamics.

## Isolated B29 GAFF2 parameterization pilot — 26 October 2026

`research/audit_nnc2215_ambertools_parameterization_pilot.py` ran AmberTools 26 Antechamber and parmchk2 on the 92-heavy-atom, source-mapped B29 linker–roof product in a declared tri-carboxylate state (formal charge −3). It reused heavy-atom coordinates from the publisher-supplied active/open model-built candidate, selected one of its 16 graph embeddings, and added explicit hydrogens with RDKit. The eight Lewis labelings of the three carboxylates produced byte-identical GAFF2 MOL2 files and byte-identical frcmod files. All 161 atoms received GAFF2 types. The rounded output partial-charge sum is −2.996005 e, leaving a +0.003995 e residual to the formal charge; no correction was applied. `parmchk2` emitted 26 additional parameter records, with no ATTN marker but two penalty scores of 41.9.

This demonstrates toolchain feasibility for one isolated, capped B29 fragment and indicates that arbitrary carboxylate oxygen labels do not change these outputs under the pinned run. It does not establish final charges or terms: the NH2 endpoint cap is not the protein-bound LysB29 NZ environment, a source-to-coordinate atom map remains ambiguous, two generated improper terms have high penalty scores, and the other protein/conjugate junction is not included. The raw parameter files are retained under `research/data/nnc2215-b29-gaff2-pilot/`; the machine-readable audit records each variant, file hash, charge sum, and parameter record. Reproduce in the locked environment with `micromamba run -p <environment-created-from-lock> python research/audit_nnc2215_ambertools_parameterization_pilot.py`.

The [AmberTools modified-residue tutorial](https://ambermd.org/tutorials/basic/tutorial5/index.php) documents Antechamber, prepgen, and parmchk2 as a candidate workflow and cautions that charges and generated terms need review. This pilot does not create a custom insulin residue library, an OpenMM system, or MD; the validated-conjugate-force-field count remains zero.

## B29 mapping sensitivity audit — 26 October 2026

The follow-on audit ran all 16 exact active/open model-coordinate mappings with all eight carboxylate Lewis labelings per mapping (128 AmberTools cases total), using four isolated workers. Within each mapping, the eight resonance labelings produced one MOL2 hash; all 128 cases shared one frcmod hash. Across mappings, there were 16 MOL2 hashes. GAFF2 atom types were constant for all mapped heavy-atom PDB serials. Partial charges varied at 19 of 92 physical heavy-atom serials, with a maximum observed range of 0.0035 e. The total partial-charge sum ranged from −2.997001 to −2.996005 e (residual +0.002999 to +0.003995 e from formal −3); every result fell within the audit's 0.01 e bookkeeping tolerance, and no correction was applied.

This shows that the generated bonded-term file and heavy-atom types did not change across these candidate mappings, while the AM1-BCC heavy-atom charges retain small mapping dependence. It does not authorize selecting one of the 16 maps, map RDKit-added hydrogens to physical PDB hydrogens, resolve the NH2 cap at LysB29, or validate the partial charges or parameters. The mapped coordinate candidate remained unchanged. The result JSON records each map, variant, output hash, and per-serial charge range; reproduce with `micromamba run -p <environment-created-from-lock> python research/audit_nnc2215_b29_mapping_sensitivity.py --max-workers 4`.

## Reproductions completed

Five Python analyses completed on the Apple M4 host:

1. **Nature 2024 Figure 3d source values:** 30 CSV values were verified against the publisher workbook's source cells and workbook SHA-256. Paired fold changes across three independent experiments are 3.1576× for 3→20 mmol/L glucose and 12.4720× for 0→20 mmol/L, matching the paper's rounded 3.2× and 12.5× values.
2. **2026 PaG pig clamp, Table 1:** the published group means, SEM, sample sizes, and p-values were transcribed. The script calculates descriptive ratios of group means: 1.1466× for low-dose PaG and 1.4226× for high-dose PaG versus saline. Individual animal values are not available in this data set, so SEM and significance tests are not recomputed.
3. **2026 PaG model, Appendix S4 Equation 4b:** the reported `CoeffPaG` value is reapplied at published human and pig concentration anchors. This reproduces only the linear insulin-sensitivity component, not receptor affinity, glucose-infusion rate, or patient response.
4. **2026 PaG binding summaries:** ratios are recalculated from the paper's reported Kd means (1.91 and 0.17 mmol/L) and zero-glucose hIR-A affinity means (0.9% and 5.4% with 70 μmol/L PaG). These are arithmetic checks on rounded published summaries; the native-MS and receptor-binding curves are not refitted.

5. **Nature 2024 Figure 3d descriptive Hill fit:** the 15 NNC2215 source-workbook values were fit with a shared apparent slope and EC50, measured experiment-specific zero-glucose baselines, and experiment-specific upper responses. The best fit is EC50 13.43 mmol/L and apparent Hill slope 1.115. In 2,000 deterministic bootstrap resamples of the three independent experiments, the EC50 2.5th percentile is 8.06 mmol/L; 64 resamples (3.2%) hit the 1,000 mmol/L optimization ceiling, censoring the upper percentile. The apparent-slope percentile range is 0.675–1.491. This fit establishes that the available points do not identify a reliable EC50; it does not establish a 6 mmol/L switch threshold. Percentile ranges are fragile with three clusters and are not formal inferential confidence intervals.

The run records and input checksums are in `data/reproduction.json`, `data/pag-clamp-reproduction.json`, `data/pag-sensitivity-reproduction.json`, `data/pag-binding-summary-reproduction.json`, and `data/hill-fit-summary.json`. The Hill-fit script is `fit_nnc2215_hill_response.py`; it regenerates `assets/hill-fit.svg`. `data/nvidia-access-status.json` records the compute-access state and source links.

## Local engine preflight

The current run record compares GROMACS CPU, OpenMM CPU, and OpenMM OpenCL using the same solvated 1UBQ ubiquitin coordinates (26,821 atoms). This Apple M4 run measured 37.420 ns/day for GROMACS CPU, 11.141 ns/day for OpenMM CPU, and 126.274 ns/day for OpenMM OpenCL. These are single-run timings, with different thermostat implementations and timing boundaries; use them as installation and rough throughput checks only. They do not compare NVIDIA hardware or reproduce any PaG/NNC2215 result.

The run record contains exact versions, setup parameters, host details, input and script hashes, and interpretation limits at `data/engine-preflight.json`. Recreate it with the pinned environment and `benchmark_engine_preflight.py`; the chart at `assets/engine-preflight.svg` is generated from that result. The input is public RCSB PDB 1UBQ. GROMACS 2026.3 was built without GPU support; OpenMM 8.6.1 ran on CPU and OpenCL. No target MD was started.

## Temporary NVIDIA T4 engine preflight

A temporary free Google Colab runtime provided a Tesla T4 (15,360 MiB, driver 580.82.07). OpenMM 8.6.1 used its OpenCL platform; no CUDA plugin was present. Three independent 100 ps timing windows on a solvated 1UBQ system (17,764 atoms; Amber14-all/TIP3P; PME; 2 fs; 300 K; single precision) measured 282.691, 270.175 and 257.885 ns/day. The median is 270.175 ns/day (sample SD 12.403). Device synchronization via `getState(getEnergy=True)` occurred inside each timed interval. The detailed run record, environment, seeds, energies, source checksum and limitations are in `data/nvidia-engine-preflight.json`; the chart is generated by `plot_nvidia_t4_preflight.py`.

This establishes that a generic protein system runs on the available T4/OpenCL stack. It does not estimate throughput for the much larger NNC2215–PaG or receptor-complex systems, validate nonstandard conjugate parameters, compare directly with the Apple M4 run, reproduce the published 500 ns calculation or measure binding. The free Colab session is temporary and is not an NVIDIA Inception or partner allocation. No project compute spend was recorded. The saved notebook is [Dia Health One — PaG GPU Reproducibility Preflight](https://colab.research.google.com/drive/1ImHF0tk_wGSWfPwtuZjC04cFN7i4fCmt?authuser=1).

## PaG-only explicit-water pilot

The arm64 chemistry environment uses RDKit 2026.03.1, AmberTools 26.0, OpenFF Toolkit 0.19.0, OpenFF Sage 2.3.0, and OpenMM 8.6.1. It builds a PaG carboxylate conformer from the PubChem CID 83944 identity record, retains the five stereocentres, assigns AM1-BCC charges, and verifies an OpenFF ligand system with net charge −1. This does not parameterize the NNC2215 macrocycle, linker, glucoside, or insulin conjugate.

Two independent OpenCL replicas completed, each with 10 ps equilibration and 100 ps production at 298.15 K and 1 atm in 0.15 M NaCl. Each solvated box contains 1,272 atoms and one PaG; its implied 0.118 M ligand concentration is a finite-box property and does not reproduce the 70 μM assay. Production mean temperatures were 298.85 K and 297.54 K, with finite energy and volume traces. This is a setup smoke test only: no NNC2215 or binding partner is present, and no affinity or published PaG simulation result is reproduced. The run manifest, raw logs, final coordinates, analysis JSON, and figure are linked from the research dashboard and indexed in `data/pag-solvation-preflight/latest.json`.

## Structure data and PaG MD method

Springer Nature Figshare publishes author-supplied open and closed NNC2215 model PDBs under CC BY 4.0. The files are included in `structures/`; source and local checksums are recorded in `data/structure-manifest.json`. These are model-built coordinates, not experimentally resolved NNC2215 structures.

`inspect_structure_models.py` verifies those checksums and records PDB connectivity metadata in `data/structure-preflight.json`: the open model has 990 coordinate records and the closed model 968, each with three SSBOND and three LINK records. This only inventories file records; it does not validate chemistry or parameters.

The follow-up audit in `data/nnc2215-model-validation.json` independently matches the expected desB30 insulin A/B sequences and coordinate sequences, checks complete backbone atoms and the three native disulfides, and finds the expected B29/B1 modification-site records. It also finds repeated atom labels in nonstandard PHX, LYX, and MC residues. The open model connects MC C31 to LYX CD6, while the closed model connects MC C31 to LYX C. These discrepancies may reflect label/representation differences, but must be reconciled against authoritative chemistry before assigning topology. Bond orders, complete conjugate stereochemistry, protonation, partial charges, and force-field parameters remain unvalidated. The PDB models do not include insulin-receptor coordinates. Recreate the checks with `python3 research/validate_nnc2215_models.py`.

The 2026 Supporting Information reports Desmond 6.6 (Schrödinger Release 2021-2), OPLS4, SPC explicit water, 0.15 mol/L NaCl, and 500 ns production. The input is an in-house crystal structure of a related glucose-bound macrocycle. The coordinates are described as in-house and are not supplied. The public NNC2215 state models are different inputs and cannot be substituted silently for an exact rerun. The recipe is recorded in `data/pag-md-protocol.json`.

A related public model source is now traced: the 2024 NNC2215 article cites Tromans et al. (2019) as the source of its macrocycle model. The 2019 Supporting Information describes Maestro MCMM models of a receptor surrogate and its glucose complex, but the reviewed supplement contains figures and methods rather than coordinate files. This narrows the search history; it does not supply the 2026 crystal seed or validate a replacement parameterization.

The current milestone report is `MILESTONE_REPORT.md`; it separates completed source-data calculations, structure-audit results, current blockers, and the evidence required for each validation gate. A proposed, not-run pilot matrix and stop criteria are in `data/target-md-pilot-protocol.json`. It uses the public NNC2215 model as a proxy input and does not claim an exact rerun of the 2026 macrocycle-analogue simulation.

The reviewed 2026 article/SI and associated public model materials do not include the original seed coordinate file. Patrik said on 17 October that he would contact the study authors; no author reply or coordinates had been received by the 18 October 10:01 UTC inbox check.

## NVIDIA compute access

- Dia Health One's NVIDIA Inception membership was approved on 2026-10-15. Membership does not allocate a GPU.
- NVIDIA's [Innovation Lab](https://www.nvidia.com/en-us/data-center/innovation-lab/) is a separate 60-day program for selected Inception members. A request must be submitted through the Inception Portal; NVIDIA reviews applications on a rolling basis and states that approval depends on availability and fit.
- NVIDIA's current Innovation Lab page explicitly lists molecular dynamics, protein folding, and drug-discovery pipelines among GPU-accelerated simulation use cases. The scientific workload fits the stated use case; selection and allocation remain outstanding.
- The authenticated Inception account contains two product records: **Dia Health One Research Dashboard** and **Dia Health computational insulin discovery platform**, both Developing / Not Accelerated. This meets the published prerequisite of at least one product record.
- The portal's **Benefits → Explore Benefits** catalog was checked after saving the product. Searching “Innovation Lab” returned zero results, so the program application tile is not available in the current account view.
- A concise support request was sent to `inceptionprogram@nvidia.com` on 2026-10-17, identifying the Dia Health One account and product and asking NVIDIA to enable the tile or provide the right application link. No reply was present at the latest inbox check. The Innovation Lab application itself has not been submitted.
- At 2026-10-17 21:40 UTC, the authenticated portal listed three pending requests: `$1,000 in Lambda Cloud Credits`, `$10,000 in AWS Cloud Credits`, and `$5,000 in Microsoft Azure Cloud Credits`. The portal states that each is pending and directs the applicant to the partner for eligibility and approval. This is not evidence of a credit award or usable compute.
- The NVIDIA portal lists the `$1,000 in Lambda Cloud Credits` request as pending and directs the applicant to complete a separate Lambda form. The form at [Lambda's NVIDIA Inception offer page](https://lambda.ai/nvidia-inception) was first opened on 2026-10-17 19:50 UTC; at that earlier check no values had been entered and no application was submitted. A later 2026-10-17 22:05 UTC review showed a redacted phone-field value on page load; this was not entered manually, the form was not submitted, and the review tab was closed. It requests applicant contact details, company and role information, countries of residence and payment card, an estimated GPU-hours-per-month value, funding, cloud alternatives, team, and traction; the offer states that a valid payment card is required and that Lambda reviews each application. No card details or personal contact information were sent to Lambda. The Lambda application and credit approval remain incomplete.
- The AWS benefit email directs the applicant to a separate AWS Activate application. No AWS application was submitted, and the company's AWS account and prior-credit eligibility have not been confirmed.
- The Microsoft Azure benefit link was opened to continue through Microsoft for Startups. Its sign-in requested multi-factor authentication for an identity different from Dia Health One's company account. No authentication code was entered and no partner application was submitted; credit approval and compute access remain unconfirmed.
- The local host is macOS 26.3 on Apple M4/arm64; no local NVIDIA device or `nvidia-smi` is available. A separate temporary free Colab runtime supplied a Tesla T4 for the generic preflight above. No NVIDIA Inception or partner allocation, paid cloud job, target-system run or compute spend is recorded. A temporary conda environment was used for the Apple M4 engine preflight (OpenMM 8.6.1 CPU/OpenCL and GROMACS 2026.3 CPU-only). Desmond is not installed.
- A read-only `gcloud` CLI check found that this shell is configured for a LabWired account and could not refresh its credentials non-interactively. No Google Cloud project APIs or compute were used; this shell does not currently provide verified Dia Health One GPU access.

The official request path is: **Inception Portal → Profile → Add product → Benefits → Explore Benefits → Innovation Lab → Request Benefit → application form.** The product step is complete; the missing tile is the current access blocker.

## Gates before any target MD run

1. Obtain Innovation Lab access or another explicitly approved compute allocation.
2. Obtain the original PaG MD starting coordinates from the authors or another attributable source. If they remain unavailable, define any public-input proxy as a separate study.
3. Validate charges and force-field parameters for PaG, the macrocycle, glucoside, linker, and insulin conjugate before interpreting energetics.
4. Repeat the OpenMM/GROMACS benchmark in the actual NVIDIA environment. The local 1UBQ preflight is only a portable engine check and is not a substitute for that measurement.
5. Run a small stability pilot with documented diagnostics and independent replicates before considering longer production simulations.

Experimental receptor-binding and pig-clamp results remain literature benchmarks. The existing CPU calculations do not predict or reproduce those assays through molecular dynamics.


## Latest evidence review — 18 October 2026, 00:05 CEST

- The authenticated NVIDIA portal still lists the $1,000 Lambda, $10,000 AWS, and $5,000 Microsoft Azure benefit requests as pending; none is a partner award or GPU allocation. `nvidia-access-status.json` now records this exact check time.
- The public NVIDIA Benefits Catalog search for “Innovation Lab” returned 0 results. The separate support request has no reply in the checked inbox. Opening the known portal application URL showed only the generic Inception sign-in/application page, not an Innovation Lab application. No additional request was submitted.
- The official Lambda form was reviewed at 00:05 CEST and later closed at 00:09 CEST. It requires a valid payment card and applicant-specific phone, residence/card countries, a six-month GPU-hours estimate, ML-team/project details, and current traction metrics. No application-specific values were entered manually and no submission was made. A redacted value appeared in the phone field on page load; it remained redacted after a clear attempt, so the form tab was closed without submission.
- Two other NVIDIA-listed GPU-cloud offers were screened. Nebius requires VC/accelerator backing and at least $5 million raised, which conflicts with the company’s recorded bootstrapped $0 funding. Scaleway lists a 3–5 year operating-history requirement plus active EU go-to-market; those eligibility facts are not established. Neither provider request was submitted.
- The Nature supplementary methods document the B29 alkyne/macrocycle click route, B1 glucoside coupling, and final formula C325H461N79O97S6. Element inventories from the model PDB residues PHX and LYX differ from the unrelated RCSB CCD entries with those same three-letter IDs. This confirms that those dictionary entries cannot be used as topology definitions for these model files. See `data/nnc2215-chemistry-source-audit.json`.
- These checks do not resolve the atom-mapped conjugate graph or the 2026 paper’s in-house seed coordinates. No target binding MD was started, no Inception/partner allocation was used, and no credit was awarded. The temporary free Colab T4 was used only for the generic benchmark documented below.

## Update — 18 October 2026, 00:45 CEST

A synchronized OpenMM 8.6.1 OpenCL preflight completed on the temporary free Colab Tesla T4. The 1UBQ benchmark used 17,764 atoms and three independent 100 ps timing windows; throughput was 282.691, 270.175 and 257.885 ns/day, with a median of 270.175 ns/day. The timing included an energy-state query to synchronize GPU work. The earlier quick notebook measurement did not include that synchronization and is superseded. The saved Colab notebook and machine-readable local manifest preserve the run; the dashboard source now references a chart generated from that manifest. This is a generic engine preflight only: no NNC2215, PaG, receptor complex, binding calculation, target simulation or paid compute was involved. The T4 came from a temporary free Colab runtime, not an Inception or partner allocation. Target conjugate topology, validated parameters, the original 2026 starting coordinates, and approved production compute remain open gates.

## Public 2026 PaG supplementary-file audit — 18 October 2026, 00:54 CEST

The ACS Figshare API record for article 32292700 was queried directly. It lists one file only: pt6c00069_si_001.pdf (2,870,120 bytes; 24 pages). The downloaded file's MD5 matches Figshare's supplied checksum (c87fafa63a891bd1f5aec832e51182ea); SHA-256 is f3d279ef9169666cef211ad9849fa0ea1ca807820eb4e3993944417b51936db2. Appendix S2 says the authors used Desmond 6.6, OPLS4, SPC water, 0.15 mol/L NaCl and 500 ns production for an in-house glucose-bound macrocycle similar to NNC2215, then used a representative glucose-complex snapshot to build a macrocycle–PaG starting configuration in Maestro and applied the same MD protocol. The public item lists no coordinate, trajectory, extracted snapshot, or Desmond project file. This establishes the missing artifacts for the exact rerun in the checked public item; it does not establish that the authors do not retain them. Machine-readable evidence is in research/data/pag-si-source-audit.json.

## NVIDIA access recheck — 18 October 2026, 00:53 CEST

The authenticated Inception Benefits page still shows the $1,000 Lambda, $10,000 AWS, and $5,000 Microsoft Azure requests as pending. The public NVIDIA Benefits Catalog search returned zero results for “Innovation Lab.” No NVIDIA partner GPU allocation was visible. The free Colab T4 remains a generic preflight only. No target MD or project compute spend is recorded.

## NNC2215 connectivity audit — 18 October 2026, 01:10 CEST

The new `audit_nnc2215_connectivity.py` audit resolves the PDB record-level evidence more precisely. Open and closed model atom serials are unique, all three named LINK endpoints resolve uniquely in each file, and duplicate atom labels can be associated with distinct serial-level CONECT neighbors. The open model has 11 duplicate-name groups and the closed model 12; the closed model additionally has an unconnected PHX H2 atom, serial 332. The macrocycle attachment differs across models: open MC C31–LYX CD6, closed MC C31–LYX C. Measured distances for both records are consistent with their PDB LINK annotations.

The wwPDB specification provides connectivity by atom serial and does not encode bond orders or distinguish formal charge states. An exploratory RDKit bond-order perception attempt on the 164-atom modified HET component did not complete within roughly one minute and was interrupted; no inferred assignments were retained. The audit therefore improves the adjacency map but does not validate bond order, stereochemistry, protonation, charges, or force-field parameters. Machine-readable findings are in `data/nnc2215-connectivity-audit.json`. Recreate with `python3 research/audit_nnc2215_connectivity.py`.

The user-provided Drive folder currently contains the progress record, company handbook, EIN letter, and incorporation certificate, but no raw molecular source data or simulation files. Research code, source extracts, structures, manifests, and generated figures currently live in the GitHub research repository. No target binding MD or NVIDIA-allocated compute was used for this audit.


## Compute access update — 18 October 2026, 02:54 CEST

The refreshed NVIDIA Inception Benefits page shows five pending catalog requests: Lambda ($1,000), AWS ($10,000), Microsoft Azure ($5,000), AI Factory France (up to 50K GPU hours + AI services), and Google Cloud startup credits ($2,000–$350,000). NVIDIA acknowledgment emails received at 00:48:56 UTC and 00:53:38 UTC confirm the last two requests were registered. Each email says the partner handles eligibility and that a separate partner application is required. These are requests, not awards; no GPU allocation or project compute spend is confirmed.

The public NVIDIA Benefits Catalog search for “Innovation Lab” returned 0 results on the latest check. The previous support request still has no reply in the point-in-time inbox search. Private communications and inbox-monitoring details are recorded separately from this public compute-status report.

### AI Factory France route

The NVIDIA request link opened AI Factory France’s services page; its computing-services section links to the DARI access portal. No application was submitted. AI Factory France’s [FAQ](https://aifactory-france.eu/faq/) says its free compute access requires the applicant’s R&D department to be based in France and the work to be open research, with a commitment to publish the scientific results. The offer materials allow a possible cross-border path only by agreement between AI Factory France and the applicant country’s AI Factory. Dia Health One’s R&D location and publication position are not yet confirmed. EuroHPC lists [HunAIFA](https://www.eurohpc-ju.europa.eu/ai-factory-antennas_en) as Hungary’s AI Factory Antenna, but referral eligibility for this specific benefit is unconfirmed.

### Hungary AI Factory support request — 18 October 2026

A support intake was submitted through [HunAIFA](https://hunaifa.eu/) asking for guidance on Hungarian and European AI Factory GPU/HPC eligibility, possible referral routes, and OpenMM/GROMACS expertise. The confirmation page stated that the responses were recorded. HunAIFA has not replied yet. This is a support request only: no AI Factory referral, compute approval, or GPU allocation is confirmed. The NVIDIA AI Factory France catalog request remains pending separately.

### Google Cloud startup credits — status update 19 October 2026

Google’s first-step email arrived at 2026-10-15 20:05:16 UTC and the company-account application was received at 20:27:42 UTC. On 19 October at 13:50:30 UTC, Google rejected the application because the company website was inaccessible or inactive when reviewed. Google explicitly invited a reply after the website is active and publicly accessible; a reply requesting re-processing was sent on 19 October. The NVIDIA $2,000–$350,000 catalog request remains pending separately; neither is a credit award. The public [eligibility page](https://cloud.google.com/startup/faq) says the award and tier are decided by Google. No credits or GPU allocation are confirmed.

### Other EuroHPC routes

- [Development Access](https://www.eurohpc-ju.europa.eu/eurohpc-ju-call-proposals-development-access_en) is a possible near-term route for code/algorithm development, testing, and optimization. The call lists a 1 November 2026 cut-off and up to 2–3 weeks from submission to access; a defensible scope here would be OpenMM/GROMACS porting and benchmarks, not target binding MD.
- [Fast Lane](https://www.eurohpc-ju.europa.eu/ai-factories/ai-factories-access-modes/fast-lane-access-ai-factories_en) advertises up to 50,000 GPU-hours over up to three months, but explicitly targets SMEs/startups performing AI activities. The current target is classical MD and has no established AI component, so fit is not assumed.
- Regular Access is open to industry and scientific fields but is intended for larger allocation needs and has a longer review-to-access path; it is a later possibility after validating the system and resource estimate.

No target NNC2215–PaG binding MD started during this access review. Full conjugate connectivity, protonation/charge and parameters remain unvalidated, and the 2026 in-house simulation seed coordinates are not present in the reviewed public SI. The temporary free Colab T4 benchmark remains a generic 1UBQ engine preflight only.


## Latest access and offer review — 18 October 2026, 03:45 UTC

The authenticated Benefits page at 03:37:56 UTC still showed the five pending partner requests listed above. The two product records satisfy the published one-product prerequisite. The Explore Benefits route opens the public catalog; its “Innovation Lab” search displayed zero results at 03:45:48 UTC, so no Innovation Lab application was submitted.

Nebius was screened out because its listed GPU-credit offer requires at least $5 million raised and VC/accelerator backing, while the current profile does not meet the funding threshold. Scaleway lists up to €42,000 in credits, but the NVIDIA listing targets a post-product scaling stage and the company’s age, account status, and other eligibility facts are not verified. Its public startup form requires verified identity. No fields or identity details were entered and no application was submitted.

No NVIDIA or partner allocation, project compute spend, or target NNC2215–PaG MD is confirmed. Current source and account evidence is in data/nvidia-access-status.json.


## Portal refresh — 18 October 2026, approximately 03:53 UTC

The authenticated NVIDIA Inception Benefits page was reloaded and now lists six pending catalog requests, adding €40,000 in Scaleway Cloud Credits to the five requests recorded at 03:37 UTC. NVIDIA’s 03:44:15 UTC receipt confirms the Scaleway catalog request and says partner processing is handled by Scaleway; the direct partner application may still be required. The separate Scaleway form has not been submitted and no identity details were transmitted. This is a pending request, not an award or GPU allocation.

The NVIDIA receipts for AI Factory France and Google Cloud also instruct applicants to complete partner processes separately. The Google for Startups Cloud Program application was separately received on 15 October and remains pending. The public Innovation Lab search still returned zero results at 03:45:48 UTC, so no Innovation Lab application was submitted. No target NNC2215–PaG MD has run.


## Referral intake status correction — 18 October 2026, 03:55 UTC

The HunAIFA support intake is submitted. A return visit to the form displays “Thank you! Your survey responses have been recorded.” The 03:15:40 UTC project email to Patrik also records the submission and confirmation. This is a request for Hungarian/European AI Factory access guidance and a possible referral; it is not an AI Factory France compute application, an approval, or an allocation.


## Current portal and analysis recheck — 18 October 2026, 07:15 UTC

A fresh authenticated NVIDIA Inception Benefits page showed the same six requests as pending: Lambda $1,000, AWS $10,000, Microsoft Azure $5,000, AI Factory France up to 50,000 GPU-hours, Google Cloud $2,000–$350,000, and Scaleway €40,000. No GPU allocation was shown. The Innovation Lab search was not repeated in this pass; its last recorded result was zero at 04:54 UTC.

The four published-data analysis scripts were rerun locally at 07:15 UTC. Nature Figure 3d paired fold changes remain 3.157629× (sample SD 0.278646) from 3 to 20 mmol/L and 12.472023× (sample SD 1.698501) from 0 to 20 mmol/L. The PaG summary ratio remains 11.235× lower reported Kd than glucose and the rounded zero-glucose receptor-affinity ratio remains 6.000× at 70 μmol/L. Pig-clamp ratios remain 1.147× (low-dose PaG, p=0.16), 1.423× (high-dose PaG, p=0.04), and 1.324× (L-glucose, p=0.01). Equation 4b changes remain +54.5%, +79.4%, and +331.6% at the reported anchors. These are source-data calculations on CPU, not target MD.

A read-only local gcloud check still resolves to the LabWired account/project rather than Dia Health One; billing-account enumeration failed at non-interactive reauthentication. No Google Cloud project or billing setting was changed and no cloud job was launched. Gmail searches through approximately 07:13 UTC found no newer Patrik reply, NVIDIA Innovation Lab response, or HunAIFA/AI Factory reply. Target NNC2215–PaG MD remains unstarted.

## NVIDIA access refresh — 18 October 2026, 07:26 UTC

The authenticated Inception Benefits page was reloaded at approximately 07:22 UTC. All six catalog requests remain pending—Lambda, AWS, Azure, AI Factory France, Google Cloud, and Scaleway—and no GPU allocation is shown. NVIDIA receipts for AI Factory France (00:48 UTC), Google Cloud (00:53 UTC), and Scaleway (03:44 UTC) confirm that each partner handles its own review and may require a separate partner-site application. The separate Google for Startups Cloud Program application is already received and pending; no credits or compute have been awarded.

NVIDIA’s official Innovation Lab page was rechecked at 07:24 UTC. It describes a two-month program for selected Inception members and lists molecular dynamics and drug-discovery pipelines among its GPU simulation use cases. The account has two product records and meets the published one-product prerequisite, but the Innovation Lab benefit tile remains unavailable. The last successfully completed catalog search returned zero results at 04:54 UTC. A follow-up search at 07:22 UTC remained on a loading indicator, so it is not recorded as a new zero-result observation. No application or allocation is confirmed.

Gmail searches through 07:26:50 UTC found no reply to the Innovation Lab support request, no newer incoming Patrik email than 17 October 21:01 UTC, and no HunAIFA/AI Factory response. Target NNC2215–PaG binding MD remains unstarted: the full conjugate graph, parameters and defensible PaG pose are still unresolved. No NVIDIA or partner GPU allocation is confirmed.


## Live dashboard reproduction and access recheck — 18 October 2026, 07:34 UTC

The deployed dashboard’s “Run reproduction” action was executed in the browser. Run `DHO-NNC2215-FIG3D-20261018073223` used Browser JavaScript arithmetic and the public Figure 3d input (SHA-256 `b7d5bddb7c6e5d2e1b3f016e98f3f18c073681654d23b9a314d6a91c4d13ba4f`). It returned 3.157629 ± 0.278646 for the paired 3→20 mmol/L fold change and 12.472023 ± 1.698501 for 0→20 mmol/L; both match the paper after rounding. The downloadable run manifest is archived at `research/data/browser-runs/DHO-NNC2215-FIG3D-20261018073223.json` and linked from the dashboard. This is a CPU-side source-data calculation, not a GPU run or target binding simulation.

The authenticated Inception page was rechecked at approximately 07:31 UTC: all six catalog requests remain pending and no GPU allocation is shown. The public Benefits Catalog search was applied at 07:33 UTC and visibly returned zero Innovation Lab results; the earlier 07:22 loading attempt is retained as incomplete historical evidence. Gmail searches through 07:34 UTC found no new Patrik message, NVIDIA Innovation Lab support reply, or HunAIFA/AI Factory response. Target NNC2215–PaG MD remains unstarted.


## NVIDIA access refresh — 18 October 2026, 08:06 UTC

The signed-in Inception Benefits page was refreshed at approximately 08:06 UTC. All six catalog benefit requests—Lambda, AWS, Azure, AI Factory France, Google Cloud, and Scaleway—remain pending; no GPU allocation is shown. The public Benefits Catalog search for “Innovation Lab” displayed “Showing 0 Results,” so no Innovation Lab application route is available in the current catalog. A Gmail search through 08:04 UTC found no reply from NVIDIA to the support request. These are partner-benefit requests, not compute awards. Target NNC2215–PaG binding MD remains unstarted.


## Published-data rerun — 18 October 2026, 08:14 UTC

The live browser reproduction completed as `DHO-NNC2215-FIG3D-20261018081424`. The Nature Figure 3d CSV checksum is `b7d5bddb7c6e5d2e1b3f016e98f3f18c073681654d23b9a314d6a91c4d13ba4f`; the browser script checksum is `8883526bcf7f853cefa07dbb2fab5e23e2f0dd15a021a4177cc096148b3b4695`. Paired experiment ratios were 3.157629 ± 0.278646 (3→20 mmol/L) and 12.472023 ± 1.698501 (0→20 mmol/L), matching the published rounded 3.2× and 12.5×. The archived manifest is `research/data/browser-runs/DHO-NNC2215-FIG3D-20261018081424.json`. This was browser-side CPU arithmetic over source data; no GPU was used and target binding MD remains unstarted.

## Local reanalysis and compute check — 18 October 2026, 08:35 UTC

Four literature-analysis scripts were rerun on local Python 3.14.4 / macOS arm64 CPU. The current Figure 3d run is `DHO-NNC2215-FIG3D-20261018T083511Z`; paired fold changes are 3.157629 ± 0.278646 (3→20 mmol/L) and 12.472023 ± 1.698501 (0→20 mmol/L), consistent with the paper's rounded 3.2× and 12.5×. The same run batch recalculated PaG summary ratios (11.235× Kd ratio and 6.000× zero-glucose receptor-affinity ratio), descriptive clamp group-mean ratios, and the SI Equation 4b component. Artifacts are the refreshed `data/reproduction.json` and `data/pag-*-reproduction.json` records. These are CPU arithmetic on published numbers and group summaries.

The current host has no NVIDIA GPU or `nvidia-smi`, `nvcc`, or GROMACS CLI; the browser inventory contained no active Colab notebook tab or new GPU runtime. No target NNC2215–PaG binding MD or free-energy calculation has run. NVIDIA/partner GPU access is not confirmed. Structure and force-field gates remain open, including unresolved B29-linked topology, stereochemistry, protonation/charges, parameters, and absent public starting coordinates for the 2026 in-house MD study.

## NVIDIA access refresh — 18 October 2026, 08:48 UTC

The authenticated Inception Benefits page was rechecked at approximately 08:44 UTC and continues to list six pending partner-benefit requests (Lambda, AWS, Microsoft Azure, AI Factory France, Google Cloud, and Scaleway); no NVIDIA GPU allocation appears. The signed-in Products page shows two company products, both Developing / Not Accelerated, so the published prerequisite of at least one product is met.

The Inception Portal’s Explore Benefits control opened the public NVIDIA Benefits Catalog. After applying the search and pressing Enter, the catalog displayed “Showing 0 Results” for “Innovation Lab.” NVIDIA’s official page says the application is made through the Innovation Lab tile in the Inception catalog, and describes the program as selective and 60 days long. It also lists molecular dynamics, protein folding, and drug-discovery pipelines among its simulation use cases. The account therefore meets the published one-product prerequisite and the research fits a named use case, but there is no available application tile or Request Benefit action in the account view. No application was submitted. The previously sent support request had no reply in the inbox check through approximately 08:48 UTC.

Patrik’s 17 October email reported an invalid short GitHub reference. The correct repository path is `research/MILESTONE_REPORT.md` in `pallagip/dia-health-one`, branch `feat/reproducibility-dashboard`. At 08:44 UTC the GitHub report page and deployed Markdown report both returned HTTP 200. At 08:45 UTC, a reply provided the complete GitHub, raw Markdown, and deployed report links and attached `MILESTONE_REPORT.md`; Patrik was both recipient and copied address.

No target NNC2215–PaG binding simulation or free-energy calculation has run. No NVIDIA or partner GPU allocation is confirmed, and no project compute spend is recorded. This portal check does not change the chemistry and parameterisation gates for target MD.

## Current source-data and GPU checkpoint — 18 October 2026, 10:01 UTC

Run `DHO-NNC2215-FIG3D-20261018T095820Z` reran the four source-data scripts on the Apple M4 CPU. All 30 Nature Figure 3d cells were verified again against a fresh download of the publisher workbook; the paired glucose-response fold changes remain 3.157629 ± 0.278646 (3→20 mmol/L) and 12.472023 ± 1.698501 (0→20 mmol/L). The run archive, source/script checksums, and output are indexed in `data/analysis-runs/` and `data/reproduction.json`. PaG ratios and Equation 4b values remain calculations from rounded or group-level published summaries.

The Benefits portal still shows six pending partner requests and no GPU allocation. The last completed Innovation Lab search returned zero results at 09:53 UTC; no application or GPU award is confirmed. Gmail checks through 10:01 UTC found no newer Patrik message and no NVIDIA support reply. A new Colab T4 connection attempt paused at the interactive reCAPTCHA; no PaG portability job started. The previous generic T4 benchmark remains the only NVIDIA-device simulation record, and it is not a target result.

The exact 2026 study rerun remains unavailable from reviewed public files because its in-house starting coordinates/system are not deposited. Patrik said he would contact the authors; no reply or coordinates had been received by this check. The NNC2215 model-built conjugate still lacks validated bond orders, stereochemistry, protonation, charges, and force-field parameters. Target binding MD remains unstarted.

## Scientist handoff preparation — 18 October 2026, 10:49 UTC

Prepared `SCIENTIST_INPUT_CHECKLIST.md` and added a concise waiting-on-study-team panel to the research dashboard. The checklist specifies the minimum files needed to reproduce the published macrocycle-analogue MD: glucose-bound crystal coordinates plus system/topology, the exact trajectory frame and Maestro/build project used to place PaG, PaG structure/protonation/parameter files, complete equilibration/production settings, and access to the reported Desmond 6.6 / OPLS4 stack. It distinguishes an exact rerun from an OpenMM/GROMACS port and from a separate proxy study using the 2024 model-built NNC2215 structures.

A targeted Gmail review through 10:49 UTC found no new author reply or source files. Patrik's latest inbound email remains 17 October 21:01 UTC and says he will contact the authors. No message was sent to the authors in this update. No target MD was started; missing source inputs and unresolved conjugate chemistry/parameter validation remain scientific gates.
