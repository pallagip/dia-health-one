# Compute and molecular-simulation readiness

**Checked:** 2026-10-17 21:40 UTC

**Decision:** the published-data analyses run reproducibly on the local CPU, a generic engine preflight ran on Apple M4 OpenCL, and a PaG-only explicit-water setup pilot completed as two short OpenMM OpenCL replicas. The PaG pilot is a small-molecule compatibility check, not an NNC2215–PaG binding study. NVIDIA Inception membership is approved and the Dia Health One Research Dashboard product profile is saved. As checked in the authenticated portal at 2026-10-17 21:40 UTC, three catalog requests are pending: $1,000 Lambda, $10,000 AWS and $5,000 Microsoft Azure. These are not awards or compute allocations. The separate Lambda and AWS partner applications have not been submitted; the Azure partner handoff stopped at a multi-factor sign-in prompt for an identity different from the company account, with no code entered or application submitted. The Innovation Lab tile is absent, and NVIDIA support was contacted on 2026-10-17. No NVIDIA GPU is allocated. The published PaG MD method is documented, but its in-house starting coordinates are not public and force-field readiness for the NNC2215 conjugate remains unresolved. No target binding MD has been made.

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

The reviewed 2026 article/SI and associated public model materials do not include the original seed coordinate file. Study authors have not yet been contacted; requesting the original files is a proposed next step.

## NVIDIA compute access

- Dia Health One's NVIDIA Inception membership was approved on 2026-10-15. Membership does not allocate a GPU.
- NVIDIA's [Innovation Lab](https://www.nvidia.com/en-us/data-center/innovation-lab/) is a separate 60-day program for selected Inception members. A request must be submitted through the Inception Portal; NVIDIA reviews applications on a rolling basis and states that approval depends on availability and fit.
- NVIDIA's current Innovation Lab page explicitly lists molecular dynamics, protein folding, and drug-discovery pipelines among GPU-accelerated simulation use cases. The scientific workload fits the stated use case; selection and allocation remain outstanding.
- The authenticated Inception account now contains the product profile **Dia Health One Research Dashboard** (Developing). The pre-existing computational-insulin product entry was left unchanged.
- The portal's **Benefits → Explore Benefits** catalog was checked after saving the product. Searching “Innovation Lab” returned zero results, so the program application tile is not available in the current account view.
- A concise support request was sent to `inceptionprogram@nvidia.com` on 2026-10-17, identifying the Dia Health One account and product and asking NVIDIA to enable the tile or provide the right application link. No reply was present at the latest inbox check. The Innovation Lab application itself has not been submitted.
- At 2026-10-17 21:40 UTC, the authenticated portal listed three pending requests: `$1,000 in Lambda Cloud Credits`, `$10,000 in AWS Cloud Credits`, and `$5,000 in Microsoft Azure Cloud Credits`. The portal states that each is pending and directs the applicant to the partner for eligibility and approval. This is not evidence of a credit award or usable compute.
- The NVIDIA portal lists the `$1,000 in Lambda Cloud Credits` request as pending and directs the applicant to complete a separate Lambda form. The form at [Lambda's NVIDIA Inception offer page](https://lambda.ai/nvidia-inception) was opened on 2026-10-17 but no values were entered and no application was submitted. It requests applicant contact details, company and role information, countries of residence and payment card, an estimated GPU-hours-per-month value, funding, cloud alternatives, team, and traction; the offer states that a valid payment card is required and that Lambda reviews each application. No card details or personal contact information were sent to Lambda. The Lambda application and credit approval remain incomplete.
- The AWS benefit email directs the applicant to a separate AWS Activate application. No AWS application was submitted, and the company's AWS account and prior-credit eligibility have not been confirmed.
- The Microsoft Azure benefit link was opened to continue through Microsoft for Startups. Its sign-in requested multi-factor authentication for an identity different from Dia Health One's company account. No authentication code was entered and no partner application was submitted; credit approval and compute access remain unconfirmed.
- The current host is macOS 26.3 on Apple M4/arm64. No NVIDIA device or `nvidia-smi` is available. A temporary conda environment was used for the generic engine preflight (OpenMM 8.6.1 CPU/OpenCL and GROMACS 2026.3 CPU-only); this does not provide NVIDIA compute. Desmond is not installed. No NVIDIA instance/job has been started and no compute spend has occurred.
- A read-only `gcloud` CLI check found that this shell is configured for a LabWired account and could not refresh its credentials non-interactively. No Google Cloud project APIs or compute were used; this shell does not currently provide verified Dia Health One GPU access.

The official request path is: **Inception Portal → Profile → Add product → Benefits → Explore Benefits → Innovation Lab → Request Benefit → application form.** The product step is complete; the missing tile is the current access blocker.

## Gates before any target MD run

1. Obtain Innovation Lab access or another explicitly approved compute allocation.
2. Obtain the original PaG MD starting coordinates from the authors or another attributable source. If they remain unavailable, define any public-input proxy as a separate study.
3. Validate charges and force-field parameters for PaG, the macrocycle, glucoside, linker, and insulin conjugate before interpreting energetics.
4. Repeat the OpenMM/GROMACS benchmark in the actual NVIDIA environment. The local 1UBQ preflight is only a portable engine check and is not a substitute for that measurement.
5. Run a small stability pilot with documented diagnostics and independent replicates before considering longer production simulations.

Experimental receptor-binding and pig-clamp results remain literature benchmarks. The existing CPU calculations do not predict or reproduce those assays through molecular dynamics.
