# Compute and molecular-simulation readiness

**Updated:** 2026-10-18 10:14 CEST (08:14 UTC; latest browser reproduction 08:14 UTC; NVIDIA portal and catalog rechecked at 08:06 UTC)

**Decision:** the published-data analyses run reproducibly on the local CPU, a generic engine preflight ran on Apple M4 OpenCL, a PaG-only explicit-water setup pilot completed as two short OpenMM OpenCL replicas, and a generic OpenMM preflight completed on a temporary free Colab Tesla T4. The PaG pilot is a small-molecule compatibility check, not an NNC2215–PaG binding study; the T4 run is a generic ubiquitin throughput check. Neither is target binding MD. NVIDIA Inception membership is approved and the Dia Health One Research Dashboard product profile is saved. The refreshed authenticated portal at 09:31 CEST / 07:31 UTC still lists six pending benefit requests: $1,000 Lambda, $10,000 AWS, $5,000 Microsoft Azure, AI Factory France up to 50K GPU hours, $2,000–$350,000 Google Cloud, and €40,000 Scaleway. The public Benefits Catalog search for Innovation Lab returned zero results at 09:33 CEST / 07:33 UTC. These requests are not awards or compute allocations. The Google for Startups Cloud Program application was submitted on 15 October and is pending review. The authenticated company-account page confirms receipt; no credits or GPU allocation are awarded. The separate Lambda, AWS, Azure, and AI Factory France partner applications have not been submitted. The Azure handoff stopped at a multi-factor prompt for a different identity; no code was entered. The Innovation Lab tile was absent in the last catalog search, and NVIDIA support was contacted on 2026-10-17. No NVIDIA Inception or partner GPU is allocated. The published PaG MD method is documented, but its in-house starting coordinates are not public and force-field readiness for the NNC2215 conjugate remains unresolved. No target binding MD has been run.

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

The reviewed 2026 article/SI and associated public model materials do not include the original seed coordinate file. Study authors have not yet been contacted; requesting the original files is a proposed next step.

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

### Google Cloud startup credits

Google’s first-step email arrived at 2026-10-15 20:05:16 UTC and directed the applicant to apply for credits. Google then confirmed receipt of the Cloud Program application at 20:27:42 UTC and said its review status would follow in a few business days. On 18 October at about 03:20 UTC, the authenticated Dia Health One company account page displayed “Thank you! We’ve received your application” and said a response is expected within 3–5 business days. Therefore the Google for Startups application is submitted and pending review. The NVIDIA $2,000–$350,000 catalog request also remains pending; neither is a credit award. An earlier read-only inspection under a mismatched personal Google account was left untouched and did not submit personal or mismatched values. The application confirmation does not expose the billing account ID it will credit. The public [eligibility page](https://cloud.google.com/startup/faq) says the award and tier are decided by Google. No credits or GPU allocation are confirmed.

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
