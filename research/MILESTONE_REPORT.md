# NNC2215–PaG reproducibility milestone report

**Date:** 2026-10-18 00:45 CEST (2026-10-17 22:45 UTC)
**Status:** Published-data reanalysis and Hill-fit identifiability check complete; a generic OpenMM throughput preflight completed on a temporary free Colab Tesla T4; NVIDIA partner-credit requests remain pending; no production GPU allocation; target binding MD not started.
**Scope:** Evidence review and readiness audit for the reported paracetamol-glucuronide (PaG) interaction with NNC2215. This is a computational research report, not a clinical model.

## Complete

- Reanalysed all 30 NNC2215 Figure 3d source-data cells against the publisher workbook. Paired fold changes across three independent experiments are 3.1576× (3 to 20 mmol/L glucose) and 12.4720× (0 to 20 mmol/L), matching the paper's rounded 3.2× and 12.5× values.
- Fit an exploratory shared-slope Hill model to the 15 NNC2215 source-workbook measurements. The best fit is EC50 13.43 mmol/L and apparent Hill slope 1.115. Resampling the three experiments 2,000 times produces a 2.5th EC50 percentile of 8.06 mmol/L; 64 resamples hit the 1,000 mmol/L fit ceiling, so the upper percentile is censored. The dataset does not identify a reliable midpoint or establish a 6 mmol/L threshold. The assumptions, residuals, raw source rows, fit record, and SVG are in `data/hill-fit-summary.json` and `assets/hill-fit.svg`.
- Recalculated the reported PaG/glucose binding-summary ratio (11.2× from rounded Kd means), the zero-glucose hIR-A affinity ratio at 70 μmol/L PaG (6.0× from rounded means), descriptive pig-clamp group-mean ratios, and the published Equation 4b coefficient component. Raw binding curves and individual animal observations are not available in the reanalysed inputs.
- Reproduced the publisher model-file hashes and checked insulin A/B sequence coverage, backbone atoms, three native disulfides, and modification-site LINK records in both public open/closed NNC2215 model PDBs.
- Completed a PaG-only chemistry and explicit-water setup preflight. PaG stereochemistry was retained and the isolated glucuronate parameterization has net charge −1. This does not validate any part of the NNC2215 conjugate.
- Completed a generic 1UBQ engine preflight on Apple M4 OpenCL/CPU. It is not an NVIDIA throughput benchmark and its performance is not transferable to the target system.
- Completed a synchronized, three-replicate, 100 ps generic 1UBQ throughput preflight using OpenMM 8.6.1 OpenCL on a temporary free Colab Tesla T4. Median throughput was 270.175 ns/day (257.885–282.691; 17,764 atoms). It is not target-system performance, a binding calculation, or an NVIDIA Inception/partner allocation. Full details are in `data/nvidia-engine-preflight.json`.

## Readiness findings and blockers

1. **The available NNC2215 PDBs are model-built, not experimentally resolved NNC2215 structures.** Their source sequences, backbone coverage, disulfide pairs, and documented modification-site records match the checks performed. The public files contain no insulin-receptor coordinates.
2. **The nonstandard conjugate topology is unresolved.** Duplicate atom names occur in the PHX, LYX, and MC nonstandard residues. The open-state model links MC atom C31 to LYX atom CD6, while the closed-state model links C31 to LYX atom C. This may be a file-labeling or chemistry issue; it must be reconciled against the source chemistry before creating a shared topology. Bond orders, stereochemistry of the complete conjugate, protonation, charges, and force-field parameters have not been validated.
3. **The original PaG MD cannot be exactly rerun from reviewed public files.** The 2026 Supporting Information reports Desmond 6.6, OPLS4, SPC water, 0.15 M NaCl, and 500 ns production using an in-house crystal structure of a related glucose-bound macrocycle. The reviewed article/SI and associated public model materials describe the structural input but do not supply those seed coordinates. The study's MD system is a related macrocycle analogue, not full NNC2215. Public 2024 NNC2215 state models are a different input and can support only a clearly labelled proxy study. The reviewed project files still lack the original seed coordinates. Patrik said he would contact the authors; the latest project status message asked for any update, but no coordinate file has been received.
4. **No target NNC2215–PaG binding simulation has run.** The authenticated NVIDIA portal lists pending $1,000 Lambda, $10,000 AWS, and $5,000 Azure credit requests; none is an award. The Lambda and AWS partner applications are unsubmitted, the Microsoft handoff requires the correct company identity, and the Innovation Lab application tile is awaiting NVIDIA support. The temporary free Colab T4 run above was a generic engine preflight only; no NVIDIA Inception or partner allocation and no project compute spend are recorded.

The record-level evidence is in [`nnc2215-model-validation.json`](data/nnc2215-model-validation.json); rerun it with [`validate_nnc2215_models.py`](validate_nnc2215_models.py). Source model hashes and provenance are in [`structure-manifest.json`](data/structure-manifest.json). The proposed pilot is in [`target-md-pilot-protocol.json`](data/target-md-pilot-protocol.json).

## Proposed next deliverable

Resolve the public-model nonstandard-residue atom/bond graph against the published chemistry and, where needed, the structure authors' source files; document an atom-mapped topology for each state; validate protonation, stereochemistry, partial charges, and force-field compatibility for the complete conjugate and PaG; and document how a PaG pose is initialized. If these checks pass, benchmark the selected engine on the allocated GPU using the prepared target system, then run the small three-condition stability pilot described in the protocol. If any topology or parameter check fails, stop before binding MD and report the unresolved evidence needed.

The two-week pilot/report estimate starts only when the necessary structural inputs, validated parameters, and compute are available. A later replicated/convergence study remains a separate four-to-eight-week planning estimate after the pilot is justified.

## Evidence required to pass the gates

| Gate | Evidence required | Current result |
|---|---|---|
| Public structure provenance | Source DOI/URL, licence, checksums, sequence/backbone/disulfide/link audit | Partial pass; model-built input, nonstandard-residue discrepancies open |
| Conjugate chemistry | Atom-mapped residue/linker graphs, bond orders, stereochemistry, protonation states, charge totals, and parameter validation | Blocked |
| Original-method reproducibility | Original macrocycle-analogue seed coordinates and compatible Desmond/OPLS4 setup, or explicitly labelled proxy inputs | Blocked for exact rerun |
| Pilot stability | Minimized valid system; finite energies; stable temperature/volume; no clashes/crashes; three independent seeds and consistent structural observables | Not run |
| GPU performance/cost | Target-system benchmark on allocated hardware with versions and measured ns/day | Generic T4 preflight complete; target-system performance/cost not measured; no production allocation |
| Binding inference | A validated free-energy/affinity method with convergence diagnostics and appropriate controls | Not part of the initial pilot |


## Update — 18 October 2026, 00:05 CEST

The Nature supplementary synthesis procedure now has a dedicated source audit at [`data/nnc2215-chemistry-source-audit.json`](data/nnc2215-chemistry-source-audit.json). It confirms the reported B29 alkyne-to-macrocycle click route, B1 glucoside coupling/deprotection, and final NNC2215 formula C325H461N79O97S6. It does not provide the model-specific atom mapping or force-field-ready topology. The model residues named PHX and LYX have coordinate-derived compositions that differ from the unrelated RCSB CCD entries using those same labels; these generic dictionary pages cannot be used to parameterize the models. The open and closed macrocycle-link atom labels still differ.

The authenticated NVIDIA portal was checked again at 00:05 CEST: Lambda ($1,000), AWS ($10,000), and Microsoft Azure ($5,000) remain catalog requests pending partner action. The Innovation Lab catalog search returns 0 results, and NVIDIA support has not replied. The official Lambda form remains incomplete. No application-specific values were entered manually and no submission was made. A redacted phone-field value appeared on page load, remained redacted after a clear attempt, and the review tab was closed without submission. Nebius was screened and is not eligible on the company’s recorded $0 raised status because its offer requires at least $5 million raised and VC/accelerator backing. Scaleway lists a 3–5 year operating-history condition and active EU go-to-market; eligibility is unverified, so no request was submitted.

No new GPU allocation or target MD run is recorded. The local dashboard source now includes the synthesis/component audit and timestamped compute-access state; these changes are pending commit and deployment.

## Update — 18 October 2026, 00:45 CEST

A synchronized generic OpenMM 8.6.1 OpenCL run completed on the free Colab Tesla T4. The solvated 1UBQ system contained 17,764 atoms; three 100 ps timed windows measured 282.691, 270.175 and 257.885 ns/day, with a 270.175 ns/day median. Each timed window included `getState(getEnergy=True)` to synchronize GPU execution. The earlier quick notebook measurement omitted that synchronization and is superseded. The saved notebook is [Dia Health One — PaG GPU Reproducibility Preflight](https://colab.research.google.com/drive/1ImHF0tk_wGSWfPwtuZjC04cFN7i4fCmt?authuser=1); the local run manifest and generated chart are linked from the dashboard.

This establishes only that a generic solvated protein can run on the temporary T4/OpenCL stack. It does not validate conjugate chemistry, estimate NNC2215–PaG throughput, reproduce the 2026 PaG MD, or predict affinity. No target system or paid compute was used. Public-model atom mapping and full conjugate parameters, original 2026 seed coordinates, and approved production compute remain gates.
