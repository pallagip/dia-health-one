# NNC2215–PaG reproducibility milestone report

**Date:** 2026-10-17 21:40 UTC
**Status:** Published-data reanalysis and Hill-fit identifiability check complete; three NVIDIA catalog credit requests pending; no GPU allocated; target binding MD not started.
**Scope:** Evidence review and readiness audit for the reported paracetamol-glucuronide (PaG) interaction with NNC2215. This is a computational research report, not a clinical model.

## Complete

- Reanalysed all 30 NNC2215 Figure 3d source-data cells against the publisher workbook. Paired fold changes across three independent experiments are 3.1576× (3 to 20 mmol/L glucose) and 12.4720× (0 to 20 mmol/L), matching the paper's rounded 3.2× and 12.5× values.
- Fit an exploratory shared-slope Hill model to the 15 NNC2215 source-workbook measurements. The best fit is EC50 13.43 mmol/L and apparent Hill slope 1.115. Resampling the three experiments 2,000 times produces a 2.5th EC50 percentile of 8.06 mmol/L; 64 resamples hit the 1,000 mmol/L fit ceiling, so the upper percentile is censored. The dataset does not identify a reliable midpoint or establish a 6 mmol/L threshold. The assumptions, residuals, raw source rows, fit record, and SVG are in `data/hill-fit-summary.json` and `assets/hill-fit.svg`.
- Recalculated the reported PaG/glucose binding-summary ratio (11.2× from rounded Kd means), the zero-glucose hIR-A affinity ratio at 70 μmol/L PaG (6.0× from rounded means), descriptive pig-clamp group-mean ratios, and the published Equation 4b coefficient component. Raw binding curves and individual animal observations are not available in the reanalysed inputs.
- Reproduced the publisher model-file hashes and checked insulin A/B sequence coverage, backbone atoms, three native disulfides, and modification-site LINK records in both public open/closed NNC2215 model PDBs.
- Completed a PaG-only chemistry and explicit-water setup preflight. PaG stereochemistry was retained and the isolated glucuronate parameterization has net charge −1. This does not validate any part of the NNC2215 conjugate.
- Completed a generic 1UBQ engine preflight on Apple M4 OpenCL/CPU. It is not an NVIDIA throughput benchmark and its performance is not transferable to the target system.

## Readiness findings and blockers

1. **The available NNC2215 PDBs are model-built, not experimentally resolved NNC2215 structures.** Their source sequences, backbone coverage, disulfide pairs, and documented modification-site records match the checks performed. The public files contain no insulin-receptor coordinates.
2. **The nonstandard conjugate topology is unresolved.** Duplicate atom names occur in the PHX, LYX, and MC nonstandard residues. The open-state model links MC atom C31 to LYX atom CD6, while the closed-state model links C31 to LYX atom C. This may be a file-labeling or chemistry issue; it must be reconciled against the source chemistry before creating a shared topology. Bond orders, stereochemistry of the complete conjugate, protonation, charges, and force-field parameters have not been validated.
3. **The original PaG MD cannot be exactly rerun from reviewed public files.** The 2026 Supporting Information reports Desmond 6.6, OPLS4, SPC water, 0.15 M NaCl, and 500 ns production using an in-house crystal structure of a related glucose-bound macrocycle. The reviewed article/SI and associated public model materials describe the structural input but do not supply those seed coordinates. The study's MD system is a related macrocycle analogue, not full NNC2215. Public 2024 NNC2215 state models are a different input and can support only a clearly labelled proxy study. No inquiry has yet been sent to the study authors; requesting the original structure files is a proposed next step.
4. **No target NNC2215–PaG binding simulation has run.** The authenticated NVIDIA portal lists pending $1,000 Lambda, $10,000 AWS, and $5,000 Azure credit requests; none is an award. The Lambda and AWS partner applications are unsubmitted, the Microsoft handoff requires the correct company identity, and the Innovation Lab application tile is awaiting NVIDIA support. No GPU is allocated and no compute spend has occurred.

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
| GPU performance/cost | Target-system benchmark on allocated hardware with versions and measured ns/day | Not run; no GPU allocated |
| Binding inference | A validated free-energy/affinity method with convergence diagnostics and appropriate controls | Not part of the initial pilot |
