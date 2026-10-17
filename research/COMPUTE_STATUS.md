# Compute and molecular-simulation readiness

**Checked:** 2026-10-17 UTC

**Decision:** source-data reanalysis is complete on CPU; atomistic GPU work remains gated.

## Reproduction completed

The Nature 2024 Figure 3d source values were extracted with workbook-cell provenance and independently verified against the publisher workbook. The Python standard-library analysis reproduces the reported paired glucose-response fold changes from three experiments. This calculation is CPU-suitable; a GPU would not improve it.

The 2026 PaG pig-clamp AUC means, SEM values, sample sizes, and reported p-values are transcribed from the main article's Table 1. A second CPU script derives treatment-to-saline ratios from those published group means. These are descriptive ratios only; individual-animal records are not available here, so the published SEM and statistical tests cannot be independently recomputed. The supporting-information Table S2 model parameters are catalogued in `data/pag-model-parameters.csv`, but they are insufficient to reconstruct the full PK/PD fit and input profiles.

A third CPU calculation reapplies Appendix S4 Equation 4b (`SI = SIbase * (1 + ConcPaG * CoeffPaG / 100)`) at the paper's 70 μmol/L human estimate and 102/426 μmol/L pig concentration anchors. It reproduces the fitted model's linear insulin-sensitivity component only; it does not predict GIR or patient outcomes.

## Compute access

- Dia Health One's NVIDIA Inception membership is approved. NVIDIA's welcome email says members can request partner benefits through the portal; membership alone does not allocate compute.
- NVIDIA's [DGX Cloud Innovation Lab](https://www.nvidia.com/en-eu/data-center/dgx-cloud/innovation-lab/) offers selected Inception members a 60-day cloud environment. A separate application is required, approval depends on fit and availability, and members do not qualify automatically. The Innovation Lab request has not been submitted.
- The NVIDIA portal currently requires account sign-in before the benefit request can be checked or submitted.
- Google's billing console confirms an active Free Trial account. Google's [Free Trial documentation](https://docs.cloud.google.com/free/docs/free-cloud-features) says GPU attachment and quota requests are unavailable while billing remains a non-billable Free Trial account. A [Google for Startups Cloud Program application](https://cloud.google.com/startup/apply) has been received; no approval is visible yet.
- No GPU instance or simulation job has been created, and no GPU spend has been incurred. No paid billing upgrade has been made.
- The local host is macOS arm64 without an NVIDIA GPU runtime. OpenMM/GROMACS GPU benchmarking is not complete.

## Independent structural gate

The 2026 PaG paper reports a 500 ns Desmond/OPLS4 simulation initialized from an in-house crystal structure of a related macrocycle. The article does not publish those starting coordinates. The NNC2215 conjugate also needs verified parameters for its non-standard switch chemistry. Therefore this project does not claim an exact MD rerun or run target-system MD from invented coordinates or unvalidated parameters.

The public article and its supplement are listed on [PubMed Central](https://pmc.ncbi.nlm.nih.gov/articles/PMC13270365/). The supplement describes the molecular-dynamics and NMR methods and publishes selected PK/PD model parameters; these values are catalogued in the research data folder. No macrocycle coordinates, target conjugate topology, or atomistic trajectory is supplied.

## Next gates

1. Sign in to the NVIDIA Inception portal and request the Innovation Lab benefit; wait for selection and access details.
2. If Google for Startups approves the application, confirm the resulting billing and GPU-quota terms before using that account for GPU work.
3. Obtain the published starting coordinates or define a separately labelled model, then validate the conjugate and PaG parameterisation.
4. Once compute and chemistry are ready, benchmark OpenMM and GROMACS on the same documented test system and protocol.
5. Run target MD only after the structure, parameters, diagnostics, and comparison observable are specified.

The source-data reanalysis, published measurements, model assumptions, and any future simulation outputs must remain separate in the dashboard and run manifests.
