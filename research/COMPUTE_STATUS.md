# Compute and molecular-simulation readiness

**Checked:** 2026-10-17 UTC

**Decision:** the published-data analyses run reproducibly on the local CPU. The NVIDIA Inception membership is approved, but the separate Innovation Lab request has not been filed and no GPU is allocated. The published PaG MD method is documented, but its in-house starting coordinates are not public and force-field readiness for the NNC2215 conjugate remains unresolved. No target molecular-dynamics run has been made.

## Reproductions completed

Three Python analyses completed on the Apple M4 host:

1. **Nature 2024 Figure 3d source values:** 30 CSV values were verified against the publisher workbook's source cells and workbook SHA-256. Paired fold changes across three independent experiments are 3.1576× for 3→20 mmol/L glucose and 12.4720× for 0→20 mmol/L, matching the paper's rounded 3.2× and 12.5× values.
2. **2026 PaG pig clamp, Table 1:** the published group means, SEM, sample sizes, and p-values were transcribed. The script calculates descriptive ratios of group means: 1.1466× for low-dose PaG and 1.4226× for high-dose PaG versus saline. Individual animal values are not available in this data set, so SEM and significance tests are not recomputed.
3. **2026 PaG model, Appendix S4 Equation 4b:** the reported `CoeffPaG` value is reapplied at published human and pig concentration anchors. This reproduces only the linear insulin-sensitivity component, not receptor affinity, glucose-infusion rate, or patient response.

The run records and input checksums are in `data/reproduction.json`, `data/pag-clamp-reproduction.json`, and `data/pag-sensitivity-reproduction.json`. `data/nvidia-access-status.json` records the compute-access state and source links.

## Structure data and PaG MD method

Springer Nature Figshare publishes author-supplied open and closed NNC2215 model PDBs under CC BY 4.0. The files are included in `structures/`; source and local checksums are recorded in `data/structure-manifest.json`. These are model-built coordinates, not experimentally resolved NNC2215 structures.

`inspect_structure_models.py` verifies those checksums and records PDB connectivity metadata in `data/structure-preflight.json`: the open model has 990 coordinate records and the closed model 968, each with three SSBOND and three LINK records. This only inventories file records; it does not validate chemistry or parameters.

The 2026 Supporting Information reports Desmond 6.6 (Schrödinger Release 2021-2), OPLS4, SPC explicit water, 0.15 mol/L NaCl, and 500 ns production. The input is an in-house crystal structure of a related glucose-bound macrocycle. The coordinates are described as in-house and are not supplied. The public NNC2215 state models are different inputs and cannot be substituted silently for an exact rerun. The recipe is recorded in `data/pag-md-protocol.json`.

## NVIDIA compute access

- Dia Health One's NVIDIA Inception membership was approved on 2026-10-15. Membership does not allocate a GPU.
- NVIDIA's [Innovation Lab](https://www.nvidia.com/en-us/data-center/innovation-lab/) is a separate 60-day program for selected Inception members. A request must be submitted through the Inception Portal; NVIDIA reviews applications on a rolling basis and states that approval depends on availability and fit.
- The official portal opens to NVIDIA sign-in. The user must complete account authentication before the benefit catalog and request status can be inspected. No Innovation Lab request has been submitted.
- The current host is macOS 26.3 on Apple M4/arm64. `nvidia-smi`, OpenMM, GROMACS, and Desmond are not installed; no NVIDIA instance/job has been started and no compute spend has occurred.

The official request path is: **Inception Portal → Profile → Add product → Benefits → Explore Benefits → Innovation Lab → Request Benefit.** Keep the final application unsubmitted until the user reviews its contents.

## Gates before any target MD run

1. Obtain Innovation Lab access or another explicitly approved compute allocation.
2. Obtain the original PaG MD starting coordinates from the authors or another attributable source. If they remain unavailable, define any public-input proxy as a separate study.
3. Validate charges and force-field parameters for PaG, the macrocycle, glucoside, linker, and insulin conjugate before interpreting energetics.
4. Benchmark OpenMM and GROMACS on the same documented test system in the actual NVIDIA environment.
5. Run a small stability pilot with documented diagnostics and independent replicates before considering longer production simulations.

Experimental receptor-binding and pig-clamp results remain literature benchmarks. The existing CPU calculations do not predict or reproduce those assays through molecular dynamics.
