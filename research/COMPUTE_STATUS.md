# Compute and molecular-simulation readiness

**Checked:** 2026-10-17 UTC

**Decision:** published-data reanalysis is complete on CPU; the 2026 PaG molecular-dynamics method is now documented, but its starting coordinates and an NVIDIA allocation are still unavailable. No target MD run has been made.

## Reproductions completed

- Nature 2024 Figure 3d source values were extracted with workbook-cell provenance and verified against the publisher workbook. The paired fold changes reproduce the reported rounded 3→20 and 0→20 mmol/L benchmarks.
- The 2026 PaG pig-clamp group means, SEM, sample sizes, and reported p-values are transcribed from Table 1. The CPU script derives treatment-to-saline ratios from group means; individual-animal records are not available here, so the published SEM and statistical tests are not recomputed.
- The PaG paper's Appendix S4 Equation 4b and Table S2 coefficient are reapplied at the published human and pig PaG concentration anchors. This reproduces the model's linear insulin-sensitivity component only, not GIR, receptor affinity, or patient outcomes.

## Structure data recovered

Springer Nature Figshare publishes author-supplied open and closed NNC2215 model PDBs under CC BY 4.0. The files are included in `structures/`; provenance and source/local checksums are recorded in `data/structure-manifest.json`.

| State model | Source PDB | Local SHA-256 |
|---|---|---|
| Active/open, free D-glucose in the B29 macrocycle | `NNC2215_active_open form.pdb` | `44c29d45b28ace4ede73a75a8e076dbca62e14b34bd3938a4b495ebd1270267a` |
| Inactive/closed, tethered glucoside in the B29 macrocycle | `NNC2215_inactive_closed form.pdb` | `a8287b859e73b274e8403ee56cd83564bc30fd0503d852829fa805fdc777e12d` |

These files are model-built coordinates, not experimentally resolved NNC2215 structures. They are not the in-house crystal structure of a related macrocycle used to initialize the 2026 PaG simulation. The files also contain non-standard chemistry that still needs validated force-field treatment before any energetic interpretation.

## Published PaG MD protocol and remaining gap

The 2026 Supporting Information is available as a separate ACS Figshare record (SI DOI `10.1021/acsptsci.6c00069.s001`, CC BY-NC 4.0). It reports:

- Desmond v6.6, Schrödinger Release 2021-2
- OPLS4 force field and SPC explicit water
- 0.15 mol/L NaCl
- 500 ns production phase
- An in-house crystal structure of a glucose-bound macrocycle similar to NNC2215 as the initial system
- A representative glucose-bound snapshot used to build the PaG configuration in Maestro, followed by the same MD protocol

The Supporting Information does not provide the starting crystal coordinates; it describes the structure as in-house and “data not shown.” This method and the public-input gap are recorded separately in `data/pag-md-protocol.json`. The local NNC2215 models can support a separately labelled proxy study, but they cannot be substituted silently and called an exact rerun of the published PaG simulation.

## Compute access

- Dia Health One's NVIDIA Inception membership is approved. Membership itself does not allocate compute.
- The separate DGX Cloud Innovation Lab request has not been submitted. The portal URL opened in Chrome on 2026-10-17 returned `net::ERR_CERT_AUTHORITY_INVALID`. No warning was bypassed and no credentials were entered on that page. A GPU instance/job has not been created and no GPU spend has occurred.
- The host is macOS arm64 without an NVIDIA GPU runtime. OpenMM/GROMACS GPU benchmarking is not complete.
- The Google Cloud account remains a Free Trial; its billing console did not allow GPU attachment/quota requests. The Google for Startups Cloud Program application was received, but no approval was visible in the last account review.

## Next gates

1. Access the NVIDIA Inception portal through a trusted browser/network, then submit the separate Innovation Lab benefit request and await selection.
2. Obtain the original macrocycle coordinates from the authors or another attributable source; otherwise define a proxy study using the published NNC2215 state models and label it separately.
3. Validate parameters for PaG and all non-standard residues/linkages in the chosen model. Do not interpret energies before this gate passes.
4. Once chemistry and compute are ready, benchmark OpenMM and GROMACS on the same documented test system, then run a small pilot with stability diagnostics.
5. Compare molecular outputs only with observables they actually predict. The published receptor-binding assay and animal clamp results remain experimental benchmarks, not outputs of the existing CPU calculations.
