# NVIDIA Innovation Lab application brief

**Archived — 25 October 2026.** NVIDIA Inception support confirmed on 19 October that Dia Health One's Inception membership is active and DGX Cloud Innovation Lab is no longer available. This brief is retained as a historical artifact and is superseded. Inception membership does not provide a GPU allocation or credit award.

**Company:** Dia Health One  
**Inception products:** Dia Health One Research Dashboard and Dia Health computational insulin discovery platform (both Developing / Not Accelerated)  
**NVIDIA account email:** andrii@diahealth.one  
**Request:** One CUDA-capable GPU environment for a 60-day molecular-dynamics research pilot and target-system benchmark.

## Project

Dia Health One is building a reproducible computational research workflow to examine published findings about glucose-responsive insulin NNC2215 and paracetamol glucuronide (PaG). The first objective is to reproduce the published experimental data and, only after chemistry and parameter checks pass, assess a clearly scoped molecular-dynamics proxy for the reported PaG interaction.

The work uses public structures and published data. It does not use patient records, make clinical predictions, recommend treatment, or provide insulin doses.

## GPU workload

- Prepare a CUDA environment for OpenMM and GROMACS and benchmark a chemically validated target system.
- If the structural and force-field gates pass, run the proposed pilot on apo, glucose-reference, and PaG conditions, with three independent 50 ns replicas per condition (450 ns total production planned).
- Record input checksums, model provenance, force-field versions, GPU/runtime details, seeds, diagnostics, and analysis outputs in versioned run manifests.
- Publish the non-clinical results and limitations through the Dia Health One reproducibility dashboard.

The initial request is for access to one suitable NVIDIA GPU, enough to benchmark a valid prepared system and begin the pilot. We are not quoting target GPU-hours yet: the measured Tesla T4 result is for a generic 17,764-atom ubiquitin system and cannot be transferred to NNC2215–PaG. After the target system is valid, we will measure its own throughput and update the compute estimate before scaling.

## Readiness and reproducibility

Completed: source-data calculations, structural inventories, a PaG-only solvent setup check, and a generic GPU engine preflight. On a temporary free Colab Tesla T4, OpenMM 8.6.1 OpenCL measured 270.175 ns/day median across three short 1UBQ windows. This is an engine check only and is not an NVIDIA allocation or target-system performance estimate.

Still required before target MD: reconcile the published open/closed model attachment records against authoritative conjugate chemistry; validate stereochemistry, protonation, charges, and force-field parameters; and obtain the original 2026 study seed coordinates or document a separate public-input proxy. The public 2026 Supporting Information record lists a 24-page PDF but no machine-readable coordinates or trajectory. Target binding MD has not started.

## Requested outcome

Access to the Innovation Lab's 60-day GPU environment and technical guidance for the CUDA molecular-dynamics workflow. A short target-specific benchmark will determine whether the proposed 450 ns pilot is feasible within that access window. The user-facing dashboard will distinguish published measurements, arithmetic reanalysis, model assumptions, simulation output, and unresolved questions.

## NVIDIA program fit and sources

NVIDIA's public Innovation Lab description lists dynamic simulations as a use case and describes GPU access as a separate, selection-based 60-day program for Inception members. The authenticated profile was rechecked on 18 October 2026 at approximately 09:59 UTC: two Developing / Not Accelerated product records are present, meeting NVIDIA’s published one-product prerequisite. The latest completed “Innovation Lab” catalog search returned 0 results at approximately 09:53 UTC; a direct page view at 10:00 UTC did not finish rendering a searchable result state. NVIDIA’s official page lists molecular dynamics, protein folding, and drug-discovery pipelines as simulation use cases and routes the application through the Innovation Lab tile in the Inception catalog. The tile and Request Benefit action remain unavailable, so no application has been submitted. A support request asking NVIDIA to enable the tile or provide the correct application route remained unanswered in the inbox check through approximately 10:01 UTC. The Benefits page still shows six pending partner requests and no GPU allocation. The €40,000 Scaleway request is still a catalog referral; its separate partner process remains incomplete.

Partner-offer screening: Nebius requires at least $5 million raised and VC/accelerator backing, which the current company profile does not meet. The Scaleway offer lists up to €42,000 but eligibility and current-customer status are unverified; its partner form requires verified identity, and no application was submitted.

- Historical NVIDIA Innovation Lab page: https://www.nvidia.com/en-us/data-center/innovation-lab/
- Inception portal: https://programs.nvidia.com/phoenix/benefits
- Research dashboard preview: https://feat-reproducibility-dashboa.dia-health-one.pages.dev/research/
- Pilot protocol: `data/target-md-pilot-protocol.json`
- Generic T4 benchmark: `data/nvidia-engine-preflight.json`
