# NVIDIA Innovation Lab application brief

**Company:** Dia Health One  
**Inception product:** Dia Health One Research Dashboard (Developing)  
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

NVIDIA's public Innovation Lab description lists dynamic simulations as a use case and describes GPU access as a separate, selection-based 60-day program for Inception members. The Inception profile and product are already set up; the Innovation Lab tile is currently absent from the accessible catalog, so this brief is prepared for the application and no application has been submitted.

- NVIDIA Innovation Lab: https://www.nvidia.com/en-us/data-center/innovation-lab/
- Inception portal: https://programs.nvidia.com/phoenix/benefits
- Research dashboard preview: https://feat-reproducibility-dashboa.dia-health-one.pages.dev/research/
- Pilot protocol: `data/target-md-pilot-protocol.json`
- Generic T4 benchmark: `data/nvidia-engine-preflight.json`
