# Compute and molecular-simulation readiness

**Checked:** 2026-10-17 UTC
**Decision:** source-data reanalysis proceeds on CPU; atomistic GPU work remains gated.

## Work completed

The source-data reanalysis ran with Python 3.14.4 and the Python standard library on a local macOS arm64 host. It processes 30 numeric source values and produces the replicate summary, paired fold-change output, input checksum, and original SVG chart. This analysis is small and CPU-suitable; using a GPU for these calculations would not improve the scientific result.

## NVIDIA/GPU environment check

- The local host has no NVIDIA GPU runtime available.
- A read-only Google Cloud Compute Engine query was attempted against a Dia candidate project. Google Cloud reported that `compute.googleapis.com` is disabled for that project, so regional accelerator quota could not be inspected.
- The API state and quota for other Dia projects are not established by that single query. NVIDIA Inception membership on its own does not confirm a GPU allocation.
- NVIDIA has approved Dia Health One for Inception. The welcome message points members to request benefits through the portal; it does not confirm that cloud credits or a GPU instance have already been allocated.
- No cloud instance was created, and no GPU compute cost was incurred.
- OpenMM/GROMACS GPU benchmarking is therefore **not reported as complete**. It requires a confirmed NVIDIA environment and a comparable, documented test system.

## Independent structural gate

The 2026 PaG paper describes a 500 ns Desmond/OPLS4 simulation initialized from an in-house crystal structure of a related macrocycle, but the article does not provide those starting coordinates. The NNC2215 conjugate also requires verified non-standard chemical parameters. Running MD from invented coordinates or unvalidated parameters would not reproduce the published calculation, so it is not included here.

## Next evidence-based gate

Before benchmarking engines or launching a simulation, establish all of the following:

1. Which Dia project is intended for research compute and whether the Inception/partner credits are actually allocated.
2. The enabled compute service, GPU quota, region, GPU model, and permitted runtime.
3. A spend cap or confirmed no-cost credit for any cloud instance.
4. Access to the published starting structure or a new, explicitly labelled model with a validated parameterisation plan.
5. The same documented benchmark system and protocol for OpenMM and GROMACS.
