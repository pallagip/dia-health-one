# Scientist handoff: inputs needed for an exact PaG simulation rerun

**Status:** prepared 18 October 2026; awaiting information from the study team. No message was sent to the authors in this update.

## What the published methods say

The 2026 study's molecular-dynamics system was a glucose-bound macrocycle analogue related to NNC2215, not the complete NNC2215 insulin conjugate. The supporting information reports a 500 ns glucose-bound trajectory, then says a representative trajectory snapshot was used in Maestro to build the macrocycle–paracetamol glucuronide (PaG) starting configuration. It reports Desmond 6.6 / Schrödinger 2021-2, OPLS4, SPC water and 0.15 M NaCl. The public Figshare record reviewed for this project lists the supporting-information PDF only; it does not list the in-house crystal coordinates, Desmond system, trajectory or extracted starting snapshot.

The source audit is recorded in [`data/pag-si-source-audit.json`](data/pag-si-source-audit.json) and [`data/pag-md-protocol.json`](data/pag-md-protocol.json). The [published supporting-information record](https://acs.figshare.com/articles/journal_contribution/Interaction_between_Paracetamol_Glucuronide_and_a_Glucose-Sensitive_Insulin_Molecule_NNC2215_Increased_Risk_of_Hypoglycemia/32292700) is the primary source.

## Request package

If the study team can share research inputs, the following would unlock an exact rerun most directly:

1. **Starting structure and system:** the in-house glucose-bound macrocycle-analogue crystal coordinates, chemical/component definitions, atom mapping, and the associated Desmond system/topology files (for example `.mae`, `.maegz`, or `.cms`). Please include any structure-preparation edits needed to get from the crystal model to the simulated system.
2. **Glucose-to-PaG handoff:** the exact glucose-trajectory frame used to seed the PaG setup (frame/time and coordinates), plus the Maestro project or build script that places PaG and constructs the PaG system. If available, include the glucose trajectory or a minimal extract containing the selected frame and its metadata.
3. **Ligand and parameter record:** the PaG molecular representation with stereochemistry and protonation state, atom names/mapping, charge/parameter files, and the parameterization method/version.
4. **Run configuration and software availability:** the Desmond input/configuration and complete equilibration and production settings, including box/solvent/ions, temperature, pressure, constraints, timestep, random seeds where recorded, and software/force-field versions. The reported stack is Desmond 6.6 / Schrödinger 2021-2 / OPLS4. Please confirm whether an executable system can be shared and whether matching software/force-field access is available. The SI gives the headline method but does not enumerate every equilibration setting.
5. **Validation outputs (if shareable):** run logs, energy/temperature diagnostics and the PaG trajectory or selected analysis frames, so a rerun can be compared against more than the reported summary.

If the full project is too large, the minimum useful first transfer is (a) the original coordinate and topology/system files, (b) the exact PaG starting frame or build project, and (c) the PaG parameter files plus the complete run settings. A checksum and a short note identifying which published figure/panel each system corresponds to would make the handoff easier to verify.

If matching Desmond/OPLS4 access is unavailable, an OpenMM or GROMACS implementation would be a port or reconstruction that needs its own parameter and observable validation; it should not be described as an exact engine-level rerun.

## Reproduction boundary

- With the original system and setup files, reproduce the published **macrocycle-analogue** simulation first and compare structural/trajectory observables.
- A study using the public 2024 NNC2215 open/closed models would be a separate **proxy study**, not an exact rerun of the 2026 PaG simulation. Those model-built files also retain open chemistry and parameterization checks in [`data/target-md-pilot-protocol.json`](data/target-md-pilot-protocol.json).
- Do not infer missing coordinates, atom mapping, charge states, parameters or a binding pose. Keep target binding MD gated until the received files and chemistry are validated.
- No simulation result, affinity estimate, clinical prediction or dosing implication is implied by this checklist.

## Current response status

The latest inbox review on 18 October 2026 at 10:49 UTC found no new study-author reply or source files. Patrik's latest inbound note remains 17 October 2026 at 21:01 UTC and says he will contact the authors. The project is prepared to inspect and checksum any files when they arrive.
