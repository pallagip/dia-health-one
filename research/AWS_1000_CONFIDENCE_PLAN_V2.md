# First-$1,000 AWS confidence plan — corrected version

**Updated:** 4 November 2026  
**Status:** executable gate plan; replaces the proposed Boltz-2 candidate-screening workflow  
**Budget source:** AWS credits  
**Objective:** buy calibrated evidence, not a large number of uncalibrated scores

## Executive decision

Do **not** spend the first $1,000 running a 50–60-candidate Boltz-2 screen or fitting Hill curves from repeated structure predictions.

Boltz-2's documented affinity module predicts the affinity of a small-molecule ligand to a **protein target**. It is not a validated affinity model for a synthetic macrocycle binding glucose or PaG. It also does not accept ligand concentration as a physical input, so repeating a prediction at labels such as 3, 5, 10 and 20 mM cannot generate a concentration-response curve or a Hill coefficient.

The valid computational route for this host–guest problem is benchmarked explicit-solvent molecular simulation with enhanced sampling, followed by prospective testing of a very small candidate set. The 2026 Chemical Science GluHUT study provides a public reference protocol and deposited data for that purpose.

## What the original proposal gets right

Keep these elements:

- strict provenance and immutable inputs;
- spot-instance cost controls and hard budget alarms;
- explicit interferent-first thinking;
- negative results retained;
- public evidence packages;
- no therapeutic efficacy claim from computation;
- no animal work before functional in-vitro evidence.

Correct these elements:

- replace Boltz-2 host–guest affinity scoring with a calibrated umbrella-sampling workflow;
- do not infer Hill coefficients from structure predictions;
- do not enumerate dozens of candidates before the parent method is validated;
- do not use short unrestrained MD as an affinity estimate;
- separate free-receptor binding from intramolecular switch behavior and insulin-receptor function;
- defer virtual-patient optimization until a measured or experimentally anchored switch curve exists.

## Scientific reference protocol

The public GluHUT free-energy study used:

- an open-cage truncated GluHUT model;
- GAFF bonded and nonbonded parameters for the host;
- RESP host charges based on HF/6-31G(d) electrostatic potentials;
- GLYCAM06 carbohydrate parameters;
- explicit TIP3P water;
- three independent unrestrained 100 ns runs per system;
- three steered unbinding paths per sugar;
- 51 umbrella windows per path;
- 0.3 Å window spacing;
- 50 kcal mol⁻¹ Å⁻² restraint;
- 10 ps restrained relaxation and 30 ns production per window;
- WHAM reconstruction and between-path uncertainty.

The published authors report qualitative agreement with the experimental sugar-affinity ranking, while also noting systematic free-energy overestimation and the limitations of a fixed-charge force field.

This protocol is a benchmark, not permission to assume that a new PaG or MC-Q1 result is correct.

## Program gates

### AWS-0 — public evidence and code freeze

**Compute spend:** $0  
**Release condition:** public confidence ledger and corrected plan are merged and accessible.

Actions:

1. Freeze exact input structures, protonation states, atom maps and hashes.
2. Publish the MC-Q1, MC-Q1D, MC-P3 and MC-P1 status with explicit evidence levels.
3. Record that MC-Q1 is a hypothesis and that the methoxy edit has a negative glucose-affinity prior from the public patent series.
4. Record the exact stop rules before running new calculations.

Deliverables:

- `research/CONFIDENCE_LEDGER.md`;
- this plan;
- public status page;
- immutable run manifest template.

### AWS-1 — reproduce deposited GluHUT analysis

**Maximum released spend:** $25  
**Purpose:** validate the analysis pipeline before generating trajectories.

Actions:

1. Obtain the authors' deposited simulation data from Zenodo DOI `10.5281/zenodo.19664426`.
2. Verify file hashes and licences.
3. Reconstruct at least one published PMF and its uncertainty from the deposited windows.
4. Reproduce the published qualitative glucose/galactose/fructose order and selected hydrogen-bond/water-occupancy summaries.
5. Record software versions, binning, equilibration discard and bootstrap/blocking choices.

Pass:

- regenerated PMF and summary statistics agree with the deposited/public record within a predefined numerical tolerance;
- analysis is deterministic from the archived inputs.

Fail:

- data are incomplete, analyses cannot be reproduced, or hidden manual choices materially affect the ranking.

No new GPU trajectories begin before AWS-1 passes.

### AWS-2 — exact chemistry and force-field gate

**Cumulative spend cap:** $100  
**Purpose:** ensure the simulated object is chemically real.

Actions:

1. Establish exact parent free-receptor graph and the relationship between the public truncated GluHUT model, the NNC2215 macrocycle reagent and the MC candidate graphs.
2. Generate and review host parameters using a source-matched RESP/GAFF workflow or a documented alternative.
3. Use GLYCAM06 for glucose and a reviewed parameterization for PaG.
4. Test at least two plausible PaG carboxylate/protonation and ion-treatment choices where chemically relevant.
5. Check energy minima, bond/angle distributions, charge totals, improper terms and structural stability.
6. Compare independent parameter routes where practical.

Pass:

- no unresolved topology or attachment ambiguity;
- stable short simulations;
- parameter alternatives do not create qualitatively incompatible host structures;
- every nonstandard parameter has provenance.

Fail:

- unresolved chemical graph, missing terms, implausible geometry, unstable host or materially contradictory parameterizations.

No candidate ranking occurs before AWS-2 passes.

### AWS-3 — benchmark trajectory reproduction

**Cumulative spend cap before review:** $350  
**Purpose:** demonstrate that our compute stack can reproduce a known host–guest result.

Run a staged benchmark:

1. **Throughput pilot:** one or two windows to measure ns/hour, storage and interruption behavior on the selected AWS instance.
2. **Single-path pilot:** one complete 51-window glucose PMF with the published spacing and production length.
3. **Full benchmark:** three independent paths only if the pilot has adequate overlap, stable free-energy plateaus and plausible bound-state structure.
4. **Hold-out sugar:** reproduce at least one lower-affinity sugar result without retuning the method.

Pre-register:

- window centers and force constants;
- equilibration discard;
- overlap threshold;
- convergence tests by trajectory fraction;
- WHAM/MBAR settings;
- bulk-region definition;
- uncertainty calculation;
- maximum acceptable deviation from the public benchmark.

Pass:

- correct qualitative sugar ranking;
- independent-path agreement within the predefined tolerance;
- stable result as later trajectory segments are added;
- adequate adjacent-window overlap;
- no hidden manual pose selection.

Fail:

- ranking changes with path, water model, reasonable protonation choice or analysis window;
- uncertainty is as large as the claimed difference;
- the benchmark cannot be reproduced without retrospective tuning.

### AWS-4 — source-matched parent glucose-versus-PaG calculation

**Additional spend release:** only after AWS-3 review  
**Purpose:** determine whether the calibrated workflow captures the known PaG liability in the parent system.

Important distinction:

- the experimental `Kd = 1.91 mM` for glucose and `0.17 mM` for PaG applies to the complete NNC2215 conjugate measured by native MS;
- the free-receptor calculation is not the same observable;
- the 2026 paper's exact in-house glucose-bound crystal structure was not made public.

Therefore the result must be labelled as a **source-matched free-receptor or model-system calculation**, not an NNC2215 affinity reproduction, unless the exact conjugate system and validated parameters are used.

Run:

- parent versus glucose;
- parent versus PaG;
- multiple unbinding paths;
- convergence and parameter sensitivity;
- explicit analysis of hydrophobic contacts, hydrogen bonds and solvent reorganisation.

Pass for ranking purposes:

- PaG-versus-glucose direction is robust across paths and reasonable parameter choices;
- uncertainty is materially smaller than the difference;
- structural interpretation is consistent with the published NMR observation of PaG aromatic shielding.

Fail:

- wrong or unstable affinity direction;
- result dominated by one path or one charge model;
- no plateau or insufficient overlap;
- parameter uncertainty exceeds the candidate effect being sought.

### AWS-5 — prospective four-compound comparison

**Release remaining credits only if AWS-4 passes.**

Compare, without protocol retuning:

1. MC-P0 parent;
2. MC-Q1 contact-site methoxy;
3. MC-Q1D distal-A methoxy control;
4. MC-P3 contact-site hydroxymethyl;
5. MC-P1 methyl size-down if the budget permits.

Pre-register the primary computational endpoint:

`DeltaDeltaG_selectivity = [G_bind(PaG) - G_bind(glucose)]_variant - [G_bind(PaG) - G_bind(glucose)]_parent`

A positive value in this convention favours improved glucose-over-PaG selectivity. Report the complete uncertainty and covariance treatment.

Advance a candidate computationally only if:

- the direction is consistent across independent paths;
- the confidence interval excludes no improvement;
- glucose binding is not predicted to deteriorate beyond the prespecified limit;
- MC-Q1 differs from MC-Q1D, supporting a local mechanism;
- the effect survives reasonable parameter alternatives.

This can move a molecule to evidence level E3. It cannot establish a hit or drug.

## Explicitly excluded from the first $1,000

Do not spend credits on:

- Boltz-2 affinity predictions for macrocycle–glucose or macrocycle–PaG host–guest binding;
- concentration-labelled Boltz runs intended to create a Hill curve;
- 50–60 unconstrained RDKit analogues;
- tandem-binder or multi-site insulin constructs before a single free binder passes;
- virtual-patient sweeps presented as molecular-discovery evidence;
- insulin–receptor complex MD before free-binder selectivity is established;
- MLIP or QM/MM refinement before the classical benchmark passes;
- large trajectory storage without a registered analysis question.

## Budget controls

1. Set AWS Budgets notifications at 10%, 25%, 50%, 75% and 90% of the credit allocation.
2. Use an automated kill condition at the agreed hard cap.
3. Record actual ns/hour and $/ns from the first throughput pilot; do not plan from generic vendor estimates.
4. Use checkpointed spot instances only after deterministic restart is verified.
5. Keep raw windows and metadata until the PMF is independently reconstructed; then apply lifecycle storage.
6. Do not release the next budget tranche merely because the previous compute completed. Release it only if the scientific gate passed.

## Confidence deliverables

Each gate produces:

- immutable input manifest;
- environment and image digest;
- exact command line;
- raw outputs and checksums;
- convergence and overlap plots;
- uncertainty calculation;
- pass/fail decision against preregistered criteria;
- public plain-language summary;
- negative results and deviations.

## Success criterion for the first $1,000

Success is **not** “50 candidates scored” or “five CIF files generated.”

The first $1,000 succeeds if it answers one of these questions defensibly:

1. Can our stack reproduce the public GluHUT host–guest benchmark?
2. Can a source-matched parent model reproduce a robust PaG-over-glucose liability?
3. Does any one of MC-Q1, MC-P3 or MC-P1 prospectively improve the computed selectivity relative to both parent and positional control?

A well-supported negative answer is useful and must stop the branch before synthesis. An uncalibrated positive score is not useful.

## External dependencies

Before AWS-4, request from the 2026 authors or platform owners:

- exact identity and coordinates of the in-house glucose-bound macrocycle analogue;
- force-field files or parameterization details used for glucose and PaG;
- initial poses and protonation states;
- any free-receptor PaG affinity measurement;
- clarification of permissible research use and collaboration/licensing path.

If these are unavailable, the computation must retain the label **public-model analogue**, and direct wet-lab binding becomes the decisive path.

## Source references

- Official Boltz documentation: affinity is supported for a small-molecule ligand against a protein target, with a recommended ligand-size range derived from training.
- Hoeg-Jensen et al., *Nature* 2024, DOI `10.1038/s41586-024-08042-3`.
- Sturis et al., *ACS Pharmacology & Translational Science* 2026, DOI `10.1021/acsptsci.6c00069`.
- Eades et al., *Chemical Science* 2026, DOI `10.1039/D6SC00738D`.
- Public GluHUT simulation data: Zenodo DOI `10.5281/zenodo.19664426`.

## Current execution decision

**AWS spend is paused at AWS-0/AWS-1 preparation.** Public evidence publication and deposited-data reproduction come first. No Boltz-2 candidate screen is authorized by this plan.
