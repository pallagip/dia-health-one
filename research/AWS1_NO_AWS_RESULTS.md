# AWS-1 no-cloud calibration results

**Updated:** 4 November 2026  
**AWS credits spent:** `$0`  
**Current stage:** public benchmark calibration and chemistry-provenance audit  
**Candidate ranking:** not authorized

## Summary

The period before AWS access was used to reproduce and stress-test the public GluHUT benchmark with GitHub Actions and independent analysis code.

The result is useful but qualified:

- the deposited WHAM workflow reproduces the published sugar-group means and qualitative sugar order;
- an independent NumPy WHAM implementation reproduces the three glucose routes within `0.000454 kcal/mol`;
- one fructofuranose route fails a prespecified short-window sensitivity test;
- all four deposited sugar systems use one exactly identical compiled 159-atom host;
- the host formula is `C60H75N15O9`, with six urea and three amide carbonyls;
- the two-atom `ROH` residue is covalently attached to the carbohydrate and completes a neutral `C6H12O6` guest; it is not a counterion;
- the source RESP/GAFF/GLYCAM generation package and standalone start coordinates are not deposited;
- no PaG or MC-Q1 ranking is yet justified.

## Public WHAM reproduction

The checksum-verified `WHAM_PMF.zip` archive was analysed using the deposited Grossfield WHAM workflow.

The reproduced sugar-group means agree with the public values within one published standard deviation. The qualitative order was reproduced without sugar-specific retuning.

A separately implemented NumPy WHAM solver was then applied to the three glucose routes. Its maximum absolute difference from the deposited workflow was:

`0.000454 kcal/mol`

This establishes that the glucose result is not dependent on the downloaded WHAM executable alone.

## Fructofuranose US-2 sensitivity failure

The deposited fructofuranose US-2 12.6 Å window contains only 205 frames, versus approximately 15,001 in neighbouring windows.

| Analysis | Binding free energy | Change from baseline |
|---|---:|---:|
| Deposited baseline | −7.521 kcal/mol | 0.000 |
| Exclude 12.6 Å window | −6.815 kcal/mol | **+0.707 kcal/mol** |

The preregistered sensitivity threshold was `0.50 kcal/mol`, so this gate failed.

The short window bridges an otherwise nearly disconnected interval:

- overlap 12.3 ↔ 12.6 Å: `0.0944`;
- overlap 12.6 ↔ 12.9 Å: `0.1371`;
- direct overlap 12.3 ↔ 12.9 Å after exclusion: `0.000133`.

The overall sugar benchmark remains useful, but fructofuranose US-2 must carry this omission sensitivity or be resampled around 12.3–12.9 Å.

## Common compiled host

The public `SMD_trajectories.zip` archive was checksum and ZIP-CRC verified and its four AMBER topologies were parsed atom by atom.

Each topology contains:

- 183 atoms;
- a 159-atom host residue labelled `MOL`;
- a 22-atom carbohydrate residue;
- a two-atom `ROH` carbohydrate-completion residue.

The compiled host is exactly identical across glucose, galactose, fructopyranose and fructofuranose:

- maximum per-atom charge difference: `0.0 e`;
- maximum mass difference: `0.0`;
- identical atom names and AMBER types;
- identical bonds, angles, dihedrals and scaling factors;
- identical exclusions;
- identical physical Lennard–Jones coefficients;
- identical compiled indexing.

Physical-host fingerprint:

`78e512d17cf1a9db672f22da49a3e59f96bb0af7dae5c059e8a4d1353de46cb1`

This supports the common-host assumption for the public sugar benchmark.

## Compiled graph chemistry

The host is one connected component with:

- formula `C60H75N15O9`;
- 84 heavy atoms;
- effectively neutral compiled charge;
- six urea carbonyls;
- three amide carbonyls;
- nine methyl groups;
- heavy-graph cycle rank 7.

The host formula and motif inventory are identical across all four systems. These values are suitable for atom-level reconciliation with the published truncated GluHUT chemical drawing. They do not reveal how the RESP charges were generated.

### Correction: `ROH` is covalent guest chemistry

The earlier residue-level charge partition could be read as a positively charged sugar balanced by a negative `ROH` counterion. Direct graph inspection disproves that interpretation.

In every system, carbohydrate atom 160 is covalently bonded to oxygen atom 183 in `ROH`:

- glucose and galactose: `C1–O1(ROH)`;
- fructopyranose and fructofuranose: `C2–O1(ROH)`;
- compiled equilibrium bond length: 1.43 Å;
- compiled force constant: 320 kcal mol⁻¹ Å⁻².

Together the two residue records form one connected neutral 24-atom guest with formula `C6H12O6`. The approximately `+0.194 e` and `−0.194 e` residue charges are an internal template partition of one molecule, not an ion pair.

This removes a false chemistry blocker. The derivation and transferability of the guest charges remain unvalidated.

## Corrected SMD archive inventory

The archive contains:

- four compiled `parm7` topologies;
- twenty `binpos` trajectories;
- three AMBER NMR-style restraint-control files whose filenames end in `.rst`;
- two protocol/documentation files;
- no standalone coordinate restart file;
- no source RESP/ESP, `mol2`, `frcmod`, `lib` or equivalent parameter-generation package.

The `.rst` files contain `&rst`, `iat`, `igr1`, `igr2`, `r2`, `rk2` and related restraint fields. They are not restart coordinates.

## Outstanding chemistry questions

The compiled host and guests can be used to reproduce the public benchmark, but the following remain unavailable:

- host RESP inputs, outputs and ESP calculation;
- source GAFF atom typing and parameter files;
- exact GLYCAM residue/template setup;
- tleap or equivalent assembly scripts;
- standalone starting coordinates;
- a source-matched PaG model;
- MC-Q1-specific parameters.

The compiled `C60H75N15O9` host must still be reconciled atom by atom with the published truncated chemical graph and synthesis.

## External requests sent

Requests were sent to:

- Adrian Mulholland and Anthony Davis for source RESP/GAFF/GLYCAM files, starting coordinates, the carbohydrate residue split and clarification of the short fructofuranose window;
- Jeppe Sturis for the exact receptor analogue, coordinates and parameters used in the NNC2215–PaG study;
- Carbometrics for platform access, rights and four-compound feasibility;
- SpiroChem for a site-defined synthesis and positional-assignment quotation;
- Domainex for orthogonal direct glucose-versus-PaG binding-method feasibility.

A conditional reply watch is active for substantive author and provider responses.

## Current decision

The public host–guest analysis pipeline has crossed an important calibration threshold, but not the threshold required for candidate ranking.

Authorized next work:

1. map the SMD restraint atom groups to exact host and guest atoms;
2. reconcile the compiled `C60H75N15O9` host with the published truncated chemical graph;
3. obtain source parameter-generation and starting-coordinate provenance;
4. prepare a source-matched parent glucose-versus-PaG calibration only after those checks pass;
5. obtain synthesis and direct-binding quotations in parallel.

Not authorized:

- transferring public host charges directly to MC-Q1;
- inventing PaG parameters without review;
- ranking MC-Q1, MC-Q1D or MC-P3 from the current benchmark;
- presenting a computed affinity as a measured hit;
- claiming a working drug.
