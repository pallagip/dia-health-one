# GluHUT SMD restraint-coordinate audit

**Date:** 4 November 2026  
**Status:** host reference mapping passed; guest-coordinate consistency gate failed  
**Scope:** deposited public GluHUT steered-MD restraint inputs only

## Question

The public GluHUT study describes a pulling coordinate based on centres of mass of selected host and sugar atoms. Before regenerating steered-MD or umbrella paths, the exact AMBER `igr1` and `igr2` atom selections were mapped onto the checksum-verified compiled topologies.

The purpose was to determine whether all four sugar systems use the same chemical definition of the reaction coordinate.

## Host reference group

The host selection is consistent across glucose, galactose, fructopyranose and fructofuranose.

It contains the same twelve `MOL` carbons in every topology:

`C35, C36, C39, C42, C45, C58, C3, C4, C5, C8, C19, C22`

The induced graph is exactly two disconnected six-member aromatic rings:

- 12 selected atoms;
- all carbon;
- all AMBER type `ca`;
- two connected components of six atoms;
- cycle rank 2.

This part of the deposited coordinate definition is internally consistent.

## Guest groups are not chemically uniform

The deposited guest selections differ by sugar.

| System | Selected chemical group |
|---|---|
| glucose | six-member pyranose ring plus exocyclic C6 |
| galactose | six-member pyranose ring plus exocyclic C6 |
| fructopyranose | six-member ring heavy atoms only |
| fructofuranose | five-member ring plus isolated hydrogen `H62` |

### Glucose and galactose

Both use seven atoms:

`C1, C2, C3, C4, C5, O5, C6`

The six ring atoms form the cyclic core. C6 is an exocyclic leaf in the selected subgraph.

### Fructopyranose

The selection contains six connected ring heavy atoms:

`C2, O6, C6, C5, C4, C3`

This is a chemically coherent six-member cyclic group.

### Fructofuranose

The selection contains:

`C2, O5, C5, H62, C4, C3`

The five heavy atoms `C2, O5, C5, C4, C3` form the furanose ring. `H62` is a hydrogen atom and is isolated from the other selected atoms in the induced group graph.

This may be intentional, a naming/indexing artifact, or a deposited input error. The archive alone cannot distinguish these possibilities.

## Restraint parameters

The mapped files use a group-centre definition with:

- `iat = -1, -1`;
- initial target `r2 = 0.0`;
- final target `r2a = 15.0`;
- `rk2 = 200.0`;
- `iresid = 0`.

The recorded parameter dictionaries are not completely identical because one input includes an `outxyz` setting not present in the others. The physical significance of that difference for the deposited production workflow has not been established.

The relationship between `rk2 = 200` in the AMBER group-restraint file and the paper’s reported restraint constant remains to be confirmed from the exact AMBER convention and author response.

## Decision

The consistency gate fails for prospective trajectory generation:

- host reference definition: **pass**;
- guest-coordinate chemical uniformity: **fail**;
- fructofuranose all-heavy-atom selection: **fail**;
- force-constant convention: **unresolved**.

Therefore:

1. reanalysis of the already deposited reaction-coordinate time series and WHAM inputs remains reproducible;
2. new SMD or umbrella trajectories must not be generated from the deposited restraint files as though they defined one uniform cross-sugar coordinate;
3. parent glucose-versus-PaG calibration is not authorized from this input set;
4. MC-Q1 or other candidate ranking remains unauthorized;
5. author clarification or a corrected fructofuranose restraint is required.

## Relationship to the fructofuranose sensitivity result

The unexplained `H62` selection occurs in the same sugar system where fructofuranose US-2 contains a 205-frame short window and fails the preregistered omission-sensitivity gate by 0.707 kcal/mol.

This is a co-occurrence, not a demonstrated causal link. The short window belongs to an umbrella route, while the current result concerns the deposited SMD group definition used to generate pulling paths. No causal claim is made without coordinates, trajectory reconstruction and author clarification.

## External clarification requested

A follow-up was sent to Adrian Mulholland and Anthony Davis asking:

- whether `H62` was intentionally included;
- which atom was intended if it is an index error;
- why glucose and galactose include exocyclic C6 while fructopyranose uses ring heavy atoms only;
- how `rk2 = 200` maps to the published restraint constant;
- whether corrected restraint inputs or exact generation scripts are available.

## Reproducibility

- Atom mapper: `research/audit_gluhut_smd_restraints.py`.
- Consistency summarizer: `research/summarize_gluhut_smd_restraint_anomaly.py`.
- Detailed mapping: `research/data/gluhut-smd-restraint-mapping.json`.
- Bounded result: `research/data/gluhut-smd-restraint-consistency.json`.
- Tests: `tests/test_audit_gluhut_smd_restraints.py` and `tests/test_summarize_gluhut_smd_restraint_anomaly.py`.
- Workflow: `.github/workflows/gluhut-smd-restraint-mapping.yml`.

## Evidence boundary

This is a chemical audit of deposited reaction-coordinate inputs. It is not a binding free-energy result, NNC2215 result, PaG result, MC-Q1 result, efficacy, safety or working-drug result.
