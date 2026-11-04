# GluHUT truncated-host architecture reconciliation

**Date:** 4 November 2026  
**Status:** compiled molecular graph and formula reconciled to the published methyl-truncated host  
**Scope:** public GluHUT benchmark host only

## Result

All four deposited AMBER systems encode the same molecular graph described by the published GluHUT architecture after truncation at the outer amide linker and replacement of each removed dendrimer by methyl.

The compiled host formula is:

`C60H75N15O9`

The primary-source architecture accounts for that formula exactly:

| Module | Multiplicity | Formula contribution |
|---|---:|---:|
| triethyl/tris(methylene) aromatic cap with three roof-side NH atoms | 2 | `2 × C15H24N3` |
| N-methyl diaminobenzamide pillar | 3 | `3 × C8H9N3O` |
| urea carbonyl | 6 | `6 × CO` |
| **Total** |  | **`C60H75N15O9`** |

## Graph checks

The deposited topology passes every preregistered architecture check:

- five independent six-member aromatic carbon rings;
- two fully substituted triethyl/tris(methylene) cap rings;
- three diaminobenzamide pillar rings;
- six urea carbonyls;
- three amide carbonyls;
- no other oxygenated carbons;
- three ethyl and three methylene-to-urea branches on each cap;
- one N-methyl benzamide on each pillar;
- six ureas forming a complete `K2,3` network between the two caps and three pillars.

Every cap is connected once to every pillar. The observed six cap–pillar pairs exactly equal the expected complete bipartite graph.

## Cross-system identity

The same architecture is present in the glucose, galactose, fructopyranose and fructofuranose topologies.

Common architecture fingerprint:

`641d6d491682d17e3f170d8340ce8087ae0eb1e9c6d57e420715c0c7a5b4f6a5`

This is stronger than formula agreement alone. It confirms the specific two-cap/three-pillar/six-urea connectivity and the methyl-truncated outer amides.

## Confidence effect

The public benchmark host is no longer an unresolved anonymous 159-atom graph. At molecular-graph and formula level, it is reconciled to the published methyl-truncated GluHUT architecture.

This supports:

- use of the exact compiled host as the common public benchmark baseline;
- atom-level comparison against future author-supplied source files;
- identification of the three N-methyl truncation sites and six urea bridges;
- precise provenance checks for charges and GAFF terms.

It does not establish:

- the original atom-name map used by the authors;
- RESP input conformers, electrostatic potentials or fitted charges;
- source `mol2`, `frcmod`, `lib` or tleap files;
- chemical correctness of the compiled numerical force field;
- transferability to PaG, MC-Q1 or complete NNC2215.

## Decision

1. Mark the public host **source-architecture reconciled at graph and formula level**.
2. Keep source parameter derivation and starting-coordinate provenance open.
3. Do not regenerate prospective PaG or MC-Q1 simulations until the charge and parameter derivation is obtained or independently rebuilt and validated.
4. Keep candidate ranking disabled.

## Reproducibility

- Script: `research/reconcile_gluhut_truncated_architecture.py`.
- Machine result: `research/data/gluhut-truncated-architecture-reconciliation.json`.
- Tests: `tests/test_reconcile_gluhut_truncated_architecture.py`.
- Workflow: `.github/workflows/gluhut-source-architecture-reconciliation.yml`.
- Public evidence PR: #39.
- Archive MD5: `9f3c6733dc608d2bd37c488c560a5443`.

## Primary sources

- Davis et al., *A biomimetic receptor for glucose*, DOI `10.1038/s41557-018-0155-z`.
- Eades et al., *Free energy simulations on a biomimetic glucose receptor*, DOI `10.1039/D6SC00738D`.

## Evidence boundary

This is a graph-provenance result for the methyl-truncated public receptor. It is not full receptor 2, NNC2215, PaG binding, MC-Q1 binding, efficacy, safety or a working-drug result.
