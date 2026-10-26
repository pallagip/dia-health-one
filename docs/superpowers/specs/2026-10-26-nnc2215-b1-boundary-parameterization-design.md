# NNC2215 B1 Local Boundary Parameterization Design

## Goal

Measure whether the pinned AmberTools GAFF2/AM1-BCC workflow can produce a charge and parameter candidate for the source-mapped B1 glucoside conjugate in a small, explicit PheB1 peptide context. This is a feasibility pilot toward the full insulin model, not a validated force field.

## Chosen scope

Extract one connected 37-heavy-atom local molecule from the existing neutral full-conjugate graph: the source-mapped B1 component (maps 2–19), PheB1, ValB2, and the AsnB3 backbone nitrogen as a neutral C-terminal amide cap. The B1 topology audit reports one exact source-to-coordinate embedding in the active/open candidate. Use that embedding for B1 atoms; map PheB1, ValB2, and the AsnB3 cap by their unique chain/residue/atom identities. The cap nitrogen keeps its active/open heavy-atom coordinate and gains hydrogens from RDKit after the downstream AsnB3 atoms are removed.

The source graph supplies the B1 stereochemistry and the pinned RDKit neutral insulin templates supply the Phe/Val stereochemistry and peptide bond orders. Reconstruct the selected graph from those audited inputs, verify that it is one connected neutral molecule with the expected atoms and bond edges, then copy active/open heavy-atom coordinates from the candidate PDB. Add hydrogens with RDKit. Do not use distances to choose atom mappings or infer bond orders.

Run Antechamber with AM1-BCC, GAFF2, declared formal charge 0, and `maxcyc=0`, followed by parmchk2 with ff14SB as a lookup set. Record input and tool hashes, per-atom names/maps/serials/types/charges, charge-sum residual, generated parameter records and penalties, and raw SDF/MOL2/frcmod/log artifacts. Check the candidate PDB hash before and after the run.

## Alternatives considered

1. **B1 linker-only fragment with an NH2 endpoint cap.** This omits the PheB1 CA bond and distorts the conjugated nitrogen environment, so it is less informative.
2. **Local PheB1/ValB2 context with an AsnB3-NH2 cap (selected).** This retains the actual PheB1 amide and adjacent peptide bond, uses the unique B1 mapping, and has a small declared C-terminal truncation. The cap and model-built geometry remain sources of uncertainty.
3. **Assign GAFF2 to the complete insulin–NNC2215 graph.** This would mix custom conjugate chemistry with unresolved B29 mapping, protonation, and full-system charge decisions, and would overstate current validation.

## Interpretation limits

The local fragment is a neutral, model-coordinate candidate. Its AM1-BCC charges depend on the chosen structure and AsnB3 cap. It does not resolve the PheB1 hydrogen-record discrepancy, define physiological microstates, validate parameter transferability or cross-boundary terms, assign the B29 fragment, parameterize a complete insulin system, establish a receptor pose, or support molecular dynamics or affinity claims. The validated-conjugate-force-field count remains zero.

## Acceptance evidence

- The extraction reproduces the expected atom inventory, formal charge, stereochemistry, and connected graph from the checked-in source records.
- All selected heavy atoms map once to the active/open model coordinates; B1 maps are the sole topology-audit embedding, and protein atoms use exact identities.
- Antechamber and parmchk2 outputs are retained with return codes, hashes, charge bookkeeping, and generated-term penalties; no charge correction is applied.
- The audit refuses stale or changed source, topology, or coordinate hashes and confirms the PDB bytes are unchanged.
- The report and dashboard state the cap, model-built coordinate scope, and limits directly.
