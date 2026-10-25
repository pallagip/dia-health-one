# NNC2215 Roof Macrocycle Graph Design

## Goal

Create a reproducible, source-attributed molecular graph for the NNC2215 roof macrocycle propyl azide so that chemistry work can progress from coordinate inventories toward an auditable working model.

## Evidence and scope

The graph is transcribed from the numbered roof propyl-azide structure in Extended Data Fig. 8i of the 2024 NNC2215 paper and checked against its Supporting Information and reported formula, C56H63N15O12. The source numbers the three-carbon roof tail C21–C23; remaining stable map IDs identify atoms in the graph and are documented by structural group.

This milestone covers the isolated macrocycle roof propyl azide only. It does not construct the insulin conjugate, resolve all source-to-PDB atom correspondences, parameterize the macrocycle, or authorize target molecular dynamics.

## Representation

Store one atom-mapped isomeric SMILES in a plain-text `.smi` file. Keep the source-drawn neutral carboxylic-acid form, and document the conventional azide resonance representation used in the SMILES. A small RDKit script parses and sanitizes the graph, verifies the formula, 83 heavy atoms, unique heavy-atom map IDs, zero formal charge, and expected ring count, then writes a provenance-rich JSON audit with the source and artifact hashes and RDKit version.

The map IDs are stable graph identifiers, not a claim that the source paper numbers every atom or that each ID has been mapped to a PDB atom. No stereocentres are assigned unless the source drawing explicitly marks one.

## Public presentation

Expose the graph file and audit from the research dashboard. Label it as a source-derived macrocycle graph candidate and display the neutral-acid/resonance-form assumptions and remaining PDB-map and force-field limits. Preserve the dashboard's existing statement that production parameters are unvalidated.

## Acceptance checks

- RDKit parses and sanitizes the mapped SMILES without errors.
- Formula is C56H63N15O12; heavy-atom count is 83; all heavy-atom map IDs are present exactly once; formal charge is zero in the depicted neutral-acid representation.
- The audit records SHA-256 for the SMILES, software version, source citations, assumptions, and explicit non-readiness limitations.
- Dashboard links resolve to the graph and audit, and copy accurately describes what has and has not been validated.
