# NNC2215–insulin force-field review handoff

**Prepared:** 26 October 2026
**Evidence baseline:** main commit `d5c5909`
**Readiness:** candidate fragments only; no validated full-conjugate force field or OpenMM target system

## Purpose

This packet collects the current structure and parameter evidence for an independent computational-chemistry review. It is intended to turn the remaining model questions into explicit decisions before anyone builds or runs a target system. The open and closed NNC2215 coordinates are publisher-supplied model-built candidates, not experimentally resolved NNC2215 structures.

NNC2215 is an existing published insulin candidate. This work advances a computational model of that candidate; it does not create a new drug or establish efficacy.

## Evidence and review decisions

| Area | Current evidence | Review needed |
|---|---|---|
| Full graph and coordinates | The source-derived neutral graph has 507 heavy atoms (398 insulin, 109 conjugate). Both corrected candidates pass heavy-atom coverage; B1 has one exact mapping and B29 retains 16 mappings per state. | Confirm source interpretation, roof stereochemistry, candidate coordinate edits, and whether B29 alternatives are physically equivalent for parameter assignment. |
| Attachment boundaries | B1, LysB29, and roof/linker joins are present in both candidates. The boundary audit reports no PDB changes. | Confirm the protein-bound amide chemistry and exact boundary atoms used for custom residues. |
| PheB1 hydrogen | Open H2 is 1.012274 Å from PheB1 N. Closed PHX H2A is 1.006376 Å away, with no N CONECT edge; the local angles are similar. | Decide whether the closed H2A belongs to PheB1 N, using chemistry/source evidence. Geometry alone does not decide this. |
| Protonation and total charge | The pH 7.4 ledger enumerates 144 group-charge combinations over the declared sites. It is group-level bookkeeping, not pKa or population evidence and not an atom-localized charge model. | Select and justify modeled microstate(s), including termini, histidines, and the three macrocycle carboxyls; report the full-system formal charge. |
| Standard protein templates | The Amber14 inventory finds literal template-name candidates for 16 of 17 observed residue labels; HIE and CYX are graph-context candidates. Modified PheB1, LysB29, and cross-boundary terms need custom treatment. | Confirm exact template, terminal, disulfide, and charge compatibility in the assembled system. |
| B1 parameter pilot | GAFF2/AM1-BCC candidate for a neutral 37-heavy-atom local fragment; charge residual +0.003002 e. The AsnB3-NH2 truncation cap and model-built geometry limit transferability. | Review charge partitioning, cap effects, missing junction terms, and compatibility with the protein force field. |
| B29 parameter pilot | GAFF2/AM1-BCC candidate for a 92-heavy-atom fragment in a declared −3 state. Across 16 coordinate mappings, total charge sums range from −2.997001 to −2.996005 e; no correction was applied. Two generated improper records have penalty score 41.9. | Review mapping-dependent charges, the endpoint cap, both high-penalty impropers, and all custom/cross-boundary terms. |
| Full-system readiness | OpenFF Sage label coverage and Amber template inventories exist, but full-system partial charges and custom parameters have not been validated. | Decide whether the evidence supports creating an OpenMM system and define acceptance criteria before any dynamics. |

## Source artifacts

- [Full conjugate graph and provenance](data/nnc2215-full-conjugate-graph.json) and [full coordinate mapping](data/nnc2215-full-coordinate-mapping.json)
- [B1/B29 boundary audit](data/nnc2215-attachment-boundary-audit.json)
- [Conditional charge-state ledger](data/nnc2215-physiological-charge-states.json)
- [OpenFF label coverage](data/nnc2215-forcefield-label-coverage.json)
- [Amber template feasibility](data/nnc2215-hybrid-feasibility.json)
- [B1 fragment pilot](data/nnc2215-b1-gaff2-boundary-pilot.json) and [B29 fragment pilot](data/nnc2215-b29-gaff2-pilot.json)
- [B29 mapping sensitivity](data/nnc2215-b29-gaff2-mapping-sensitivity.json)
- [Proposed target-MD protocol and stop criteria](data/target-md-pilot-protocol.json)
- [Source and coordinate provenance](data/structure-manifest.json)

## Gate to an initial target system

1. Record an expert-reviewed, atom-mapped chemical graph, stereochemistry, selected microstate, and total formal charge for the complete conjugate.
2. Resolve the B1 and B29 custom-residue boundaries and review charges, bonded terms, impropers, and cross-boundary interactions. Preserve alternatives when the evidence does not select one mapping.
3. Assemble the chosen open or closed candidate in OpenMM. Verify atom and bond coverage, parameter completeness, charge totals, clashes, and finite energy before dynamics.
4. After a valid system and a PaG-pose method are reviewed, run restrained minimization/equilibration and a short stability check. Scale to the planned independent-seed pilot only if those checks pass.

The exact 2026 PaG simulation remains unavailable for reproduction because its in-house starting coordinates/system were not present in the reviewed public materials. Any use of the public NNC2215 candidates is a separately labelled proxy study. No target-binding MD, affinity estimate, or clinical conclusion has been produced.

## Reviewer decision record

- Reviewer / relevant expertise:
- Date and evidence baseline:
- Graph and stereochemistry decision:
- Selected microstate and formal charge:
- B1/B29 boundary and hydrogen decisions:
- Parameter approval, required changes, or rejection:
- Conditions required before OpenMM system creation:
