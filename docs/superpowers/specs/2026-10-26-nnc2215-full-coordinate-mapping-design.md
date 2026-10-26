# NNC2215 Full Conjugate Coordinate Mapping Design

## Goal

Create a reproducible, source-linked ledger that maps every heavy atom in the 507-atom neutral insulin–NNC2215 graph to the active/open and inactive/closed model-built coordinate candidates. Preserve every exact B29 graph-to-coordinate assignment as an ensemble so the ledger can support later topology work without selecting an arbitrary atom map.

This closes a coordinate-coverage and atom-identity gate. It does not validate the model-built conformations or produce a simulation-ready model.

## Current evidence and design decision

The source-derived neutral graph contains 507 heavy atoms: 398 atoms from the insulin residue templates and 109 mapped conjugate atoms. The corrected open candidate contains one additional, separate 12-heavy-atom BGL residue; after excluding that unlinked free glucose, both candidate structures have the same 507-heavy-atom elemental inventory as the graph. The existing full-graph audit validates sequences, three disulfide records, three attachment LINK records, and source hashes, but deliberately does not map graph atoms to coordinates.

The B1 topology audit has one exact element-labeled mapping per state. The B29 topology audit has 16 exact mappings per state. The selected approach composes these mappings with exact chain, residue, residue-name, atom-name, and element matches for the 398 insulin-template atoms. It checks the B1 PheB1 N and B29 LysB29 NZ endpoint identities and emits 16 full-graph mappings per candidate state. It never uses coordinates or distances to select, rank, or collapse those mappings.

## Chosen scope

Add `research/audit_nnc2215_full_coordinate_mapping.py` and its deterministic report `research/data/nnc2215-full-coordinate-mapping.json`. The audit will:

1. Load the checked-in full neutral graph, corrected open and closed candidate PDBs, coordinate reconciliation record, B1 topology audit, B29 topology audit, and attachment-boundary audit. Verify the relevant input hashes against the upstream audit records before mapping and confirm candidate PDB bytes are unchanged afterward.
2. Map each `insulin_chain_template` graph atom by its exact chain, residue number, residue name, atom name, and element identity. Require one match for each of the 398 atoms in each candidate.
3. Compose the unique B1 mapping and each of the 16 B29 mappings with the protein mapping. Map B1 source map 1 onto the existing PheB1 N graph node and B29 source map 84 onto the existing LysB29 NZ graph node; do not duplicate either attachment endpoint in the full graph. Require every full-graph node to map once and every mapped assignment to use 507 distinct candidate heavy-atom records.
4. Identify the open candidate's separate BGL residue by residue identity and atom inventory. Confirm it has no adjacency edge to the 507-atom graph, exclude it from graph mapping, and report its 12 heavy atoms as an unlinked extra component. Do not alter or rewrite either candidate PDB.
5. Store, for each state and each of its 16 full mappings, the mapping ID, full-graph node ID, atom provenance/source map, candidate PDB serial and identity, element, and observed XYZ coordinates in angstroms. Include per-state input hashes, mapping counts, atom coverage, uniqueness, and unmatched or extra atom inventories.
6. Record source-graph heavy-bond distances for each full mapping. For C/N/O pairs, apply the already disclosed broad heuristic screen from the B29 topology audit and report every out-of-range edge and its endpoint atoms. For the three S–S disulfide edges, record distances without assigning a threshold. These are geometry observations only; they cannot validate bond order, atom identity, mapping, stereochemistry, or chemical correctness.
7. Record the existing PheB1 hydrogen-record discrepancy and other known hydrogen-record differences as unresolved context. PDB hydrogens are not mapped to graph nodes and do not determine protonation or hydrogen ownership.

Update `research/README.md`, `research/COMPUTE_STATUS.md`, `research/NNC2215_CHEMISTRY_RECONCILIATION.md`, `research/MILESTONE_REPORT.md`, and `research/index.html` to link the report and rerun command and to describe the mapping coverage and limitations. The page must continue to state that coordinates are model-built and that no complete force field, OpenMM target system, target MD, affinity, or clinical result has been established.

## Error handling and integrity

The audit must stop with a clear failed or partial report if an upstream hash is stale, a required graph node or candidate atom is absent or duplicated, an element or endpoint identity disagrees, the expected 1/16 component mapping counts change, the BGL exclusion is not exact, or a full assignment does not cover 507 unique candidate heavy-atom serials. A failed report must not publish incomplete atom-to-coordinate assignments as complete.

Coordinates are copied as observations only. No coordinates, atom names, source graphs, candidate structures, or prior audit outputs are modified. Geometry flags are retained in the report and never used to reject, score, or choose among the 16 B29 mappings.

## Alternatives considered

1. **Begin integrated Amber topology construction now.** This is more direct toward simulation but depends on unresolved B29 mapping, whole-system microstate selection, and custom boundary terms. It could produce a superficially complete topology before atom-to-coordinate coverage is established.
2. **Expand the pH 7.4 charge-state ledger and place hydrogens.** This would add conditional chemical states, but would not close the missing full-graph coordinate mapping and would still lack boundary parameters.
3. **Full graph-to-coordinate mapping ledger (selected).** It closes a needed identity and coverage gate, retains mapping ambiguity explicitly, and produces a reusable atom-indexed scaffold record without implying chemical or force-field validation.

## Acceptance evidence

- Each candidate state maps all 398 insulin-template atoms by unique exact identity.
- The B1 topology mapping count is one and the B29 mapping count is sixteen in each state; the combined report contains sixteen complete 507-node assignments per state.
- In every assignment, each full-graph heavy atom maps once, each candidate graph atom serial is unique, the B1 and B29 attachment endpoints resolve to the existing protein nodes, and all source graph atoms have finite XYZ coordinates.
- The candidate heavy-atom inventory equals the graph inventory after the explicitly validated BGL exclusion in the open state; all unmatched and extra records are listed.
- The report records geometry observations and out-of-range heuristic screens without using them to select or rank mappings.
- Candidate PDB hashes match their upstream records before and after the audit. No candidate PDB or upstream audit file is modified.
- Documentation and dashboard links resolve locally and state that this coordinate ledger is not a parameterized or simulation-ready model.
- The audit is reproducible in the existing pinned research environment. No separate software test suite is added or run.

## Limits

Both PDB candidates are model-built derivatives, not experimentally resolved NNC2215 structures. The 16 B29 assignments remain unresolved; coordinate geometry is not evidence that selects one. This work does not assign bond orders from PDB records, settle the PheB1 hydrogen discrepancy, assign a full physiological microstate or total partial charges, construct complete parameters or cross-boundary terms, create an OpenMM system, or establish target-binding MD, affinity, pharmacological efficacy, or clinical readiness. The validated conjugate force-field count remains zero.
