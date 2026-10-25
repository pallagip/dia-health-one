# NNC2215 Roof Graph to Coordinate Candidate Mapping Design

## Goal

Measure how the source-derived roof macrocycle graph fits the existing open and closed corrected-coordinate candidates, and publish the result without presenting an inferred atom map as authoritative.

## Method

Read each corrected candidate PDB, extract heavy atoms and adjacency from CONECT plus LINK records, then use an element-labeled NetworkX subgraph monomorphism against the RDKit graph. Apply only the previously documented C21–C23 and azide-tail candidate anchors from the linker-mapping audit. Preserve all valid mappings and report PDB atoms that remain outside the macrocycle graph.

The check compares adjacency only. It ignores bond orders, formal charges, protonation, force-field atom types, and coordinates beyond the existing source-corrected candidate PDBs. It does not edit the PDB files or promote their candidate status.

## Output

Write one JSON audit with input hashes, source-tail anchor provenance, per-state formula/count inventory, all matching graph maps, per-map ambiguity, unmatched linker atoms, and explicit limitations. Update the dashboard and chemistry notes to report the exact match count and the still-ambiguous atom assignments.

## Acceptance checks

- Both coordinate candidates expose the expected 91-heavy-atom B29-linked fragment (C62N15O14).
- Each candidate contains the complete 83-atom source graph as an element-matched, edge-preserving subgraph after applying the documented tail anchors.
- The audit records every valid candidate embedding, all unmatched atoms, and the source of each fixed anchor.
- No candidate PDB coordinate or atom record is changed by the audit.
- Public copy explicitly says the graph-to-PDB assignment remains ambiguous and does not validate bond order or force-field chemistry.
