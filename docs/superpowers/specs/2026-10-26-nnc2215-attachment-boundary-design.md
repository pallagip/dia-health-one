# NNC2215 attachment-boundary closure audit

**Date:** 26 October 2026

**Status:** The audit scope was approved; implementation awaits review of this written specification.

## Goal

Audit how the source-derived B1 glucoside and B29 linker–roof components join to the local insulin residue topology in the corrected open and closed NNC2215 candidate coordinates. Establish the heavy-atom boundary edges and their local geometry, and explicitly report the PheB1 nitrogen hydrogen-record discrepancy. Preserve the coordinate candidates and do not infer or repair hydrogen ownership.

## Current evidence and gap

- The source-mapped B1 graph has one exact adjacency mapping in each corrected candidate. Its boundary is B1 map 1 (PheB1 N) to map 2 (the B1 linker carbonyl carbon), observed as PheB1 N–PHX C2 in each candidate.
- The source-mapped B29 graph matches 16 exact adjacency mappings per state. It includes LysB29 NZ but intentionally excludes LysB29 CE, the native residue context.
- The B29 source component maps linker carbonyl C3 as map 84 and LysB29 NZ as map 92; their candidate link is LYS B29 NZ–LYX B30 CD1. The roof-tail edge joins source maps 78–79 (C21–C22) across the MC/LYX residue boundary. Its candidate PDB atom name is LYX CD6 in the open model and LYX C in the closed model.
- The corrected candidates record all three expected heavy-atom LINK edges: PheB1 N–PHX C2, LysB29 NZ–LYX CD1, and roof MC C31–LYX CD6/C.
- The current global structure and connectivity audits describe the original publisher PDB records. The component audits use the corrected candidates, so the local attachment boundaries have not yet been reconciled in one report.
- PheB1 N has a CONECT hydrogen neighbor named H2 in the open candidate. The closed candidate has no CONECT hydrogen neighbor at that N, while PHX H2A lies 1.006 Å from PheB1 N without a CONECT edge. These records do not establish the chemically correct hydrogen ownership.

## Proposed design

Add a focused audit script that reads the two corrected candidate PDBs, the source B1/B29 graph and map files, and the existing topology-audit JSONs. For each state, it will resolve the three expected heavy-atom joins by full atom identity, check that the corresponding source-boundary edges and native protein bonds are present, and calculate distances without using them to rank or select mappings.

The report will include a per-state boundary ledger with source map IDs, candidate PDB atom identities, edge provenance, and measured distances. It will also record the local neighbors and CONECT hydrogen records at PheB1 N and LysB29 NZ. The expected local neutral boundary patterns are PheB1 N bonded to CA, the B1 carbonyl carbon, and one hydrogen; LysB29 NZ bonded to CE, the B29 carbamate carbonyl carbon, and one hydrogen. Any mismatch in the observed hydrogen records will be surfaced as an open item and will not be repaired or silently attributed to another atom. For the open-state PheB1 H2 and the closed-state PHX H2A coordinate, the audit reports N–H distances, local angles, and nearest heavy-atom identities as a descriptive comparison. These measurements must not create a connectivity edge or assign hydrogen ownership.

## Invariants and output

The audit must assert that both candidate input hashes match the reconciliation record and remain unchanged during the run; each expected heavy-atom boundary is unique; source component graphs and their existing exact mapping counts remain unchanged; and the PheB1 N–CA and LysB29 CE–NZ native residue bonds resolve uniquely. A failure of a required heavy-atom join or identity check must prevent a “boundary closure” pass. Hydrogen-record discrepancies are reported separately and must not invalidate the heavy-atom adjacency result.

Write `research/data/nnc2215-attachment-boundary-audit.json` with source paths and hashes, per-state boundary rows, coordinate distances, local heavy-atom and hydrogen neighbors, explicit discrepancies, and a scoped assessment. Add a rerun command and download link to the research documentation and dashboard. Do not alter candidate PDBs or the B1/B29 source graphs.

## Scope limits

This audit closes only the local heavy-atom attachment boundaries and characterizes their recorded hydrogen neighbors. It does not determine a unique full B29 atom mapping, assign B29 carboxyl oxygen names, establish complete-conjugate bond orders or stereochemistry, compute a full insulin protonation/charge ledger, validate force-field parameters, run target-binding MD, or establish a clinical conclusion. The corrected files remain model-built candidates rather than experimentally solved NNC2215 structures.

## Acceptance criteria

1. Open and closed candidates each resolve the B1 attachment link, B29 Lys attachment link, and roof-tail inter-residue link by unique atom identities.
2. Source boundary edges and native PheB1 N–CA / LysB29 CE–NZ context are shown explicitly and cross-checked against the existing B1/B29 audits.
3. Geometry is reported as observation only; no distance is used to choose a graph mapping or hydrogen assignment.
4. The PheB1 H2/H2A discrepancy and LysB29 NZ hydrogen neighbors are reported per state; the coordinate-geometry comparison remains descriptive and does not assign an owner; no atom or coordinate is changed.
5. Existing B1 and B29 audit results, input hashes, and B29's 16 mappings per state are preserved.
6. The dashboard and README distinguish verified heavy-atom boundary closure from unresolved hydrogen chemistry and full-model readiness.
