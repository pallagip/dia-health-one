# NNC2215 source-derived full conjugate graph

**Date:** 26 October 2026  
**Status:** Scoped continuation of the authorized NNC2215 modeling research.

## Goal

Construct a source-derived neutral constitutional graph for the full desB30 insulin–NNC2215 conjugate and audit its sequence, connectivity, formula, atom valence, provenance, and local attachment boundaries. The graph is a reproducible chemistry representation, not an experimentally solved structure or simulation-ready model.

## Approaches considered

1. **Build full coordinates and select one B29 mapping.** This would use geometric choices that are not supported by the existing 16 exact B29 mappings per state. It is out of scope.
2. **Check only formula arithmetic.** This repeats existing source-formula work and would not close the open full-graph connectivity gate.
3. **Assemble the source-derived constitutional graph with RDKit FASTA residue templates and the mapped conjugate fragments.** This is selected: it adds auditable residue-level provenance and valence checks without using coordinates to resolve mapping ambiguity.

## Inputs and construction

- Read A- and B-chain sequences from the existing formula audit. Construct each chain with the pinned environment's RDKit `MolFromFASTA`, retaining RDKit residue and atom-name metadata.
- Join the three cysteine sulfur pairs specified by the insulin source sequence and cross-check them against both corrected candidate coordinate files: A6–A11, A7–B7, and A20–B19.
- Parse the mapped B1 glucoside graph. Remove map 1, the PheB1 N boundary placeholder, and connect its map-2 carbonyl carbon to the PheB1 N atom in the B-chain graph.
- Parse the mapped B29 linker–roof product graph. Remove map 92, the LysB29 NZ terminal cap, and connect its map-84 carbonyl carbon to LysB29 NZ in the B-chain graph. Preserve all remaining source map IDs on the appended nodes.
- Keep hydrogen atoms implicit and report each atom's implicit hydrogen count. This avoids reusing disputed coordinate-file hydrogen ownership. Apply RDKit sanitization and explicit valence checks to the assembled graph.
- Record per-node origin: A/B chain residue and template atom name, B1 source map, or B29 source map. Record source file hashes and the pinned RDKit version.

## Outputs

- Add `research/audit_nnc2215_full_conjugate_graph.py` to assemble and audit the graph.
- Add `research/data/nnc2215-full-conjugate-graph.json` with the atom and bond ledger, provenance, source/input hashes, formula and element counts, connected components, cycle rank, valence summary, state-by-state PDB sequence/LINK/disulfide cross-checks, retained B1/B29 mapping counts, and scoped assessment.
- Update `research/README.md`, `research/NNC2215_CHEMISTRY_RECONCILIATION.md`, and `research/index.html` with the audit result, rerun command, and JSON download link.

## Required checks

1. Sequence inputs equal A `GIVEQCCTSICSLYQLENYCN` and B `FVNQHLCGSHLVEALYLVCGERGFFYTPK` and match both corrected candidate PDB chain sequences.
2. Exactly three source disulfides are present: A6–A11, A7–B7, and A20–B19.
3. The assembled graph has one connected component, passes RDKit sanitization/valence checks, and reports its calculated cycle rank from `E − V + 1`.
4. The neutral graph formula equals `C325H461N79O97S6`, matching the published reported formula and the existing source formula audit. Report this as a graph/formula consistency check, not independent experimental validation.
5. B1 map 1 is replaced only by PheB1 N, and B29 map 92 only by LysB29 NZ; the source boundary atoms and corresponding candidate PDB LINK records agree in open and closed states.
6. Existing audits retain one exact B1 mapping and sixteen exact B29 mappings per state. Do not rank or choose among the B29 mappings.
7. The JSON identifies the distinct conditional tri-carboxylate audit as a separate graph representation and does not merge its charge state into this neutral graph.

## Limits and unresolved gates

This graph does not resolve PheB1 hydrogen ownership, protonation microstates, total charge at physiological pH, histidine tautomer choice beyond the pinned RDKit neutral template, roof stereochemistry, B29 coordinate-map ambiguity, neutral B29 oxygen-label transfers, force-field parameters, structural relaxation, target-binding MD, drug efficacy, or clinical readiness. The corrected PDB files are used only for sequence and named connectivity cross-checks; no coordinates or atom maps are changed.
