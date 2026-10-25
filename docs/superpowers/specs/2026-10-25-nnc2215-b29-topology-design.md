# NNC2215 B29 Topology Audit Design

**Date:** 2026-10-25  
**Status:** User approved this scope; review this design before implementation.

## Goal

Create an explicit, source-mapped graph for the NNC2215 B29 linker and roof macrocycle after the azide–alkyne click reaction, then check whether the graph is consistent with both corrected open and closed coordinate candidates. This is the next chemistry gate before full-conjugate topology construction, parameterization, or molecular dynamics.

## Current evidence

- The source-drawn roof propyl-azide graph has formula C56H63N15O12, 83 heavy atoms, 89 heavy-atom edges, and cycle rank 7.
- The corrected candidate PDB fragments containing LYX and MC have 91 heavy atoms (heavy-element inventory C62N15O14), 98 edges, one connected component, and cycle rank 8. Including the LysB29 NZ attachment atom and its LINK edge gives 92 heavy atoms (heavy-element inventory C62N16O14), 99 edges, and cycle rank 8.
- The eight atoms outside the roof graph are assigned to a C6O2 linker under the existing candidate partition. The count change is consistent with a connected, acyclic eight-atom linker plus two click-product connections to the azide, creating the triazole cycle. Counts alone do not establish bond order or regioisomer.
- A read-only symmetry calculation found that the 16 existing element/adjacency embeddings per state use the same 83 PDB atoms and equal the anchored automorphisms of the source graph. Preserving the source graph's bond labels leaves two source-label permutations. This preliminary result must be reproduced by the audit. It does not assign target PDB bond orders.
- The PDBs are author-supplied model-built candidates, not experimentally resolved NNC2215 structures. PDB CONECT/LINK records encode adjacency, not bond order or formal charge.

Primary chemistry references are the [2024 NNC2215 article](https://www.nature.com/articles/s41586-024-08042-3), its [Supporting Information](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-024-08042-3/MediaObjects/41586_2024_8042_MOESM1_ESM.pdf), and the [WO2023144240A1 patent](https://patents.google.com/patent/WO2023144240A1/en). The public model coordinates and their provenance are in the [Figshare dataset](https://springernature.figshare.com/articles/dataset/3D_models_of_NNC2215/26526460).

## Proposed design

1. Build a versioned atom-mapped product graph from the published B29 carbamate linker, roof macrocycle, and product drawing for the azide–alkyne coupling. Preserve source atom identifiers where supported; record every transformation from the isolated azide and alkyne reactants to the product. Treat the existing tail anchors and the two-carbon coordinate correction as candidate assumptions, not source-authoritative atom names.
2. Include the LysB29 NZ attachment node and bond so the checked component covers the full B29 attachment. Keep the graph's protonation state explicit as a representation choice; defer pH microstates, partial charges, and force-field atom types.
3. Match the complete product graph one-to-one against both corrected PDB components by element and adjacency, including the LysB29 NZ endpoint. Transfer bond-order labels from the source-derived product graph only after that mapping; the PDB records themselves do not supply those labels. Represent symmetry-equivalent source maps as equivalence classes and do not use adjacency-only monomorphism to claim a complete product match.
4. Audit atom and edge counts, formula, connectedness, cycle rank, valence, carbamate connectivity, triazole ring membership/regiochemistry, and the B29 attachment. Record coordinate distances for each proposed heavy-atom bond in both states against documented broad bond-class plausibility intervals. Distances are a geometry screen; they cannot independently establish bond order.
5. Compare the open and closed graphs, hashes, and hydrogen-neighbor records. Preserve all input PDBs byte-for-byte. Report mismatches as unresolved instead of repairing coordinates or silently selecting a candidate.

## Deliverables

- A mapped product graph artifact for the B29 linker–macrocycle component, with an atom-map legend and explicit source references.
- A reproducible audit script and JSON record containing the source graph, mappings up to symmetry, candidate PDB hashes, graph and valence summaries, bond-distance observations, and unresolved assumptions.
- A concise update to the chemistry reconciliation report and research dashboard only after the graph audit has a defensible result.

## Acceptance criteria

- The product graph has the expected heavy-atom inventory C62N16O14 when the LysB29 NZ attachment is included, 92 heavy atoms, 99 edges, one connected component, and cycle rank 8; any deviation is investigated and reported.
- The published product drawing supports the represented click-product regiochemistry. If it does not resolve the regioisomer, preserve both chemically plausible graphs and stop short of a unique assignment.
- Every graph atom maps one-to-one onto the intended atoms in each candidate PDB; no atom or edge is silently dropped. Symmetry-equivalent mappings are reported as such.
- All assigned valences are valid for the documented bond-order representation. Bond-distance outliers and open/closed connectivity or hydrogen differences are surfaced without changing the PDBs.
- Provenance and PDB hashes match the existing reconciliation record, and the run can be reproduced from checked-in inputs.

## Out of scope

- Building or validating the B1 glucoside component or the entire insulin conjugate.
- Assigning physiological protonation microstates, partial charges, or force-field parameters.
- Generating a production topology, modifying coordinates, minimizing a system, or running molecular dynamics.
- Inferring glucose or PaG affinity, receptor activation, clinical efficacy, or safety from these coordinates.
- Reproducing the 2026 PaG simulation, whose reported in-house seed coordinates are not in the reviewed public input set.

## Alternatives considered

- **Proceed with the B29 topology audit (selected):** resolves a bounded, chemistry-critical part of the coordinate model and establishes whether later parameterization is justified.
- **Attempt full-conjugate force-field preparation now:** could appear to advance quickly, but it would combine the unresolved B1 graph, B29 bond orders, protonation, and parameters into a system whose energies or dynamics could not be interpreted.
- **Wait for the original PaG study package:** would be needed for an exact reproduction, but the reviewed public paper and SI do not provide its in-house seed coordinates or simulation files. It would not independently resolve the public NNC2215 candidate topology.

## Limits

Even a successful audit would establish a source-supported candidate topology for one conjugate component, not an experimental structure, a validated force field, target-binding simulation readiness, or a drug candidate. The later gates remain full B1/insulin topology, protonation and charge review, custom parameter validation, system stability, and a defensible PaG pose.
