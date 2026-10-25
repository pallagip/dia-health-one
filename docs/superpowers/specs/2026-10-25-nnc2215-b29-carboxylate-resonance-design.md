# NNC2215 B29 carboxylate resonance audit

## Goal

Report the source-to-coordinate B29 bond-label ambiguity for a declared pH 7.4 charge state, while preserving the existing neutral-acid graph audit.

## Context

The current neutral source graph has three macrocycle carboxylic acids at map groups `[52, 53, 54]`, `[60, 61, 62]`, and `[68, 69, 70]`. The eight current bond-label transfer classes differ only by swapping the two oxygen map assignments within each group. The two exact mappings inside each class are related by a source graph automorphism. This is a mapping ambiguity, not eight different linker or click topologies.

The available source record does not provide site-specific pKa measurements for these three groups. Therefore pH 7.4 is a declared modeling scenario, and the fully deprotonated state is an explicit assumption to be checked later during parameterization.

## Design

- Keep the neutral source SMILES, its 16 exact graph mappings per coordinate state, and its eight neutral bond-label classes unchanged.
- Add a separate resonance-normalized graph comparison for a scenario with all three macrocycle carboxyl groups deprotonated. In this comparison only, the two C–O edges in each carboxylate are treated as equivalent; all other source bond labels remain unchanged.
- Enumerate the same exact coordinate mappings and group them under the normalized graph automorphisms. Record the normalized class count and the atom-map ambiguity separately from the neutral results.
- Record the assumed local charge change as three removed acid protons and net charge −3 for the isolated B29 component, while retaining the LysB29 NZ NH2 boundary cap. Do not claim this capped component is the full conjugate charge.
- Do not modify either PDB or select an oxygen mapping from coordinate distances. Do not create force-field charges or parameter files.
- Keep the dashboard's overall full-chemistry gate open. Explain that the oxygen-label classes collapse for the declared carboxylate state while full atom naming, complete-conjugate protonation, and force-field readiness remain open.

## Acceptance criteria

1. The audit verifies exactly the three declared C/O map groups and checks the neutral reference formula before applying the stated hydrogen/charge delta.
2. For each state, the neutral analysis still reports 16 exact mappings and eight bond-label transfer classes.
3. For each state, the resonance-normalized comparison reports 16 exact mappings and one chemically equivalent bond-label class.
4. The output explicitly labels pH 7.4 and full deprotonation as a modeling assumption, not measured pKa evidence.
5. The output keeps atom-label correspondence, complete-conjugate protonation, charges, and force-field parameters open; coordinate hashes remain unchanged.
6. The research page and chemistry note distinguish neutral mapping classes from the conditional carboxylate equivalence result.

## Risks and boundaries

The normalization is a graph-equivalence analysis, not a validated protonation calculation. A future full conjugate model still needs a site-specific protonation inventory and parameterization. The PDBs are model-built structures, and their oxygen names are not authoritative source labels.
