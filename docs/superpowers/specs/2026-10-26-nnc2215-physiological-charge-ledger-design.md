# NNC2215 conditional pH 7.4 charge-state ledger

**Date:** 26 October 2026

**Status:** Authorized next modeling gate; assumptions remain explicit and reversible.

## Goal

Turn the source-derived neutral full-conjugate graph into a machine-readable inventory of plausible formal charge and hydrogen-count states around pH 7.4. Provide a clearly labeled canonical candidate for downstream model preparation while enumerating the uncertain macrocycle carboxyl, histidine, and A-chain N-terminal states.

This is a charge-state bookkeeping model. It does not predict pKa values or microstate populations, localize partial charges, assign force-field atom types, or claim that any one state is measured.

## Selected approach

Use the merged full-conjugate atom ledger as the sole graph input. Resolve ionizable sites by residue identity and mapped-source IDs, count proton changes relative to the neutral graph, and enumerate combinations at the group level. Keep carboxylate resonance and histidine tautomer ambiguity explicit instead of selecting PDB oxygen names or claiming atom-specific force-field charges.

At pH 7.4, provide one conventional candidate: six insulin carboxyl groups (four Glu sidechains plus two chain C termini) deprotonated, three B29 macrocycle carboxyl groups deprotonated, ArgB22 protonated, A-chain GlyA1 N terminus protonated, and both histidines neutral in their pinned RDKit-template tautomers. This candidate has net formal charge −7 and formula C325H454N79O97S6. It is a bookkeeping starting point, not a measured or probability-ranked microstate.

Enumerate 144 states from three binary macrocycle-carboxyl states, two A1 amino states, and three states for each histidine (neutral template tautomer, neutral alternate tautomer, protonated). Fix conventional ionization assumptions for the four Glu groups, both chain C termini, and ArgB22; disclose that local pKa shifts are not modeled. Do not enumerate tyrosine deprotonation, which is outside the usual pH 7.4 baseline used here; record it as an excluded sensitivity.

## Inputs and outputs

- Input: `research/data/nnc2215-full-conjugate-graph.json` and its provenance fields.
- Create `research/audit_nnc2215_physiological_charge_states.py` to validate the source atom identities and emit the site ledger, explicit assumptions, state enumeration, formula/hydrogen accounting, and scope limits.
- Create `research/data/nnc2215-physiological-charge-states.json` with canonical candidate state, 144 enumerated group-state combinations, charge/formula ranges, input hash, and checks.
- Update the research README, chemistry reconciliation, and dashboard with the result and rerun command.

## Required checks

1. The input graph hash and formula are recorded; no coordinate data are read.
2. Exactly four insulin Glu sidechains, two chain C-terminal carboxyls, three mapped B29 macrocycle carboxyls, two histidines, A-chain N terminus, and ArgB22 are found by provenance.
3. PheB1 N and LysB29 NZ remain neutral conjugate amides/carbamates in this ledger; the disputed PheB1 hydrogen ownership is not resolved.
4. The canonical candidate's net charge is −7 and hydrogen count is 454 relative to the neutral source graph formula C325H461N79O97S6.
5. Enumeration yields exactly 144 combinations without assigning state probabilities; total charge spans −8 to −2 and hydrogen count spans 453 to 459 under the declared fixed-group assumptions.
6. The ledger states that group formal charges are not atom-localized partial charges and do not make the conjugate force-field or MD ready.

## Open gates

Local pKa shifts, measured protonation, partial atomic charges, histidine sidechain/tautomer energetics, tyrosine sensitivity, force-field parameters for the conjugate, coordinate mapping and hydrogen ownership remain open. The conditional tri-carboxylate graph stays separate from the neutral constitutional graph.
