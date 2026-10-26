# NNC2215 hybrid force-field feasibility inventory

**Date:** 26 October 2026  
**Status:** Approved direction; scoped feasibility audit, not a parameterization.

## Goal

Determine whether a candidate hybrid route can retain standard Amber14SB residue templates for the insulin sequence while separately identifying the NNC2215 conjugate atoms and the modified insulin attachment sites that need custom parameters.

## Approaches considered

1. **Assign NAGL charges to the entire 507-heavy-atom graph.** The pinned Sage force field names the NAGL AM1-BCC model, but the model documentation reports no test molecules above 89 heavy atoms. A whole-graph charge vector would be a large size extrapolation, so this is not the route to a simulation-ready model.
2. **Use Amber14SB templates for the entire graph.** The Amber force field contains templates for the standard insulin amino-acid residue names, but the B1 glucoside, B29 carbamate-linker, and macrocycle are not ordinary amino-acid templates. The two modified attachment residues also cannot be treated as unchanged PHE and LYS templates.
3. **Inventory a hybrid Amber plus custom-conjugate route.** This is selected. It gives a concrete coverage map and identifies the custom terms and chemistry inputs required before parameterization.

## Inputs and method

- Read the source-derived neutral full-graph JSON, the existing Sage label-coverage JSON, and the pinned environment lock.
- Load the installed OpenMM `amber14-all.xml` include tree. Record every included XML file hash and the set of residue-template names it defines.
- Group graph atoms by insulin chain, residue number, and residue name using their recorded provenance. Check which sequence residue names occur in the Amber template inventory.
- Identify graph edges that cross from an insulin-template atom to a mapped conjugate atom. Report the B1 Phe N to mapped B1 linker carbonyl join and the LysB29 NZ to mapped B29 carbamate carbonyl join by source node and map identifiers.
- Mark PHEB1 and LYSB29 as modified sites even when ordinary PHE/LYS templates exist. Mark the mapped B1 and B29 atoms as custom-conjugate scope.
- Carry forward the NAGL model's published size-domain evidence and current open graph, protonation, stereochemistry, atom-mapping, and charge-assignment gates as provenance, not as computed parameter results.

## Outputs and acceptance

- Add `research/audit_nnc2215_hybrid_feasibility.py` and generate `research/data/nnc2215-hybrid-feasibility.json`.
- Report formula, heavy-atom and bond counts; the 50 insulin residue identities and their template-name availability; modified PHEB1 and LYSB29 sites; both cross-origin conjugate joins; mapped-conjugate atom count; hashes and package versions; and a scoped go/no-go assessment.
- The audit passes its inventory checks only if every graph residue is classified, all ordinary insulin residue names are present in the Amber template inventory, and exactly the two expected insulin-to-conjugate joins resolve by provenance.
- Use status `STANDARD_RESIDUE_TEMPLATE_NAMES_PRESENT; CUSTOM_CONJUGATE_PARAMETERS_REQUIRED` when those checks pass. This means a hybrid route is structurally identifiable, not that any force-field parameters have been generated or validated.
- Update the research README, compute status, chemistry reconciliation, milestone report, and dashboard with a rerun command and JSON link.

## Scope limits

This inventory does not call OpenMM's system builder, generate partial charges, choose a protonation microstate, create or validate custom bonded/nonbonded parameters, use candidate coordinates, build a solvated system, run MD, estimate binding, or support a clinical inference. Amber template-name presence does not prove an exact residue-template match. PHEB1 hydrogen ownership, the unique B29 coordinate map, roof stereochemistry, local protonation, and compatible cross-boundary parameters remain open. The output must keep the validated conjugate-force-field count at zero.
