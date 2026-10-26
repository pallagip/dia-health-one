# NNC2215 hybrid force-field feasibility inventory

**Date:** 26 October 2026  
**Status:** Approved direction; scoped feasibility audit, not a parameterization.

## Goal

Determine whether a candidate hybrid route can retain standard Amber14SB residue templates for the insulin sequence while separately identifying the NNC2215 conjugate atoms and the modified insulin attachment sites that need custom parameters.

## Approaches considered

1. **Assign NAGL charges to the entire 507-heavy-atom graph.** The pinned Sage force field names the NAGL AM1-BCC model, but the model documentation reports no test molecules above 89 heavy atoms. A whole-graph charge vector would be a large size extrapolation, so this is not the route to a simulation-ready model.
2. **Use Amber14SB templates for the entire graph.** The Amber force field contains templates for the standard insulin amino-acid residue names, but the B1 glucoside, B29 carbamate-linker, and macrocycle are not ordinary amino-acid templates. The two modified attachment residues also cannot be treated as unchanged PHE and LYS templates.
3. **Inventory Amber14SB standard-residue templates plus an AmberTools custom-polymer route for modified residues and conjugate fragments.** This is selected. Amber's documented modified-residue workflow is a useful precedent, while its generated charges/types still require manual review. OpenFF's documented Sage plus ff14SB Interchange example concerns separate protein and ligand components and warns the combination is experimental; it does not validate a covalently attached conjugate route.

## Inputs and method

- Read the source-derived neutral full-graph JSON, the existing Sage label-coverage JSON, and the pinned environment lock.
- Load the installed OpenMM `amber14-all.xml` include tree. Record every included XML file hash and the set of residue-template names it defines.
- Group graph atoms by insulin chain, residue number, and residue name using their recorded provenance. Distinguish literal residue-name matches from chemically relevant Amber variants: histidine HIE/HID/HIP forms and disulfide cysteine CYX. Record that neutral-graph termini and physiological-pH charge-state choices are not exact Amber terminal-template matches.
- Record the documented AmberTools custom-polymer workflow (Antechamber, prepgen, parmchk2, LEaP) as a candidate method and the existing isolated PaG GAFF2 result as toolchain precedent only. Keep the OpenFF Sage/ff14SB component-combination example labeled experimental and separate-component; do not claim it supports the covalent NNC2215 joins.
- Identify graph edges that cross from an insulin-template atom to a mapped conjugate atom. Report the B1 Phe N to mapped B1 linker carbonyl join and the LysB29 NZ to mapped B29 carbamate carbonyl join by source node and map identifiers.
- Mark PHEB1 and LYSB29 as modified sites even when ordinary PHE/LYS templates exist. Mark the mapped B1 and B29 atoms as custom-conjugate scope.
- Carry forward the NAGL model's published size-domain evidence and current open graph, protonation, stereochemistry, atom-mapping, and charge-assignment gates as provenance, not as computed parameter results.

## Outputs and acceptance

- Add `research/audit_nnc2215_hybrid_feasibility.py` and generate `research/data/nnc2215-hybrid-feasibility.json`.
- Report formula, heavy-atom and bond counts; the 50 insulin residue identities; literal template-name and variant candidate availability; the six disulfide cysteines; modified PHEB1 and LYSB29 sites; both cross-origin conjugate joins; mapped-conjugate atom count; hashes and package versions; and a scoped go/no-go assessment.
- The audit passes its inventory checks only if every graph residue is classified: 16 observed residue names have literal Amber template names, the two neutral-template histidines can be identified with the available HIE/HID/HIP variants from their hydrogen ledger, all six disulfide cysteines are identified as CYX candidates, and exactly the two expected insulin-to-conjugate joins resolve by provenance. Do not claim exact template matching or charge-state compatibility.
- Use status `AMBER_TEMPLATE_CANDIDATES_IDENTIFIED; CUSTOM_CONJUGATE_PARAMETERS_REQUIRED` when those checks pass. This means the base/variant categories and custom regions are located, not that any force-field parameters have been generated or validated.
- Update the research README, compute status, chemistry reconciliation, milestone report, and dashboard with a rerun command and JSON link.

## Scope limits

This inventory does not call OpenMM's system builder, generate partial charges, choose a protonation microstate, create or validate custom bonded/nonbonded parameters, use candidate coordinates, build a solvated system, run MD, estimate binding, or support a clinical inference. Template or variant candidate presence does not prove an exact residue-template match or charge compatibility. PHEB1 hydrogen ownership, the unique B29 coordinate map, roof stereochemistry, local protonation, terminal-state matching, and compatible cross-boundary parameters remain open. AmberTools' custom polymer-residue workflow is a candidate precedent; it is not parameterization evidence for NNC2215. The output must keep the validated conjugate-force-field count at zero.

## Primary method references

- OpenFF NAGL AM1-BCC GNN v1.0 model domain and data summary: <https://docs.openforcefield.org/projects/nagl-models/en/stable/models/openff-gnn-am1bcc-1.0.0/>.
- AmberTools custom polymer-residue workflow: <https://ambermd.org/tutorials/basic/tutorial5/index.php>.
- OpenFF Sage plus ff14SB Interchange example and combination warning: <https://docs.openforcefield.org/en/latest/examples/openforcefield/openff-interchange/protein_ligand/protein_ligand.html>.
