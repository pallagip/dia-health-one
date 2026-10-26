# NNC2215 Full Conjugate Coordinate Mapping Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use `superpowers:executing-plans` to implement this plan task-by-task in the current isolated worktree. The report is the scientific audit artifact; do not add or run a software test suite.

**Goal:** Map all 507 heavy atoms in the source-derived insulin–NNC2215 graph to both corrected model-built coordinate candidates while retaining all 16 B29 assignments per state.

**Architecture:** A standard-library Python audit reads the full graph, candidate PDBs, and existing B1/B29 mapping audits. It maps insulin atoms by exact PDB identity, composes each existing conjugate mapping with that protein map, validates coverage and hashes, and writes a deterministic JSON ledger with coordinates and geometry observations. Documentation and the dashboard expose the bounded result without describing it as a parameterized model.

**Tech Stack:** Python 3 standard library, JSON, fixed-column PDB parsing, SHA-256, static HTML.

---

## Files

- Create `research/audit_nnc2215_full_coordinate_mapping.py`: input integrity, PDB atom parsing, full graph mapping, geometry observations, fail-closed report generation.
- Create `research/data/nnc2215-full-coordinate-mapping.json`: mapping ledger for both coordinate candidates and all exact full-graph mappings.
- Modify `research/README.md`, `research/COMPUTE_STATUS.md`, `research/NNC2215_CHEMISTRY_RECONCILIATION.md`, `research/MILESTONE_REPORT.md`, and `research/index.html`: publish report link, rerun command, result counts, and scope limits.

## Task 1: Implement input integrity and PDB parsing

**Files:** Create `research/audit_nnc2215_full_coordinate_mapping.py`.

- [x] **Step 1: Define audit paths and stable atom records.** Use repository-relative paths for the full graph JSON, both candidate PDBs, reconciliation record, B1/B29 topology audits, boundary audit, and output JSON. Add a frozen record for serial, record type, chain, residue number, insertion code, residue name, atom name, altloc, element, and XYZ.

```python
@dataclass(frozen=True)
class PdbAtom:
    serial: int
    record_type: str
    chain: str
    residue_number: int
    insertion_code: str
    residue_name: str
    atom_name: str
    altloc: str
    element: str
    xyz: tuple[float, float, float]
```

- [x] **Step 2: Add deterministic file and PDB readers.** Implement `sha256_file(path)`, `load_json(path)`, and `parse_pdb(path)`. Parse `ATOM`, `HETATM`, `CONECT`, and `LINK`; reject duplicate serials, absent/invalid coordinates, missing elements, malformed atom records, and unresolved connectivity endpoints. Resolve LINK records by full chain/residue/insertion/altloc/name identity and return undirected adjacency edges alongside the atom table. Keep all records in the parser and expose heavy atoms by excluding H and D.
- [x] **Step 3: Verify upstream provenance before mapping.** Check candidate PDB hashes against both the reconciliation record and the B1/B29 topology state records. Check the current full graph and audit inputs against their recorded digests. Snapshot candidate hashes before parsing and compare again after the audit.
- [x] **Step 4: Run the audit's input preflight.** Run `python3 research/audit_nnc2215_full_coordinate_mapping.py --preflight` after implementing the CLI preflight. Expected output: both candidate hashes match; source graph has 507 heavy atoms; upstream B1/B29 mappings report 1 and 16 exact embeddings per state.

## Task 2: Compose the complete graph-to-coordinate mappings

**Files:** Modify `research/audit_nnc2215_full_coordinate_mapping.py`.

- [x] **Step 1: Map the insulin-template atoms by exact identity.** Index heavy PDB atoms by `(chain, residue_number, insertion_code, residue_name, atom_name, element)`. For each of the 398 full-graph atoms with `provenance.origin == "insulin_chain_template"`, require exactly one candidate record with the same identity.
- [x] **Step 2: Adapt the B1 and B29 upstream maps.** Read B1 serials from `candidate_states[].all_exact_mappings[].map_id_to_target[str(map_id)].target_atom.serial`. Read B29 serials from `states[].complete_product_exact_graph_comparison.exact_mappings[].map_id_to_pdb_serial[str(map_id)]`. Normalize graph provenance's integer map IDs with `str(map_id)` for JSON lookups; preserve each upstream mapping ID and all map-to-serial assignments. Reject duplicate states, mapping IDs, and non-exact map-key sets before composing coordinates.
- [x] **Step 3: Merge appendage atoms onto their full-graph nodes.** Index `mapped_source_fragment` graph nodes by `(component, int(source_map_id))`. Map B1 map 1 to the existing protein PheB1 N node and B29 map 92 to the existing protein LysB29 NZ node; verify the upstream serials equal those protein identities. Keep B29 map 84 as the linker carbamate carbonyl carbon and map it to its own conjugate node. Map all other source maps to their corresponding graph nodes without creating endpoint duplicates.
- [x] **Step 4: Enumerate full assignments without coordinate ranking.** For each state, combine its one B1 mapping with each B29 mapping. Require exactly 16 combined mappings, each covering all 507 graph node IDs exactly once and all 507 distinct candidate heavy-atom serials. Never sort or filter mappings by XYZ or bond distance.
- [x] **Step 5: Validate the candidate atom inventory.** Require the graph's 507 heavy-atom element inventory to equal the candidate heavy-atom inventory after excluding only the open state's separate 12-heavy-atom BGL residue. Verify neither CONECT nor LINK has an edge from BGL to the mapped 507-atom graph, list BGL records separately, and list any other unmatched or extra records.

```python
def validate_complete_assignment(
    graph_node_ids: set[str],
    node_to_atom: dict[str, PdbAtom],
) -> None:
    if set(node_to_atom) != graph_node_ids:
        raise ValueError("Full graph node coverage is incomplete")
    serials = [atom.serial for atom in node_to_atom.values()]
    if len(serials) != 507 or len(set(serials)) != 507:
        raise ValueError("Full assignment must use 507 unique heavy-atom serials")
```

## Task 3: Record coordinate and geometry evidence

**Files:** Modify `research/audit_nnc2215_full_coordinate_mapping.py`.

- [x] **Step 1: Serialize complete per-mapping atom records.** Store shared graph provenance and source-bond ledgers, each state's candidate atom identity/XYZ table, and one graph-index-to-PDB-serial vector per mapping. These lossless joins reconstruct each complete atom record without repeating immutable PDB data. Record exact coverage and uniqueness counts per assignment.
- [x] **Step 2: Compute source-edge distances after mapping.** For each of the 527 full-graph heavy-atom edges, store its graph endpoints, bond type/order, and screen interval once in `source_bonds`; retain a 527-value distance vector per mapping aligned by source-bond list order. Apply the documented C/N/O heuristic intervals and list out-of-range edge IDs. Mark all nine sulfur-containing edges (six C–S and three S–S) as unscreened and retain their distances without a threshold. Keep every mapping and every edge even when outside a screen.
- [x] **Step 3: Retain known hydrogen-record issues as unresolved context.** Include the PheB1 H discrepancy and existing LYX hydrogen-neighbor differences from the attachment-boundary audit. Do not map PDB hydrogens, repair records, or infer protonation.
- [x] **Step 4: Implement fail-closed report behavior.** On stale hashes, missing or duplicate identities, changed mapping counts, endpoint disagreement, invalid BGL exclusion, incomplete graph coverage, or non-finite coordinates, emit a `failed` or `partial` report with diagnostics and no complete assignment claim. Include a deterministic assessment that keeps force-field, whole-system charge, OpenMM, and MD readiness false.

## Task 4: Publish the bounded result in project documentation

**Files:** Modify the five documentation/dashboard files listed above.

- [x] **Step 1: Add the README research entry.** Link the JSON and audit script; state 507 graph atoms, 398 protein atoms, 109 conjugate atoms, 1 B1 mapping, 16 B29 mappings per state, and the separate open BGL exclusion. Include `python3 research/audit_nnc2215_full_coordinate_mapping.py` as the rerun command.
- [x] **Step 2: Update compute status and chemistry reconciliation.** Record mapping coverage, hash checks, geometry-screen limits, B29 ambiguity, model-built coordinate caveat, and unresolved PheB1 hydrogens, microstates, parameters, and force-field readiness.
- [x] **Step 3: Add the milestone report entry and website card.** Link the machine-readable report from `research/MILESTONE_REPORT.md` and `research/index.html`. State that this is a coordinate-coverage ledger and that it does not establish bond orders, a validated force field, OpenMM system, MD, binding, efficacy, or clinical readiness.
- [x] **Step 4: Inspect every new link and run `git diff --check`.** Expected result: all local artifact paths exist, and `git diff --check` emits no whitespace errors.

## Task 5: Reproduce, review, merge, and deliver the milestone

**Files:** Audit outputs, documentation, and branch history.

- [x] **Step 1: Generate the report in the repository environment.** Run `python3 research/audit_nnc2215_full_coordinate_mapping.py`. Expected assessment: both states pass input integrity; each has 16 full mappings; each mapping covers 507 unique graph nodes and candidate heavy serials; only the open chain D, residue 1 BGL is excluded as a separate 12-heavy-atom component; all source-edge distances are retained.
- [x] **Step 2: Repeat the audit and compare output hashes.** Run the same command a second time and verify the JSON is byte-identical, the candidate PDB hashes are unchanged, and no existing source graph, PDB, or upstream audit file changed.
- [x] **Step 3: Review the diff and commit the implementation.** Run `git diff --check`, inspect the report and each modified page, and commit the complete milestone with a focused message.
- [ ] **Step 4: Open and merge a PR to `main` under the standing authorization.** Confirm repository checks succeed and the merged tree includes the script, report, and documentation.
- [ ] **Step 5: Verify production publication.** Confirm Cloudflare Pages and GitHub Pages deploy the merge commit. Fetch the live research page and report JSON; verify the JSON bytes match the merged artifact and every linked artifact resolves.
- [ ] **Step 6: Send Patrick the milestone update.** Reply in the existing project thread with the merged PR, live page/report links, mapping counts, geometry-screen scope, and explicit unresolved model gates. Verify the sent reply.

## Completion boundary

This milestone establishes a source-linked heavy-atom coordinate mapping ledger for two model-built candidates and retains all B29 mapping ambiguity. It does not select an authoritative B29 assignment, resolve bond labels or hydrogen ownership, assign a validated physiological microstate, parameterize the complete conjugate, build a complete OpenMM system, or produce target MD or affinity evidence.
