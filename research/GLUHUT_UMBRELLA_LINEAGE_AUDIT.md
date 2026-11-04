# GluHUT umbrella-trajectory lineage audit

**Date:** 4 November 2026  
**Status:** route/window naming reconciled; full histogram lineage and regeneration failed  
**Scope:** public `US_Trajectories.zip` versus deposited WHAM inputs

## Question

The deposited umbrella archive contains path-labelled trajectory files, while the WHAM package contains the reaction-coordinate time series used for the published PMFs. This audit asks whether the trajectory archive is complete enough to regenerate those histogram inputs and whether every nominal route and window represents distinct data.

It does not calculate a new PMF or change the already reproduced deposited WHAM results.

## Verified sources

### Umbrella trajectory archive

- file: `US_Trajectories.zip`;
- size: 123,756,824 bytes;
- MD5: `c460e6a1ee97ac531fc284d8f07e381e`;
- SHA-256: `1e68f57ef15c3c010a767ab3c8181b5eb8df106d1ee96066f6afa4d54cb8f517`;
- ZIP CRCs: verified.

### WHAM package

- file: `WHAM_PMF.zip`;
- size: 41,560,830 bytes;
- MD5: `dbc6723b8adcb85deb8230f63b6a51da`;
- SHA-256: `e0e7b7c04326b8fe86cf71e2193315f34fb865168cf11117b3caf9253119c6eb`.

## Archive layout

The umbrella archive contains:

- 714 Scripps BINPOS trajectory files;
- 14 raw route directories;
- 51 path-labelled windows per route;
- path labels from 0.0 to 15.0 Å in 0.3 Å increments;
- 11 topology files, three restraint inputs and two protocol files;
- no source force-field generation package;
- no authoritative route-to-WHAM lineage manifest.

The 14 raw route directories are:

- glucose runs 1–3;
- fructopyranose runs 1–3;
- fructofuranose runs 1–3;
- galactose path-B runs 1–3;
- galactose path-A runs 4–5.

The deposited WHAM package contains 13 route records. The count difference is explained by an exact complete-route duplicate in the trajectory archive.

## Complete galactose path-A duplication

Galactose run4(a) and run5(a) are exact copies across all 51 windows.

They have:

- the same ordered set of window centres;
- the same file SHA-256 at every corresponding window;
- the same complete-route content signature;
- one mapped WHAM route: `galactose/US-routeA`.

Therefore:

- archive route directories: 14;
- unique complete-route content signatures: 13;
- deposited WHAM routes: 13.

Galactose path-A has two directory copies but only one unique trajectory route. The two copies must not be counted as independent simulations.

## Duplicate windows inside galactose path-B routes

Two path-B routes contain byte-identical files under different window labels:

- galactose run1(b): 0.0 and 0.3 Å files are identical;
- galactose run2(b): 0.0 and 0.3 Å files are identical.

The path names therefore indicate 51 centres, but those two routes each contain only 50 unique trajectory files.

This is a file-provenance failure. It does not establish which file is wrong or whether the corresponding full-resolution trajectory was also duplicated.

## WHAM-only 2.25 Å galactose windows

The trajectory archive path labels for each galactose path-B run contain 51 centres:

`0.0, 0.3, 0.6, ... 15.0 Å`

Each corresponding WHAM route contains 52 centres. The additional centre is:

`2.25 Å`

This occurs in:

- galactose run1(b) → `galactose/US-1`;
- galactose run2(b) → `galactose/US-2`;
- galactose run3(b) → `galactose/US-3`.

No 2.25 Å trajectory file exists in the umbrella archive. Conversely, no trajectory-only centre exists in those routes.

Therefore route and centre mapping is not one-to-one. The 2.25 Å WHAM inputs came from data that are not represented by a named trajectory file in the public umbrella archive.

## Fructofuranose frame-count failures

The filenames contain `100-frames`, but this is not a reliable record count.

### Fructofuranose run 2

- 50 windows contain 10 frames each;
- the 12.6 Å window contains only **one frame**;
- unique recorded frame counts: `[1, 10]`.

The one-frame trajectory is:

`US_Trajectories/fructofuranose/run 2/prod-12-6-100-frames.binpos`

The same 12.6 Å centre is the known WHAM anomaly containing only 205 reaction-coordinate samples rather than approximately 15,001. This localizes the data-provenance problem to the same route and centre.

It does not prove that the one-frame visualization file directly generated the 205-sample WHAM input.

### Fructofuranose run 3

All 51 files contain **10 frames**, not 100.

### Other routes

The remaining route files contain 100 recorded frames per window.

## Visualization versus WHAM sampling

For ordinary routes:

- visualization trajectory: 100 frames per window;
- deposited WHAM input: approximately 15,001 samples per window;
- ratio: approximately 150:1.

For fructofuranose runs 2 and 3, the gap is substantially larger because the public visualization files contain 10 frames, and one run-2 file contains only one.

The public umbrella archive therefore does not contain the full reaction-coordinate time series used to build the WHAM histograms.

Missing provenance includes:

- the full-resolution umbrella trajectories;
- the exact frame-selection/downsampling script;
- the time mapping of the retained visualization frames;
- the source of the galactose 2.25 Å windows;
- corrected or complete fructofuranose run-2/run-3 files;
- an authoritative route-to-WHAM lineage manifest.

## What remains reproducible

The following remains reproducible:

- checksum verification of both archives;
- deterministic parsing of all deposited WHAM inputs;
- reproduction of the deposited WHAM PMFs and qualitative sugar ranking;
- exact route/window/file inventory of the visualization archive;
- duplicate-file and duplicate-route detection;
- identification of mismatched centres and record counts.

The following is **not** reproducible from the umbrella archive:

- generation of the deposited histogram time series;
- independent recomputation of those histograms from full trajectories;
- verification of equilibration discard and statistical inefficiency from the public trajectory files;
- prospective regeneration of the original umbrella protocol.

## Decision

1. Continue to treat the deposited WHAM inputs as an auditable published data product.
2. Do not describe the umbrella visualization archive as sufficient to regenerate the PMFs.
3. Do not count galactose run4(a) and run5(a) as independent routes.
4. Do not treat the 0.0 and 0.3 Å files in galactose run1(b)/run2(b) as independent window trajectories.
5. Do not ignore the missing 2.25 Å galactose path-B lineage.
6. Do not infer direct lineage between the one-frame fructofuranose 12.6 Å file and the 205-sample WHAM input without author confirmation.
7. Keep new PMF calculation, parent glucose-versus-PaG calibration and MC-Q1 ranking blocked until the missing full-resolution data and lineage are obtained or a new benchmark is generated prospectively.

## Author clarification requested

The requested package now includes:

- the original full-resolution umbrella trajectories;
- the exact downsampling or frame-selection script;
- a window-by-window route-to-WHAM manifest;
- the missing galactose 2.25 Å trajectory data;
- corrected galactose 0.0/0.3 files if the duplicate is accidental;
- an explanation for duplicate path-A run4/run5 directories;
- complete fructofuranose run2/run3 visualization or production files;
- the source reaction-coordinate time series and equilibration/discard rules.

## Reproducibility

- Archive inventory: `research/data/gluhut-us-archive-inventory.json`.
- Route lineage: `research/data/gluhut-us-route-lineage.json`.
- Compact anomaly ledger: `research/data/gluhut-us-lineage-anomalies.json`.
- Inventory tool: `research/audit_gluhut_us_archive.py`.
- Route reconciler: `research/summarize_gluhut_us_route_lineage.py`.
- Anomaly summarizer: `research/summarize_gluhut_us_lineage_anomalies.py`.
- Public evidence PRs: #52, #55 and #56.

## Evidence boundary

This is a provenance audit of public umbrella trajectory and WHAM files. It is not a new PMF, a binding affinity, an NNC2215 result, a PaG result, an MC-Q1 result, efficacy, safety or a working-drug result.

## Delivery status

**Delivered:** Gmail confirmed the umbrella-specific provenance addendum was sent on 4 November 2026 to Wayne Powell (CC: Anthony Davis) in the existing GluHUT reproducibility thread. The earlier 4 November Gmail rate limit was transient and the saved draft was successfully sent without editing it.

- Gmail sent-message ID: `1a1a8f2999166f97`.
- Gmail thread ID: `1a1a7ce37fa8d5f5`.
- Gmail-reported message timestamp: `2026-11-04T23:43:53` (provider response has no timezone suffix).
- Delivery evidence: Gmail returned message ID, matching thread ID and `SENT` label, and the sent message was independently read back.
- Original draft ID `r-943222025130158032` is no longer pending; the `draft_created` status field records its historical creation, not an unsent message.

This is evidence of a successful **Gmail send**, not of recipient delivery confirmation, reading, author acknowledgement, or a scientific answer. All source-lineage and parameterization gaps remain open until independently resolved.

