# Dia Health One public-repository agent contract

This repository contains the **public evidence record and website** for the Dia Health One glucose-responsive-insulin program.

The canonical private technical and decision archive is `w1ne/dia-health-one-research-archive`. Scientific work must remain synchronized with that archive. Public agents must not create a second, conflicting scientific state.

## Read before editing

1. Read `research/status/status.json`.
2. Read:
   - `research/AWS1_GATE_DECISION.md`;
   - `research/AWS1_NO_AWS_RESULTS.md`;
   - `research/CONFIDENCE_LEDGER.md`;
   - `research/MCQ1_EVIDENCE_RECORD.md`.
3. Search open pull requests before starting work.
4. Check whether the change is:
   - website/presentation only;
   - public evidence synchronization;
   - a new scientific result.
5. For scientific changes, verify that a canonical private issue/PR exists and reference it in the public PR.

## Repository roles

- `w1ne/dia-health-one-research-archive` — canonical private science, procurement, provider communication and decision log.
- `pallagip/dia-health-one` — public bounded evidence, machine-readable status and website.

Do not publish confidential provider correspondence, private quotations, licence discussions, unpublished structures not approved for release, or internal procurement data.

## Current public scientific state

- AWS-1: **PASS for deposited WHAM-analysis reproduction**, with mandatory qualifications.
- Fructofuranose US-2 route robustness: **FAIL**.
- Raw trajectory-to-WHAM lineage: **FAIL**.
- New SMD/umbrella generation from deposited inputs: not authorized.
- AWS-2 source-controlled parent-calibration readiness: active in the private archive.
- AWS spend: `$0`.
- MC-Q1: evidence level **E2**, unsynthesized and unmeasured.
- No Dia Health One compound is a measured hit, working switch, insulin candidate or drug.

The first physical E4 gate remains synthesis and two-method direct binding of MC-P0, MC-Q1, MC-Q1D and MC-P3.

## Latest bounded science records

- Public PR #60 / `cb083036fb5410f94412781531de89db94b2ba5e` — bounded AWS-1 PASS.
- Public PR #38 / `75be8df449ccda9c1ae3c1d84581fe2018eea880` — bounded fructofuranose restraint alternatives.
- Private PR #85 / `7d8e51e6828cbd4c760414671416f117f83684c4` — canonical multi-agent coordination contract and MC-P0 graph outputs.

The public `main` branch may include later design or website-only commits. Do not interpret a newer website commit as a new scientific result.

## Public-PR requirements

Every scientific public PR must state:

- the canonical private issue or PR being mirrored;
- what result is public;
- the exact pass/fail/block decision;
- what remains unknown;
- whether machine-readable status changes;
- why the change does not overstate a hit, efficacy, safety or drug claim.

A website-only PR must not silently alter scientific values or evidence levels.

## Evidence boundaries

Do not:

- describe MC-Q1 as a hit, lead, working switch or drug;
- authorize candidate ranking before parent calibration;
- present Boltz protein–ligand affinity as macrocycle–sugar affinity;
- infer a Hill curve from concentration-labelled structure predictions;
- claim that public raw trajectories regenerate the deposited WHAM series;
- hide the failed fructofuranose route or SMD-coordinate gates;
- transfer public truncated-host parameters directly to PaG or MC-Q1;
- equate a free-receptor model with complete NNC2215 experimental binding;
- publish confidential provider or licensing content.

## Coordination before work

Before a scientific edit, search the canonical private repository for:

- active issue #82 — AWS-2 parent-calibration readiness;
- issue #58 — physical MC-Q1 discovery track;
- `AGENTS.md`;
- `research/private-data/agent-coordination-state.json`.

If private access is unavailable, keep the public task presentation-only and state that canonical science synchronization is required before merge.

## Merge discipline

- One bounded question per branch.
- Preserve failed gates and negative results.
- Use deterministic outputs and stable hashes when possible.
- Verify public links and machine-readable JSON before merge.
- Do not claim a workflow passed when no runner or no step executed.
- Do not claim a generated file exists until it is present on `main`.
- Keep `research/status/status.json` consistent with the human-readable records.

## Current public task boundary

Public work may currently:

- improve discoverability and presentation of already merged evidence;
- synchronize bounded conclusions from the private canonical archive;
- validate links, metadata, machine-readable consistency and deployment;
- make evidence boundaries clearer.

Public work may not initiate or rank new drug candidates without the corresponding private canonical gate and evidence package.
