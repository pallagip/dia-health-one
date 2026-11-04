# Dia Health One confidence ledger

**Updated:** 4 November 2026  
**Purpose:** state exactly what is known, what is inferred, what is untested, and which result would increase or decrease confidence in the current glucose-responsive insulin program.

## Bottom line

Dia Health One has **not discovered a working drug or an experimentally confirmed hit**.

The project has established a real therapeutic problem, reproduced key public evidence, defined a compact and falsifiable candidate experiment, and eliminated several weak computational directions. The current primary rescue hypothesis, **MC-Q1**, remains a structure-defined but unsynthesized molecule. Its glucose affinity and PaG selectivity are unknown.

The next confidence-changing result is not another static score. It is either:

1. an independently reproducible host–guest free-energy benchmark followed by a source-matched glucose-versus-PaG calculation; or
2. preferably, authenticated MC-P0/MC-Q1/MC-Q1D/MC-P3 material with direct binding measurements.

## Evidence levels

| Level | Meaning | Required evidence |
|---|---|---|
| **E0 — concept** | An idea or design specification | Written hypothesis only |
| **E1 — defined** | Exact graph, sequence, protocol or model input exists | Reproducible structure/protocol and provenance |
| **E2 — computationally stress-tested** | Result survives prespecified sensitivity and negative-control tests | Reproducible code, outputs, alternative assumptions and falsification controls |
| **E3 — calibrated computation** | Method reproduces an external experimental benchmark before predicting a new case | Independent benchmark, convergence, uncertainty and hold-out prediction |
| **E4 — directly measured** | Authenticated material gives reproducible direct physical data | Identity, purity, raw data, independent runs and orthogonal method |
| **E5 — functionally validated** | Complete switch or insulin construct passes functional biological gates | Human-cell signaling, selectivity, reversibility and PaG/glucuronide resistance |
| **E6 — preclinical candidate** | Manufacturable candidate passes animal PK/PD and safety-oriented gates | Reproducible in-vivo efficacy, PK/PD, formulation and risk package |

A project statement must not use a higher evidence level than its supporting record.

## Current confidence by claim

| Claim | Current level | Confidence assessment | Evidence | Missing result |
|---|---:|---|---|---|
| NNC2215 is glucose-responsive | E5 for the published molecule | High | Published receptor, cellular and animal evidence | Independent replication is desirable but not necessary to establish the public result |
| PaG is a material NNC2215 liability | E5 for the published molecule | High | Direct native-MS binding, receptor assay, pig clamp and PK/PD modeling | None for the conclusion that NNC2215 itself is not viable |
| The archived NNC2215 chemical graph is sufficiently defined for bookkeeping and local audits | E2 | Moderate | Source reconciliation, graph audits and mapping sensitivity analyses | Authoritative atom mapping and complete experimentally validated force field |
| One roof site is more exposed to the modeled PaG aromatic appendage than two distal sites | E2 | Moderate | 1,728- and 13,824-placement rigid-grid audits plus distal controls | Experimental PaG-bound structure or calibrated explicit-solvent free-energy result |
| MC-Q1 is a useful causal experiment | E2 | Moderate | Exact mapped graph; MC-Q1D positional control; MC-P3 constitutional control; rotor and grid stress tests | Synthesis feasibility, material identity and direct binding |
| MC-Q1 preserves glucose affinity | E0 | Low | No measurement; public full-roof methoxy series is a negative prior | Direct glucose binding versus same-run MC-P0 |
| MC-Q1 rejects PaG | E0 | Unknown | No direct or calibrated free-energy result | Direct PaG binding and selectivity ratio |
| MC-Q1 can be synthesized and positionally assigned | E0–E1 | Low to moderate | Plausible precursor concepts and an RFI exist | Written provider acceptance, route, price and analytical assignment plan |
| MC-Q1 can form a useful reversible switch | E0 | Unknown | No reporter construct | Non-insulin reporter dose response, reversibility, HSA and interferent panel |
| MC-Q1 insulin conjugate is potent and safe enough to advance | E0 | Unknown | No conjugate exists | hIR-A/hIR-B, pIR/pAKT, IGF-1R, potency, Emax, reversibility and PaG data |
| DHO-GC1 independent binder concept can produce an owned platform | E0–E1 | Low but strategically important | Binding specification and provider-ready discovery design | Sequence-defined hit, native-glucose binding and glucuronide rejection |

## Current candidate set

### Minimum information package

1. **MC-P0** — parent reference.
2. **MC-Q1** — contact-site methoxy hypothesis.
3. **MC-Q1D** — same methoxy edit at distal roof site A; positional-isomer control.
4. **MC-P3** — contact-site hydroxymethyl; constitutional-isomer control.

### First expansion control

5. **MC-P1** — contact-site methyl size-down control.

The public patent series makes MC-P1 more important than lower-information analogues: a receptor with a full methoxy roof had substantially weaker glucose affinity than a related methyl-roof receptor. This is not direct MC-Q1 evidence, but it removes any assumption that methoxy is intrinsically glucose-preserving.

## What the existing computation supports

The current computational work supports the following bounded claims:

- a source-derived chemical graph and conditional coordinate mapping exist;
- the modeled PaG aromatic appendage can enter one local roof corridor in the rigid placement screen;
- the contact-versus-distal distinction survives an eightfold denser torsion grid;
- distal roof site A is a cleaner MC-Q1D control than distal site B in the idealized glucose-clearance scan;
- a generic C6-side donor or charge edit is poorly selective in the current geometry because PaG and glucose overlap at that end;
- MC-Q1 is **not** supported as a simple steric blocker; any benefit must arise from local polarity, electronics, hydration, solvation or conformational effects;
- a four-compound physical experiment can distinguish the main causal explanations.

It does **not** support a binding affinity, selectivity ratio, Hill coefficient, receptor response, efficacy or safety claim for MC-Q1.

## Confidence-changing gates

### Gate C0 — public evidence and provenance

Pass when all public claims, scripts, generated outputs and limitations are accessible from the public research site, with stable hashes and no result described as private-only.

### Gate C1 — method-domain validity

Pass only when the chosen computational method is designed for a small-molecule host–guest system. Protein–ligand affinity models do not pass this gate for GluHUT–glucose or GluHUT–PaG.

### Gate C2 — independent benchmark reproduction

Before predicting MC-Q1, reproduce the published GluHUT–glucose free-energy workflow or reanalyse its deposited trajectories:

- source-matched host and guest structures;
- GAFF/RESP host treatment and GLYCAM06 carbohydrate treatment, or a justified equivalent;
- explicit water;
- three independent paths or an equivalently justified sampling design;
- convergence and overlap diagnostics;
- uncertainty that includes between-path variation;
- qualitative reproduction of the published glucose-versus-galactose/fructose ranking.

A method that cannot reproduce the public benchmark must not rank new variants.

### Gate C3 — parent glucose-versus-PaG calibration

Apply the validated host–guest workflow to a source-matched parent receptor and both glucose and PaG. The result must be compared with the experimental direction and magnitude available for NNC2215, while clearly separating free-receptor and full-conjugate observables.

Fail if:

- the parent trend has the wrong direction;
- convergence is path-dependent beyond the prespecified tolerance;
- force-field or protonation alternatives reverse the result;
- the exact parent structure cannot be established.

### Gate C4 — prospective computational hold-out

Only after C2 and C3 pass, compare MC-Q1, MC-Q1D, MC-P3 and preferably MC-P1 using the same protocol. Pre-register the ranking rule and do not redesign the compounds after seeing the result without creating a new prospective round.

A computational result can move a candidate to **E3**, not E4.

### Gate C5 — synthesis feasibility

Pass when a provider accepts exact mapped structures and provides:

- route and regioisomer strategy;
- analytical position assignment;
- deliverable mass and purity;
- price, timing and stop points;
- rights and data-ownership terms.

### Gate C6 — direct free-binder experiment

Pass MC-Q1 only if:

- `Kd(PaG) / Kd(glucose) >= 5`;
- the lower 95% confidence bound is above 1;
- glucose affinity is no worse than fivefold versus same-run MC-P0;
- an orthogonal method confirms the selectivity direction;
- MC-Q1 outperforms MC-Q1D, supporting a local-site mechanism;
- identity, positional assignment, solubility and aggregation controls are acceptable.

This is the first gate that can move MC-Q1 to E4.

### Gate C7 — non-insulin reporter

Require a physiologically relevant glucose transition, low-to-high-to-low reversibility, albumin robustness and a full glucuronide/sugar counter-screen before insulin conjugation.

### Gate C8 — insulin function

Require human-cell functional signaling, practical potency and Emax, IR-A/IR-B balance, IGF-1R selectivity, PaG resistance and reversibility before animal work.

## Current stop rules

Stop or redirect Track A if any of the following occurs:

- the source-matched parent cannot be parameterized without unresolved chemistry;
- the benchmark method does not reproduce the published host–guest trend;
- MC-Q1 and MC-Q1D behave similarly, eliminating the proposed local mechanism;
- MC-Q1 loses more than fivefold glucose affinity;
- no variant reaches `Kd(PaG) / Kd(glucose) > 1`;
- rights or provider terms do not permit a developable asset;
- synthesis and positional assignment are economically unrealistic.

If Track A stops, preserve the negative result and move resources to the independent DHO-GC1 campaign rather than expanding the same macrocycle library indefinitely.

## Current status in one sentence

**The program has a credible problem, a disciplined four-compound experiment and moderate confidence in where to test, but essentially no confidence yet that MC-Q1 is a successful binder or drug; that confidence must come from calibrated host–guest computation and authenticated direct-binding data.**
