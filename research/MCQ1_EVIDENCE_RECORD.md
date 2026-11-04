# MC-Q1 evidence record

**Updated:** 4 November 2026  
**Evidence level:** E2 — defined and computationally stress-tested  
**Experimental status:** no compound in this record has been synthesized or measured by Dia Health One

## Problem being addressed

NNC2215 established that an insulin molecule can change activity with glucose. It also exposed a decisive liability: paracetamol glucuronide (PaG) binds the NNC2215 glucose-sensing system more strongly than glucose and can increase activity at low glucose. The published point estimates are approximately:

- glucose `Kd = 1.91 ± 0.10 mM`;
- PaG `Kd = 0.17 ± 0.01 mM`;
- point-estimate PaG affinity preference: 11.24-fold.

The rescue question is therefore not merely “can the receptor bind glucose?” It is “can a receptor retain useful glucose recognition while rejecting PaG and a broader glucuronide panel?”

## Primary hypothesis

**MC-Q1** replaces one model-local roof ethyl substituent with methoxy at the site most exposed to the modeled PaG aromatic appendage.

The current mechanism label is:

> contact-site polarity, electronics, hydration, solvation or rotamer hypothesis

MC-Q1 is **not** supported as a simple steric blocker. In the denser fixed-geometry screen, the idealized ether oxygen retained almost the same PaG-contact coverage as the parent proximal carbon.

## Causal control set

| Compound | Role | Question answered |
|---|---|---|
| **MC-P0** | Parent reference | What are glucose and PaG binding in the same assay block? |
| **MC-Q1** | Contact-site methoxy | Does the local edit improve glucose-over-PaG selectivity? |
| **MC-Q1D** | Distal-A methoxy positional control | Is any effect specific to the proposed contact site? |
| **MC-P3** | Contact-site hydroxymethyl constitutional control | Does topology matter beyond adding polarity? |
| **MC-P1** | First expansion control; methyl size-down | Does smaller nonpolar substitution preserve glucose better than methoxy? |

MC-Q1 and MC-Q1D have the same formula and scalar descriptors but different positions relative to the anchored linker. MC-Q1 and MC-P3 are constitutional isomers.

## Geometry stress tests

The PaG aromatic appendage was sampled on the same three-torsion rigid grid at two resolutions.

### Denser 15-degree grid

- total placements: 13,824;
- contact-site proximal carbon placements below 4 Å: **5,376**;
- distal-A placements below 4 Å: **0**;
- distal-B placements below 4 Å: **0**.

The contact-versus-distal distinction therefore survived an eightfold increase in grid size.

### Methoxy rotor clearance against modeled glucose

- contact MC-Q1: **72/72** idealized rotamers free of an immediate van der Waals clash;
- distal-A MC-Q1D: **72/72**;
- distal-B control: **52/72**.

Distal A is therefore the defined MC-Q1D position.

These are fixed-coordinate clearance and site-selection results, not binding predictions or conformational probabilities.

## Rejected design direction

A generic C6-side donor, charge or polarity edit is not advanced. In the archived lowest-overlap PaG placement, the PaG glucuronate carboxylate nearly overlays the modeled glucose C6 hydroxymethyl region and faces the same modeled urea donors. Such an edit is likely to weaken desired glucose recognition together with PaG binding.

An apparent nearby carbon in the author-deposited model was also rejected as a design anchor because source-chemistry reconciliation removes it as an excess methylene.

## Negative public prior for methoxy

The public macrocycle patent reports a related full-roof comparison:

- mixed trimethoxy/triethyl receptor: glucose `Ka = 1,310 ± 33 M⁻¹`;
- mixed trimethyl/triethyl receptor: glucose `Ka = 14,926 ± 1,566 M⁻¹`.

The point-estimate difference is 11.39-fold, approximately 1.44 kcal/mol at 298.15 K, in favour of the methyl-containing receptor.

This does not directly test the single-site MC-Q1 edit and contains no PaG result. It does mean that methoxy cannot be assumed to preserve glucose affinity. MC-P1 is therefore the first expansion control.

## Quantitative direct-binding gate

MC-Q1 advances from hypothesis to a measured hit only if authenticated material satisfies all of the following:

1. `Kd(PaG) / Kd(glucose) >= 5`;
2. the lower 95% confidence bound is above 1;
3. glucose affinity is no worse than fivefold relative to same-run MC-P0;
4. an orthogonal direct-binding method confirms the selectivity direction;
5. MC-Q1 outperforms MC-Q1D, supporting a local-site mechanism;
6. identity, positional assignment, solubility and aggregation controls are acceptable.

Stretch target: `Kd(PaG) / Kd(glucose) >= 10`.

No insulin conjugation or animal work is justified before this free-binder gate and a subsequent non-insulin reversible reporter gate pass.

## What would falsify MC-Q1

Stop or redesign if:

- MC-Q1 and MC-Q1D behave similarly;
- MC-Q1 loses more than fivefold glucose affinity;
- no variant reaches `Kd(PaG) / Kd(glucose) > 1`;
- direct binding fails orthogonal confirmation;
- material cannot be synthesized and positionally assigned as one defined species;
- rights or provider terms do not support a developable asset.

A failed result must be retained as a program decision rather than hidden or replaced by a larger uncalibrated library.

## Current confidence statement

There is moderate confidence that the four-compound set is a useful experiment. There is essentially no direct confidence yet that MC-Q1 is a successful binder, switch or drug. That confidence requires calibrated host–guest computation and authenticated physical measurements.

## Primary references

- Hoeg-Jensen et al., glucose-sensitive insulin NNC2215, *Nature* 2024, DOI: `10.1038/s41586-024-08042-3`.
- Sturis et al., PaG interaction with NNC2215, *ACS Pharmacology & Translational Science* 2026, DOI: `10.1021/acsptsci.6c00069`.
- Eades et al., GluHUT explicit-solvent free-energy simulations, *Chemical Science* 2026, DOI: `10.1039/D6SC00738D`.
- Macrocycle patent family: `WO2018167503` / `US10800747`.
