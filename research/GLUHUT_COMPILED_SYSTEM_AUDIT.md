# Compiled GluHUT system audit

**Date:** 4 November 2026  
**Status:** public benchmark-system audit complete; source parameter-generation provenance remains missing  
**Scope:** four deposited GluHUT–sugar AMBER topologies in `SMD_trajectories.zip`

## Question

Did the deposited glucose, galactose, fructopyranose and fructofuranose simulations use one common compiled GluHUT host, or were sugar-specific host parameters hidden inside the four AMBER topology files?

This matters because a benchmark comparison is interpretable only if the host itself is held constant across the sugar series.

## Verified archive

- file: `SMD_trajectories.zip`;
- size: 37,931,874 bytes;
- MD5: `9f3c6733dc608d2bd37c488c560a5443`;
- SHA-256: `04c58666328764aeb47bac0d9f8607546f9763892cb9a333b108630df70e3943`;
- ZIP CRC: verified.

The corrected content-aware inventory contains:

- four compiled AMBER `parm7` topologies;
- twenty `binpos` trajectories;
- three AMBER NMR-style restraint-control files with `.rst` suffix;
- one SMD control file and one galactose note;
- no standalone coordinate restart file;
- no source RESP, electrostatic-potential, `mol2`, `frcmod`, `lib` or equivalent parameter-generation file.

The `.rst` suffix was initially ambiguous. Inspection of the file contents shows `&rst`, `iat`, `igr1`, `igr2`, `r2`, `rk2` and related restraint fields. These are restraint inputs, not restart coordinates.

## Topology structure

Each of the four compiled systems has:

- 183 atoms;
- three residues;
- a 159-atom synthetic host residue labelled `MOL`;
- a 22-atom sugar residue;
- a two-atom residue labelled `ROH` that balances the sugar charge;
- total compiled system charge effectively zero within floating-point representation.

The sugar residue labels are:

| System | Sugar label | Sugar charge | ROH charge |
|---|---|---:|---:|
| glucose | `0GB` | +0.194 e | −0.194 e |
| galactose | `0LB` | +0.194 e | −0.194 e |
| fructopyranose | `0CB` | +0.194 e | −0.194 e |
| fructofuranose | `0CU` | +0.194 e | −0.194 e |

The non-integer sugar and `ROH` charges are part of the deposited compiled systems and require source-level review before chemical interpretation. This audit records them; it does not validate them.

## Main result

The 159-atom host is compiled identically in all four topologies.

Exact equality was confirmed for:

- atom count and residue label;
- atom names and atomic numbers;
- per-atom partial charges;
- atomic masses;
- AMBER atom types and type indices;
- bonds and their force constants/equilibrium distances;
- angles and their force constants/equilibrium values;
- dihedrals, periodicities, phases and scaling factors;
- excluded-atom lists;
- physical Lennard–Jones A/B coefficients;
- compiled indexing metadata.

For every sugar-versus-reference comparison:

- maximum absolute host charge difference: `0.0 e`;
- maximum absolute host mass difference: `0.0`;
- every physical component hash matched;
- every compiled-indexing component hash matched.

Common physical host fingerprint:

`78e512d17cf1a9db672f22da49a3e59f96bb0af7dae5c059e8a4d1353de46cb1`

Common compiled-indexing fingerprint:

`33acc851a118914732d8a89f21ec2a93808e1cae2948fc024dd7430026cf10ad`

## Interpretation

This materially increases confidence in the public sugar benchmark:

- differences among the deposited sugar simulations were not created by hidden changes to the compiled host charges, bonded terms, exclusions or Lennard–Jones parameters;
- the same compiled 159-atom host can be treated as the common baseline for reproducing the public GluHUT sugar series;
- the previously reproduced sugar ordering is therefore a comparison against one internally consistent compiled host model.

It does **not** establish that the host parameterization is chemically correct. The archive does not contain the source RESP calculation, quantum-chemistry inputs, atom-typing history, original GAFF/GLYCAM files or standalone start coordinates.

## Decision

1. Use this exact compiled host only for reproducing and stress-testing the public GluHUT benchmark.
2. Do not transfer its charges or parameters directly to PaG, MC-Q1, MC-Q1D, MC-P3 or the complete NNC2215 conjugate.
3. Request the missing source RESP/GAFF/GLYCAM generation package from the authors.
4. Independently reconcile the 159-atom compiled host with an exact chemical graph before new ligand calculations.
5. Review the deposited sugar and `ROH` charge construction before treating the guest models as transferable.
6. Keep candidate ranking and PaG parameterization disabled.

## Confidence effect

- Public analysis reproduction: supported.
- Common-host assumption across the four sugar systems: **supported exactly at the compiled topology level**.
- Source parameter derivation: not reproducible from the deposit.
- Chemical correctness of the compiled host: not established.
- PaG calibration: not performed.
- MC-Q1 ranking: not authorized.
- Drug discovery claim: none.

## Reproducibility

- Parser and audit: `research/audit_gluhut_compiled_systems.py`.
- Corrected archive inventory: `research/audit_gluhut_smd_archive.py`.
- Machine-readable result: `research/data/gluhut-compiled-system-audit.json`.
- Corrected machine-readable inventory: `research/data/gluhut-smd-archive-inventory.json`.
- Tests: `tests/test_audit_gluhut_compiled_systems.py` and `tests/test_audit_gluhut_smd_archive.py`.
- Workflow: `.github/workflows/gluhut-compiled-system-audit.yml`.

## Evidence boundary

This is an atom-level audit of deposited compiled AMBER benchmark systems. It is not an affinity calculation, an NNC2215 result, a PaG result, an MC-Q1 result, efficacy, safety or a working-drug result.
