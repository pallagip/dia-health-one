# Dia Health One research reproducibility demo

This directory contains a small, inspectable reanalysis of source values from the NNC2215 glucose-response study, plus a dashboard that presents the published PaG findings separately. It is a computational research demonstration, not a clinical model.

## Run the analysis

From the repository root:

```sh
python3 research/run_reproduction.py
python3 research/reproduce_pag_clamp.py
```

The script uses only the Python standard library. It reads `data/source-data.csv`, calculates the mean and sample standard deviation at each glucose concentration, calculates paired fold changes within each independent experiment, and writes:

- `data/reproduction.json` — machine-readable result and input SHA-256.
- `assets/glucose-response.svg` — original vector chart from the same values.

The JSON result includes the exact source rows used, with experiment IDs and workbook cell references, so the result can be reviewed without reconstructing the input from the chart.

The second script processes the published group-level summaries for the 2026 pig clamp study. It writes `data/pag-clamp-reproduction.json` with the source CSV checksum and descriptive treatment-to-saline ratios. These ratios divide published group means; they are not paired animal-level effects. The study's individual animal measurements and hypothesis tests are not reanalysed by this script.

The dashboard also re-runs the paired calculations in the browser from the CSV and offers a timestamped JSON run manifest. Serve the repository over HTTP to use the browser run button; browsers do not allow `fetch()` of a sibling CSV from a `file://` page.

## Source and extraction

- Høeg-Jensen et al., “Glucose-sensitive insulin with attenuation of hypoglycaemia,” *Nature* (2024), [article and methods](https://www.nature.com/articles/s41586-024-08042-3), DOI `10.1038/s41586-024-08042-3`.
- The publisher's [Source Data workbook](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-024-08042-3/MediaObjects/41586_2024_8042_MOESM5_ESM.xlsx), worksheet `Figure 3d`.
- Workbook downloaded 2026-10-17 UTC; SHA-256 `09a9b6aa4103f9d61af7df03886d59f16d6b09f69f2fe768260548061c7dc81d`.
- Source assay label: hIR-A affinity relative to human insulin, with 1.5% HSA. Three independent experiments, five D-glucose concentrations (0, 3, 5, 10, 20 mmol/L), and two series (NNC2215 and insulin degludec).
- `source-data.csv` contains the workbook's numeric values, not values digitised from a plot. Each row records the source worksheet and source cell. The workbook itself is linked above and is not copied into this repository.

The analysis reports the mean and sample SD across the three experiment-level values at each concentration. To compare with the paper's reported fold changes, it first divides the high-glucose value by the corresponding low-glucose value within each experiment, then reports the mean and sample SD across those three paired ratios. This gives 3.1576× for 3→20 mmol/L and 12.4720× for 0→20 mmol/L, which round to the reported 3.2× and 12.5×. Dividing the group means instead gives different values (3.1379× and 12.2725×); both calculations are preserved in the output to make the analysis choice visible.

To verify that the checked-in CSV still matches the publisher workbook, download the workbook and run:

```sh
python3 research/verify_source_extract.py /path/to/41586_2024_8042_MOESM5_ESM.xlsx
```

The verifier checks the workbook checksum, worksheet name, source-cell references, and all 30 numeric values.

## PaG evidence boundary

Sturis et al., [“Interaction between Paracetamol Glucuronide and a Glucose-Sensitive Insulin Molecule (NNC2215): Increased Risk of Hypoglycemia”](https://pubs.acs.org/doi/10.1021/acsptsci.6c00069), *ACS Pharmacology & Translational Science* (2026), DOI `10.1021/acsptsci.6c00069`, is represented on the dashboard through explicitly labelled published summary values. The reported Kd means are 1.91 mM for D-glucose and 0.17 mM for PaG; their ratio (about 11.2) is calculated from those means. The paper reports 0.9% and 5.4% hIR-A relative affinity at zero glucose without PaG and with 70 μM PaG. These rounded summaries are not independently reanalysed raw values.

The PaG study describes a 500 ns Desmond/OPLS4 simulation from an in-house crystal structure of a related macrocycle, but the article does not supply those starting coordinates. Therefore this demo does not claim an exact molecular-dynamics rerun. PDB 6PXV and 2VJZ are reference structures used in the 2024 paper's modelling; they are not solved atomic structures of every NNC2215 switch state.

The pig-clamp panel reports the article's GIR AUC from 300 to 420 minutes at a 3 mmol/L glucose clamp. Means, SEM, sample sizes, and p-values are transcribed from the main article's Table 1; ratios to saline are derived from the reported means. `data/pag-model-parameters.csv` transcribes the minimal-model values in Supporting Information Table S2. The equations and these selected parameter estimates are not enough to recreate the authors' full PK/PD fit: individual-animal time series and some model inputs are not published with the article.

## Interpretation limits

- No Hill curve or 6 mmol/L threshold is fitted or inferred.
- No molecular dynamics, docking, binding free-energy calculation, or patient-response model is run.
- GPU acceleration is not useful for this small source-data reanalysis. Any later GPU simulation is gated on published structure availability, validated parameters for the conjugate and ligands, and a verified compute allocation.
- All source data, derived calculations, model-built structures, and future simulation outputs must remain separately labelled.
