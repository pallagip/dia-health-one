"use strict";

const sourceUrl = "data/source-data.csv";
const resultUrl = "data/reproduction.json";
let latestManifest = null;

const mean = values => values.reduce((total, value) => total + value, 0) / values.length;
const sampleSd = values => {
  if (values.length < 2) return 0;
  const average = mean(values);
  return Math.sqrt(values.reduce((sum, value) => sum + (value - average) ** 2, 0) / (values.length - 1));
};
const fmt = (value, digits = 2) => Number(value).toFixed(digits);

async function loadPublishedSummary() {
  const response = await fetch(resultUrl, { cache: "no-cache" });
  if (!response.ok) throw new Error("Could not load the published source-data summary.");
  const result = await response.json();
  const rows = result.series.NNC2215.map((point, index) => {
    const control = result.series.Insulin_degludec[index];
    return `<tr><td>${point.glucose_mmol_L}</td><td>${fmt(point.mean)} ± ${fmt(point.sd)}</td><td>${fmt(control.mean)}</td></tr>`;
  }).join("");
  document.getElementById("data-rows").innerHTML = rows;
  document.getElementById("fold-3-20").textContent = `${fmt(result.fold_changes[0].mean_within_experiment_ratio)}×`;
  document.getElementById("fold-0-20").innerHTML = `${fmt(result.fold_changes[1].mean_within_experiment_ratio)}<span>×</span>`;
}

async function loadHillFitSummary() {
  const response = await fetch("data/hill-fit-summary.json", { cache: "no-cache" });
  if (!response.ok) throw new Error("Could not load the exploratory Hill-fit record.");
  const record = await response.json();
  const point = record.derived.point_estimate;
  const bootstrap = record.derived.cluster_bootstrap;
  const ranges = bootstrap.percentile_ranges_not_formal_confidence_intervals;
  const bounds = bootstrap.parameter_bound_diagnostics;

  document.getElementById("hill-fit-slope").textContent = point.apparent_hill_coefficient.toFixed(2);
  document.getElementById("hill-fit-slope-range").textContent =
    `Bootstrap percentiles ${ranges.apparent_hill_coefficient.lower_2_5_percent.toFixed(2)}–${ranges.apparent_hill_coefficient.upper_97_5_percent.toFixed(2)}`;
  document.getElementById("hill-fit-ec50").textContent = point.ec50_mmol_L.toFixed(1);
  document.getElementById("hill-fit-ec50-range").textContent = bounds.upper_percentile_is_censored
    ? `Upper tail hits ${bounds.ec50_upper_bound_mmol_L.toLocaleString("en-US")} mmol/L fit ceiling in ${bounds.ec50_upper_bound_hit_count}/${bootstrap.resamples} resamples`
    : `Bootstrap percentiles ${ranges.ec50_mmol_L.lower_2_5_percent.toFixed(1)}–${ranges.ec50_mmol_L.upper_97_5_percent.toFixed(1)} mmol/L`;
  document.getElementById("hill-fit-note").textContent =
    `${record.interpretation.summary} ${record.interpretation.limits.join(" ")}`;
}

function parseCsv(text) {
  const lines = text.trim().split(/\r?\n/);
  const headers = lines.shift().split(",");
  return lines.map(line => {
    const values = line.split(",");
    const row = Object.fromEntries(headers.map((header, index) => [header, values[index]]));
    row.experiment = Number(row.experiment);
    row.glucose_mmol_L = Number(row.glucose_mmol_L);
    row.relative_affinity_percent = Number(row.relative_affinity_percent);
    return row;
  });
}

function summarize(rows, analyte) {
  return [0, 3, 5, 10, 20].map(glucose => {
    const values = rows.filter(row => row.analyte === analyte && row.glucose_mmol_L === glucose).map(row => row.relative_affinity_percent);
    return { glucose_mmol_L: glucose, values, n: values.length, mean: mean(values), sd: sampleSd(values) };
  });
}

function pairedFold(rows, low, high, reported) {
  const experiments = new Map();
  rows.filter(row => row.analyte === "NNC2215" && (row.glucose_mmol_L === low || row.glucose_mmol_L === high)).forEach(row => {
    if (!experiments.has(row.experiment)) experiments.set(row.experiment, {});
    experiments.get(row.experiment)[row.glucose_mmol_L] = row.relative_affinity_percent;
  });
  const ratios = [...experiments.values()].map(pair => pair[high] / pair[low]);
  const average = mean(ratios);
  return {
    from_glucose_mmol_L: low,
    to_glucose_mmol_L: high,
    replicate_ratios: ratios,
    mean_within_experiment_ratio: average,
    sd_within_experiment_ratio: sampleSd(ratios),
    published_reported_fold: reported,
    matches_published_rounding: Number(average.toFixed(1)) === reported
  };
}

async function sha256(text) {
  const digest = await crypto.subtle.digest("SHA-256", new TextEncoder().encode(text));
  return [...new Uint8Array(digest)].map(value => value.toString(16).padStart(2, "0")).join("");
}

function renderOutput(manifest) {
  const folds = manifest.fold_changes;
  const checks = [
    { title: "Glucose response · 3 → 20 mmol/L", fold: folds[0], note: "Mean of within-experiment ratios; paired by source experiment." },
    { title: "Glucose response · 0 → 20 mmol/L", fold: folds[1], note: "Mean of within-experiment ratios; paired by source experiment." }
  ];
  document.getElementById("output-grid").innerHTML = checks.map(item => {
    const fold = item.fold;
    return `<article class="output-card"><div class="output-card-head"><strong>${item.title}</strong><span class="output-badge">${fold.matches_published_rounding ? "ROUNDING MATCH" : "REVIEW"}</span></div><p><code>${fmt(fold.mean_within_experiment_ratio, 3)}× ± ${fmt(fold.sd_within_experiment_ratio, 3)}</code> mean ± sample SD of paired ratios across three experiments; paper reports ${fold.published_reported_fold.toFixed(1)}×. Replicate ratios: ${fold.replicate_ratios.map(value => fmt(value, 2)).join(", ")}×. ${item.note}</p></article>`;
  }).join("");
  document.getElementById("input-checksum").textContent = `${manifest.run.input_sha256.slice(0, 16)}…`;
  document.getElementById("output-meta").textContent = `RUN ${manifest.run.id} · ${manifest.run.created_at_utc} · ${manifest.run.engine} · Input SHA-256 ${manifest.run.input_sha256} · Source: Nature (2024), Figure 3d. No NNC2215–PaG binding MD or patient-response model was run.`;
  document.getElementById("run-output").hidden = false;
  document.getElementById("run-output").scrollIntoView({ behavior: "smooth", block: "nearest" });
}

document.getElementById("run-button").addEventListener("click", async event => {
  const button = event.currentTarget;
  button.disabled = true;
  button.innerHTML = "Calculating… <span>•••</span>";
  try {
    const [response, scriptResponse] = await Promise.all([
      fetch(sourceUrl, { cache: "no-cache" }),
      fetch("research.js", { cache: "no-cache" })
    ]);
    if (!response.ok || !scriptResponse.ok) throw new Error("Analysis inputs could not be loaded.");
    const [csv, scriptSource] = await Promise.all([response.text(), scriptResponse.text()]);
    const rows = parseCsv(csv);
    const folds = [pairedFold(rows, 3, 20, 3.2), pairedFold(rows, 0, 20, 12.5)];
    const createdAt = new Date().toISOString();
    const manifest = {
      schema: "dia-health-one.reproduction.v2",
      run: {
        id: `DHO-NNC2215-FIG3D-${createdAt.replace(/[-:.TZ]/g, "").slice(0, 14)}`,
        created_at_utc: createdAt,
        engine: "Browser JavaScript arithmetic",
        engine_version: "ES2022",
        analysis: "source-data reanalysis; group means, sample SD, paired experiment fold changes",
        input_sha256: await sha256(csv),
        script_sha256: await sha256(scriptSource),
        gpu_simulation: false,
        clinical_use: false
      },
      source: {
        doi: "10.1038/s41586-024-08042-3",
        article: "https://www.nature.com/articles/s41586-024-08042-3",
        source_data: "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-024-08042-3/MediaObjects/41586_2024_8042_MOESM5_ESM.xlsx",
        source_data_sha256: "09a9b6aa4103f9d61af7df03886d59f16d6b09f69f2fe768260548061c7dc81d",
        accessed_utc_date: "2026-10-17",
        sheet: "Figure 3d",
        assay: "hIR-A affinity relative to human insulin; 1.5% HSA",
        source_rows: rows.length,
        replicate_design: "Three independent experiments"
      },
      inputs: rows,
      series: { NNC2215: summarize(rows, "NNC2215"), Insulin_degludec: summarize(rows, "Insulin degludec") },
      fold_changes: folds,
      interpretation: {
        hill_fit: "A separate SciPy shared-slope Hill fit is available in data/hill-fit-summary.json. It is exploratory; its EC50 bootstrap upper tail is censored at the 1000 mmol/L fit ceiling and does not establish a threshold.",
        scope: "Published in-vitro data reanalysis only; no molecular dynamics, binding free-energy, or patient-response model was run."
      }
    };
    latestManifest = manifest;
    renderOutput(manifest);
  } catch (error) {
    document.getElementById("output-grid").innerHTML = `<article class="output-card"><div class="output-card-head"><strong>Source analysis could not run</strong><span class="output-badge">INPUT ERROR</span></div><p>${error.message} Please use a local web server or the hosted research page so the source CSV can be loaded.</p></article>`;
    document.getElementById("run-output").hidden = false;
  } finally {
    button.disabled = false;
    button.innerHTML = "Run again <span>↻</span>";
  }
});

document.getElementById("manifest-button").addEventListener("click", () => {
  if (!latestManifest) return;
  const blob = new Blob([JSON.stringify(latestManifest, null, 2) + "\n"], { type: "application/json" });
  const url = URL.createObjectURL(blob);
  const anchor = document.createElement("a");
  anchor.href = url;
  anchor.download = `${latestManifest.run.id.toLowerCase()}.json`;
  anchor.click();
  URL.revokeObjectURL(url);
});

loadPublishedSummary().catch(error => {
  document.getElementById("data-rows").innerHTML = `<tr><td colspan="3">${error.message} The CSV and analysis script remain available below.</td></tr>`;
});

loadHillFitSummary().catch(error => {
  document.getElementById("hill-fit-note").textContent = `${error.message} Use the downloadable fit record and script for the full result.`;
});
