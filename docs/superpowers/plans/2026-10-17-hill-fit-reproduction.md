# Hill Fit Reproduction Implementation Plan

> **For agentic workers:** Implement inline under the user's existing authorization to reproduce the published source data and improve the research dashboard.

**Goal:** Add an auditable exploratory Hill fit to the Nature 2024 NNC2215 receptor-affinity source data and show why its midpoint is not reliably identified.

**Architecture:** A pinned scientific Python environment reads the checked-in source-data CSV, fits a shared apparent Hill slope and EC50 with experiment-specific upper responses and measured zero-glucose baselines, then writes a checksummed JSON record and SVG. The static dashboard embeds the result with its assumptions and uncertainty limits; the existing browser reanalysis remains unchanged.

**Tech Stack:** Python, NumPy, SciPy, Matplotlib, existing CSV/JSON/SVG assets, static HTML/CSS.

---

### Task 1: Add the reproducible Hill fit

**Files:**
- Create: `research/fit_nnc2215_hill_response.py`
- Create: `research/data/hill-fit-summary.json`
- Modify: `research/README.md`

- [x] Read `research/data/source-data.csv` and fit the exact 15 NNC2215 measurements across five glucose concentrations and three independent experiments.
- [x] Fix each experiment's measured 0 mmol/L response as its baseline; estimate experiment-specific upper response, with EC50 and Hill slope shared across experiments.
- [x] Use deterministic 2,000-resample cluster bootstrap over the three independent experiments and record median, 2.5th/97.5th percentiles, fit residuals, bounds, software versions, input/script hashes, source, assumptions, and limitations.
- [x] Label the model descriptive; do not interpret the apparent Hill slope as site count or infer an activation threshold.

### Task 2: Generate the visual evidence

**Files:**
- Create: `research/assets/hill-fit.svg`
- Modify: `research/fit_nnc2215_hill_response.py`

- [x] Plot the three experiment-level source series, group mean and sample SD, shared-slope fit, and bootstrap prediction interval over the tested 0–20 mmol/L range.
- [x] Use Dia Health One's existing cream, purple, green, and mono-label palette; annotate that the EC50 interval extends far beyond tested concentrations.

### Task 3: Integrate the result into the dashboard

**Files:**
- Modify: `research/index.html`
- Modify: `research/research.css`
- Modify: `research/research.js`

- [x] Add a compact result panel below the source-data chart with Hill slope and EC50 uncertainty, fit links, and a clear non-identifiability statement.
- [x] Update captions so the source measurements remain primary and the fitted curve is explicitly exploratory.
- [x] Keep the long-term 6 mmol/L design target separate from the published NNC2215 response.

### Task 4: Reproduce and review the artifact

**Files:**
- Review: `research/data/hill-fit-summary.json`
- Review: `research/assets/hill-fit.svg`
- Review: hosted or local `research/index.html`

- [x] Run the analysis script twice with the pinned local chemistry environment; both runs produced the same point estimates, bootstrap ranges, and boundary-hit count.
- [x] Confirm the recorded source and script SHA-256 values match the current files.
- [x] Inspect the desktop chart and desktop/mobile panel layout; verify uncertainty and caveats remain readable.
- [x] Record the scientific interpretation and caveats in `research/README.md`, `research/COMPUTE_STATUS.md`, `research/MILESTONE_REPORT.md`, and the compute-access JSON.
