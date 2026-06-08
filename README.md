# RSV Cost-Effectiveness Model — British Columbia

An individual-level microsimulation model evaluating the cost-effectiveness of
seven RSV immunisation strategies for infants in British Columbia, Canada.

---

## Overview

The model simulates monthly health-state transitions for a cohort of 10,000
infants over a 24-month follow-up period and computes costs (vaccine,
inpatient, outpatient) and quality-adjusted life years (QALYs) for each
immunisation strategy. Results feed into a cost-effectiveness analysis (CEA)
with one-way and two-way deterministic sensitivity analyses.

This repository contains the simulation code and supporting data for the paper:
**"Cost-effectiveness of infant and maternal RSV immunization strategies in
British Columbia, Canada."**

**Strategies evaluated**

| # | Strategy |
|---|----------|
| 1 | No Intervention |
| 2 | Palivizumab (historical standard of care) |
| 3 | Nirsevimab — high- and moderate-risk infants (with catch-up) |
| 4 | RSVpreF maternal vaccine only |
| 5 | RSVpreF + nirsevimab (high-risk infants) |
| 6 | RSVpreF + nirsevimab (high- and moderate-risk infants) |
| 7 | Nirsevimab (All, Low In-Season) |
| 8 | RSVpreF + Nirsevimab (All, Low In-Season) |

**Key findings:** Nirsevimab for high- and moderate-risk infants (Strategy 3)
is the most cost-effective option at a $50,000/QALY threshold. The combined
RSVpreF + nirsevimab strategy (Strategy 6) offers greater health benefits but
at an incremental cost of ~$119K/QALY vs. Strategy 3; when compared directly
to palivizumab, Strategy 6 is cost-effective at ~$18K/QALY.

---

## Model structure

**Health states:** Newborn (NB) → Healthy (H) ↔ RSV infection (RSV) → Death (D)

**Treatment types:** No treatment (NT), Doctor visit (DV), Emergency department
(ED), Paediatric ward (PW), ICU

**Hospitalisation probability** is estimated from a logistic regression model
fitted to BC surveillance data, with predictors for age (polynomial in
month-of-life), RSV seasonality (calendar month), risk level, and palivizumab
use.

**Vaccine efficacy** is modelled as constant for the first 5 months
post-administration, then declining linearly to zero by month 10.

---

## Project structure

```
RSV-Cost-effectiveness-model/
├── main.R                          # Entry point — run this to reproduce all results
├── R/
│   ├── 01_packages.R               # Package loading (CRAN + GitHub)
│   ├── 02_parameters.R             # Fixed constants: time horizon, state names, efficacy scalars
│   ├── 03_data_loading.R           # build_params(): reads CSVs, returns canonical l_params
│   ├── 04_vaccine_efficacy.R       # vaccine_efficacy(), build_vaccine_efficacy_list()
│   ├── 05_cohort.R                 # create_cohort(): individual-level characteristics
│   ├── 06_strategies.R             # Vaccination strategy functions + apply_strategy() dispatch
│   ├── 07_model_functions.R        # Hospitalisation probability, state transitions,
│   │                               #   vaccine modifiers, treatment assignment, costs, QALYs
│   ├── 08_microsimulation.R        # MicroSim(): the core simulation loop
│   ├── 09_run_simulation.R         # run_base_case(): runs all strategies × replications
│   ├── 10_cea.R                    # run_cea(): ICER table and CE-plane plot
│   └── 11_sensitivity_analysis.R   # calculate_ce_out(), run_owsa(), run_twsa()
├── Data/
│   ├── RSV_AllMetrics_withCI.csv               # Epidemiological parameters with 95% CIs
│   └── RSV_Hosp_model_coefficients.csv         # Logistic regression coefficients
├── Output/                         # Generated CSV results and PNG plots
├── report/
│   └── analysis_report.Rmd         # Self-contained R Markdown report
└── README.md
```

---

## Requirements

- **R** ≥ 4.0
- Packages are installed automatically on first run via `pacman`:
  `dplyr`, `ggplot2`, `scales`, `dampack`, `matrixStats`, `truncnorm`,
  `reshape2`, `stringr`, `knitr`, `markdown`, `devtools`, and
  `DARTH-git/darthtools` (from GitHub)

---

## How to run

### Full analysis

Open R (or RStudio) in the project root and run:

```r
source("main.R")
```

This will, in order:

1. Load all packages
2. Build the parameter list from the data files
3. Run the base-case simulation (10 iterations × 7 strategies)
4. Produce the ICER table and CE-plane plot
5. Run the one-way sensitivity analysis (tornado and curve plots)
6. Run the two-way sensitivity analysis (RSVpreF × nirsevimab price heatmap)

All outputs (CSV and PNG files) are written to `./Output/`.

### Rendered report

Knit `report/analysis_report.Rmd` from RStudio, or run:

```r
rmarkdown::render("report/analysis_report.Rmd", output_dir = "Output")
```

### Partial runs

Source `main.R` up through the parameter-building step, then call individual
functions as needed:

```r
source("R/01_packages.R")
source("R/02_parameters.R")
source("R/03_data_loading.R")
source("R/04_vaccine_efficacy.R")
# ... source remaining R/ scripts ...

l_params              <- build_params()
l_params              <- add_rsvma_vectors(l_params)
vaccine_efficacy_list <- build_vaccine_efficacy_list(l_params)

# Run only the sensitivity analysis at a different WTP threshold
run_owsa(l_params, n_wtp = 100000)
```

---

## Extending the model

### Adding a new immunisation strategy

1. Write a strategy function in `R/06_strategies.R` following the same pattern
   as `Palivizumab()` or `Nirsevimab_HighRisk_WithCatchup()`.
2. Register it in the `strategy_fns` list inside `apply_strategy()`.
3. Add the strategy name to `v_names_str` in `R/02_parameters.R`.

### Changing model parameters

All base-case parameters are defined in `build_params()` (`R/03_data_loading.R`)
and the fixed constants in `R/02_parameters.R`. Sensitivity-analysis ranges are
set in `build_sa_ranges()` (`R/11_sensitivity_analysis.R`).

To run a quick scenario without modifying any files:

```r
l_alt <- l_params
l_alt$c_nirsevimab <- 300          # lower nirsevimab price
l_alt <- add_rsvma_vectors(l_alt)  # rebuild RSVMA efficacy vectors
run_base_case(l_alt, vaccine_efficacy_list)
```

---

## Data sources

| Parameter | Source |
|-----------|--------|
| ICU and ward unit costs | Vadlamudi et al. 2025 |
| ED and doctor-visit costs | Shoukat et al. 2023 |
| Palivizumab efficacy | NACI guidance |
| RSVpreF (Abrysvo) efficacy | Kampmann et al. 2023 |
| Nirsevimab (Beyfortus) efficacy | Simões et al. 2023 |
| Hospitalisation rates, ICU probabilities, length of stay | BC provincial surveillance data |

---

## Citation

If you use this model in your work, please cite:

> Taleshi, Javad, et al. "Cost-effectiveness of infant and maternal RSV immunization strategies, in British Columbia, Canada." Vaccine 68 (2025): 127936.

## Contact

For questions about the model or manuscript, please contact the corresponding
author listed in the paper.
