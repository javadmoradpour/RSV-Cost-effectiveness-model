# Cost-effectiveness of infant and maternal RSV immunization strategies in British Columbia, Canada

This repository contains the simulation code and supporting data for the paper **"Cost-effectiveness of infant and maternal RSV immunization strategies, in British Columbia, Canada."** The analysis evaluates recently available preventive options (nirsevimab and RSVpreF) alongside the historical palivizumab standard of care using a discrete-event simulation informed by provincial epidemiological and costing data.

## Abstract
Background: Respiratory syncytial virus (RSV) is a leading cause of lower respiratory tract infections in young children and results in significant healthcare burden and costs. To reduce the impact of RSV in this population, the monoclonal antibody palivizumab has historically been used. Recently, new preventive options have become available, including a longer-acting monoclonal antibody (nirsevimab) and a maternal vaccine (RSVpreF).

Methods: We developed a discrete-event simulation model using epidemiological and cost data from British Columbia, Canada, and published efficacy estimates. The model simulated a cohort of 100,000 newborns and followed them up to 24 months. We conducted the analysis from a healthcare system perspective, evaluating five immunization strategies: (1) the historical palivizumab standard of care for high-risk children; (2) nirsevimab for high- and moderate-risk children; (3) in-season maternal RSVpreF vaccination combined with nirsevimab for high-risk children; (4) in-season maternal RSVpreF plus nirsevimab for high- and moderate-risk children; and (5) nirsevimab for all infants. We conducted a sequential cost-effectiveness analysis, ordering strategies by cost, excluding dominated or extendedly dominated options, and evaluating the remaining strategies stepwise. To support policy interpretation, we also performed a pairwise analysis comparing each strategy directly with the historical standard of care.

Results: In the sequential analysis, strategy 2 was the most cost-effective option. Strategy 4 provided additional health gains but was not cost-effective incrementally (ICER ≈ $119,000 per QALY vs strategy 2). Strategy 5 offered the greatest overall health benefits but was the least cost-effective option. When compared directly with the historical standard of care, however, strategy 4 was cost-effective (ICER ≈ $18,000 per QALY).

Interpretation: These findings support policy recommendations to prioritize nirsevimab for high- and moderate-risk infants as the most cost-effective strategy. Maternal RSVpreF vaccination offers added health benefits and is cost-effective compared with the historical standard of care, though not when considered incrementally.

## Repository structure
- `Microsim.Rmd`: R Markdown document implementing the discrete-event simulation, including parameter definitions, model structure, and economic evaluation workflow.
- `Data/RSV_AllMetrics_withCI.csv`: Summary of key epidemiological and program parameters (means and confidence intervals) used to populate the model.
- `Data/RSV_Hosp_model_coefficients.csv`: Logistic regression coefficients for hospitalization risk as a function of month of life, with standard errors and confidence limits.

## Requirements
- R (version 4.0+ recommended)
- R packages: `pacman`, `devtools`, `dplyr`, `scales`, `ellipse`, `ggplot2`, `lazyeval`, `igraph`, `truncnorm`, `ggraph`, `reshape2`, `knitr`, `markdown`, `stringr`, `dampack`, `matrixStats`, and `darthtools` (from GitHub: `DARTH-git/darthtools`). The R Markdown document installs any missing packages via `pacman` when rendered.

## Running the analysis
1. Open the project in R or RStudio.
2. Ensure the working directory is the repository root so the `Data/` files are discoverable.
3. Knit `Microsim.Rmd` to HTML or PDF (via the **Knit** button in RStudio or by running `rmarkdown::render("Microsim.Rmd")`).
4. The document will run the discrete-event simulation, generate figures, and summarize cost-effectiveness outcomes for the five immunization strategies.

## Reproducibility notes
- The simulation seeds and number of iterations are defined in the `l_params` list within `Microsim.Rmd`. Adjust these values if you wish to explore alternative cohort sizes or probabilistic sensitivity analyses.
- Cost inputs, efficacy assumptions, and treatment pathways are parameterized near the top of `Microsim.Rmd` for easy updating as new evidence becomes available.

## Contact
For questions about the model or manuscript, please contact the corresponding author listed in the paper.
