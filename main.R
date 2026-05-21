# =============================================================================
# main.R
# Entry point for the RSV cost-effectiveness microsimulation.
#
# Run this script from the project root to reproduce the full analysis:
#   source("main.R")          # interactive session
#   Rscript main.R            # command line
#
# What it does (in order):
#   1. Loads packages
#   2. Defines fixed constants (time horizon, strategy names, …)
#   3. Reads data and builds the canonical parameter list (l_params)
#   4. Computes vaccine-efficacy profiles
#   5. Runs the base-case microsimulation
#   6. Performs cost-effectiveness analysis (ICER table + CE-plane plot)
#   7. Runs one-way sensitivity analysis (tornado and curve plots)
#   8. Runs two-way sensitivity analysis (price heatmap)
#
# All outputs (CSV, PNG) are written to ./Output/.
#
# To run a partial analysis, source main.R through step 4 and then call
# individual functions (run_owsa, run_twsa, etc.) as needed.
# =============================================================================

rm(list = ls())
gc()

# --- 1. Source all module scripts (order matters) ----------------------------
source("R/01_packages.R")
source("R/02_parameters.R")
source("R/03_data_loading.R")
source("R/04_vaccine_efficacy.R")
source("R/05_cohort.R")
source("R/06_strategies.R")
source("R/07_model_functions.R")
source("R/08_microsimulation.R")
source("R/09_run_simulation.R")
source("R/10_cea.R")
source("R/11_sensitivity_analysis.R")

# --- 2. Build parameters and efficacy objects --------------------------------
l_params              <- build_params(data_path = "./Data")
l_params              <- add_rsvma_vectors(l_params)
vaccine_efficacy_list <- build_vaccine_efficacy_list(l_params)

# --- 3. Base-case simulation -------------------------------------------------
results_table <- run_base_case(
  l_params              = l_params,
  vaccine_efficacy_list = vaccine_efficacy_list,
  output_dir            = "./Output"
)
print(results_table)

# --- 4. Cost-effectiveness analysis -----------------------------------------
cea_out <- run_cea(results_table, output_dir = "./Output")
print(cea_out$table_cea)

# --- 5. One-way sensitivity analysis ----------------------------------------
owsa_nmb <- run_owsa(
  l_params   = l_params,
  output_dir = "./Output",
  nsamp      = 20,
  n_wtp      = 50000
)

# --- 6. Two-way sensitivity analysis ----------------------------------------
twsa_nmb <- run_twsa(
  l_params   = l_params,
  output_dir = "./Output",
  nsamp      = 20,
  n_wtp      = 50000
)
