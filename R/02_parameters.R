# =============================================================================
# 02_parameters.R
# Fixed model constants — values that never change across scenarios or
# sensitivity analyses.  All time-varying or data-loaded parameters live
# in l_params (built in 03_data_loading.R).
# =============================================================================

options(scipen = 999)

# --- Time horizon and cohort size ---
cycle_length <- 1 / 12  # one month, expressed in years
n_cycles     <- 24      # follow-up length (monthly cycles)
n_i          <- 10000   # cohort size (births per year)

v_names_cycles <- paste("cycle", 0:n_cycles)
v_names_states <- c("NB", "H", "RSV", "D")
n_states       <- length(v_names_states)

v_names_treatments <- c("NT", "DV", "ED", "ICU", "PW")
n_treatments       <- length(v_names_treatments)

# --- Discount rates ---
d_c <- 0.00
d_e <- 0.00
v_dwc <- 1 / ((1 + d_e * cycle_length) ^ (0:n_cycles))
v_dwe <- 1 / ((1 + d_c * cycle_length) ^ (0:n_cycles))

# --- Immunisation strategies evaluated ---
v_names_str <- c(
  "No Intervention",
  "Palivizumab",
  "Nirsevimab (High & Mod)",
  "RSVpreF",
  "RSVpreF + Nirsevimab (High)",
  "RSVpreF + Nirsevimab (High & Mod)",
  "Nirsevimab (Universal)"
)
n_str <- length(v_names_str)

# --- Simulation replications ---
n_sim <- 10   # Monte Carlo iterations for base case
seed  <- 100

# --- RSV season timing (calendar months) ---
sim_start_month        <- 5   # simulation starts in May
rsv_season_start_month <- 11  # RSV season begins in November
sim_to_rsv_month_gap   <- rsv_season_start_month - sim_start_month  # = 6

# --- Outpatient visit probabilities by age group (proportion of RSV cases) ---
# First 6 months, months 7-12, months 13-24
prob_RSVMA_Inp <- c(rep(0.083, 6), rep(0.021, 6), rep(0.027, n_cycles - 12))
prob_Outp_ED   <- c(rep(0.29,  6), rep(0.24,  6), rep(0.33,  n_cycles - 12))

# --- Utility parameters ---
u_H     <- 1      # healthy-state utility weight
u_Total <- 45.3   # discounted life expectancy at birth (years)
u_rs    <- u_Total - (n_cycles - 6) * cycle_length  # residual life lost at death

# --- Vaccine duration parameters ---
valid_months_RSVpreF    <- 5   # months of constant peak efficacy
valid_months_nirsevimab <- 5
zero_month              <- 10  # month at which efficacy declines to zero

# --- Vaccine efficacy constants (scalars used throughout) ---
eff_RSVMA_palivizumab <- 0.70   # NACI estimate
eff_inpat_palivizumab <- 0.82
eff_ICU_palivizumab   <- 0.82

efficacy_RSVMA_RSVpreF <- 0.525  # Kampmann 2023
efficacy_inpat_RSVpreF <- 0.564
efficacy_ICU_RSVpreF   <- 0.709

efficacy_RSVMA_nirsevimab <- 0.795  # Simões 2023
efficacy_inpat_nirsevimab <- 0.773
efficacy_ICU_nirsevimab   <- 0.860

# Nirsevimab coverage for high-risk infants is assumed to be 100 %
coverage_nirsevimab_HighRisk <- 1.0
