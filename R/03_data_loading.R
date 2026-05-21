# =============================================================================
# 03_data_loading.R
# Builds the canonical model parameter list (l_params) by combining
# hard-coded inputs with values read from the two data CSV files.
#
# l_params is the single source of truth for all model parameters and is
# passed to every model function.  Sensitivity analyses create modified
# copies of l_params without touching any global state.
# =============================================================================

#' Build the canonical model parameter list.
#'
#' @param data_path Directory containing the two CSV data files.
#' @return A named list of all model parameters.
#' @details
#'   Two external files are required:
#'   \itemize{
#'     \item \code{RSV_AllMetrics_withCI.csv} — epidemiological estimates with
#'       confidence intervals.
#'     \item \code{RSV_Hosp_model_coefficients.csv} — logistic-regression
#'       coefficients for monthly hospitalisation probability.
#'   }
build_params <- function(data_path = "./Data") {

  # --- Helper: pull one value from the metrics table -----------------------
  all_metrics <- read.csv(file.path(data_path, "RSV_AllMetrics_withCI.csv"))
  get_metric  <- function(var, col = "Mean") {
    all_metrics[all_metrics$Variable == var, col]
  }

  # --- Hospitalisation model coefficients ----------------------------------
  hosp_coef_raw   <- read.csv(
    file.path(data_path, "RSV_Hosp_model_coefficients.csv"),
    stringsAsFactors = FALSE
  )
  hosp_model_beta <- setNames(hosp_coef_raw$estimate, hosp_coef_raw$term)

  # --- Population risk distribution ----------------------------------------
  prob_high <- get_metric("prob_HighRisk")
  prob_low  <- get_metric("prob_LowRisk")
  v_risk_prob <- c(
    HighRisk     = prob_high,
    LowRisk      = prob_low,
    ModerateRisk = 1 - prob_high - prob_low
  )

  # -------------------------------------------------------------------------
  list(
    # Simulation settings
    base_str = "No Intervention",
    n_sim    = n_sim,

    # Population
    v_risk_prob             = v_risk_prob,
    coverage_palivizumab_y2 = get_metric("pvz_coverage_second_season"),

    # Mortality
    p_mort     = get_metric("monthly_mortality"),
    p_mort_Inp = get_metric("monthly_mortality"),

    # ICU probability conditional on hospitalisation
    Prob_ICU_LowRisk          = get_metric("Prob_ICU_LowRisk"),
    Prob_ICU_ModerateRisk     = get_metric("Prob_ICU_ModerateRisk"),
    Prob_ICU_HighRisk         = get_metric("Prob_ICU_HighRisk"),
    RR_Prob_ICU_Second_season = get_metric("RR_Prob_ICU_Second"),

    # Length of stay (days) — ward and ICU by risk level
    LOS_hosp_LowRisk      = get_metric("LOS_hosp_LowRisk"),
    LOS_hosp_ModerateRisk = get_metric("LOS_hosp_ModerateRisk"),
    LOS_hosp_HighRisk     = get_metric("LOS_hosp_HighRisk"),
    LOS_icu_LowRisk       = get_metric("LOS_icu_LowRisk"),
    LOS_icu_ModerateRisk  = get_metric("LOS_icu_ModerateRisk"),
    LOS_icu_HighRisk      = get_metric("LOS_icu_HighRisk"),

    # Unit costs
    PD_ICU              = 4615,  # per day, ICU (Vadlamudi 2025)
    PD_PW               = 2283,  # per day, paediatric ward (Vadlamudi 2025)
    c_ED                = 342,   # ED visit (Shoukat 2023)
    c_DV                = 229,   # doctor visit (Shoukat 2023)
    c_RSVpreF           = 125,   # RSVpreF per dose
    c_nirsevimab        = 450,   # nirsevimab per dose
    c_Admin_vac         = 15,    # administration per vaccine dose
    c_palivizumab       = 3600,  # palivizumab per RSV season
    c_Admin_palivizumab = 150,   # administration per palivizumab season

    # Disutilities (infant and caregiver components retained separately
    # so analysts can toggle caregiver disutility easily)
    Du_ICU_Infants    = 0.0245,
    Du_ICU_Caregivers = 0.0097,
    Du_PW_Infants     = 0.0169,
    Du_PW_Caregivers  = 0.0067,
    Du_ED_Infants     = 0.0135,
    Du_ED_Caregivers  = 0.00675,
    Du_DV_Infants     = 0.00845,
    Du_DV_Caregivers  = 0.00423,

    # Vaccine coverage
    coverage_RSVpreF                 = 0.65,
    coverage_nirsevimab_ModerateRisk = 0.70,

    # Vaccine efficacy scalars — varied in sensitivity analyses
    eff_RSVMA_palivizumab     = eff_RSVMA_palivizumab,
    eff_inpat_palivizumab     = eff_inpat_palivizumab,
    eff_ICU_palivizumab       = eff_ICU_palivizumab,
    efficacy_RSVMA_RSVpreF    = efficacy_RSVMA_RSVpreF,
    efficacy_inpat_RSVpreF    = efficacy_inpat_RSVpreF,
    efficacy_ICU_RSVpreF      = efficacy_ICU_RSVpreF,
    efficacy_RSVMA_nirsevimab = efficacy_RSVMA_nirsevimab,
    efficacy_inpat_nirsevimab = efficacy_inpat_nirsevimab,
    efficacy_ICU_nirsevimab   = efficacy_ICU_nirsevimab,

    # Hospitalisation model (used in compute_hosp_prob())
    hosp_model_beta = hosp_model_beta,

    # Simulation settings carried inside l_params so SA wrappers can
    # override them (e.g. TWSA uses a larger cohort) without touching globals.
    n_i  = n_i,
    seed = seed,

    # Raw metrics table — retained so sensitivity-analysis bounds can be
    # looked up from confidence intervals without re-reading the file.
    all_metrics = all_metrics
  )
}
