# =============================================================================
# 11_sensitivity_analysis.R
# Deterministic sensitivity analysis (DSA):
#   - One-way sensitivity analysis (OWSA) with tornado plot
#   - Two-way sensitivity analysis (TWSA) with heatmap
#
# The sensitivity analyses re-use l_params as the base case.  Modified copies
# (never the original list) are passed to calculate_ce_out() so that all
# global state remains untouched.
#
# Note: vaccine_efficacy_list must be available as a global because the
# dampack interface (run_owsa_det / run_twsa_det) only passes l_params to
# the user-supplied function.  Since inpat/ICU efficacies are not varied in
# the SA, the global vaccine_efficacy_list is always valid.
# =============================================================================

# --- SA wrapper --------------------------------------------------------------

#' Run all strategies for one parameter configuration and return NMB.
#'
#' This function signature is required by \code{dampack::run_owsa_det()} and
#' \code{run_twsa_det()}: it must accept \code{l_params} as its first argument
#' followed by any extra keyword arguments (here \code{n_wtp}).
#'
#' @param l_params Named list of model parameters (possibly SA-modified).
#' @param n_wtp    Willingness-to-pay threshold per QALY ($/QALY).
#' @return Data frame with columns \code{Strategy}, \code{Cost},
#'   \code{Effect}, and \code{NMB}.
calculate_ce_out <- function(l_params, n_wtp = 50000) {

  # Rebuild RSVMA efficacy vectors whenever scalar efficacy changes
  l_params_mod <- add_rsvma_vectors(l_params)

  n_sim_run  <- l_params_mod$n_sim
  base_str   <- l_params_mod$base_str
  seed_run   <- l_params_mod$seed

  v_tot_cost <- v_tot_qaly <- setNames(numeric(n_str), v_names_str)

  for (j in seq_len(n_sim_run)) {
    set.seed(seed_run + j)
    df_X <- create_cohort(l_params_mod)

    for (i in seq_len(n_str)) {
      out <- MicroSim(l_params_mod, vaccine_efficacy_list, df_X,
                      v_names_str[i], seed_run + j * 10L,
                      full_output = FALSE)
      inv_n          <- 1 / n_sim_run
      v_tot_cost[i]  <- v_tot_cost[i] + sum(out$tc) * inv_n
      v_tot_qaly[i]  <- v_tot_qaly[i] + sum(out$te) * inv_n
    }
  }

  base_cost <- v_tot_cost[base_str]
  base_qaly <- v_tot_qaly[base_str]

  v_nmb <- ((v_tot_qaly - base_qaly) * n_wtp -
               (v_tot_cost - base_cost)) / 1000

  data.frame(
    Strategy = v_names_str,
    Cost     = v_tot_cost,
    Effect   = v_tot_qaly,
    NMB      = v_nmb,
    stringsAsFactors = FALSE
  )
}


# --- SA parameter ranges -----------------------------------------------------

#' Build lower- and upper-bound parameter lists for sensitivity analysis.
#'
#' @param l_params Base-case parameter list.
#' @return A list with \code{prm_l} (lower bounds) and \code{prm_u} (upper
#'   bounds), each a named list suitable for passing to
#'   \code{dampack::run_owsa_det()}.
build_sa_ranges <- function(l_params) {
  rng25   <- function(x) list(l = x * 0.75, u = x * 1.25)
  get_ci  <- function(var, col) {
    l_params$all_metrics[l_params$all_metrics$Variable == var, col]
  }

  prm_l <- list(
    PD_ICU  = rng25(l_params$PD_ICU)$l,
    PD_PW   = rng25(l_params$PD_PW)$l,
    c_ED    = rng25(l_params$c_ED)$l,
    c_DV    = rng25(l_params$c_DV)$l,

    c_RSVpreF           = 100,
    c_nirsevimab        = 400,
    c_Admin_vac         = 0,
    c_palivizumab       = rng25(l_params$c_palivizumab)$l,
    c_Admin_palivizumab = 0,

    Du_ICU_Infants = 0.0145,
    Du_PW_Infants  = 0.01,
    Du_ED_Infants  = 0.008,
    Du_DV_Infants  = 0.005,

    Prob_ICU_LowRisk      = get_ci("Prob_ICU_LowRisk",      "LCL"),
    Prob_ICU_ModerateRisk = get_ci("Prob_ICU_ModerateRisk", "LCL"),
    Prob_ICU_HighRisk     = get_ci("Prob_ICU_HighRisk",     "LCL"),

    LOS_hosp_LowRisk      = get_ci("LOS_hosp_LowRisk",      "LCL"),
    LOS_hosp_ModerateRisk = get_ci("LOS_hosp_ModerateRisk", "LCL"),
    LOS_hosp_HighRisk     = get_ci("LOS_hosp_HighRisk",     "LCL"),
    LOS_icu_LowRisk       = get_ci("LOS_icu_LowRisk",       "LCL"),
    LOS_icu_ModerateRisk  = get_ci("LOS_icu_ModerateRisk",  "LCL"),
    LOS_icu_HighRisk      = get_ci("LOS_icu_HighRisk",      "LCL"),

    coverage_RSVpreF                 = 0.5,
    coverage_nirsevimab_ModerateRisk = 0.6,

    eff_RSVMA_palivizumab     = 0.19,
    efficacy_RSVMA_RSVpreF    = 0.287,
    efficacy_RSVMA_nirsevimab = 0.659,

    p_mort_Inp = get_ci("monthly_mortality", "LCL")
  )

  prm_u <- list(
    PD_ICU  = rng25(l_params$PD_ICU)$u,
    PD_PW   = rng25(l_params$PD_PW)$u,
    c_ED    = rng25(l_params$c_ED)$u,
    c_DV    = rng25(l_params$c_DV)$u,

    c_RSVpreF           = 250,
    c_nirsevimab        = 950,
    c_Admin_vac         = 30,
    c_palivizumab       = rng25(l_params$c_palivizumab)$u,
    c_Admin_palivizumab = rng25(l_params$c_Admin_palivizumab)$u,

    Du_ICU_Infants = 0.1053,
    Du_PW_Infants  = 0.0726,
    Du_ED_Infants  = 0.0454,
    Du_DV_Infants  = 0.0454,

    Prob_ICU_LowRisk      = get_ci("Prob_ICU_LowRisk",      "UCL"),
    Prob_ICU_ModerateRisk = get_ci("Prob_ICU_ModerateRisk", "UCL"),
    Prob_ICU_HighRisk     = get_ci("Prob_ICU_HighRisk",     "UCL"),

    LOS_hosp_LowRisk      = get_ci("LOS_hosp_LowRisk",      "UCL"),
    LOS_hosp_ModerateRisk = get_ci("LOS_hosp_ModerateRisk", "UCL"),
    LOS_hosp_HighRisk     = get_ci("LOS_hosp_HighRisk",     "UCL"),
    LOS_icu_LowRisk       = get_ci("LOS_icu_LowRisk",       "UCL"),
    LOS_icu_ModerateRisk  = get_ci("LOS_icu_ModerateRisk",  "UCL"),
    LOS_icu_HighRisk      = get_ci("LOS_icu_HighRisk",      "UCL"),

    coverage_RSVpreF                 = 0.8,
    coverage_nirsevimab_ModerateRisk = 0.9,

    eff_RSVMA_palivizumab     = 0.90,
    efficacy_RSVMA_RSVpreF    = 0.689,
    efficacy_RSVMA_nirsevimab = 0.877,

    p_mort_Inp = 0.00192  # upper assumption (NACI)
  )

  list(prm_l = prm_l, prm_u = prm_u)
}


# --- OWSA --------------------------------------------------------------------

#' Run one-way sensitivity analysis and produce tornado and curve plots.
#'
#' @param l_params    Base-case parameter list.
#' @param output_dir  Directory for PNG output files.
#' @param nsamp       Number of parameter values sampled per parameter.
#' @param n_wtp       Willingness-to-pay threshold ($/QALY).
#' @param n_sim_sa    Replications used for each OWSA model evaluation.
#'   Defaults to \code{l_params$n_sim_sa} (separate from the base-case
#'   \code{n_sim}).  Increase if SA curves still appear noisy.
#' @param n_i_sa      Cohort size used for each OWSA evaluation.  Defaults to
#'   \code{l_params$n_i}.
#' @param analytic_mort Logical (default \code{TRUE}).  When \code{TRUE},
#'   replaces the stochastic end-of-simulation death QALY penalty with an
#'   expected-value calculation at every cycle.  This eliminates Monte-Carlo
#'   noise from the rare inpatient-mortality event and produces smooth OWSA
#'   curves for the inpatient-mortality parameter without any change to the
#'   state-transition logic.  See \code{Effs()} for details.
#' @return The \code{owsa_nmb} object returned by \code{run_owsa_det()}.
run_owsa <- function(l_params, output_dir = "./Output", nsamp = 20,
                     n_wtp = 50000,
                     n_sim_sa      = l_params$n_sim_sa,
                     n_i_sa        = l_params$n_i_sa,
                     analytic_mort = TRUE) {

  # Build a modified parameter list for SA runs
  l_params_sa               <- l_params
  l_params_sa$n_sim         <- n_sim_sa
  l_params_sa$n_i           <- n_i_sa
  l_params_sa$analytic_mort <- analytic_mort

  sa_ranges <- build_sa_ranges(l_params_sa)
  prm_l <- sa_ranges$prm_l
  prm_u <- sa_ranges$prm_u

  owsa_pars <- c(
    "c_RSVpreF", "c_nirsevimab",
    "PD_ICU", "PD_PW", "c_ED", "c_DV",
    "coverage_RSVpreF", "coverage_nirsevimab_ModerateRisk",
    "eff_RSVMA_palivizumab", "efficacy_RSVMA_RSVpreF",
    "efficacy_RSVMA_nirsevimab", "p_mort_Inp"
  )

  df_params_owsa <- data.frame(
    pars = owsa_pars,
    min  = unlist(prm_l[owsa_pars]),
    max  = unlist(prm_u[owsa_pars]),
    stringsAsFactors = FALSE
  )

  owsa_nmb <- run_owsa_det(
    params_range    = df_params_owsa,
    params_basecase = l_params_sa,
    nsamp           = nsamp,
    FUN             = calculate_ce_out,
    outcomes        = "NMB",
    strategies      = v_names_str,
    n_wtp           = n_wtp
  )

  nice_params <- c(
    c_ED                              = "ED cost",
    c_DV                              = "DV cost",
    c_nirsevimab                      = "Nirsevimab price",
    c_RSVpreF                         = "RSVpreF price",
    coverage_RSVpreF                  = "RSVpreF coverage",
    coverage_nirsevimab_ModerateRisk  = "Nirsevimab coverage (mod)",
    eff_RSVMA_palivizumab             = "Palivizumab efficacy",
    efficacy_RSVMA_nirsevimab         = "Nirsevimab efficacy",
    efficacy_RSVMA_RSVpreF            = "RSVpreF efficacy",
    p_mort_Inp                        = "Inpatient mortality",
    PD_ICU                            = "ICU cost per day",
    PD_PW                             = "Ward cost per day"
  )

  nice_strats <- setNames(v_names_str, make.names(v_names_str))

  owsa_nice <- owsa_nmb %>%
    mutate(
      parameter = factor(parameter,
                         levels = names(nice_params), labels = nice_params),
      strategy  = factor(strategy,
                         levels = names(nice_strats),  labels = nice_strats)
    )

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  p_tornado <- owsa_tornado(owsa = owsa_nice)
  ggsave(file.path(output_dir, "owsa_tornado.png"),
         p_tornado, width = 8, height = 6, dpi = 300)

  p_owsa_all <- plot(owsa_nice, txtsize = 11, n_x_ticks = 5,
                     facet_scales = "free") +
    labs(y = "E[NMB]", x = "Parameter value") +
    theme(legend.text = element_text(size = 11),
          axis.text.x = element_text(size = 7))
  ggsave(file.path(output_dir, "owsa_all_strategies.png"),
         p_owsa_all, width = 12, height = 8, dpi = 500)

  p_owsa_excl <- owsa_nice %>%
    filter(strategy != "Nirsevimab (Universal)") %>%
    droplevels() %>%
    plot(txtsize = 11, n_x_ticks = 5, facet_scales = "free") +
    labs(y = "E[NMB]", x = "Parameter value") +
    theme(legend.text = element_text(size = 11),
          axis.text.x = element_text(size = 7))
  ggsave(file.path(output_dir, "owsa_excl_universal.png"),
         p_owsa_excl, width = 12, height = 8, dpi = 500)

  cat("OWSA plots saved to:", output_dir, "\n")
  owsa_nmb
}


# --- TWSA --------------------------------------------------------------------

#' Run two-way sensitivity analysis (RSVpreF price × nirsevimab price).
#'
#' @param l_params      Base-case parameter list.
#' @param n_i_twsa      Cohort size for TWSA runs. Defaults to
#'   \code{l_params$n_i_sa}; override to a larger value for extra stability.
#' @param output_dir    Directory for output files.
#' @param nsamp         Grid points per price axis.
#' @param n_wtp         Willingness-to-pay threshold ($/QALY).
#' @param analytic_mort Logical (default \code{TRUE}).  See \code{run_owsa()}
#'   for details.
#' @return The \code{twsa_nmb} object returned by \code{run_twsa_det()}.
run_twsa <- function(l_params, n_i_twsa = l_params$n_i_sa, output_dir = "./Output",
                     nsamp = 20, n_wtp = 50000, analytic_mort = TRUE) {

  # Override cohort size and mort mode without mutating any global variable.
  l_params_twsa               <- l_params
  l_params_twsa$n_i           <- n_i_twsa
  l_params_twsa$analytic_mort <- analytic_mort

  df_params_twsa <- data.frame(
    pars = c("c_RSVpreF", "c_nirsevimab"),
    min  = c(0,   0),
    max  = c(300, 1000),
    stringsAsFactors = FALSE
  )

  twsa_nmb <- run_twsa_det(
    params_range    = df_params_twsa,
    params_basecase = l_params_twsa,
    nsamp           = nsamp,
    FUN             = calculate_ce_out,
    outcomes        = "NMB",
    strategies      = v_names_str,
    n_wtp           = n_wtp
  )

  nice_strats_twsa <- setNames(v_names_str, make.names(v_names_str))

  twsa_nice <- twsa_nmb %>%
    mutate(strategy = factor(strategy,
                             levels = names(nice_strats_twsa),
                             labels = nice_strats_twsa)) %>%
    rename(`RSVpreF price` = c_RSVpreF, `Nirsevimab price` = c_nirsevimab)

  p_twsa <- plot(twsa_nice, txtsize = 12, facet_scales = "free") +
    scale_fill_manual(values = c(
      "No Intervention"                   = "#999999",
      "Palivizumab"                       = "#000000",
      "RSVpreF"                           = "#0072B2",
      "Nirsevimab (High & Mod)"           = "#E69F00",
      "RSVpreF + Nirsevimab (High)"       = "#56B4E9",
      "RSVpreF + Nirsevimab (High & Mod)" = "#009E73",
      "Nirsevimab (Universal)"            = "#F0E442"
    )) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
    theme(legend.text = element_text(size = 10))

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(output_dir, "twsa_nmb.png"),
         p_twsa, width = 8, height = 6, dpi = 500)
  cat("TWSA plot saved to:", output_dir, "\n")

  twsa_nmb
}
