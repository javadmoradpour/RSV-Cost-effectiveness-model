# =============================================================================
# 07_model_functions.R
# Core model functions: hospitalisation probability, state transitions,
# vaccine effect modifiers, treatment assignment, costs, and utilities.
#
# Performance notes
# -----------------
# * precompute_lp_static() and precompute_lp_calendar() extract the parts
#   of the hospitalisation log-odds that are constant across cycles (or
#   across individuals).  MicroSim() calls these once before the cycle loop
#   and passes the pre-computed arrays to state_transition().
# * Vaccine_effect_modifier() is called exactly once per cycle in
#   treatment_type() by consolidating the RSVMA and inpat/ICU computations
#   that were previously two separate calls.
# =============================================================================

# --- Hospitalisation probability helpers -------------------------------------

#' Pre-compute the individual-level static component of the hospitalisation LP.
#'
#' @param df_X          Cohort data frame.
#' @param hosp_model_beta Named numeric vector of logistic-regression
#'   coefficients.
#' @return Numeric vector of length \code{nrow(df_X)} containing the intercept
#'   plus risk-level and pvz contributions for each infant.  This vector does
#'   not change during the simulation.
precompute_lp_static <- function(df_X, hosp_model_beta) {
  lp_risk <- numeric(nrow(df_X))
  lp_risk[df_X$Risk == "ModerateRisk"] <- hosp_model_beta["riskModerateRisk"]
  lp_risk[df_X$Risk == "HighRisk"]     <- hosp_model_beta["riskHighRisk"]

  lp_pvz <- hosp_model_beta["pvz"] * df_X$pvz

  as.numeric(hosp_model_beta["(Intercept)"] + lp_pvz + lp_risk)
}


#' Pre-compute the calendar-month contribution to the hospitalisation LP for
#' every cycle.
#'
#' @param hosp_model_beta Named numeric vector of logistic-regression
#'   coefficients.
#' @return Numeric vector of length \code{n_cycles}.  August (the reference
#'   month in the logistic model) contributes 0.
precompute_lp_calendar <- function(hosp_model_beta) {
  sapply(seq_len(n_cycles), function(t) {
    cal_month <- ((t - 1L) + (sim_start_month - 1L)) %% 12L + 1L
    if (cal_month == 8L) {
      0
    } else {
      key <- paste0("calendar_month", cal_month)
      hosp_model_beta[key]
    }
  })
}


#' Compute per-infant monthly hospitalisation probability for cycle t.
#'
#' @param t             Integer cycle index (1-based).
#' @param df_X          Cohort data frame.
#' @param lp_static     Pre-computed static LP vector (from
#'   \code{precompute_lp_static()}).
#' @param lp_calendar   Pre-computed calendar LP vector (from
#'   \code{precompute_lp_calendar()}).
#' @param l_params      Model parameter list (needs \code{hosp_model_beta}
#'   and \code{eff_inpat_palivizumab}).
#' @param da            Small offset added to month-of-life to avoid division
#'   by zero at age 0 (default 0.25, matching the original model fit).
#' @return Numeric vector of length \code{nrow(df_X)} with hospitalisation
#'   probabilities.  Infants not yet born receive probability 0.
#' @details
#'   For palivizumab recipients, the model coefficient already reflects their
#'   reduced risk.  The function reverses that effect (multiplies by
#'   \code{1/(1-eff_inpat_palivizumab)}) to recover the counterfactual
#'   unvaccinated probability; vaccine effects are then applied uniformly in
#'   \code{Vaccine_effect_modifier()}.
compute_hosp_prob <- function(t, df_X, lp_static, lp_calendar, l_params,
                               da = 0.25) {
  born <- t >= df_X$MonthBorn
  prob <- numeric(nrow(df_X))

  if (any(born)) {
    beta <- l_params$hosp_model_beta
    mol  <- t - df_X$MonthBorn[born] + da  # month-of-life (offset)

    lp_mol <- beta["I(1/month_of_life)"] / mol +
              beta["month_of_life"]       * mol +
              beta["I(month_of_life^2)"]  * mol^2 +
              beta["I(month_of_life^3)"]  * mol^3

    lp <- lp_static[born] + lp_calendar[t] + as.numeric(lp_mol)
    prob[born] <- 1 / (1 + exp(-lp))

    # Reverse the pvz coefficient so subsequent vaccine modifiers are applied
    # on the unprotected baseline
    pvz_born <- born & df_X$pvz == 1L
    if (any(pvz_born)) {
      adj <- 1 / (1 - l_params$eff_inpat_palivizumab)
      prob[pvz_born] <- prob[pvz_born] * adj
    }
  }

  prob
}


# --- State transition ---------------------------------------------------------

#' Simulate monthly health-state transitions for all infants.
#'
#' @param l_params    Model parameter list.
#' @param M_t         Character vector of current states (\code{"NB"},
#'   \code{"H"}, \code{"RSV"}, \code{"D"}).
#' @param Tr_t        Character vector of current treatments (used to set
#'   inpatient mortality for RSV state).
#' @param t           Integer cycle index.
#' @param df_X        Cohort data frame.
#' @param lp_static   Pre-computed static LP vector.
#' @param lp_calendar Pre-computed calendar LP vector.
#' @return Character vector of next-cycle states.
state_transition <- function(l_params, M_t, Tr_t, t, df_X,
                              lp_static, lp_calendar) {
  p_NB <- as.numeric(df_X$MonthBorn == t)

  p_H_D   <- l_params$p_mort
  p_RSV_D <- ifelse(Tr_t %in% c("ICU", "PW"), l_params$p_mort_Inp,
                    l_params$p_mort)

  p_H_Inp <- compute_hosp_prob(t, df_X, lp_static, lp_calendar, l_params)

  # Convert inpatient probability to overall RSV-episode probability
  month_of_life <- ifelse(t >= df_X$MonthBorn, t - df_X$MonthBorn, NA_real_)
  denom         <- ifelse(!is.na(month_of_life),
                          prob_RSVMA_Inp[month_of_life + 1L], 1)
  p_H_RSV <- pmin(ifelse(!is.na(month_of_life), p_H_Inp / denom, 0), 1)

  m_p_t <- matrix(0, nrow = n_states, ncol = length(M_t),
                  dimnames = list(v_names_states, NULL))

  idx <- M_t == "NB"
  m_p_t["NB",  idx] <- 1 - p_NB[idx]
  m_p_t["H",   idx] <- p_NB[idx] * (1 - p_H_RSV[idx])
  m_p_t["RSV", idx] <- p_NB[idx] *      p_H_RSV[idx]

  idx <- M_t == "H"
  m_p_t["H",   idx] <- (1 - p_H_D) * (1 - p_H_RSV[idx])
  m_p_t["RSV", idx] <- (1 - p_H_D) *      p_H_RSV[idx]
  m_p_t["D",   idx] <- p_H_D

  idx <- M_t == "RSV"
  m_p_t["H", idx] <- 1 - p_RSV_D[idx]
  m_p_t["D", idx] <-     p_RSV_D[idx]

  idx <- M_t == "D"
  m_p_t["D", idx] <- 1

  samplev(t(m_p_t), 1)
}


# --- Vaccine effect modifier -------------------------------------------------

#' Compute per-infant probability modifiers for RSV risk, hospitalisation, and
#' ICU admission based on current vaccination status.
#'
#' @param vaccine_efficacy_list Named list of efficacy vectors/scalars (from
#'   \code{build_vaccine_efficacy_list()}).  Provides inpatient and ICU
#'   components.
#' @param l_params              Model parameter list.  Provides RSVMA
#'   components (which can vary in sensitivity analyses).
#' @param t                     Integer cycle index.
#' @param df_X                  Cohort data frame.
#' @return A list with three numeric vectors of length \code{nrow(df_X)}:
#'   \describe{
#'     \item{pmf_RSVMA}{Multiplier on overall RSV-episode risk.}
#'     \item{pmf_inpat}{Multiplier on hospitalisation probability given RSV.}
#'     \item{pmf_ICU}{Multiplier on ICU probability given hospitalisation.}
#'   }
#'   A value of 1 means no protection; values < 1 indicate partial protection.
#' @details
#'   When multiple vaccines overlap, the function takes the maximum efficacy
#'   (best protection) component-wise before converting to a probability
#'   modifier via \code{(1 - efficacy)}.
Vaccine_effect_modifier <- function(vaccine_efficacy_list, l_params, t, df_X) {
  n <- nrow(df_X)
  eff_RSVMA <- numeric(n)
  eff_inpat  <- numeric(n)
  eff_ICU    <- numeric(n)

  # Palivizumab — year 1
  idx_pvz1 <- !is.na(df_X$palivizumab_Month_y1) &
    df_X$palivizumab_Month_y1 <= t & t < 13L
  if (any(idx_pvz1)) {
    eff_RSVMA[idx_pvz1] <- pmax(eff_RSVMA[idx_pvz1], l_params$eff_RSVMA_palivizumab)
    eff_inpat[idx_pvz1] <- pmax(eff_inpat[idx_pvz1], vaccine_efficacy_list$eff_inpat_palivizumab)
    eff_ICU  [idx_pvz1] <- pmax(eff_ICU  [idx_pvz1], vaccine_efficacy_list$eff_ICU_palivizumab)
  }

  # Palivizumab — year 2
  idx_pvz2 <- !is.na(df_X$palivizumab_Month_y2) &
    df_X$palivizumab_Month_y2 <= t & t < 25L
  if (any(idx_pvz2)) {
    eff_RSVMA[idx_pvz2] <- pmax(eff_RSVMA[idx_pvz2], l_params$eff_RSVMA_palivizumab)
    eff_inpat[idx_pvz2] <- pmax(eff_inpat[idx_pvz2], vaccine_efficacy_list$eff_inpat_palivizumab)
    eff_ICU  [idx_pvz2] <- pmax(eff_ICU  [idx_pvz2], vaccine_efficacy_list$eff_ICU_palivizumab)
  }

  # RSVpreF
  idx_rsvf <- !is.na(df_X$RSVpreF_Month) & df_X$RSVpreF_Month <= t
  if (any(idx_rsvf)) {
    d <- as.integer(1 + t - df_X$RSVpreF_Month[idx_rsvf])
    eff_RSVMA[idx_rsvf] <- pmax(eff_RSVMA[idx_rsvf], l_params$eff_RSVMA_RSVpreF[d])
    eff_inpat[idx_rsvf] <- pmax(eff_inpat[idx_rsvf], vaccine_efficacy_list$eff_inpat_RSVpreF[d])
    eff_ICU  [idx_rsvf] <- pmax(eff_ICU  [idx_rsvf], vaccine_efficacy_list$eff_ICU_RSVpreF[d])
  }

  # Nirsevimab — year 1
  idx_nirs1 <- !is.na(df_X$nirsevimab_Month) & df_X$nirsevimab_Month <= t
  if (any(idx_nirs1)) {
    d <- as.integer(1 + t - df_X$nirsevimab_Month[idx_nirs1])
    eff_RSVMA[idx_nirs1] <- pmax(eff_RSVMA[idx_nirs1], l_params$eff_RSVMA_nirsevimab[d])
    eff_inpat[idx_nirs1] <- pmax(eff_inpat[idx_nirs1], vaccine_efficacy_list$eff_inpat_nirsevimab[d])
    eff_ICU  [idx_nirs1] <- pmax(eff_ICU  [idx_nirs1], vaccine_efficacy_list$eff_ICU_nirsevimab[d])
  }

  # Nirsevimab — year 2
  idx_nirs2 <- !is.na(df_X$nirsevimab_Month_y2) &
    df_X$nirsevimab_Month_y2 <= t & t < 25L
  if (any(idx_nirs2)) {
    d <- as.integer(1 + t - df_X$nirsevimab_Month_y2[idx_nirs2])
    eff_RSVMA[idx_nirs2] <- pmax(eff_RSVMA[idx_nirs2], l_params$eff_RSVMA_nirsevimab[d])
    eff_inpat[idx_nirs2] <- pmax(eff_inpat[idx_nirs2], vaccine_efficacy_list$eff_inpat_nirsevimab[d])
    eff_ICU  [idx_nirs2] <- pmax(eff_ICU  [idx_nirs2], vaccine_efficacy_list$eff_ICU_nirsevimab[d])
  }

  # Convert efficacies to probability multipliers.
  # pmf_inpat and pmf_ICU are conditional on the preceding level being avoided,
  # so they are bounded at 1.
  pmf_RSVMA <- 1 - eff_RSVMA
  pmf_inpat  <- pmin((1 - eff_inpat) / pmax(1 - eff_RSVMA, 1e-9), 1)
  pmf_ICU    <- pmin((1 - eff_ICU)   / pmax(1 - eff_inpat,  1e-9), 1)

  list(pmf_RSVMA = pmf_RSVMA, pmf_inpat = pmf_inpat, pmf_ICU = pmf_ICU)
}


# --- Treatment type ----------------------------------------------------------

#' Assign a treatment type to every infant in the current cycle.
#'
#' @param vaccine_efficacy_list Efficacy list from \code{build_vaccine_efficacy_list()}.
#' @param l_params              Model parameter list.
#' @param M_t                   Character vector of current health states.
#' @param t                     Integer cycle index.
#' @param df_X                  Cohort data frame.
#' @return Character vector of treatment types: one of \code{"NT"} (no
#'   treatment), \code{"ICU"}, \code{"PW"} (paediatric ward), \code{"ED"}
#'   (emergency department), or \code{"DV"} (doctor visit).
treatment_type <- function(vaccine_efficacy_list, l_params, M_t, t, df_X) {
  # Single call using l_params for RSVMA (SA-aware) and vaccine_efficacy_list
  # for inpat/ICU (not varied in SA)
  pmf <- Vaccine_effect_modifier(vaccine_efficacy_list, l_params, t, df_X)
  pmf_RSVMA <- pmf$pmf_RSVMA
  pmf_inpat  <- pmf$pmf_inpat
  pmf_ICU    <- pmf$pmf_ICU

  month_of_life <- ifelse(t >= df_X$MonthBorn, t - df_X$MonthBorn, NA_real_)

  Inp_prob <- pmin(
    ifelse(!is.na(month_of_life), prob_RSVMA_Inp[month_of_life + 1L], 0) * pmf_inpat,
    1
  )

  ICU_prob <- with(l_params, {
    p <- case_when(
      df_X$Risk == "HighRisk"     ~ Prob_ICU_HighRisk,
      df_X$Risk == "ModerateRisk" ~ Prob_ICU_ModerateRisk,
      TRUE                        ~ Prob_ICU_LowRisk
    )
    if (t > 12L) p <- p * RR_Prob_ICU_Second_season
    p * pmf_ICU
  })

  ED_prob <- ifelse(!is.na(month_of_life), prob_Outp_ED[month_of_life + 1L], 0)

  m_p_tr <- matrix(0, nrow = n_treatments, ncol = length(M_t),
                   dimnames = list(v_names_treatments, NULL))

  # Non-RSV states always receive no treatment
  for (s in c("NB", "H", "D")) {
    idx <- M_t == s
    m_p_tr["NT", idx] <- 1
  }

  idx <- M_t == "RSV"
  m_p_tr["NT",  idx] <- 1 - pmf_RSVMA[idx]
  m_p_tr["ICU", idx] <- pmf_RSVMA[idx] *      Inp_prob[idx] *      ICU_prob[idx]
  m_p_tr["PW",  idx] <- pmf_RSVMA[idx] *      Inp_prob[idx] * (1 - ICU_prob[idx])
  m_p_tr["ED",  idx] <- pmf_RSVMA[idx] * (1 - Inp_prob[idx]) *     ED_prob[idx]
  m_p_tr["DV",  idx] <- pmf_RSVMA[idx] * (1 - Inp_prob[idx]) * (1 - ED_prob[idx])

  samplev(t(m_p_tr), 1)
}


# --- Costs -------------------------------------------------------------------

#' Compute per-infant costs in the current cycle.
#'
#' @param l_params Model parameter list.
#' @param Tr_t     Character vector of treatment types.
#' @param t        Integer cycle index (used to identify vaccine administration
#'   months).
#' @param df_X     Cohort data frame.
#' @return A list with numeric vectors of length \code{nrow(df_X)}:
#'   \code{c_t} (total), \code{c_inpat} (inpatient only), \code{c_vac}
#'   (all vaccine costs combined), \code{c_vac_pvz} (palivizumab),
#'   \code{c_vac_nirs} (nirsevimab), \code{c_vac_rsvf} (RSVpreF).
Costs <- function(l_params, Tr_t, t, df_X) {
  with(as.list(l_params), {
    n <- length(Tr_t)
    c_inpat <- numeric(n)
    c_outpat <- numeric(n)

    if (any(Tr_t == "ICU")) {
      ICU_mask <- Tr_t == "ICU"
      LoS_ICU  <- case_when(
        df_X$Risk == "HighRisk"     ~ LOS_icu_HighRisk,
        df_X$Risk == "ModerateRisk" ~ LOS_icu_ModerateRisk,
        TRUE                        ~ LOS_icu_LowRisk
      )
      c_inpat[ICU_mask] <- (PD_ICU * LoS_ICU)[ICU_mask]
    }

    if (any(Tr_t == "PW")) {
      PW_mask <- Tr_t == "PW"
      LoS_PW  <- case_when(
        df_X$Risk == "HighRisk"     ~ LOS_hosp_HighRisk,
        df_X$Risk == "ModerateRisk" ~ LOS_hosp_ModerateRisk,
        TRUE                        ~ LOS_hosp_LowRisk
      )
      c_inpat[PW_mask] <- (PD_PW * LoS_PW)[PW_mask]
    }

    c_outpat <- ifelse(Tr_t == "ED", c_ED,
                ifelse(Tr_t == "DV", c_DV, 0))

    # Vaccine costs: charge in the cycle when each vaccine is administered
    vac_pvz1  <- ifelse(!is.na(df_X$palivizumab_Month_y1) &
                          df_X$palivizumab_Month_y1 == t,
                        c_palivizumab + c_Admin_palivizumab, 0)
    vac_pvz2  <- ifelse(!is.na(df_X$palivizumab_Month_y2) &
                          df_X$palivizumab_Month_y2 == t,
                        c_palivizumab + c_Admin_palivizumab, 0)
    vac_nirs1 <- ifelse(!is.na(df_X$nirsevimab_Month) &
                          df_X$nirsevimab_Month == t,
                        c_nirsevimab + c_Admin_vac, 0)
    vac_nirs2 <- ifelse(!is.na(df_X$nirsevimab_Month_y2) &
                          df_X$nirsevimab_Month_y2 == t,
                        c_nirsevimab + c_Admin_vac, 0)
    vac_rsvf  <- ifelse(!is.na(df_X$RSVpreF_Month) &
                          df_X$RSVpreF_Month == t,
                        c_RSVpreF + c_Admin_vac, 0)
    c_vac_pvz  <- vac_pvz1  + vac_pvz2
    c_vac_nirs <- vac_nirs1 + vac_nirs2
    c_vac_rsvf <- vac_rsvf
    c_vac      <- c_vac_pvz + c_vac_nirs + c_vac_rsvf

    list(c_t       = c_inpat + c_outpat + c_vac,
         c_inpat   = c_inpat,
         c_vac     = c_vac,
         c_vac_pvz = c_vac_pvz,
         c_vac_nirs= c_vac_nirs,
         c_vac_rsvf= c_vac_rsvf)
  })
}


# --- Utilities ---------------------------------------------------------------

#' Compute per-infant QALYs accrued in the current cycle.
#'
#' @param l_params Model parameter list (needs disutility values and
#'   \code{analytic_mort}).
#' @param M_t      Character vector of current health states.
#' @param Tr_t     Character vector of treatment types.
#' @param df_X     Cohort data frame (unused here but kept for API consistency).
#' @return Numeric vector of per-infant QALYs (always >= 0).
#' @details
#'   Base utility is 1 QALY/year × \code{cycle_length} for alive (H and RSV)
#'   states.  Treatment-specific disutilities are subtracted.  Caregiver
#'   disutility columns are retained in \code{l_params} for optional
#'   activation.
#'
#'   When \code{l_params$analytic_mort = TRUE}, the expected QALY loss from
#'   death is also subtracted analytically at each cycle:
#'   \itemize{
#'     \item Non-inpatient individuals: \code{p_mort × u_rs}
#'     \item ICU/PW individuals:        \code{p_mort_Inp × u_rs}
#'   }
#'   This replaces the end-of-simulation stochastic death penalty in
#'   \code{MicroSim()} and eliminates all Monte-Carlo noise from this rare
#'   event, producing smooth sensitivity-analysis curves.
Effs <- function(l_params, M_t, Tr_t, df_X) {
  with(as.list(l_params), {
    u_t <- numeric(length(M_t))
    u_t[M_t %in% c("H", "RSV")] <- u_H * cycle_length

    u_t[Tr_t == "ICU"] <- u_t[Tr_t == "ICU"] - Du_ICU_Infants
    u_t[Tr_t == "PW"]  <- u_t[Tr_t == "PW"]  - Du_PW_Infants
    u_t[Tr_t == "ED"]  <- u_t[Tr_t == "ED"]  - Du_ED_Infants
    u_t[Tr_t == "DV"]  <- u_t[Tr_t == "DV"]  - Du_DV_Infants

    if (isTRUE(analytic_mort)) {
      # Inpatient (ICU/PW) individuals face elevated mortality
      inpat <- M_t %in% c("H", "RSV") & Tr_t %in% c("ICU", "PW")
      u_t[inpat] <- u_t[inpat] - p_mort_Inp * u_rs

      # All other alive individuals face baseline mortality
      non_inpat <- M_t %in% c("H", "RSV") & !Tr_t %in% c("ICU", "PW")
      u_t[non_inpat] <- u_t[non_inpat] - p_mort * u_rs
    }

    pmax(u_t, 0)
  })
}
