# =============================================================================
# 06_strategies.R
# Functions that assign vaccination months to a cohort data frame according
# to each immunisation strategy.  Each function receives the cohort data frame
# and returns it with the appropriate *_Month columns filled in.
#
# apply_strategy() is the single dispatch point used by MicroSim().
# To add a new strategy: (1) write a helper function below, (2) register it
# in the strategy_fns list inside apply_strategy().
# =============================================================================

# --- Individual strategy functions -------------------------------------------

#' Palivizumab: high-risk infants in season 1; a fraction in season 2.
#'
#' @param df_X          Cohort data frame.
#' @param coverage_pvz_y2 Proportion of season-1 recipients also vaccinated in
#'   the second RSV season.
#' @return Modified cohort data frame.
Palivizumab <- function(df_X, coverage_pvz_y2) {
  eligible_y1 <- df_X$pvz == 1L
  df_X$palivizumab_Month_y1[eligible_y1] <- df_X$MonthBorn[eligible_y1]

  eligible_y2 <- !is.na(df_X$palivizumab_Month_y1) &
    runif(nrow(df_X)) <= coverage_pvz_y2
  df_X$palivizumab_Month_y2[eligible_y2] <- 12L + (sim_to_rsv_month_gap + 1L)

  df_X
}


#' Nirsevimab replaces palivizumab for palivizumab-eligible infants.
#'
#' @param df_X            Cohort data frame.
#' @param coverage_pvz_y2 Second-season coverage probability.
#' @return Modified cohort data frame.
Nirsevimab_Replacement <- function(df_X, coverage_pvz_y2) {
  eligible_y1 <- df_X$pvz == 1L
  df_X$nirsevimab_Month[eligible_y1] <- df_X$MonthBorn[eligible_y1]

  eligible_y2 <- !is.na(df_X$nirsevimab_Month) &
    runif(nrow(df_X)) <= coverage_pvz_y2
  df_X$nirsevimab_Month_y2[eligible_y2] <- 12L + (sim_to_rsv_month_gap + 1L)

  df_X
}


#' Nirsevimab for all high-risk infants with catch-up at RSV-season start.
#'
#' Infants born before the RSV season receive nirsevimab at season start;
#' those born during the season receive it at birth.
#'
#' @param df_X                        Cohort data frame.
#' @param coverage_nirsevimab_HighRisk Uptake proportion (default 1).
#' @param coverage_pvz_y2             Second-season palivizumab-replacement
#'   coverage (passed to \code{Nirsevimab_Replacement}).
#' @return Modified cohort data frame.
Nirsevimab_HighRisk_WithCatchup <- function(df_X, coverage_nirsevimab_HighRisk,
                                            coverage_pvz_y2) {
  df_X <- Nirsevimab_Replacement(df_X, coverage_pvz_y2)

  eligible <- df_X$Risk == "HighRisk" &
    is.na(df_X$nirsevimab_Month) &
    runif(nrow(df_X)) <= coverage_nirsevimab_HighRisk

  df_X$nirsevimab_Month[eligible] <- ifelse(
    df_X$MonthBorn[eligible] <= sim_to_rsv_month_gap,
    sim_to_rsv_month_gap + 1L,
    df_X$MonthBorn[eligible]
  )
  df_X
}


#' Nirsevimab for moderate-risk infants with catch-up at RSV-season start.
#'
#' Infants born before the RSV season and not already protected by RSVpreF
#' receive nirsevimab at season start; those born during the season receive it
#' at birth.
#'
#' @param df_X                           Cohort data frame.
#' @param coverage_nirsevimab_ModerateRisk Uptake proportion.
#' @return Modified cohort data frame.
Nirsevimab_ModerateRisk_WithCatchup <- function(df_X,
                                                 coverage_nirsevimab_ModerateRisk) {
  eligible <- df_X$Risk == "ModerateRisk" &
    is.na(df_X$RSVpreF_Month) &
    runif(nrow(df_X)) <= coverage_nirsevimab_ModerateRisk

  df_X$nirsevimab_Month[eligible] <- ifelse(
    df_X$MonthBorn[eligible] <= sim_to_rsv_month_gap,
    sim_to_rsv_month_gap + 1L,
    df_X$MonthBorn[eligible]
  )
  df_X
}


#' Nirsevimab for low-risk infants born during the RSV season (no catch-up).
#'
#' Only low-risk infants born after the RSV season starts are vaccinated at
#' birth.  Moderate-risk infants are handled separately by
#' \code{Nirsevimab_ModerateRisk_WithCatchup()}, which includes catch-up.
#'
#' @param df_X                           Cohort data frame.
#' @param coverage_nirsevimab_ModerateRisk Uptake proportion.
#' @return Modified cohort data frame.
Nirsevimab_LowRisk_InSeason <- function(df_X, coverage_nirsevimab_ModerateRisk) {
  eligible <- df_X$Risk == "LowRisk" &
    df_X$MonthBorn > sim_to_rsv_month_gap &
    runif(nrow(df_X)) <= coverage_nirsevimab_ModerateRisk

  df_X$nirsevimab_Month[eligible] <- df_X$MonthBorn[eligible]
  df_X
}


#' RSVpreF maternal vaccine for infants born on or after the RSVpreF offer date.
#'
#' The offer window is controlled by \code{rsvpref_gap} (cycles from simulation
#' start to the first month RSVpreF is offered).  By default this equals
#' \code{sim_to_rsv_month_gap} (RSV season start), but it can differ when
#' \code{l_params$rsvpref_start_month} is changed independently.
#'
#' @param df_X            Cohort data frame.
#' @param coverage_RSVpreF Uptake proportion.
#' @param rsvpref_gap     Number of cycles from simulation start before RSVpreF
#'   is offered (\code{rsvpref_start_month - sim_start_month}).
#' @return Modified cohort data frame.
RSVpreF_RSVSeason <- function(df_X, coverage_RSVpreF, rsvpref_gap) {
  eligible <- df_X$Risk %in% c("ModerateRisk", "LowRisk") &
    df_X$MonthBorn > rsvpref_gap &
    runif(nrow(df_X)) <= coverage_RSVpreF

  df_X$RSVpreF_Month[eligible] <- df_X$MonthBorn[eligible]
  df_X
}


# --- Strategy dispatch -------------------------------------------------------

#' Apply a named immunisation strategy to a cohort.
#'
#' @param df_X     Cohort data frame produced by \code{create_cohort()}.
#' @param Str      Strategy name — must match one of \code{v_names_str}.
#' @param l_params Model parameter list (provides coverage values).
#' @return Cohort data frame with vaccination-month columns filled in.
#' @details
#'   All coverage parameters are extracted from \code{l_params} so that
#'   sensitivity-analysis runs with modified parameters are handled correctly
#'   without relying on global state.  Add new strategies by extending the
#'   \code{strategy_fns} list.
apply_strategy <- function(df_X, Str, l_params) {
  pvz_y2      <- l_params$coverage_palivizumab_y2
  nirs_mod    <- l_params$coverage_nirsevimab_ModerateRisk
  rsvf_cov    <- l_params$coverage_RSVpreF
  rsvpref_gap <- l_params$rsvpref_start_month - sim_start_month

  strategy_fns <- list(
    "No Intervention" = function(df) df,

    "Palivizumab" = function(df) {
      Palivizumab(df, pvz_y2)
    },

    "Nirsevimab (High & Mod)" = function(df) {
      df <- Nirsevimab_HighRisk_WithCatchup(df, coverage_nirsevimab_HighRisk,
                                            pvz_y2)
      Nirsevimab_ModerateRisk_WithCatchup(df, nirs_mod)
    },

    "RSVpreF" = function(df) {
      RSVpreF_RSVSeason(df, rsvf_cov, rsvpref_gap)
    },

    "RSVpreF + Nirsevimab (High)" = function(df) {
      df <- RSVpreF_RSVSeason(df, rsvf_cov, rsvpref_gap)
      Nirsevimab_HighRisk_WithCatchup(df, coverage_nirsevimab_HighRisk, pvz_y2)
    },

    "RSVpreF + Nirsevimab (High & Mod)" = function(df) {
      df <- RSVpreF_RSVSeason(df, rsvf_cov, rsvpref_gap)
      df <- Nirsevimab_HighRisk_WithCatchup(df, coverage_nirsevimab_HighRisk,
                                            pvz_y2)
      Nirsevimab_ModerateRisk_WithCatchup(df, nirs_mod)
    },

    "Nirsevimab (Universal)" = function(df) {
      df <- Nirsevimab_HighRisk_WithCatchup(df, coverage_nirsevimab_HighRisk,
                                            pvz_y2)
      df <- Nirsevimab_ModerateRisk_WithCatchup(df, nirs_mod)
      Nirsevimab_LowRisk_InSeason(df, nirs_mod)
    }
  )

  if (!Str %in% names(strategy_fns)) {
    stop(sprintf("Unknown strategy: '%s'.\nAvailable: %s",
                 Str, paste(names(strategy_fns), collapse = ", ")))
  }
  strategy_fns[[Str]](df_X)
}
