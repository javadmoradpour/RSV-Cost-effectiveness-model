# =============================================================================
# 05_cohort.R
# Creates a cohort of individual-level characteristics used throughout the
# simulation.  Each row is one infant tracked over n_cycles monthly cycles.
# =============================================================================

#' Generate a cohort of infants with individual-level attributes.
#'
#' @param l_params Named list of model parameters.  Must contain
#'   \code{n_i} (cohort size) and \code{v_risk_prob} (named probability
#'   vector for risk levels).
#' @return A data frame with one row per infant and columns:
#'   \describe{
#'     \item{ID}{Integer index.}
#'     \item{MonthBorn}{Calendar-cycle in which the infant enters the model
#'       (1–12, sampled uniformly).}
#'     \item{M_init}{Initial health state — always \code{"NB"} (newborn).}
#'     \item{Risk}{Risk category: \code{"HighRisk"}, \code{"ModerateRisk"},
#'       or \code{"LowRisk"}.}
#'     \item{pvz}{1 if the infant is high-risk and born after the RSV season
#'       has started (palivizumab-eligible under current guidelines), 0
#'       otherwise.}
#'     \item{RSVpreF_Month, nirsevimab_Month, nirsevimab_Month_y2,
#'       palivizumab_Month_y1, palivizumab_Month_y2}{Cycle in which each
#'       vaccine is administered; \code{NA} if not vaccinated.}
#'   }
create_cohort <- function(l_params) {
  n       <- l_params$n_i
  v_monthborn <- sample(1:12, size = n, replace = TRUE)
  v_risk      <- sample(
    names(l_params$v_risk_prob), size = n, replace = TRUE,
    prob = l_params$v_risk_prob
  )
  v_pvz <- as.integer(v_risk == "HighRisk" & v_monthborn > sim_to_rsv_month_gap)

  data.frame(
    ID                   = seq_len(n),
    MonthBorn            = v_monthborn,
    M_init               = "NB",
    Risk                 = v_risk,
    pvz                  = v_pvz,
    RSVpreF_Month        = NA_real_,
    nirsevimab_Month     = NA_real_,
    nirsevimab_Month_y2  = NA_real_,
    palivizumab_Month_y1 = NA_real_,
    palivizumab_Month_y2 = NA_real_,
    stringsAsFactors     = FALSE
  )
}
