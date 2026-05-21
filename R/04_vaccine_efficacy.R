# =============================================================================
# 04_vaccine_efficacy.R
# Functions to compute time-varying vaccine efficacy profiles and to build
# the efficacy list used by Vaccine_effect_modifier().
# =============================================================================

#' Compute a time-varying vaccine efficacy vector over the simulation horizon.
#'
#' @param efficacy    Scalar peak efficacy (0–1).
#' @param valid_months Integer number of months at peak efficacy before decline
#'   begins.
#' @param zero_month  Integer month at which efficacy reaches zero.
#' @return Numeric vector of length \code{n_cycles} with per-cycle efficacy.
#' @details
#'   Efficacy is constant at \code{efficacy} for the first \code{valid_months}
#'   cycles, then declines linearly to zero at \code{zero_month}, and remains
#'   zero thereafter.
vaccine_efficacy <- function(efficacy, valid_months, zero_month = zero_month) {
  eff_vec <- numeric(n_cycles)

  end_flat <- min(valid_months, n_cycles)
  eff_vec[seq_len(end_flat)] <- efficacy

  if (n_cycles > valid_months && valid_months < zero_month) {
    decline_end    <- min(zero_month, n_cycles)
    decline_months <- (valid_months + 1):decline_end
    eff_vec[decline_months] <- efficacy *
      (1 - (decline_months - valid_months) / (zero_month - valid_months))
  }

  eff_vec
}


#' Build the list of time-varying efficacy vectors passed to
#' \code{Vaccine_effect_modifier()}.
#'
#' @param l_params Named list of model parameters (from \code{build_params()}).
#' @return A named list with per-cycle efficacy vectors for RSVpreF and
#'   nirsevimab (inpatient and ICU components) plus scalar constants for
#'   palivizumab.  RSVMA vectors are stored in \code{l_params} rather than
#'   here because they are varied in sensitivity analyses.
build_vaccine_efficacy_list <- function(l_params) {
  list(
    # Palivizumab: constant efficacy (assumed throughout the RSV season)
    eff_inpat_palivizumab = l_params$eff_inpat_palivizumab,
    eff_ICU_palivizumab   = l_params$eff_ICU_palivizumab,

    # RSVpreF: time-varying inpatient and ICU efficacy
    eff_inpat_RSVpreF = vaccine_efficacy(l_params$efficacy_inpat_RSVpreF,
                                         valid_months_RSVpreF, zero_month),
    eff_ICU_RSVpreF   = vaccine_efficacy(l_params$efficacy_ICU_RSVpreF,
                                         valid_months_RSVpreF, zero_month),

    # Nirsevimab: time-varying inpatient and ICU efficacy
    eff_inpat_nirsevimab = vaccine_efficacy(l_params$efficacy_inpat_nirsevimab,
                                             valid_months_nirsevimab, zero_month),
    eff_ICU_nirsevimab   = vaccine_efficacy(l_params$efficacy_ICU_nirsevimab,
                                             valid_months_nirsevimab, zero_month)
  )
}


#' Add time-varying RSVMA efficacy vectors to l_params.
#'
#' @param l_params Named list of model parameters.
#' @return The same list with \code{eff_RSVMA_RSVpreF} and
#'   \code{eff_RSVMA_nirsevimab} vectors appended.  These are kept in
#'   \code{l_params} (not in \code{vaccine_efficacy_list}) so that
#'   sensitivity-analysis code can modify the underlying scalar and
#'   regenerate the vector without rebuilding the full efficacy list.
add_rsvma_vectors <- function(l_params) {
  l_params$eff_RSVMA_RSVpreF <- vaccine_efficacy(
    l_params$efficacy_RSVMA_RSVpreF, valid_months_RSVpreF, zero_month
  )
  l_params$eff_RSVMA_nirsevimab <- vaccine_efficacy(
    l_params$efficacy_RSVMA_nirsevimab, valid_months_nirsevimab, zero_month
  )
  l_params
}
