# =============================================================================
# 08_microsimulation.R
# The core individual-level simulation loop.
#
# MicroSim() runs one strategy for one Monte Carlo replication.  It accepts a
# `full_output` flag so that the deterministic sensitivity-analysis wrapper
# (calculate_ce_out()) can skip tracking treatment counts and separate cost
# components, reducing both computation and memory during the many SA runs.
# =============================================================================

#' Run the microsimulation for a single immunisation strategy.
#'
#' @param l_params              Model parameter list (from \code{build_params()}
#'   augmented by \code{add_rsvma_vectors()}).
#' @param vaccine_efficacy_list Time-varying efficacy list (from
#'   \code{build_vaccine_efficacy_list()}).
#' @param df_X                  Cohort data frame (from \code{create_cohort()}).
#' @param Str                   Strategy name (must match \code{v_names_str}).
#' @param seed                  Random seed for reproducibility.
#' @param full_output           Logical.  If \code{TRUE} (default) also returns
#'   disaggregated vaccine/inpatient costs and treatment counts.  Set to
#'   \code{FALSE} for sensitivity-analysis runs where only total costs and
#'   QALYs are needed.
#' @return A named list with:
#'   \describe{
#'     \item{tc}{Numeric vector (length \code{nrow(df_X)}) of total per-person
#'       costs.}
#'     \item{te}{Numeric vector of per-person QALYs.}
#'     \item{tc_vac, tc_inpat}{Disaggregated cost vectors (only when
#'       \code{full_output = TRUE}).}
#'     \item{ICU_count, PW_count, ED_count, DV_count, D_count}{Aggregate
#'       treatment/death counts across the cohort and all cycles (only when
#'       \code{full_output = TRUE}).}
#'   }
MicroSim <- function(l_params, vaccine_efficacy_list, df_X, Str, seed,
                     full_output = TRUE) {
  set.seed(seed)

  n_i  <- nrow(df_X)
  df_X <- apply_strategy(df_X, Str, l_params)

  # Pre-compute the parts of the hospitalisation LP that are fixed throughout
  # the simulation (individual risk/pvz) or fixed per cycle (calendar month).
  lp_static   <- precompute_lp_static(df_X, l_params$hosp_model_beta)
  lp_calendar <- precompute_lp_calendar(l_params$hosp_model_beta)

  # Pre-allocate result matrices (n_i × n_cycles+1)
  col_names <- paste("cycle", 0:n_cycles)
  row_names <- paste("ind",   seq_len(n_i))

  m_M  <- matrix(NA_character_, n_i, n_cycles + 1L,
                 dimnames = list(row_names, col_names))
  m_Tr <- matrix(NA_character_, n_i, n_cycles + 1L,
                 dimnames = list(row_names, col_names))
  m_C  <- matrix(0,            n_i, n_cycles + 1L,
                 dimnames = list(row_names, col_names))
  m_E  <- matrix(0,            n_i, n_cycles + 1L,
                 dimnames = list(row_names, col_names))

  if (full_output) {
    m_C_vac   <- matrix(0, n_i, n_cycles + 1L,
                        dimnames = list(row_names, col_names))
    m_C_inpat <- matrix(0, n_i, n_cycles + 1L,
                        dimnames = list(row_names, col_names))
  }

  # Initialise cycle 0
  m_M[, 1L]  <- as.character(df_X$M_init)
  m_Tr[, 1L] <- treatment_type(vaccine_efficacy_list, l_params, m_M[, 1L], 0L, df_X)

  cost_0     <- Costs(l_params, m_Tr[, 1L], 0L, df_X)
  m_C[, 1L]  <- cost_0$c_t
  m_E[, 1L]  <- Effs(l_params, m_M[, 1L], m_Tr[, 1L], df_X)
  if (full_output) {
    m_C_vac[, 1L]   <- cost_0$c_vac
    m_C_inpat[, 1L] <- cost_0$c_inpat
  }

  set.seed(seed + 1L)

  for (t in seq_len(n_cycles)) {
    m_M[, t + 1L]  <- state_transition(l_params, m_M[, t], m_Tr[, t], t,
                                        df_X, lp_static, lp_calendar)
    m_Tr[, t + 1L] <- treatment_type(vaccine_efficacy_list, l_params,
                                      m_M[, t + 1L], t, df_X)

    cost_t         <- Costs(l_params, m_Tr[, t + 1L], t, df_X)
    m_C[, t + 1L]  <- cost_t$c_t
    m_E[, t + 1L]  <- Effs(l_params, m_M[, t + 1L], m_Tr[, t + 1L], df_X)

    if (full_output) {
      m_C_vac[, t + 1L]   <- cost_t$c_vac
      m_C_inpat[, t + 1L] <- cost_t$c_inpat
    }
  }

  tc <- rowSums(m_C)
  te <- rowSums(m_E)

  # Subtract quality-of-life loss from premature death.
  # When analytic_mort = TRUE the penalty has already been applied as an
  # expected value in every cycle inside Effs(), so we skip the stochastic
  # end-of-simulation adjustment to avoid double-counting.
  if (!isTRUE(l_params$analytic_mort)) {
    v_du_D <- ifelse(m_M[, n_cycles + 1L] == "D", u_rs, 0)
    te     <- te - v_du_D
  }

  if (!full_output) {
    return(list(tc = tc, te = te))
  }

  # Aggregate treatment counts across all cycles (excluding cycle-0 column)
  tr_counts <- table(factor(as.vector(m_Tr[, -1L]),
                            levels = v_names_treatments))

  list(
    tc        = tc,
    tc_vac    = rowSums(m_C_vac),
    tc_inpat  = rowSums(m_C_inpat),
    te        = te,
    ICU_count = as.integer(tr_counts["ICU"]),
    PW_count  = as.integer(tr_counts["PW"]),
    ED_count  = as.integer(tr_counts["ED"]),
    DV_count  = as.integer(tr_counts["DV"]),
    D_count   = sum(m_M[, n_cycles + 1L] == "D")
  )
}
