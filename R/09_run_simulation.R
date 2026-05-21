# =============================================================================
# 09_run_simulation.R
# Runs the base-case microsimulation across all strategies and replications,
# then computes incremental costs, QALYs, and naive ICERs vs. No Intervention.
#
# Outputs:
#   results_table  — data frame, also saved to Output/RSV_results.csv
# =============================================================================

#' Run the base-case simulation and return a results data frame.
#'
#' @param l_params              Model parameter list (must include
#'   \code{n_i}, \code{n_sim}, \code{seed}, and \code{base_str}).
#' @param vaccine_efficacy_list Time-varying efficacy list.
#' @param output_dir            Directory for CSV output.
#' @return Data frame with one row per strategy containing costs, QALYs,
#'   treatment counts, incremental values, and naive ICERs vs. the base
#'   strategy.
run_base_case <- function(l_params, vaccine_efficacy_list,
                          output_dir = "./Output") {

  n_i_run  <- l_params$n_i
  n_sim_run <- l_params$n_sim
  seed_run  <- l_params$seed

  acc <- function() setNames(numeric(n_str), v_names_str)
  v_tot_cost   <- acc(); v_vac_cost   <- acc(); v_inpat_cost <- acc()
  v_tot_qaly   <- acc(); v_ave_cost   <- acc(); v_ave_qaly   <- acc()
  v_ICU_count  <- acc(); v_PW_count   <- acc(); v_ED_count   <- acc()
  v_DV_count   <- acc(); v_D_count    <- acc()

  cat("Starting base-case simulation...\n")
  sim_time <- system.time({

    for (j in seq_len(n_sim_run)) {
      if (j %% max(1L, n_sim_run %/% 10L) == 0L)
        cat(sprintf("\r%d%% done", round(j / n_sim_run * 100)))

      set.seed(seed_run + j)
      df_X <- create_cohort(l_params)

      for (i in seq_len(n_str)) {
        out  <- MicroSim(l_params, vaccine_efficacy_list, df_X,
                         v_names_str[i], seed_run + j * 10L,
                         full_output = TRUE)
        inv_n <- 1 / n_sim_run
        v_tot_cost[i]  <- v_tot_cost[i]   + sum(out$tc)       * inv_n
        v_vac_cost[i]  <- v_vac_cost[i]   + sum(out$tc_vac)   * inv_n
        v_inpat_cost[i]<- v_inpat_cost[i] + sum(out$tc_inpat) * inv_n
        v_tot_qaly[i]  <- v_tot_qaly[i]   + sum(out$te)       * inv_n
        v_ave_cost[i]  <- v_ave_cost[i]   + mean(out$tc)      * inv_n
        v_ave_qaly[i]  <- v_ave_qaly[i]   + mean(out$te)      * inv_n
        v_ICU_count[i] <- v_ICU_count[i]  + out$ICU_count     * inv_n
        v_PW_count[i]  <- v_PW_count[i]   + out$PW_count      * inv_n
        v_ED_count[i]  <- v_ED_count[i]   + out$ED_count      * inv_n
        v_DV_count[i]  <- v_DV_count[i]   + out$DV_count      * inv_n
        v_D_count[i]   <- v_D_count[i]    + out$D_count       * inv_n
      }
    }
  })

  cat(sprintf("\nSimulation complete in %.1f seconds.\n", sim_time["elapsed"]))

  base_cost <- v_tot_cost[l_params$base_str]
  base_qaly <- v_tot_qaly[l_params$base_str]

  results_table <- data.frame(
    Strategy        = v_names_str,
    Total_Cost      = as.numeric(v_tot_cost),
    Vaccine_Cost    = as.numeric(v_vac_cost),
    Inpatient_Cost  = as.numeric(v_inpat_cost),
    ICU_Count       = as.numeric(v_ICU_count),
    PW_Count        = as.numeric(v_PW_count),
    ED_Count        = as.numeric(v_ED_count),
    DV_Count        = as.numeric(v_DV_count),
    D_Count         = as.numeric(v_D_count),
    Total_QALY      = as.numeric(v_tot_qaly),
    stringsAsFactors = FALSE
  )
  results_table$Incremental_Cost <- results_table$Total_Cost - base_cost
  results_table$Incremental_QALY <- results_table$Total_QALY - base_qaly
  results_table$ICER             <- results_table$Incremental_Cost /
                                    results_table$Incremental_QALY

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  out_path <- file.path(output_dir, "RSV_results.csv")
  write.csv(results_table, out_path, row.names = FALSE)
  cat("Results saved to:", out_path, "\n")

  results_table
}
