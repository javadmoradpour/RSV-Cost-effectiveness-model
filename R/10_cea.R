# =============================================================================
# 10_cea.R
# Cost-effectiveness analysis using the dampack package.
# Computes ICERs on the efficiency frontier, formats results, and produces
# the CE-plane plot.
# =============================================================================

#' Compute ICERs and plot the cost-effectiveness frontier.
#'
#' @param results_table Data frame returned by \code{run_base_case()}.
#' @param output_dir    Directory for PNG output.
#' @return A list with \code{df_cea} (raw dampack ICER table) and
#'   \code{table_cea} (formatted table).
run_cea <- function(results_table, output_dir = "./Output") {
  df_cea <- calculate_icers(
    cost       = results_table$Total_Cost,
    effect     = results_table$Total_QALY,
    strategies = results_table$Strategy
  )

  table_cea <- format_table_cea(df_cea)

  icer_plot <- plot(df_cea, label = "all", txtsize = 14) +
    expand_limits(x = max(table_cea$QALYs) + 0.1) +
    theme(legend.position = c(0.15, 0.85))

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  plot_path <- file.path(output_dir, "RSV_ICER_Plot.png")
  ggsave(filename = plot_path, plot = icer_plot, width = 8, height = 6, dpi = 300)
  cat("ICER plot saved to:", plot_path, "\n")

  cea_path <- file.path(output_dir, "RSV_CEA_table.csv")
  write.csv(table_cea, cea_path, row.names = FALSE)
  cat("CEA table saved to:", cea_path, "\n")

  list(df_cea = df_cea, table_cea = table_cea, icer_plot = icer_plot)
}
