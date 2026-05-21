# =============================================================================
# 01_packages.R
# Load all required packages, installing from CRAN and GitHub if needed.
# =============================================================================

if (!require("pacman")) install.packages("pacman")
library(pacman)

p_load(
  "devtools", "dplyr", "scales", "ellipse", "ggplot2", "lazyeval",
  "igraph", "truncnorm", "ggraph", "reshape2", "knitr", "markdown",
  "stringr", "dampack", "matrixStats"
)
p_load_gh("DARTH-git/darthtools")
