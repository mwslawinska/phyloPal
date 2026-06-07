#' @importFrom magrittr %>%
#' @importFrom rlang .data :=
#' @importFrom stats formula
#' @importFrom utils combn modifyList
#' @importFrom ggplot2 after_stat
NULL

# Suppress R CMD check notes for variables used in dplyr/ggplot2 pipelines
utils::globalVariables(c(
  ".", "RA", "ASV_ID", "Abundance", "abundance_status",
  "tax_val", "tax_type", "tax_color_plot", "stratum",
  "label", "x", "xend", "y", "yend",
  ".orig", ".denom", ".group", ".tax", ".lab",
  "tax_color"
))