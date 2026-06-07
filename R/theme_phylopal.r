#' A clean minimal theme for phyloPal plots
#'
#' A lightweight ggplot2 theme with clean axes and no background grid,
#' suitable for microbiome composition plots.
#'
#' @param base_size Base font size (default: 11)
#' @param base_family Base font family (default: "")
#'
#' @return A ggplot2 theme object
#'
#' @examples
#' \dontrun{
#' p + theme_phylopal()
#' }
#'
#' @export
theme_phylopal <- function(base_size = 11, base_family = "") {
  ggplot2::theme_bw(base_size = base_size, base_family = base_family) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(color = "black", linewidth = 0.4),
    axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.4),
    strip.background = ggplot2::element_rect(fill = NA, color = "black", linewidth = 0.4),
    legend.background = ggplot2::element_blank(),
    legend.key = ggplot2::element_blank(),
    legend.key.height = ggplot2::unit(10, "pt"),
    legend.key.width = ggplot2::unit(12, "pt")
  )
}