#' Create an alluvial (Sankey) plot of taxonomic composition across groups
#'
#' Draws an alluvial plot where strata represent taxa (or taxon labels) within each group.
#' By default, strata are drawn for each \code{tax_val}. Optionally, strata can be collapsed into
#' one bucket per group (e.g. for \code{"unique low abundant"} taxa) to improve readability and
#' avoid numerical issues when many extremely small strata are present.
#' 
#' @param data A data frame (typically output of \code{\link{classify_taxa_patterns}}) containing:
#'   \itemize{
#'     \item \code{RA}: relative abundance per \code{group_col} (0--1; typically sums to 1 within group),
#'     \item \code{tax_val}: displayed stratum label (often the taxon name),
#'     \item \code{tax_color}: palette key used for color mapping (may equal \code{tax_val} for most taxa),
#'     \item the grouping column named by \code{group_col}.
#'   }
#' @param custom_palette Named character vector of colors. Names must match palette keys in \code{data$tax_color}.
#'   The order of \code{names(custom_palette)} defines legend order and (if \code{order_strata_by_palette=TRUE})
#'   influences stratum stacking order.
#' @param tax_level Character string used for the legend title (e.g. \code{"Family"}). Does not affect aggregation.
#' @param group_col Character string naming the grouping column used on the x-axis.
#' @param alpha Numeric transparency for flows (default: 0.7).
#' @param stratum_width Numeric width of strata passed to \code{ggalluvial::geom_stratum()} (default: 1/3).
#' @param flow_width Numeric width of flows passed to \code{ggalluvial::geom_flow()} (default: 1/3).
#' @param add_labels Logical; if \code{TRUE}, add stratum labels (default: FALSE).
#' @param label_size Numeric text size for labels (default: 3).
#' @param theme_obj A ggplot2 theme object (default: \code{ggplot2::theme_minimal()}).
#' @param line_width Numeric line width for stratum borders (default: 0.2).
#' @param legend_position Legend position (default: \code{"right"}).
#' @param y_axis_label y-axis label (default: \code{"Relative Abundance"}).
#' @param x_axis_label x-axis label (default: empty string).
#' @param order_strata_by_palette Logical; if \code{TRUE} (default), strata are ordered by palette-key order
#'   (via \code{tax_color}). This makes stacking and legend consistent across plots.
#' @param within_key_sort Character; when \code{order_strata_by_palette=TRUE}, controls how \code{tax_val} are ordered
#'   within each palette key. \code{"alphabetical"} (default) sorts taxa within each key; \code{"none"} keeps taxa in the
#'   order encountered within each key after key-ordering is applied.
#' @param collapse_strata Logical; if TRUE, collapse selected palette-key bins into a single stratum per group.
#'   This is useful when many taxa share one palette key (e.g. \code{"unique low abundant"}) and would otherwise
#'   generate many microscopic strata. Default: FALSE.
#' @param collapse_keys Character vector of palette keys to collapse when \code{collapse_strata=TRUE}.
#'   Default: \code{c("unique low abundant","shared low abundant","unknown")}.
#' @param bucket_small Logical; if TRUE, taxa with \code{RA < bucket_min_ra} are collapsed into one bucket per group.
#'   This is independent of classification thresholds and is intended as a plotting stabilization/readability option.
#'   Default: FALSE.
#' @param bucket_min_ra Numeric in (0,1); taxa with \code{RA < bucket_min_ra} are put into \code{bucket_label}.
#'   Default: 0.001 (0.1%).
#' @param bucket_label Character label used for the bucket created by \code{bucket_small}. Must be present in the palette
#'   (or will be appended as an extra key and drawn with \code{na.value}). Default: \code{"very low abundant"}.
#' @param bucket_has_flows Logical; if FALSE (default), the bucket created by \code{bucket_small} will not have flows.
#'   If TRUE, flows are drawn for that bucket like any other shared stratum.
#'
#' @section Ordering behavior (palette-driven workflow):
#' \itemize{
#'   \item \strong{Legend order} follows \code{names(custom_palette)}.
#'   \item If \code{order_strata_by_palette=TRUE}, \strong{stacking order} is driven by palette-key order
#'         (\code{tax_color}). This is useful when multiple taxa share the same palette key (e.g. low-abundance bins).
#' }
#'
#' @section Flow behavior:
#' Flows are drawn only for taxa (\code{tax_val}) that occur in more than one group
#' (i.e. present in \eqn{>1} distinct \code{group_col} levels). Taxa unique to a single group are shown as strata only.
#'
#' @section Palette coverage:
#' \code{custom_palette} should ideally include a color for every value of \code{data$tax_color}.
#' Missing palette keys are displayed using \code{na.value} (default \code{"grey50"}).
#'
#' @section Stratum identity vs. color identity:
#' \itemize{
#'   \item \code{tax_val} is the taxon label used for strata by default (one stratum per taxon).
#'   \item \code{tax_color} is the palette key controlling fill colors and legend order (multiple taxa may share one key).
#'   \item If \code{collapse_strata=TRUE} and/or \code{bucket_small=TRUE}, the plot uses a derived stratum ID so that
#'         selected bins (e.g. \code{"unique low abundant"}) or tiny taxa are drawn as a single bucket per group.
#' }
#' 
#' @return A \code{ggplot} object.
#'
#' @seealso \code{\link{prepare_alluvial_data}}, \code{\link{classify_taxa_patterns}},
#'   \code{\link{generate_alluvial_palette}}, \code{\link{create_alluvial_plot}}
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' library(ggplot2)
#'
#' toy <- tibble::tibble(
#'   Group = rep(c("G1","G2"), each = 3),
#'   tax_val   = rep(c("A","B","C"), times = 2),
#'   RA        = c(0.7, 0.2, 0.1,  0.6, 0.02, 0.38),
#'   tax_color = c("A","B","C", "A","unique low abundant","C")
#' )
#'
#' pal <- c(
#'   "unique low abundant" = "#A5A5A5",
#'   "A" = "#1b9e77",
#'   "B" = "#d95f02",
#'   "C" = "#7570b3"
#' )
#'
#' plot_alluvial(toy, pal, tax_level = "Class", group_col = "Group")
#'}
#' @export
plot_alluvial <- function(data,
                          custom_palette,
                          tax_level,
                          group_col = "Group",
                          alpha = 0.7,
                          stratum_width = 1/3,
                          flow_width = 1/3,
                          add_labels = FALSE,
                          label_size = 3,
                          theme_obj = ggplot2::theme_minimal(),
                          line_width = 0.2,
                          legend_position = "right",
                          y_axis_label = "Relative Abundance",
                          x_axis_label = "",
                          order_strata_by_palette = TRUE,
                          within_key_sort = c("alphabetical", "none"),
                          collapse_strata = FALSE,
                          collapse_keys = c("unique low abundant", "shared low abundant", "unknown"),
                          bucket_small = FALSE,
                          bucket_min_ra = 0.001,
                          bucket_label = "very low abundant",
                          bucket_has_flows = FALSE) {

  within_key_sort <- match.arg(within_key_sort)

  # ---------------------------
  # Input validation
  # ---------------------------
  if (!"tax_color" %in% colnames(data)) stop("Data must contain 'tax_color' column (from classify_taxa_patterns)")
  if (!"tax_val" %in% colnames(data)) stop("Data must contain 'tax_val' column")
  if (!"RA" %in% colnames(data)) stop("Data must contain 'RA' column")
  if (!group_col %in% colnames(data)) stop("group_col '", group_col, "' not found in data")
  if (!is.numeric(bucket_min_ra) || length(bucket_min_ra) != 1 || bucket_min_ra <= 0 || bucket_min_ra >= 1) {
    stop("bucket_min_ra must be a single numeric value in (0,1).")
  }

  palette_names <- names(custom_palette)
  if (is.null(palette_names) || anyNA(palette_names) || any(palette_names == "")) {
    stop("custom_palette must be a *named* vector (names are palette keys).")
  }

  # ---------------------------
  # Palette-driven legend order (tax_color)
  # ---------------------------
  tax_color_chr <- as.character(data$tax_color)
  extras <- setdiff(unique(tax_color_chr), palette_names)

  # If we introduce bucket_label, ensure it appears as a palette key level (color may come from na.value if missing)
  if (isTRUE(bucket_small)) extras <- unique(c(extras, bucket_label))

  palette_levels <- c(palette_names, setdiff(extras, palette_names))
  data$tax_color <- factor(tax_color_chr, levels = palette_levels)

  # X order
  data[[group_col]] <- factor(data[[group_col]], levels = unique(data[[group_col]]))

  # ---------------------------
  # NEW: define plotting stratum and fill keys (and optionally collapse)
  # ---------------------------
  tax_val_chr <- as.character(data$tax_val)
  tax_col_chr <- as.character(data$tax_color)

  # (A) collapse selected palette-key bins into one stratum
  stratum_plot <- tax_val_chr
  if (isTRUE(collapse_strata)) {
    stratum_plot <- ifelse(tax_col_chr %in% collapse_keys, tax_col_chr, tax_val_chr)
  }

  # (B) bucket very small taxa by RA cutoff
  bucket_mask <- rep(FALSE, nrow(data))
  if (isTRUE(bucket_small)) {
    bucket_mask <- data$RA < bucket_min_ra
    stratum_plot <- ifelse(bucket_mask, bucket_label, stratum_plot)
    tax_col_chr  <- ifelse(bucket_mask, bucket_label, tax_col_chr)
  }

  data$stratum_plot <- stratum_plot
  data$tax_color_plot <- factor(tax_col_chr, levels = palette_levels)

  # ---------------------------
  # If we collapsed strata / bucketed, we must aggregate to unique (x, stratum)
  # ---------------------------
  class_col <- if ("category" %in% names(data)) "category" else if ("tax_type" %in% names(data)) "tax_type" else NA_character_
  if (is.na(class_col)) {
    stop("plot_alluvial() needs a 'category' or 'tax_type' column to decide which taxa should have flows.")
  }

  data_plot <- data %>%
    dplyr::group_by(.data[[group_col]], stratum_plot, tax_color_plot, .data[[class_col]]) %>%
    dplyr::summarise(RA = sum(RA, na.rm = TRUE), .groups = "drop") %>%
    dplyr::group_by(.data[[group_col]]) %>%
    dplyr::mutate(RA = RA / sum(RA, na.rm = TRUE)) %>%  # keep bars at 1 after aggregation
    dplyr::ungroup()

  # ---------------------------
  # Palette-driven stacking order (stratum_plot) via tax_color_plot
  # ---------------------------
  if (isTRUE(order_strata_by_palette)) {
    tmp <- data_plot %>%
      dplyr::distinct(tax_color_plot = as.character(.data$tax_color_plot),
                      stratum_plot = as.character(.data$stratum_plot)) %>%
      dplyr::mutate(tax_color_plot = factor(tax_color_plot, levels = palette_levels))

    if (within_key_sort == "alphabetical") {
      tmp <- tmp %>% dplyr::arrange(tax_color_plot, stratum_plot)
    } else {
      tmp <- tmp %>% dplyr::arrange(tax_color_plot)
    }

    stratum_levels <- unique(tmp$stratum_plot)
    data_plot$stratum_plot <- factor(as.character(data_plot$stratum_plot), levels = stratum_levels)
  }

  # ---------------------------
  # Shared strata (flows): keep your semantic rule, but apply to stratum_plot
  # ---------------------------
  shared_levels <- c("shared abundant", "shared low abundant", "shared mixed abundance")

  shared_data <- data_plot %>%
    dplyr::filter(
      .data[[class_col]] %in% shared_levels,
      !is.na(stratum_plot),
      as.character(stratum_plot) != "unknown"
    )

  # By default: do not draw flows for the bucket created by bucket_small
  if (isTRUE(bucket_small) && !isTRUE(bucket_has_flows)) {
    shared_data <- shared_data %>% dplyr::filter(as.character(stratum_plot) != bucket_label)
  }

  # ---------------------------
  # Plot
  # ---------------------------
  p <- ggplot2::ggplot(data_plot) +

    ggalluvial::geom_stratum(
      ggplot2::aes(
        x = .data[[group_col]],
        stratum = stratum_plot,
        y = RA,
        fill = tax_color_plot
      ),
      width = stratum_width,
      color = "black",
      linewidth = line_width,
      alpha = 1
    ) +

    {
      if (nrow(shared_data) > 0)
        ggalluvial::geom_flow(
          data = shared_data,
          ggplot2::aes(
            x = .data[[group_col]],
            stratum = stratum_plot,
            alluvium = stratum_plot,
            y = RA,
            fill = tax_color_plot
          ),
          width = flow_width,
          alpha = alpha
        )
    } +

    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::scale_fill_manual(
      values = custom_palette,
      breaks = palette_levels,
      drop = TRUE,
      na.value = "grey50"
    ) +
    ggplot2::labs(
      y = y_axis_label,
      x = x_axis_label,
      fill = tax_level
    ) +
    theme_obj +
    ggplot2::theme(legend.position = legend_position)

  if (isTRUE(add_labels)) {
    p <- p +
      ggplot2::geom_text(
        stat = "stratum",
        ggplot2::aes(label = after_stat(stratum)),
        size = label_size
      )
  }

  return(p)
}


#' Create an alluvial plot via the full workflow wrapper
#'
#' Runs the alluvial workflow end-to-end:
#' \enumerate{
#'   \item \code{\link{prepare_alluvial_data}}: aggregate ASV-level \code{RA} to \code{tax_level} and compute group means,
#'   \item \code{\link{classify_taxa_patterns}}: assign abundance-pattern classes and create \code{tax_color},
#'   \item \code{\link{generate_alluvial_palette}}: generate a palette for \code{tax_color} (skipped if \code{custom_palette} is provided),
#'   \item \code{\link{plot_alluvial}}: draw the plot.
#' }
#'
#' Additional arguments are routed to the underlying functions via the \code{*_args} lists.
#'
#' @param data Raw long-format data with ASV-level \code{RA} and taxonomy columns.
#' @param tax_level Character string naming the taxonomy rank to plot (e.g. \code{"Class"}, \code{"Family"}).
#' @param group_col Character string naming the grouping column used on the x-axis.
#' @param SampleID_col Character string naming the sample ID column (default: \code{"SampleID"}).
#' @param groups Optional character vector restricting groups included. If \code{NULL}, all groups are used.
#' @param low_abundance_threshold Numeric threshold passed to \code{\link{classify_taxa_patterns}} (default: 0.01).
#' @param palette_list Optional vector of HCL palette names passed to \code{\link{generate_alluvial_palette}}.
#' @param custom_palette Optional named palette mapping \code{tax_color} keys to colors. If provided, palette generation is skipped.
#' @param prepare_args Named list of extra arguments passed to \code{\link{prepare_alluvial_data}}.
#' @param classify_args Named list of extra arguments passed to \code{\link{classify_taxa_patterns}}.
#' @param palette_args Named list of extra arguments passed to \code{\link{generate_alluvial_palette}}
#'   (only used if \code{custom_palette} is \code{NULL}).
#' @param plot_args Named list of extra arguments passed to \code{\link{plot_alluvial}}.
#' @param return_all Logical; if \code{TRUE}, returns a list with intermediate objects:
#'   \code{data_prepared}, \code{data_classified}, \code{palette}, \code{plot}.
#'
#' @return A \code{ggplot} object, or a named list if \code{return_all = TRUE}.
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#'
#' toy <- tibble::tibble(
#'   SampleID = c("S1","S1","S2","S2"),
#'   Group    = c("G1","G1","G1","G1"),
#'   Class    = c("A","B","A","C"),
#'   RA       = c(0.7, 0.3, 0.6, 0.4)
#' )
#'
#' res <- create_alluvial_plot(
#'   data = toy,
#'   tax_level = "Class",
#'   group_col = "Group",
#'   prepare_args = list(complete_zero = TRUE, clean_taxonomy = FALSE),
#'   classify_args = list(detect_mixed_abundance = TRUE),
#'   plot_args = list(alpha = 0.7),
#'   return_all = TRUE
#' )
#'
#' res$plot
#'}
#' @export
create_alluvial_plot <- function(data,
                                 tax_level,
                                 group_col,
                                 SampleID_col = "SampleID",
                                 groups = NULL,
                                 low_abundance_threshold = 0.01,
                                 palette_list = NULL,
                                 custom_palette = NULL,
                                 prepare_args = list(),
                                 classify_args = list(),
                                 palette_args = list(),
                                 plot_args = list(),
                                 return_all = FALSE) {

  # ---- 1) Prepare ----
  prep_call <- c(
    list(
      data = data,
      tax_level = tax_level,
      group_col = group_col,
      SampleID_col = SampleID_col,
      groups = groups
    ),
    prepare_args
  )
  allu_data <- do.call(prepare_alluvial_data, prep_call)

  # ---- 2) Classify ----
  class_call <- c(
    list(
      data = allu_data,
      tax_level = tax_level,
      group_col = group_col,
      groups = groups,
      low_abundance_threshold = low_abundance_threshold
    ),
    classify_args
  )
  allu_classified <- do.call(classify_taxa_patterns, class_call)

  # ---- 3) Palette ----
  if (is.null(custom_palette)) {
    pal_call <- c(
      list(
        data = allu_classified,
        tax_level = tax_level,
        palette_list = palette_list
      ),
      palette_args
    )
    allu_palette <- do.call(generate_alluvial_palette, pal_call)
  } else {
    allu_palette <- custom_palette
  }

  # ---- 4) Plot ----
  plot_call <- c(
    list(
      data = allu_classified,
      custom_palette = allu_palette,
      tax_level = tax_level,
      group_col = group_col
    ),
    plot_args
  )
  p <- do.call(plot_alluvial, plot_call)

  if (isTRUE(return_all)) {
    return(list(
      data_prepared = allu_data,
      data_classified = allu_classified,
      palette = allu_palette,
      plot = p
    ))
  }

  p
}
