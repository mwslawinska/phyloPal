#' Create grouped ASV matrix by averaging within groups
#'
#' Collapses a sample-by-taxon matrix into a group-by-taxon matrix by taking the mean
#' abundance across samples within each group. This is useful for building dendrograms
#' on group centroids (e.g. habitat-level mean communities).
#'
#' @param asv_matrix Numeric matrix with samples as columns and taxa/ASVs as rows.
#'   Column names must match \code{metadata[[sample_col]]}.
#' @param metadata Data frame with one row per sample (or at least one row per sample ID).
#' @param sample_col Column name in \code{metadata} containing sample IDs (default: \code{"SampleID"}).
#' @param group_col Column name in \code{metadata} defining groups to average within.
#' @param group_order How to order group columns in the output:
#'   \itemize{
#'     \item \code{"alphabetical"}: sort group names A--Z (default)
#'     \item \code{"metadata"}: order by first appearance in \code{metadata[[group_col]]}
#'     \item \code{"custom"}: use \code{group_levels} (and append unseen groups at the end)
#'   }
#' @param group_levels Character vector giving desired group order when \code{group_order="custom"}.
#' @param drop_unmapped_samples Logical; if \code{TRUE}, drop samples missing a group assignment.
#'   If \code{FALSE} (default), error if any sample in \code{asv_matrix} is missing \code{group_col} in metadata.
#'
#' @return Numeric matrix with groups as columns and taxa/ASVs as rows (mean abundance per group).
#'
#' @examples
#' mat <- matrix(c(1, 2, 3, 4), nrow = 2)
#' rownames(mat) <- c("ASV1", "ASV2")
#' colnames(mat) <- c("S1", "S2")
#' meta <- data.frame(SampleID = c("S1", "S2"), Group = c("G1", "G1"))
#' create_grouped_matrix(mat, meta, sample_col = "SampleID", group_col = "Group")
#'
#' @export
create_grouped_matrix <- function(asv_matrix,
                                  metadata,
                                  sample_col = "SampleID",
                                  group_col,
                                  group_order = c("alphabetical", "metadata", "custom"),
                                  group_levels = NULL,
                                  drop_unmapped_samples = FALSE) {
  group_order <- match.arg(group_order)

  if (!sample_col %in% colnames(metadata)) {
    stop("sample_col '", sample_col, "' not found in metadata")
  }
  if (!group_col %in% colnames(metadata)) {
    stop("group_col '", group_col, "' not found in metadata")
  }

  asv_df <- as.data.frame(asv_matrix, check.names = FALSE)
  asv_df$ASV_ID <- rownames(asv_df)

  asv_long <- asv_df %>%
    tidyr::pivot_longer(
      cols = -ASV_ID,
      names_to = sample_col,
      values_to = "Abundance"
    ) %>%
    dplyr::left_join(metadata, by = sample_col)

  # samples missing group assignment
  missing_group_samples <- asv_long %>%
    dplyr::filter(is.na(.data[[group_col]])) %>%
    dplyr::pull(.data[[sample_col]]) %>%
    unique()

  if (length(missing_group_samples) > 0) {
    if (isTRUE(drop_unmapped_samples)) {
      asv_long <- asv_long %>% dplyr::filter(!is.na(.data[[group_col]]))
    } else {
      stop(
        "Some samples in asv_matrix are missing '", group_col, "' in metadata: ",
        paste(missing_group_samples, collapse = ", ")
      )
    }
  }

  asv_grouped <- asv_long %>%
    dplyr::group_by(ASV_ID, .data[[group_col]]) %>%
    dplyr::summarise(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop")

  # Determine group levels (column order)
  groups_present <- as.character(unique(asv_grouped[[group_col]]))

  if (group_order == "alphabetical") {
    group_levels_final <- sort(groups_present)
  } else if (group_order == "metadata") {
    # First-appearance order in metadata (stable and intuitive)
    group_levels_final <- metadata %>%
      dplyr::pull(.data[[group_col]]) %>%
      as.character() %>%
      (\(x) x[!is.na(x)])() %>%
      unique()
    group_levels_final <- intersect(group_levels_final, groups_present)
  } else { # custom
    if (is.null(group_levels) || length(group_levels) == 0) {
      stop("group_levels must be provided when group_order='custom'.")
    }
    group_levels_final <- intersect(group_levels, groups_present)
    missing_levels <- setdiff(groups_present, group_levels_final)
    if (length(missing_levels) > 0) {
      warning(
        "Some groups are not in group_levels and will be appended at the end: ",
        paste(missing_levels, collapse = ", ")
      )
      group_levels_final <- c(group_levels_final, sort(missing_levels))
    }
  }

  asv_wide <- asv_grouped %>%
    dplyr::mutate(!!group_col := factor(.data[[group_col]], levels = group_levels_final)) %>%
    tidyr::pivot_wider(
      names_from = !!rlang::sym(group_col),
      values_from = Abundance,
      values_fill = 0
    )

  asv_matrix_grouped <- as.matrix(asv_wide[, -1, drop = FALSE])
  rownames(asv_matrix_grouped) <- asv_wide$ASV_ID

  asv_matrix_grouped
}

#' Build a dendrogram (and leaf order) from an ASV matrix
#'
#' Computes a dissimilarity matrix between columns (samples or groups) using
#' \code{vegan::vegdist()} on \code{t(mat)}, clusters with \code{stats::hclust()},
#' and returns a dendrogram with an optional \code{hang}.
#'
#' @param mat Numeric matrix with taxa/ASVs as rows and units-to-cluster as columns
#'   (samples or groups). Columns will be clustered.
#' @param distance_method Distance method passed to \code{vegan::vegdist()} (default: \code{"bray"}).
#' @param cluster_method Linkage method passed to \code{stats::hclust()}
#'   (default: \code{"ward.D2"}).
#' @param hang Numeric hang parameter passed to \code{dendextend::hang.dendrogram()}
#'   (default: -1).
#'
#' @return A list with:
#' \describe{
#'   \item{dendrogram}{A \code{dendrogram} object.}
#'   \item{order}{Character vector of leaf labels in plotting order.}
#'   \item{hc}{The \code{hclust} object.}
#'   \item{dist}{The \code{dist} object from \code{vegan::vegdist}.}
#' }
#'
#' @examples
#' if (requireNamespace("vegan", quietly = TRUE) &&
#'   requireNamespace("dendextend", quietly = TRUE)) {
#'   mat <- matrix(runif(20), nrow = 5)
#'   colnames(mat) <- c("G1", "G2", "G3", "G4")
#'   res <- build_dendrogram(mat)
#'   res$order
#' }
#'
#' @export
build_dendrogram <- function(mat,
                             distance_method = "bray",
                             cluster_method = "ward.D2",
                             hang = -1) {
  if (!requireNamespace("vegan", quietly = TRUE)) {
    stop("Package 'vegan' required. Install with: install.packages('vegan')")
  }
  if (!requireNamespace("dendextend", quietly = TRUE)) {
    stop("Package 'dendextend' required. Install with: install.packages('dendextend')")
  }

  d <- vegan::vegdist(t(mat), method = distance_method)
  hc <- stats::hclust(d, method = cluster_method)
  dend <- stats::as.dendrogram(hc)
  dend <- dendextend::hang.dendrogram(dend, hang = hang)

  list(
    dendrogram = dend,
    order = labels(dend),
    hc = hc,
    dist = d
  )
}

#' Plot a dendrogram as ggplot with optional tip aesthetics (color/fill/shape/size) and labels
#'
#' Converts a \code{dendrogram} into ggplot layers using \code{ggdendro}. Draws branches
#' as line segments and overlays tip (leaf) points. Tip aesthetics can be mapped to
#' outline color (\code{color}), interior fill (\code{fill}; visible for shapes 21-25),
#' shape (\code{shape}), and size (\code{size}).
#'
#' Leaf order (left-to-right in the rendered dendrogram layout *before orientation transforms*)
#' is stored on the returned plot as \code{attr(plot, "leaf_order")}. This is useful to order
#' the alluvial axis to match the dendrogram.
#'
#' @param dend A \code{dendrogram} object, or a list from \code{\link{build_dendrogram}} containing \code{$dendrogram}.
#' @param metadata Optional data frame used to build \code{tip_meta} automatically (recommended).
#' @param label_from Column in \code{metadata} that matches the dendrogram leaf labels.
#'
#' @param color_by Optional column mapped to point outline color.
#' @param color_palette Optional named vector for \code{color_by} levels.
#' @param fill_by Optional column mapped to point interior fill.
#' @param fill_palette Optional named vector for \code{fill_by} levels.
#' @param shape_by Optional column mapped to point shape.
#' @param shape_palette Optional values passed to \code{ggplot2::scale_shape_manual()}.
#' @param size_by Optional column mapped to point size.
#' @param size_palette Optional named numeric vector passed to \code{ggplot2::scale_size_manual()}.
#' @param tip_meta Optional precomputed tip metadata. Must contain \code{label_col} and any requested mapping columns.
#' @param label_col Name of the label column inside \code{tip_meta} (default: \code{"label"}).
#' @param point_size Numeric tip point size (default: 2). Ignored when \code{size_by} is used.
#' @param point_stroke Numeric outline width for points (\code{stroke}); only applies to shapes 21--25.
#'   If \code{NULL}, uses ggplot2 default.
#' @param default_shape Integer shape used when \code{fill_by} is set but \code{shape_by} is \code{NULL}
#'   (default: 21). This ensures mapped fills are visible.
#' @param orientation Dendrogram orientation:
#'   \describe{
#'     \item{\code{"top"}}{Tips at bottom (default ggdendro layout).}
#'     \item{\code{"bottom"}}{Tips at top (vertical mirror).}
#'     \item{\code{"right"}}{Tips on left, branches extend to the right (rotated 90°).}
#'     \item{\code{"left"}}{Tips on right, branches extend to the left (rotated 90° + mirrored).}
#'   }
#' @param first_branch For \code{orientation="left"} or \code{"right"}, controls whether the first leaf
#'   in the layout appears at the top (\code{"top"}) or bottom (\code{"bottom"}). This is implemented by
#'   reversing the (post-flip) vertical ordering via \code{scale_x_reverse()}.
#'   For \code{orientation="top"} or \code{"bottom"}, this parameter is currently ignored.
#' @param add_tip_labels Logical; if \code{TRUE}, draw tip labels.
#' @param tip_label_size Numeric text size for tip labels.
#' @param leaf_order Optional character vector giving the desired leaf order (matching dendrogram labels).
#' @param rotate_to_leaf_order Logical; if \code{TRUE} (default), rotate dendrogram to match \code{leaf_order}
#'   using \code{dendextend::rotate()} (requires dendextend).
#' @param theme_obj ggplot2 theme to apply (default: \code{ggplot2::theme_void()}).
#' @param guides_obj Optional \code{ggplot2::guides(...)} added to the plot.
#' @param legend_position Legend position passed to \code{theme(legend.position=...)}.
#' @param legend_position_inside Numeric length-2 vector (x, y) between 0 and 1 used when \code{legend_position="inside"}.
#'
#' @return A \code{ggplot} object. Leaf order (left-to-right) is stored as \code{attr(plot, "leaf_order")}.
#' @export
plot_dendrogram <- function(dend,
                            metadata = NULL,
                            label_from = NULL,
                            color_by = NULL,
                            color_palette = NULL,
                            fill_by = NULL,
                            fill_palette = NULL,
                            shape_by = NULL,
                            shape_palette = NULL,
                            size_by = NULL,
                            size_palette = NULL,
                            tip_meta = NULL,
                            label_col = "label",
                            point_size = 2,
                            point_stroke = NULL, # NULL = ggplot default
                            default_shape = 21,
                            orientation = c("top", "right", "bottom", "left"),
                            first_branch = c("bottom", "top"), # relevant for left/right (post-flip vertical order)
                            add_tip_labels = FALSE,
                            tip_label_size = 3,
                            leaf_order = NULL,
                            rotate_to_leaf_order = TRUE,
                            theme_obj = ggplot2::theme_void(),
                            guides_obj = NULL,
                            legend_position = "right",
                            legend_position_inside = NULL) {

  orientation <- match.arg(orientation)
  first_branch <- match.arg(first_branch)

  if (!requireNamespace("ggdendro", quietly = TRUE)) {
    stop("Package 'ggdendro' required. Install with: install.packages('ggdendro')")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' required.")
  if (!requireNamespace("rlang", quietly = TRUE)) stop("Package 'rlang' required.")

  # Accept list output from build_dendrogram()
  dend_obj <- dend
  if (is.list(dend) && "dendrogram" %in% names(dend)) dend_obj <- dend$dendrogram

  # ---- Optional: rotate dendrogram to a desired tip order (freeze layout) ----
  if (!is.null(leaf_order) && isTRUE(rotate_to_leaf_order)) {
    if (!requireNamespace("dendextend", quietly = TRUE)) {
      stop("Package 'dendextend' required for rotate_to_leaf_order=TRUE. Install with: install.packages('dendextend')")
    }
    leaf_order2 <- intersect(as.character(leaf_order), labels(dend_obj))
    if (length(leaf_order2) == 0) stop("leaf_order does not match any labels in the dendrogram.")
    dend_obj <- dendextend::rotate(dend_obj, order = leaf_order2)
  }

  dd <- ggdendro::dendro_data(dend_obj)
  seg <- ggdendro::segment(dd)
  lab <- ggdendro::label(dd) # x, y, label
  lab <- lab %>% dplyr::rename(!!label_col := label)

  # Leaf order (left-to-right in base ggdendro layout, before orientation transforms)
  leaf_order_from_plot <- lab %>%
    dplyr::arrange(x) %>%
    dplyr::pull(.data[[label_col]]) %>%
    as.character()

  # ---- Build tip_meta automatically if not provided ----
  if (is.null(tip_meta)) {
    if (is.null(metadata)) stop("Provide either tip_meta or metadata.")
    if (is.null(label_from)) stop("Provide label_from (column in metadata matching dendrogram labels).")
    if (!label_from %in% names(metadata)) stop("label_from not found in metadata: ", label_from)

    needed <- c(label_from, color_by, fill_by, shape_by, size_by)
    needed <- needed[!is.null(needed)]
    needed <- unique(needed)

    missing_needed <- setdiff(needed, names(metadata))
    if (length(missing_needed) > 0) stop("Missing columns in metadata: ", paste(missing_needed, collapse = ", "))

    tip_meta <- metadata %>%
      dplyr::select(dplyr::all_of(needed)) %>%
      dplyr::mutate(.lab = as.character(.data[[label_from]])) %>%
      dplyr::distinct(.lab, .keep_all = TRUE) %>%
      dplyr::rename(!!label_col := .lab)
  } else {
    if (!label_col %in% names(tip_meta)) stop("tip_meta must contain '", label_col, "'.")
    tip_meta[[label_col]] <- as.character(tip_meta[[label_col]])
  }

  lab2 <- lab %>% dplyr::left_join(tip_meta, by = label_col)

  # ---- Palette coverage checks ----
  if (!is.null(color_by)) {
    present <- unique(stats::na.omit(as.character(lab2[[color_by]])))
    if (is.null(color_palette)) stop("color_palette must be provided when color_by is used.")
    missing_cols <- setdiff(present, names(color_palette))
    if (length(missing_cols) > 0) stop("color_palette missing: ", paste(missing_cols, collapse = ", "))
  }
  if (!is.null(fill_by)) {
    present <- unique(stats::na.omit(as.character(lab2[[fill_by]])))
    if (is.null(fill_palette)) stop("fill_palette must be provided when fill_by is used.")
    missing_fills <- setdiff(present, names(fill_palette))
    if (length(missing_fills) > 0) stop("fill_palette missing: ", paste(missing_fills, collapse = ", "))
  }
  if (!is.null(size_by)) {
    present <- unique(stats::na.omit(as.character(lab2[[size_by]])))
    if (is.null(size_palette)) stop("size_palette must be provided when size_by is used.")
    missing_sizes <- setdiff(present, names(size_palette))
    if (length(missing_sizes) > 0) stop("size_palette missing: ", paste(missing_sizes, collapse = ", "))
  }

  # ---- Base plot ----
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = seg,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      linewidth = 0.3
    )

  # Tip point aesthetics (outline=color, interior=fill like ggplot)
  aes_args <- list(x = quote(x), y = quote(y))
  if (!is.null(color_by)) aes_args$color <- rlang::sym(color_by)
  if (!is.null(fill_by))  aes_args$fill  <- rlang::sym(fill_by)
  if (!is.null(shape_by)) aes_args$shape <- rlang::sym(shape_by)
  if (!is.null(size_by))  aes_args$size  <- rlang::sym(size_by)

  point_args <- list(
    data = lab2,
    mapping = do.call(ggplot2::aes, aes_args)
  )

  # fixed aesthetics (only when not mapped)
  if (is.null(size_by)) point_args$size <- point_size
  if (!is.null(point_stroke)) point_args$stroke <- point_stroke

  # IMPORTANT: if fill is mapped but shape is not, force a fillable shape
  if (!is.null(fill_by) && is.null(shape_by)) {
    point_args$shape <- default_shape
  }

  p <- p + do.call(ggplot2::geom_point, point_args)

  # Optional tip labels
  if (isTRUE(add_tip_labels)) {
    p <- p + ggplot2::geom_text(
      data = lab2,
      ggplot2::aes(x = x, y = y, label = .data[[label_col]]),
      size = tip_label_size,
      hjust = -0.1,
      vjust = 0.5
    )
  }

  # ---- Scales ----
  if (!is.null(color_by)) {
    p <- p + ggplot2::scale_color_manual(values = color_palette, breaks = names(color_palette))
  }
  if (!is.null(fill_by)) {
    p <- p + ggplot2::scale_fill_manual(values = fill_palette, breaks = names(fill_palette))
  }
  if (!is.null(shape_by) && !is.null(shape_palette)) {
    p <- p + ggplot2::scale_shape_manual(values = shape_palette, breaks = names(shape_palette))
  }
  if (!is.null(size_by)) {
    p <- p + ggplot2::scale_size_manual(values = size_palette, breaks = names(size_palette))
  }

  # ---- Orientation ----
  if (orientation == "right") {
    p <- p + ggplot2::coord_flip(clip = "off")
    if (first_branch == "top") p <- p + ggplot2::scale_x_reverse()
  } else if (orientation == "left") {
    p <- p + ggplot2::coord_flip(clip = "off") + ggplot2::scale_y_reverse()
    if (first_branch == "top") p <- p + ggplot2::scale_x_reverse()
  } else if (orientation == "bottom") {
    p <- p + ggplot2::coord_cartesian(clip = "off") + ggplot2::scale_y_reverse()
  } else {
    p <- p + ggplot2::coord_cartesian(clip = "off")
  }

  # ---- Theme ----
  p <- p +
    theme_obj +
    ggplot2::theme(
      legend.position = legend_position,
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank()
    )

  # Optional inside legend position (used by some ggplot2 extensions)
  if (!is.null(legend_position_inside)) {
    p <- p + ggplot2::theme(legend.position.inside = legend_position_inside)
  }

  if (!is.null(guides_obj)) p <- p + guides_obj

  attr(p, "leaf_order") <- leaf_order_from_plot
  p
}


#' Extract a legend grob from a ggplot
#'
#' Convenience helper around \code{cowplot::get_legend()}.
#'
#' @param p A \code{ggplot} object.
#'
#' @return A grob containing the legend, or \code{NULL} if the plot has no legend.
#'
#' @export
extract_legend_grob <- function(p) {
  if (!requireNamespace("cowplot", quietly = TRUE)) {
    stop("Package 'cowplot' required. Install with install.packages('cowplot').")
  }
  cowplot::get_legend(p)
}


#' Combine dendrogram and alluvial plot with alignment and flexible legend placement
#'
#' Compose a dendrogram panel and an alluvial panel into one figure with stable alignment.
#' Supports vertical stacking (top/bottom) and horizontal composition (left/right).
#'
#' @param alluvial_plot A plain ggplot (alluvial panel). Must not be an aplot composite.
#' @param dendrogram_plot A plain ggplot (dendrogram panel). Must not be an aplot composite.
#'
#' @param dend_position "top","bottom","left","right". Default "top".
#' @param dend_height Relative size of dendrogram panel. For top/bottom = fraction of height;
#'   for left/right = fraction of width. Default 0.15.
#' @param strip_alluvial_x If TRUE, remove alluvial x axis title/text/ticks (useful for top/bottom stacks).
#' @param legend "separate","omit","together"
#' @param legend_source "alluvial","dendrogram","both"
#' @param legend_position "right","bottom" (only used when legend="separate")
#' @param legend_rel_width Relative legend column width (legend_position="right")
#' @param legend_rel_height Relative legend row height (legend_position="bottom")
#' @param alluvial_margins Plot margin for alluvial panel (ggplot2::margin()).
#' @param dendrogram_margins Plot margin for dendrogram panel (ggplot2::margin()).
#' @param outer_margins Outer margin around full grob. Accepts ggplot2::margin(), unit length,
#'   or numeric c(t,r,b,l) in pt.
#' @param alluvial_legend_overrides Optional list of ggplot additions applied only during legend extraction.
#' @param dendrogram_legend_overrides Optional list of ggplot additions applied only during legend extraction.
#' @param legend_box "vertical" or "horizontal" when legend_source="both".
#' @param align "panel" (recommended) or "full_grob".
#'   "panel" equalizes only the panel region; "full_grob" uses unit.pmax on all slots.
#' @param x_expand_zero If TRUE, set expand=0 on x scales where no x scale is defined yet.
#' @param leaf_order Character vector giving leaf order. If NULL, uses attr(dendrogram_plot,"leaf_order").
#' @param align_x_centers For top/bottom: lock x-centers between dend tips and alluvial strata.
#'   Uses scale_x_discrete(limits=leaf_order) on alluvial, and coord_cartesian(xlim=...) on dendrogram.
#' @param overwrite_x_scales If TRUE, overwrite existing x scales when applying align_x_centers.
#' @param x_labels Optional relabeling for the alluvial x axis when
#'   \code{align_x_centers = TRUE}. Does not change the underlying factor
#'   levels used for alignment. Accepts:
#'   \itemize{
#'     \item A named character vector mapping \code{leaf_order} values
#'       (e.g. Bio_SampleID) to display labels.
#'     \item A function passed to \code{scale_x_discrete(labels = ...)}.
#'   }
#'   If \code{NULL} (default), original x scale labels are used.
#'   When \code{overwrite_x_scales = TRUE}, any existing x scale on the
#'   alluvial plot will be replaced and \code{x_labels} will be applied.
#' @param x_drop Logical; passed to \code{scale_x_discrete(drop = ...)}
#'   when enforcing alignment with \code{leaf_order}. Default \code{FALSE}.
#'   Set to \code{TRUE} to drop unused factor levels.
#' @param align_y_centers For left/right: lock y-centers. Intended for dendrogram plotted with
#'   plot_dendrogram(orientation="left"/"right") and alluvial with coord_flip().
#'   IMPORTANT: even with coord_flip(), group order is still controlled by scale_x_discrete on alluvial.
#'   Uses coord_cartesian(ylim=...) on dendrogram to avoid dropping segments.
#' @param overwrite_y_scales If TRUE, overwrite existing scales used for align_y_centers.
#' @param dend_limits_left,dend_limits_right “sweet spot” for dendrogram x view:
#'   dend xlim = c(dend_limits_left, n + dend_limits_right). Use e.g. -0.75 and 0.2.
#' @param dend_limits_bottom,dend_limits_top “sweet spot” for dendrogram y view:
#'   dend ylim = c(dend_limits_bottom, n + dend_limits_top).
#'
#' #' @details
#' When \code{align_x_centers = TRUE}, the function enforces identical
#' x positions between dendrogram tips and alluvial strata using
#' \code{leaf_order}. The underlying factor levels are not modified.
#' Use \code{x_labels} to change only the displayed axis text
#' (e.g. mapping Bio_SampleID to Subsite_BID) without breaking alignment.
#' To place the dendrogram ideally above the alluvial plot (leaves in the centre of bars) 
#' try different combinations of \code{dend_limits_left},\code{dend_limits_right} 
#' or \code{dend_limits_bottom},\code{dend_limits_top}.
#' 
#' 
#' @return A grob (gridExtra::arrangeGrob) suitable for ggsave().
#' @export
combine_dendrogram_alluvial <- function(alluvial_plot,
                                        dendrogram_plot,
                                        dend_position = c("top", "bottom", "left", "right"),
                                        dend_height = 0.15,
                                        strip_alluvial_x = TRUE,
                                        legend = c("separate", "omit", "together"),
                                        legend_source = c("alluvial", "dendrogram", "both"),
                                        legend_position = c("right", "bottom"),
                                        legend_rel_width = 0.28,
                                        legend_rel_height = 0.18,
                                        alluvial_margins   = ggplot2::margin(0, 0, 0, 0),
                                        dendrogram_margins = ggplot2::margin(0, 0, 0, 0),
                                        outer_margins      = ggplot2::margin(0, 0, 0, 0),
                                        alluvial_legend_overrides = NULL,
                                        dendrogram_legend_overrides = NULL,
                                        legend_box = c("vertical", "horizontal"),
                                        align = c("panel", "full_grob"),
                                        x_expand_zero = FALSE,
                                        leaf_order = NULL,
                                        align_x_centers = FALSE,
                                        overwrite_x_scales = FALSE,
                                        x_labels = NULL,
                                        x_drop = FALSE, 
                                        dend_limits_left = 0.5,
                                        dend_limits_right = 0.5,
                                        align_y_centers = FALSE,
                                        overwrite_y_scales = FALSE,
                                        dend_limits_bottom = 0.5,
                                        dend_limits_top = 0.5) {
  dend_position <- match.arg(dend_position)
  legend <- match.arg(legend)
  legend_source <- match.arg(legend_source)
  legend_position <- match.arg(legend_position)
  legend_box <- match.arg(legend_box)
  align <- match.arg(align)

  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required.")
  if (!requireNamespace("cowplot", quietly = TRUE)) stop("cowplot required.")
  if (!requireNamespace("grid", quietly = TRUE)) stop("grid required.")
  if (!requireNamespace("gridExtra", quietly = TRUE)) stop("gridExtra required.")

  # ---- Guard: aplot objects will break ggplot_build() ----
  if (inherits(alluvial_plot, "aplot")) {
    stop("alluvial_plot is an 'aplot' object. Pass the original ggplot (before insert_*).")
  }
  if (inherits(dendrogram_plot, "aplot")) {
    stop("dendrogram_plot is an 'aplot' object. Pass the original ggplot (before insert_*).")
  }
  if (!inherits(alluvial_plot, "ggplot")) stop("alluvial_plot must be a ggplot.")
  if (!inherits(dendrogram_plot, "ggplot")) stop("dendrogram_plot must be a ggplot.")

  # ---- Helpers ----
  add_overrides <- function(p, overrides) {
    if (is.null(overrides)) return(p)
    for (obj in overrides) p <- p + obj
    p
  }

  has_scale <- function(p, axis = c("x", "y")) {
    axis <- match.arg(axis)
    any(vapply(p$scales$scales, function(s) axis %in% s$aesthetics, logical(1)))
  }

  add_outer_margins <- function(g, m, default_unit = "pt") {
    as_unit1 <- function(x) {
      if (inherits(x, "unit")) return(x)
      if (is.numeric(x) && length(x) == 1) return(grid::unit(x, default_unit))
      stop("Unsupported margin component type: ", paste(class(x), collapse = "/"))
    }

    if (inherits(m, "margin")) {
      top    <- as_unit1(m[["t"]])
      right  <- as_unit1(m[["r"]])
      bottom <- as_unit1(m[["b"]])
      left   <- as_unit1(m[["l"]])
    } else if (inherits(m, "unit") && length(m) == 4) {
      top <- m[[1]]; right <- m[[2]]; bottom <- m[[3]]; left <- m[[4]]
    } else if (is.numeric(m) && length(m) == 4) {
      top    <- grid::unit(m[[1]], default_unit)
      right  <- grid::unit(m[[2]], default_unit)
      bottom <- grid::unit(m[[3]], default_unit)
      left   <- grid::unit(m[[4]], default_unit)
    } else {
      stop("outer_margins must be ggplot2::margin(), grid::unit length-4, or numeric c(t,r,b,l).")
    }

    gridExtra::arrangeGrob(
      grobs = list(
        grid::nullGrob(), grid::nullGrob(), grid::nullGrob(),
        grid::nullGrob(), g,                 grid::nullGrob(),
        grid::nullGrob(), grid::nullGrob(), grid::nullGrob()
      ),
      ncol = 3,
      widths  = grid::unit.c(left, grid::unit(1, "null"), right),
      heights = grid::unit.c(top,  grid::unit(1, "null"), bottom)
    )
  }

  panel_cols <- function(g) {
    idx <- which(grepl("^panel", g$layout$name))
    sort(unique(unlist(Map(seq, g$layout$l[idx], g$layout$r[idx]))))
  }
  panel_rows <- function(g) {
    idx <- which(grepl("^panel", g$layout$name))
    sort(unique(unlist(Map(seq, g$layout$t[idx], g$layout$b[idx]))))
  }

  align_panel_widths <- function(g1, g2) {
    c1 <- panel_cols(g1)
    c2 <- panel_cols(g2)
    if (length(c1) == 0 || length(c2) == 0) {
      w <- grid::unit.pmax(g1$widths, g2$widths)
      g1$widths <- w; g2$widths <- w
      return(list(g1 = g1, g2 = g2))
    }
    w_panel <- grid::unit.pmax(sum(g1$widths[c1]), sum(g2$widths[c2]))
    g1$widths[c1] <- grid::unit(rep(0, length(c1)), "pt")
    g2$widths[c2] <- grid::unit(rep(0, length(c2)), "pt")
    g1$widths[c1[1]] <- w_panel
    g2$widths[c2[1]] <- w_panel
    list(g1 = g1, g2 = g2)
  }

  align_panel_heights <- function(g1, g2) {
    r1 <- panel_rows(g1)
    r2 <- panel_rows(g2)
    if (length(r1) == 0 || length(r2) == 0) {
      h <- grid::unit.pmax(g1$heights, g2$heights)
      g1$heights <- h; g2$heights <- h
      return(list(g1 = g1, g2 = g2))
    }
    h_panel <- grid::unit.pmax(sum(g1$heights[r1]), sum(g2$heights[r2]))
    g1$heights[r1] <- grid::unit(rep(0, length(r1)), "pt")
    g2$heights[r2] <- grid::unit(rep(0, length(r2)), "pt")
    g1$heights[r1[1]] <- h_panel
    g2$heights[r2[1]] <- h_panel
    list(g1 = g1, g2 = g2)
  }
  normalize_x_labels <- function(x_labels, leaf_order) {
  if (is.null(x_labels)) return(NULL)

  # function labels are fine
  if (is.function(x_labels)) return(x_labels)

  # named vector mapping is best
  if (is.atomic(x_labels) && !is.null(names(x_labels))) {
    # ensure we can index in leaf order
    lab <- x_labels[leaf_order]
    # if any NA, fall back to original leaf values for those
    lab[is.na(lab)] <- leaf_order[is.na(lab)]
    return(lab)
  }

  stop("x_labels must be NULL, a function, or a named vector with names matching leaf_order.")
}

  # ---- Base plots ----
  p_allu <- alluvial_plot   + ggplot2::theme(plot.margin = alluvial_margins)
  p_dend <- dendrogram_plot + ggplot2::theme(plot.margin = dendrogram_margins)

  if (isTRUE(strip_alluvial_x)) {
    p_allu <- p_allu + ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.text.x  = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )
  }

  # Optional: reduce drift from expansion (only if axis has no scale yet)
  if (isTRUE(x_expand_zero)) {
    if (!has_scale(p_allu, "x")) {
      p_allu <- p_allu + ggplot2::scale_x_discrete(expand = ggplot2::expansion(mult = 0, add = 0))
    }
    if (!has_scale(p_dend, "x")) {
      p_dend <- p_dend + ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0, add = 0))
    }
    if (!has_scale(p_dend, "y")) {
      p_dend <- p_dend + ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0, add = 0))
    }
  }

  # Leaf order for center locking
  if (is.null(leaf_order)) leaf_order <- attr(dendrogram_plot, "leaf_order")
  if ((isTRUE(align_x_centers) || isTRUE(align_y_centers)) && is.null(leaf_order)) {
    stop("align_*_centers=TRUE requires leaf_order (or dendrogram_plot must have attr(,'leaf_order')).")
  }
  if (!is.null(leaf_order)) {
    leaf_order <- as.character(leaf_order)
    n <- length(leaf_order)
  }

  # ---- Center lock for TOP/BOTTOM (tips vary along x) ----
  if (isTRUE(align_x_centers)) {
    can_set_allu <- overwrite_x_scales || !has_scale(p_allu, "x")
    can_set_dend <- overwrite_x_scales || !has_scale(p_dend, "x")

    if (can_set_allu) {
  x_lab <- normalize_x_labels(x_labels, leaf_order)

  labels_arg <- if (is.null(x_lab)) ggplot2::waiver() else x_lab

p_allu <- p_allu +
  ggplot2::scale_x_discrete(
    limits = leaf_order,
    labels = labels_arg,
    drop   = x_drop,
    expand = ggplot2::expansion(mult = 0, add = 0)
  )
}

    if (can_set_dend) {
      # ZOOM (no dropping) -> avoids "Removed rows..." warnings
      p_dend <- p_dend +
        ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0, add = 0)) +
        ggplot2::coord_cartesian(
          xlim = c(dend_limits_left, n + dend_limits_right),
          clip = "off"
        )
    }
  }

  # ---- Center lock for LEFT/RIGHT (tips vary along y) ----
  if (isTRUE(align_y_centers)) {
    # IMPORTANT:
    # even if alluvial uses coord_flip(), its group order is still controlled by scale_x_discrete()
    can_set_allu <- overwrite_y_scales || !has_scale(p_allu, "x")
    can_set_dend <- overwrite_y_scales || !has_scale(p_dend, "y")

    if (can_set_allu) {
      p_allu <- p_allu +
        ggplot2::scale_x_discrete(
          limits = leaf_order,
          expand = ggplot2::expansion(mult = 0, add = 0)
        )
    }

    if (can_set_dend) {
      p_dend <- p_dend +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0, add = 0)) +
        ggplot2::coord_cartesian(
          ylim = c(dend_limits_bottom, n + dend_limits_top),
          clip = "off"
        )
    }
  }

  # ---- Legend extraction ----
  legend_grobs <- list()

  if (legend == "separate") {
    if (legend_source %in% c("alluvial", "both")) {
      p_leg <- add_overrides(p_allu, alluvial_legend_overrides)
      lg <- cowplot::get_legend(p_leg)
      if (!is.null(lg)) legend_grobs$alluvial <- lg
      p_allu <- p_allu + ggplot2::theme(legend.position = "none")
    }
    if (legend_source %in% c("dendrogram", "both")) {
      p_leg <- add_overrides(p_dend, dendrogram_legend_overrides)
      lg <- cowplot::get_legend(p_leg)
      if (!is.null(lg)) legend_grobs$dendrogram <- lg
      p_dend <- p_dend + ggplot2::theme(legend.position = "none")
    }
  } else if (legend == "omit") {
    p_allu <- p_allu + ggplot2::theme(legend.position = "none")
    p_dend <- p_dend + ggplot2::theme(legend.position = "none")
  }

  combined_legend <- NULL
  if (legend == "separate" && length(legend_grobs) > 0) {
    lg_list <- legend_grobs[c("dendrogram", "alluvial")]
    lg_list <- lg_list[!vapply(lg_list, is.null, logical(1))]

    combined_legend <- if (length(lg_list) == 1) {
      lg_list[[1]]
    } else if (legend_box == "vertical") {
      gridExtra::arrangeGrob(grobs = lg_list, ncol = 1)
    } else {
      gridExtra::arrangeGrob(grobs = lg_list, nrow = 1)
    }
  }

  # ---- Convert to grobs + align ----
  g_allu <- ggplot2::ggplotGrob(p_allu)
  g_dend <- ggplot2::ggplotGrob(p_dend)

  # ---- Compose ----
  out <- NULL

  if (dend_position %in% c("top", "bottom")) {
    if (align == "panel") {
      res <- align_panel_widths(g_allu, g_dend)
      g_allu <- res$g1; g_dend <- res$g2
    } else {
      w <- grid::unit.pmax(g_allu$widths, g_dend$widths)
      g_allu$widths <- w; g_dend$widths <- w
    }

    rel_heights <- c(dend_height, 1 - dend_height)
    stacked <- if (dend_position == "top") {
      gridExtra::arrangeGrob(g_dend, g_allu, ncol = 1, heights = rel_heights)
    } else {
      gridExtra::arrangeGrob(g_allu, g_dend, ncol = 1, heights = rel_heights)
    }
    out <- stacked

  } else { # left/right
    if (align == "panel") {
      res <- align_panel_heights(g_allu, g_dend)
      g_allu <- res$g1; g_dend <- res$g2
    } else {
      h <- grid::unit.pmax(g_allu$heights, g_dend$heights)
      g_allu$heights <- h; g_dend$heights <- h
    }

    rel_widths <- c(dend_height, 1 - dend_height)
    stacked <- if (dend_position == "left") {
      gridExtra::arrangeGrob(g_dend, g_allu, nrow = 1, widths = rel_widths)
    } else {
      gridExtra::arrangeGrob(g_allu, g_dend, nrow = 1, widths = rel_widths)
    }
    out <- stacked
  }

  # ---- Place legend outside ----
  if (legend == "separate" && !is.null(combined_legend)) {
    out <- if (legend_position == "right") {
      gridExtra::arrangeGrob(
        out, combined_legend,
        nrow = 1,
        widths = grid::unit.c(grid::unit(1, "null"), grid::unit(legend_rel_width, "null"))
      )
    } else {
      gridExtra::arrangeGrob(
        out, combined_legend,
        ncol = 1,
        heights = grid::unit.c(grid::unit(1, "null"), grid::unit(legend_rel_height, "null"))
      )
    }
  }

  # ---- Outer margins ----
  if (!is.null(outer_margins)) out <- add_outer_margins(out, outer_margins)
  out
}



#' Create an aligned alluvial + dendrogram figure (full workflow wrapper)
#'
#' High-level convenience wrapper that builds a composite figure combining:
#' \enumerate{
#'   \item A grouped ASV matrix (mean abundance per \code{group_col}),
#'   \item A hierarchical clustering dendrogram built on that grouped matrix,
#'   \item An alluvial (Sankey) plot of taxonomic composition across the same groups,
#'   \item A final aligned composite figure with optional external legend.
#' }
#'
#' This function is designed for reproducible figures where
#' clustering structure (top panel) and compositional structure (bottom panel)
#' share identical group ordering and horizontal alignment.
#'
#' @param asv_matrix Numeric matrix with samples as columns and taxa/ASVs as rows.
#'   Column names must match \code{metadata[[sample_col]]}.
#' @param metadata Data frame with one row per sample. Must contain
#'   \code{sample_col} and \code{group_col}.
#' @param sample_col Character string specifying the sample ID column in \code{metadata}
#'   (default: \code{"SampleID"}).
#' @param group_col Character string specifying the grouping variable used both for:
#'   \itemize{
#'     \item averaging samples in the ASV matrix (clustering),
#'     \item defining the x-axis of the alluvial plot.
#'   }
#' @param alluvial_data Long-format data used for the alluvial workflow
#'   (typically ASV-level relative abundance data).
#' @param tax_level Taxonomic column used for aggregation in the alluvial plot
#'   (e.g. \code{"Family"}, \code{"Class"}).
#' @param group_order Determines group ordering in both dendrogram and alluvial plot.
#'   One of:
#'   \itemize{
#'     \item \code{"metadata"} – first-appearance order in metadata,
#'     \item \code{"alphabetical"} – sorted group names,
#'     \item a character vector of explicit group names.
#'   }
#' @param low_abundance_threshold Numeric threshold passed to
#'   \code{\link{classify_taxa_patterns}} via \code{\link{create_alluvial_plot}}
#'   (default: 0.01).
#' @param palette_list Character vector of palette names passed to
#'   \code{\link{generate_alluvial_palette}} (default: \code{NULL}).
#' @param custom_palette Optional named color vector. If provided,
#'   automatic palette generation is skipped.
#' @param dend_color_by Column in \code{metadata} used to color dendrogram tips.
#'   Defaults to \code{group_col}.
#' @param dend_color_palette Named vector mapping levels of \code{dend_color_by}
#'   to colors.
#' @param dend_shape_by Optional column in \code{metadata} used to set
#'   dendrogram tip shapes (default: \code{NULL}).
#' @param dend_shape_palette Optional vector passed to
#'   \code{ggplot2::scale_shape_manual()} (default: \code{NULL}).
#' @param theme_alluvial ggplot2 theme added to the alluvial plot
#'   (default: \code{ggplot2::theme_minimal()}).
#' @param theme_dendrogram ggplot2 theme added to the dendrogram plot
#'   (default: \code{ggplot2::theme_void()}).
#' @param build_dendrogram_args Named list forwarded to
#'   \code{\link{build_dendrogram}} (e.g. \code{distance_method}, \code{cluster_method}).
#' @param plot_dendrogram_args Named list forwarded to
#'   \code{\link{plot_dendrogram}} (e.g. \code{orientation}, \code{point_size}, \code{guides_obj}).
#' @param alluvial_args Named list forwarded to
#'   \code{\link{create_alluvial_plot}}. May include routed argument lists
#'   (\code{prepare_args}, \code{classify_args}, \code{palette_args}, \code{plot_args})
#'   and \code{return_all}.
#' @param combine_args Named list forwarded to
#'   \code{\link{combine_dendrogram_alluvial}}
#'   (e.g. \code{dend_height}, \code{legend}, \code{legend_position}).
#' @param post_plot_guides Modify legend guides of the plot (default: NULL).
#' @param dend_limits_left,dend_limits_right “sweet spot” for dendrogram x view:
#'   dend xlim = c(dend_limits_left, n + dend_limits_right). Use e.g. -0.75 and 0.2.
#' @param dend_limits_bottom,dend_limits_top “sweet spot” for dendrogram y view:
#'   dend ylim = c(dend_limits_bottom, n + dend_limits_top).
#'
#' #' @section Workflow overview:
#' Internally, the function performs:
#' \enumerate{
#'   \item \code{\link{create_grouped_matrix}} – average ASV abundances within \code{group_col}.
#'   \item \code{\link{build_dendrogram}} – compute distance + hierarchical clustering.
#'   \item \code{\link{plot_dendrogram}} – render dendrogram as a ggplot.
#'   \item \code{\link{create_alluvial_plot}} – run the complete alluvial pipeline
#'         (prepare → classify → palette → plot).
#'   \item \code{\link{combine_dendrogram_alluvial}} – align and combine both panels.
#' }
#'
#' @section Argument ownership and routing:
#' This wrapper \strong{owns} the structural arguments passed to
#' \code{\link{create_alluvial_plot}}:
#' \code{data}, \code{tax_level}, \code{group_col}, \code{SampleID_col},
#' \code{groups}, \code{low_abundance_threshold}, \code{palette_list},
#' and \code{custom_palette}.
#'
#' Any duplicated arguments provided inside \code{alluvial_args}
#' (including nested \code{prepare_args}, \code{classify_args},
#' \code{palette_args}, or \code{plot_args}) are automatically removed
#' to prevent argument-matching conflicts.
#'
#' @return A named list containing:
#' \describe{
#'   \item{grouped_matrix}{Matrix used for clustering (rows = taxa, columns = groups).}
#'   \item{group_levels}{Final group ordering applied to both panels.}
#'   \item{dendrogram}{List returned by \code{\link{build_dendrogram}}.}
#'   \item{dendrogram_plot}{Dendrogram rendered as a ggplot.}
#'   \item{alluvial}{Output of \code{\link{create_alluvial_plot}}
#'         (list if \code{return_all=TRUE}).}
#'   \item{combined_plot}{Final composite grob suitable for \code{ggsave()}.}
#' }
#'
#' @details
#' This wrapper guarantees that:
#' \itemize{
#'   \item The dendrogram is computed on the same group structure used
#'         for the alluvial x-axis.
#'   \item Group ordering is identical across clustering and composition panels.
#'   \item Horizontal panel alignment is preserved in the final composite.
#' }
#'
#' @examples
#' # See the package vignette for a full reproducible example including
#' # ASV matrix, metadata, and long-format alluvial input.
#'
#' @export
create_alluvial_dendrogram_plot <- function(asv_matrix,
                                            metadata,
                                            sample_col = "SampleID",
                                            group_col,
                                            alluvial_data,
                                            tax_level,
                                            group_order = c("metadata", "alphabetical"),
                                            low_abundance_threshold = 0.01,
                                            palette_list = NULL,
                                            custom_palette = NULL,
                                            dend_color_by = NULL,
                                            dend_color_palette,
                                            dend_shape_by = NULL,
                                            dend_shape_palette = NULL,
                                            theme_alluvial = ggplot2::theme_minimal(),
                                            theme_dendrogram = ggplot2::theme_void(),
                                            build_dendrogram_args = list(),
                                            plot_dendrogram_args = list(),
                                            alluvial_args = list(),
                                            combine_args = list(),
                                            post_plot_guides = NULL,
                                            dend_limits_left = NULL,
                                            dend_limits_right = NULL,
                                            dend_limits_top = NULL,
                                            dend_limits_bottom = NULL) {
  # ---------------------------
  # Helpers
  # ---------------------------
  .trim_chr <- function(x) {
    x <- as.character(x)
    x <- trimws(x)
    x[!is.na(x)]
  }

  .drop_names_recursive <- function(x, drop) {
    if (!is.list(x)) {
      return(x)
    }
    for (nm in drop) {
      if (!is.null(x[[nm]])) x[[nm]] <- NULL
    }
    # also drop in nested routed args (if they exist)
    routed <- c("prepare_args", "classify_args", "palette_args", "plot_args")
    for (r in routed) {
      if (is.list(x[[r]])) {
        for (nm in drop) {
          if (!is.null(x[[r]][[nm]])) x[[r]][[nm]] <- NULL
        }
      }
    }
    x
  }

  # ---------------------------
  # Validate basics
  # ---------------------------
  if (!is.matrix(asv_matrix)) stop("asv_matrix must be a matrix.")
  if (!sample_col %in% names(metadata)) stop("sample_col not found in metadata: ", sample_col)
  if (!group_col %in% names(metadata)) stop("group_col not found in metadata: ", group_col)

  if (is.null(dend_color_by)) dend_color_by <- group_col
  if (!dend_color_by %in% names(metadata)) stop("dend_color_by not found in metadata: ", dend_color_by)

  # ---------------------------
  # 1) Group ASV matrix -> grouped matrix
  # ---------------------------
  gp_mat <- create_grouped_matrix(
    asv_matrix = asv_matrix,
    metadata = metadata,
    sample_col = sample_col,
    group_col = group_col
  )

  if (ncol(gp_mat) == 0) stop("Grouped matrix has 0 columns. Check group_col/sample_col mapping.")

  # ---------------------------
  # 1b) Determine & apply group order
  # ---------------------------
  # Handle group_order being passed as c("metadata","alphabetical") by taking the first
  if (is.character(group_order) && length(group_order) > 1 &&
    all(group_order %in% c("metadata", "alphabetical"))) {
    group_order <- group_order[1]
  }

  mat_groups <- .trim_chr(colnames(gp_mat))

  if (is.character(group_order) && length(group_order) == 1 && group_order %in% c("metadata", "alphabetical")) {
    if (group_order == "metadata") {
      ord <- .trim_chr(unique(metadata[[group_col]]))
      group_levels <- intersect(ord, mat_groups)
      # fallback: if metadata-derived order doesn't intersect, keep matrix col order
      if (length(group_levels) == 0) group_levels <- mat_groups
    } else {
      group_levels <- sort(mat_groups)
    }
  } else if (is.character(group_order) && length(group_order) >= 1) {
    group_levels <- intersect(.trim_chr(group_order), mat_groups)
    # fallback: if explicit vector doesn't intersect, keep matrix col order
    if (length(group_levels) == 0) group_levels <- mat_groups
  } else {
    stop("group_order must be 'metadata', 'alphabetical', or a character vector of group names.")
  }

  # apply ordering to matrix (use original colnames matching)
  # we match trimmed -> original to avoid whitespace issues
  map_trim_to_orig <- stats::setNames(colnames(gp_mat), .trim_chr(colnames(gp_mat)))
  gp_mat <- gp_mat[, map_trim_to_orig[group_levels], drop = FALSE]

  # ---------------------------
  # 2) Dendrogram
  # ---------------------------
  dend_call <- modifyList(
    list(mat = gp_mat, distance_method = "bray", cluster_method = "ward.D2", hang = -1),
    build_dendrogram_args
  )
  dend <- do.call(build_dendrogram, dend_call)

  # ---------------------------
  # 3) Dendrogram plot (labels are group names)
  # ---------------------------
  dend_plot_call <- modifyList(
    list(
      dend = dend,
      metadata = metadata,
      label_from = group_col,
      color_by = dend_color_by,
      color_palette = dend_color_palette,
      shape_by = dend_shape_by,
      shape_palette = dend_shape_palette,
      orientation = "top"
    ),
    plot_dendrogram_args
  )
  dend_plot <- do.call(plot_dendrogram, dend_plot_call) + theme_dendrogram

  # ---------------------------
  # 4) Alluvial plot via create_alluvial_plot()
  #    Avoid duplicated args by stripping wrapper-owned names from alluvial_args (incl nested)
  # ---------------------------
  if (is.null(alluvial_args$return_all)) alluvial_args$return_all <- TRUE

  drop_structural <- c(
    "data", "tax_level", "group_col", "SampleID_col", "groups",
    "low_abundance_threshold", "palette_list", "custom_palette"
  )
  alluvial_args <- .drop_names_recursive(alluvial_args, drop_structural)

  allu_call <- c(
    list(
      data = alluvial_data,
      tax_level = tax_level,
      group_col = group_col,
      SampleID_col = sample_col, # keep consistent naming
      groups = group_levels,
      low_abundance_threshold = low_abundance_threshold,
      palette_list = palette_list,
      custom_palette = custom_palette
    ),
    alluvial_args
  )

  allu_res <- do.call(create_alluvial_plot, allu_call)

  # theme on the plot itself (keeps routed-args clean)
  if (is.list(allu_res) && !is.null(allu_res$plot)) {
    p_allu <- allu_res$plot + theme_alluvial
    allu_res$plot <- p_allu
  } else {
    p_allu <- allu_res + theme_alluvial
    allu_res <- list(plot = p_allu)
  }
  if (!is.null(post_plot_guides)) {
    p_allu <- p_allu + do.call(ggplot2::guides, post_plot_guides)
  }
  # ---------------------------
  # 5) Combine
  # ---------------------------
combine_defaults <- list(
    dend_position      = "top",
    dend_height        = 0.15,
    strip_alluvial_x   = FALSE,
    legend             = "separate",
    legend_source      = "both",
    legend_position    = "right",
    legend_rel_width   = 0.32,
    dendrogram_margins = ggplot2::margin(0.25, 0.25, 0.25, 0.25),
    align_x_centers    = TRUE,
    overwrite_x_scales = TRUE,
    leaf_order         = dend$order
  )
  if (!is.null(dend_limits_left))  combine_defaults$dend_limits_left  <- dend_limits_left
  if (!is.null(dend_limits_right)) combine_defaults$dend_limits_right <- dend_limits_right
  combine_call <- modifyList(combine_defaults, combine_args)

  combined <- do.call(
    combine_dendrogram_alluvial,
    c(list(alluvial_plot = p_allu, dendrogram_plot = dend_plot), combine_call)
  )

  list(
    grouped_matrix = gp_mat,
    group_levels = group_levels,
    dendrogram = dend,
    dendrogram_plot = dend_plot,
    alluvial = allu_res,
    combined_plot = combined
  )
}
