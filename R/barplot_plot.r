#' Create stacked taxonomic barplots with optional faceting and strip theming
#'
#' Builds a stacked barplot of taxonomic composition (relative abundance) from
#' long-format data (typically produced by \code{\link{process_barplot_data}}).
#' The function maps \code{RA} to bar height and \code{tax_level} to fill color,
#' applies a user-supplied palette, and optionally facets the plot. If
#' \code{facet_strip_colors} is provided, facet-strip background colors are applied
#' via \pkg{ggh4x} (if installed).
#'
#' @param data A data frame in long format that must contain:
#'   \itemize{
#'     \item \code{RA}: relative abundance values (ideally in 0--1).
#'     \item the column named by \code{tax_level}: taxon labels.
#'     \item the column named by \code{x_axis_var} (or the first column if \code{x_axis_var = NULL}).
#'   }
#'   This is typically the output of \code{\link{process_barplot_data}}.
#' @param tax_level Character string giving the taxonomy column to plot as fill
#'   (e.g. \code{"Family"}, \code{"Class"}).
#' @param palette Named character vector of colors. Names must match the taxon
#'   labels in \code{data[[tax_level]]}. The order of \code{names(palette)} is used
#'   to set the factor levels for stacking and legend order.
#' @param x_axis_var Character string giving the x-axis column name. If \code{NULL},
#'   the first column of \code{data} is used (a message is printed).
#' @param facet_by Optional faceting specification. Either:
#'   \itemize{
#'     \item a single column name as a string (e.g. \code{"SampleType"}), interpreted as \code{. ~ SampleType}, or
#'     \item a formula string (e.g. \code{"Site ~ SampleType"} or \code{". ~ SampleType"}).
#'   }
#' @param facet_scales Facet scaling, passed to \code{facet_grid()} or
#'   \code{ggh4x::facet_grid2()}. Common values: \code{"fixed"}, \code{"free"},
#'   \code{"free_x"}, \code{"free_y"}.
#' @param facet_space Facet spacing, passed to \code{facet_grid()} or
#'   \code{ggh4x::facet_grid2()}. Common values: \code{"fixed"}, \code{"free"},
#'   \code{"free_x"}, \code{"free_y"}.
#' @param facet_strip_colors Optional named character vector of colors used to fill
#'   facet strip backgrounds. Requires \pkg{ggh4x}; if \pkg{ggh4x} is not installed,
#'   a warning is issued and default faceting is used. Names should correspond to
#'   facet levels (values of the faceting variable). Compatible with output from \code{\link{generate_grouped_palette}}.
#' @param show_x_labels Logical; if \code{FALSE} (default), x-axis tick labels and
#'   ticks are removed (useful when there are many samples). If \code{TRUE}, labels
#'   are shown and rotated using \code{x_label_angle}.
#' @param x_label_angle Numeric angle used for x-axis labels when
#'   \code{show_x_labels = TRUE} (default: 90).
#' @param add_black_outline Logical; if \code{TRUE} (default), draws a black outline
#'   around bars via \code{geom_col(color = "black")}.
#' @param line_width Numeric line width for bar outlines (default: 0.2). Only used
#'   when \code{add_black_outline = TRUE}.
#' @param theme_obj A ggplot2 theme object applied to the plot
#'   (default: \code{ggplot2::theme_minimal()}).
#' @param legend_title Legend title for the fill scale. If \code{NULL} (default),
#'   uses the value of \code{tax_level}.
#' @param y_label y-axis label (default: \code{"Relative abundance"}).
#' @param x_label x-axis label (default: empty string).
#'
#' @return A \code{ggplot} object.
#' 
#' @seealso \code{\link{process_barplot_data}}, \code{\link{generate_palette_hcl}},
#' \code{\link{generate_grouped_palette}}
#' 
#' @examples
#' library(ggplot2)
#' library(dplyr)
#'
#' # Minimal toy long-format dataset
#' toy <- tibble::tibble(
#'   SampleID = rep(c("S1","S2"), each = 3),
#'   Family   = rep(c("A","B","unknown"), times = 2),
#'   RA       = c(0.6, 0.3, 0.1,
#'               0.2, 0.7, 0.1)
#' )
#'
#' pal <- c(A = "#1b9e77", B = "#d95f02", unknown = "#000000")
#'
#' p <- plot_taxonomic_barplot(
#'   data = toy,
#'   tax_level = "Family",
#'   palette = pal,
#'   x_axis_var = "SampleID",
#'   show_x_labels = TRUE
#' )
#' p
#'
#' # Faceting (string column name -> interpreted as ". ~ SampleType")
#' toy2 <- toy %>% mutate(SampleType = rep(c("Soil","Water"), each = 3))
#' p2 <- plot_taxonomic_barplot(
#'   data = toy2,
#'   tax_level = "Family",
#'   palette = pal,
#'   x_axis_var = "SampleID",
#'   facet_by = "SampleType"
#' )
#' p2
#'
#' @export
plot_taxonomic_barplot <- function(data,
                                   tax_level,
                                   palette,
                                   x_axis_var = NULL,
                                   facet_by = NULL,
                                   facet_scales = "free",
                                   facet_space = "free_x",
                                   facet_strip_colors = NULL,
                                   show_x_labels = FALSE,
                                   x_label_angle = 90,
                                   add_black_outline = TRUE,
                                   line_width = 0.2,
                                   theme_obj = ggplot2::theme_minimal(),
                                   legend_title = NULL,
                                   y_label = "Relative abundance",
                                   x_label = "") {
  
  # Input validation
  if (!(tax_level %in% colnames(data))) {
    stop("tax_level '", tax_level, "' not found in data")
  }
  
  if (!"RA" %in% colnames(data)) {
    stop("Data must contain 'RA' column (relative abundance)")
  }
  
  # Auto-detect x_axis_var if not provided
  if (is.null(x_axis_var)) {
    x_axis_var <- colnames(data)[1]
    message("Using '", x_axis_var, "' for x-axis")
  }
  
  if (!(x_axis_var %in% colnames(data))) {
    stop("x_axis_var '", x_axis_var, "' not found in data")
  }
  
  # Set legend title
  if (is.null(legend_title)) {
    legend_title <- tax_level
  }
  
  # Set taxonomic levels according to palette
  tax_names <- names(palette)
  data[[tax_level]] <- factor(data[[tax_level]], levels = tax_names)
  
  # Create base plot
  p <- ggplot2::ggplot(data, ggplot2::aes(x = .data[[x_axis_var]], 
                        y = RA, 
                        fill = .data[[tax_level]])) +
    ggplot2::geom_col(color = if(add_black_outline) "black" else NA,
                      linewidth = if(add_black_outline) line_width else 0) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::scale_fill_manual(values = palette, limits = force, na.value = "black") +
    ggplot2::labs(y = y_label, x = x_label, fill = legend_title) +
    theme_obj
  
  # Handle x-axis labels
  if (!show_x_labels) {
    p <- p +
      ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      )
  } else {
    p <- p +
      ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(angle = x_label_angle))
  }
  
  # Add faceting if requested
  if (!is.null(facet_by)) {
    # Check if facet_by is a formula or column name(s)
    if (is.character(facet_by) && !grepl("~", facet_by)) {
      # Simple column name - create formula
      facet_formula <- formula(paste(". ~", facet_by))
      facet_var <- facet_by
    } else if (is.character(facet_by) && grepl("~", facet_by)) {
      # Already a formula string
      facet_formula <- formula(facet_by)
      # Extract variable name(s) for strip colors
      facet_var <- gsub(".*~ *", "", facet_by)
      facet_var <- trimws(strsplit(facet_var, "\\+")[[1]][1])
    } else {
      stop("facet_by must be a column name or formula string")
    }
    
    # Check if ggh4x is available for custom strip colors
    if (!is.null(facet_strip_colors)) {
      if (!requireNamespace("ggh4x", quietly = TRUE)) {
        warning("Package 'ggh4x' required for custom facet strip colors. Install with: install.packages('ggh4x'). Using default faceting instead.")
        facet_strip_colors <- NULL
      }
    }
    
    # Add faceting with or without custom strip colors
    if (!is.null(facet_strip_colors)) {
      # Match strip colors to facet levels
      facet_levels <- unique(data[[facet_var]])
      
      # Create ordered strip colors matching facet levels
      if (is.null(names(facet_strip_colors))) {
        strip_colors_ordered <- facet_strip_colors
      } else {
        strip_colors_ordered <- facet_strip_colors[as.character(facet_levels)]
        strip_colors_ordered <- strip_colors_ordered[!is.na(strip_colors_ordered)]
      }
      
      # With custom strip colors (requires ggh4x)
      p <- p +
        ggh4x::facet_grid2(
          facet_formula,
          scales = facet_scales,
          space = facet_space,
          strip = ggh4x::strip_themed(
            background_x = ggh4x::elem_list_rect(fill = strip_colors_ordered)
          )
        ) +
        ggplot2::theme(strip.text = ggplot2::element_text(color = "black"))
    } else {
      # Standard faceting
      p <- p +
        ggplot2::facet_grid(
          facet_formula,
          scales = facet_scales,
          space = facet_space,
          labeller = ggplot2::labeller(.default = ggplot2::label_wrap_gen(width = 12))
        )
    }
  }
  
  return(p)
}