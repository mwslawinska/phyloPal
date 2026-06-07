#' Add alpha transparency to hex colors
#'
#' @param hex Character vector of hex colors
#' @param alpha Numeric value between 0 (transparent) and 1 (opaque)
#'
#' @return Character vector of hex colors with alpha channel
#'
#' @examples
#' add_alpha("#FF0000", 0.5)  # Semi-transparent red
#' add_alpha(c("#FF0000", "#00FF00"), 0.7)
#'
#' @export
add_alpha <- function(hex, alpha = 1) {
  if (alpha < 0 || alpha > 1) {
    stop("alpha must be between 0 and 1")
  }
  
  # Convert hex to RGB
  rgb_vals <- grDevices::col2rgb(hex)
  
  # Add alpha and convert back to hex
  grDevices::rgb(
    rgb_vals[1, ], 
    rgb_vals[2, ], 
    rgb_vals[3, ],
    alpha = alpha * 255,
    maxColorValue = 255
  )
}

#' Generate HCL color palette for taxonomic visualization
#'
#' Generates a named HCL-based color palette for taxa at a specified
#' taxonomic level. Optionally supports hierarchical grouping, where
#' colors are assigned within higher taxonomic groups (e.g. coloring
#' families grouped by phylum).
#'
#' @param data A data frame containing at least \code{tax_level}. If
#'   \code{group_by_higher_tax} is used, the corresponding higher-taxonomy
#'   column must also be present.
#' @param tax_level Character string specifying the taxonomic level to
#'   generate colors for (e.g. \code{"Family"}, \code{"Genus"}).
#' @param group_by_higher_tax Optional higher taxonomic level used to group
#'   taxa before color assignment (e.g. \code{"Phylum"}, \code{"Class"}).
#'   Recommended only for small datasets (e.g. synthetic communities with
#'   fewer than ~10 higher-level taxa). For complex natural communities,
#'   the ungrouped palette is typically more interpretable.
#' @param group_palette_map Optional named list mapping higher taxa to palette
#'   specifications (required if \code{group_by_higher_tax} is used). Each
#'   element can be either:
#'   \itemize{
#'     \item a single palette name (character), or
#'     \item a list with elements \code{palette} and \code{side} for diverging palettes.
#'   }
#' @param group_default_side Default side for grouped diverging palettes:
#'   \code{"both"}, \code{"left"}, or \code{"right"} (default: \code{"both"}).
#' @param order_by_higher_tax Logical; if \code{TRUE}, taxa are ordered
#'   deterministically by higher taxonomy and then by \code{tax_level}
#'   (default: \code{TRUE} when \code{group_by_higher_tax} is used).
#' @param order_groups How to order higher-taxonomy groups when
#'   \code{order_by_higher_tax = TRUE}. Options:
#'   \code{"alphabetical"} (default) or \code{"data"}.
#' @param order_within_groups How to order taxa within each higher-taxonomy
#'   group when \code{order_by_higher_tax = TRUE}. Options:
#'   \code{"alphabetical"} (default) or \code{"data"}.
#' @param fixed_colors_enabled Logical; if \code{TRUE}, apply fixed colors
#'   to specific taxa (default: \code{TRUE}).
#' @param fixed_colors Named character vector mapping specific taxa to
#'   fixed colors (default: \code{c("unknown" = "#000000",
#'   "low abundant" = "#E5E5E5")}).
#' @param fixed_colors_position Position of fixed-color taxa in the
#'   resulting palette: \code{"end"} (default), \code{"beginning"},
#'   or \code{"alphabetical"}.
#' @param palette_list Character vector of HCL palette names used for
#'   dynamically generated colors. If \code{NULL}, a predefined set is used.
#' @param cmax Maximum chroma value for generated HCL colors (default: 100).
#' @param luminance Luminance value (single number) or range (length-2 numeric)
#'   used in HCL generation (default: 65).
#' @param power Power transformation applied in HCL color interpolation
#'   (default: 1.2).
#' @param shuffle Logical; if \code{TRUE}, dynamically generated colors are
#'   shuffled to reduce adjacency similarity (default: \code{TRUE}).
#' @param seed Numeric seed used when \code{shuffle = TRUE} (default: 42).
#'
#' @details
#' \strong{Color generation workflow:}
#' \enumerate{
#'   \item Identify unique taxa at \code{tax_level}.
#'   \item Optionally group taxa by \code{group_by_higher_tax}.
#'   \item Generate HCL-based colors either globally (ungrouped) or within
#'         each higher-taxonomy group.
#'   \item Apply ordering rules (if enabled).
#'   \item Insert fixed colors according to \code{fixed_colors_position}.
#' }
#'
#' When hierarchical grouping is enabled, each higher-taxonomy group
#' receives its own palette (defined in \code{group_palette_map}). This
#' produces visually coherent sub-palettes but may reduce distinguishability
#' in highly diverse datasets.
#'
#' @return A named character vector where names correspond to taxa at
#'   \code{tax_level} and values are hexadecimal color codes.
#'
#' @seealso \code{\link{process_barplot_data}}, \code{\link{plot_taxonomic_barplot}}
#' 
#' @examples
#' \dontrun{
#' library(dplyr)
#'
#' toy <- tibble::tibble(
#'   Phylum = c("P1","P1","P2","P2"),
#'   Family = c("F1","F2","F3","F4")
#' )
#'
#' # Simple ungrouped palette
#' pal_simple <- generate_tax_palette(
#'   data = toy,
#'   tax_level = "Family"
#' )
#'
#' # Hierarchical grouping by Phylum
#' pal_grouped <- generate_tax_palette(
#'   data = toy,
#'   tax_level = "Family",
#'   group_by_higher_tax = "Phylum",
#'   group_palette_map = list(
#'     P1 = "Blues",
#'     P2 = "Reds"
#'   )
#' )
#'
#' pal_simple
#' pal_grouped
#'}
#' @export
generate_palette_hcl <- function(data,
                                 tax_level,
                                 group_by_higher_tax = NULL,
                                 group_palette_map = NULL,
                                 group_default_side = "both",
                                 order_by_higher_tax = FALSE,
                                 order_groups = c("alphabetical", "data"),
                                 order_within_groups = c("alphabetical", "data"),
                                 fixed_colors_enabled = TRUE,
                                 fixed_colors = c("unknown" = "#000000", "low abundant" = "#E5E5E5"),
                                 fixed_colors_position = "end",
                                 palette_list = NULL,
                                 cmax = 100,
                                 luminance = 65,
                                 power = 1.2,
                                 shuffle = TRUE,
                                 seed = 42) {

  order_groups <- match.arg(order_groups)
  order_within_groups <- match.arg(order_within_groups)

  # -----------------------------
  # Input validation
  # -----------------------------
  required_cols <- c(tax_level)
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  if (!(fixed_colors_position %in% c("end", "beginning", "alphabetical"))) {
    stop("Invalid fixed_colors_position. Must be one of 'end', 'beginning', or 'alphabetical'.")
  }

  if (!is.null(group_by_higher_tax) && is.null(group_palette_map)) {
    stop("group_palette_map must be provided when using group_by_higher_tax")
  }

  if (!is.null(group_by_higher_tax)) {
    if (!(group_by_higher_tax %in% colnames(data))) {
      stop("group_by_higher_tax '", group_by_higher_tax, "' not found in data")
    }

    n_higher_taxa <- dplyr::n_distinct(data[[group_by_higher_tax]], na.rm = TRUE)
    if (n_higher_taxa > 10) {
      warning("You have ", n_higher_taxa, " unique values in '", group_by_higher_tax,
              "'. Hierarchical grouping is recommended for <10 higher taxa. ",
              "Consider using ungrouped palette for better color distinction.")
    }
  }

  # Auto-enable order_by_higher_tax if group_by_higher_tax is used
  if (!is.null(group_by_higher_tax) && !order_by_higher_tax) {
    order_by_higher_tax <- TRUE
  }

  # Get unique taxa names
  tax_names <- unique(data[[tax_level]])
  tax_names[is.na(tax_names)] <- "unknown"

  # ==================================================================
  # HIERARCHICAL GROUPING PATH
  # ==================================================================
  if (!is.null(group_by_higher_tax) && !is.null(group_palette_map)) {

    # Determine which taxa get fixed colors vs dynamic colors
    if (fixed_colors_enabled) {
      dynamic_tax_names <- tax_names[!(tax_names %in% names(fixed_colors))]
      fixed_tax <- names(fixed_colors)[names(fixed_colors) %in% tax_names]
      fixed_colors_used <- fixed_colors[fixed_tax]
    } else {
      dynamic_tax_names <- tax_names
      fixed_tax <- character(0)
      fixed_colors_used <- character(0)
    }

    # Data subset used for grouping + ordering
    data_for_grouping <- data %>%
      dplyr::filter(.data[[tax_level]] %in% dynamic_tax_names) %>%
      dplyr::select(dplyr::all_of(c(group_by_higher_tax, tax_level))) %>%
      dplyr::distinct()

    # Remove rows where group_by_higher_tax is also a fixed color key (rare but defensive)
    if (fixed_colors_enabled) {
      data_for_grouping <- data_for_grouping %>%
        dplyr::filter(!(.data[[group_by_higher_tax]] %in% names(fixed_colors)))
    }

    # Check palette assignments
    higher_taxa_in_data <- unique(data_for_grouping[[group_by_higher_tax]])
    missing_assignments <- setdiff(higher_taxa_in_data, names(group_palette_map))
    if (length(missing_assignments) > 0) {
      stop("Missing palette assignments for: ", paste(missing_assignments, collapse = ", "),
           ". Add them to group_palette_map.")
    }

    # Generate grouped palette for dynamic taxa
    dynamic_colors <- generate_grouped_palette(
      data = data_for_grouping,
      group_col = group_by_higher_tax,
      item_col = tax_level,
      palette_map = group_palette_map,
      default_side = group_default_side,
      cmax = cmax,
      luminance = luminance,
      power = power,
      shuffle = shuffle,
      seed = seed
    )

    # ---- NEW: deterministic ordering of the PALETTE names ----
    if (isTRUE(order_by_higher_tax)) {

      ord_df <- data_for_grouping %>%
        dplyr::mutate(
          .group = as.character(.data[[group_by_higher_tax]]),
          .tax   = as.character(.data[[tax_level]])
        ) %>%
        dplyr::distinct(.group, .tax)

      # Order groups
      if (order_groups == "alphabetical") {
        ord_df <- ord_df %>% dplyr::mutate(.group = factor(.group, levels = sort(unique(.group))))
      } else {
        ord_df <- ord_df %>% dplyr::mutate(.group = factor(.group, levels = unique(.group)))
      }

      # Order taxa within groups
      if (order_within_groups == "alphabetical") {
        ord_df <- ord_df %>% dplyr::arrange(.group, .tax)
      } else {
        # "data" order: keep first appearance order within each group
        ord_df <- ord_df %>%
          dplyr::group_by(.group) %>%
          dplyr::mutate(.tax = factor(.tax, levels = unique(.tax))) %>%
          dplyr::ungroup() %>%
          dplyr::arrange(.group, .tax)
      }

      ordered_tax <- ord_df$.tax %>% unique()

      # Reorder dynamic_colors to match the computed order (keep any extras at the end)
      present <- intersect(ordered_tax, names(dynamic_colors))
      extra <- setdiff(names(dynamic_colors), present)
      dynamic_colors <- dynamic_colors[c(present, extra)]
    }

    # Combine with fixed colors based on position
    if (!fixed_colors_enabled) {
      col_tax <- dynamic_colors
    } else {
      if (fixed_colors_position == "end") {
        col_tax <- c(dynamic_colors, fixed_colors_used)
      } else if (fixed_colors_position == "beginning") {
        col_tax <- c(fixed_colors_used, dynamic_colors)
      } else if (fixed_colors_position == "alphabetical") {
        col_tax <- c(dynamic_colors, fixed_colors_used)
        col_tax <- col_tax[sort(names(col_tax))]
      }
    }

    return(col_tax)
  }

  # ==================================================================
  # STANDARD UNGROUPED PATH
  # ==================================================================

  # Default palettes
  default_palettes <- c("RedYel", "Emrld", "Red-Blue", "YlOrRd", "Purp", "Peach", "Blues")
  if (is.null(palette_list) || length(palette_list) == 0) {
    palette_list <- default_palettes
  }

  # Validate palette names
  available_palettes <- grDevices::hcl.pals()
  wrong_palettes <- setdiff(palette_list, available_palettes)
  if (length(wrong_palettes) > 0) {
    stop("Invalid palette names: ", paste(wrong_palettes, collapse = ", "),
         ". Available palettes are: ", paste(available_palettes, collapse = ", "))
  }

  # Determine which taxa get fixed colors vs dynamic colors
  if (fixed_colors_enabled) {
    dynamic_tax <- tax_names[!(tax_names %in% names(fixed_colors))]
    fixed_tax <- names(fixed_colors)[names(fixed_colors) %in% tax_names]
    fixed_colors_used <- fixed_colors[fixed_tax]
  } else {
    dynamic_tax <- tax_names
    fixed_tax <- character(0)
    fixed_colors_used <- character(0)
  }

  # Generate dynamic colors if needed
  if (length(dynamic_tax) > 0) {
    n_dyn <- length(dynamic_tax)
    n_pal <- length(palette_list)
    base_n <- floor(n_dyn / n_pal)
    extra <- n_dyn %% n_pal

    if (shuffle) set.seed(seed)

    dynamic_colors <- unlist(lapply(seq_along(palette_list), function(i) {
      n_i <- base_n + ifelse(i <= extra, 1, 0)
      if (n_i > 0) {
        cols <- colorspace::sequential_hcl(
          n = n_i,
          palette = palette_list[i],
          cmax = cmax,
          l = luminance,
          power = power
        )
        if (shuffle) sample(cols) else cols
      } else NULL
    }))

    if (shuffle) dynamic_colors <- sample(dynamic_colors)
    names(dynamic_colors) <- sort(dynamic_tax)
  } else {
    dynamic_colors <- character(0)
  }

  # Combine colors based on position preference
  if (!fixed_colors_enabled) {
    col_tax <- dynamic_colors
  } else if (fixed_colors_position == "end") {
    col_tax <- c(dynamic_colors, fixed_colors_used)
  } else if (fixed_colors_position == "beginning") {
    col_tax <- c(fixed_colors_used, dynamic_colors)
  } else if (fixed_colors_position == "alphabetical") {
    col_tax <- c(dynamic_colors, fixed_colors_used)
    col_tax <- col_tax[sort(names(col_tax))]
  }

  return(col_tax)
}

#' Generate grouped palettes using HCL colors
#'
#' Generates a named vector of colors for items grouped by a categorical variable.
#' Each group is assigned a color "family" based on an HCL palette name, and colors
#' for individual items within that group are generated by varying hue, chroma,
#' and luminance around a base hue derived from the palette.
#'
#' @param data A data frame containing grouping and item columns.
#' @param group_col Character string specifying the column containing group labels.
#' @param item_col Character string specifying the column containing item names
#'   (each item receives a unique color).
#' @param palette_map Named list mapping group names to HCL palette names
#'   (see \code{grDevices::hcl.pals()}). Example:
#'   \code{list(GroupA = "YlOrBr", GroupB = "BuPu", GroupC = "Greens")}.
#' @param default_side Character string specifying which side of a diverging palette
#'   to use when determining the base hue. One of \code{"both"}, \code{"left"},
#'   \code{"right"}, \code{"low"}, or \code{"high"}. (default: "both")
#' @param cmax Numeric value controlling maximum chroma when sampling the base palette
#'   to determine a representative hue. (default: 80)
#' @param luminance Numeric value or length-2 vector controlling luminance.
#'   If length 2, defines the luminance range used to generate colors
#'   (e.g. \code{c(45, 80)}). (default: c(40, 80))
#' @param power Numeric value controlling the shape of sequential/diverging
#'   HCL palettes when deriving the base hue. (default: 1.2)
#' @param shuffle Logical; if \code{TRUE}, randomly shuffles colors within each group. (default: FALSE)
#' @param seed Integer seed used for reproducible shuffling. (default: 42)
#' @param increase_separation Logical indicating whether to increase perceptual
#'   separation between colours within a group by spreading hues and alternating
#'   luminance/chroma values. (default: FALSE)
#' @param hue_span Numeric controlling the width (in degrees) of the hue band
#'   around the palette's base hue used to generate colors within a group. (default: 25)
#' @param c_range Numeric vector of length 2 specifying the chroma range used
#'   when generating colors (default values typically between 20 and 80 give
#'   visually pleasing results). (default: c(25, 85))
#' @param l_range Numeric vector of length 2 specifying the luminance range used
#'   when generating colors. Restricting luminance (e.g. \code{c(45, 80)}) helps
#'   avoid colors that are too dark or too close to white. (default: c(40, 80))
#'
#' @return A named character vector of hex colors, where names correspond to
#'   the unique values of \code{item_col}.
#'
#' @details
#' Colors are generated in the HCL (Hue–Chroma–Luminance) color space using
#' \code{\link[colorspace]{polarLUV}}. Compared with RGB or HSV color models,
#' HCL is designed to better reflect human color perception, meaning that
#' changes in luminance and chroma correspond more closely to perceived
#' differences between colors. This makes HCL particularly suitable for
#' data visualisation where many categories must be distinguishable.
#'
#' For each group, a base hue is derived from the specified HCL palette.
#' Colors for individual items are then generated within a narrow hue band
#' around this base hue while varying luminance and chroma. When
#' \code{increase_separation = TRUE}, the function spreads hues slightly and
#' alternates luminance and chroma values to increase perceptual separation
#' between neighbouring colors.
#'
#' Colors are generated in HCL space using \code{\link[colorspace]{polarLUV}} and
#' converted to hex values. If a color falls outside the displayable sRGB gamut,
#' chroma is reduced automatically until a valid color is obtained.
#' 
#' When \code{increase_separation = TRUE}, the function increases perceptual
#' differences between neighbouring colors by spreading hues slightly and
#' alternating luminance and chroma values. This typically improves readability
#' when many items must be distinguished within the same group.
#
#' @examples
#' \dontrun{
#' example_data <- data.frame(
#'   Type = rep(c("soil", "moss"), each = 5),
#'   Description = paste0("item_", 1:10)
#' )
#'
#' pal <- generate_grouped_palette(
#'   data = example_data,
#'   group_col = "Type",
#'   item_col = "Description",
#'   palette_map = list(
#'     soil = "YlOrBr",
#'     moss = "BuPu"
#'   ),
#'   increase_separation = TRUE,
#'   hue_span = 30,
#'   c_range = c(25, 70),
#'   luminance = c(45, 80)
#' )
#'
#' print(pal)
#' }
#' @export
generate_grouped_palette <- function(data,
                                     group_col,
                                     item_col,
                                     palette_map,
                                     default_side = "both",
                                     cmax = 80,
                                     luminance = c(40, 80),
                                     power = 1.2,
                                     shuffle = FALSE,
                                     seed = 42,
                                     increase_separation = FALSE,
                                     hue_span = 25,
                                     c_range = c(25, 85),
                                     l_range = c(40, 80)) {

  # --- validation (light) ---
  if (!group_col %in% names(data)) stop("group_col '", group_col, "' not found in data")
  if (!item_col %in% names(data)) stop("item_col '", item_col, "' not found in data")
  if (!is.list(palette_map)) stop("palette_map must be a named list")
  if (!default_side %in% c("both", "left", "right", "low", "high")) {
    stop("default_side must be one of: 'both', 'left', 'right', 'low', 'high'")
  }
  if (!is.numeric(luminance) || !length(luminance) %in% c(1, 2)) {
    stop("luminance must be numeric length 1 or 2 (e.g. 65 or c(40,80))")
  }

  # derive a base hue from the *middle* of the named HCL palette
  palette_base_hue <- function(palette_name, default_side) {
    diverging <- grDevices::hcl.pals(type = "diverging")
    sequential <- grDevices::hcl.pals(type = "sequential")
    qualitative <- grDevices::hcl.pals(type = "qualitative")

    if (palette_name %in% diverging) {
      pal <- colorspace::diverging_hcl(n = 9, palette = palette_name, cmax = 100, l = 60, power = power)
      # choose side-ish base if requested
      mid <- if (default_side %in% c("left", "low")) pal[3] else if (default_side %in% c("right", "high")) pal[7] else pal[5]
    } else if (palette_name %in% sequential) {
      pal <- colorspace::sequential_hcl(n = 9, palette = palette_name, cmax = 100, l = 60, power = power)
      mid <- pal[5]
    } else if (palette_name %in% qualitative) {
      pal <- colorspace::qualitative_hcl(n = 9, palette = palette_name, l = 60)
      mid <- pal[1]
    } else {
      stop("Palette '", palette_name, "' not found in grDevices::hcl.pals()")
    }

    hcl <- methods::as(colorspace::hex2RGB(mid), "polarLUV")
    h <- as.numeric(hcl@coords[1, "H"])
    if (is.na(h)) h <- 0
    h
  }

  # create a “nice” set of H/C/L values for n items
 make_hcl_set <- function(n, base_h, hue_span, c_range, l_range, increase_separation) {
  if (n == 1) {
    H <- base_h %% 360
    C <- mean(c_range)
    L <- mean(l_range)
    return(colorspace::hex(colorspace::polarLUV(L = L, C = C, H = H), fixup = TRUE))
  }

  span <- if (increase_separation) hue_span else max(5, hue_span / 2)
  H <- seq(base_h - span, base_h + span, length.out = n) %% 360

  Lseq <- seq(l_range[1], l_range[2], length.out = n)
  Cseq <- seq(c_range[2], c_range[1], length.out = n)  # opposite direction

  zigzag <- function(x) {
    odd <- x[seq(1, length(x), by = 2)]
    even <- rev(x[seq(2, length(x), by = 2)])
    c(odd, even)[seq_len(length(x))]
  }

  if (increase_separation) {
    L <- zigzag(Lseq)
    C <- zigzag(Cseq)
  } else {
    L <- Lseq
    C <- Cseq
  }

  # Try conversion; if any NA remain, progressively reduce chroma
  for (k in 0:6) {
    Ck <- C * (0.90^k)  # reduce chroma stepwise if needed
    cols <- colorspace::hex(colorspace::polarLUV(L = L, C = Ck, H = H), fixup = TRUE)
    if (!anyNA(cols)) return(cols)
  }

  # Last resort: force any remaining NA to a neutral grey (should be rare)
  cols[is.na(cols)] <- "#808080"
  cols
}

  groups <- unique(data[[group_col]])
  missing_groups <- setdiff(groups, names(palette_map))
  if (length(missing_groups) > 0) {
    stop("Missing palette assignments for groups: ", paste(missing_groups, collapse = ", "))
  }

  if (shuffle) set.seed(seed)

  all_colors <- list()

  for (g in groups) {
    items <- data |>
      dplyr::filter(.data[[group_col]] == g) |>
      dplyr::pull(.data[[item_col]]) |>
      unique()

    n <- length(items)
    if (n == 0) next

    pal_spec <- palette_map[[g]]
    if (is.list(pal_spec)) {
      pal_name <- pal_spec$palette
      side <- if (!is.null(pal_spec$side)) pal_spec$side else default_side
    } else {
      pal_name <- pal_spec
      side <- default_side
    }

    base_h <- palette_base_hue(pal_name, side)

    cols <- make_hcl_set(
      n = n,
      base_h = base_h,
      hue_span = hue_span,
      c_range = c_range,
      l_range = if (length(luminance) == 2) luminance else l_range,
      increase_separation = increase_separation
    )

    if (shuffle) cols <- sample(cols)

    names(cols) <- items
    all_colors[[g]] <- cols
  }

  out <- unlist(all_colors, use.names = FALSE)
  names(out) <- unlist(lapply(all_colors, names), use.names = FALSE)
  out
}

#' Generate a color palette for alluvial plots from tax_color keys
#'
#' Creates a named color palette for alluvial/Sankey plots where the palette keys
#' are taken from \code{data$tax_color}. Three keys are treated as fixed categories
#' by default (\code{"unique low abundant"}, \code{"unknown"}, \code{"shared low abundant"}),
#' taxa listed in \code{special_taxa} receive their own colors, and all remaining
#' keys receive automatically generated HCL colors.
#'
#' @param data A data frame that must contain a \code{tax_color} column.
#'   The returned palette names correspond to the unique values of \code{data$tax_color}.
#' @param tax_level Character string (currently unused). Included for API compatibility
#'   with higher-level workflow wrappers.
#' @param group_by_higher_tax Optional higher taxonomic level used to group
#'   taxa before color assignment (e.g. \code{"Phylum"}, \code{"Class"}).
#'   Recommended only for small datasets (e.g. synthetic communities with
#'   fewer than ~10 higher-level taxa). For complex natural communities,
#'   the ungrouped palette is typically more interpretable.
#' @param group_palette_map Optional named list mapping higher taxa to palette
#'   specifications (required if \code{group_by_higher_tax} is used). Each
#'   element can be either:
#'   \itemize{
#'     \item a single palette name (character), or
#'     \item a list with elements \code{palette} and \code{side} for diverging palettes.
#'   }
#' @param group_default_side Default side for grouped diverging palettes:
#'   \code{"both"}, \code{"left"}, or \code{"right"} (default: \code{"both"}).
#' @param special_taxa Character vector of palette keys (i.e., values of \code{tax_color})
#'   that should receive their own colors (drawn from \code{palette = "Dark 3"}) and be
#'   excluded from the automatically generated pool.
#' @param special_side Controls ordering of \code{special_taxa}. Currently only
#'   \code{"alphabetical"} has an effect (special taxa are sorted); other values are
#'   accepted for forward compatibility.
#' @param palette_list Character vector of \pkg{colorspace} sequential HCL palette names
#'   used to generate colors for non-fixed, non-special keys. If \code{NULL} or empty,
#'   defaults to \code{c("RedYel","Emrld","Red-Blue","YlOrRd","Purp","Peach","Blues")}.
#' @param custom_palette Optional named character vector of colors. If provided, the
#'   function validates that it covers all unique values of \code{data$tax_color} and
#'   returns it unchanged.
#' @param fixed_colors Named character vector of fixed colors for reserved keys.
#'   Names must match values used in \code{data$tax_color}. Defaults to:
#'   \code{c("unique low abundant"="#A5A5A5","unknown"="#000000","shared low abundant"="#E5E5E5")}.
#' @param cmax Maximum chroma for generated HCL colors (passed to
#'   \code{colorspace::sequential_hcl()}, default: 100).
#' @param luminance Luminance (argument \code{l}) for generated HCL colors. Can be a single
#'   value or a length-2 range (default: 65).
#' @param power Power transformation for generated HCL colors (default: 1.2).
#' @param shuffle Logical; if \code{TRUE}, shuffle colors within each palette chunk
#'   (default: \code{FALSE}).
#' @param seed Seed used when \code{shuffle = TRUE} (default: 42).
#'
#' @details
#' The palette is returned in a deterministic name order:
#' \enumerate{
#'   \item \code{"unique low abundant"}
#'   \item \code{"unknown"}
#'   \item dynamically colored keys (all remaining \code{tax_color} values excluding fixed and special)
#'   \item \code{special_taxa} keys
#'   \item \code{"shared low abundant"}
#' }
#' Only keys that are present in \code{names(col_tax)} are kept in the final output.
#'
#' @return A named character vector of colors. Names are palette keys (values of \code{data$tax_color}).
#'
#' #' @seealso \code{\link{classify_taxa_patterns}}, \code{\link{plot_alluvial}},
#'   \code{\link{generate_palette_hcl}}
#' 
#' @examples
#' \dontrun{
#' # Minimal example using tax_color keys
#' df <- tibble::tibble(
#'   tax_color = c("unknown","unique low abundant","shared low abundant",
#' "Flavobacteriaceae","Bacillaceae"),
#'   RA = c(0.1, 0.05, 0.2, 0.4, 0.25)
#' )
#'
#' pal <- generate_alluvial_palette(df, tax_level = "Family")
#' pal
#'
#' # With special taxa that keep distinct colors
#' pal2 <- generate_alluvial_palette(
#'   df,
#'   tax_level = "Family",
#'   special_taxa = c("Flavobacteriaceae")
#' )
#' pal2
#'}
#' @export
generate_alluvial_palette <- function(data,
                                      tax_level,
                                      group_by_higher_tax = NULL,     
                                      group_palette_map = NULL,
                                      group_default_side = "both",
                                      special_taxa = character(0),
                                      special_side = c("both", "left", "right", "alphabetical"),
                                      palette_list = NULL,
                                      custom_palette = NULL,
                                      fixed_colors = c(
                                        "unique low abundant" = "#A5A5A5",
                                        "unknown" = "#000000",
                                        "shared low abundant" = "#E5E5E5"
                                      ),
                                      cmax = 100,
                                      luminance = 65,
                                      power = 1.2,
                                      shuffle = FALSE,
                                      seed = 42) {
  
  special_side <- match.arg(special_side)
  
  if (!is.null(group_by_higher_tax) && is.null(group_palette_map)) {
    stop("group_palette_map must be provided when using group_by_higher_tax")
  }
  if (!is.null(group_by_higher_tax) && !group_by_higher_tax %in% colnames(data)) {
    stop("group_by_higher_tax '", group_by_higher_tax, "' not found in data")
  }

  # Return user-defined palette if provided
  if (!is.null(custom_palette)) {
    missing <- setdiff(unique(data$tax_color), names(custom_palette))
    if (length(missing) > 0) stop("Custom palette missing colors for: ", paste(missing, collapse = ", "))
    return(custom_palette)
  }
  
  # Default palettes
  if (is.null(palette_list) || length(palette_list) == 0) {
    palette_list <- c("RedYel", "Emrld", "Red-Blue", "YlOrRd", "Purp", "Peach", "Blues")
  }
  
  # Generating grouped palette
  if (!is.null(group_by_higher_tax) && !is.null(group_palette_map)) {
    
    all_taxa <- unique(data$tax_color)
    dynamic_taxa <- setdiff(all_taxa, c(names(fixed_colors), special_taxa))
    
    # Build grouping lookup from data
    data_for_grouping <- data %>%
      dplyr::filter(tax_color %in% dynamic_taxa) %>%
      dplyr::filter(!(!!rlang::sym(group_by_higher_tax) %in% names(fixed_colors))) %>%
      dplyr::select(!!rlang::sym(group_by_higher_tax), tax_color) %>%
      dplyr::distinct()
    
    # Check all higher taxa have palette assignments
    higher_taxa_in_data <- unique(data_for_grouping[[group_by_higher_tax]])
    missing_assignments <- setdiff(higher_taxa_in_data, names(group_palette_map))
    if (length(missing_assignments) > 0) {
      stop("Missing palette assignments for: ", 
           paste(missing_assignments, collapse = ", "),
           ". Add them to group_palette_map.")
    }
    
    # Generate grouped colors
    dynamic_colors <- generate_grouped_palette(
      data        = data_for_grouping,
      group_col   = group_by_higher_tax,
      item_col    = "tax_color",
      palette_map = group_palette_map,
      default_side = group_default_side,
      cmax        = cmax,
      luminance   = luminance,
      power       = power,
      shuffle     = shuffle,
      seed        = seed
    )
    
    # Special taxa colors
    if (length(special_taxa) > 0) {
      special_taxa_in_data <- intersect(special_taxa, all_taxa)
      if (special_side == "alphabetical") special_taxa_in_data <- sort(special_taxa_in_data)
      special_colors <- colorspace::qualitative_hcl(
        n = length(special_taxa_in_data), palette = "Dark 3"
      )
      names(special_colors) <- special_taxa_in_data
    } else {
      special_colors <- character(0)
    }
    
    # Order and return
    sorted_tax_names <- c(
      "unique low abundant",
      "unknown",
      names(dynamic_colors),
      names(special_colors),
      "shared low abundant"
    )
    col_tax <- c(dynamic_colors, special_colors, fixed_colors)
    col_tax <- col_tax[intersect(sorted_tax_names, names(col_tax))]
    return(col_tax)
  }

  # Identify taxa categories
  all_taxa <- unique(data$tax_color)
  dynamic_taxa <- setdiff(all_taxa, names(fixed_colors))
  dynamic_taxa <- setdiff(dynamic_taxa, special_taxa)
  
  # Generate dynamic colors
  n_dyn <- length(dynamic_taxa)
  n_pal <- length(palette_list)
  base_n <- floor(n_dyn / n_pal)
  extra <- n_dyn %% n_pal
  if (shuffle) set.seed(seed)
  dynamic_colors <- unlist(lapply(seq_along(palette_list), function(i) {
    n_i <- base_n + ifelse(i <= extra, 1, 0)
    if (n_i > 0) {
      cols <- colorspace::sequential_hcl(palette = palette_list[i], n = n_i,
                                         cmax = cmax, l = luminance, power = power)
      if (shuffle) sample(cols) else cols
    } else NULL
  }))
  names(dynamic_colors) <- dynamic_taxa
  
  # Special taxa ordering
  if (length(special_taxa) > 0) {
    special_taxa_in_data <- intersect(special_taxa, all_taxa)
    if (special_side == "alphabetical") special_taxa_in_data <- sort(special_taxa_in_data)
    special_colors <- colorspace::qualitative_hcl(n = length(special_taxa_in_data), palette = "Dark 3")
    names(special_colors) <- special_taxa_in_data
  } else {
    special_colors <- character(0)
  }
  
  # Final stacking order: bottom → top visually
  # Base: shared low abundant, then dynamic, then special, then unique low, top unknown
  sorted_tax_names <- c(
    "unique low abundant",
    "unknown",
    names(dynamic_colors),
    names(special_colors),
    "shared low abundant"
  )
  
  # Combine all colors
  col_tax <- c(dynamic_colors, special_colors, fixed_colors)
  
  # Ensure final order matches stacking
  col_tax <- col_tax[intersect(sorted_tax_names, names(col_tax))]
  
  return(col_tax)
}