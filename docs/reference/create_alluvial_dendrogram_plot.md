# Create an aligned alluvial + dendrogram figure (full workflow wrapper)

High-level convenience wrapper that builds a composite figure combining:

1.  A grouped ASV matrix (mean abundance per `group_col`),

2.  A hierarchical clustering dendrogram built on that grouped matrix,

3.  An alluvial (Sankey) plot of taxonomic composition across the same
    groups,

4.  A final aligned composite figure with optional external legend.

## Usage

``` r
create_alluvial_dendrogram_plot(
  asv_matrix,
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
  dend_limits_bottom = NULL
)
```

## Arguments

- asv_matrix:

  Numeric matrix with samples as columns and taxa/ASVs as rows. Column
  names must match `metadata[[sample_col]]`.

- metadata:

  Data frame with one row per sample. Must contain `sample_col` and
  `group_col`.

- sample_col:

  Character string specifying the sample ID column in `metadata`
  (default: `"SampleID"`).

- group_col:

  Character string specifying the grouping variable used both for:

  - averaging samples in the ASV matrix (clustering),

  - defining the x-axis of the alluvial plot.

- alluvial_data:

  Long-format data used for the alluvial workflow (typically ASV-level
  relative abundance data).

- tax_level:

  Taxonomic column used for aggregation in the alluvial plot (e.g.
  `"Family"`, `"Class"`).

- group_order:

  Determines group ordering in both dendrogram and alluvial plot. One
  of:

  - `"metadata"` – first-appearance order in metadata,

  - `"alphabetical"` – sorted group names,

  - a character vector of explicit group names.

- low_abundance_threshold:

  Numeric threshold passed to
  [`classify_taxa_patterns`](classify_taxa_patterns.md) via
  [`create_alluvial_plot`](create_alluvial_plot.md) (default: 0.01).

- palette_list:

  Character vector of palette names passed to
  [`generate_alluvial_palette`](generate_alluvial_palette.md) (default:
  `NULL`).

- custom_palette:

  Optional named color vector. If provided, automatic palette generation
  is skipped.

- dend_color_by:

  Column in `metadata` used to color dendrogram tips. Defaults to
  `group_col`.

- dend_color_palette:

  Named vector mapping levels of `dend_color_by` to colors.

- dend_shape_by:

  Optional column in `metadata` used to set dendrogram tip shapes
  (default: `NULL`).

- dend_shape_palette:

  Optional vector passed to
  [`ggplot2::scale_shape_manual()`](https://ggplot2.tidyverse.org/reference/scale_manual.html)
  (default: `NULL`).

- theme_alluvial:

  ggplot2 theme added to the alluvial plot (default:
  [`ggplot2::theme_minimal()`](https://ggplot2.tidyverse.org/reference/ggtheme.html)).

- theme_dendrogram:

  ggplot2 theme added to the dendrogram plot (default:
  [`ggplot2::theme_void()`](https://ggplot2.tidyverse.org/reference/ggtheme.html)).

- build_dendrogram_args:

  Named list forwarded to [`build_dendrogram`](build_dendrogram.md)
  (e.g. `distance_method`, `cluster_method`).

- plot_dendrogram_args:

  Named list forwarded to [`plot_dendrogram`](plot_dendrogram.md) (e.g.
  `orientation`, `point_size`, `guides_obj`).

- alluvial_args:

  Named list forwarded to
  [`create_alluvial_plot`](create_alluvial_plot.md). May include routed
  argument lists (`prepare_args`, `classify_args`, `palette_args`,
  `plot_args`) and `return_all`.

- combine_args:

  Named list forwarded to
  [`combine_dendrogram_alluvial`](combine_dendrogram_alluvial.md) (e.g.
  `dend_height`, `legend`, `legend_position`).

- post_plot_guides:

  Modify legend guides of the plot (default: NULL).

- dend_limits_left, dend_limits_right:

  “sweet spot” for dendrogram x view: dend xlim = c(dend_limits_left,
  n + dend_limits_right). Use e.g. -0.75 and 0.2.

- dend_limits_bottom, dend_limits_top:

  “sweet spot” for dendrogram y view: dend ylim = c(dend_limits_bottom,
  n + dend_limits_top).

  \#' @section Workflow overview: Internally, the function performs:

  1.  [`create_grouped_matrix`](create_grouped_matrix.md) – average ASV
      abundances within `group_col`.

  2.  [`build_dendrogram`](build_dendrogram.md) – compute distance +
      hierarchical clustering.

  3.  [`plot_dendrogram`](plot_dendrogram.md) – render dendrogram as a
      ggplot.

  4.  [`create_alluvial_plot`](create_alluvial_plot.md) – run the
      complete alluvial pipeline (prepare → classify → palette → plot).

  5.  [`combine_dendrogram_alluvial`](combine_dendrogram_alluvial.md) –
      align and combine both panels.

## Value

A named list containing:

- grouped_matrix:

  Matrix used for clustering (rows = taxa, columns = groups).

- group_levels:

  Final group ordering applied to both panels.

- dendrogram:

  List returned by [`build_dendrogram`](build_dendrogram.md).

- dendrogram_plot:

  Dendrogram rendered as a ggplot.

- alluvial:

  Output of [`create_alluvial_plot`](create_alluvial_plot.md) (list if
  `return_all=TRUE`).

- combined_plot:

  Final composite grob suitable for `ggsave()`.

## Details

This function is designed for reproducible figures where clustering
structure (top panel) and compositional structure (bottom panel) share
identical group ordering and horizontal alignment.

This wrapper guarantees that:

- The dendrogram is computed on the same group structure used for the
  alluvial x-axis.

- Group ordering is identical across clustering and composition panels.

- Horizontal panel alignment is preserved in the final composite.

## Argument ownership and routing

This wrapper **owns** the structural arguments passed to
[`create_alluvial_plot`](create_alluvial_plot.md): `data`, `tax_level`,
`group_col`, `SampleID_col`, `groups`, `low_abundance_threshold`,
`palette_list`, and `custom_palette`.

Any duplicated arguments provided inside `alluvial_args` (including
nested `prepare_args`, `classify_args`, `palette_args`, or `plot_args`)
are automatically removed to prevent argument-matching conflicts.

## Examples

``` r
# See the package vignette for a full reproducible example including
# ASV matrix, metadata, and long-format alluvial input.
```
