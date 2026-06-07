# Create an alluvial (Sankey) plot of taxonomic composition across groups

Draws an alluvial plot where strata represent taxa (or taxon labels)
within each group. By default, strata are drawn for each `tax_val`.
Optionally, strata can be collapsed into one bucket per group (e.g. for
`"unique low abundant"` taxa) to improve readability and avoid numerical
issues when many extremely small strata are present.

## Usage

``` r
plot_alluvial(
  data,
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
  bucket_has_flows = FALSE
)
```

## Arguments

- data:

  A data frame (typically output of
  [`classify_taxa_patterns`](classify_taxa_patterns.md)) containing:

  - `RA`: relative abundance per `group_col` (0–1; typically sums to 1
    within group),

  - `tax_val`: displayed stratum label (often the taxon name),

  - `tax_color`: palette key used for color mapping (may equal `tax_val`
    for most taxa),

  - the grouping column named by `group_col`.

- custom_palette:

  Named character vector of colors. Names must match palette keys in
  `data$tax_color`. The order of `names(custom_palette)` defines legend
  order and (if `order_strata_by_palette=TRUE`) influences stratum
  stacking order.

- tax_level:

  Character string used for the legend title (e.g. `"Family"`). Does not
  affect aggregation.

- group_col:

  Character string naming the grouping column used on the x-axis.

- alpha:

  Numeric transparency for flows (default: 0.7).

- stratum_width:

  Numeric width of strata passed to
  [`ggalluvial::geom_stratum()`](http://corybrunson.github.io/ggalluvial/reference/geom_stratum.md)
  (default: 1/3).

- flow_width:

  Numeric width of flows passed to
  [`ggalluvial::geom_flow()`](http://corybrunson.github.io/ggalluvial/reference/geom_flow.md)
  (default: 1/3).

- add_labels:

  Logical; if `TRUE`, add stratum labels (default: FALSE).

- label_size:

  Numeric text size for labels (default: 3).

- theme_obj:

  A ggplot2 theme object (default:
  [`ggplot2::theme_minimal()`](https://ggplot2.tidyverse.org/reference/ggtheme.html)).

- line_width:

  Numeric line width for stratum borders (default: 0.2).

- legend_position:

  Legend position (default: `"right"`).

- y_axis_label:

  y-axis label (default: `"Relative Abundance"`).

- x_axis_label:

  x-axis label (default: empty string).

- order_strata_by_palette:

  Logical; if `TRUE` (default), strata are ordered by palette-key order
  (via `tax_color`). This makes stacking and legend consistent across
  plots.

- within_key_sort:

  Character; when `order_strata_by_palette=TRUE`, controls how `tax_val`
  are ordered within each palette key. `"alphabetical"` (default) sorts
  taxa within each key; `"none"` keeps taxa in the order encountered
  within each key after key-ordering is applied.

- collapse_strata:

  Logical; if TRUE, collapse selected palette-key bins into a single
  stratum per group. This is useful when many taxa share one palette key
  (e.g. `"unique low abundant"`) and would otherwise generate many
  microscopic strata. Default: FALSE.

- collapse_keys:

  Character vector of palette keys to collapse when
  `collapse_strata=TRUE`. Default:
  `c("unique low abundant","shared low abundant","unknown")`.

- bucket_small:

  Logical; if TRUE, taxa with `RA < bucket_min_ra` are collapsed into
  one bucket per group. This is independent of classification thresholds
  and is intended as a plotting stabilization/readability option.
  Default: FALSE.

- bucket_min_ra:

  Numeric in (0,1); taxa with `RA < bucket_min_ra` are put into
  `bucket_label`. Default: 0.001 (0.1%).

- bucket_label:

  Character label used for the bucket created by `bucket_small`. Must be
  present in the palette (or will be appended as an extra key and drawn
  with `na.value`). Default: `"very low abundant"`.

- bucket_has_flows:

  Logical; if FALSE (default), the bucket created by `bucket_small` will
  not have flows. If TRUE, flows are drawn for that bucket like any
  other shared stratum.

## Value

A `ggplot` object.

## Ordering behavior (palette-driven workflow)

- **Legend order** follows `names(custom_palette)`.

- If `order_strata_by_palette=TRUE`, **stacking order** is driven by
  palette-key order (`tax_color`). This is useful when multiple taxa
  share the same palette key (e.g. low-abundance bins).

## Flow behavior

Flows are drawn only for taxa (`tax_val`) that occur in more than one
group (i.e. present in \\\>1\\ distinct `group_col` levels). Taxa unique
to a single group are shown as strata only.

## Palette coverage

`custom_palette` should ideally include a color for every value of
`data$tax_color`. Missing palette keys are displayed using `na.value`
(default `"grey50"`).

## Stratum identity vs. color identity

- `tax_val` is the taxon label used for strata by default (one stratum
  per taxon).

- `tax_color` is the palette key controlling fill colors and legend
  order (multiple taxa may share one key).

- If `collapse_strata=TRUE` and/or `bucket_small=TRUE`, the plot uses a
  derived stratum ID so that selected bins (e.g.
  `"unique low abundant"`) or tiny taxa are drawn as a single bucket per
  group.

## See also

[`prepare_alluvial_data`](prepare_alluvial_data.md),
[`classify_taxa_patterns`](classify_taxa_patterns.md),
[`generate_alluvial_palette`](generate_alluvial_palette.md),
[`create_alluvial_plot`](create_alluvial_plot.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(dplyr)
library(ggplot2)

toy <- tibble::tibble(
  Group = rep(c("G1","G2"), each = 3),
  tax_val   = rep(c("A","B","C"), times = 2),
  RA        = c(0.7, 0.2, 0.1,  0.6, 0.02, 0.38),
  tax_color = c("A","B","C", "A","unique low abundant","C")
)

pal <- c(
  "unique low abundant" = "#A5A5A5",
  "A" = "#1b9e77",
  "B" = "#d95f02",
  "C" = "#7570b3"
)

plot_alluvial(toy, pal, tax_level = "Class", group_col = "Group")
} # }
```
