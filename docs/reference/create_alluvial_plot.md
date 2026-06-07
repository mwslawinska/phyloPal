# Create an alluvial plot via the full workflow wrapper

Runs the alluvial workflow end-to-end:

1.  [`prepare_alluvial_data`](prepare_alluvial_data.md): aggregate
    ASV-level `RA` to `tax_level` and compute group means,

2.  [`classify_taxa_patterns`](classify_taxa_patterns.md): assign
    abundance-pattern classes and create `tax_color`,

3.  [`generate_alluvial_palette`](generate_alluvial_palette.md):
    generate a palette for `tax_color` (skipped if `custom_palette` is
    provided),

4.  [`plot_alluvial`](plot_alluvial.md): draw the plot.

## Usage

``` r
create_alluvial_plot(
  data,
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
  return_all = FALSE
)
```

## Arguments

- data:

  Raw long-format data with ASV-level `RA` and taxonomy columns.

- tax_level:

  Character string naming the taxonomy rank to plot (e.g. `"Class"`,
  `"Family"`).

- group_col:

  Character string naming the grouping column used on the x-axis.

- SampleID_col:

  Character string naming the sample ID column (default: `"SampleID"`).

- groups:

  Optional character vector restricting groups included. If `NULL`, all
  groups are used.

- low_abundance_threshold:

  Numeric threshold passed to
  [`classify_taxa_patterns`](classify_taxa_patterns.md) (default: 0.01).

- palette_list:

  Optional vector of HCL palette names passed to
  [`generate_alluvial_palette`](generate_alluvial_palette.md).

- custom_palette:

  Optional named palette mapping `tax_color` keys to colors. If
  provided, palette generation is skipped.

- prepare_args:

  Named list of extra arguments passed to
  [`prepare_alluvial_data`](prepare_alluvial_data.md).

- classify_args:

  Named list of extra arguments passed to
  [`classify_taxa_patterns`](classify_taxa_patterns.md).

- palette_args:

  Named list of extra arguments passed to
  [`generate_alluvial_palette`](generate_alluvial_palette.md) (only used
  if `custom_palette` is `NULL`).

- plot_args:

  Named list of extra arguments passed to
  [`plot_alluvial`](plot_alluvial.md).

- return_all:

  Logical; if `TRUE`, returns a list with intermediate objects:
  `data_prepared`, `data_classified`, `palette`, `plot`.

## Value

A `ggplot` object, or a named list if `return_all = TRUE`.

## Details

Additional arguments are routed to the underlying functions via the
`*_args` lists.

## Examples

``` r
if (FALSE) { # \dontrun{
library(dplyr)

toy <- tibble::tibble(
  SampleID = c("S1","S1","S2","S2"),
  Group    = c("G1","G1","G1","G1"),
  Class    = c("A","B","A","C"),
  RA       = c(0.7, 0.3, 0.6, 0.4)
)

res <- create_alluvial_plot(
  data = toy,
  tax_level = "Class",
  group_col = "Group",
  prepare_args = list(complete_zero = TRUE, clean_taxonomy = FALSE),
  classify_args = list(detect_mixed_abundance = TRUE),
  plot_args = list(alpha = 0.7),
  return_all = TRUE
)

res$plot
} # }
```
