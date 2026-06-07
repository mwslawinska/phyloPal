# Generate HCL color palette for taxonomic visualization

Generates a named HCL-based color palette for taxa at a specified
taxonomic level. Optionally supports hierarchical grouping, where colors
are assigned within higher taxonomic groups (e.g. coloring families
grouped by phylum).

## Usage

``` r
generate_palette_hcl(
  data,
  tax_level,
  group_by_higher_tax = NULL,
  group_palette_map = NULL,
  group_default_side = "both",
  order_by_higher_tax = FALSE,
  order_groups = c("alphabetical", "data"),
  order_within_groups = c("alphabetical", "data"),
  fixed_colors_enabled = TRUE,
  fixed_colors = c(unknown = "#000000", `low abundant` = "#E5E5E5"),
  fixed_colors_position = "end",
  palette_list = NULL,
  cmax = 100,
  luminance = 65,
  power = 1.2,
  shuffle = TRUE,
  seed = 42
)
```

## Arguments

- data:

  A data frame containing at least `tax_level`. If `group_by_higher_tax`
  is used, the corresponding higher-taxonomy column must also be
  present.

- tax_level:

  Character string specifying the taxonomic level to generate colors for
  (e.g. `"Family"`, `"Genus"`).

- group_by_higher_tax:

  Optional higher taxonomic level used to group taxa before color
  assignment (e.g. `"Phylum"`, `"Class"`). Recommended only for small
  datasets (e.g. synthetic communities with fewer than ~10 higher-level
  taxa). For complex natural communities, the ungrouped palette is
  typically more interpretable.

- group_palette_map:

  Optional named list mapping higher taxa to palette specifications
  (required if `group_by_higher_tax` is used). Each element can be
  either:

  - a single palette name (character), or

  - a list with elements `palette` and `side` for diverging palettes.

- group_default_side:

  Default side for grouped diverging palettes: `"both"`, `"left"`, or
  `"right"` (default: `"both"`).

- order_by_higher_tax:

  Logical; if `TRUE`, taxa are ordered deterministically by higher
  taxonomy and then by `tax_level` (default: `TRUE` when
  `group_by_higher_tax` is used).

- order_groups:

  How to order higher-taxonomy groups when `order_by_higher_tax = TRUE`.
  Options: `"alphabetical"` (default) or `"data"`.

- order_within_groups:

  How to order taxa within each higher-taxonomy group when
  `order_by_higher_tax = TRUE`. Options: `"alphabetical"` (default) or
  `"data"`.

- fixed_colors_enabled:

  Logical; if `TRUE`, apply fixed colors to specific taxa (default:
  `TRUE`).

- fixed_colors:

  Named character vector mapping specific taxa to fixed colors (default:
  `c("unknown" = "#000000", "low abundant" = "#E5E5E5")`).

- fixed_colors_position:

  Position of fixed-color taxa in the resulting palette: `"end"`
  (default), `"beginning"`, or `"alphabetical"`.

- palette_list:

  Character vector of HCL palette names used for dynamically generated
  colors. If `NULL`, a predefined set is used.

- cmax:

  Maximum chroma value for generated HCL colors (default: 100).

- luminance:

  Luminance value (single number) or range (length-2 numeric) used in
  HCL generation (default: 65).

- power:

  Power transformation applied in HCL color interpolation (default:
  1.2).

- shuffle:

  Logical; if `TRUE`, dynamically generated colors are shuffled to
  reduce adjacency similarity (default: `TRUE`).

- seed:

  Numeric seed used when `shuffle = TRUE` (default: 42).

## Value

A named character vector where names correspond to taxa at `tax_level`
and values are hexadecimal color codes.

## Details

**Color generation workflow:**

1.  Identify unique taxa at `tax_level`.

2.  Optionally group taxa by `group_by_higher_tax`.

3.  Generate HCL-based colors either globally (ungrouped) or within each
    higher-taxonomy group.

4.  Apply ordering rules (if enabled).

5.  Insert fixed colors according to `fixed_colors_position`.

When hierarchical grouping is enabled, each higher-taxonomy group
receives its own palette (defined in `group_palette_map`). This produces
visually coherent sub-palettes but may reduce distinguishability in
highly diverse datasets.

## See also

[`process_barplot_data`](process_barplot_data.md),
[`plot_taxonomic_barplot`](plot_taxonomic_barplot.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(dplyr)

toy <- tibble::tibble(
  Phylum = c("P1","P1","P2","P2"),
  Family = c("F1","F2","F3","F4")
)

# Simple ungrouped palette
pal_simple <- generate_tax_palette(
  data = toy,
  tax_level = "Family"
)

# Hierarchical grouping by Phylum
pal_grouped <- generate_tax_palette(
  data = toy,
  tax_level = "Family",
  group_by_higher_tax = "Phylum",
  group_palette_map = list(
    P1 = "Blues",
    P2 = "Reds"
  )
)

pal_simple
pal_grouped
} # }
```
