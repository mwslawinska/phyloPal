# Generate a color palette for alluvial plots from tax_color keys

Creates a named color palette for alluvial/Sankey plots where the
palette keys are taken from `data$tax_color`. Three keys are treated as
fixed categories by default (`"unique low abundant"`, `"unknown"`,
`"shared low abundant"`), taxa listed in `special_taxa` receive their
own colors, and all remaining keys receive automatically generated HCL
colors.

## Usage

``` r
generate_alluvial_palette(
  data,
  tax_level,
  group_by_higher_tax = NULL,
  group_palette_map = NULL,
  group_default_side = "both",
  special_taxa = character(0),
  special_side = c("both", "left", "right", "alphabetical"),
  palette_list = NULL,
  custom_palette = NULL,
  fixed_colors = c(`unique low abundant` = "#A5A5A5", unknown = "#000000",
    `shared low abundant` = "#E5E5E5"),
  cmax = 100,
  luminance = 65,
  power = 1.2,
  shuffle = FALSE,
  seed = 42
)
```

## Arguments

- data:

  A data frame that must contain a `tax_color` column. The returned
  palette names correspond to the unique values of `data$tax_color`.

- tax_level:

  Character string (currently unused). Included for API compatibility
  with higher-level workflow wrappers.

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

- special_taxa:

  Character vector of palette keys (i.e., values of `tax_color`) that
  should receive their own colors (drawn from `palette = "Dark 3"`) and
  be excluded from the automatically generated pool.

- special_side:

  Controls ordering of `special_taxa`. Currently only `"alphabetical"`
  has an effect (special taxa are sorted); other values are accepted for
  forward compatibility.

- palette_list:

  Character vector of colorspace sequential HCL palette names used to
  generate colors for non-fixed, non-special keys. If `NULL` or empty,
  defaults to
  `c("RedYel","Emrld","Red-Blue","YlOrRd","Purp","Peach","Blues")`.

- custom_palette:

  Optional named character vector of colors. If provided, the function
  validates that it covers all unique values of `data$tax_color` and
  returns it unchanged.

- fixed_colors:

  Named character vector of fixed colors for reserved keys. Names must
  match values used in `data$tax_color`. Defaults to:
  `c("unique low abundant"="#A5A5A5","unknown"="#000000","shared low abundant"="#E5E5E5")`.

- cmax:

  Maximum chroma for generated HCL colors (passed to
  [`colorspace::sequential_hcl()`](https://colorspace.R-Forge.R-project.org/reference/hcl_palettes.html),
  default: 100).

- luminance:

  Luminance (argument `l`) for generated HCL colors. Can be a single
  value or a length-2 range (default: 65).

- power:

  Power transformation for generated HCL colors (default: 1.2).

- shuffle:

  Logical; if `TRUE`, shuffle colors within each palette chunk (default:
  `FALSE`).

- seed:

  Seed used when `shuffle = TRUE` (default: 42).

## Value

A named character vector of colors. Names are palette keys (values of
`data$tax_color`).

\#' @seealso [`classify_taxa_patterns`](classify_taxa_patterns.md),
[`plot_alluvial`](plot_alluvial.md),
[`generate_palette_hcl`](generate_palette_hcl.md)

## Details

The palette is returned in a deterministic name order:

1.  `"unique low abundant"`

2.  `"unknown"`

3.  dynamically colored keys (all remaining `tax_color` values excluding
    fixed and special)

4.  `special_taxa` keys

5.  `"shared low abundant"`

Only keys that are present in `names(col_tax)` are kept in the final
output.

## Examples

``` r
if (FALSE) { # \dontrun{
# Minimal example using tax_color keys
df <- tibble::tibble(
  tax_color = c("unknown","unique low abundant","shared low abundant",
"Flavobacteriaceae","Bacillaceae"),
  RA = c(0.1, 0.05, 0.2, 0.4, 0.25)
)

pal <- generate_alluvial_palette(df, tax_level = "Family")
pal

# With special taxa that keep distinct colors
pal2 <- generate_alluvial_palette(
  df,
  tax_level = "Family",
  special_taxa = c("Flavobacteriaceae")
)
pal2
} # }
```
