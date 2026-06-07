# Extract a legend grob from a ggplot

Convenience helper around
[`cowplot::get_legend()`](https://wilkelab.org/cowplot/reference/get_legend.html).

## Usage

``` r
extract_legend_grob(p)
```

## Arguments

- p:

  A `ggplot` object.

## Value

A grob containing the legend, or `NULL` if the plot has no legend.
