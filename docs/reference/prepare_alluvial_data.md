# Aggregate microbiome data for alluvial plot visualization

Aggregates ASV-level data with pre-calculated relative abundance (`RA`)
to a chosen taxonomic level and computes mean group compositions
(sample-weighted). Optionally completes missing taxa with zeros to
ensure a rectangular structure required by alluvial/Sankey plots.

## Usage

``` r
prepare_alluvial_data(
  data,
  tax_level,
  group_col = "Group",
  SampleID_col = "SampleID",
  groups = NULL,
  complete_zero = TRUE,
  clean_taxonomy = TRUE,
  preserve_higher_taxonomy = FALSE,
  hierarchy = c("Domain", "Phylum", "Class", "Order", "Family", "Genus"),
  clean_levels = c("Phylum", "Class", "Order", "Family", "Genus")
)
```

## Arguments

- data:

  Data frame with `RA` already calculated at ASV level.

- tax_level:

  Character string specifying the taxonomic level to aggregate to (e.g.
  `"Family"`).

- group_col:

  Character string specifying the grouping column (e.g. `"SampleType"`).

- SampleID_col:

  Character string specifying the sample ID column (default:
  `"SampleID"`).

- groups:

  Character vector of groups to include (NULL = all groups present in
  `group_col`).

- complete_zero:

  Logical; if TRUE, completes missing taxa with `RA = 0` within each
  sample before summarizing (default: TRUE).

- clean_taxonomy:

  Logical; if TRUE, clean taxonomy using `replace_incertae_sedis_NAs`
  (default: TRUE).

- preserve_higher_taxonomy:

  Logical; if TRUE, keep higher taxonomic ranks up to `tax_level`
  (default: FALSE).

- hierarchy:

  Taxonomic hierarchy for cleaning (default:
  `c("Domain","Phylum","Class","Order","Family","Genus")`).

- clean_levels:

  Taxonomic levels to clean (default:
  `c("Phylum","Class","Order","Family","Genus")`).

## Value

A data frame with mean relative abundance per `group_col` and
`tax_level`.

- RA:

  Mean relative abundance per group after normalization (0–1).

- `<group_col>`:

  Grouping variable used on the alluvial x-axis (the column named by
  `group_col`).

- `<tax_level>`:

  Taxon identifier at the aggregation level (the column named by
  `tax_level`).

- Higher taxonomy columns:

  Only when `preserve_higher_taxonomy = TRUE`: higher ranks up to
  `tax_level`.

## `complete_zero` behavior

- `TRUE`:

  Completes missing taxa within each sample with `RA = 0`, so taxa are
  averaged over all samples.

- `FALSE`:

  Taxa absent from some samples are averaged only over samples where
  they appear (can inflate sporadic taxa).

## `preserve_higher_taxonomy` behavior

- `TRUE`:

  Keeps higher ranks up to `tax_level`. Note: rows created by the final
  completion may have missing higher-rank labels unless those columns
  are also completed/filled.

- `FALSE`:

  Only returns `tax_level` plus grouping columns.

## Data processing in the function step-by-step

1.  Optionally clean taxonomy using `replace_incertae_sedis_NAs`.

2.  Aggregate ASV-level `RA` to `tax_level` within each sample using
    [`sum()`](https://rdrr.io/r/base/sum.html).

3.  If `complete_zero = TRUE`, add missing taxa per sample with `RA = 0`
    (taxa are taken from all values observed in `data[[tax_level]]`).

4.  Compute the mean `RA` per `group_col` (mean across samples), then
    normalize within each group so `sum(RA) = 1`.

5.  Ensure a complete `group_col` \\\times\\ `tax_level` grid (missing
    combinations are set to `RA = 0`).

## Examples

``` r
library(dplyr)
library(tidyr)

# Toy example showing how complete_zero changes group means:
# Taxon B is absent from sample S1 (no row), but present in S2.
toy <- tibble::tibble(
  SampleID = c("S1","S2","S2"),
  Group    = c("G1","G1","G1"),
  Class    = c("A","A","B"),
  RA       = c(1.0, 0.5, 0.5)
)

# Without completion: B is averaged only over samples where it appears (inflated)
out_no0 <- prepare_alluvial_data(
  data = toy,
  tax_level = "Class",
  group_col = "Group",
  SampleID_col = "SampleID",
  complete_zero = FALSE,
  clean_taxonomy = FALSE
)

# With completion: missing taxa count as 0 in those samples
out_0 <- prepare_alluvial_data(
  data = toy,
  tax_level = "Class",
  group_col = "Group",
  SampleID_col = "SampleID",
  complete_zero = TRUE,
  clean_taxonomy = FALSE
)

out_no0 %>% arrange(Class)
#> # A tibble: 2 × 3
#>   Group Class    RA
#>   <chr> <chr> <dbl>
#> 1 G1    A       0.6
#> 2 G1    B       0.4
out_0   %>% arrange(Class)
#> # A tibble: 2 × 3
#>   Group Class    RA
#>   <chr> <chr> <dbl>
#> 1 G1    A      0.75
#> 2 G1    B      0.25
```
