# Process barplot data at specified taxonomic level

Prepares long-format data for barplot creation at a specified taxonomic
level. Handles unknown taxonomy, groups low-abundance taxa, and
aggregates relative abundances (RA). Works with data that already has RA
(RA in the 0-1 range) calculated at ASV level.

## Usage

``` r
process_barplot_data(
  data,
  tax_level,
  group_vars,
  SampleID_col = "SampleID",
  low_abundance_threshold = 0.01,
  low_abundance_basis = c("post_aggregation", "per_sample"),
  keep_ratype = c("collapse", "separate"),
  tax_original_suffix = "_original",
  preserve_higher_taxonomy = FALSE,
  clean_taxonomy = TRUE,
  hierarchy = c("Domain", "Phylum", "Class", "Order", "Family", "Genus"),
  clean_levels = c("Phylum", "Class", "Order", "Family", "Genus"),
  agg_fun = c("sum", "mean"),
  normalize_by = NULL,
  drop_zero = TRUE,
  unknown_label = "unknown",
  low_label = "low abundant"
)
```

## Arguments

- data:

  A data frame with RA (relative abundance) already calculated per-ASV.

- tax_level:

  Character string specifying the taxonomic level to plot (e.g.,
  "Family", "Genus").

- group_vars:

  Character vector of grouping variables (e.g. SampleID, SampleType).

- SampleID_col:

  Character string specifying the sample ID column (default: "SampleID")

- low_abundance_threshold:

  Numeric threshold below which taxa are considered low abundance
  (default: 0.01).

- low_abundance_basis:

  When to identify low-abundance taxa:

  "per_sample"

  :   Mark taxa as low abundant at individual sample level BEFORE
      aggregation

  "post_aggregation"

  :   Mark taxa as low abundant AFTER aggregating across samples
      (default)

- keep_ratype:

  Controls whether low-abundance taxa are collapsed or kept separate:

  "collapse"

  :   Relabel taxa below threshold as "low abundant" and collapse into
      one bin (default)

  "separate"

  :   Create `<tax_level>_original`; keep low/unknown as plot labels but
      do not merge originals

- tax_original_suffix:

  Suffix used to create the original-taxon column when
  `keep_ratype="separate"` (default: "\_original"). Example:
  tax_level="Class" -\> "Class_original".

- preserve_higher_taxonomy:

  Logical, whether to keep higher taxonomic levels (default: FALSE).

- clean_taxonomy:

  Logical, whether to clean taxonomy using replace_incertae_sedis_NAs
  (default: TRUE).

- hierarchy:

  Taxonomic hierarchy for cleaning (default: c("Domain", "Phylum",
  "Class", "Order", "Family", "Genus")).

- clean_levels:

  Levels to clean (default: c("Phylum", "Class", "Order", "Family",
  "Genus")). Uses `replace_incertae_sedis_NAs`

- agg_fun:

  Aggregation function for RA when summarizing groups: `"sum"` or
  `"mean"` (default: "sum").

- normalize_by:

  Character vector (or single string) giving grouping columns within
  which RA should sum to 1. If NULL: normalize within `SampleID_col` if
  present, else within `group_vars` excluding taxonomy columns.

- drop_zero:

  Logical; if TRUE, drop rows with `RA <= 0` at the end (default: TRUE).
  Useful to avoid inflating plotting keys with structural zeros.

- unknown_label:

  Label used for unknown taxonomy (default: "unknown").

- low_label:

  Label used for low-abundance bin (default: "low abundant").

## Value

A data frame with processed relative abundances at specified taxonomic
level.

- RA:

  Relative abundance after aggregation and normalization (0–1).

- `<tax_level>`:

  Taxon labels used for plotting (may include
  `low_label`/`unknown_label`).

- `<tax_level>_original`:

  Only when `keep_ratype="separate"`: original taxon identity.

- `<Higher>_true`:

  Only when `preserve_higher_taxonomy=TRUE` and
  `keep_ratype="separate"`: true higher ranks.

## Low abundance handling (controlled by `keep_ratype`)

- collapse:

  Relabel taxa below the threshold as "low abundant" and collapse into a
  single bin.

- separate:

  Keep each original taxon in `<tax_level>_original`, while
  `<tax_level>` is replaced by `low_label`/`unknown_label` for flagged
  taxa.

## Higher taxonomy handling (when `preserve_higher_taxonomy = TRUE`)

- collapse:

  Higher taxonomy columns are also relabeled to "low
  abundant"/"unknown".

- separate:

  True higher taxonomy is re-attached via `<tax_level>_original` into
  `<Higher>_true` columns (e.g. `Phylum_true`).

## Examples

``` r
library(dplyr)

# Minimal toy dataset (already aggregated at ASV level with RA in 0–1 range)
toy <- tibble::tibble(
  SampleID = c("S1","S1","S2","S2"),
  Group    = c("A","A","A","A"),
  Phylum   = c("P1","P2","P1","P2"),
  Class    = c("C1","C2","C1","C2"),
  RA       = c(0.9, 0.1, 0.6, 0.4)
)

# Collapse low-abundance taxa into a single bin
out_collapse <- process_barplot_data(
  data = toy,
  tax_level = "Class",
  group_vars = "SampleID",
  low_abundance_threshold = 0.2,
  keep_ratype = "collapse",
  clean_taxonomy = FALSE
)

# Keep low-abundance taxa separate but flagged
out_separate <- process_barplot_data(
  data = toy,
  tax_level = "Class",
  group_vars = "SampleID",
  low_abundance_threshold = 0.2,
  keep_ratype = "separate",
  clean_taxonomy = FALSE
)

# In separate mode, an additional column "Class_original" appears
names(out_separate)
#> [1] "SampleID"       "Class"          "Class_original" "RA"            
```
