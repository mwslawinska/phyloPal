# Normalize Incertae sedis and fill missing taxonomy ranks

Cleans hierarchical taxonomy columns by standardizing common "Incertae
sedis" variants, filling missing child ranks from stable parent ranks,
and replacing empty/uninformative entries with `"unknown"`.

## Usage

``` r
replace_incertae_sedis_NAs(
  df,
  hierarchy = c("Domain", "Phylum", "Class", "Order", "Family", "Genus"),
  clean_levels = c("Phylum", "Class", "Order", "Family", "Genus")
)
```

## Arguments

- df:

  A data.frame or tibble containing taxonomy columns.

- hierarchy:

  Character vector of column names in rank order from highest (parent)
  to lowest (child). Only columns present in `df` are used.

- clean_levels:

  Character vector of ranks to normalize and fill. Only ranks present in
  both `df` and `hierarchy` are used.

## Value

An object of the same class as `df` (tibble in, tibble out) with cleaned
taxonomy columns.

## Details

The function is designed for rank-ordered taxonomy (e.g., Domain →
Phylum → Class → ...), and performs a downward pass through the
hierarchy to propagate parent information.

The following rules are applied:

- "Incertae sedis" variants (case-insensitive; underscores/spaces
  tolerated) are converted to `"inc. sed."`.

- Empty strings (`""`) in `clean_levels` are treated as missing (`NA`).

- A parent is considered *unstable* if it is missing (`NA`) or matches
  `"unknown"`, `"unclassified"`, or `"inc. sed."` (case-insensitive).

- If a child is `"inc. sed."` and the parent is stable, the child
  becomes `"<parent> inc. sed."`; otherwise it becomes `"unknown"`.

- If a child is missing and the parent is stable, the child becomes
  `"<parent>, unclassified"`; otherwise it becomes `"unknown"`.

- Underscores are converted to spaces and values are trimmed.

## Examples

``` r
tax_df <- data.frame(
  Domain = c("Bacteria", "Bacteria", NA),
  Phylum = c("Proteobacteria", "IncertaE_Sedis", "Actinobacteriota"),
  Class  = c(NA, "Alphaproteobacteria", ""),
  Order  = c("", "Rhizobiales", "inc. sed."),
  Family = c("Rhizobiaceae", NA, "unknown"),
  Genus  = c("", "inc. sed.", NA)
)

replace_incertae_sedis_NAs(
  df = tax_df,
  hierarchy = c("Domain", "Phylum", "Class", "Order", "Family", "Genus"),
  clean_levels = c("Phylum", "Class", "Order", "Family", "Genus")
)
#>     Domain             Phylum                          Class       Order
#> 1 Bacteria     Proteobacteria   Proteobacteria, unclassified     unknown
#> 2 Bacteria Bacteria inc. sed.            Alphaproteobacteria Rhizobiales
#> 3     <NA>   Actinobacteriota Actinobacteriota, unclassified     unknown
#>                      Family                      Genus
#> 1              Rhizobiaceae Rhizobiaceae, unclassified
#> 2 Rhizobiales, unclassified                    unknown
#> 3                   unknown                    unknown
```
