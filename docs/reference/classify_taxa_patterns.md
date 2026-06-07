# Classify taxa by abundance patterns across groups

Classifies taxa as shared/unique and abundant/low abundant based on
their relative abundances across multiple groups. Optionally detects
“mixed abundance” taxa that are abundant in some groups but low abundant
in others.

## Usage

``` r
classify_taxa_patterns(
  data,
  tax_level,
  group_col = "Group",
  groups = NULL,
  low_abundance_threshold = 0.01,
  special_taxa = character(0),
  detect_mixed_abundance = TRUE,
  unknown_label = "unknown"
)
```

## Arguments

- data:

  A data frame containing at least `tax_level`, `group_col`, and `RA`
  (typically the output of
  [`prepare_alluvial_data()`](prepare_alluvial_data.md)).

- tax_level:

  Character string specifying the taxonomic level (e.g. `"Family"`,
  `"Class"`).

- group_col:

  Character string specifying the grouping column (default: `"Group"`).

- groups:

  Character vector of groups to compare. If NULL, uses all unique values
  in `data[[group_col]]`.

- low_abundance_threshold:

  Numeric threshold below which taxa are considered low abundant
  (default: 0.01).

- special_taxa:

  Character vector of taxa that should never be marked as low abundant
  (default: `character(0)`).

- detect_mixed_abundance:

  Logical; if TRUE, detect taxa that are abundant in some groups and low
  abundant in others and label as `"shared mixed abundance"` (default:
  TRUE). These taxa keep their own palette key (`tax_color = <taxon>`)
  rather than being mapped to a low-abundance bin.

- unknown_label:

  Label used for unknown taxa (default: `"unknown"`).

## Value

The input data with additional classification columns.

- tax_type:

  Primary abundance pattern classification (e.g. shared/unique;
  abundant/low; and optionally mixed).

- category:

  Final category used for summaries and plotting (may match `tax_type`
  depending on workflow).

- tax_color:

  A label used for color mapping (e.g. low-abundant bins may share a
  common color label).

&nbsp;

- tax_val:

  Character copy of `<tax_level>` used internally for classification.

Other columns (including higher taxonomy ranks) are preserved if present
in `data`.

## Details

This function does not set plotting order (legend or stratum stacking).
In the palette-driven workflow, ordering is handled in
[`plot_alluvial()`](plot_alluvial.md) via the order of
`names(custom_palette)`.

## Examples

``` r
library(dplyr)
#> 
#> Attaching package: ‘dplyr’
#> The following objects are masked from ‘package:stats’:
#> 
#>     filter, lag
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, setdiff, setequal, union

# Two groups, three taxa:
# - A is abundant in both groups (shared abundant)
# - B is abundant in G1 but low in G2 (mixed abundance if enabled)
# - C appears only in G2 (unique)
toy <- tibble::tibble(
  Group = c("G1", "G1", "G1", "G2", "G2", "G2"),
  Class = c("A", "B", "unknown", "A", "B", "C"),
  RA    = c(0.80, 0.15, 0.05, 0.80, 0.005, 0.195)
)

out <- classify_taxa_patterns(
  data = toy,
  tax_level = "Class",
  group_col = "Group",
  low_abundance_threshold = 0.01,
  detect_mixed_abundance = TRUE
)

out %>% arrange(Group, Class)
#> # A tibble: 6 × 7
#>   Group Class      RA tax_val tax_type               category          tax_color
#>   <chr> <chr>   <dbl> <chr>   <chr>                  <chr>             <chr>    
#> 1 G1    A       0.8   A       shared abundant        shared abundant   A        
#> 2 G1    B       0.15  B       shared mixed abundance shared mixed abu… B        
#> 3 G1    unknown 0.05  unknown unknown                unknown           unknown  
#> 4 G2    A       0.8   A       shared abundant        shared abundant   A        
#> 5 G2    B       0.005 B       shared mixed abundance shared mixed abu… B        
#> 6 G2    C       0.195 C       unique abundant        unique abundant   C        

out_mixed <- classify_taxa_patterns(toy, "Class", "Group", detect_mixed_abundance = TRUE)
out_nomix <- classify_taxa_patterns(toy, "Class", "Group", detect_mixed_abundance = FALSE)
out_mixed %>%
  filter(Class == "B") %>%
  dplyr::select(Group, Class, RA, tax_type, tax_color)
#> # A tibble: 2 × 5
#>   Group Class    RA tax_type               tax_color
#>   <chr> <chr> <dbl> <chr>                  <chr>    
#> 1 G1    B     0.15  shared mixed abundance B        
#> 2 G2    B     0.005 shared mixed abundance B        
out_nomix %>%
  filter(Class == "B") %>%
  dplyr::select(Group, Class, RA, tax_type, tax_color)
#> # A tibble: 2 × 5
#>   Group Class    RA tax_type            tax_color          
#>   <chr> <chr> <dbl> <chr>               <chr>              
#> 1 G1    B     0.15  unique abundant     B                  
#> 2 G2    B     0.005 unique low abundant unique low abundant
```
