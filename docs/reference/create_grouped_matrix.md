# Create grouped ASV matrix by averaging within groups

Collapses a sample-by-taxon matrix into a group-by-taxon matrix by
taking the mean abundance across samples within each group. This is
useful for building dendrograms on group centroids (e.g. habitat-level
mean communities).

## Usage

``` r
create_grouped_matrix(
  asv_matrix,
  metadata,
  sample_col = "SampleID",
  group_col,
  group_order = c("alphabetical", "metadata", "custom"),
  group_levels = NULL,
  drop_unmapped_samples = FALSE
)
```

## Arguments

- asv_matrix:

  Numeric matrix with samples as columns and taxa/ASVs as rows. Column
  names must match `metadata[[sample_col]]`.

- metadata:

  Data frame with one row per sample (or at least one row per sample
  ID).

- sample_col:

  Column name in `metadata` containing sample IDs (default:
  `"SampleID"`).

- group_col:

  Column name in `metadata` defining groups to average within.

- group_order:

  How to order group columns in the output:

  - `"alphabetical"`: sort group names A–Z (default)

  - `"metadata"`: order by first appearance in `metadata[[group_col]]`

  - `"custom"`: use `group_levels` (and append unseen groups at the end)

- group_levels:

  Character vector giving desired group order when
  `group_order="custom"`.

- drop_unmapped_samples:

  Logical; if `TRUE`, drop samples missing a group assignment. If
  `FALSE` (default), error if any sample in `asv_matrix` is missing
  `group_col` in metadata.

## Value

Numeric matrix with groups as columns and taxa/ASVs as rows (mean
abundance per group).

## Examples

``` r
mat <- matrix(c(1, 2, 3, 4), nrow = 2)
rownames(mat) <- c("ASV1", "ASV2")
colnames(mat) <- c("S1", "S2")
meta <- data.frame(SampleID = c("S1", "S2"), Group = c("G1", "G1"))
create_grouped_matrix(mat, meta, sample_col = "SampleID", group_col = "Group")
#>      G1
#> ASV1  2
#> ASV2  3
```
