# Build a dendrogram (and leaf order) from an ASV matrix

Computes a dissimilarity matrix between columns (samples or groups)
using
[`vegan::vegdist()`](https://vegandevs.github.io/vegan/reference/vegdist.html)
on `t(mat)`, clusters with
[`stats::hclust()`](https://rdrr.io/r/stats/hclust.html), and returns a
dendrogram with an optional `hang`.

## Usage

``` r
build_dendrogram(
  mat,
  distance_method = "bray",
  cluster_method = "ward.D2",
  hang = -1
)
```

## Arguments

- mat:

  Numeric matrix with taxa/ASVs as rows and units-to-cluster as columns
  (samples or groups). Columns will be clustered.

- distance_method:

  Distance method passed to
  [`vegan::vegdist()`](https://vegandevs.github.io/vegan/reference/vegdist.html)
  (default: `"bray"`).

- cluster_method:

  Linkage method passed to
  [`stats::hclust()`](https://rdrr.io/r/stats/hclust.html) (default:
  `"ward.D2"`).

- hang:

  Numeric hang parameter passed to
  [`dendextend::hang.dendrogram()`](https://talgalili.github.io/dendextend/reference/hang.dendrogram.html)
  (default: -1).

## Value

A list with:

- dendrogram:

  A `dendrogram` object.

- order:

  Character vector of leaf labels in plotting order.

- hc:

  The `hclust` object.

- dist:

  The `dist` object from
  [`vegan::vegdist`](https://vegandevs.github.io/vegan/reference/vegdist.html).

## Examples

``` r
if (requireNamespace("vegan", quietly = TRUE) &&
  requireNamespace("dendextend", quietly = TRUE)) {
  mat <- matrix(runif(20), nrow = 5)
  colnames(mat) <- c("G1", "G2", "G3", "G4")
  res <- build_dendrogram(mat)
  res$order
}
#> Registered S3 method overwritten by 'dendextend':
#>   method     from 
#>   rev.hclust vegan
#> [1] "G1" "G3" "G2" "G4"
```
