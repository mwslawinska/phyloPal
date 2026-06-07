# Example OTU matrix

A plain matrix of OTU/ASV counts with ASVs as rows and samples as
columns, subset to the samples present in `example_microbiome`. Used as
input for [`create_grouped_matrix()`](create_grouped_matrix.md) and
[`build_dendrogram()`](build_dendrogram.md).

## Usage

``` r
em_otu
```

## Format

A numeric matrix with ASV identifiers as row names and sample
identifiers as column names.

@source Caporaso, J.G., et al. (2011). Global patterns of 16S rRNA
diversity at a depth of millions of sequences per sample. *PNAS*, 108,
4516–4522. Available via
[`phyloseq::GlobalPatterns`](https://rdrr.io/pkg/phyloseq/man/data-GlobalPatterns.html).
