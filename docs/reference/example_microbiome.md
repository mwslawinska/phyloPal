# Example microbiome dataset

A subset of the `GlobalPatterns` dataset from the phyloseq package,
filtered to five habitat types (Soil, Ocean, Freshwater, Freshwater
creek, Sediment estuary) and reshaped to long format. Relative abundance
(RA) is pre-calculated per ASV per sample.

## Usage

``` r
example_microbiome
```

## Format

A data frame with columns:

- SampleID:

  Sample identifier

- OTU:

  ASV/OTU identifier

- Counts:

  Raw read counts

- Depth:

  Total sequencing depth per sample

- RA:

  Relative abundance (Counts / Depth)

- Habitat:

  Broad habitat category (Terrestrial, Oceanic, Freshwater, Brackish)

- Description:

  Sample description)

- SampleType:

  Original sample type label

- Kingdom:

  Taxonomic kingdom

- Phylum:

  Taxonomic phylum

- Class:

  Taxonomic class

- Order:

  Taxonomic order

- Family:

  Taxonomic family

- Genus:

  Taxonomic genus

@source Caporaso, J.G., et al. (2011). Global patterns of 16S rRNA
diversity at a depth of millions of sequences per sample. *PNAS*, 108,
4516–4522. Available via
[`phyloseq::GlobalPatterns`](https://rdrr.io/pkg/phyloseq/man/data-GlobalPatterns.html).
