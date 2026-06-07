#' Example microbiome dataset
#'
#' A subset of the \code{GlobalPatterns} dataset from the \pkg{phyloseq}
#' package, filtered to five habitat types (Soil, Ocean, Freshwater,
#' Freshwater creek, Sediment estuary) and reshaped to long format.
#' Relative abundance (RA) is pre-calculated per ASV per sample.
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{SampleID}{Sample identifier}
#'   \item{OTU}{ASV/OTU identifier}
#'   \item{Counts}{Raw read counts}
#'   \item{Depth}{Total sequencing depth per sample}
#'   \item{RA}{Relative abundance (Counts / Depth)}
#'   \item{Habitat}{Broad habitat category (Terrestrial, Oceanic, Freshwater, Brackish)}
#'   \item{Description}{Sample description)}
#'   \item{SampleType}{Original sample type label}
#'   \item{Kingdom}{Taxonomic kingdom}
#'   \item{Phylum}{Taxonomic phylum}
#'   \item{Class}{Taxonomic class}
#'   \item{Order}{Taxonomic order}
#'   \item{Family}{Taxonomic family}
#'   \item{Genus}{Taxonomic genus}
#' }
#'
#'  @source Caporaso, J.G., et al. (2011). Global patterns of 16S rRNA
#'   diversity at a depth of millions of sequences per sample.
#'   \emph{PNAS}, 108, 4516--4522.
#'   Available via \code{phyloseq::GlobalPatterns}.
"example_microbiome"

#' Example sample metadata
#'
#' Sample metadata for the \code{example_microbiome} dataset,
#' containing sample identifiers, sample types, and habitat classifications.
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{SampleID}{Sample identifier}
#'   \item{SampleType}{Original sample type label}
#'   \item{Habitat}{Broad habitat category (Terrestrial, Oceanic, Freshwater, Brackish)}
#'   \item{Description}{Sample description)}
#' }
#'
#'  @source Caporaso, J.G., et al. (2011). Global patterns of 16S rRNA
#'   diversity at a depth of millions of sequences per sample.
#'   \emph{PNAS}, 108, 4516--4522.
#'   Available via \code{phyloseq::GlobalPatterns}.
"em_metadata"

#' Example OTU matrix
#'
#' A plain matrix of OTU/ASV counts with ASVs as rows and samples as columns,
#' subset to the samples present in \code{example_microbiome}.
#' Used as input for \code{create_grouped_matrix()} and \code{build_dendrogram()}.
#'
#' @format A numeric matrix with ASV identifiers as row names and
#'   sample identifiers as column names.
#'
#'  @source Caporaso, J.G., et al. (2011). Global patterns of 16S rRNA
#'   diversity at a depth of millions of sequences per sample.
#'   \emph{PNAS}, 108, 4516--4522.
#'   Available via \code{phyloseq::GlobalPatterns}.
"em_otu"