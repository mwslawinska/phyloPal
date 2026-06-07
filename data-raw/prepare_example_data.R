library(phyloseq)   # for GlobalPatterns
library(dplyr)
library(tidyr)
library(tibble)
library(phyloPal)

gp <- data(GlobalPatterns)

gp_otu <- otu_table(GlobalPatterns)
gp_tax <- tax_table(GlobalPatterns)

gp_sample <- sample_data(GlobalPatterns) %>%
    as.data.frame() %>%
    as_tibble() %>%
    dplyr::rename("SampleID" = "X.SampleID") %>%
    filter(SampleType %in% c("Soil", "Ocean", "Freshwater", "Freshwater (creek)",
    "Sediment (estuary)")) %>%
    mutate(Habitat = case_when(
        SampleType == "Soil" ~ "Terrestrial",
        SampleType == "Ocean" ~ "Oceanic",
        SampleType %in% c("Freshwater", "Freshwater (creek)") ~ "Freshwater",
        SampleType == "Sediment (estuary)" ~ "Brackish",
        TRUE ~ "Other"
    ))


gp_otu_l <- gp_otu %>%
    as.data.frame() %>%
    rownames_to_column(var = "OTU") %>%
    as_tibble() %>%
    pivot_longer(cols = -OTU, names_to = "Sample", values_to = "Counts") %>%
    filter(Sample %in% gp_sample$SampleID) %>%
    dplyr::rename("SampleID" = "Sample")

gp_tax_l <- gp_tax %>%
    as.data.frame() %>%
    rownames_to_column(var = "OTU") %>%
    as_tibble() %>%
    filter(OTU %in% gp_otu_l$OTU)

#calculate depth
depth_df <- gp_otu_l %>%
    group_by(SampleID) %>%
    dplyr::summarise(Depth = sum(Counts)) %>%
    ungroup()

# merge all together
example_microbiome <- gp_otu_l %>%
    inner_join(gp_tax_l, by = "OTU") %>%
    inner_join(gp_sample, by = "SampleID") %>%
    inner_join(depth_df, by = "SampleID") %>%
    mutate(RA = Counts / Depth) %>%
    dplyr::select(SampleID, OTU, Counts, Depth, RA,  
                SampleType, Habitat, Description,
                Kingdom, Phylum, Class, Order, Family, Genus)

em_otu <- gp_otu[, colnames(gp_otu) %in% example_microbiome$SampleID]
em_otu <- as(em_otu, "matrix")

em_metadata <- gp_sample %>%
  dplyr::select(SampleID, SampleType, Habitat, Description)

usethis::use_data(example_microbiome, em_metadata, em_otu, overwrite = TRUE)
