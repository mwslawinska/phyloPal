# Generates figures for README.md
# Run once manually: source("data-raw/prepare_readme_figures.r")
# Requires: devtools::install() first

options(bitmapType = "cairo")  # for HPC/no X11

library(phyloPal)
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggplotify)

dir.create("man/figures", recursive = TRUE, showWarnings = FALSE)

data(example_microbiome)
data(em_metadata)
data(em_otu)

### this is going to be the package example data
# clean taxonomy
em_cleaned <- example_microbiome %>%
replace_incertae_sedis_NAs()

###########
# add alpha
cols <- c(
  add_alpha("#FF0000", 0.2),
  add_alpha("#FF0000", 0.5),
  add_alpha("#FF0000", 0.8),
  add_alpha("#FF0000", 1)
)

df_alpha <- data.frame(
  x = factor(c("alpha 0.2", "alpha 0.5", "alpha 0.8", "alpha 1.0"),
             levels = c("alpha 0.2", "alpha 0.5", "alpha 0.8", "alpha 1.0")),
  y = 1,
  col = cols
)

p <- ggplot(df_alpha, aes(x = x, y = y, fill = col)) +
  geom_col(width = 0.9, color = NA) +
  scale_fill_identity() +
  labs(
    title = "Effect of alpha transparency",
    x = NULL,
    y = NULL
  ) +
  theme_phylopal() +
  scale_y_continuous(expand = c(0,0)) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    aspect.ratio = 2
  ) 

ggsave(
  filename = paste0("man/figures/alpha_barplot.png"),
  plot = p,
  width = 8,
  height = 6,
  units = "in",
  dpi = 150
)

###########
# Barplot workflow
# Process data for barplot
em_barplot_processed <- process_barplot_data(
  em_cleaned,
  tax_level = "Class",
  group_vars = c("SampleType", "SampleID", "Habitat"),
  low_abundance_basis = "per_sample",
  low_abundance_threshold = 0.01,
  agg_fun = "sum",
  keep_ratype = "separate",
  clean_taxonomy = FALSE
)


em_barplot_processed2 <- process_barplot_data(
  example_microbiome,
tax_level = "Class",
group_vars = c("SampleType", "SampleID", "Habitat"),
low_abundance_threshold = 0.01,
keep_ratype = "collapse",
clean_taxonomy = TRUE,
preserve_higher_taxonomy = T
)

barplot_pal <- generate_palette_hcl(
  data = em_barplot_processed,
  tax_level = "Class",
  fixed_colors_enabled = TRUE,
  fixed_colors_position = "end",
  palette_list = c("Reds", "Purples", "BrwnYl", "Blues", "TealGrn"),
  cmax = 65,
  luminance = c(20,90),
    power = 1.2,
    shuffle = FALSE)

habitat_palette <- generate_grouped_palette(
  data = em_cleaned,
  group_col = "Habitat",
  item_col = "SampleType",
  palette_map = list(
    "Terrestrial" = "BrwnYl",
    "Oceanic" = "Blues",
    "Freshwater" = "Greens",
    "Brackish" = "PuRd"
  ),
  luminance = 65,
  power = 1.2
)
p_barplot <- plot_taxonomic_barplot(
  data = em_barplot_processed,
  tax_level = "Class",
  palette = barplot_pal,
  x_axis_var = "SampleID",
  facet_by = "SampleType",
  facet_strip_colors = habitat_palette,
  theme_obj = theme_phylopal()
) + 
  guides(
    fill = guide_legend(
      ncol = 1
    )
  ) 

ggsave(
  filename = paste0("man/figures/em_barplot.png"),
  plot = p_barplot,
  width = 8,
  height = 6,
  units = "in",
  dpi = 150
)

p_barplot2 <- plot_taxonomic_barplot(
  data = em_barplot_processed2,
  tax_level = "Class",
  palette = barplot_pal,
  x_axis_var = "SampleID",
  facet_by = "SampleType",
  facet_strip_colors = habitat_palette,
  theme_obj = theme_phylopal()
) + 
  guides(
    fill = guide_legend(
      ncol = 1
    )
  ) 

ggsave(
  filename = paste0("man/figures/em_barplot2.png"),
  plot = p_barplot2,
  width = 8,
  height = 6,
  units = "in",
  dpi = 150
)


em_processed_grouped <- process_barplot_data(
  example_microbiome,
  tax_level = "Class",
  group_vars = c("SampleType", "Habitat"),
  normalize_by = c("SampleType", "Habitat"),
  low_abundance_threshold = 0.1,
  preserve_higher_taxonomy = TRUE,
  low_abundance_basis = "per_sample",
  agg_fun = "sum",
  keep_ratype = "separate"
)

barplot_pal_grouped <- generate_palette_hcl(
  data = em_processed_grouped,
  tax_level = "Class",
   group_by_higher_tax = "Phylum",
   order_by_higher_tax = TRUE,
   group_palette_map = list(
     "Actinobacteria" = "Blues",
     "Proteobacteria" = "Greens",
     "Cyanobacteria" = list(palette = "Purple-Orange", side = "right"),
     "Acidobacteria" = list(palette = "Purple-Orange", side = "right"),
     "Bacteroidetes" = "Burg",
     "Verrucomicrobia" = "BrwnYl"),
  fixed_colors_enabled = TRUE,
  fixed_colors_position = "end",
  order_groups = "alphabetical",
  order_within_groups = "alphabetical",
  cmax = 65,
  luminance = c(20,90),
    power = 1.2,
    shuffle = FALSE)

p_barplot_grouped <- plot_taxonomic_barplot(
  data = em_processed_grouped,
  tax_level = "Class",
  palette = barplot_pal_grouped,
  x_axis_var = "SampleType",
  facet_by = "Habitat",
  theme_obj = theme_phylopal() 
) + 
  ggplot2::guides(
    fill = guide_legend(
      ncol = 1
    )
  )  + theme(axis.text.x = element_text(size =11, angle = 45, hjust = 1, vjust = 1),
    axis.ticks.x = ggplot2::element_line(color = "black", linewidth = 0.4))

ggsave(
  filename = paste0("man/figures/em_barplot_grouped.png"),
  plot = p_barplot_grouped,
  width = 8,
  height = 8,
  units = "in",
  dpi = 150
)

# grouped, but palette which isn't grouped by Phylum
em_processed_grouped2 <- process_barplot_data(
  example_microbiome,
  tax_level = "Class",
  group_vars = c("SampleType", "Habitat"),
  normalize_by = c("SampleType", "Habitat"),
  low_abundance_threshold = 0.01,
  preserve_higher_taxonomy = TRUE,
  low_abundance_basis = "per_sample",
  agg_fun = "sum",
  keep_ratype = "separate"
)

p_barplot_grouped2 <- plot_taxonomic_barplot(
  data = em_processed_grouped2,
  tax_level = "Class",
  palette = barplot_pal,
  x_axis_var = "SampleType",
  facet_by = "Habitat",
  theme_obj = theme_phylopal() 
) + 
  ggplot2::guides(
    fill = guide_legend(
      ncol = 1
    )
  )  + theme(axis.text.x = element_text(size =11, angle = 45, hjust = 1, vjust = 1),
    axis.ticks.x = ggplot2::element_line(color = "black", linewidth = 0.4))

ggsave(
  filename = paste0("man/figures/em_barplot_grouped2.png"),
  plot = p_barplot_grouped2,
  width = 8,
  height = 8,
  units = "in",
  dpi = 150
)


###########################
# alluvial plot workflow

example_microbiome$SampleType <- factor(example_microbiome$SampleType, 
levels = unique(example_microbiome$SampleType))

em_allu <- prepare_alluvial_data(example_microbiome,
tax_level = "Class",
group_col = c("SampleType"),
clean_taxonomy = TRUE
)

em_allu_classified <- classify_taxa_patterns(
  data = em_allu,
  tax_level = "Class",
  group_col = c("SampleType")
)

allu_pal <- generate_alluvial_palette(
    data = em_allu_classified,
  palette_list = c("Reds", "Purples", "BrwnYl", "Blues", "TealGrn"),
  cmax = 65,
  luminance = c(20,90),
    power = 1.2,
    )

p_allu <- plot_alluvial(em_allu_classified, 
custom_palette = allu_pal,
tax_level = "Class", 
group_col = "SampleType",
theme_obj = theme_phylopal(),
line_width = 0.2,
x_axis_label = "Sample Type"
) +
theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(
    fill = guide_legend(
      ncol = 1
    )
  )
ggsave(
  filename = paste0("man/figures/em_alluvial.png"),
  plot = p_allu,
  width = 8,
  height = 6,
  units = "in",
  dpi = 150
)

# wrapper to avoid doing everything on your own
em_allu_wrapper <- create_alluvial_plot(
  data = example_microbiome,
  tax_level = "Class",
  group_col = "SampleType",
  prepare_args = list(clean_taxonomy = TRUE),
  palette_list = c("Reds", "Purples", "BrwnYl", "Blues", "TealGrn"),
  palette_args = list(
    cmax = 65,
    luminance = c(20, 90),
    power = 1.2
  ),
  plot_args = list(
    theme_obj = theme_phylopal(),
    line_width = 0.2,
    x_axis_label = "Sample Type"
  )
) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1)) +
  ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1))

  ggsave(
  filename = paste0("man/figures/em_alluvial_wrapper.png"),
  plot = em_allu_wrapper,
  width = 8,
  height = 6,
  units = "in",
  dpi = 150
)



em_allu2 <- prepare_alluvial_data(example_microbiome,
tax_level = "Order",
group_col = c("SampleType"),
clean_taxonomy = TRUE,
preserve_higher_taxonomy = T
)

em_allu_classified2 <- classify_taxa_patterns(
  data = em_allu2,
  tax_level = "Order",
  group_col = c("SampleType"),
  low_abundance_threshold = 0.1
)

allu_pal2 <- generate_alluvial_palette(
  data = em_allu_classified2,
  tax_level = "Order",
  group_by_higher_tax = "Class",  
  group_palette_map = list(         
    Actinobacteria = "Blues",
    Nostocophycideae = list(palette = "Purple-Brown", side = "right"),
    Chloroplast = "Greens",
    Flavobacteria = list(palette = "Blue-Yellow 3", side = "right"),
    Gammaproteobacteria = "Reds",
    Deltaproteobacteria = list(palette = "Purple-Brown", side = "left"),
    Spartobacteria = "Burg",
    Betaproteobacteria = "Heat"
  )
)

p_allu_group <- plot_alluvial(em_allu_classified2, 
custom_palette = allu_pal2,
tax_level = "Order", 
group_col = "SampleType",
theme_obj = theme_phylopal(),
line_width = 0.2,
x_axis_label = "Sample Type"
) +
ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  guides(
    fill = guide_legend(
      ncol = 1
    )
  )
ggsave(
  filename = paste0("man/figures/em_alluvial_grouped.png"),
  plot = p_allu_group,
  width = 8,
  height = 6,
  units = "in",
  dpi = 150
)

#############
# alluvial with dendrogram


em_otu_grouped <- create_grouped_matrix(
asv_matrix = em_otu,
metadata = em_metadata,
sample_col = "SampleID",
group_col= "SampleType",
group_order = "metadata"
)

em_dendrogram <- build_dendrogram(
  mat = em_otu_grouped,
  distance_method = "bray",
  cluster_method = "ward.D2"
)

em_dendrogram_plot <- plot_dendrogram(
  dend = em_dendrogram,
  metadata = em_metadata,
  label_from = "SampleType",      
  color_by = "SampleType",
  color_palette = habitat_palette,
  point_size = 2,
  orientation = "top",
  shape_by = "Habitat",
  theme_obj = theme_void() + theme(text = element_text(size = 7, color = "black"),
  legend.title = element_text(size = 7, color = "black"),)
)


ggsave(
  filename = paste0("man/figures/em_dend.png"),
  plot = em_dendrogram_plot,
  width = 8,
  height = 4,
  units = "in",
  dpi = 150
)

p_allu4dend <- create_alluvial_plot(
  data = example_microbiome,
  tax_level = "Class",
  group_col = "SampleType",
  prepare_args = list(clean_taxonomy = TRUE),
  palette_list = c("Reds", "Purples", "BrwnYl", "Blues", "TealGrn"),
  palette_args = list(
    cmax = 65,
    luminance = c(20, 90),
    power = 1.2
  ),
  plot_args = list(
    theme_obj = theme_phylopal(),
    line_width = 0.2,
    x_axis_label = "Sample Type"
  )
) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
  ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1))


# dendrogram and alluvial with legends
dendrogram_alluplot <- combine_dendrogram_alluvial(
  alluvial_plot   = p_allu4dend +
  scale_y_continuous(expand = c(0,0), breaks = seq(0,1,0.1), limits = c(0,1))+
  ggplot2::guides(fill = guide_legend(ncol =1, title = "Class")),
  dendrogram_plot = em_dendrogram_plot +
  ggplot2::guides(color = guide_legend(ncol = 2, title = "Sample Type"), shape = guide_legend(ncol = 2)),
  dend_position   = "top",
  dend_height     = 0.15,
  strip_alluvial_x = FALSE,
  legend          = "separate",
  legend_source   = "both",       
  legend_position = "right",
  legend_rel_width = 0.75,            
  alluvial_margins    = ggplot2::margin(0, 0, 0, 0, unit = "cm"),
  dendrogram_margins    = ggplot2::margin(0, 0, 0.15, 0, unit = "cm"),
  outer_margins    = ggplot2::margin(0.2, 0.2, 0.2, 0.2, unit = "cm"),
  align = "panel",
  x_expand_zero = TRUE,
  align_x_centers = TRUE,
  leaf_order = em_dendrogram$order,
  overwrite_x_scales = TRUE,
  dend_limits_left = 0.4,  
  dend_limits_right = 0.18
) 


ggsave( 
  filename = paste0("man/figures/em_dend_allu_legend.png"),
  plot = dendrogram_alluplot,
  width = 8, height = 9, 
  units = "in", dpi = 150
)

# dendrogram + alluvial without legends
dendrogram_alluplot <- combine_dendrogram_alluvial(
  alluvial_plot   = p_allu4dend,
  dendrogram_plot = em_dendrogram_plot,
  dend_position   = "top",
  dend_height     = 0.15,
  strip_alluvial_x = FALSE,
  legend          = "omit",
  legend_source   = "both",       
  legend_position = "right",
  legend_rel_width = 0.75,            
  alluvial_margins    = ggplot2::margin(0, 0, 0, 0, unit = "cm"),
  dendrogram_margins    = ggplot2::margin(0, 0, 0.15, 0, unit = "cm"),
  outer_margins    = ggplot2::margin(0.2, 0.2, 0.2, 0.2, unit = "cm"),
  align = "panel",
  x_expand_zero = TRUE,
  align_x_centers = TRUE,
  leaf_order = em_dendrogram$order,
  overwrite_x_scales = TRUE,
  dend_limits_left = 0.54,  
  dend_limits_right = 0.16
) 


ggsave( 
  filename = paste0("man/figures/em_dend_allu.png"),
  plot = dendrogram_alluplot,
  width = 8, height = 6, units = "in", dpi = 150
)

# If you want to have more control over the legends

dendrogram_alluplot_control <- combine_dendrogram_alluvial(
  alluvial_plot   = p_allu4dend,
  dendrogram_plot = em_dendrogram_plot,
  dend_position   = "top",
  dend_height     = 0.15,
  strip_alluvial_x = TRUE,
  legend          = "omit",
  legend_source   = "both",       
  legend_position = "right",
  legend_rel_width = 0.75,            
  alluvial_margins    = ggplot2::margin(0, 0, 0, 0, unit = "cm"),
  dendrogram_margins    = ggplot2::margin(0, 0, 0.15, 0, unit = "cm"),
  outer_margins    = ggplot2::margin(0.2, 0.2, 0.2, 0.2, unit = "cm"),
  align = "panel",
  x_expand_zero = TRUE,
  align_x_centers = TRUE,
  leaf_order = em_dendrogram$order,
  overwrite_x_scales = TRUE,
  dend_limits_left = 0.3,  
  dend_limits_right = 0.18
) 


legend_alluplot <- ggpubr::get_legend(p_allu4dend + 
                                            guides(
                                              fill = guide_legend(ncol =1, title = "Class")))


legend_dendrogram <- ggpubr::get_legend(
  em_dendrogram_plot + 
  ggplot2::guides(color = guide_legend(ncol = 1, title = "Sample Type"), shape = guide_legend(ncol = 1))+
    ggplot2::theme(legend.position = "right", 
          legend.box = "vertical",
          legend.title.position = "top",
          plot.margin = margin(0,0,0,0))
)

p_allu_full <- plot_grid(
  plot_grid(
    as.ggplot(dendrogram_alluplot_control),
    plot_grid(
      legend_dendrogram,
      legend_alluplot,
      rel_heights = c(0.6, 1),
      rel_widths = c(1, 1),
      ncol = 1,
      align = "hv", axis = "tblr"
    ),
    # rel_heights = c(1,1),
    rel_widths = c(1,0.6),
    ncol = 2,
    align = "hv", axis = "tblr"
  )
)


ggsave( 
  filename = paste0("man/figures/em_dend_allu_control.png"),
  plot = p_allu_full,
  width = 8, height = 9, units = "in", dpi = 150
)

res <- create_alluvial_dendrogram_plot(
  asv_matrix = em_otu,
  metadata = em_metadata,
  sample_col = "SampleID",
  group_col  = "SampleType",
  alluvial_data = example_microbiome,
  tax_level = "Class",
  dend_color_palette = habitat_palette,
  dend_shape_by = "Habitat",
  theme_alluvial = theme_phylopal(),
  theme_dendrogram = ggplot2::theme_void(),
  alluvial_args = list(
    return_all = TRUE,
    prepare_args = list(clean_taxonomy = TRUE),
    classify_args = list(low_abundance_threshold = 0.01),
    palette_args = list(
      palette_list = c("Reds", "Purples", "BrwnYl", "Blues", "TealGrn"),
      cmax = 65,
      luminance = c(20, 90),
      power = 1.2
    ),
    plot_args = list(
      line_width = 0.2,
      x_axis_label = "Sample Type"
    )
  ),
  post_plot_guides   = list(      # guides applied to alluvial
    fill = ggplot2::guide_legend(ncol = 1, title = "Class")
  ),    
  dend_limits_left = 0.4,  
  dend_limits_right = 0.18, 
  combine_args = list(
    legend_rel_width = 0.5,
    strip_alluvial_x = TRUE,  
    alluvial_margins = ggplot2::margin(0, 0, 0, 0, unit = "cm"),
    outer_margins    = ggplot2::margin(0.2, 0.5, 0.2, 0.2, unit = "cm") 
  )
)

p_em_alluvial_dend_wrapper <- res$combined_plot

ggsave(
  filename = paste0("man/figures/em_alluvial_dend_wrapper.png"),
  plot = p_em_alluvial_dend_wrapper,
  width = 8,
  height = 10,
  units = "in",
  dpi = 150
)
