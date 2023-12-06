
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("4_Analysis0_file_loading_and_prep.R")


vir.phyl <- tax_glom(PSV.no.realm, "Genus", NArm = FALSE)
ps0 <- transform_sample_counts(vir.phyl, function(x) x / sum(x))

#Load phyloseq files to ampvis2 format


PSVamp <- amp_load(ps0)
#Bacteriome
PSVamp$metadata$Treatment <- factor(PSVamp$metadata$Treatment, levels = groups)
PSVamp$metadata$Sample_time_point <- factor(PSVamp$metadata$Sample_time_point, levels = c("Arrival","1w_after_FVT","Termination"))

p_abundance_heatmap_vir <-amp_heatmap(PSVamp,
                                  group_by = "Treatment",
                                  facet_by = "Sample_time_point",
                                  plot_values = FALSE,
                                  tax_show = 15,
                                  tax_aggregate = "Genus" ,
                                  tax_empty = "best",
                                  showRemainingTaxa = TRUE,
                                  normalise = TRUE,
                                  color_vector = c("white", "blue"),
                                  plot_colorscale = "sqrt",
                                  plot_legendbreaks = c(1, 10, 25, 50)
) +
  scale_x_discrete(breaks=groups,
                   labels=labels_groups_2rows)+
  theme_classic() +
  mytheme_abundance
p_abundance_heatmap_vir

