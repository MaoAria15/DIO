
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("Analysis0_file_loading_and_prep.R")


vir.phyl <- tax_glom(PSB, "Species", NArm = FALSE)
ps0 <- transform_sample_counts(vir.phyl, function(x) x / sum(x))

#Load phyloseq files to ampvis2 format


PSBamp <- amp_load(ps0)
#Bacteriome
PSBamp$metadata$Treatment <- factor(PSBamp$metadata$Treatment, levels = groups)
PSBamp$metadata$Sample_time_point <- factor(PSBamp$metadata$Sample_time_point, levels = c("Arrival","Before_1st_FVT","1w_after_FVT","Termination"))

p_abundance_heatmap <-amp_heatmap(PSBamp,
                          group_by = "Treatment",
                          facet_by = "Sample_time_point",
                          plot_values = FALSE,
                          tax_show = 20,
                          tax_aggregate = "Species" ,
                          tax_empty = "best",
                          showRemainingTaxa = TRUE,
                          normalise = TRUE,
                          color_vector = c("white", "red"),
                          plot_colorscale = "sqrt",
                          plot_legendbreaks = c(1, 10, 25, 50)
              ) +
  scale_x_discrete(breaks=groups,
                   labels=labels_groups_2rows)+
                  theme_classic() +
                  mytheme_abundance

p_abundance_heatmap

# all_variables <- ls()
# variable_to_keep <- "p_abundance_heatmap"
# variables_to_remove <- setdiff(all_variables, variable_to_keep)
# rm(list = variables_to_remove)
# rm(variable_to_keep,variables_to_remove,all_variables)