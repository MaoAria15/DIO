
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("Bacteriome/p_alpha_diversity.R")
source("Bacteriome/p_unifrac_diversity.R")
source("Bacteriome/p_abundance_heatmap.R")
source("Bacteriome/p_DESeq_heatmap.R")
source("Correlation_analysis/p_node_network_bac_fatfacs.R")
source("Correlation_analysis/p_node_network_bac_pheno.R")
p1<-ggarrange(p_alpha_bac,p_unifac_stat_bac,
          ncol = 2,
          nrow = 1,
          heights = c(1,1),
          widths = c(1,1),
          common.legend = TRUE,
          labels= c("A","B"),
          legend = "top")

p_fig3_1<-ggarrange(p1, p_abundance_heatmap,
          ncol = 1,
          nrow = 2,
          heights = c(1,1.5),
          labels = c(NA,"C"))


p_fig3_3<-ggarrange(p_bac_pheno_node,p_bac_fats_node,
          ncol = 2,
          nrow = 1,
          common.legend = TRUE,
          labels  = c("E","F"),
          legend = "top")


ggarrange(p_fig3_1,p_fig3_3,
          nrow = 2,
          ncol = 1,
          height = c(3,2),
          common.legend = FALSE)

p_Deseq_heatmap #p_fig3_2


p_arrival_stat  #Manually added
p_beforefvt_stat #Manually added
p_afterfvt_stat #Manually added
p_termination_stat# Manually added 