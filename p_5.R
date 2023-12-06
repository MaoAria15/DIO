setwd("Virome_2")

source("p_vir_alpha.R")
source("p_vir_beta_bray.R")
source("p_host_abundance_heatmap.R")
source("Virome_2/p_vir_abundance_heatmap.R")
p_vir_corr_heatmap_sim
p_abundance_heatmap_vir
p_abundance_heatmap_host
p_vir_Deseq_heatmap_simple
p_alpha_vir
p_beta_vir
# p_arrival_stat,
# p_afterfvt_stat,
# p_termination_stat,
ggarrange(p_alpha_vir,
          p_beta_vir,
          # p_abundance_heatmap_vir,
          # p_abundance_heatmap_host,
          labels = c("A","B","C","D"),
          ncol = 1,
          nrow = 4,
          heights = c(3,3,3,3),
          commen.legend= FALSE
          )
