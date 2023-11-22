setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("p_rel_abundance_bac.R")
source("p_FvsB_ratio.R")
source("p_abundance_barplot_termination.R")


p_supplementary_fig3<-ggarrange(ggarrange(p_Allobaculum,p_Bacteroides,p_Bifidobacterium,p_Clostridium,
                    p_Dehalobacterium,p_Lactobacillus,p_Lactococcus,p_Oscillospira,
                    ncol = 4,
                    nrow = 2,
                    labels = c("A","B","C","D","E","F","G","H"),
                    common.legend = TRUE),
          ggarrange(p_BvsF,p_abundance_bar_termination,
                    ncol = 2,
                    nrow = 1,
                    widths = c(2,3),
                    labels = c("I","J")
                    ),
          nrow = 2,
          heights = c(2,1),
          widths = c(4,4))
