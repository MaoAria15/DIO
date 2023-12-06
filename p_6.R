setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("Phenotype/sub_test1.R")
source("Bacteriome/sub_beta_diversity .R")
source("Bacteriome/sub_rel_abundance_bac.R")


p_termination_stat # add manually



ggarrange(ggarrange(p_histo,  p_ewat_mg,p_Weight_gain_perc,
                    p_termination,p_AUC1_w18, p_ogtt_w18,
          common.legend = TRUE,
          legend = "top",
          ncol = 3,
          nrow = 2,
          labels = c("A","B","C","D","E","F")
          ),
          ggarrange(p_Allobaculum,p_Bacteroides,p_Bifidobacterium,p_Clostridium,
          p_Dehalobacterium,p_Lactobacillus,p_Lactococcus,p_Oscillospira,
          ncol = 4,
          nrow = 2,
          labels = c("G","H","I","J","K","L","M","N"),
          common.legend = TRUE),
          ncol = 1,
          nrow = 2
          )
