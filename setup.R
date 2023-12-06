library(ggplot2)
library(MESS)
library(forcats)
library(tidyverse)
library(psycho)
library(ggpubr)
library(ggsci)
library(rstatix)
library(ampvis2)
library(DESeq2)
library(plyr)
library(cowplot)
library(RVAideMemoire)
library(data.table)
library(microbiome)
library(forestmangr)
library(writexl)
library(viridis)
library(readxl)
library(phyloseq)      # necessary to import the data from Excel file
library(dplyr)        # filter and reformat data frames
library(stringr)
library(vegan)
library(metagenomeSeq)
library(tidyr)
library(RColorBrewer)
library(reshape2)
library(writexl)
library(directlabels)
library(glue)
library(ggpmisc)
library(ComplexHeatmap)
library(magick)
library(colorRamp2)
library(DESeq2)
library(circlize)
library(microDecon)
library(psych)
library(pheatmap)
library(igraph)
library(ggraph)
library(influential)
library(showtext)
library(gridExtra)
cols <- c("#F27970","#BB9727","#54B345","#32B897","#05B9E2","#8983BF")

cols_his <- c("grey","grey51","grey12")

group_col <- c(FVT_ChP="#F27970", FVT_SDT="#BB9727", FVT_PyT="#54B345",Unmodified_FVT="#32B897",Obese_control="#05B9E2", Lean_control="#8983BF")
cytokine_cols <- c("#F27970","#32B897","#05B9E2","#8983BF")
Obse_compare <-list(c("Obese_control","FVT_ChP"),
                    c("Obese_control","FVT_SDT"),
                    c("Obese_control","FVT_PyT"),
                    c("Obese_control","Unmodified_FVT"))
cytokine_Obse_compare <-list(c("Obese_control","FVT_ChP"),
                             c("Obese_control","Unmodified_FVT"))
neg_compare <- list(c("Obese_control","Lean_control"))



cytokine_groups <- c("FVT_ChP", "Unmodified_FVT", "Obese_control","Lean_control")
groups <- c("FVT_ChP", "FVT_SDT", "FVT_PyT","Unmodified_FVT", "Obese_control","Lean_control")
labels_groups <- c("FVT-ChP", "FVT-SDT", "FVT-PyT","Unmodified-FVT", "Obese-control","Lean-control")
labels_groups_2rows <- c("FVT\nChP", "FVT\nSDT", "FVT\nPyT","Unmodified\nFVT", "Obese\ncontrol","Lean\ncontrol")


#Subtest 
liver_compare_pos <- list(c("Unmodified_FVT_Response","Obese_control_Non_Response"),
                          c("Unmodified_FVT_Non_Response","Obese_control_Non_Response"))
liver_compare_neg <- list(c("Obese_control_Non_Response","Lean_control_Response"))


liver_group <- c("Unmodified_FVT_Response","Unmodified_FVT_Non_Response","Obese_control_Non_Response","Lean_control_Response")
liver_cols <- c("#32A897","#32C897","#05B9E2","#8983BF")






mytheme_with_x <- theme(text = element_text(size = 8, colour = "Black",family = "Arial"),
                 axis.line=element_line(size=0.5),
                 #panel.border = element_blank(),
                 axis.text=element_text(size = 8, colour = "Black"),
                 axis.ticks=element_line(size=1, colour = "Black"),
                 strip.background = element_rect(colour = "white", fill = "white"),
                 axis.text.x=element_text(size= 8, angle = 0,vjust = 0.6),
                 axis.title = element_text(size = 8, face = "bold"),
                 strip.text.x = element_text(angle = 30, size=8, face = "bold"),
                 legend.text = element_text(size=8),
                 legend.key.size = unit(8, "pt"),
                 legend.title = element_text(size = 8,face = "bold"),
                 title = element_text(size =8, face = "bold")
)
mytheme <- theme(text = element_text(size = 8, colour = "Black",family = "Arial"),
                 axis.line=element_line(size=0.5),
                 #panel.border = element_blank(),
                 axis.text=element_text(size = 8, colour = "Black"),
                 axis.ticks=element_line(size=1, colour = "Black"),
                 axis.ticks.x = element_blank(),
                 strip.background = element_rect(colour = "white", fill = "white"),
                 axis.text.x=element_blank(),
                 axis.title = element_text(size = 8, face = "bold"),
                 strip.text.x = element_text(angle = 30, size=8, face = "bold"),
                 legend.text = element_text(size=8),
                 legend.key.size = unit(8, "pt"),
                 legend.title = element_text(size = 8,face = "bold"),
                 title = element_text(size =8, face = "bold")
)
mytheme_alpha <- theme(text = element_text(size = 8, colour = "Black",family = "Arial"),
                 axis.line=element_line(size=0.5),
                 #panel.border = element_blank(),
                 axis.text=element_text(size = 8, colour = "Black"),
                 axis.ticks=element_line(size=1, colour = "Black"),
                 axis.ticks.x = element_blank(),
                 strip.background = element_rect(colour = "white", fill = "white"),
                 # axis.text.x=element_blank(),
                 axis.title = element_text(size = 8, face = "bold"),
                 strip.text.x = element_text(angle = 30, size=8, face = "bold"),
                 legend.text = element_text(size=8),
                 legend.key.size = unit(8, "pt"),
                 legend.title = element_text(size = 8,face = "bold"),
                 title = element_text(size =8, face = "bold")
)
mytheme_beta <- theme(text = element_text(size = 8, colour = "Black",family = "Arial"),
                       axis.line=element_line(size=0.5),
                       #panel.border = element_blank(),
                       axis.text=element_text(size = 8, colour = "Black"),
                       axis.ticks=element_line(size=1, colour = "Black"),
                       strip.background = element_rect(colour = "white", fill = "white"),
                       axis.title = element_text(size = 8, face = "bold"),
                       strip.text.x = element_text(angle = 30, size=8, face = "bold"),
                       legend.text = element_text(size=8),
                       legend.key.size = unit(8, "pt"),
                       legend.title = element_text(size = 8,face = "bold"),
                       title = element_text(size =8, face = "bold")
)
mytheme_abundance <- theme(text = element_text(size = 8, colour = "Black",family = "Arial"),
                      axis.line=element_line(size=0.5),
                      #panel.border = element_blank(),
                      axis.text=element_text(size = 8, colour = "Black"),
                      axis.text.x=element_text(size = 8,angle = 0,colour = "Black"),
                      axis.text.y=element_text(size = 8,face = "italic",colour = "Black"),
                      axis.ticks=element_line(size=1, colour = "Black"),
                      strip.background = element_rect(colour = "white", fill = "white"),
                      axis.title = element_text(size = 8, face = "bold"),
                      strip.text.x = element_text(angle = 0, size=8, face = "bold"),
                      legend.text = element_text(size=8),
                      legend.key.size = unit(8, "pt"),
                      legend.title = element_text(size = 8,face = "bold"),
                      title = element_text(size =8, face = "bold")
)
mytheme_his <-  theme(text = element_text(size = 8, colour = "Black",family = "Arial"),
                      axis.line=element_line(size=0.5),
                      #panel.border = element_blank(),
                      axis.text=element_text(size = 8, colour = "Black"),
                      axis.ticks=element_line(size=1, colour = "Black"),
                      strip.background = element_rect(colour = "white", fill = "white"),
                      axis.text.x=element_text(size= 8, angle = 0,vjust = 0.6),
                      axis.title = element_text(size = 8, face = "bold"),
                      strip.text.x = element_text(angle = 30, size=8, face = "bold"),
                      legend.text = element_text(size=8),
                      legend.key.size = unit(8, "pt"),
                      legend.title = element_text(size = 8,face = "bold"),
                      title = element_text(size =8, face = "bold"),
                      legend.background = element_rect(color = "black",linewidth = 1))
  Data <- read.delim('data/p_data/phenotype_metadata.txt', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
Data$Group <- fct_relevel(Data$Group, groups) #Releveling to control


filter_and_replace <- function(data, threshold = 0.05) {
  data %>%
    filter(p < threshold) %>%
    mutate(p_signif = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE      ~ ""
    )) %>%
    select(-p)  # Remove the original p-value column
}


