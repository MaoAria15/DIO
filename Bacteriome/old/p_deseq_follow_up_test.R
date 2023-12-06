setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("4_Analysis0_file_loading_and_prep.R")

  tab_CHP_all <- read.table("stat_result/desep2_genus_top10///FVT_ChP_HFD.tsv", sep = "\t",header = T, row.names = 1)
tab_SDT_all <- read.table("stat_result/desep2_genus_top10/FVT_SDT_HFD.tsv", sep = "\t", header = T, row.names = 1)
tab_PYT_all <- read.table("stat_result/desep2_genus_top10/FVT_PyT_HFD.tsv", sep = "\t",header = T, row.names = 1)
tab_FVT_all <- read.table("stat_result/desep2_genus_top10/Unmodified_FVT_HFD.tsv", sep = "\t", header = T, row.names = 1)
tab_LFD_all <- read.table("stat_result/desep2_genus_top10/Lean_control_HFD.tsv", sep = "\t",header = T, row.names = 1)


#Turn result into a table

# tab_CHP_all$tax <- paste(rownames(tab_CHP_all),tab_CHP_all$Genus, sep = "_")
# tab_SDT_all$tax <- paste(rownames(tab_SDT_all),tab_SDT_all$Genus, sep = "_")
# tab_PYT_all$tax <- paste(rownames(tab_PYT_all),tab_PYT_all$Genus, sep = "_")
# tab_FVT_all$tax <- paste(rownames(tab_FVT_all),tab_FVT_all$Genus, sep = "_")
# tab_LFD_all$tax <- paste(rownames(tab_LFD_all),tab_LFD_all$Genus, sep = "_")







##############################Volcano plot#####################################
# #Load packages needed for making volcano plot
library(EnhancedVolcano)
# #Cut off padj=0.05, log2FoldChange=0.6
p_vol_CHP<-  EnhancedVolcano(tab_CHP_all,
                                  x="log2FoldChange",
                                  y="pvalue",
                                  # lab = tab_CHP_all$tax,
                                  lab = tab_CHP_all$Genus,
                                  pCutoff = 0.05,
                                  FCcutoff = 0.6,
                                  pointSize = 4.0,
                                  labSize = 3.0,
                                  # boxedLabels = TRUE,
                                  title = 'Obese-control vs FVT-ChP')
p_vol_SDT<-  EnhancedVolcano(tab_SDT_all,
                                  x="log2FoldChange",
                                  y="pvalue",
                                  # lab = tab_SDT_all$tax,
                             lab = tab_CHP_all$Genus,
                                  pCutoff = 0.05,
                                  FCcutoff = 0.6,
                                  pointSize = 4.0,
                                  labSize = 3.0,
                                  # boxedLabels = TRUE,
                                  title = 'Obese-control vs FVT-SDT')
p_vol_PYT<-  EnhancedVolcano(tab_PYT_all,
                                  x="log2FoldChange",
                                  y="pvalue",
                                  # lab = tab_PYT_all$tax,
                             lab = tab_CHP_all$Genus,
                                  pCutoff = 0.05,
                                  FCcutoff = 0.6,
                                  pointSize = 4.0,
                                  labSize = 3.0,
                                  # boxedLabels = TRUE,
                                  title = 'Obese-control vs FVT-PyT')
p_vol_FVT<-  EnhancedVolcano(tab_FVT_all,
                                  x="log2FoldChange",
                                  y="pvalue",
                                  # lab = tab_FVT_all$tax,
                             lab = tab_CHP_all$Genus,
                                  pCutoff = 0.05,
                                  FCcutoff = 0.6,
                                  pointSize = 4.0,
                                  labSize = 3.0,
                                  # boxedLabels = TRUE,
                                  title = 'Obese-control vs Unmodified-FVT')
p_vol_LFD<-  EnhancedVolcano(tab_LFD_all,
                                  x="log2FoldChange",
                                  y="pvalue",
                                  # lab =tab_LFD_all$tax,
                             lab = tab_CHP_all$Genus,
                                  pCutoff = 0.05,
                                  FCcutoff = 0.6,
                                  pointSize = 4.0,
                                  labSize = 3.0,
                                  # boxedLabels = TRUE,
                                  title = 'Obese-control vs Lean-control')



ggarrange(p_vol_CHP,p_vol_SDT+rremove("ylab"),p_vol_PYT+rremove("ylab"),
          p_vol_FVT,p_vol_LFD+rremove("ylab"),
          nrow=2, ncol = 3,
          # labels = c("A","B","C","D","E","F"),
          font.label = list(size = 20),
          common.legend = TRUE,
          legend = "bottom")
#
##############################Log2FoldChange###################################




#Make log2foldchange barplot

# Create a placeholder colors vector (replace with your desired color values)

p_CHP <-   ggplot(tab_CHP_all, aes(y = Genus, x = log2FoldChange, color = Phylum)) +
                        geom_vline(xintercept = 0.0, color = "orange", size = 0.5, lty = 2) +
                        geom_point(size = 3) +
                        theme(legend.text = element_text(),
                              title = element_text(size=8),
                              axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5, size = 8),
                              axis.text.y = element_text(size = 8, face = "italic", colour = "black")) +
                        ggtitle("Obese-control vs FVT-ChP") +
                        geom_text(mapping = aes(x = -25, y = 0.75, label = "FDR < 0.01, |LFC| > 2"), color = "red") +
                        scale_x_continuous(limits = c(-10, 10), n.breaks = 10) +
                        scale_y_discrete(expand = c(0.00005, 0.8)) +
                        scale_color_manual(values = c("red","green","yellow","blue","black","white","pink","cyan","gray","orange","brown","purple")) +
                        ylab("Genus")

p_SDT <-   ggplot(tab_SDT_all, aes(y = Genus, x = log2FoldChange, color = Phylum)) +
  geom_vline(xintercept = 0.0, color = "orange", size = 0.5, lty = 2) +
  geom_point(size = 3) +
  theme(legend.text = element_text(face = "bold"),
        title = element_text(size=8),
        axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5, size = 8),
        axis.text.y = element_blank()) +
  ggtitle("vs FVT-SDT") +
  geom_text(mapping = aes(x = -25, y = 0.75, label = "FDR < 0.01, |LFC| > 2"), color = "red") +
  scale_x_continuous(limits = c(-10, 10), n.breaks = 10) +
  scale_y_discrete(expand = c(0.00005, 0.8)) +
  scale_color_manual(values = c("red","green","yellow","blue","black","white","pink","cyan","gray","orange","brown","purple")) +
  ylab("Genus")
p_PYT <-   ggplot(tab_PYT_all, aes(y = Genus, x = log2FoldChange, color = Phylum)) +
  geom_vline(xintercept = 0.0, color = "orange", size = 0.5, lty = 2) +
  geom_point(size = 3) +
  theme(legend.text = element_text(face = "bold"),
        title = element_text(size=8),
        axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5, size = 8),
        axis.text.y = element_blank()) +
  ggtitle("vs FVT-PyT") +
  geom_text(mapping = aes(x = -25, y = 0.75, label = "FDR < 0.01, |LFC| > 2"), color = "red") +
  scale_x_continuous(limits = c(-10, 10), n.breaks = 10) +
  scale_y_discrete(expand = c(0.00005, 0.8)) +
  scale_color_manual(values = c("red","green","yellow","blue","black","white","pink","cyan","gray","orange","brown","purple")) +
  ylab("Genus")
p_FVT <-   ggplot(tab_FVT_all, aes(y = Genus, x = log2FoldChange, color = Phylum)) +
  geom_vline(xintercept = 0.0, color = "orange", size = 0.5, lty = 2) +
  geom_point(size = 3) +
  theme(legend.text = element_text(face = "bold"),
        title = element_text(size=8),
        axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5, size = 8),
        axis.text.y = element_blank()) +
  ggtitle("vs Unmodified-FVT") +
  geom_text(mapping = aes(x = -25, y = 0.75, label = "FDR < 0.01, |LFC| > 2"), color = "red") +
  scale_x_continuous(limits = c(-10, 10), n.breaks = 10) +
  scale_y_discrete(expand = c(0.00005, 0.8)) +
  scale_color_manual(values = c("red","green","yellow","blue","black","white","pink","cyan","gray","orange","brown","purple")) +
  ylab("Genus")
p_LFD <-   ggplot(tab_LFD_all, aes(y = Genus, x = log2FoldChange, color = Phylum)) +
  geom_vline(xintercept = 0.0, color = "orange", size = 0.5, lty = 2) +
  geom_point(size = 3) +
  theme(legend.text = element_text(face = "bold"),
        title = element_text(size=8),
        axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5, size = 8),
        axis.text.y = element_blank()) +
  ggtitle("vs Lean-control") +
  geom_text(mapping = aes(x = -25, y = 0.75, label = "FDR < 0.01, |LFC| > 2"), color = "red") +
  scale_x_continuous(limits = c(-10, 10), n.breaks = 10) +
  scale_y_discrete(expand = c(0.00005, 0.8)) +
  scale_color_manual(values = c("red","green","yellow","blue","black","white","pink","cyan","gray","orange","brown","purple")) +
  ylab("Genus")


ggarrange(p_CHP,p_SDT+rremove("ylab"),p_PYT+rremove("ylab"),
          p_FVT+rremove("ylab"),p_LFD+rremove("ylab"),
          nrow=1, ncol = 5,
          # labels = c("A","B","C","D","E","F"),
          font.label = list(size = 20),
          widths = c(2,1,1,1,1),
          common.legend = TRUE,
          legend = "bottom")

