setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("4_Analysis0_file_loading_and_prep.R")

ps <- subset_samples(PSV.no.realm, Treatment %in% c("Unmodified_FVT","Obese_control","Lean_control"))
ps@sam_data$Liver_response <- paste(ps@sam_data$Treatment,ps@sam_data$Liver_response, sep = "_")
ps <- subset_samples(ps, Sample_time_point %in% "Termination")
ps <- tax_glom(ps, "Genus", NArm = FALSE) #select a level to compare
### loop for all the grouping compared with Obese control
Liver_response <- unique(sample_data(ps)$Liver_response)

Liver_response
Liver_response <- Liver_response[Liver_response!="Obese_control_Response"]
path_table <- "stat_result/desep2-sub/"
dir.create(path_table)

for (group in Liver_response){
  ps.sub <- prune_samples(sample_data(ps)$Liver_response %in% c(group, "Obese_control_Response"), ps)
  ps.sub
  # remove all error taxa
  ps.ds <- phyloseq_to_deseq2(ps.sub, ~Liver_response)
  # solve rows without a zero, deseq need to calculate the geometric zero, 
  cts <- counts(ps.ds)
  geoMeans <- apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
  dds <- estimateSizeFactors(ps.ds, geoMeans=geoMeans)
  ps.ds <-  DESeq(dds, test="Wald", fitType="parametric")
  # result
  res = results(ps.ds, cooksCutoff = FALSE)
  #alpha = 0.0001
  sigtab = res
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
  head(sigtab)
  theme_set(theme_bw())
  scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
  }
  # Phylum order
  # x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
  # x = sort(x, TRUE)
  # sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
  # # Speciea order
  # x = tapply(sigtab$log2FoldChange, sigtab$Species, function(x) max(x))
  # x = sort(x, TRUE)
  # sigtab$Species = factor(as.character(sigtab$Species), levels=names(x))
  write.table(data.frame(sigtab), paste0(path_table,str_replace(group," ",""),"_HFD.tsv"), sep="\t", col.names = NA)
}






















tab_FVTN_all <- read.table("stat_result/desep2-sub/Unmodified_FVT_Non_Response_HFD.tsv", sep = "\t",header = T, row.names = 1)
tab_FVT_all <- read.table("stat_result/desep2-sub/Unmodified_FVT_Response_HFD.tsv", sep = "\t", header = T, row.names = 1)
tab_HFDN_all <- read.table("stat_result/desep2-sub/Obese_control_Non_Response_HFD.tsv", sep = "\t",header = T, row.names = 1)
tab_LFD_all <- read.table("stat_result/desep2-sub/Lean_control_Non_Response_HFD.tsv", sep = "\t", header = T, row.names = 1)





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
p_vol_FVTN<-  EnhancedVolcano(tab_FVTN_all,
                                  x="log2FoldChange",
                                  y="pvalue",
                                  lab = tab_FVTN_all$Genus,
                                  # lab = tab_FVTN_all$Species,
                                  pCutoff = 0.05,
                                  FCcutoff = 0.6,
                                  pointSize = 4.0,
                                  labSize = 3.0,
                                  # boxedLabels = TRUE,
                                  title = 'Obese-control vs Unmodified-FVT non response')
p_vol_FVT<-  EnhancedVolcano(tab_FVT_all,
                                  x="log2FoldChange",
                                  y="pvalue",
                                  lab = tab_FVT_all$Genus,
                             # lab = tab_CHP_all$Species,
                                  pCutoff = 0.05,
                                  FCcutoff = 0.6,
                                  pointSize = 4.0,
                                  labSize = 3.0,
                                  # boxedLabels = TRUE,
                                  title = 'Obese-control vs Unmodified-FVT response')
p_vol_HFDN<-  EnhancedVolcano(tab_HFDN_all,
                                  x="log2FoldChange",
                                  y="pvalue",
                                  lab = tab_HFDN_all$Genus,
                             # lab = tab_CHP_all$Species,
                                  pCutoff = 0.05,
                                  FCcutoff = 0.6,
                                  pointSize = 4.0,
                                  labSize = 3.0,
                                  # boxedLabels = TRUE,
                                  title = 'Obese-control vs Obese-control non response')
p_vol_LFD<-  EnhancedVolcano(tab_LFD_all,
                                  x="log2FoldChange",
                                  y="pvalue",
                                  lab = tab_LFD_all$Genus,
                             # lab = tab_CHP_all$Species,
                                  pCutoff = 0.05,
                                  FCcutoff = 0.6,
                                  pointSize = 4.0,
                                  labSize = 3.0,
                                  # boxedLabels = TRUE,
                                  title = 'Obese-control vs Lean-control')




ggarrange(p_vol_FVTN,p_vol_FVT+rremove("ylab"),p_vol_HFDN,
          p_vol_LFD+rremove("ylab"),
          nrow=2, ncol = 2,
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

