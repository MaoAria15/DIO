setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("4_Analysis0_file_loading_and_prep.R")
ps <- subset_samples(PSV.no.realm, Sample_time_point %in% "Termination")
#################################################################################
#draw heatmap 

#Include deseq2 p.adj  taxonomy for rows
#Include Litter Information and Group NEC Information for colors
#load two list of deseq2
tab_CHP_all <- read.table("stat_result/desep2/FVT_ChP_HFD.tsv", sep = "\t",header = T, row.names = 1)
tab_SDT_all <- read.table("stat_result/desep2/FVT_SDT_HFD.tsv", sep = "\t", header = T, row.names = 1)
tab_PYT_all <- read.table("stat_result/desep2/FVT_PyT_HFD.tsv", sep = "\t",header = T, row.names = 1)
tab_FVT_all <- read.table("stat_result/desep2/Unmodified_FVT_HFD.tsv", sep = "\t", header = T, row.names = 1)
tab_LFD_all <- read.table("stat_result/desep2/Lean_control_HFD.tsv", sep = "\t",header = T, row.names = 1)



tab_CHP <- subset(tab_CHP_all, padj < 0.05)
tab_SDT <-subset(tab_SDT_all, padj < 0.05)
tab_PYT <- subset(tab_PYT_all, padj < 0.05)
tab_FVT <-subset(tab_FVT_all, padj < 0.05)
tab_LFD <- subset(tab_LFD_all, padj < 0.05)

OTU <- unique(c(rownames(tab_CHP),
                rownames(tab_SDT),
                rownames(tab_PYT),
                rownames(tab_FVT),
                rownames(tab_LFD)
))


ps.rel <- transform_sample_counts(ps, function(x) x/sum(x)*100)
ps.rel.sig <- prune_taxa(rownames(otu_table(ps.rel)) %in% OTU,ps.rel)

#select the rel-abun > 0.1%

# at least 1% relative abundance appearance in 10% samples
# mat <- as.matrix(otu_table(ps.rel.sig))
# species2keep <- rownames(mat)[rowSums(mat>=1)/length(colnames(mat))> 0.1]
# ps.rel.sig <- prune_taxa(species2keep,ps.rel.sig)

otu_abun_select <- data.frame(otu_table(ps.rel.sig), check.names = F)

#import relavant metadata
metadata <- data.frame(sample_data(ps.rel.sig))
tax.clean <- data.frame(tax_table(ps.rel.sig))


# create a variable to define the subgroup
# order the matrix by the subgroup
metadata$Treatment <- factor(metadata$Treatment, levels = groups)
meta_order <- metadata[order(metadata$Treatment),]

# re_order the col
mat <- otu_abun_select
mat <- as.matrix(mat[,rownames(meta_order)])

base_mean = rowMeans(mat)
mat_scaled = t(scale(t(mat)))

# calculate heatmap annotation
tax_heatmap <- tax.clean[rownames(mat_scaled),]
tax_heatmap$CHP <- sapply(rownames(tax_heatmap), function(x) ifelse(x %in% rownames(tab_CHP),"*","ns"))
tax_heatmap$SDT <- sapply(rownames(tax_heatmap), function(x) ifelse(x %in% rownames(tab_SDT),"*","ns"))
tax_heatmap$PYT <- sapply(rownames(tax_heatmap), function(x) ifelse(x %in% rownames(tab_PYT),"*","ns"))
tax_heatmap$FVT <- sapply(rownames(tax_heatmap), function(x) ifelse(x %in% rownames(tab_FVT),"*","ns"))
tax_heatmap$LFD <- sapply(rownames(tax_heatmap), function(x) ifelse(x %in% rownames(tab_LFD),"*","ns"))

index <- match(rownames(tax_heatmap), rownames(tab_CHP_all))
index
tax_heatmap$CHP_p <- tab_CHP_all$padj[index]


index <- match(rownames(tax_heatmap), rownames(tab_SDT_all))
index
tax_heatmap$SDT_p <- tab_SDT_all$padj[index]

index <- match(rownames(tax_heatmap), rownames(tab_PYT_all))
index
tax_heatmap$PYT_p <- tab_PYT_all$padj[index]

index <- match(rownames(tax_heatmap), rownames(tab_FVT_all))
index
tax_heatmap$FVT_p <- tab_FVT_all$padj[index]

index <- match(rownames(tax_heatmap), rownames(tab_LFD_all))
index
tax_heatmap$LFD_p <- tab_LFD_all$padj[index]

tax_heatmap <- tax_heatmap[order(tax_heatmap$CHP_p,tax_heatmap$SDT_p, tax_heatmap$PYT_p, tax_heatmap$FVT_p, tax_heatmap$LFD_p),]

max(c(-log10(tax_heatmap$CHP_p), -log10(tax_heatmap$SDT_p),-log10(tax_heatmap$PYT_p), -log10(tax_heatmap$FVT_p),-log10(tax_heatmap$LFD_p)))
mean(c(-log10(tax_heatmap$CHP_p), -log10(tax_heatmap$SDT_p),-log10(tax_heatmap$PYT_p), -log10(tax_heatmap$FVT_p),-log10(tax_heatmap$LFD_p)))


my_palette <- c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "khaki2", "firebrick")

# adjust tax_heatmap genus and species, pasteurella
tax_heatmap$Family <- as.character(tax_heatmap$Family)
# tax_heatmap$Species <- as.character(tax_heatmap$Species)
tax_heatmap <- tax_heatmap[order(rownames(mat_scaled)),]
common_rows <- intersect(rownames(tax_heatmap), rownames(mat_scaled))
tax_heatmap <- tax_heatmap[common_rows, ]
mat_scaled <- mat_scaled[common_rows, ]

rownames(tax_heatmap)<- paste(rownames(tax_heatmap),tax_heatmap$Family,sep = "_")
rownames(mat_scaled) <-rownames(tax_heatmap)
plot <- mat_scaled

family <- unique(as.character(tax_heatmap$Family))
family_col <- colorRampPalette(my_palette)(length(family))
names(family_col) <- family


pvalue_col_fun = colorRamp2(c(1,0.1,0.05), c("red", "white", "lightseagreen"))
ha_row <- HeatmapAnnotation(ObesevsChP=anno_simple(-log10(tax_heatmap$CHP_p),col = pvalue_col_fun, pch = na_if(tax_heatmap$CHP,"ns"), gp = gpar(fontsize(1))),
                            ObesevsSDT=anno_simple(-log10(tax_heatmap$SDT_p),col = pvalue_col_fun, pch = na_if(tax_heatmap$SDT,"ns"), gp = gpar(fontsize(1))),
                            ObesevsPyT=anno_simple(-log10(tax_heatmap$PYT_p),col = pvalue_col_fun, pch = na_if(tax_heatmap$PYT,"ns"), gp = gpar(fontsize(1))),
                            ObesevsUnmodifiedFVT=anno_simple(-log10(tax_heatmap$FVT_p),col = pvalue_col_fun, pch = na_if(tax_heatmap$FVT,"ns"), gp = gpar(fontsize(1))),
                            ObesevsLFD=anno_simple(-log10(tax_heatmap$LFD_p),col = pvalue_col_fun, pch = na_if(tax_heatmap$LFD,"ns"), gp = gpar(fontsize(1))),
                            Family=anno_simple(tax_heatmap$Family, col = family_col),
                            which = "row")

ha_row_txt <- rowAnnotation(labels = anno_text(rownames(tax_heatmap),which = "row",gp=gpar(fontsize=2, face= "italic")))

ha_col = HeatmapAnnotation(Treatment=meta_order$Treatment,
                           col = list(Treatment=c("FVT_ChP"="#F27970",
                                                  "FVT_SDT"="#BB9727",
                                                  "FVT_PyT"="#54B345",
                                                  "Unmodified_FVT"="#32B897",
                                                  "Obese_control"="#05B9E2",
                                                  "Lean_control"="#8983BF"))
)




Hist <- ComplexHeatmap::pheatmap(plot, 
                                 cluster_cols = FALSE, cluster_rows = TRUE,
                                 name="Z-score", col=colorRamp2(c(-2, 0, 2), c("dodgerblue4", "white","deeppink3")),
                                 top_annotation = ha_col,
                                 left_annotation = ha_row,
                                 right_annotation = ha_row_txt,
                                 show_colnames = FALSE, 
                                 show_rownames = FALSE,
                                 heatmap_legend_param = list(legend_direction = "horizontal")
)
Hist

# define the two legend
lgd_genus = Legend(title = "Family", legend_gp = gpar(fill = family_col),labels = family, ncol = 2)

lgd_sig = Legend(title= " ", pch = "*", type = "points", labels = "p < 0.05")


lgd_pvalue = Legend(title = "p value",
                    col_fun  = pvalue_col_fun,
                    at = c(0, 1, 2),
                    labels = c("1","0.1","0.05"),
                    direction = "horizontal")

p_vir_Deseq_heatmap <-draw(Hist,
                       heatmap_legend_list=list(lgd_genus,lgd_pvalue,lgd_sig),
                       heatmap_legend_side = "bottom", annotation_legend_side = "bottom")



##Simple one in figure4#########################################################################



ha_row <- HeatmapAnnotation(Family=anno_simple(tax_heatmap$Family, col = family_col, gp=gpar(fontsize=4)),
                            which = "row")

ha_col = HeatmapAnnotation(Treatment=meta_order$Treatment,
                           col = list(Treatment=c("FVT_ChP"="#F27970",
                                                  "FVT_SDT"="#BB9727",
                                                  "FVT_PyT"="#54B345",
                                                  "Unmodified_FVT"="#32B897",
                                                  "Obese_control"="#05B9E2",
                                                  "Lean_control"="#8983BF")),
                                      show_legend=FALSE, gp=gpar(fontsize=4)
)



plot <- mat_scaled[rownames(tax_heatmap),]
p_vir_Deseq_heatmap_simple <- ComplexHeatmap::pheatmap(plot, 
                                 cluster_cols = FALSE, cluster_rows = TRUE,
                                 name="Z-score", col=colorRamp2(c(-2, 0, 2), c("dodgerblue4", "white","deeppink3")),
                                 top_annotation = ha_col,
                                 left_annotation = ha_row,
                                 show_colnames = FALSE, 
                                 show_rownames = FALSE,
                                 
                                 heatmap_legend_param = list(legend_direction = "vertical")
)

p_vir_Deseq_heatmap_simple
























