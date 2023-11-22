setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("Analysis0_file_loading_and_prep.R")

ps <- subset_samples(PSB, Sample_time_point %in% "Termination")
ps <- tax_glom(ps, "Species", NArm = FALSE) #select a level to compare
### loop for all the grouping compared with Obese control
Treatments <- unique(sample_data(ps)$Treatment)

Treatments
Treatments <- Treatments[Treatments!="Obese_control"]
path_table <- "stat_result/desep2/"
dir.create(path_table)

for (group in Treatments){
  ps.sub <- prune_samples(sample_data(ps)$Treatment %in% c(group, "Obese_control"), ps)
  ps.sub
  # remove all error taxa
  ps.ds <- phyloseq_to_deseq2(ps.sub, ~Treatment)
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