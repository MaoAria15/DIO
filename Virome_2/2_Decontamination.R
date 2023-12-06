
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##############################Load sequencing data#######################################
##Import files

source("1_Phyloseq_import_vir.R")



#Remove the contaminants from phyloseqs
#########For PSV##################################################
ps <- PSV
ps_control <- subset_samples(ps, Sample_or_Control %in% "Neg_control")
ps_sample<- subset_samples(ps, Sample_or_Control %in% "TRUE_sample")
otu_control <- as.data.frame(otu_table(ps_control))
otu_sample<- as.data.frame(otu_table(ps_sample))
n_control <- ncol(otu_control)
n_sample <- ncol(otu_sample)
tax<- as.data.frame(phyloseq::tax_table(ps))
# Combine the columns into a single "Taxonomy" column
tax$Taxonomy <- paste(
  tax$Domain,
  tax$Realm,
  tax$Kingdom,
  tax$Phylum,
  tax$Class,
  tax$Order,
  tax$Family,
  tax$Genus,
  sep = "; "
)
decon_table <- data.frame(row.names = NULL, OTU_ID = rownames(otu_control) , otu_control,otu_sample,Taxonomy=tax$Taxonomy)

result <- decon(decon_table,numb.blanks = n_control,numb.ind = n_sample,taxa = T, runs = 1, regression = 1)
result_decon <- as.data.frame(result$decon.table)
otu<- result_decon[,colnames(result_decon)!=c("OTU_ID","Mean.blank","Taxonomy")]
otu<- otu[,colnames(otu)!=c("Taxonomy")]
rownames(otu) <- result_decon$OTU_ID
metadata<- sample_data(ps)
metadata<-metadata[intersect(colnames(otu),rownames(metadata)),]
tax <- phyloseq::tax_table(ps)
tax<- tax[intersect(rownames(otu),rownames(tax)),]
otu <- otu_table(otu,taxa_are_rows=T)
tax <- phyloseq::tax_table(tax)
metadata <- sample_data(metadata)
PSV <- phyloseq(otu,tax,metadata)



#########For PSV.HOST##################################################
ps <- PSV.HOST
ps_control <- subset_samples(ps, Sample_or_Control %in% "Neg_control")
ps_sample<- subset_samples(ps, Sample_or_Control %in% "TRUE_sample")
otu_control <- as.data.frame(otu_table(ps_control))
otu_sample<- as.data.frame(otu_table(ps_sample))
n_control <- ncol(otu_control)
n_sample <- ncol(otu_sample)
tax<- as.data.frame(phyloseq::tax_table(ps))
# Combine the columns into a single "Taxonomy" column
tax$Taxonomy <- paste(
  tax$Kingdom,
  tax$Phylum,
  tax$Class,
  tax$Order,
  tax$Family,
  tax$Genus,
  tax$Species,
  sep = "; "
)
decon_table <- data.frame(row.names = NULL, OTU_ID = rownames(otu_control) , otu_control,otu_sample,Taxonomy=tax$Taxonomy)

result <- decon(decon_table,numb.blanks = n_control,numb.ind = n_sample,taxa = T, runs = 1, regression = 1)
result_decon <- as.data.frame(result$decon.table)
otu<- result_decon[,colnames(result_decon)!=c("OTU_ID","Mean.blank","Taxonomy")]
otu<- otu[,colnames(otu)!=c("Taxonomy")]
rownames(otu) <- result_decon$OTU_ID
metadata<- sample_data(ps)
metadata<-metadata[intersect(colnames(otu),rownames(metadata)),]
tax <- phyloseq::tax_table(ps)
tax<- tax[intersect(rownames(otu),rownames(tax)),]
otu <- otu_table(otu,taxa_are_rows=T)
tax <- phyloseq::tax_table(tax)
metadata <- sample_data(metadata)
PSV.HOST <- phyloseq(otu,tax,metadata)

##############################################################
rm(ps,ps_control,result_decon,result,decon_table,metadata,otu_control,otu_sample,ps_sample,otu,tax,n_control,n_sample)
