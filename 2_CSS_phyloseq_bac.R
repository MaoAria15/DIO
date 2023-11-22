###CSS Nomralization###

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#######################Virome#######################

#Load Phyloseq

source("1_Phyloseq_import_bac.R")
OTU_table <- as.data.frame(otu_table(PSB))
taxonomy <- as.data.frame(tax_table(PSB))
mapping <- as.data.frame(sample_data(PSB))
tree<- phyloseq::phy_tree(PSB)
dim(OTU_table)
dim(taxonomy)
dim(mapping)

#CSS normalization 
data.metagenomeSeq = newMRexperiment(OTU_table)
new_count = filterData(data.metagenomeSeq, depth = 1)
new_count[rowSums(new_count)>0,] %>% dim                   #some samples do not have reads after clean up steps.

data.cumnorm <- cumNorm(obj =new_count, p = cumNormStat(new_count))
OTU_read_count_CSS = data.frame(MRcounts(data.cumnorm, norm=TRUE, log=TRUE))
otu_mat = MRcounts(data.cumnorm, norm=TRUE, log=TRUE)
dim(otu_mat) # sample control was removed during normalization

# Filter the mapping file to keep only the samples present in the OTU table
sample_id_otu <- colnames(otu_mat)
map_mat <- mapping[rownames(mapping) %in% sample_id_otu, ]
dim(map_mat)




##Construct individual tables is phyloseq format
otu_mat <- as.matrix(otu_mat)
taxonomy <- as.matrix(taxonomy)
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)
TAX <- phyloseq::tax_table(taxonomy)
samples <- sample_data(map_mat)

##Combine in phyloseq object

PSB.CSS <- phyloseq::phyloseq(OTU, TAX, samples,tree)

#Clean up

rm(data.cumnorm, data.metagenomeSeq, OTU_table, taxonomy, rownames, map_mat, otu_mat, samples, OTU, TAX, OTU_read_count_CSS, mapping, new_count )

