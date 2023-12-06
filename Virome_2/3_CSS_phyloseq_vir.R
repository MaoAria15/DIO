##CSS Nomralization###

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#######################Virome#######################

#Load Phyloseq : make sure if we want to have decontam ps, if yes choose source2, no choose source1

# source("1_Phyloseq_import_vir.R")

source("2_Decontamination.R")

#########################PSV table without Realm - for Ampvis2

PSV.no.realm <- subset_taxa(PSV, select = -c(Realm,Domain)) 

#######################Remove low read samples#####################

PSV.no.realm = prune_samples(sample_sums(PSV) >= 15000, PSV.no.realm)


# Preparing PSV.no.Realm.CSS####################################################################
ps<- PSV.no.realm
OTU_table <- as.data.frame(otu_table(ps))
taxonomy <- as.data.frame(tax_table(ps))
mapping <- as.data.frame(sample_data(ps))
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

PSV.no.Realm.CSS <- phyloseq::phyloseq(OTU, TAX, samples)


#Clean up

rm(data.cumnorm, data.metagenomeSeq, OTU_table, taxonomy, 
    map_mat, otu_mat, samples, OTU, TAX, OTU_read_count_CSS, 
   mapping, new_count,ps,sample_id_otu,PSV )

