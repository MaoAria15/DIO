###CSS Nomralization###

#Use after running perl script (run perl on TPM normalized table for virome analysis)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(metagenomeSeq)
library(dplyr)
library(stringr)
library(phyloseq)
library(reshape2)


#######################Virome#######################

#Load OTU_table

OTU_table <- read.delim("data/OTU_UMI16S_192_PrePhage1.txt",stringsAsFactors = F,header= T)

#######Do CSS norm
dim(OTU_table)
rownames <- OTU_table[,1]

taxonomy <- OTU_table["Taxonomy"]
#Remove otu names
zOTU <- OTU_table[-1]
#Remove taxonomy
zOTU <- within(zOTU, rm("Taxonomy"))
dim(zOTU)
class(zOTU)
OTU_read_count = as.data.frame(zOTU)


data.metagenomeSeq = newMRexperiment(zOTU)
                                     
data.cumnorm = cumNorm(data.metagenomeSeq, p=cumNormStatFast(data.metagenomeSeq))
OTU_read_count_CSS = data.frame(MRcounts(data.cumnorm, norm=TRUE, log=TRUE))
#data.cumnorm
otu_mat = MRcounts(data.cumnorm, norm=TRUE, log=TRUE)

row.names(otu_mat) <- rownames
#Add taxonomy back
otu_mat <- cbind(otu_mat,taxonomy)

#########Phyloseq import

#generate taxonomy table from tax table
tax_mat <- read.delim('data/taxonomy.txt', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)


# Remove "k__", "p__" etc from taxonomy column
tax_mat[,1] <- str_remove_all(tax_mat[,1], str_c(c("k__", "p__", "c__", "o__", "f__", "g__", "s__"), collapse="|"))
tax <- tax_mat %>%
  dplyr::select(Taxonomy) %>%
  separate(Taxonomy, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), ";")

tax.clean <- data.frame(row.names = row.names(tax),
                        Kingdom = str_replace(tax[,1], "k__",""),
                        Phylum = str_replace(tax[,2], "p__",""),
                        Class = str_replace(tax[,3], "c__",""),
                        Order = str_replace(tax[,4], "o__",""),
                        Family = str_replace(tax[,5], "f__",""),
                        Genus = str_replace(tax[,6], "g__",""),
                        Species = str_replace(tax[,7], "s__",""),
                        stringsAsFactors = FALSE)

tax.clean <- subset(tax.clean, Kingdom != "Unassigned")
tax.clean[is.na(tax.clean)] <- ""

for (i in 1:nrow(tax.clean)){
  for(j in 1:7) {
    if (tax.clean[i,j] == "uncultured_bacterium") {
      tax.clean[i,j] <- "";
    }
  }
}

for (i in 1:nrow(tax.clean)){
  for(j in 1:7) {
    if (tax.clean[i,j] == "uncultured") {
      tax.clean[i,j] <- "";
    }
  }
}

tax.clean[tax.clean=="__"] <- ""

for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,7] != ""){
    tax.clean$Species[i] <- paste(tax.clean$Genus[i], tax.clean$Species[i], sep = "_")
  } else if (tax.clean[i,2] == ""){
    kingdom <- paste("Unclassified", tax.clean[i,1], sep = "_")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("Unclassified", tax.clean[i,2], sep = "_")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("Unclassified", tax.clean[i,3], sep = "_")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("Unclassified", tax.clean[i,4], sep = "_")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste("Unclassified", tax.clean[i,5], sep = "_")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste("Unclassified ",tax.clean$Genus[i], sep = "_")
  }
}

#Convert to dataframe

tax_mat <- as.matrix.data.frame(tax.clean)

#Remove taxonomy from otu table
otu_mat <- within(otu_mat,rm("Taxonomy"))

otu_mat <- as.matrix(otu_mat)
map_mat <- read.delim('data/Mapping.txt', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)


##Construct individual tables is phyloseq format

OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)
TAX <- tax_table(tax_mat)
samples <- sample_data(map_mat)

##Combine in phyloseq object

PSB.CSS <- phyloseq(OTU, TAX, samples)

#Clean up

rm(data.cumnorm, data.metagenomeSeq, OTU_table, taxonomy, zOTU, p, rownames, tax_mat, map_mat, otu_mat, samples, split, OTU, TAX)

