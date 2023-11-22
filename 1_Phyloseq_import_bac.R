

setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))
source("setup.R")
######################################Construct Phyloseq objects##########################################


###########Virome

##Load raw data 

#OTU table
otu_mat <- read.delim('data/16S_data/otu_result.tsv', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)

#04/11-2021 taxonomy from Josue
tax_mat <- read.delim('data/16S_data/greengenes_taxa.txt', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)


tax <- tax_mat %>%
  dplyr::select(taxonomy) %>%
  separate(taxonomy, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), ";")

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
tax.clean[tax.clean=="__"] <- ""

for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,7] != ""){
    tax.clean$Species[i] <- paste(tax.clean$Genus[i], tax.clean$Species[i], sep = " ")
  } else if (tax.clean[i,2] == ""){
    kingdom <- paste("k", tax.clean[i,1], sep = "_")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("p", tax.clean[i,2], sep = "_")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("c", tax.clean[i,3], sep = "_")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("o", tax.clean[i,4], sep = "_")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste("f", tax.clean[i,5], sep = "_")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste("g ",tax.clean$Genus[i], sep = "_")
  }
}

# 
tax_mat <- as.matrix.data.frame(tax.clean)

#Remove taxonomy from otu table
otu_mat <- within(otu_mat,rm("taxonomy"))
#Round all values
rownames <- rownames(otu_mat)
otu_mat <- mutate_all(otu_mat, round) %>% mutate_all(as.numeric)
rownames(otu_mat) <- rownames

map_mat <- read.delim('data/16S_data/DIO_mapping.txt', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)

# map_mat <- dplyr::mutate_at(map_mat,1:9,funs(factor))
# map_mat <- dplyr::mutate_if(map_mat, is.character, as.numeric)
TREE = read_tree("data/16S_data/zOTU.tree")
##Construct individual tables is phyloseq format

OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)
TAX <- tax_table(tax_mat)
samples <- sample_data(map_mat)
##Combine in phyloseq object

PSB <- phyloseq(OTU, TAX, samples,TREE)

#Make metadata files for easy overview

meta.bac <- as.data.frame(as.matrix(sample_data(PSB)))

#############Remove left over objects

rm(tax_mat, map_mat, otu_mat, samples, OTU, TAX,tax,tax.clean,meta.bac,TREE)



