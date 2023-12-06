
setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))
source("setup.R")

#OTU table
host_mat <- read.csv('data/Virome_data_2/Host_prediction.csv')

otu <- read.delim('data/Virome_data_2/tpm_vOTU_table.tsv', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)
# Group the data by the 'out' column and filter for rows with the highest 'score'
host_mat <- host_mat %>%
  group_by(Virus) %>%
  filter(Confidence_score == max(Confidence_score))


host_mat$Confidence_score1 <- ifelse(is.na(host_mat$Confidence_score1), "", host_mat$Confidence_score1)

host_mat <- host_mat %>%
  group_by(Virus) %>%
  filter(Confidence_score1 == max(Confidence_score1))

host_mat$Confidence_score2 <- ifelse(is.na(host_mat$Confidence_score2), "", host_mat$Confidence_score2)
host_mat <- host_mat %>%
  group_by(Virus) %>%
  filter(Confidence_score2 == max(Confidence_score2))
host_mat <- host_mat %>%
  distinct(Virus, .keep_all = TRUE)

# Assuming your data frame is called 'host_mat'
duplicates <- host_mat$Virus[duplicated(host_mat$Virus)]

taxonomy <- subset(host_mat[,c("Virus","Taxonomy")])

otu_classified <- intersect(rownames(otu),taxonomy$Virus)
otu_all <- rownames(otu)

missing_votus <- otu_all[!otu_all %in% otu_classified]
Unknown_tax <- data.frame(Virus= missing_votus,
                          Taxonomy = rep("Unknown;;;;;", length(missing_votus)))
taxonomy <- rbind(Unknown_tax,taxonomy)

rownames <- taxonomy$Virus

taxonomy <- subset(taxonomy, select = -Virus)
rownames(taxonomy) <- rownames

taxonomy <- taxonomy %>%
  dplyr::select(Taxonomy) %>%
  separate(Taxonomy, c("Domain","Phylum","Class","Order","Family","Genus"), ";")

taxonomy[] <- lapply(taxonomy, function(x) gsub("(f_|d_|c_|o_|p_|g_)", "", x))

for (i in 1:nrow(taxonomy)){
  if (taxonomy[i,6] != ""){
    taxonomy$Genus[i] <- taxonomy$Genus[i]
  }  else if (taxonomy[i,2] == ""){
    domain <- paste("Unclassified", taxonomy[i,1], sep = "_")
    taxonomy[i, 2:6] <- domain
  } else if (taxonomy[i,3] == ""){
    phylum <- paste("Unclassified", taxonomy[i,2], sep = "_")
    taxonomy[i, 3:6] <- phylum
  } else if (taxonomy[i,4] == ""){
    class <- paste("Unclassified", taxonomy[i,3], sep = "_")
    taxonomy[i, 4:6] <- class
  } else if (taxonomy[i,5] == ""){
    order <- paste("Unclassified", taxonomy[i,4], sep = "_")
    taxonomy[i, 5:6] <- order
  } else if (taxonomy[i,6] == ""){
    taxonomy$Genus[i] <- paste("Unclassified ",taxonomy$Family[i], sep = "_")
  }
}
taxonomy[] <- lapply(taxonomy, function(x) gsub("Unclassified_Unknown", "Unknown", x))


write.csv(taxonomy,"data/Virome_data_2/host_taxonomy.csv")

