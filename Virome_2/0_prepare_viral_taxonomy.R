
setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))
source("setup.R")

#OTU table
otu_mat <- read.csv('data/Virome_data_2/viral_taxonomy.csv')

otu <- read.delim('data/Virome_data_2/tpm_vOTU_table.tsv', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)

otu_classified <- intersect(rownames(otu),otu_mat$Votu)
otu_all <- rownames(otu)
missing_votus <- otu_all[!otu_all %in% otu_classified]
Unknown_tax <- data.frame(Votu= missing_votus,
                          Taxonomy = rep("Unknown;;;;;;;;", length(missing_votus)))
taxonomy <- rbind(Unknown_tax,otu_mat)
rownames <- taxonomy$Votu
taxonomy <- subset(taxonomy, select = -Votu)
rownames(taxonomy) <- rownames


tax_clean <- taxonomy %>%
  dplyr::select(Taxonomy) %>%
  separate(Taxonomy, c("Domain","Realm","Kingdom","Phylum","Class","Order","Family","Subfamily","Genus"), ";")
tax_clean[] <- lapply(tax_clean, function(x) gsub("Unclassified", "", x))


for (i in 1:nrow(tax_clean)){
  if (tax_clean[i,9] != ""){
    tax_clean$Genus[i] <- paste("",tax_clean$Genus[i], sep = "")
  }  else if (tax_clean[i,2] == ""){
    domain <- paste("Unclassified", tax_clean[i,1], sep = "_")
    tax_clean[i, 2:9] <- domain
  } else if (tax_clean[i,3] == ""){
    realm <- paste("Unclassified", tax_clean[i,2], sep = "_")
    tax_clean[i, 3:9] <- realm
  } else if (tax_clean[i,4] == ""){
    phylum <- paste("Unclassified", tax_clean[i,3], sep = "_")
    tax_clean[i, 4:9] <- phylum
  } else if (tax_clean[i,5] == ""){
    class <- paste("Unclassified", tax_clean[i,4], sep = "_")
    tax_clean[i, 5:9] <- class
  } else if (tax_clean[i,6] == ""){
    order <- paste("Unclassified", tax_clean[i,5], sep = "_")
    tax_clean[i, 6:9] <- order
  } else if (tax_clean[i,7] == ""){
    family <- paste("Unclassified", tax_clean[i,6], sep = "_")
    tax_clean[i, 7:9] <- family
  } else if (tax_clean[i,8] == ""){
    subfamily <- paste("Unclassified", tax_clean[i,7], sep = "_")
    tax_clean[i, 8:9] <- subfamily
  } else if (tax_clean[i,9] == ""){
    tax_clean$Genus[i] <- paste("Unclassified ",tax_clean$Subfamily[i], sep = "_")
  }
}
tax_clean[] <- lapply(tax_clean, function(x) gsub("Unclassified_Unknown", "Unknown", x))


for (i in 1:nrow(tax_clean)){
  if (tax_clean[i,8] != ""){
    tax_clean$Subfamily[i] <- paste("",tax_clean$Subfamily[i], sep = "")
  }   else if (tax_clean[i,5] == ""){
    class <- paste("Unclassified", tax_clean[i,4], sep = "_")
    tax_clean[i, 5:9] <- class
  } else if (tax_clean[i,6] == ""){
    order <- paste("Unclassified", tax_clean[i,5], sep = "_")
    tax_clean[i, 6:9] <- order
  } else if (tax_clean[i,7] == ""){
    family <- paste("Unclassified", tax_clean[i,6], sep = "_")
    tax_clean[i, 7:9] <- family
  } else if (tax_clean[i,9] == ""){
    tax_clean$Subfamily[i] <- paste("Unclassified ",tax_clean$Family[i], sep = "_")
  }
}

for (i in 1:nrow(tax_clean)){
  if (tax_clean[i,7] != ""){
    tax_clean$Family[i] <- paste("",tax_clean$Family[i], sep = "")
  }   else if (tax_clean[i,5] == ""){
    class <- paste("Unclassified", tax_clean[i,4], sep = "_")
    tax_clean[i, 5:9] <- class
  } else if (tax_clean[i,6] == ""){
    order <- paste("Unclassified", tax_clean[i,5], sep = "_")
    tax_clean[i, 6:9] <- order
  } else if (tax_clean[i,9] == ""){
    tax_clean$Family[i] <- paste("Unclassified ",tax_clean$Order[i], sep = "_")
  }
}


for (i in 1:nrow(tax_clean)){
  if (tax_clean[i,6] != ""){
    tax_clean$Order[i] <- paste("",tax_clean$Order[i], sep = "")
  }   else if (tax_clean[i,5] == ""){
    class <- paste("Unclassified", tax_clean[i,4], sep = "_")
    tax_clean[i, 5:9] <- class
  } else if (tax_clean[i,9] == ""){
    tax_clean$Order[i] <- paste("Unclassified ",tax_clean$Class[i], sep = "_")
  }
}


# Still few level need to be maunally fit

write.csv(tax_clean,"data/Virome_data_2/Viral_taxonomy.csv")

