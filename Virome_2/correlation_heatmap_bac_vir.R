setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))
source("setup.R")




# Get bacteria data prepared###############################################################################
ps <- readRDS("data/corr_data/ps_bacterium.rds")
ps<- subset_samples(ps, Sample_time_point %in% c("Termination"))



#merge OTU having same genus

merge_ps <- transform_sample_counts(ps, function(x) x / sum(x))



# separate dataset
otu_mat <- as.data.frame(phyloseq::otu_table(merge_ps))
taxonomy <- as.data.frame(phyloseq::tax_table(merge_ps))
mapping <- as.data.frame(phyloseq::sample_data(merge_ps))

#Prepare data set for next step
bac_otu <- otu_mat
rownames(bac_otu) <- paste(rownames(bac_otu), taxonomy$Species, sep ="_")
mouse_id <- mapping$Mouse_ID
colnames(bac_otu) <- mouse_id


#Getting the top 30 genus
rown<-rownames(bac_otu)
bac_otu$tax <- rown

bac_otu <- bac_otu %>% 
  select(-tax)
bac_otu$mean <- rowMeans(bac_otu, na.rm = TRUE)
bac_otu<-bac_otu %>%
  filter(mean > 0.001)%>%
  select(-mean)


bac_otu <-t(bac_otu)


# rm(merge_ps,otu_mat,ps,taxonomy,top,t_bac_otu,mapping)


# Get virome data prepared#########################################################################################################
ps <- readRDS("data/corr_data/ps_virome.rds")
ps<- subset_samples(ps, Sample_time_point %in% c("Termination"))



merge_ps <- transform_sample_counts(ps, function(x) x / sum(x))

# separate dataset

mapping_vir <- as.data.frame(sample_data(merge_ps))
otu_mat <- as.data.frame(otu_table(merge_ps))
taxonomy <- as.data.frame(tax_table(merge_ps))
#Prepare data set for next step
vir_otu <- otu_mat
vir_otu$tax <- paste(rownames(vir_otu),taxonomy$Genus, sep = "_")
# Merge duplicated tax names and calculate the sum of corresponding values
vir_otu <- aggregate(. ~ tax, data = vir_otu, FUN = sum)
tax <- vir_otu$tax
rownames(vir_otu) <- tax
vir_otu <- vir_otu %>% select(-tax)
mouse_id <- mapping_vir$Mouse_ID
colnames(vir_otu) <- mouse_id


vir_otu$mean <- rowMeans(vir_otu, na.rm = TRUE)
vir_otu<-vir_otu %>%
  filter(mean > 0.001)%>%
  select(-mean)


vir_otu <-t(vir_otu)
#only keep top 20 bac genus in bac_otu

rm(taxonomy,mapping, merge_ps,otu_mat,ps,merged_t_vir_otu,t_vir_otu,top,mapping_vir,class,fam_vec,family,genus_vec,i,j,kingdom,order,phylum,pig_id,pig_vir_id,rown,sample_id_otu,singletones,tax)


corr1 <- bac_otu
corr2 <- vir_otu

# Identify common column names
common_rows <- intersect(rownames(corr1), rownames(corr2))
#Get the number of sample compared
num_sample <- length(common_rows)

# Subset the data frames to keep only the common columns
corr1 <- corr1[common_rows, ]
corr2 <- corr2[common_rows, ]


#Correlation test

res <- corr.test(x=corr2,y=corr1, use = "pairwise", method = "spearman", adjust = "fdr", alpha= 0.5, ci= FALSE, minlength = 5)
vir_bac <- cor_tab(corr2,corr1)
write.csv(cor.data,"Virome_2/stat_result/bac_vir_spearman_corr.csv")
r<- as.data.frame(res$r)
# Remove rows that only have NA values
r <- r[!apply(is.na(r) | r == "", 1, all), ]
# Remove columns that only have NA values
r <- r[, colSums(is.na(r)) != nrow(r)]

# write.csv(r,"Virome_2/stat_result/vir_bac_spearman_corr_high_resolution.csv")
#make heatamp
p_vir_corr_heatmap <- Heatmap(r,
                                name = "Correlation",
                                column_title = "Spearman correlation between 16s rRNA and virome",
                                column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                                row_names_gp = gpar(fontsize = 6, fontface = "italic"),
                                column_names_gp = gpar(fontsize = 6, fontface = "italic"))

