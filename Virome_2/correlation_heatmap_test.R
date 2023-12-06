setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))
source("setup.R")


# Get bacteria data prepared###############################################################################
ps <- readRDS("data/corr_data/ps_bacterium.rds")
ps<- subset_samples(ps, Sample_time_point %in% c("Termination"))

#merge OTU having same genus
ps <- phyloseq::tax_glom(ps, "Species")
ps <- phyloseq::prune_taxa(!(rownames(otu_table(ps)) %in% "OTU_656"), ps)# Remove the specific OTU by name for non-unique value when setting 'row.names': ??uncultured_rumen??

# merge_ps <- merge_samples(ps,"Group")
merge_ps <- transform_sample_counts(ps, function(x) x / sum(x))



# separate dataset
otu_mat <- as.data.frame(phyloseq::otu_table(merge_ps))
taxonomy <- as.data.frame(phyloseq::tax_table(merge_ps))
mapping <- as.data.frame(phyloseq::sample_data(merge_ps))

#Prepare data set for next step
bac_otu <- otu_mat
genus_vec <- taxonomy[,7]    # for species

rownames(bac_otu) <- genus_vec
mouse_id <- mapping$Mouse_ID
colnames(bac_otu) <- mouse_id


#Getting the top 30 genus
rown<-rownames(bac_otu)
bac_otu$tax <- rown

bac_otu <- bac_otu %>% 
  select(-tax)
bac_otu$mean <- rowMeans(bac_otu, na.rm = TRUE)
bac_otu<-bac_otu %>%
  filter(mean > 0.01)%>%
  select(-mean)


bac_otu <-t(bac_otu)

# rm(merge_ps,otu_mat,ps,taxonomy,top,t_bac_otu,mapping)

# Get virome data prepared#########################################################################################################
ps <- readRDS("data/corr_data/ps_virome.rds")
ps<- subset_samples(ps, Sample_time_point %in% c("Termination"))
# merge vOTU having same family
ps <- tax_glom(ps, "Genus")
# merge_ps <- merge_samples(ps,"Group")
merge_ps <- transform_sample_counts(ps, function(x) x / sum(x))
# separate dataset
mapping_vir <- as.data.frame(sample_data(merge_ps))
otu_mat <- as.data.frame(otu_table(merge_ps))
taxonomy <- as.data.frame(tax_table(merge_ps))
#Prepare data set for next step
vir_otu <- otu_mat
fam_vec <- taxonomy[,7]
vir_otu$tax <- fam_vec
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




corr1 <- vir_otu
corr2 <- bac_otu



common_rows <- intersect(rownames(corr1), rownames(corr2))# Identify common column names
num_sample <- length(common_rows)#Get the number of sample compared

# Subset the data frames to keep only the common columns
corr1 <- corr1[common_rows, ]
corr2 <- corr2[common_rows, ]

#Correlation test

res <- corr.test(x=corr2,y=corr1, use = "pairwise", method = "spearman", adjust = "fdr", alpha= 0.5, ci= FALSE, minlength = 5)
r<- as.data.frame(res$r)
r <- r[!apply(is.na(r) | r == "", 1, all), ]# Remove rows that only have NA values
r <- r[, colSums(is.na(r)) != nrow(r)]# Remove columns that only have NA values
write.csv(r,"Virome_2/stat_result/vir_bac_spearman_corr.csv")
p_vir_corr_heatmap_sim <- Heatmap(r,
                                    name = "Correlation",
                                    column_title = "Spearman correlation between 16s rRNA and virome",
                                    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                                    row_names_gp = gpar(fontsize = 10, fontface = "italic"),
                                    column_names_gp = gpar(fontsize = 10, fontface = "italic"))

