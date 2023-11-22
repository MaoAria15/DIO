
setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))
source("setup.R")

# Get clinical data prepared###################################################################################################
clinical_data <- Data

rownames(clinical_data) <- clinical_data$Mouse_ID
clinical_data <- clinical_data[,!(colnames(clinical_data) %in% c("Mouse_ID","Group","Liver_response",
                                                                 "Steatosis","Inflammation","Hepatocyte_injury",
                                                                 "OGTT_t0_w13","OGTT_t15_w13","OGTT_t30_w13","OGTT_t60_w13",
                                                                 "OGTT_t90_w13","OGTT_t120_w13",
                                                                 "OGTT_t0_w18","OGTT_t15_w18","OGTT_t30_w18","OGTT_t60_w18",
                                                                 "OGTT_t90_w18","OGTT_t120_w18"))]

# clinical_data <- t(clinical_data)


clinical_data <- clinical_data[,(colnames(clinical_data) %in% c("fat_Dendritic_cells","fat_Macrophages","fat_M1_macrophages",
                                                                "fat_B_cells","fat_Cytotoxic_T_cells","fat_Th_cells","fat_Central_memory_T_cells_CD8",
                                                                "fat_Effector_memory_T_cells_CD8","fat_Naive_T_cells_CD8","fat_Central_memory_T_cells_CD4",
                                                                "fat_Effector_memory_T_cells_CD4","fat_Naive_T_cells_CD4" ))]
# clinical_data <- clinical_data[,(colnames(clinical_data) %in% c("mln_Dendritic_cells","mln_Macrophages","mln_M1_macrophages",         
#                                                                 "mln_B_cells","mln_Cytotoxic_T_cells","mln_Th_cells","mln_Central_memory_T_cells_CD8",
#                                                                 "mln_Effector_memory_T_cells_CD8","mln_Naive_T_cells_CD8","mln_Central_memory_T_cells_CD4",
#                                                                 "mln_Effector_memory_T_cells_CD4","mln_Naive_T_cells_CD4" ))]

clinical_data <- as.data.frame(clinical_data)

rm(merge_ps,mapping,mapping_vir,otu_mat,ps,S.PSB.Colon,taxonomy,t_clinical_data)





# Get bacteria data prepared###############################################################################
ps <- readRDS("data/corr_data/ps_bacterium.rds")
ps<- subset_samples(ps, Sample_time_point %in% c("Termination"))



#merge OTU having same genus


# ps <- tax_glom(ps, "Genus")
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
# genus_vec <- taxonomy[,6]
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








#Find common pig 

common_row_names <- intersect(rownames(bac_otu), rownames(clinical_data))
# Keep only rows with common row names in all three dataframes
bac_otu <- subset(bac_otu, rownames(bac_otu) %in% common_row_names)
clinical_data <- subset(clinical_data, rownames(clinical_data) %in% common_row_names)


# data <- bind_rows(data_vir_bac,t_clinical_data)

#PreParing type form
node_names <- rownames(t(bac_otu))
type_bac <- data.frame(node=node_names,type="Bacteria")


node_names <- rownames(t(clinical_data))
type_clinical <- data.frame(node=node_names,type="Clinical_result")

type <- bind_rows(type_bac,type_clinical)



### function to calculate respective scaled spearman correlation
cor_tab <- function(dataframe1,dataframe2){
  cor <-corr.test(dataframe1,dataframe2, use = "pairwise",method="spearman",adjust="fdr", alpha=.05,ci=FALSE)
  cor.r <- data.frame(cor$r) # æå–Rå€?
  cor.p <- data.frame(cor$p) # æå–på€?
  cor.r$from <- rownames(cor.r)
  cor.p$from <- rownames(cor.p)
  p <- cor.p %>%
    gather(key = "to", value = "p", -from) %>%
    data.frame()
  cor.data <- cor.r %>%
    gather(key = "to",value = "r", -from) %>%
    data.frame() %>%
    left_join(p, by=c("from","to")) %>%
    filter(p<=0.05, from !=to) %>%
    mutate(
      linecolor=ifelse(r>0, "positive","negative"),  #set the color and line type
      linesize=abs(r)  #set line wildth
    )
  return(cor.data)
}


#Calculate spearman correlation among all clinial and microbiome data

# t_data_filtered <-t(data_filtered)
bac_clinical <- cor_tab(bac_otu,clinical_data)


#merge correlationdata together
cor.data <- rbind(bac_clinical)



#Construct NET
#calculate the degree of each node


c(as.character(cor.data$from),as.character(cor.data$to)) %>%
  dplyr::as_tibble() %>%
  dplyr::group_by(value) %>%
  dplyr::summarize(n=n()) -> vertices
colnames(vertices) <- c("node", "n")  



#### add type character
vertices %>%
  select(-n) %>% # delete n because it is not accurate
  left_join(type,by="node") -> vertices 

# The nodes in the network diagram will be drawn sequentially according to 
# the order of the node attribute files. In order to make the positions of 
# variables of the same type close, the nodes are sorted according to the node attributes.
vertices$type = factor(vertices$type,levels = unique(vertices$type))
vertices <- arrange(vertices,type)
write.csv(vertices,"clinical_microbiome_vertices.csv")
head(vertices)



colnames(cor.data)[3] <- "Coefficient" 
# construct graph net structure data
g <- graph_from_data_frame(cor.data, vertices = vertices, directed = FALSE )
g
vcount(g) # number of nodesï¼?18
ecount(g) # number of link: 16
#get.vertex.attribute(g) # View the node properties contained in the graph
#get.edge.attribute(g) # View the link properties contained in the graph


### making simple graph net
is.simple(g) # For non-simple graphs, the number of links will be high, so it needs to be converted to a simple graph.

g



p_bac_fats_node <- ggraph(g, layout = 'linear', circular=TRUE) +
  geom_edge_link(aes(edge_alpha = abs(Coefficient), edge_width = abs(Coefficient), color = Coefficient)) +
  guides(edge_alpha = "none", edge_width = "none") +
  scale_edge_colour_gradientn(limits = c(-1, 1), colors = c("#f1a340","#998ec3")) +
  geom_node_point(color = "#43a2ca", size = 5) +
  geom_node_text(aes(label = name), repel = TRUE) +
  theme_graph(base_family="sans") 



# write.csv(cor.data,"Result-16s/Simone/Updated correlation 13102023/spearman-correlation(clinical-16s-virome).csv")





