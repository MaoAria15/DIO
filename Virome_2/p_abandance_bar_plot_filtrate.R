setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("4_Analysis0_file_loading_and_prep.R")

ps <- PSV.Filtrate
vir.phyl <- tax_glom(ps, "Family", NArm = FALSE)

ps0 <- transform_sample_counts(vir.phyl, function(x) x *100/ sum(x))

ps1 <- merge_samples(ps0, "Treatment")

ps1 <- transform_sample_counts(ps1, function(x) x / sum(x))


#Create melted dataframe
df <- psmelt(ps1)

#Select last non-empty taxonomic rank
df[df==""] <- NA

df$tax <- apply(df, 1, function(x) tail(na.omit(x), 1))


top<-df %>%
  # filter(level=="Phylum")%>%
  group_by(tax)%>%
  dplyr::summarize(mean_abund=mean(Abundance), .groups = "drop")%>%
  arrange(desc(mean_abund))

#Find top 30 tax
top

top30 <- top$tax[1:18] #for species
# top30 <- top$tax[1:5] #for phylum
df0 <- df %>%
  mutate(tax = fct_other(df$tax, keep=c(as.matrix(top30))))%>%
  arrange(desc(tax))

#Set order for samples
df0$Sample <- factor(df0$Sample, levels = groups[1:4])


df0$Abundance <- df0$Abundance*100
top30 <-c(top30,"Other")
df0$tax <- factor(df0$tax,top30)
#join all qualitative palettes
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

p_abundance_bar_termination <- ggplot(df0, aes(Sample, Abundance, fill = fct_reorder(tax, Abundance, .fun = mean))) + 
  geom_col(width = 0.8) +
  theme_classic() +
  scale_fill_jco(name = "Taxonomy") +
  scale_fill_manual(values=rep(col_vector,10),name="Family")+
  mytheme_beta+
  theme(legend.text = element_text(face = "italic"))+
  scale_x_discrete(breaks=groups[1:4],
                   labels=c(glue("FVT\nChP"),
                            glue("FVT\nSDT"),
                            glue("FVT\nPyT"),
                            glue("Unmodified\nFVT")
                   )
  ) +
  labs(x= "",y="Mean Relative abundance (%)", title = "Termination") 

p_abundance_bar_termination



