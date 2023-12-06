setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("Analysis0_file_loading_and_prep.R")

ps <- subset_samples(PSB, Sample_time_point %in% c("Termination"))
ps.rel <- transform_sample_counts(ps, function(x) x / sum(x))


# Firmicutes
ps0 <- tax_glom(subset_taxa(ps.rel, Phylum=="Firmicutes"), "Phylum", NArm = FALSE)
#Select desired rank
df_Firmicutes <- psmelt(ps0)
colnames(df_Firmicutes)[colnames(df_Firmicutes) == "Abundance"] <- "Firmicutes"


# Bacteroidetes
ps0 <- tax_glom(subset_taxa(ps.rel, Phylum=="Bacteroidetes"), "Phylum", NArm = FALSE)
#Select desired rank
df_Bacteroidetes <- psmelt(ps0)
colnames(df_Bacteroidetes)[colnames(df_Bacteroidetes) == "Abundance"] <- "Bacteroidetes"

#merge 2 
merged_df <- merge(df_Bacteroidetes, df_Firmicutes, by = "ProjectID", all = TRUE)%>%
  select(ProjectID, Bacteroidetes, Firmicutes,Treatment.x,Sample.x)
merged_df$ratio <- merged_df$Bacteroidetes/merged_df$Firmicutes

colnames(merged_df)[colnames(merged_df) == "Treatment.x"] <- "Treatment"

# kt <- kruskal.test(data = merged_df,ratio~Treatment)


stat.test <- merged_df %>%
  wilcox_test(ratio~Treatment,
              comparisons = Obse_compare,
              p.adjust.method = "fdr",
              paired = FALSE,
              alternative = "two.sided",
              detailed = TRUE) 
stat.test
stat_neg<- merged_df %>%
  wilcox_test(ratio~Treatment,
              comparisons = neg_compare,
              p.adjust.method = "fdr",
              paired = FALSE,
              alternative = "two.sided",
              detailed = TRUE) 
stat.test[5,] <-stat_neg


# count numbers
chp_n <-  count(merged_df$Treatment == "FVT_ChP")[[2,2]]
sdt_n <-  count(merged_df$Treatment == "FVT_SDT")[[2,2]]
pty_n <-  count(merged_df$Treatment == "FVT_PyT")[[2,2]]
fvt_n <-  count(merged_df$Treatment == "Unmodified_FVT")[[2,2]]
obese_n <-  count(merged_df$Treatment == "Obese_control")[[2,2]]
lean_n <-  count(merged_df$Treatment == "Lean_control")[[2,2]]
merged_df$Treatment <- factor(merged_df$Treatment,levels = groups)
p_BvsF <- ggplot(merged_df, aes(x= Treatment, y= ratio)) +
  stat_boxplot(geom ='errorbar', linetype=1, width=0.2) +
  geom_boxplot(outlier.shape = NA,alpha = 1,aes(fill=Treatment), coef=1.5, show.legend = FALSE) +
  geom_jitter(show.legend=FALSE, width=0.25, shape=21, fill="black") +
  # stat_summary(fun=mean, show.legend=FALSE, geom="crossbar", linetype=3,
  #              color="black",width=0.75, size=0.4)+
  scale_x_discrete(breaks=groups,
                   labels=c(glue("(N={chp_n})"),
                            glue("(N={sdt_n})"),
                            glue("(N={pty_n})"),
                            glue("(N={fvt_n})"),
                            glue("(N={obese_n})"),
                            glue("(N={lean_n})")
                   )
  ) +
  scale_fill_manual(values = cols) +
  labs(x="", y="Ration", title="Ration of Bacteroidetes/Firmicutes at Termination") +
  stat_pvalue_manual(stat.test,label = "p.adj.signif", tip.length = 0, size = 4,
                     y.position = c(0.6,NA,NA,NA,0.5))+ 
  theme_classic() +
  mytheme_alpha
p_BvsF

write_xlsx(stat.test, "Bacteriome/stat_result/Bacteroidets vs Firmicutes abundance.xlsx")


