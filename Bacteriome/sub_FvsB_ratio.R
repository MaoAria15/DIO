setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("Analysis0_file_loading_and_prep.R")
ps <- subset_samples(PSB, Treatment %in% c("Unmodified_FVT","Obese_control","Lean_control"))
ps@sam_data$Liver_response <- paste(ps@sam_data$Treatment,ps@sam_data$Liver_response, sep = "_")
ps <- subset_samples(ps, Sample_time_point %in% c("Termination"))
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
  select(ProjectID, Bacteroidetes, Firmicutes,Liver_response.x,Sample.x)
merged_df$ratio <- merged_df$Bacteroidetes/merged_df$Firmicutes

colnames(merged_df)[colnames(merged_df) == "Liver_response.x"] <- "Liver_response"


# Defined the respsonse group and non response group
merged_df <- subset(merged_df, Liver_response %in% c("Unmodified_FVT_Non_Response",
                                                     "Unmodified_FVT_Response",
                                                     "Obese_control_Response",
                                                     "Lean_control_Non_Response"))
merged_df$Liver_response <- factor(merged_df$Liver_response,level = c("Unmodified_FVT_Non_Response",
                                                                            "Unmodified_FVT_Response",
                                                                            "Obese_control_Response",
                                                                            "Lean_control_Non_Response"))    


# kt <- kruskal.test(data = merged_df,ratio~Treatment)


stat.test <- merged_df %>%
  wilcox_test(ratio~Liver_response,
              comparisons = liver_compare_pos,
              p.adjust.method = "fdr",
              paired = FALSE,
              alternative = "two.sided",
              detailed = TRUE) 
stat.test
stat_neg<- merged_df %>%
  wilcox_test(ratio~Liver_response,
              comparisons = liver_compare_neg,
              p.adjust.method = "fdr",
              paired = FALSE,
              alternative = "two.sided",
              detailed = TRUE) 
stat.test[3,] <-stat_neg


# count numbers
n1 <-  count(merged_df$Liver_response == "Unmodified_FVT_Non_Response")[[2,2]]
n2 <-  count(merged_df$Liver_response == "Unmodified_FVT_Response")[[2,2]]
n3 <-  count(merged_df$Liver_response == "Obese_control_Non_Response")[[2,2]]
n4 <-  count(merged_df$Liver_response == "Obese_control_Response")[[2,2]]
n5 <-  count(merged_df$Liver_response == "Lean_control_Non_Response")[[2,2]]

merged_df$Liver_response <- factor(merged_df$Liver_response,levels = liver_group)
p_BvsF <- ggplot(merged_df, aes(x= Liver_response, y= ratio)) +
  stat_boxplot(geom ='errorbar', linetype=1, width=0.2) +
  geom_boxplot(outlier.shape = NA,alpha = 1,aes(fill=Liver_response), coef=1.5, show.legend = FALSE) +
  geom_jitter(show.legend=FALSE, width=0.25, shape=21, fill="black") +
  # stat_summary(fun=mean, show.legend=FALSE, geom="crossbar", linetype=3,
  #              color="black",width=0.75, size=0.4)+
  scale_x_discrete(breaks=liver_group,
                   labels=c(glue("(N={n1})"),
                            glue("(N={n2})"),
                            # glue("(N={n3})"),
                            glue("(N={n4})"),
                            glue("(N={n5})")
                   )
  ) +
  scale_fill_manual(values = liver_cols) +
  labs(x="", y="Ration", title="Ration of Bacteroidetes/Firmicutes at Termination") +
  stat_pvalue_manual(stat.test,label = "p.adj.signif", tip.length = 0, size = 4,
                     y.position = c(NA,NA,0.8))+ 
  theme_classic() +
  mytheme_alpha
p_BvsF

write_xlsx(stat.test, "stat_result/Bacteroidets vs Firmicutes abundance.xlsx")


