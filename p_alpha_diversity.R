
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("Analysis0_file_loading_and_prep.R")


Alpha_obs <- c("Shannon")
# Alpha_obs <- c("Observe")
set.seed(111) # keep result reproductive
ps.rarefied = rarefy_even_depth(PSB, rngseed=1, sample.size=2000, replace=F)



#Normalize to mean read count
#Preparing data sheet for further differential analysis
rich <- estimate_richness(ps.rarefied)
tab <- subset(rich, select = Alpha_obs)
index <- match(rownames(rich), rownames(sample_data(ps.rarefied)))
tab$Treatment <-sample_data(ps.rarefied)$Treatment[index]
tab$Sample_time_point <- sample_data(ps.rarefied)$Sample_time_point[index]

tab$Sample <- rownames(rich)

###############################################################################
#Divide into 4 timepoints
tab

tab_arrival <- tab[tab$Sample_time_point == "Arrival",]
tab_beforefvt <- tab[tab$Sample_time_point == "Before_1st_FVT",]
tab_afterfvt <- tab[tab$Sample_time_point == "1w_after_FVT",]
tab_termination <- tab[tab$Sample_time_point == "Termination",]


# Test overall statistical difference
kt_arrival <- kruskal.test(data = tab_arrival,Shannon~Treatment)
kt_beforefvt <- kruskal.test(data = tab_beforefvt,Shannon~Treatment)
kt_afterfvt <- kruskal.test(data = tab_afterfvt,Shannon~Treatment)
kt_termination <- kruskal.test(data = tab_termination,Shannon~Treatment)

kt_chi <- c(kt_arrival$statistic[[1]],kt_beforefvt$statistic[[1]],kt_afterfvt$statistic[[1]],kt_termination$statistic[[1]])
kt_p <- c(kt_arrival$p.value[[1]],kt_beforefvt$p.value[[1]],kt_afterfvt$p.value[[1]],kt_termination$p.value[[1]])
kruskal_sum <- data.frame(
  Time_point = c("Arrival","Before_1st_FVT","1w_after_FVT","Termination"),
  chi_squared = kt_chi,
  p_value = kt_p
)

kruskal_sum

###############################################################################



#1 arrival



stat_arrival_median <- tab_arrival %>%
  group_by(Treatment) %>%
  get_summary_stats(Shannon,type = c("five_number"))
stat_arrival_median

stat_arrival_neg<- tab_arrival %>%
  wilcox_test(Shannon~Treatment,
              comparisons = neg_compare,
              p.adjust.method = "fdr",
              paired = FALSE,
              alternative = "two.sided",
              detailed = TRUE) 
stat_arrival_neg

stat_arrival <- tab_arrival %>%
  wilcox_test(Shannon~Treatment,
              comparisons = Obse_compare,
              p.adjust.method = "fdr",
              paired = FALSE,
              alternative = "two.sided",
              detailed = TRUE) 
stat_arrival


# 2 before fvt
stat_beforefvt_median <- tab_beforefvt %>%
  group_by(Treatment) %>%
  get_summary_stats(Shannon,type = "five_number")
stat_beforefvt_median

stat_beforefvt <- tab_beforefvt %>%
  wilcox_test(Shannon~Treatment,
              comparisons = Obse_compare,
              p.adjust.method = "fdr",
              paired = FALSE,
              alternative = "two.sided",
              detailed = TRUE) 
stat_beforefvt
stat_beforefvt_neg<- tab_beforefvt %>%
  wilcox_test(Shannon~Treatment,
              comparisons = neg_compare,
              p.adjust.method = "fdr",
              paired = FALSE,
              alternative = "two.sided",
              detailed = TRUE) 
stat_beforefvt_neg

# 3 afterfvt
stat_afterfvt_median <- tab_afterfvt %>%
  group_by(Treatment) %>%
  get_summary_stats(Shannon,type = "five_number")
stat_afterfvt_median

stat_afterfvt <- tab_afterfvt %>%
  wilcox_test(Shannon~Treatment,
              comparisons = Obse_compare,
              p.adjust.method = "fdr",
              paired = FALSE,
              alternative = "two.sided",
              detailed = TRUE) 
stat_afterfvt
stat_afterfvt_neg<- tab_afterfvt %>%
  wilcox_test(Shannon~Treatment,
              comparisons = neg_compare,
              p.adjust.method = "fdr",
              paired = FALSE,
              alternative = "two.sided",
              detailed = TRUE) 
stat_afterfvt_neg

# 4 termination
stat_termination_median <- tab_termination %>%
  group_by(Treatment) %>%
  get_summary_stats(Shannon,type = "five_number")
stat_termination_median

stat_termination <- tab_termination %>%
  wilcox_test(Shannon~Treatment,
              comparisons = Obse_compare,
              p.adjust.method = "fdr",
              paired = FALSE,
              alternative = "two.sided",
              detailed = TRUE) 
stat_termination
stat_termination_neg<- tab_termination %>%
  wilcox_test(Shannon~Treatment,
              comparisons = neg_compare,
              p.adjust.method = "fdr",
              paired = FALSE,
              alternative = "two.sided",
              detailed = TRUE) 
stat_termination_neg


# sheets <- list("Kruskal_test" = kruskal_sum,
# "Arrival_median" = stat_arrival_median, "Arrival" = stat_arrival,"Arrival_neg" = stat_arrival_neg,
# "Beforefvt_median" = stat_beforefvt_median, "Beforefvt" = stat_beforefvt, "Beforefvt_neg" = stat_beforefvt_neg,
# "Afterfvt_median" = stat_afterfvt_median, "Afterfvt" = stat_afterfvt,"Afterfvt_neg" = stat_afterfvt_neg,
# "Termination_median" = stat_termination_median,"Termination" = stat_termination, "Termination_neg" = stat_termination_neg
# )
# 
# write_xlsx(sheets, path = "Bacteriome/stat_result/alpha_diversity_stat.xlsx")




###############################################################################

# Get the sample size 
chp_n <-  count(tab_arrival$Treatment == "FVT_ChP")[[2,2]]
sdt_n <-  count(tab_arrival$Treatment == "FVT_SDT")[[2,2]]
pty_n <-  count(tab_arrival$Treatment == "FVT_PyT")[[2,2]]
fvt_n <-  count(tab_arrival$Treatment == "Unmodified_FVT")[[2,2]]
obese_n <-  count(tab_arrival$Treatment == "Obese_control")[[2,2]]
lean_n <-  count(tab_arrival$Treatment == "Lean_control")[[2,2]]



#no significant exist in kruskal_test

tab_arrival$Treatment <- factor(tab_arrival$Treatment,levels = groups)
p_arrival <- ggplot(tab_arrival, aes(x= Treatment, y= Shannon)) +
  stat_boxplot(geom ='errorbar', linetype=1, width=0.2) +
  geom_boxplot(outlier.shape = NA,alpha = 1,aes(fill=Treatment), coef=1.5) +
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
  labs(x="", y="Shannon diversity index", title="Arrival") +
  # stat_pvalue_manual(stat.test,label = "p.signif", tip.length = 0, size = 6)+ 
  theme_classic() +
  mytheme_alpha
p_arrival



###########################################################################

# Get the sample size 
chp_n <-  count(tab_beforefvt$Treatment == "FVT_ChP")[[2,2]]
sdt_n <-  count(tab_beforefvt$Treatment == "FVT_SDT")[[2,2]]
pty_n <-  count(tab_beforefvt$Treatment == "FVT_PyT")[[2,2]]
fvt_n <-  count(tab_beforefvt$Treatment == "Unmodified_FVT")[[2,2]]
obese_n <-  count(tab_beforefvt$Treatment == "Obese_control")[[2,2]]
lean_n <-  count(tab_beforefvt$Treatment == "Lean_control")[[2,2]]


stat.test <- stat_beforefvt
stat.test[5,] <-stat_beforefvt_neg
tab_beforefvt$Treatment <- factor(tab_beforefvt$Treatment,levels = groups)




p_beforefvt <- ggplot(tab_beforefvt, aes(x= Treatment, y= Shannon)) +
  stat_boxplot(geom ='errorbar', linetype=1, width=0.2) +
  geom_boxplot(outlier.shape = NA,alpha = 1,aes(fill=Treatment), coef=1.5) +
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
  labs(x="", y="Shannon diversity index", title="Before 1st FVT") +
  stat_pvalue_manual(stat.test,label = "p.adj.signif", tip.length = 0, size = 6,
                     y.position = c(NA,NA,NA,NA,4.7))+ 
  theme_classic() +
  mytheme_alpha
p_beforefvt
###############################################################################

# Get the sample size 
chp_n <-  count(tab_afterfvt$Treatment == "FVT_ChP")[[2,2]]
sdt_n <-  count(tab_afterfvt$Treatment == "FVT_SDT")[[2,2]]
pty_n <-  count(tab_afterfvt$Treatment == "FVT_PyT")[[2,2]]
fvt_n <-  count(tab_afterfvt$Treatment == "Unmodified_FVT")[[2,2]]
obese_n <-  count(tab_afterfvt$Treatment == "Obese_control")[[2,2]]
lean_n <-  count(tab_afterfvt$Treatment == "Lean_control")[[2,2]]

stat.test <- stat_afterfvt
stat.test[5,] <-stat_afterfvt_neg

tab_afterfvt$Treatment <- factor(tab_afterfvt$Treatment,levels = groups)




p_afterfvt <- ggplot(tab_afterfvt, aes(x= Treatment, y= Shannon)) +
  stat_boxplot(geom ='errorbar', linetype=1, width=0.2) +
  geom_boxplot(outlier.shape = NA,alpha = 1,aes(fill=Treatment), coef=1.5) +
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
  labs(x="", y="Shannon diversity index", title="1 week after 2nd FVT") +
  stat_pvalue_manual(stat.test,label = "p.adj.signif", tip.length = 0, size = 6,
                     y.position = c(NA,NA,NA,NA,4.7))+ 
  theme_classic() +
  mytheme_alpha
p_afterfvt
###############################################################################



# Get the sample size 
chp_n <-  count(tab_termination$Treatment == "FVT_ChP")[[2,2]]
sdt_n <-  count(tab_termination$Treatment == "FVT_SDT")[[2,2]]
pty_n <-  count(tab_termination$Treatment == "FVT_PyT")[[2,2]]
fvt_n <-  count(tab_termination$Treatment == "Unmodified_FVT")[[2,2]]
obese_n <-  count(tab_termination$Treatment == "Obese_control")[[2,2]]
lean_n <-  count(tab_termination$Treatment == "Lean_control")[[2,2]]

stat.test <- stat_termination
stat.test[5,] <-stat_termination_neg

tab_termination$Treatment <- factor(tab_termination$Treatment,levels = groups)




p_termination <- ggplot(tab_termination, aes(x= Treatment, y= Shannon)) +
  stat_boxplot(geom ='errorbar', linetype=1, width=0.2) +
  geom_boxplot(outlier.shape = NA,alpha = 1,aes(fill=Treatment), coef=1.5) +
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
  labs(x="", y="Shannon diversity index", title="Termination") +
  stat_pvalue_manual(stat.test,label = "p.adj.signif", tip.length = 0, size = 6,
                     y.position = c(4.9,NA,NA,NA,5.0))+ 
  theme_classic() +
  mytheme_alpha
p_termination


########### merge plot#############################

p_alpha_bac <-ggarrange(p_arrival,
                        p_beforefvt+rremove("ylab"),
                        p_afterfvt+rremove("ylab"),
                        p_termination+rremove("ylab"),
          nrow=1, ncol = 4,
          font.label = list(size = 20),
          common.legend = TRUE,
          legend = "bottom")



