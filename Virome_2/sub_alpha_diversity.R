
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("4_Analysis0_file_loading_and_prep.R")


Alpha_obs <- c("Shannon")
# Alpha_obs <- c("Observe")
set.seed(111) # keep result reproductive
ps <- subset_samples(PSV.no.realm, Treatment %in% c("Unmodified_FVT","Obese_control","Lean_control"))
ps@sam_data$Liver_response <- paste(ps@sam_data$Treatment,ps@sam_data$Liver_response, sep = "_")
ps.rarefied = rarefy_even_depth(ps, rngseed=1, sample.size=2000, replace=F)



#Normalize to mean read count
#Preparing data sheet for further differential analysis
rich <- estimate_richness(ps.rarefied)
tab <- subset(rich, select = Alpha_obs)
index <- match(rownames(rich), rownames(sample_data(ps.rarefied)))
tab$Liver_response <-sample_data(ps.rarefied)$Liver_response[index]
tab$Sample_time_point <- sample_data(ps.rarefied)$Sample_time_point[index]

tab$Sample <- rownames(rich)

###############################################################################
#Divide into 4 timepoints
tab

tab_arrival <- tab[tab$Sample_time_point == "Arrival",]

tab_afterfvt <- tab[tab$Sample_time_point == "1w_after_FVT",]
tab_termination <- tab[tab$Sample_time_point == "Termination",]


# Test overall statistical difference
kt_arrival <- kruskal.test(data = tab_arrival,Shannon~Liver_response)

kt_afterfvt <- kruskal.test(data = tab_afterfvt,Shannon~Liver_response)
kt_termination <- kruskal.test(data = tab_termination,Shannon~Liver_response)

kt_chi <- c(kt_arrival$statistic[[1]],kt_afterfvt$statistic[[1]],kt_termination$statistic[[1]])
kt_p <- c(kt_arrival$p.value[[1]],kt_afterfvt$p.value[[1]],kt_termination$p.value[[1]])
kruskal_sum <- data.frame(
  Time_point = c("Arrival","1w_after_FVT","Termination"),
  chi_squared = kt_chi,
  p_value = kt_p
)

kruskal_sum


###############################################################################



#1 arrival



stat_arrival_median <- tab_arrival %>%
  group_by(Liver_response) %>%
  get_summary_stats(Shannon,type = c("five_number"))
stat_arrival_median

stat_arrival_neg<- tab_arrival %>%
  wilcox_test(Shannon~Liver_response,
              comparisons = liver_compare_pos,
              p.adjust.method = "fdr",
              paired = FALSE,
              alternative = "two.sided",
              detailed = TRUE) 
stat_arrival_neg

stat_arrival <- tab_arrival %>%
  wilcox_test(Shannon~Liver_response,
              comparisons = liver_compare_neg,
              p.adjust.method = "fdr",
              paired = FALSE,
              alternative = "two.sided",
              detailed = TRUE) 
stat_arrival




# 3 afterfvt
stat_afterfvt_median <- tab_afterfvt %>%
  group_by(Liver_response) %>%
  get_summary_stats(Shannon,type = "five_number")
stat_afterfvt_median

stat_afterfvt <- tab_afterfvt %>%
  wilcox_test(Shannon~Liver_response,
              comparisons = liver_compare_pos,
              p.adjust.method = "fdr",
              paired = FALSE,
              alternative = "two.sided",
              detailed = TRUE) 
stat_afterfvt
stat_afterfvt_neg<- tab_afterfvt %>%
  wilcox_test(Shannon~Liver_response,
              comparisons = liver_compare_neg,
              p.adjust.method = "fdr",
              paired = FALSE,
              alternative = "two.sided",
              detailed = TRUE) 
stat_afterfvt_neg

# 4 termination
stat_termination_median <- tab_termination %>%
  group_by(Liver_response) %>%
  get_summary_stats(Shannon,type = "five_number")
stat_termination_median

stat_termination <- tab_termination %>%
  wilcox_test(Shannon~Liver_response,
              comparisons = liver_compare_pos,
              p.adjust.method = "fdr",
              paired = FALSE,
              alternative = "two.sided",
              detailed = TRUE) 
stat_termination
stat_termination_neg<- tab_termination %>%
  wilcox_test(Shannon~Liver_response,
              comparisons = liver_compare_neg,
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




###############################################################################

# Get the sample size 
n1 <-  count(tab_arrival$Liver_response == "Unmodified_FVT_Non_Response")[[2,2]]
n2 <-  count(tab_arrival$Liver_response == "Unmodified_FVT_Response")[[2,2]]
n3 <-  count(tab_arrival$Liver_response == "Obese_control_Non_Response")[[2,2]]
n4 <-  count(tab_arrival$Liver_response == "Obese_control_Response")[[2,2]]
n5 <-  count(tab_arrival$Liver_response == "Lean_control_Non_Response")[[2,2]]




#no significant exist in kruskal_test

tab_arrival$Liver_response <- factor(tab_arrival$Liver_response,levels = liver_group)
p_arrival <- ggplot(tab_arrival, aes(x= Liver_response, y= Shannon)) +
  stat_boxplot(geom ='errorbar', linetype=1, width=0.2) +
  geom_boxplot(outlier.shape = NA,alpha = 1,aes(fill=Liver_response), coef=1.5) +
  geom_jitter(show.legend=FALSE, width=0.25, shape=21, fill="black") +
  stat_summary(fun=mean, show.legend=FALSE, geom="crossbar", linetype=3,
               color="black",width=0.75, size=0.4)+
  scale_x_discrete(breaks=liver_group,
                   labels=c(glue("(N={n1})"),
                            glue("(N={n2})"),
                            glue("(N={n3})"),
                            glue("(N={n4})"),
                            glue("(N={n5})")
                   )
  ) +
  scale_fill_manual(values = liver_cols) +
  labs(x="", y="Shannon diversity index", title="Arrival") +
  # stat_pvalue_manual(stat.test,label = "p.signif", tip.length = 0, size = 6)+ 
  theme_classic() +
  mytheme_alpha
p_arrival




###############################################################################

# Get the sample size 
n1 <-  count(tab_afterfvt$Liver_response == "Unmodified_FVT_Non_Response")[[2,2]]
n2 <-  count(tab_afterfvt$Liver_response == "Unmodified_FVT_Response")[[2,2]]
n3 <-  count(tab_afterfvt$Liver_response == "Obese_control_Non_Response")[[2,2]]
n4 <-  count(tab_afterfvt$Liver_response == "Obese_control_Response")[[2,2]]
n5 <-  count(tab_arrival$Liver_response == "Lean_control_Non_Response")[[2,2]]




#no significant exist in kruskal_test

tab_afterfvt$Liver_response <- factor(tab_afterfvt$Liver_response,levels = liver_group)
p_afterfvt <- ggplot(tab_afterfvt, aes(x= Liver_response, y= Shannon)) +
  stat_boxplot(geom ='errorbar', linetype=1, width=0.2) +
  geom_boxplot(outlier.shape = NA,alpha = 1,aes(fill=Liver_response), coef=1.5) +
  geom_jitter(show.legend=FALSE, width=0.25, shape=21, fill="black") +
  stat_summary(fun=mean, show.legend=FALSE, geom="crossbar", linetype=3,
               color="black",width=0.75, size=0.4)+
  scale_x_discrete(breaks=liver_group,
                   labels=c(glue("(N={n1})"),
                            glue("(N={n2})"),
                            glue("(N={n3})"),
                            glue("(N={n4})"),
                            glue("(N={n5})")
                   )
  ) +
  scale_fill_manual(values = liver_cols) +
  labs(x="", y="Shannon diversity index", title="1 week after FVT") +
  # stat_pvalue_manual(stat.test,label = "p.signif", tip.length = 0, size = 6)+ 
  theme_classic() +
  mytheme_alpha
p_afterfvt
###############################################################################



n1 <-  count(tab_arrival$Liver_response == "Unmodified_FVT_Non_Response")[[2,2]]
n2 <-  count(tab_arrival$Liver_response == "Unmodified_FVT_Response")[[2,2]]
n3 <-  count(tab_arrival$Liver_response == "Obese_control_Non_Response")[[2,2]]
n4 <-  count(tab_arrival$Liver_response == "Obese_control_Response")[[2,2]]
n5 <-  count(tab_arrival$Liver_response == "Lean_control_Non_Response")[[2,2]]




#no significant exist in kruskal_test

tab_termination$Liver_response <- factor(tab_termination$Liver_response,levels = liver_group)
p_termination <- ggplot(tab_termination, aes(x= Liver_response, y= Shannon)) +
  stat_boxplot(geom ='errorbar', linetype=1, width=0.2) +
  geom_boxplot(outlier.shape = NA,alpha = 1,aes(fill=Liver_response), coef=1.5) +
  geom_jitter(show.legend=FALSE, width=0.25, shape=21, fill="black") +
  stat_summary(fun=mean, show.legend=FALSE, geom="crossbar", linetype=3,
               color="black",width=0.75, size=0.4)+
  scale_x_discrete(breaks=liver_group,
                   labels=c(glue("(N={n1})"),
                            glue("(N={n2})"),
                            glue("(N={n3})"),
                            glue("(N={n4})"),
                            glue("(N={n5})")
                   )
  ) +
  scale_fill_manual(values = liver_cols) +
  labs(x="", y="Shannon diversity index", title="Termination") +
  # stat_pvalue_manual(stat.test,label = "p.signif", tip.length = 0, size = 6)+ 
  theme_classic() +
  mytheme_alpha

p_termination


########### merge plot#############################

p_alpha_bac <-ggarrange(p_arrival,
                        p_afterfvt+rremove("ylab"),
                        p_termination+rremove("ylab"),
          nrow=1, ncol =3,
          font.label = list(size = 20),
          common.legend = TRUE,
          legend = "bottom")



