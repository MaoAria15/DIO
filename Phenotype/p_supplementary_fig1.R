
setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))
source("setup.R")


#Count numbers

count_number <- function(df)
{
  vector<-c()
  vector[1] <-  count(df$Group == "FVT_ChP")[[2,2]]
  vector[2] <-  count(df$Group == "FVT_SDT")[[2,2]]
  vector[3] <-  count(df$Group == "FVT_PyT")[[2,2]]
  vector[4] <-  count(df$Group == "Unmodified_FVT")[[2,2]]
  vector[5] <-  count(df$Group == "Obese_control")[[2,2]]
  vector[6] <-  count(df$Group == "Lean_control")[[2,2]]
  return(vector)
}



#statistical function


stat_single_compare <- function(df){
  stat_bac<- df %>%
    wilcox_test(Parameter~Group,
                comparisons = Obse_compare,
                p.adjust.method = "fdr",
                paired = FALSE, 
                alternative = "two.sided",
                detailed = TRUE) 
  stat_bac
  stat_bac_neg<- df %>%
    wilcox_test(Parameter~Group,
                comparisons = neg_compare,
                p.adjust.method = "fdr",
                paired = FALSE,
                alternative = "two.sided",
                detailed = TRUE) 
  stat_bac_neg
  
  stat.test <- stat_bac
  stat.test[5,] <-stat_bac_neg
  return(stat.test)
}
#plot function

ggplot_boxplot_abundance <- function(df)
{
  ggplot(df, aes(x= Group, y= Parameter))+ 
    scale_fill_manual(values = cols) +
    stat_boxplot(geom ='errorbar', linetype=1, width=0.2) + 
    geom_boxplot(outlier.shape = NA, aes(fill=Group)) +
    geom_jitter(show.legend=FALSE, width=0.25, shape=21, fill="black") +
    scale_x_discrete(breaks=groups,
                     labels=c(glue("(N={counts[1]})"),
                              glue("(N={counts[2]})"),
                              glue("(N={counts[3]})"),
                              glue("(N={counts[4]})"),
                              glue("(N={counts[5]})"),
                              glue("(N={counts[6]})"))) +
    # stat_summary(fun=mean, show.legend=FALSE, geom="crossbar", linetype=3,
    #              color="black",width=0.75, linewidth=0.4)+
    theme_classic() +
    mytheme_beta
}

ggplot_violinplot_abundance <- function(df)
{ ggviolin(df, x = "Group", y = "Parameter", add = c("mean_sd"), color ="black", shape = "dose",palette = "jco", fill = "Group")+
    scale_fill_manual(values = cols) +
    # stat_boxplot(geom ='errorbar', linetype=1, width=0.2) + 
    # geom_boxplot(outlier.shape = NA, aes(fill=Group)) +
    geom_jitter(show.legend=FALSE, width=0.25, shape=21, fill="black") +
    scale_x_discrete(breaks=groups,
                     labels=c(glue("(N={counts[1]})"),
                              glue("(N={counts[2]})"),
                              glue("(N={counts[3]})"),
                              glue("(N={counts[4]})"),
                              glue("(N={counts[5]})"),
                              glue("(N={counts[6]})"))) +
    # stat_summary(fun=mean, show.legend=FALSE, geom="crossbar", linetype=3,
    #              color="black",width=0.75, linewidth=0.4)+
    theme_classic() +
    theme(legend.position = "none")+
    mytheme_beta
}

#tAUC OGTT week 13###########################################################################
sheet <- list()
variable <- "AUC1_w13"
df<-Data[,c("Mouse_ID","Group",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat<- stat_single_compare(df)
stat
counts <- count_number(df)
p_AUC1_w13  <- ggplot_boxplot_abundance(df)+
  labs(x="", y="mg/dL", title="tAUC OGTT at week 13") +
  stat_pvalue_manual(stat,label = "p.adj.signif", tip.length = 0, size =4,y.position = c(NA,2000,NA,NA,2100))

p_AUC1_w13

sheet[["AUC1_w13"]] <- stat

write_xlsx(sheet, "Phenotype/stat_result/tAUC_W13.xlsx")

#OGTT W13########################################################################
#stat
sheet <- list()
variable <- "OGTT_t0_w13"
df<-Data[,c("Mouse_ID","Group",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat<- stat_single_compare(df)
stat
sheet[["OGTT_t0_w13"]] <- stat

variable <- "OGTT_t15_w13"
df<-Data[,c("Mouse_ID","Group",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat<- stat_single_compare(df)
stat
sheet[["OGTT_t15_w13"]] <- stat

variable <- "OGTT_t30_w13"
df<-Data[,c("Mouse_ID","Group",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat<- stat_single_compare(df)
stat
sheet[["OGTT_t30_w13"]] <- stat

variable <- "OGTT_t60_w13"
df<-Data[,c("Mouse_ID","Group",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat<- stat_single_compare(df)
stat
sheet[["OGTT_t60_w13"]] <- stat

variable <- "OGTT_t90_w13"
df<-Data[,c("Mouse_ID","Group",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat<- stat_single_compare(df)
stat
sheet[["OGTT_t90_w13"]] <- stat

variable <- "OGTT_t120_w13"
df<-Data[,c("Mouse_ID","Group",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat<- stat_single_compare(df)
stat
sheet[["OGTT_t120_w13"]] <- stat

write_xlsx(sheet, "Phenotype/stat_result/OGTT_W13.xlsx")




#plot
Data_ogtt_w13 <- Data[, c("Group", "OGTT_t0_w13",  "OGTT_t15_w13", "OGTT_t30_w13",
                          "OGTT_t60_w13", "OGTT_t90_w13", "OGTT_t120_w13")]

long_data_ogtt_w13 <- pivot_longer(Data_ogtt_w13,
                                   cols = -Group,  # Columns to keep as-is (non-pivoted)
                                   names_to = "Time",  # New column for time points
                                   values_to = "Blood_glu" ) # New column for values


long_data_ogtt_w13$Time[long_data_ogtt_w13$Time == "OGTT_t0_w13"] <-"0"
long_data_ogtt_w13$Time[long_data_ogtt_w13$Time == "OGTT_t15_w13"] <-"15"
long_data_ogtt_w13$Time[long_data_ogtt_w13$Time == "OGTT_t30_w13"] <-"30"
long_data_ogtt_w13$Time[long_data_ogtt_w13$Time == "OGTT_t60_w13"] <-"60"
long_data_ogtt_w13$Time[long_data_ogtt_w13$Time == "OGTT_t90_w13"] <-"90"
long_data_ogtt_w13$Time[long_data_ogtt_w13$Time == "OGTT_t120_w13"] <-"120"
long_data_ogtt_w13$Time <- as.numeric(long_data_ogtt_w13$Time)
# long_data_ogtt_w13$Time <- fct_relevel(long_data_ogtt_w13$Time, c("0", "15", "30",
#                                                                   "3", "OGTT_t90_w13", "OGTT_t120_w13"))
str(long_data_ogtt_w13)

p_ogtt_w13 <- ggline(long_data_ogtt_w13, x = "Time", y = "Blood_glu",
                     color = "Group", add = ("mean_sd"),
                     numeric.x.axis = TRUE,
                     position = position_dodge(0.2),
                     size = 0.8, shape = "Group",
                     add.params = list(size = 0.2, width = 0.5)) +
  # stat_line(geom ='errorbar', linetype=1, width=0.2) +
  scale_color_manual(values = cols) +
  scale_y_continuous(limits = c(5, 25)) +
  scale_x_continuous(breaks = c(0, 15, 30, 60, 90, 120), labels = c("0", "15", "30", "60", "90", "120")) +
  labs(x = "Time (min)", y = "Blood glucose (mM)", title = "OGTT at Week 13") +
  theme_classic() +
  mytheme_with_x

p_ogtt_w13






#####Sub-categories of pathological score evaluating hepatocyte injury, inflammation and steatosis###########################
sheet <- list()
# #Cumulative score
# variable <- "Cumulative_score"
# df<-Data[,c("Mouse_ID","Group",variable)]
# colnames(df)[3] <- "Parameter"
# df <- subset(df, Parameter != "")
# stat <- stat_single_compare(df)
# stat
# counts <- count_number(df)
# p_Cumulative_score  <- ggplot_violinplot_abundance(df)+
#   labs(x="", y="score", title="Cumulative score") +
# stat_pvalue_manual(stat,label = "p.adj.signif", tip.length = 0, size =4,y.position = c(NA,NA,NA,NA,8))
# 
# p_Cumulative_score
# 
# sheet[["Cumulative_score"]] <- stat


#Hepatocyte ingury
variable <- "Hepatocyte_injury"
df<-Data[,c("Mouse_ID","Group",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat <- stat_single_compare(df)
stat
counts <- count_number(df)
p_Hepatocyte_injury  <- ggplot_violinplot_abundance(df)+
  labs(x="", y="score", title="Hepatocyte injury") +
  stat_pvalue_manual(stat,label = "p.adj.signif", tip.length = 0, size =4,y.position = c(NA,NA,NA,4,5))

p_Hepatocyte_injury

sheet[["Hepatocyte_injury"]] <- stat



#Inflammation
variable <- "Inflammation"
df<-Data[,c("Mouse_ID","Group",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat <- stat_single_compare(df)
stat
counts <- count_number(df)
p_Inflammation  <- ggplot_violinplot_abundance(df)+
  labs(x="", y="score", title="Inflammation") 
  # stat_pvalue_manual(stat,label = "p.adj.signif", tip.length = 0, size =4,y.position = c(NA,NA,NA,NA,8))

p_Inflammation

sheet[["Inflammation"]] <- stat



#Steatosis
variable <- "Steatosis"
df<-Data[,c("Mouse_ID","Group",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat <- stat_single_compare(df)
stat
counts <- count_number(df)
p_Steatosis  <- ggplot_violinplot_abundance(df)+
  labs(x="", y="score", title="Steatosis") +
stat_pvalue_manual(stat,label = "p.adj.signif", tip.length = 0, size =4,y.position = c(NA,NA,NA,NA,8))

p_Steatosis

sheet[["Steatosis"]] <- stat
write_xlsx(sheet, "Phenotype/stat_result/histology.xlsx")
########### arrange p_cumulative socre#####################

p_histology<- ggarrange(         p_Hepatocyte_injury,
                                 p_Inflammation+rremove("ylab"),
                                 p_Steatosis+rremove("ylab"),
                                 ncol = 3, nrow = 1,
                                 labels = c("C","D","E"),
                                 common.legend = FALSE,
                                 font.label = list(size = 15))

p_histology

#merge ogtt13 facs

p_blood_Glucose_w13 <- ggarrange(p_AUC1_w13,p_ogtt_w13,
                        labels = c("A", "B"),
                        ncol = 2, nrow = 1,
                        widths = c(2,3),
                        common.legend = TRUE,
                        font.label = list(size = 15))
p_blood_Glucose_w13

#fig s1

p_supplementary_fig1 <- ggarrange(p_blood_Glucose_w13,p_histology,
                                  ncol = 1,
                                  nrow = 2,
                                  common.legend = TRUE,
                                  font.label = list(size=15))
p_supplementary_fig1




