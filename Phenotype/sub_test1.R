setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))
source("setup.R")

#Take the subgroup needed in the analysis : Unmodified_FVT, Obese_control, Lean_control

Sub_df <- subset(Data, !(Group %in% c("FVT_PyT", "FVT_SDT")))
#Take the subgroup needed in the analysis : Unmodified_FVT, Obese_control, Lean_control

Sub_df <- subset(Sub_df, Group %in% c("Unmodified_FVT","Obese_control","Lean_control"))
# Defined the respsonse group and non response group

Sub_df$Liver_NADFL_Response <- paste(Sub_df$Group, Sub_df$Liver_response, sep = "_" )
# Defined the respsonse group and non response group
Sub_df <- subset(Sub_df, Liver_NADFL_Response %in% c("Unmodified_FVT_Response",
                                                     "Unmodified_FVT_Non_Response",
                                                     "Obese_control_Non_Response",
                                                     "Lean_control_Response"))
Sub_df$Liver_NADFL_Response <- factor(Sub_df$Liver_NADFL_Response,level = c("Unmodified_FVT_Response",
                                                                            "Unmodified_FVT_Non_Response",
                                                                            "Obese_control_Non_Response",
                                                                            "Lean_control_Response"))          
#Count numbers

#Count numbers

count_number <- function(df)
{
  vector<-c()
  vector[1] <-  count(df$Liver_NADFL_Response == "Unmodified_FVT_Response")[[2,2]]
  vector[2] <-  count(df$Liver_NADFL_Response == "Unmodified_FVT_Non_Response")[[2,2]]
  # vector[3] <-  count(df$Liver_NADFL_Response == "Obese_control_Non_Response")[[2,2]]
  vector[3] <-  count(df$Liver_NADFL_Response == "Obese_control_Non_Response")[[2,2]]
  vector[4] <-  count(df$Liver_NADFL_Response == "Lean_control_Response")[[2,2]]
  
  return(vector)
}



#statistical function





stat_single_compare <- function(df){
  stat_bac<- df %>%
    wilcox_test(Parameter~Liver_NADFL_Response,
                comparisons = liver_compare_pos,
                p.adjust.method = "fdr",
                paired = FALSE, 
                alternative = "two.sided",
                detailed = TRUE) 
  
  stat_bac_neg<- df %>%
    wilcox_test(Parameter~Liver_NADFL_Response,
                comparisons = liver_compare_neg,
                p.adjust.method = "fdr",
                paired = FALSE,
                alternative = "two.sided",
                detailed = TRUE) 
  stat_bac_neg
  
  stat.test <- stat_bac
  stat.test[3,] <-stat_bac_neg
  
  return(stat.test)
}
#plot function
#plot function

ggplot_boxplot_abundance <- function(df)
{
  ggplot(df, aes(x= Liver_NADFL_Response, y= Parameter))+ 
    scale_fill_manual(values = liver_cols) +
    stat_boxplot(geom ='errorbar', linetype=1, width=0.2) + 
    geom_boxplot(outlier.shape = NA, aes(fill=Liver_NADFL_Response)) +
    geom_jitter(show.legend=FALSE, width=0.25, shape=21, fill="black") +
    scale_x_discrete(breaks=liver_group,
                     labels=c(glue("(N={counts[1]})"),
                              glue("(N={counts[2]})"),
                              glue("(N={counts[3]})"),
                              glue("(N={counts[4]})"))) +
    # stat_summary(fun=mean, show.legend=FALSE, geom="crossbar", linetype=3,
    #              color="black",width=0.75, linewidth=0.4)+
    theme_classic() +
    mytheme_beta
}

#eWat#####################################################################################


sheet <- list()
variable <- "ewat_mg"
df<-Sub_df[,c("Mouse_ID","Liver_NADFL_Response",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat<- stat_single_compare(df)
stat
counts <- count_number(df)
p_ewat_mg  <- ggplot_boxplot_abundance(df)+
  labs(x="", y="mg", title="Epididymal White Adipose Tissue") +
  stat_pvalue_manual(stat,label = "p.adj.signif", tip.length = 0, size =4,y.position = c(NA,NA,3000))

p_ewat_mg

sheet[["ewat_mg"]] <- stat

write_xlsx(sheet, "Phenotype/stat_result/Sub_group/ewat_mg.xlsx")



#Weight gain#####################################################################################


sheet <- list()
variable <- "Weight_gain_perc"
df<-Sub_df[,c("Mouse_ID","Liver_NADFL_Response",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat<- stat_single_compare(df)
stat
counts <- count_number(df)
p_Weight_gain_perc  <- ggplot_boxplot_abundance(df)+
  labs(x="", y="%", title="Weight Gain Percentage") +
  stat_pvalue_manual(stat,label = "p.adj.signif", tip.length = 0, size =4,y.position = c(NA,NA,350))

p_Weight_gain_perc

sheet[["Weight_gain_perc"]] <- stat

write_xlsx(sheet, "Phenotype/stat_result/Sub_group/Weight_gain_perc.xlsx")

#tAUC OGTT week 18###########################################################################

sheet <- list()
variable <- "AUC1_w18"
df<-Sub_df[,c("Mouse_ID","Liver_NADFL_Response",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat<- stat_single_compare(df)
stat
counts <- count_number(df)
p_AUC1_w18  <- ggplot_boxplot_abundance(df)+
  labs(x="", y="mg/gL", title="tAUC OGTT at week 18") +
  stat_pvalue_manual(stat,label = "p.adj.signif", tip.length = 0, size =4,y.position = c(NA,NA,2000))

p_AUC1_w18

sheet[["AUC1_w18"]] <- stat

write_xlsx(sheet, "Phenotype/stat_result/Sub_group/AUC1_w18.xlsx")

#OGTT W18########################################################################


#stat
sheet <- list()
variable <- "OGTT_t0_w18"
df<-Sub_df[,c("Mouse_ID","Liver_NADFL_Response",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat<- stat_single_compare(df)
stat
sheet[["OGTT_t0_w18"]] <- stat

variable <- "OGTT_t15_w18"
df<-Sub_df[,c("Mouse_ID","Liver_NADFL_Response",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat<- stat_single_compare(df)
stat
sheet[["OGTT_t15_w18"]] <- stat

variable <- "OGTT_t30_w18"
df<-Sub_df[,c("Mouse_ID","Liver_NADFL_Response",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat<- stat_single_compare(df)
stat
sheet[["OGTT_t30_w18"]] <- stat

variable <- "OGTT_t60_w18"
df<-Sub_df[,c("Mouse_ID","Liver_NADFL_Response",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat<- stat_single_compare(df)
stat
sheet[["OGTT_t60_w18"]] <- stat

variable <- "OGTT_t90_w18"
df<-Sub_df[,c("Mouse_ID","Liver_NADFL_Response",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat<- stat_single_compare(df)
stat
sheet[["OGTT_t90_w18"]] <- stat

variable <- "OGTT_t120_w18"
df<-Sub_df[,c("Mouse_ID","Liver_NADFL_Response",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat<- stat_single_compare(df)
stat
sheet[["OGTT_t120_w18"]] <- stat

  write_xlsx(sheet, "Phenotype/stat_result/Sub_group/OGTT_w18.xlsx")



Data_ogtt_w18 <- Sub_df[, c("Liver_NADFL_Response", "OGTT_t0_w18",  "OGTT_t15_w18", "OGTT_t30_w18",
                          "OGTT_t60_w18", "OGTT_t90_w18", "OGTT_t120_w18")]

long_data_ogtt_w18 <- pivot_longer(Data_ogtt_w18,
                                   cols = -Liver_NADFL_Response,  # Columns to keep as-is (non-pivoted)
                                   names_to = "Time",  # New column for time points
                                   values_to = "Blood_glu" ) # New column for values


long_data_ogtt_w18$Time[long_data_ogtt_w18$Time == "OGTT_t0_w18"] <-"0"
long_data_ogtt_w18$Time[long_data_ogtt_w18$Time == "OGTT_t15_w18"] <-"15"
long_data_ogtt_w18$Time[long_data_ogtt_w18$Time == "OGTT_t30_w18"] <-"30"
long_data_ogtt_w18$Time[long_data_ogtt_w18$Time == "OGTT_t60_w18"] <-"60"
long_data_ogtt_w18$Time[long_data_ogtt_w18$Time == "OGTT_t90_w18"] <-"90"
long_data_ogtt_w18$Time[long_data_ogtt_w18$Time == "OGTT_t120_w18"] <-"120"
long_data_ogtt_w18$Time <- as.numeric(long_data_ogtt_w18$Time)
# long_data_ogtt_w18$Time <- fct_relevel(long_data_ogtt_w18$Time, c("0", "15", "30",
#                                                                   "3", "OGTT_t90_w18", "OGTT_t120_w18"))
str(long_data_ogtt_w18)

p_ogtt_w18 <- ggline(long_data_ogtt_w18, x = "Time", y = "Blood_glu",
                     color = "Liver_NADFL_Response", add = ("mean_sd"),
                     numeric.x.axis = TRUE,
                     position = position_dodge(0.2),
                     size = 0.8, shape = "Liver_NADFL_Response",
                     add.params = list(size = 0.2, width = 0.5)) +
  scale_color_manual(values = liver_cols) +
  scale_y_continuous(limits = c(5, 25)) +
  scale_x_continuous(breaks = c(0, 15, 30, 60, 90, 120), labels = c("0", "15", "30", "60", "90", "120")) +
  labs(x = "Time (min)", y = "Blood glucose (mM)", title = "OGTT at Week 18") +
  theme_classic() +
  theme(legend.position = "none")+
  mytheme_with_x

p_ogtt_w18

#Liver Histo######################################################################

Data_his <- Sub_df[, c("Liver_NADFL_Response","Steatosis", "Inflammation",  "Hepatocyte_injury")]

long_data_his <- pivot_longer(Data_his,
                              cols = -Liver_NADFL_Response,  # Columns to keep as-is (non-pivoted)
                              names_to = "His",  # New column for time points
                              values_to = "Count" ) # New column for values



variable <- "Cumulative_score"
df<-Sub_df[,c("Mouse_ID","Liver_NADFL_Response",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat<- stat_single_compare(df)
stat

counts <- count_number(df)



write_xlsx(stat, "Phenotype/stat_result/Sub_group/Cumulative_score.xlsx")

p_histo <- ggplot(long_data_his, aes(fill=His, y=Count/8, x=Liver_NADFL_Response)) + 
  geom_bar(position="stack", stat="identity")+
  # scale_fill_viridis(discrete = T) +
  scale_fill_manual(values = cols_his) +
  scale_y_continuous(limits = c(0, 7)) +
  labs(x = "", y = "Cumulative score", title = "Cumulative score",fill="") +
  # stat_pvalue_manual(stat,label = "p.adj.signif", tip.length = 0, size =4,y.position = c(NA,NA,NA,NA,1900))+
  theme_classic() +
  scale_x_discrete(breaks=liver_group,
                   labels=c(
                     glue("Unmodified\nFVT\nResponse"),
                     glue("Unmodified\nFVT\nNon-response"),
                     glue("Obese\ncontrol"),
                     glue("Lean\ncontrol")
                   ))+
  mytheme_his

p_histo


ggarrange(p_ewat_mg,
          p_Weight_gain_perc,
          p_AUC1_w18,
          p_ogtt_w18,
          common.legend = TRUE,
          legend = "right")
