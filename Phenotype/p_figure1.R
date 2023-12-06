

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

#eWat#####################################################################################


sheet <- list()
variable <- "ewat_mg"
df<-Data[,c("Mouse_ID","Group",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat<- stat_single_compare(df)
stat
counts <- count_number(df)
p_ewat_mg  <- ggplot_boxplot_abundance(df)+
  labs(x="", y="mg", title="Epididymal White Adipose Tissue") +
  stat_pvalue_manual(stat,label = "p.adj.signif", tip.length = 0, size =4,y.position = c(NA,NA,NA,NA,3200))

p_ewat_mg

sheet[["ewat_mg"]] <- stat

write_xlsx(sheet, "Phenotype/stat_result/ewat_mg.xlsx")



#Weight gain#####################################################################################


sheet <- list()
variable <- "Weight_gain_perc"
df<-Data[,c("Mouse_ID","Group",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat<- stat_single_compare(df)
stat
counts <- count_number(df)
p_Weight_gain_perc  <- ggplot_boxplot_abundance(df)+
  labs(x="", y="%", title="Weight Gain Percentage") +
  stat_pvalue_manual(stat,label = "p.adj.signif", tip.length = 0, size =4,y.position = c(NA,NA,NA,NA,300))

p_Weight_gain_perc

sheet[["Weight_gain_perc"]] <- stat

write_xlsx(sheet, "Phenotype/stat_result/Weight_gain_perc.xlsx")

#tAUC OGTT week 18###########################################################################

sheet <- list()
variable <- "AUC1_w18"
df<-Data[,c("Mouse_ID","Group",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat<- stat_single_compare(df)
stat
counts <- count_number(df)
p_AUC1_w18  <- ggplot_boxplot_abundance(df)+
  labs(x="", y="mg/gL", title="tAUC OGTT at week 18") +
  stat_pvalue_manual(stat,label = "p.adj.signif", tip.length = 0, size =4,y.position = c(2000,NA,NA,NA,1900))

p_AUC1_w18

sheet[["AUC1_w18"]] <- stat

write_xlsx(sheet, "Phenotype/stat_result/AUC1_w18.xlsx")

#OGTT W18########################################################################


#stat
sheet <- list()
variable <- "OGTT_t0_w18"
df<-Data[,c("Mouse_ID","Group",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat<- stat_single_compare(df)
stat
sheet[["OGTT_t0_w18"]] <- stat

variable <- "OGTT_t15_w18"
df<-Data[,c("Mouse_ID","Group",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat<- stat_single_compare(df)
stat
sheet[["OGTT_t15_w18"]] <- stat

variable <- "OGTT_t30_w18"
df<-Data[,c("Mouse_ID","Group",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat<- stat_single_compare(df)
stat
sheet[["OGTT_t30_w18"]] <- stat

variable <- "OGTT_t60_w18"
df<-Data[,c("Mouse_ID","Group",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat<- stat_single_compare(df)
stat
sheet[["OGTT_t60_w18"]] <- stat

variable <- "OGTT_t90_w18"
df<-Data[,c("Mouse_ID","Group",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat<- stat_single_compare(df)
stat
sheet[["OGTT_t90_w18"]] <- stat

variable <- "OGTT_t120_w18"
df<-Data[,c("Mouse_ID","Group",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat<- stat_single_compare(df)
stat
sheet[["OGTT_t120_w18"]] <- stat

write_xlsx(sheet, "Phenotype/stat_result/OGTT_w18.xlsx")



Data_ogtt_w18 <- Data[, c("Group", "OGTT_t0_w18",  "OGTT_t15_w18", "OGTT_t30_w18",
                                    "OGTT_t60_w18", "OGTT_t90_w18", "OGTT_t120_w18")]

long_data_ogtt_w18 <- pivot_longer(Data_ogtt_w18,
                                   cols = -Group,  # Columns to keep as-is (non-pivoted)
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
                      color = "Group", add = ("mean_sd"),
                      numeric.x.axis = TRUE,
                      position = position_dodge(0.2),
                      size = 0.8, shape = "Group",
                      add.params = list(size = 0.2, width = 0.5)) +
  scale_color_manual(values = cols) +
  scale_y_continuous(limits = c(5, 25)) +
  scale_x_continuous(breaks = c(0, 15, 30, 60, 90, 120), labels = c("0", "15", "30", "60", "90", "120")) +
  labs(x = "Time (min)", y = "Blood glucose (mM)", title = "OGTT at Week 18") +
  theme_classic() +
  mytheme_with_x

p_ogtt_w18

#Liver Histo######################################################################

Data_his <- Data[, c("Group","Steatosis", "Inflammation",  "Hepatocyte_injury")]

long_data_his <- pivot_longer(Data_his,
                                   cols = -Group,  # Columns to keep as-is (non-pivoted)
                                   names_to = "His",  # New column for time points
                                   values_to = "Count" ) # New column for values



variable <- "Cumulative_score"
df<-Data[,c("Mouse_ID","Group",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat<- stat_single_compare(df)
stat
counts <- count_number(df)

sheet[["Cumulative_score"]] <- stat

write_xlsx(sheet, "Phenotype/stat_result/Cumulative_score.xlsx")


p_histo <- ggplot(long_data_his, aes(fill=His, y=Count/8, x=Group)) + 
  geom_bar(position="stack", stat="identity")+
  # scale_fill_viridis(discrete = T) +
  scale_fill_manual(values = cols_his) +
  scale_y_continuous(limits = c(0, 7)) +
  labs(x = "", y = "Cumulative score", title = "Cumulative score",fill="") +
  # stat_pvalue_manual(stat,label = "p.adj.signif", tip.length = 0, size =4,y.position = c(NA,NA,NA,NA,1900))+
  theme_classic() +
  scale_x_discrete(breaks=groups,
                   labels=c(
                     glue("FVT\nChP"),
                     glue("FVT\nSDT"),
                     glue("FVT\nPyT"),
                     glue("Unmodified\nFVT"),
                     glue("Obese\ncontrol"),
                     glue("Lean\ncontrol")
                   ))+
  mytheme_his

p_histo








#merge fig1


p_fig1<- ggarrange(ggarrange(p_Weight_gain_perc,
                   p_ewat_mg,
                   p_ogtt_w18,
                   p_AUC1_w18,
                   labels = c("A", "B","C","D"),
                   ncol = 2, nrow = 2,
                   common.legend = TRUE,
                   font.label = list(size = 15)),
                 ggarrange( p_histo,p_histo,
                  labels= c("E","F"),
                  ncol = 2,nrow=1,
                  common.legend=TRUE,
                  font.label = list(size=15)),
                 ncol=1,nrow=2,
                 common.legend=FALSE,
                 font.label = list(size=15),
                 heights = c(4,2),
                 widths = c(4,4))
p_fig1

