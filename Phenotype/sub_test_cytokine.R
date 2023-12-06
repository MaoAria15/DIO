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
    labs(x="", y="pg/ml", title=variable) +
    theme_classic() +
    mytheme_beta
}


sheet <-list()
#1. GM_CSF_p#################################################################

variable <- "GM_CSF"
df<-Sub_df[,c("Mouse_ID","Liver_NADFL_Response",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat_cytokine <- stat_single_compare(df)
stat_cytokine
counts <- count_number(df)
p_GM_CSF  <- ggplot_boxplot_abundance(df)
# stat_pvalue_manual(stat_cytokine,label = "p.adj.signif", tip.length = 0, size =4,y.position = c(5.2,4.8,4.4,4,3.6))

p_GM_CSF

sheet[["GM_CSF"]] <- stat_cytokine

#2. IFN_g#################################################################


variable <- "IFN_g"
df<-Sub_df[,c("Mouse_ID","Liver_NADFL_Response",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat_cytokine <- stat_single_compare(df)
stat_cytokine
counts <- count_number(df)
p_IFN_g  <- ggplot_boxplot_abundance(df)
# stat_pvalue_manual(stat_cytokine,label = "p.adj.signif", tip.length = 0, size =4,y.position = c(5.2,4.8,4.4,4,3.6))

p_IFN_g

sheet[["IFN_g"]] <- stat_cytokine

#3. IL_10#################################################################

variable <- "IL_10"
df<-Sub_df[,c("Mouse_ID","Liver_NADFL_Response",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat_cytokine <- stat_single_compare(df)
stat_cytokine
counts <- count_number(df)
p_IL_10  <- ggplot_boxplot_abundance(df)
# stat_pvalue_manual(stat_cytokine,label = "p.adj.signif", tip.length = 0, size =4,y.position = c(5.2,4.8,4.4,4,3.6))

p_IL_10

sheet[["IL_10"]] <- stat_cytokine
#4. IL_15#################################################################

variable <- "IL_15"
df<-Sub_df[,c("Mouse_ID","Liver_NADFL_Response",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat_cytokine <- stat_single_compare(df)
stat_cytokine
counts <- count_number(df)
p_IL_15  <- ggplot_boxplot_abundance(df)
  # stat_pvalue_manual(stat_cytokine,label = "p.adj.signif", tip.length = 0, size =4,y.position = c(42,NA,NA))

p_IL_15

sheet[["IL_15"]] <- stat_cytokine

#5. IL_17#################################################################

variable <- "IL_17"
df<-Sub_df[,c("Mouse_ID","Liver_NADFL_Response",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat_cytokine <- stat_single_compare(df)
stat_cytokine
# counts <- count_number(df)
p_IL_17  <- ggplot_boxplot_abundance(df)
# stat_pvalue_manual(stat_cytokine,label = "p.adj.signif", tip.length = 0, size =4,y.position = c(5.2,4.8,4.4,4,3.6))

p_IL_17

sheet[["IL_17"]] <- stat_cytokine

#6. IL_6#################################################################

variable <- "IL_6"
df<-Sub_df[,c("Mouse_ID","Liver_NADFL_Response",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat_cytokine <- stat_single_compare(df)
stat_cytokine
# counts <- count_number(df)
p_IL_6 <- ggplot_boxplot_abundance(df)
# stat_pvalue_manual(stat_cytokine,label = "p.adj.signif", tip.length = 0, size =4,y.position = c(5.2,4.8,4.4,4,3.6))

p_IL_6

sheet[["IL_6"]] <- stat_cytokine

#7. KC_GRO#################################################################
variable <- "KC_GRO"
df<-Sub_df[,c("Mouse_ID","Liver_NADFL_Response",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat_cytokine <- stat_single_compare(df)
stat_cytokine
counts <- count_number(df)
p_KC_GRO  <- ggplot_boxplot_abundance(df)
# stat_pvalue_manual(stat_cytokine,label = "p.adj.signif", tip.length = 0, size =4,y.position = c(5.2,4.8,4.4,4,3.6))

p_KC_GRO

sheet[["KC_GRO"]] <- stat_cytokine
#8. MIP_2 #################################################################
variable <- "MIP_2"
df<-Sub_df[,c("Mouse_ID","Liver_NADFL_Response",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat_cytokine <- stat_single_compare(df)
stat_cytokine
counts <- count_number(df)
p_MIP_2  <- ggplot_boxplot_abundance(df)+
  stat_pvalue_manual(stat_cytokine,label = "p.adj.signif", tip.length = 0, size =4,y.position = c(64,NA,60))

p_MIP_2

sheet[["MIP_2"]] <- stat_cytokine
#9. TNF_a#################################################################

variable <- "TNF_a"
df<-Sub_df[,c("Mouse_ID","Liver_NADFL_Response",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat_cytokine <- stat_single_compare(df)
stat_cytokine
counts <- count_number(df)
p_TNF_a  <- ggplot_boxplot_abundance(df)
  # stat_pvalue_manual(stat_cytokine,label = "p.adj.signif", tip.length = 0, size =4,y.position = c(20,NA,NA))

p_TNF_a

sheet[["TNF_a"]] <- stat_cytokine

#10. IL_22#################################################################

variable <- "IL_22"
df<-Sub_df[,c("Mouse_ID","Liver_NADFL_Response",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat_cytokine <- stat_single_compare(df)
stat_cytokine
# counts <- count_number(df)
p_IL_22  <- ggplot_boxplot_abundance(df)
# stat_pvalue_manual(stat_GM_CSF_p,label = "p.adj.signif", tip.length = 0, size =4,y.position = c(5.2,4.8,4.4,4,3.6))

p_IL_22

sheet[["IL_22"]] <- stat_cytokine



#merge mln facs###################################################################

p_cytokine <- ggarrange(p_GM_CSF,
                        p_IFN_g,
                        p_IL_10,
                        p_IL_15,
                        p_IL_17,
                        p_IL_6,
                        p_KC_GRO,
                        p_MIP_2,
                        p_TNF_a,
                        p_IL_22,
                        labels = c("A", "B", "C", "D", "E", "F",
                                   "G", "H","I","J"),
                        ncol = 5, nrow = 2,
                        common.legend = TRUE,
                        legend = "top",
                        font.label = list(size = 15))
p_cytokine
###################################################################################
write_xlsx(sheet, "Phenotype/stat_result/Sub_group/cytokines.xlsx")
