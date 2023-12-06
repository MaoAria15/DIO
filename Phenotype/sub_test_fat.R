setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))
source("setup.R")

# Change scale %
fat_columns <- grep("^fat_", names(Data), value = TRUE) 
Data[, fat_columns] <- Data[, fat_columns] * 100

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
    # stat_summary(fun=mean, show.legend=FALSE, geom="crossbar", linetype=3,
    #              color="black",width=0.75, linewidth=0.4)+
    theme_classic() +
    mytheme_beta
}
sheet <- list()
#1. p_fat_Dendritic_cells#################################################################

variable <- "fat_Dendritic_cells"
df<-Sub_df[,c("Mouse_ID","Liver_NADFL_Response",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat_mln <- stat_single_compare(df)
stat_mln

counts <- count_number(df)

p_fat_Dendritic_cells  <- ggplot_boxplot_abundance(df)+
  labs(x="", y="%CD11c+ of CD45+", title="Dendritic Cells") +
stat_pvalue_manual(stat_mln,label = "p.adj.signif", tip.length = 0, size =4,y.position = c(3.5,NA,3.2))

p_fat_Dendritic_cells

sheet[["fat_Dendritic_cells"]] <- stat_mln



#2. p_fat_Macrophages#################################################################
variable <- "fat_Macrophages"
df<-Sub_df[,c("Mouse_ID","Liver_NADFL_Response",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat_mln <- stat_single_compare(df)
stat_mln

counts <- count_number(df)
p_fat_Macrophages  <- ggplot_boxplot_abundance(df)+
  labs(x="", y="%F4.80+ of CD45+", title="Macrophages")  
# stat_pvalue_manual(stat_mln,label = "p.adj.signif", tip.length = 0, size =4,y.position = c(5.2,4.8,4.4,4,3.6))

p_fat_Macrophages

sheet[["p_fat_Macrophages"]] <- stat_mln
#3. fat_M1_macrophages#################################################################

variable <- "fat_M1_macrophages"
df<-Sub_df[,c("Mouse_ID","Liver_NADFL_Response",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat_mln <- stat_single_compare(df)
stat_mln

counts <- count_number(df)
p_fat_M1_macrophages  <- ggplot_boxplot_abundance(df)+
  labs(x="", y="%CD11c+ of F4.80+", title="M1 Macrophages") +
  stat_pvalue_manual(stat_mln,label = "p.adj.signif", tip.length = 0, size =4,y.position = c(8,NA,7.8))

p_fat_M1_macrophages

sheet[["fat_M1_macrophages"]] <- stat_mln

#4. fat_B_cells#################################################################

variable <- "fat_B_cells"
df<-Sub_df[,c("Mouse_ID","Liver_NADFL_Response",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat_mln <- stat_single_compare(df)
stat_mln

counts <- count_number(df)
p_fat_B_cells  <- ggplot_boxplot_abundance(df)+
  labs(x="", y="%CD19c+ of CD45+", title="B Cells") 
# stat_pvalue_manual(stat_mln,label = "p.adj.signif", tip.length = 0, size =4,y.position = c(5.2,4.8,4.4,4,3.6))

p_fat_B_cells

sheet[["fat_B_cells"]] <- stat_mln

#5. fat_Th_cells#################################################################

variable <- "fat_Th_cells"
df<-Sub_df[,c("Mouse_ID","Liver_NADFL_Response",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat_mln <- stat_single_compare(df)
stat_mln

counts <- count_number(df)
p_fat_Th_cells  <- ggplot_boxplot_abundance(df)+
  labs(x="", y="%CD4+ of TCRab+", title="T Helper Cells") +
  stat_pvalue_manual(stat_mln,label = "p.adj.signif", tip.length = 0, size =4,y.position = c(12,13,NA))
p_fat_Th_cells

sheet[["fat_Th_cells"]] <- stat_mln

#6. fat_Cytotoxic_T_cells#################################################################

variable <- "fat_Cytotoxic_T_cells"
df<-Sub_df[,c("Mouse_ID","Liver_NADFL_Response",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat_mln <- stat_single_compare(df)
stat_mln

counts <- count_number(df)
p_fat_Cytotoxic_T_cells  <- ggplot_boxplot_abundance(df)+
  labs(x="", y="%CD8a+ of TCRab+", title="Cytotoxic T cells") 
  # stat_pvalue_manual(stat_mln,label = "p.adj.signif", tip.length = 0, size =4,y.position = c(NA,NA,30))

p_fat_Cytotoxic_T_cells

sheet[["fat_Cytotoxic_T_cells"]] <- stat_mln

#7. fat_Central_memory_T_cells_CD8#################################################################

variable <- "fat_Central_memory_T_cells_CD8"
df<-Sub_df[,c("Mouse_ID","Liver_NADFL_Response",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat_mln <- stat_single_compare(df)
stat_mln

counts <- count_number(df)
p_fat_Central_memory_T_cells_CD8  <- ggplot_boxplot_abundance(df)+
  labs(x="", y="%CD44+CD62l+ of CD8+", title="Central memory CD8+ T cells") +
  stat_pvalue_manual(stat_mln,label = "p.adj.signif", tip.length = 0, size =4,y.position = c(50,NA,48))

p_fat_Central_memory_T_cells_CD8

sheet[["fat_Central_memory_T_cells_CD8"]] <- stat_mln

#8. fat_Effector_memory_T_cells_CD8 #################################################################

variable <- "fat_Effector_memory_T_cells_CD8"
df<-Sub_df[,c("Mouse_ID","Liver_NADFL_Response",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat_mln <- stat_single_compare(df)
stat_mln
counts <- count_number(df)
p_fat_Effector_memory_T_cells_CD8 <- ggplot_boxplot_abundance(df)+
  labs(x="", y="%CD44+CD62l- of CD8+", title="Effector memory CD8+ T cells")  
# stat_pvalue_manual(stat_mln,label = "p.adj.signif", tip.length = 0, size =4,y.position = c(5.2,4.8,4.4,4,3.6))

p_fat_Effector_memory_T_cells_CD8

sheet[["fat_Effector_memory_T_cells_CD8"]] <- stat_mln

#9. fat_Naive_T_cells_CD8#################################################################

variable <- "fat_Naive_T_cells_CD8"
df<-Sub_df[,c("Mouse_ID","Liver_NADFL_Response",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat_mln <- stat_single_compare(df)
stat_mln
counts <- count_number(df)
p_fat_Naive_T_cells_CD8  <- ggplot_boxplot_abundance(df)+
  labs(x="", y="%CD44-CD62l+ of CD8+", title="Naïve CD8+ T cells") +
stat_pvalue_manual(stat_mln,label = "p.adj.signif", tip.length = 0, size =4,y.position = c(NA,NA,30))

p_fat_Naive_T_cells_CD8

sheet[["fat_Naive_T_cells_CD8"]] <- stat_mln

#10. fat_Central_memory_T_cells_CD4################################################################

variable <- "fat_Central_memory_T_cells_CD4"
df<-Sub_df[,c("Mouse_ID","Liver_NADFL_Response",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat_mln <- stat_single_compare(df)
stat_mln
counts <- count_number(df)
p_fat_Central_memory_T_cells_CD4  <- ggplot_boxplot_abundance(df)+
  labs(x="", y="%CD44+CD62l+ of CD4+", title="Central memory CD4+ T cells") 
# stat_pvalue_manual(stat_mln,label = "p.adj", tip.length = 0, size =2,y.position = c(NA,77,NA,NA,NA))

p_fat_Central_memory_T_cells_CD4

sheet[["fat_Central_memory_T_cells_CD4"]] <- stat_mln

#11. fat_Effector_memory_T_cells_CD4#################################################################

variable <- "fat_Effector_memory_T_cells_CD4"
df<-Sub_df[,c("Mouse_ID","Liver_NADFL_Response",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat_mln <- stat_single_compare(df)
stat_mln
counts <- count_number(df)
p_fat_Effector_memory_T_cells_CD4  <- ggplot_boxplot_abundance(df)+
  labs(x="", y="%CD44+CD62l- of CD4+", title="Effector memory CD4+ T cells") 
# stat_pvalue_manual(stat_mln,label = "p.adj.signif", tip.length = 0, size =4,y.position = c(5.2,4.8,4.4,4,3.6))

p_fat_Effector_memory_T_cells_CD4

sheet[["fat_Effector_memory_T_cells_CD4"]] <- stat_mln

#12. fat_Naive_T_cells_CD4#################################################################


variable <- "fat_Naive_T_cells_CD4"
df<-Sub_df[,c("Mouse_ID","Liver_NADFL_Response",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat_mln <- stat_single_compare(df)
stat_mln
counts <- count_number(df)
p_fat_Naive_T_cells_CD4  <- ggplot_boxplot_abundance(df)+
  labs(x="", y="%CD44-CD62l+ of CD4+", title="Naïve CD4+ T cells")  
# stat_pvalue_manual(stat_mln,label = "p.adj", tip.length = 0, size =2,y.position = c(NA,20,NA,NA,NA))

p_fat_Naive_T_cells_CD4

sheet[["fat_Naive_T_cells_CD4"]] <- stat_mln


#merge mln facs

p_fat_facs <- ggarrange(p_fat_Dendritic_cells,
                        p_fat_Macrophages,
                        p_fat_M1_macrophages,
                        p_fat_B_cells,
                        p_fat_Th_cells,
                        p_fat_Cytotoxic_T_cells,
                        p_fat_Central_memory_T_cells_CD8,
                        p_fat_Effector_memory_T_cells_CD8,
                        p_fat_Naive_T_cells_CD8,
                        p_fat_Central_memory_T_cells_CD4,
                        p_fat_Effector_memory_T_cells_CD4,
                        p_fat_Naive_T_cells_CD4,
                        labels = c("A", "B", "C", "D", "E", "F",
                                   "G", "H","I","J","K","L"),
                        ncol = 3, nrow = 4,
                        common.legend = TRUE,
                        legend = "top",
                        font.label = list(size = 15))
p_fat_facs


###################################################################################
write_xlsx(sheet, "Phenotype/stat_result/Sub_group/FACs_fat.xlsx")
