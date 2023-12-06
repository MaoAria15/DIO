setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))
source("setup.R")

# Change scale %
mln_columns <- grep("^mln_", names(Data), value = TRUE) 
Data[, mln_columns] <- Data[, mln_columns] * 100


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

sheet <- list()

#1. p_mln_Dendritic_cells#################################################################

variable <- "mln_Dendritic_cells"
df<-Data[,c("Mouse_ID","Group",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat_mln <- stat_single_compare(df)
stat_mln

counts <- count_number(df)

p_mln_Dendritic_cells  <- ggplot_boxplot_abundance(df)+
  labs(x="", y="%CD11c+ of CD45+", title="Dendritic Cells") 
# stat_pvalue_manual(stat_mln,label = "p.adj.signif", tip.length = 0, size =4,y.position = c(5.2,4.8,4.4,4,3.6))

p_mln_Dendritic_cells

sheet[["mln_Dendritic_cells"]] <- stat_mln



#2. p_mln_Macrophages#################################################################
variable <- "mln_Macrophages"
df<-Data[,c("Mouse_ID","Group",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat_mln <- stat_single_compare(df)
stat_mln

counts <- count_number(df)
p_mln_Macrophages  <- ggplot_boxplot_abundance(df)+
  labs(x="", y="%F4.80+ of CD45+", title="Macrophages")  
# stat_pvalue_manual(stat_mln,label = "p.adj.signif", tip.length = 0, size =4,y.position = c(5.2,4.8,4.4,4,3.6))

p_mln_Macrophages

sheet[["p_mln_Macrophages"]] <- stat_mln
#3. mln_M1_macrophages#################################################################

variable <- "mln_M1_macrophages"
df<-Data[,c("Mouse_ID","Group",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat_mln <- stat_single_compare(df)
stat_mln

counts <- count_number(df)
p_mln_M1_macrophages  <- ggplot_boxplot_abundance(df)+
  labs(x="", y="%CD11c+ of F4.80+", title="M1 Macrophages") +
stat_pvalue_manual(stat_mln,label = "p.adj.signif", tip.length = 0, size =4,y.position = c(NA,NA,NA,NA,30))

p_mln_M1_macrophages

sheet[["mln_M1_macrophages"]] <- stat_mln

#4. mln_B_cells#################################################################

variable <- "mln_B_cells"
df<-Data[,c("Mouse_ID","Group",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat_mln <- stat_single_compare(df)
stat_mln

counts <- count_number(df)
p_mln_B_cells  <- ggplot_boxplot_abundance(df)+
  labs(x="", y="%CD19c+ of CD45+", title="B Cells") 
# stat_pvalue_manual(stat_mln,label = "p.adj.signif", tip.length = 0, size =4,y.position = c(5.2,4.8,4.4,4,3.6))

p_mln_B_cells

sheet[["mln_B_cells"]] <- stat_mln

#5. mln_Th_cells#################################################################

variable <- "mln_Th_cells"
df<-Data[,c("Mouse_ID","Group",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat_mln <- stat_single_compare(df)
stat_mln

counts <- count_number(df)
p_mln_Th_cells  <- ggplot_boxplot_abundance(df)+
  labs(x="", y="%CD4+ of TCRab+", title="T Helper Cells") +
  stat_pvalue_manual(stat_mln,label = "p.adj.signif", tip.length = 0, size =4,y.position = c(NA,NA,NA,NA,12))
p_mln_Th_cells

sheet[["mln_Th_cells"]] <- stat_mln

#6. mln_Cytotoxic_T_cells#################################################################

variable <- "mln_Cytotoxic_T_cells"
df<-Data[,c("Mouse_ID","Group",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat_mln <- stat_single_compare(df)
stat_mln

counts <- count_number(df)
p_mln_Cytotoxic_T_cells  <- ggplot_boxplot_abundance(df)+
  labs(x="", y="%CD8a+ of TCRab+", title="Cytotoxic T cells") +
stat_pvalue_manual(stat_mln,label = "p.adj.signif", tip.length = 0, size =4,y.position = c(NA,NA,NA,NA,30))

p_mln_Cytotoxic_T_cells

sheet[["mln_Cytotoxic_T_cells"]] <- stat_mln

#7. mln_Central_memory_T_cells_CD8#################################################################

variable <- "mln_Central_memory_T_cells_CD8"
df<-Data[,c("Mouse_ID","Group",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat_mln <- stat_single_compare(df)
stat_mln

counts <- count_number(df)
p_mln_Central_memory_T_cells_CD8  <- ggplot_boxplot_abundance(df)+
  labs(x="", y="%CD44+CD62l+ of CD8+", title="Central memory CD8+ T cells") +
stat_pvalue_manual(stat_mln,label = "p.adj", tip.length = 0, size =2,y.position = c(NA,65,NA,NA,NA))

p_mln_Central_memory_T_cells_CD8

sheet[["mln_Central_memory_T_cells_CD8"]] <- stat_mln

#8. mln_Effector_memory_T_cells_CD8 #################################################################

variable <- "mln_Effector_memory_T_cells_CD8"
df<-Data[,c("Mouse_ID","Group",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat_mln <- stat_single_compare(df)
stat_mln
counts <- count_number(df)
p_mln_Effector_memory_T_cells_CD8 <- ggplot_boxplot_abundance(df)+
  labs(x="", y="%CD44+CD62l- of CD8+", title="Effector memory CD8+ T cells")  
# stat_pvalue_manual(stat_mln,label = "p.adj.signif", tip.length = 0, size =4,y.position = c(5.2,4.8,4.4,4,3.6))

p_mln_Effector_memory_T_cells_CD8

sheet[["mln_Effector_memory_T_cells_CD8"]] <- stat_mln

#9. mln_Naive_T_cells_CD8#################################################################

variable <- "mln_Naive_T_cells_CD8"
df<-Data[,c("Mouse_ID","Group",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat_mln <- stat_single_compare(df)
stat_mln
counts <- count_number(df)
p_mln_Naive_T_cells_CD8  <- ggplot_boxplot_abundance(df)+
  labs(x="", y="%CD44-CD62l+ of CD8+", title="Naïve CD8+ T cells") 
# stat_pvalue_manual(stat_mln,label = "p.adj.signif", tip.length = 0, size =4,y.position = c(5.2,4.8,4.4,4,3.6))

p_mln_Naive_T_cells_CD8

sheet[["mln_Naive_T_cells_CD8"]] <- stat_mln

#10. mln_Central_memory_T_cells_CD4################################################################

variable <- "mln_Central_memory_T_cells_CD4"
df<-Data[,c("Mouse_ID","Group",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat_mln <- stat_single_compare(df)
stat_mln
counts <- count_number(df)
p_mln_Central_memory_T_cells_CD4  <- ggplot_boxplot_abundance(df)+
  labs(x="", y="%CD44+CD62l+ of CD4+", title="Central memory CD4+ T cells") +
stat_pvalue_manual(stat_mln,label = "p.adj", tip.length = 0, size =2,y.position = c(NA,77,NA,NA,NA))

p_mln_Central_memory_T_cells_CD4

sheet[["mln_Central_memory_T_cells_CD4"]] <- stat_mln

#11. mln_Effector_memory_T_cells_CD4#################################################################

variable <- "mln_Effector_memory_T_cells_CD4"
df<-Data[,c("Mouse_ID","Group",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat_mln <- stat_single_compare(df)
stat_mln
counts <- count_number(df)
p_mln_Effector_memory_T_cells_CD4  <- ggplot_boxplot_abundance(df)+
  labs(x="", y="%CD44+CD62l- of CD4+", title="Effector memory CD4+ T cells") 
# stat_pvalue_manual(stat_mln,label = "p.adj.signif", tip.length = 0, size =4,y.position = c(5.2,4.8,4.4,4,3.6))

p_mln_Effector_memory_T_cells_CD4

sheet[["mln_Effector_memory_T_cells_CD4"]] <- stat_mln

#12. mln_Naive_T_cells_CD4#################################################################


variable <- "mln_Naive_T_cells_CD4"
df<-Data[,c("Mouse_ID","Group",variable)]
colnames(df)[3] <- "Parameter"
df <- subset(df, Parameter != "")
stat_mln <- stat_single_compare(df)
stat_mln
counts <- count_number(df)
p_mln_Naive_T_cells_CD4  <- ggplot_boxplot_abundance(df)+
  labs(x="", y="%CD44-CD62l+ of CD4+", title="Naïve CD4+ T cells")  +
stat_pvalue_manual(stat_mln,label = "p.adj", tip.length = 0, size =2,y.position = c(NA,20,NA,NA,NA))

p_mln_Naive_T_cells_CD4

sheet[["mln_Naive_T_cells_CD4"]] <- stat_mln


#merge mln facs

p_mln_facs <- ggarrange(p_mln_Dendritic_cells,
                        p_mln_Macrophages,
                        p_mln_M1_macrophages,
                        p_mln_B_cells,
                        p_mln_Th_cells,
                        p_mln_Cytotoxic_T_cells,
                        p_mln_Central_memory_T_cells_CD8,
                        p_mln_Effector_memory_T_cells_CD8,
                        p_mln_Naive_T_cells_CD8,
                        p_mln_Central_memory_T_cells_CD4,
                        p_mln_Effector_memory_T_cells_CD4,
                        p_mln_Naive_T_cells_CD4,
                        labels = c("A", "B", "C", "D", "E", "F",
                                   "G", "H","I","J","K","L"),
                        ncol = 3, nrow = 4,
                        common.legend = TRUE,
                        font.label = list(size = 15))
p_mln_facs


###################################################################################
write_xlsx(sheet, "Phenotype/stat_result/FACs_mln.xlsx")
