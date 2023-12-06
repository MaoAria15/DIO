setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("Analysis0_file_loading_and_prep.R")

ps <- subset_samples(PSB, Sample_time_point %in% c("Termination"))
ps.rel <- transform_sample_counts(ps, function(x) x / sum(x))
ps.rel@sam_data$Liver_response <- paste(ps.rel@sam_data$Treatment,ps.rel@sam_data$Liver_response, sep = "_")
ps.rel<- subset_samples(ps.rel, Liver_response %in% c("Unmodified_FVT_Response",
                                                      "Unmodified_FVT_Non_Response",
                                                      "Obese_control_Non_Response",
                                                      "Lean_control_Response"))
ps.rel@sam_data$Liver_response <- factor(ps.rel@sam_data$Liver_response,level = c("Unmodified_FVT_Response",
                                                                                  "Unmodified_FVT_Non_Response",
                                                                                  "Obese_control_Non_Response",
                                                                                  "Lean_control_Response"))    



#####setting####################################################################
df <- data.frame(sample_data(ps.rel))

# count numbers
n1 <-  count(df$Liver_response == "Unmodified_FVT_Response")[[2,2]]
n2 <-  count(df$Liver_response == "Unmodified_FVT_Non_Response")[[2,2]]
n3 <-  count(df$Liver_response == "Obese_control_Non_Response")[[2,2]]
n4 <-  count(df$Liver_response == "Lean_control_Response")[[2,2]]


counts<-c(glue("(N={n1})"),
                glue("(N={n2})"),
                glue("(N={n3})"),
                glue("(N={n4})"))
#statistical function
stat_single_compare <- function(df){
  stat_bac<- df %>%
    wilcox_test(Abundance~Liver_response,
                comparisons = liver_compare_pos,
                p.adjust.method = "fdr",
                paired = FALSE, 
                alternative = "two.sided",
                detailed = TRUE) 
  stat_bac
  stat_bac_neg<- df %>%
    wilcox_test(Abundance~Liver_response,
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
ggplot_boxplot_abundance <- function(df)
{
  ggplot(df_bac,aes(x=Liver_response,y=Abundance)) +
    stat_boxplot(geom ='errorbar', linetype=1, width=0.2) +
    geom_boxplot(outlier.shape = NA,alpha = 1,aes(fill=Liver_response), coef=1.5) +
    geom_jitter(show.legend=FALSE, width=0.25, shape=21, fill="black") +
    scale_x_discrete(breaks=liver_group,
                     labels=counts) +
    scale_fill_manual(values = liver_cols) +
    ylab("Mean Relative Abundance (%)")+
    theme_classic() +
    mytheme_alpha
}

sheet <- list()

#Allobaculum##########################################################################
ps0 <- tax_glom(subset_taxa(ps.rel, Genus=="Allobaculum"), "Genus", NArm = FALSE)
df <- psmelt(ps0)
df[df==""] <- NA
#Select last non-empty taxonomic rank
df$tax <- apply(df, 1, function(x) tail(stats::na.omit(x), 1))
#Add "unknown" to taxonomy columns where genus is not classified 
df$tax[is.na(df$Genus)] <- paste0(df$tax[is.na(df$Genus)],"unknown genus")

df_bac<-df
df_bac$Liver_response <- factor(df_bac$Liver_response, levels = liver_group)
df_bac$Abundance <- df_bac$Abundance*100

#stat
stat_bac<- stat_single_compare(df_bac)
df_bac$Liver_response <- factor(df_bac$Liver_response, levels = liver_group)



p_Allobaculum <- ggplot_boxplot_abundance(df_bac)+
  stat_pvalue_manual(stat_bac,label = "p.adj.signif", tip.length = 0, size = 4,
                     y.position = c(15,17,19))+
  ggtitle("Allobaculum") 
  
p_Allobaculum 
  
sheet <- list("Allobaculum" = stat_bac)

#Bacteroides##########################################################################
ps0 <- tax_glom(subset_taxa(ps.rel, Genus=="Bacteroides"), "Genus", NArm = FALSE)
df <- psmelt(ps0)
df[df==""] <- NA
#Select last non-empty taxonomic rank
df$tax <- apply(df, 1, function(x) tail(stats::na.omit(x), 1))
#Add "unknown" to taxonomy columns where genus is not classified 
df$tax[is.na(df$Genus)] <- paste0(df$tax[is.na(df$Genus)],"unknown genus")

df_bac<-df
df_bac$Liver_response <- factor(df_bac$Liver_response, levels = liver_group)
df_bac$Abundance <- df_bac$Abundance*100

#stat
stat_bac<- stat_single_compare(df_bac)
df_bac$Liver_response <- factor(df_bac$Liver_response, levels = liver_group)


p_Bacteroides <- ggplot_boxplot_abundance(df_bac)+
  # stat_pvalue_manual(stat_bac,label = "p.adj", tip.length = 0, size = 3,
  #                    y.position = c(15,17,NA))+
  ggtitle("Bacteroides") 

p_Bacteroides

sheet[["Bacteroides"]] <- stat_bac

#Bifidobacterium##########################################################################
ps0 <- tax_glom(subset_taxa(ps.rel, Genus=="Bifidobacterium"), "Genus", NArm = FALSE)
df <- psmelt(ps0)
df[df==""] <- NA
#Select last non-empty taxonomic rank
df$tax <- apply(df, 1, function(x) tail(stats::na.omit(x), 1))
#Add "unknown" to taxonomy columns where genus is not classified 
df$tax[is.na(df$Genus)] <- paste0(df$tax[is.na(df$Genus)],"unknown genus")

df_bac<-df
df_bac$Liver_response <- factor(df_bac$Liver_response, levels = liver_group)
df_bac$Abundance <- df_bac$Abundance*100

#stat
  stat_bac<- df_bac %>%
    wilcox_test(Abundance~Liver_response,
                comparisons = liver_compare_neg,
                p.adjust.method = "fdr",
                paired = FALSE,
                alternative = "two.sided",
                detailed = TRUE) 
# stat_bac<- stat_single_compare(df_bac)
df_bac$Liver_response <- factor(df_bac$Liver_response, levels = liver_group)



p_Bifidobacterium <- ggplot_boxplot_abundance(df_bac)+
  stat_pvalue_manual(stat_bac,label = "p.adj.signif", tip.length = 0, size = 4,
                     y.position = c(2.8))+
  ggtitle("Bifidobacterium") 

p_Bifidobacterium

sheet[["Bifidobacterium"]] <- stat_bac

#Clostridium##########################################################################
ps0 <- tax_glom(subset_taxa(ps.rel, Genus=="Clostridium"), "Genus", NArm = FALSE)
df <- psmelt(ps0)
df[df==""] <- NA
#Select last non-empty taxonomic rank
df$tax <- apply(df, 1, function(x) tail(stats::na.omit(x), 1))
#Add "unknown" to taxonomy columns where genus is not classified 
df$tax[is.na(df$Genus)] <- paste0(df$tax[is.na(df$Genus)],"unknown genus")
df_bac<-df
df_bac$Liver_response <- factor(df_bac$Liver_response, levels = liver_group)
df_bac$Abundance <- df_bac$Abundance*100

#stat
stat_bac<- stat_single_compare(df_bac)
df_bac$Liver_response <- factor(df_bac$Liver_response, levels = liver_group)


p_Clostridium <- ggplot_boxplot_abundance(df_bac)+
  # stat_pvalue_manual(stat_bac,label = "p.adj.signif", tip.length = 0, size = 6,
  #                    y.position = c(NA,NA,15,17,19))+
  ggtitle("Clostridium") 

p_Clostridium 

sheet[["Clostridium"]] <- stat_bac


#Dehalobacterium##########################################################################
ps0 <- tax_glom(subset_taxa(ps.rel, Genus=="Dehalobacterium"), "Genus", NArm = FALSE)
df <- psmelt(ps0)
df[df==""] <- NA
#Select last non-empty taxonomic rank
df$tax <- apply(df, 1, function(x) tail(stats::na.omit(x), 1))
#Add "unknown" to taxonomy columns where genus is not classified 
df$tax[is.na(df$Genus)] <- paste0(df$tax[is.na(df$Genus)],"unknown genus")
df_bac<-df
df_bac$Liver_response <- factor(df_bac$Liver_response, levels = liver_group)
df_bac$Abundance <- df_bac$Abundance*100

#stat
stat_bac<- stat_single_compare(df_bac)
df_bac$Liver_response <- factor(df_bac$Liver_response, levels = liver_group)


p_Dehalobacterium <- ggplot_boxplot_abundance(df_bac)+
  # stat_pvalue_manual(stat_bac,label = "p.adj.signif", tip.length = 0, size = 6,
  #                    y.position = c(NA,NA,15,17,19))+
  ggtitle("Dehalobacterium")

p_Dehalobacterium 

sheet[["Dehalobacterium"]] <- stat_bac


#Lactobacillus##########################################################################
ps0 <- tax_glom(subset_taxa(ps.rel, Genus=="Lactobacillus"), "Genus", NArm = FALSE)
df <- psmelt(ps0)
df[df==""] <- NA
#Select last non-empty taxonomic rank
df$tax <- apply(df, 1, function(x) tail(stats::na.omit(x), 1))
#Add "unknown" to taxonomy columns where genus is not classified 
df$tax[is.na(df$Genus)] <- paste0(df$tax[is.na(df$Genus)],"unknown genus")
df_bac<-df
df_bac$Liver_response <- factor(df_bac$Liver_response, levels = liver_group)
df_bac$Abundance <- df_bac$Abundance*100

#stat
stat_bac<- stat_single_compare(df_bac)
df_bac$Liver_response <- factor(df_bac$Liver_response, levels = liver_group)


p_Lactobacillus <- ggplot_boxplot_abundance(df_bac)+
  stat_pvalue_manual(stat_bac,label = "p.adj.signif", tip.length = 0, size = 3,
                     y.position = c(26,NA,25))+
  ggtitle("lactobacilli")

p_Lactobacillus 

sheet[["Lactobacillus"]] <- stat_bac
#Lactococcus##########################################################################
ps0 <- tax_glom(subset_taxa(ps.rel, Genus=="Lactococcus"), "Genus", NArm = FALSE)
df <- psmelt(ps0)
df[df==""] <- NA
#Select last non-empty taxonomic rank
df$tax <- apply(df, 1, function(x) tail(stats::na.omit(x), 1))
#Add "unknown" to taxonomy columns where genus is not classified 
df$tax[is.na(df$Genus)] <- paste0(df$tax[is.na(df$Genus)],"unknown genus")
df_bac<-df
df_bac$Liver_response <- factor(df_bac$Liver_response, levels = liver_group)
df_bac$Abundance <- df_bac$Abundance*100

#stat
stat_bac<- stat_single_compare(df_bac)
df_bac$Liver_response <- factor(df_bac$Liver_response, levels = liver_group)


p_Lactococcus <- ggplot_boxplot_abundance(df_bac)+
  stat_pvalue_manual(stat_bac,label = "p.adj.signif", tip.length = 0, size = 4,
                     y.position = c(NA,NA,48))+
  ggtitle("Lactococcus")

p_Lactococcus

sheet[["Lactococcus"]] <- stat_bac

#Oscillospira##########################################################################
ps0 <- tax_glom(subset_taxa(ps.rel, Genus=="Oscillospira"), "Genus", NArm = FALSE)
df <- psmelt(ps0)
df[df==""] <- NA
#Select last non-empty taxonomic rank
df$tax <- apply(df, 1, function(x) tail(stats::na.omit(x), 1))
#Add "unknown" to taxonomy columns where genus is not classified 
df$tax[is.na(df$Genus)] <- paste0(df$tax[is.na(df$Genus)],"unknown genus")
df_bac<-df
df_bac$Liver_response <- factor(df_bac$Liver_response, levels = liver_group)
df_bac$Abundance <- df_bac$Abundance*100

#stat
stat_bac<- stat_single_compare(df_bac)
df_bac$Liver_response <- factor(df_bac$Liver_response, levels = liver_group)


p_Oscillospira <- ggplot_boxplot_abundance(df_bac)+
  stat_pvalue_manual(stat_bac,label = "p.adj.signif", tip.length = 0, size = 3,
                     y.position = c(4,NA,3.8))+
  ggtitle("Oscillospira")

p_Oscillospira 

sheet[["Oscillospira"]] <- stat_bac

###############################################################################
write_xlsx(sheet, "Bacteriome/stat_result/sub/single_genus_abundance.xlsx")
all_variables <- ls()
variable_to_keep <- c("p_Allobaculum","p_Oscillospira","p_Bacteroides","p_Bifidobacterium",
                      "p_Clostridium","p_Dehalobacterium","p_Lactobacillus","p_Lactococcus")
# sub_supplementary_fig3<-ggarrange(p_Allobaculum,p_Bacteroides,p_Bifidobacterium,p_Clostridium,
#                                           p_Dehalobacterium,p_Lactobacillus,p_Lactococcus,p_Oscillospira,
#                                           ncol = 4,
#                                           nrow = 2,
#                                           labels = c("A","B","C","D","E","F","G","H"),
#                                           common.legend = TRUE)

