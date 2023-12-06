setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("Analysis0_file_loading_and_prep.R")

ps <- subset_samples(PSB, Sample_time_point %in% c("Termination"))
ps.rel <- transform_sample_counts(ps, function(x) x / sum(x))

#####setting####################################################################
df <- data.frame(sample_data(ps.rel))

# count numbers
chp_n <-  count(df$Treatment == "FVT_ChP")[[2,2]]
sdt_n <-  count(df$Treatment == "FVT_SDT")[[2,2]]
pty_n <-  count(df$Treatment == "FVT_PyT")[[2,2]]
fvt_n <-  count(df$Treatment == "Unmodified_FVT")[[2,2]]
obese_n <-  count(df$Treatment == "Obese_control")[[2,2]]
lean_n <-  count(df$Treatment == "Lean_control")[[2,2]]

counts<-c(glue("(N={chp_n})"),
                glue("(N={sdt_n})"),
                glue("(N={pty_n})"),
                glue("(N={fvt_n})"),
                glue("(N={obese_n})"),
                glue("(N={lean_n})"))
#statistical function
stat_single_compare <- function(df){
  stat_bac<- df %>%
    wilcox_test(Abundance~Treatment,
                comparisons = Obse_compare,
                p.adjust.method = "fdr",
                paired = FALSE, 
                alternative = "two.sided",
                detailed = TRUE) 
  stat_bac
  stat_bac_neg<- df %>%
    wilcox_test(Abundance~Treatment,
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
ggplot_boxplot_abundance <- function(df)
{
  ggplot(df_bac,aes(x=Treatment,y=Abundance)) +
    stat_boxplot(geom ='errorbar', linetype=1, width=0.2) +
    geom_boxplot(outlier.shape = NA,alpha = 1,aes(fill=Treatment), coef=1.5) +
    geom_jitter(show.legend=FALSE, width=0.25, shape=21, fill="black") +
    # stat_summary(fun=mean, show.legend=FALSE, geom="crossbar", linetype=3,
    #              color="black",width=0.75, size=0.4)+
    scale_x_discrete(breaks=groups,
                     labels=counts) +
    scale_fill_manual(values = cols) +
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
df_bac$Treatment <- factor(df_bac$Treatment, levels = groups)
df_bac$Abundance <- df_bac$Abundance*100

#stat
stat_bac<- stat_single_compare(df_bac)
df_bac$Treatment <- factor(df_bac$Treatment,levels = groups)


p_Allobaculum <- ggplot_boxplot_abundance(df_bac)+
  stat_pvalue_manual(stat_bac,label = "p.adj.signif", tip.length = 0, size = 4,
                     y.position = c(NA,NA,15,17,19))+
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
df_bac$Treatment <- factor(df_bac$Treatment, levels = groups)
df_bac$Abundance <- df_bac$Abundance*100

#stat
stat_bac<- stat_single_compare(df_bac)
df_bac$Treatment <- factor(df_bac$Treatment,levels = groups)


p_Bacteroides <- ggplot_boxplot_abundance(df_bac)+
  stat_pvalue_manual(stat_bac,label = "p.adj", tip.length = 0, size = 3,
                     y.position = c(15,17,NA,NA,NA))+
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
df_bac$Treatment <- factor(df_bac$Treatment, levels = groups)
df_bac$Abundance <- df_bac$Abundance*100

#stat
stat_bac<- df_bac %>%
  wilcox_test(Abundance~Treatment,
              comparisons = neg_compare,
              p.adjust.method = "fdr",
              paired = FALSE,
              alternative = "two.sided",
              detailed = TRUE) 
stat_bac

stat.test <- stat_bac
df_bac$Treatment <- factor(df_bac$Treatment,levels = groups)


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
df_bac$Treatment <- factor(df_bac$Treatment, levels = groups)
df_bac$Abundance <- df_bac$Abundance*100

#stat
stat_bac<- stat_single_compare(df_bac)
df_bac$Treatment <- factor(df_bac$Treatment,levels = groups)


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
df_bac$Treatment <- factor(df_bac$Treatment, levels = groups)
df_bac$Abundance <- df_bac$Abundance*100

#stat
stat_bac<- stat_single_compare(df_bac)
df_bac$Treatment <- factor(df_bac$Treatment,levels = groups)


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
df_bac$Treatment <- factor(df_bac$Treatment, levels = groups)
df_bac$Abundance <- df_bac$Abundance*100

#stat
stat_bac<- stat_single_compare(df_bac)
df_bac$Treatment <- factor(df_bac$Treatment,levels = groups)


p_Lactobacillus <- ggplot_boxplot_abundance(df_bac)+
  stat_pvalue_manual(stat_bac,label = "p.adj.signif", tip.length = 0, size = 3,
                     y.position = c(40,NA,NA,NA,25))+
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
df_bac$Treatment <- factor(df_bac$Treatment, levels = groups)
df_bac$Abundance <- df_bac$Abundance*100

#stat
stat_bac<- stat_single_compare(df_bac)
df_bac$Treatment <- factor(df_bac$Treatment,levels = groups)


p_Lactococcus <- ggplot_boxplot_abundance(df_bac)+
  stat_pvalue_manual(stat_bac,label = "p.adj.signif", tip.length = 0, size = 4,
                     y.position = c(NA,NA,NA,NA,48))+
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
df_bac$Treatment <- factor(df_bac$Treatment, levels = groups)
df_bac$Abundance <- df_bac$Abundance*100

#stat
stat_bac<- stat_single_compare(df_bac)
df_bac$Treatment <- factor(df_bac$Treatment,levels = groups)


p_Oscillospira <- ggplot_boxplot_abundance(df_bac)+
  stat_pvalue_manual(stat_bac,label = "p.adj", tip.length = 0, size = 3,
                     y.position = c(NA,NA,NA,NA,3.8))+
  ggtitle("Oscillospira")

p_Oscillospira 

sheet[["Oscillospira"]] <- stat_bac

###############################################################################
write_xlsx(sheet, "Bacteriome/stat_result/single_genus_abundance.xlsx")
all_variables <- ls()
variable_to_keep <- c("p_Allobaculum","p_Oscillospira","p_Bacteroides","p_Bifidobacterium",
                      "p_Clostridium","p_Dehalobacterium","p_Lactobacillus","p_Lactococcus")

p_rel_abundance<- ggarrange(p_Allobaculum,p_Bacteroides,p_Bifidobacterium,p_Clostridium,
          p_Dehalobacterium,p_Lactobacillus,p_Lactococcus,p_Oscillospira,
          ncol = 4,
          nrow = 2,
          labels = c("F","G","H","I","J","K","L","M"),
          common.legend = TRUE,
          legend = "top")


# ggarrange(p_supplementary_fig1,p_rel_abundance,
#           ncol = 1,
#           nrow = 2,
#           heights = c(5,4),
#           common.legend = TRUE,
#           legend = "top")
