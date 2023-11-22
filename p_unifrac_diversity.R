
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("Analysis0_file_loading_and_prep.R")

set.seed(19950915)
ps.arrival <- subset_samples(PSB.CSS, Sample_time_point %in% "Arrival")
ps.beforefvt <- subset_samples(PSB.CSS, Sample_time_point %in% "Before_1st_FVT")
ps.afterfvt <- subset_samples(PSB.CSS, Sample_time_point %in% "1w_after_FVT")
ps.termination <- subset_samples(PSB.CSS, Sample_time_point %in% "Termination")

desired_order <- c("FVT_ChP", "FVT_SDT", "FVT_PyT", "Unmodified_FVT", "Lean_control")
Method="unifrac"
#############################stat###########################################
#1 Arrival


ps <- ps.arrival
GP.ord <- ordinate(ps, "PCoA", Method)
new_colnames <- gsub("Axis.", "PCoA", colnames(GP.ord$vectors))
colnames(GP.ord$vectors) <- new_colnames
bray.PSB <- phyloseq::distance(ps, method = Method)
# make a data frame from the sample_data
sampledf.PSB <- data.frame(sample_data(ps))
adonis_sub <- adonis2(bray.PSB ~ Treatment, method = Method, data = sampledf.PSB, permutations = 999)

adonis_R2 <- adonis_sub[1,3]
adonis_p <- adonis_sub[1,5]

ps@sam_data$Treatment <- factor(ps@sam_data$Treatment, levels = groups)
p_arrival <- phyloseq::plot_ordination(ps, GP.ord, color="Treatment") + 
  stat_ellipse(geom = "polygon", level = 0.95, fill = NA, linewidth = 1,show.legend = FALSE) +
  geom_point(alpha=.8, size=2,aes(color=Treatment)) +
  scale_color_manual(values= cols)+
  ggtitle(paste(" Arrival\n",
                "R2=", format(adonis_R2[1],digits = 2),", P=",format(adonis_p[1], digits = 1),"\n")) +
  theme_classic()+
  mytheme_beta 



#Calculate pairwise p value result
metadata <- data.frame(sample_data(ps))
cbn <- combn(x=unique(metadata$Treatment), m = 2)
cbn <- cbn[,1:5]

p <- c()
r <- c()
for(i in 1:ncol(cbn)){
  ps.subs <- subset_samples(ps, Treatment %in% cbn[,i])
  metadata_sub <- data.frame(sample_data(ps.subs))
  permanova_pairwise <- adonis2(phyloseq::distance(ps.subs, method = Method) ~ Treatment, 
                                data = metadata_sub)
  p <- c(p, permanova_pairwise$`Pr(>F)`[1])
  r <- c(r, permanova_pairwise$R2[1])
}

p_control <- p[1]
p_treatment <- p[2:5]
p_treatment.adj <- round(p.adjust(p_treatment, method = "fdr"),digits=3)

p_treatment.adj <- c("_",p_treatment.adj)
t_arrival <- data.frame(Group1=c( "Obese-control","Obese-control","Obese-control","Obese-control","Obese-control"),
                        Group2=c("Lean-control","FVT-SDT","FVT-ChP","Unmodified-FVT","FVT-PyT"),
                        R2= round(r,digits = 3),
                        p=round(p,digits = 3),
                        p.adj=p_treatment.adj)







p_arrival_stat<-ggplot()+
  theme_void() +  # Remove background
  theme(
    axis.title = element_blank(),  # Remove axis titles
    axis.text = element_blank(),   # Remove axis text
    axis.ticks = element_blank()   # Remove axis ticks
  )+
  annotate(geom="table",
           x=0.5,
           y=0.5,
           label = list(t_arrival))
p_arrival_stat
adonis_sub_arrival <- adonis_sub

#############################stat###########################################
#2  before fvt


ps <- ps.beforefvt

GP.ord <- ordinate(ps, "PCoA", Method)

new_colnames <- gsub("Axis.", "PCoA", colnames(GP.ord$vectors))
colnames(GP.ord$vectors) <- new_colnames

bray.PSB <- phyloseq::distance(ps, method = Method)
# make a data frame from the sample_data
sampledf.PSB <- data.frame(sample_data(ps))

adonis_sub <- adonis2(bray.PSB ~ Treatment, method = Method, data = sampledf.PSB, permutations = 999)
adonis_sub
adonis_R2 <- adonis_sub[1,3]
adonis_p <- adonis_sub[1,5]

ps@sam_data$Treatment <- factor(ps@sam_data$Treatment, levels = groups)
p_beforefvt <- phyloseq::plot_ordination(ps, GP.ord, color="Treatment") + 
  stat_ellipse(geom = "polygon", level = 0.95, fill = NA, linewidth = 1,show.legend = FALSE) +
  geom_point(alpha=.8, size=2,aes(color=Treatment)) +
  scale_color_manual(values= cols)+
  ggtitle(paste( " Before 1st FVT\n",
                 "R2=", format(adonis_R2[1],digits = 2),", P=",format(adonis_p[1], digits = 1),"\n")) +
  theme_classic()+
  mytheme_beta 



#Calculate pairwise p value result
metadata <- data.frame(sample_data(ps))
cbn <- combn(x=unique(metadata$Treatment), m = 2)
cbn <- cbn[,1:5]

p <- c()
r <- c()
for(i in 1:ncol(cbn)){
  ps.subs <- subset_samples(ps, Treatment %in% cbn[,i])
  metadata_sub <- data.frame(sample_data(ps.subs))
  permanova_pairwise <- adonis2(phyloseq::distance(ps.subs, method = Method) ~ Treatment, 
                                data = metadata_sub)
  p <- c(p, permanova_pairwise$`Pr(>F)`[1])
  r <- c(r, permanova_pairwise$R2[1])
}

p_control <- p[1]
p_treatment <- p[2:5]
p_treatment.adj <- round(p.adjust(p_treatment, method = "fdr"),digits=3)

p_treatment.adj <- c("_",p_treatment.adj)
t_beforefvt <- data.frame(Group1=c( "Obese-control","Obese-control","Obese-control","Obese-control","Obese-control"),
                        Group2=c("Lean-control","FVT-SDT","FVT-ChP","Unmodified-FVT","FVT-PyT"),
                        R2= round(r,digits = 3),
                        p=round(p,digits = 3),
                        p.adj=p_treatment.adj)








p_beforefvt_stat<-ggplot()+
  theme_void() +  # Remove background
  theme(
    axis.title = element_blank(),  # Remove axis titles
    axis.text = element_blank(),   # Remove axis text
    axis.ticks = element_blank()   # Remove axis ticks
  )+
  annotate(geom="table",
           x=0.5,
           y=0.5,
           label = list(t_beforefvt))
p_beforefvt_stat

adonis_sub_beforefvt <- adonis_sub

#############################stat###########################################
#3  after fvt


ps <- ps.afterfvt

GP.ord <- ordinate(ps, "PCoA", Method)

new_colnames <- gsub("Axis.", "PCoA", colnames(GP.ord$vectors))
colnames(GP.ord$vectors) <- new_colnames

bray.PSB <- phyloseq::distance(ps, method = Method)
# make a data frame from the sample_data
sampledf.PSB <- data.frame(sample_data(ps))

adonis_sub <- adonis2(bray.PSB ~ Treatment, method = Method, data = sampledf.PSB, permutations = 999)

adonis_R2 <- adonis_sub[1,3]
adonis_p <- adonis_sub[1,5]

ps@sam_data$Treatment <- factor(ps@sam_data$Treatment, levels = groups)
p_afterfvt <- phyloseq::plot_ordination(ps, GP.ord, color="Treatment") + 
  stat_ellipse(geom = "polygon", level = 0.95, fill = NA, linewidth = 1,show.legend = FALSE) +
  geom_point(alpha=.8, size=2,aes(color=Treatment)) +
  scale_color_manual(values= cols)+
  ggtitle(paste(" 1 week after 2nd FVT\n",
                "R2=", format(adonis_R2[1],digits = 2),", P=",format(adonis_p[1], digits = 1),"\n")) +
  theme_classic()+
  mytheme_beta 



#Calculate pairwise p value result
metadata <- data.frame(sample_data(ps))
cbn <- combn(x=unique(metadata$Treatment), m = 2)
cbn <- cbn[,1:5]

p <- c()
r <- c()
for(i in 1:ncol(cbn)){
  ps.subs <- subset_samples(ps, Treatment %in% cbn[,i])
  metadata_sub <- data.frame(sample_data(ps.subs))
  permanova_pairwise <- adonis2(phyloseq::distance(ps.subs, method = Method) ~ Treatment, 
                                data = metadata_sub)
  p <- c(p, permanova_pairwise$`Pr(>F)`[1])
  r <- c(r, permanova_pairwise$R2[1])
}

p_control <- p[1]
p_treatment <- p[2:5]
p_treatment.adj <- round(p.adjust(p_treatment, method = "fdr"),digits=3)

p_treatment.adj <- c("_",p_treatment.adj)
t_afterfvt <- data.frame(Group1=c( "Obese-control","Obese-control","Obese-control","Obese-control","Obese-control"),
                        Group2=c("Lean-control","FVT-SDT","FVT-ChP","Unmodified-FVT","FVT-PyT"),
                        R2= round(r,digits = 3),
                        p=round(p,digits = 3),
                        p.adj=p_treatment.adj)








p_afterfvt_stat<-ggplot()+
  theme_void() +  # Remove background
  theme(
    axis.title = element_blank(),  # Remove axis titles
    axis.text = element_blank(),   # Remove axis text
    axis.ticks = element_blank()   # Remove axis ticks
  ) + 
  annotate(geom="table",
           x=0.5,
           y=0.5,
           label = list(t_afterfvt))
p_afterfvt_stat


adonis_sub_afterfvt <-adonis_sub
#############################stat###########################################
#4 termination


ps <- ps.termination

GP.ord <- ordinate(ps, "PCoA", Method)

new_colnames <- gsub("Axis.", "PCoA", colnames(GP.ord$vectors))
colnames(GP.ord$vectors) <- new_colnames

bray.PSB <- phyloseq::distance(ps, method = Method)
# make a data frame from the sample_data
sampledf.PSB <- data.frame(sample_data(ps))

adonis_sub <- adonis2(bray.PSB ~ Treatment, method = Method, data = sampledf.PSB, permutations = 999)

adonis_R2 <- adonis_sub[1,3]
adonis_p <- adonis_sub[1,5]

ps@sam_data$Treatment <- factor(ps@sam_data$Treatment, levels = groups)
p_termination <- phyloseq::plot_ordination(ps, GP.ord, color="Treatment") + 
  stat_ellipse(geom = "polygon", level = 0.95, fill = NA, linewidth = 1,show.legend = FALSE) +
  geom_point(alpha=.8, size=2,aes(color=Treatment)) +
  scale_color_manual(values= cols)+
  ggtitle(paste( " Termination\n",
                 "R2=", format(adonis_R2[1],digits = 2),", P=",format(adonis_p[1], digits = 1),"\n")) +
  theme_classic()+
  mytheme_beta 



#Calculate pairwise p value result
metadata <- data.frame(sample_data(ps))
cbn <- combn(x=unique(metadata$Treatment), m = 2)
cbn <- cbn[,1:5]

p <- c()
r <- c()
for(i in 1:ncol(cbn)){
  ps.subs <- subset_samples(ps, Treatment %in% cbn[,i])
  metadata_sub <- data.frame(sample_data(ps.subs))
  permanova_pairwise <- adonis2(phyloseq::distance(ps.subs, method = Method) ~ Treatment, 
                                data = metadata_sub)
  p <- c(p, permanova_pairwise$`Pr(>F)`[1])
  r <- c(r, permanova_pairwise$R2[1])
}

p_control <- p[1]
p_treatment <- p[2:5]
p_treatment.adj <- round(p.adjust(p_treatment, method = "fdr"),digits=3)

p_treatment.adj <- c("_",p_treatment.adj)
t_termination <- data.frame(Group1=c( "Obese-control","Obese-control","Obese-control","Obese-control","Obese-control"),
                        Group2=c("Lean-control","FVT-SDT","FVT-ChP","Unmodified-FVT","FVT-PyT"),
                        R2= round(r,digits = 3),
                        p=round(p,digits = 3),
                        p.adj=p_treatment.adj)






p_termination_stat<-ggplot()+
  theme_void() +  # Remove background
  theme(
    axis.title = element_blank(),  # Remove axis titles
    axis.text = element_blank(),   # Remove axis text
    axis.ticks = element_blank(),
    # Remove axis ticks
  )+
  annotate(geom="table",
           x=0.5,
           y=0.5,
           label = list(t_termination))
p_termination_stat

adonis_sub_termination <- adonis_sub
################################################################################

# p_unifac_bac<- ggarrange(p_arrival,
#                        p_beforefvt,
#                        p_afterfvt,
#                        p_termination,
#                        nrow=1, ncol = 4,
#                        font.label = list(size = 20),
#                        common.legend = TRUE,
#                        legend = "bottom")
p_unifac_stat_bac <- ggarrange(p_arrival,
                             p_beforefvt,
                             p_afterfvt,
                             p_termination,
                             # p_arrival_stat,
                             # p_beforefvt_stat,
                             # p_afterfvt_stat,
                             # p_termination_stat,
                             nrow=1, ncol = 4,
                             font.label = list(size = 10),
                             common.legend = TRUE,
                             legend = "bottom")
p_unifac_stat_bac

################################################################################
#  sheets <- list(
#  "adonis_sub_arrival" =  adonis_sub_arrival,
# "adonis_sub_beforefvt" = adonis_sub_beforefvt,
#  "adonis_sub_afterfvt" =  adonis_sub_afterfvt,
#  "adonis_sub_termination" = adonis_sub_termination
#  )
# 
#  write_xlsx(sheets, "stat_result/unifac_permanova_all.xlsx")
# 
#  sheet2 <- list(
#    "adonis_sub_arrival" =  t_arrival,
#    "adonis_sub_beforefvt" = t_beforefvt,
#    "adonis_sub_afterfvt" =  t_afterfvt,
#    "adonis_sub_termination" = t_termination
#  )
# 
#  write_xlsx(sheet2, "stat_result/unifac_permanova_pairs.xlsx")


# all_variables <- ls()
# variable_to_keep <- c("p_unifac_bac","p_unifac_stat_bac")
# variables_to_remove <- setdiff(all_variables, variable_to_keep)
# rm(list = variables_to_remove)
# rm(variable_to_keep,variables_to_remove,all_variables)
