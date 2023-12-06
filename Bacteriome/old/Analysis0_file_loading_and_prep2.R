###################DIO study data analysis###################
library(tidyverse)
library(ggpubr)
library(ggsci)
library(rstatix)
library(ampvis2)
library(plyr)
library(cowplot)
library(RVAideMemoire)
library(data.table)
library(DESeq2)
library(viridis)
library(vegan)
library(directlabels)
library(usdm)
library("varhandle") 
#library(Rarefy)
#library(EnhancedVolcano)
#install.packages("remotes")
#remotes::install_github("schuyler-smith/phylosmith")
library(remotes)
#library(phylosmith)

#Load phyloseq to amvis2 converter
devtools::source_gist("8d0ca4206a66be7ff6d76fc4ab8e66c6")

#Coulour pallettes
col_fil <- pal_jco("default")(10)

col_scale <- scale_color_jco()

#Set working directory to script directory

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##############################Load sequencing data#######################################

##Import files

#Bacteria

source("Phyloseq_import_bac_vir.R")

source("CSS_phyloseq_bac_vir.R")


#######################################Remove low read samples###############################################################

PSB.CSS = prune_samples(sample_sums(PSB) >= 100, PSB.CSS)
PSB = prune_samples(sample_sums(PSB) >= 100, PSB)

#####################Choose the samples e.g. timepoint like Termination######
##Ensure that all data is not overwrited
PSB.all <- PSB  
PSB.CSS.all <- PSB.CSS

PSB <- subset_samples(PSB.all,Sample_time_point == "Termination")
PSB.CSS <- subset_samples(PSB.CSS.all,Sample_time_point == "Termination")


#You can also do it more selective by the Remove column in the mapping file. 
#PSB <- subset_samples(PSB.all,Remove == "s")
#PSB.CSS <- subset_samples(PSB.CSS.all,Remove == "s")



