

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##############################Load sequencing data#######################################

##Import files

#Bacteria

source("2_CSS_phyloseq_bac.R")


############Set data categories - and order



PSB.CSS = prune_samples(sample_sums(PSB) >= 2000, PSB.CSS)
PSB = prune_samples(sample_sums(PSB) >= 2000, PSB)





PSB <- subset_samples(PSB, Sample_time_point %in% c("Arrival","Before_1st_FVT","1w_after_FVT","Termination"))
PSB.CSS <- subset_samples(PSB.CSS, Sample_time_point %in% c("Arrival","Before_1st_FVT","1w_after_FVT","Termination"))
PSB.CSS@sam_data$Sample_time_point <- fct_relevel(PSB.CSS@sam_data$Sample_time_point, "Arrival","Before_1st_FVT","1w_after_FVT","Termination")
PSB@sam_data$Sample_time_point <- fct_relevel(PSB@sam_data$Sample_time_point, "Arrival","Before_1st_FVT","1w_after_FVT","Termination")




# saveRDS(PSB, file = "data/corr_data/ps_bacterium.rds")

