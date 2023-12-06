

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##############################Load sequencing data#######################################

##Import files


source("3_CSS_phyloseq_vir.R")



#take the filtrate out
PSV.Filtrate <- subset_samples(PSV.no.realm, ProjectID %in% c("FVT_ChP", "FVT_SDT", "FVT_PyT", "Unmodified_FVT"))
PSV.Filtrate.CSS <- subset_samples(PSV.no.Realm.CSS, ProjectID %in% c("FVT_ChP", "FVT_SDT", "FVT_PyT", "Unmodified_FVT"))
PSV.Filtrate.HOST <- subset_samples(PSV.HOST, ProjectID %in% c("FVT_ChP", "FVT_SDT", "FVT_PyT", "Unmodified_FVT"))



#only recipents sample
PSV.no.realm <- subset_samples(PSV.no.realm, Sample_time_point %in% c("Arrival","1w_after_FVT","Termination"))
PSV.no.Realm.CSS <- subset_samples(PSV.no.Realm.CSS, Sample_time_point %in% c("Arrival","1w_after_FVT","Termination"))
PSV.HOST <- subset_samples(PSV.HOST, Sample_time_point %in% c("Arrival","1w_after_FVT","Termination"))



# saveRDS(PSV.no.realm,file = "data/corr_data/ps_virome.rds")
