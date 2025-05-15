setwd("/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/TRF-metaT/")

library(tidyverse)
library(data.table)

####################################################
norm_data<-"data/pfam_metaG/pfam-TPM_clean_noNT.txt"
dat_path<-"data/pfam_metaG/"
####################################################
#PF04563.18 is RPOB

#read in TPM table
pfamTPM<-fread(norm_data)%>%
  column_to_rownames("#OTU ID")

for (col in colnames(pfamTPM)) {
  pfamTPM[[col]] <- pfamTPM[[col]] / pfamTPM[[col]][491]
}

pfamTPM<-pfamTPM%>%rownames_to_column("FeatureID")

write.table(pfamTPM,paste0(dat_path,"pfam_clean_noNT_TPM_normRPOB.txt"),sep = "\t",row.names = FALSE, quote=FALSE)



