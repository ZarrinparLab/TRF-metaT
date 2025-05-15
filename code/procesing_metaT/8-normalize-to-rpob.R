setwd("/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/TRF-metaT/")

library(tidyverse)
library(data.table)
####################################################
norm_data1<-"data/pfam_metaT/pfam-TPM_clean_noNT.tsv"
dat_path1<-"data/pfam_metaT/"
norm_data2<-"data/DE_analysis/g-diversity-core-metrics114K/rarefied_table/pfam-144k_TPM.tsv"
dat_path2<-"data/DE_analysis/g-diversity-core-metrics114K/rarefied_table/"
####################################################
#for non-rarefied metaT data
#PF04563.18 is RPOB

#read in TPM table
pfamTPM<-fread(norm_data1)%>%
  column_to_rownames("FeatureID")

for (col in colnames(pfamTPM)) {
  pfamTPM[[col]] <- pfamTPM[[col]] / pfamTPM[[col]][849]
}

pfamTPM<-pfamTPM%>%rownames_to_column("FeatureID")

write.table(pfamTPM,paste0(dat_path1,"pfam-TPM_clean_noNT_normRPOB.txt"),sep = "\t",row.names = FALSE, quote=FALSE)


####################################################
#for rarefied metaT data
#PF04563.18 is RPOB

#read in TPM table
pfamTPM<-fread(norm_data2)%>%
  column_to_rownames("#FeatureID")%>%dplyr::select(-Name)

for (col in colnames(pfamTPM)) {
  pfamTPM[[col]] <- pfamTPM[[col]] / pfamTPM[[col]][239]
}

pfamTPM<-pfamTPM%>%rownames_to_column("FeatureID")

write.table(pfamTPM,paste0(dat_path2,"pfam-144k_clean_noNT_normRPOB.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

