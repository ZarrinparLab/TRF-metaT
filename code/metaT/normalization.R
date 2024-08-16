setwd("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_m_results/")

library(tidyverse)
library(data.table)
####################################################
#PF04563.18 is RPOB


#read in TPM table
pfamTPM<-fread("pfam/pfam-TPM_clean_noNT.tsv")%>%
  column_to_rownames("FeatureID")

for (col in colnames(pfamTPM)) {
  pfamTPM[[col]] <- pfamTPM[[col]] / pfamTPM[[col]][849]
}

pfamTPM<-pfamTPM%>%rownames_to_column("FeatureID")

write.table(pfamTPM,"pfam/pfam-TPM_clean_noNT_normRPOB.txt",sep = "\t",row.names = FALSE, quote=FALSE)


####################################################
#PF04563.18 is RPOB


#read in TPM table
pfamTPM<-fread("g-diversity-core-metrics114K/rarefied_table/pfam-144k_TPM.tsv")%>%
  column_to_rownames("#FeatureID")%>%dplyr::select(-Name)

for (col in colnames(pfamTPM)) {
  pfamTPM[[col]] <- pfamTPM[[col]] / pfamTPM[[col]][239]
}

pfamTPM<-pfamTPM%>%rownames_to_column("FeatureID")

write.table(pfamTPM,"g-diversity-core-metrics114K/rarefied_table/pfam-144k_clean_noNT_normRPOB.txt",sep = "\t",row.names = FALSE, quote=FALSE)



