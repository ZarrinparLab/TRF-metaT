setwd(setwd("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metagenomic/woltka2_results/filtered_metaG"))

library(tidyverse)
library(data.table)
####################################################
#PF04563.18 is RPOB


#read in TPM table
pfamTPM<-fread("pfam_TPM/pfam_clean_noNT_TPM.txt")%>%
  column_to_rownames("#OTU ID")

for (col in colnames(pfamTPM)) {
  pfamTPM[[col]] <- pfamTPM[[col]] / pfamTPM[[col]][491]
}

pfamTPM<-pfamTPM%>%rownames_to_column("FeatureID")

write.table(pfamTPM,"pfam_TPM/pfam_clean_noNT_TPM_normRPOB.txt",sep = "\t",row.names = FALSE, quote=FALSE)



