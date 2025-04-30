setwd("~/Notebooks/sfloresr/TRF-metaT/")

library(data.table)
library(stringr)
library(tibble)
library(dplyr)
library(vegan)

#####################################################################################################
#paths
dat_library<-"data/ENB_invivo_metab/GNPS_fecal_mtb/GNPS_outputs/library/merged_results_with_gnps.tsv"
dat_feattab<-"data/ENB_invivo_metab/GNPS_fecal_mtb/GNPS_outputs/clustering/featuretable_reformated.csv"
dat_path<-"data/ENB_invivo_metab/GNPS_fecal_mtb/GNPS_outputs/"
#####################################################################################################
#filtering for just bile acids (BA)

annot<-fread(dat_library)%>%
  filter(Organism=="BILELIB19" | Organism=="GNPS-BILE-ACID-MODIFICATIONS" )
write.table(annot,paste0(dat_path,"library/merged_results_with_gnps_justBA.tsv"),sep = "\t",row.names = FALSE, quote=FALSE)  

mtb<-fread(dat_feattab)%>%
  filter(`row ID` %in% annot$`#Scan#`)
           write_csv(mtb, paste0(dat_path,"clustering/featuretable_reformated_justBA.csv"))
         
#####################################################################################################
#Cleanup BA feature quant table from MZmine 

norm <- read.csv(paste0(dat_path,"clustering/featuretable_reformated_justBA.csv"))%>%
  dplyr::select(-"Unnamed..162")
norm <- norm %>% rename_with(~str_replace_all(., ".mzML.Peak.area", ".mzML Peak area"))
norm <- norm %>% dplyr::select(-(row.m.z:ccs), -(annotation.network.number:neutral.M.mass))
norm <- norm %>% mutate(across(row.ID:correlation.group.ID, as.character))
norm <- norm %>% rename_with(~str_replace_all(., "ut.", "ut-"))

## Merge ions from the same molecule (iimn)
iimn <- norm %>%
  group_by(correlation.group.ID) %>%
  summarise(across(where(is.numeric), ~sum(., na.rm = TRUE)),
            collapserowID = paste(row.ID, collapse = '|')) %>%
  filter(!is.na(correlation.group.ID)) %>%
  dplyr::rename(row.ID = correlation.group.ID)%>%
  dplyr::select(row.ID,collapserowID,everything())
labels <- norm %>% dplyr::select(row.ID, correlation.group.ID)
iimn$row.ID <- paste0(iimn$row.ID, "_i")
norm_iimn <- norm %>% filter(is.na(correlation.group.ID)) %>% dplyr::select(-correlation.group.ID)%>%
  mutate(collapserowID=NA)%>%
  dplyr::select(row.ID,collapserowID,everything())
norm_iimn <- rbind(norm_iimn, iimn)%>%
  dplyr::select(-collapserowID)

write_csv(norm_iimn, paste0(dat_path,"clustering/featuretable_reformated_justBA_summcollapse.csv"))

norm_iimn_fornorm<-norm_iimn%>%
  column_to_rownames("row.ID")%>%as.matrix()

norm_iimn_rclr <- norm_iimn_fornorm %>% decostand(method = "rclr")%>%
  as.data.frame()%>%rownames_to_column("row.ID")
write_csv(norm_iimn_rclr, paste0(dat_path,"clustering/featuretable_reformated_justBA_summcollapse_cln_rclr.csv"))

## Read in FBMN library IDs from GNPS
library_ID <- read.delim(paste0(dat_path,"library/merged_results_with_gnps_justBA.tsv"))%>%
  dplyr::rename(row.ID=X.Scan.)
library_ID$row.ID<- as.character(library_ID$row.ID)
library_ID <- right_join(labels, library_ID, by = "row.ID", multiple = "all")

result_list_iin <- library_ID %>%
  filter(!is.na(correlation.group.ID)) %>%
  group_by(correlation.group.ID) %>%
  summarise(
    mn_mz = mean(SpecMZ, na.rm = TRUE),
    across(-SpecMZ, ~paste(., collapse = "|")))%>%
  dplyr::rename(FeatureID=correlation.group.ID)

result_list_iin$FeatureID <- paste0(result_list_iin$FeatureID, "_i")

result_list <- library_ID %>%
  filter(is.na(correlation.group.ID)) %>%
  dplyr::select(c(1,12,2:11,13:49))%>%
  dplyr::rename(FeatureID=row.ID,x=correlation.group.ID, mn_mz=SpecMZ)%>%
  dplyr::rename(row.ID=x)

result_list<-rbind(result_list_iin,result_list)
write.table(result_list,paste0(dat_path,"library/merged_results_with_gnps_justBA_collapse.tsv"),sep = "\t",row.names = FALSE, quote=FALSE)  
      