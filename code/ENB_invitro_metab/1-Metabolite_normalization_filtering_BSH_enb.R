setwd("~/Notebooks/sfloresr/TRF-metaT/")

library(data.table)
library(dplyr)
library(stringr)
library(vegan)
#####################################################################################################
#paths
dat_library<-"data/ENB_invitro_metab/GNPS/GNPS_outputs/networking/library-results-merged_results_with_gnps.tsv"
dat_feattab<-"data/ENB_invitro_metab/GNPS/GNPS_outputs/clustering/featuretable_reformated.csv"
dat_path<-"data/ENB_invitro_metab/GNPS/GNPS_outputs/"

#####################################################################################################
#combine the additional plate that was run to table

##original run
annot<-fread(dat_library)%>%
  filter(Organism=="BILELIB19" | Organism=="GNPS-BILE-ACID-MODIFICATIONS" )
write.table(annot,paste0(dat_path,"networking/library-results-merged_results_with_gnps_justBA.tsv"),sep = "\t",row.names = FALSE, quote=FALSE)  

mtb<-fread(dat_feattab)%>%
  filter(`row ID` %in% annot$`#Scan#`)
write_csv(mtb, paste0(dat_path,"clustering/featuretable_reformated_justBA.csv"))

## Cleanup feature quant table from MZmine ##########################################################
norm <- read.csv(paste0(dat_path,"clustering/featuretable_reformated_justBA.csv"))%>%
  dplyr::select(-"Unnamed..655")
norm <- norm %>% rename_with(~str_replace_all(., ".mzML.Peak.area", ".mzML Peak area"))
norm <- norm %>% rename_with(~str_replace_all(., "EcAZ.1.cat", "EcAZ-1-cat"))
norm <- norm %>% rename_with(~str_replace_all(., "AZ.52", "AZ-52"))
norm <- norm %>% rename_with(~str_replace_all(., "LCAG.95", "LCAG-95"))
norm <- norm %>% dplyr::select(-(row.m.z:ccs), -(annotation.network.number:neutral.M.mass))
norm <- norm %>% mutate(across(row.ID:correlation.group.ID, as.character))

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
write_csv(norm_iimn, paste0(dat_path,"clustering/featuretable_reformated_justBA_summcollapse_cln.csv"))

norm_iimn_fornorm<-norm_iimn%>%
  column_to_rownames("row.ID")%>%as.matrix()

norm_iimn_rclr <- norm_iimn_fornorm %>% decostand(method = "rclr")%>%
  as.data.frame()%>%rownames_to_column("row.ID")
write_csv(norm_iimn_rclr, paste0(dat_path,"clustering/featuretable_reformated_justBA_summcollapse_cln_rclr.csv"))

## Read in FBMN library IDs from GNPS
library_ID <- read.delim(paste0(dat_path,"networking/library-results-merged_results_with_gnps_justBA.tsv"))%>%
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
  dplyr::select(c(1,12,2:11,13:46))%>%
  dplyr::rename(FeatureID=row.ID,x=correlation.group.ID, mn_mz=SpecMZ)%>%
  dplyr::rename(row.ID=x)

result_list<-rbind(result_list_iin,result_list)
write.table(result_list,paste0(dat_path,"networking/library-results-merged_results_with_gnps_justBA_collapse.tsv"),sep = "\t",row.names = FALSE, quote=FALSE)  
