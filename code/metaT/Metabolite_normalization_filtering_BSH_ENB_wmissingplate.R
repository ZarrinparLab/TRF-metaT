setwd("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_m_results/")
library(tidyverse)
library(data.table)
library(stringr)
library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)
library(fuzzyjoin)

#####################################################################################################
#combine the additional plate that was run to table

##original run
annot<-fread("BSH/BSH_ENB/GNPS/nf_output/networking/library-results-merged_results_with_gnps.tsv")%>%
  filter(Organism=="BILELIB19" | Organism=="GNPS-BILE-ACID-MODIFICATIONS" )
write.table(annot,"BSH/BSH_ENB/GNPS/nf_output/networking/library-results-merged_results_with_gnps_justBA.tsv",sep = "\t",row.names = FALSE, quote=FALSE)  

mtb<-fread("BSH/BSH_ENB/GNPS/nf_output/clustering/featuretable_reformated.csv")%>%
  filter(`row ID` %in% annot$`#Scan#`)
write_csv(mtb, "BSH/BSH_ENB/GNPS/nf_output/clustering/featuretable_reformated_justBA.csv")

md<-fread("BSH/BSH_ENB/GNPS/metadata_filename/ZT_BSH_metadata.tsv")%>%
  filter(!(plate=="P1"|plate=="P2"))

mtb_rmp12<-mtb%>%
  dplyr::select(c(1:14),any_of(md$sample_name))
write_csv(mtb, "BSH/BSH_ENB/GNPS/nf_output/clustering/featuretable_reformated_justBA_rmp12.csv")

##missing plates
annot2<-fread("BSH/BSH_ENB/GNPS_missingplates/nf_output/networking/library-results-merged_results_with_gnps.tsv")%>%
filter(Organism=="BILELIB19" | Organism=="GNPS-BILE-ACID-MODIFICATIONS" )
write.table(annot,"BSH/BSH_ENB/GNPS_missingplates/nf_output/networking/library-results-merged_results_with_gnps_justBA.tsv",sep = "\t",row.names = FALSE, quote=FALSE)  

mtb2<-fread("BSH/BSH_ENB/GNPS_missingplates/nf_output/clustering/featuretable_reformated.csv")%>%
  filter(`row ID` %in% annot2$`#Scan#`)
write_csv(mtb, "BSH/BSH_ENB/GNPS_missingplates/nf_output/clustering/featuretable_reformated_justBA.csv")

#combined based on annot
annot_sub<-annot%>%
  dplyr::select(`#Scan#`,Compound_Name)%>%
  dplyr::rename(`row ID`=`#Scan#`)

mtb_rmp12_sub<-mtb_rmp12%>%
  dplyr::select(`row ID`,`row m/z`,`row retention time`)%>%
  left_join(.,annot_sub,by="row ID")

annot2_sub<-annot2%>%
  dplyr::select(`#Scan#`,Compound_Name)%>%
  dplyr::rename(`row ID`=`#Scan#`)

mtb2_sub<-mtb2%>%
  dplyr::select(`row ID`,`row m/z`,`row retention time`)%>%
  left_join(.,annot2_sub,by="row ID")

initial_join <- left_join(mtb_rmp12_sub, mtb2_sub, by = "Compound_Name")
combannot <- initial_join %>%
  group_by(Compound_Name) %>%
  #nest()
  filter(abs(`row retention time.x` - `row retention time.y`) <= 0.05)
  
  #subset annotations only in plate 1/2 for original run 

## Cleanup feature quant table from MZmine ##########################################################
norm <- read.csv("BSH/BSH_ENB/GNPS/nf_output/clustering/featuretable_reformated_justBA_rmp12.csv")%>%
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

write_csv(norm_iimn, "BSH/BSH_ENB/GNPS/nf_output/clustering/featuretable_reformated_justBA_summcollapse_rmp12.csv")
write_csv(norm_iimn, "BSH/BSH_ENB/GNPS/nf_output/clustering/featuretable_reformated_justBA_summcollapse_rmp12_cln.csv")

norm_iimn_fornorm<-norm_iimn%>%
  column_to_rownames("row.ID")%>%as.matrix()

norm_iimn_rclr <- norm_iimn_fornorm %>% decostand(method = "rclr")%>%
  as.data.frame()%>%rownames_to_column("row.ID")
write_csv(norm_iimn_rclr, "BSH/BSH_ENB/GNPS/nf_output/clustering/featuretable_reformated_justBA_summcollapse_rmp12_cln_rclr.csv")

# norm_iimn_alr <- norm_iimn_fornorm %>% decostand(method = "alr")%>%
#   as.data.frame()%>%rownames_to_column("row.ID")
# 
# norm_iimn_clr <- norm_iimn_fornorm %>% decostand(method = "clr",pseudocount = 1)%>%
#   as.data.frame()%>%rownames_to_column("row.ID")
# 
# norm_iimn_log <- norm_iimn_fornorm %>% decostand(method = "log")%>%
#   as.data.frame()%>%rownames_to_column("row.ID")
# write_csv(norm_iimn_log, "data_files/FBMN_simone/gnps_quant_collapseiin_sum_corrgroup_log.csv")

## Read in FBMN library IDs from GNPS
library_ID <- read.delim("BSH/BSH_ENB/GNPS/nf_output/networking/library-results-merged_results_with_gnps_justBA.tsv")%>%
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
write.table(result_list,"BSH/BSH_ENB/GNPS/nf_output/networking/library-results-merged_results_with_gnps_justBA_collapse_rmp12.tsv",sep = "\t",row.names = FALSE, quote=FALSE)  
