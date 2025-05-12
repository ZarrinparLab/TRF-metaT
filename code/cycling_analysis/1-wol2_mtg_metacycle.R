setwd("~/Notebooks/sfloresr/TRF-metaT/")

library(MetaCycle)
library(dplyr)
library(data.table)
####################################################
#paths
dat_metadata<-"data/metaG_metadata_noNT.txt"
rna_nn<-"data/pfam_metaG/pfam_clean_noNT_TPM.txt"
rna_rpob<-"data/pfam_metaG/pfam_clean_noNT_TPM_normRPOB.txt"
dat_path<-"data/cycling_analysis/mgx_cyclic_analysis/"
####################################################
#not norm to rpob

sampledata<- fread(dat_metadata)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  arrange(zt_time)

dna <-fread(rna_nn)%>%
  dplyr::rename(FeatureID=`#OTU ID`)

filter_dna_FA<-dna%>%dplyr::select(1,8:13,2:7)
write.table(filter_dna_FA,paste0(dat_path,"wol2_dna_pipeline/filtered_dna_FA.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

meta2d(infile=paste0(dat_path,"wol2_dna_pipeline/filtered_dna_FA.txt"), cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir=paste0(dat_path,"wol2_dna_pipeline/FA_metacycle"),timepoints=rep(seq(1, 21, by=4),each=2),minper=20,maxper=24)

filter_dna_FT<-dna%>%dplyr::select(1,20:25,14:19)
write.table(filter_dna_FT,paste0(dat_path,"wol2_dna_pipeline/filtered_dna_FT.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

meta2d(infile=paste0(dat_path,"wol2_dna_pipeline/filtered_dna_FT.txt"), cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir=paste0(dat_path,"wol2_dna_pipeline/FT_metacycle"),timepoints=rep(seq(1, 21, by=4),each=2),minper=20,maxper=24)

filter_dna_NA<-dna%>%dplyr::select(1,35:43,26:34)
write.table(filter_dna_NA,paste0(dat_path,"wol2_dna_pipeline/filtered_dna_NA.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

meta2d(infile=paste0(dat_path,"wol2_dna_pipeline/filtered_dna_NA.txt"), cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir=paste0(dat_path,"wol2_dna_pipeline/NA_metacycle"),timepoints=rep(seq(1, 21, by=4),each=3),minper=20,maxper=24)

#rpob is not cycling
####################################################
#norm to rpob

sampledata<- fread(dat_metadata)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  arrange(zt_time)

dna <-fread(rna_rpob)

filter_dna_FA<-dna%>%dplyr::select(1,8:13,2:7)
write.table(filter_dna_FA,paste0(dat_path,"wol2_dna_pipeline_rpob/filtered_dna_FA.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

meta2d(infile=paste0(dat_path,"wol2_dna_pipeline_rpob/filtered_dna_FA.txt"), cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir=paste0(dat_path,"wol2_dna_pipeline_rpob/FA_metacycle"),timepoints=rep(seq(1, 21, by=4),each=2),minper=20,maxper=24)

filter_dna_FT<-dna%>%dplyr::select(1,20:25,14:19)
write.table(filter_dna_FT,paste0(dat_path,"wol2_dna_pipeline_rpob/filtered_dna_FT.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

meta2d(infile=paste0(dat_path,"wol2_dna_pipeline_rpob/filtered_dna_FT.txt"), cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir=paste0(dat_path,"wol2_dna_pipeline_rpob/FT_metacycle"),timepoints=rep(seq(1, 21, by=4),each=2),minper=20,maxper=24)

filter_dna_NA<-dna%>%dplyr::select(1,35:43,26:34)
write.table(filter_dna_NA,paste0(dat_path,"wol2_dna_pipeline_rpob/filtered_dna_NA.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

meta2d(infile=paste0(dat_path,"wol2_dna_pipeline_rpob/filtered_dna_NA.txt"), cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir=paste0(dat_path,"wol2_dna_pipeline_rpob/NA_metacycle"),timepoints=rep(seq(1, 21, by=4),each=3),minper=20,maxper=24)
