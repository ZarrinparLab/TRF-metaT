setwd("~/Notebooks/sfloresr/TRF-metaT/")

library(MetaCycle)
library(dplyr)
library(data.table)
####################################################
#paths
dat_metadata<-"data/cycling_analysis/metaT_metadata_ztcat_noNT.txt"
rna_nn<-"data/pfam_metaT/pfam-TPM_clean_noNT.tsv"
rna_rpob<-"data/pfam_metaT/pfam-TPM_clean_noNT_normRPOB.txt"
rna_144k<-"data/diversity_analysis/g-diversity-core-metrics114K/rarefied_table/pfam-144k_clean_noNT_normRPOB.txt"
dat_path<-"data/cycling_analysis/mtx_cyclic_analysis/"
####################################################
#not norm to rpob

sampledata<- fread(dat_metadata)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  arrange(zt_time)

rna <-fread(rna_nn)

filter_rna_FA<-rna%>%dplyr::select(1,2:13)
write.table(filter_rna_FA,paste0(dat_path,"wol2_rna_pipeline/filtered_rna_FA.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

meta2d(infile=paste0(dat_path,"wol2_rna_pipeline/filtered_rna_FA.txt"), cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir=paste0(dat_path,"wol2_rna_pipeline/FA_metacycle"),timepoints=rep(seq(1, 21, by=4),each=2),minper=20,maxper=24)

filter_rna_FT<-rna%>%dplyr::select(1,14:25)
write.table(filter_rna_FT,paste0(dat_path,"wol2_rna_pipeline/filtered_rna_FT.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

meta2d(infile=paste0(dat_path,"wol2_rna_pipeline/filtered_rna_FT.txt"), cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir=paste0(dat_path,"wol2_rna_pipeline/FT_metacycle"),timepoints=rep(seq(1, 21, by=4),each=2),minper=20,maxper=24)

filter_rna_NA<-rna%>%dplyr::select(1,26:43)
write.table(filter_rna_NA,paste0(dat_path,"wol2_rna_pipeline/filtered_rna_NA.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

meta2d(infile=paste0(dat_path,"wol2_rna_pipeline/filtered_rna_NA.txt"), cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir=paste0(dat_path,"wol2_rna_pipeline/NA_metacycle"),timepoints=rep(seq(1, 21, by=4),each=3),minper=20,maxper=24)

#rpob is not cycling

####################################################
#norm to rpob

sampledata<- fread(dat_metadata)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  arrange(zt_time)

rna <-fread(rna_rpob)

filter_rna_FA<-rna%>%dplyr::select(1,2:13)
write.table(filter_rna_FA,paste0(dat_path,"wol2_rna_pipeline_rpob/filtered_rna_FA.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

meta2d(infile=paste0(dat_path,"wol2_rna_pipeline_rpob/filtered_rna_FA.txt"), cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir=paste0(dat_path,"wol2_rna_pipeline_rpob/FA_metacycle"),timepoints=rep(seq(1, 21, by=4),each=2),minper=20,maxper=24)

filter_rna_FT<-rna%>%dplyr::select(1,14:25)
write.table(filter_rna_FT,paste0(dat_path,"wol2_rna_pipeline_rpob/filtered_rna_FT.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

meta2d(infile=paste0(dat_path,"wol2_rna_pipeline_rpob/filtered_rna_FT.txt"), cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir=paste0(dat_path,"wol2_rna_pipeline_rpob/FT_metacycle"),timepoints=rep(seq(1, 21, by=4),each=2),minper=20,maxper=24)

filter_rna_NA<-rna%>%dplyr::select(1,26:43)
write.table(filter_rna_NA,paste0(dat_path,"wol2_rna_pipeline_rpob/filtered_rna_NA.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

meta2d(infile=paste0(dat_path,"wol2_rna_pipeline_rpob/filtered_rna_NA.txt"), cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir=paste0(dat_path,"wol2_rna_pipeline_rpob/NA_metacycle"),timepoints=rep(seq(1, 21, by=4),each=3),minper=20,maxper=24)


####################################################
#rarefied norm to rpob
sampledata<- fread(dat_metadata)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  arrange(zt_time)

rna <-fread(rna_144k)

filter_rna_FA<-rna%>%dplyr::select(1,2:13)
write.table(filter_rna_FA,paste0(dat_path,"wol2_rna_pipeline_144k_rpob/filtered_rna_FA.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

meta2d(infile=paste0(dat_path,"wol2_rna_pipeline_144k_rpob/filtered_rna_FA.txt"), cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir=paste0(dat_path,"wol2_rna_pipeline_144k_rpob/FA_metacycle"),timepoints=rep(seq(1, 21, by=4),each=2),minper=20,maxper=24)

filter_rna_FT<-rna%>%dplyr::select(1,14:25)
write.table(filter_rna_FT,paste0(dat_path,"wol2_rna_pipeline_144k_rpob/filtered_rna_FT.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

meta2d(infile=paste0(dat_path,"wol2_rna_pipeline_144k_rpob/filtered_rna_FT.txt"), cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir=paste0(dat_path,"wol2_rna_pipeline_144k_rpob/FT_metacycle"),timepoints=rep(seq(1, 21, by=4),each=2),minper=20,maxper=24)

filter_rna_NA<-rna%>%dplyr::select(1,26:43)
write.table(filter_rna_NA,paste0(dat_path,"wol2_rna_pipeline_144k_rpob/filtered_rna_NA.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

meta2d(infile=paste0(dat_path,"wol2_rna_pipeline_144k_rpob/filtered_rna_NA.txt"), cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir=paste0(dat_path,"wol2_rna_pipeline_144k_rpob/NA_metacycle"),timepoints=rep(seq(1, 21, by=4),each=3),minper=20,maxper=24)

