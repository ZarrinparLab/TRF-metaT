#metatranscriptomics cycling 

setwd("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_m_results/pfam/cyclic_analysis/")
library(MetaCycle)
library(tidyverse)
library(data.table)
####################################################
#not norm to rpob

sampledata<- fread("~/scratch/TRF_multiomics/metatranscript/woltka2_m_results/metaT_metadata_ztcat_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  arrange(zt_time)

rna <-fread("~/scratch/TRF_multiomics/metatranscript/woltka2_m_results/pfam/pfam-TPM_clean_noNT.tsv")

filter_rna_FA<-rna%>%dplyr::select(1,2:13)
write.table(filter_rna_FA,"wol2_rna_pipeline/filtered_rna_FA.txt",sep = "\t",row.names = FALSE, quote=FALSE)

meta2d(infile="wol2_rna_pipeline/filtered_rna_FA.txt", cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir="wol2_rna_pipeline/FA_metacycle",timepoints=rep(seq(1, 21, by=4),each=2),minper=20,maxper=24)

filter_rna_FT<-rna%>%dplyr::select(1,14:25)
write.table(filter_rna_FT,"wol2_rna_pipeline/filtered_rna_FT.txt",sep = "\t",row.names = FALSE, quote=FALSE)

meta2d(infile="wol2_rna_pipeline/filtered_rna_FT.txt", cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir="wol2_rna_pipeline/FT_metacycle",timepoints=rep(seq(1, 21, by=4),each=2),minper=20,maxper=24)

filter_rna_NA<-rna%>%dplyr::select(1,26:43)
write.table(filter_rna_NA,"wol2_rna_pipeline/filtered_rna_NA.txt",sep = "\t",row.names = FALSE, quote=FALSE)

meta2d(infile="wol2_rna_pipeline/filtered_rna_NA.txt", cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir="wol2_rna_pipeline/NA_metacycle",timepoints=rep(seq(1, 21, by=4),each=3),minper=20,maxper=24)

#rpob is not cycling
####################################################

#norm to rpob
sampledata<- fread("~/scratch/TRF_multiomics/metatranscript/woltka2_m_results/metaT_metadata_ztcat_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  arrange(zt_time)

rna <-fread("~/scratch/TRF_multiomics/metatranscript/woltka2_m_results/pfam/pfam-TPM_clean_noNT_normRPOB.txt")

filter_rna_FA<-rna%>%dplyr::select(1,2:13)
write.table(filter_rna_FA,"wol2_rna_pipeline_rpob/filtered_rna_FA.txt",sep = "\t",row.names = FALSE, quote=FALSE)

meta2d(infile="wol2_rna_pipeline_rpob/filtered_rna_FA.txt", cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir="wol2_rna_pipeline_rpob/FA_metacycle",timepoints=rep(seq(1, 21, by=4),each=2),minper=20,maxper=24)

filter_rna_FT<-rna%>%dplyr::select(1,14:25)
write.table(filter_rna_FT,"wol2_rna_pipeline_rpob/filtered_rna_FT.txt",sep = "\t",row.names = FALSE, quote=FALSE)

meta2d(infile="wol2_rna_pipeline_rpob/filtered_rna_FT.txt", cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir="wol2_rna_pipeline_rpob/FT_metacycle",timepoints=rep(seq(1, 21, by=4),each=2),minper=20,maxper=24)

filter_rna_NA<-rna%>%dplyr::select(1,26:43)
write.table(filter_rna_NA,"wol2_rna_pipeline_rpob/filtered_rna_NA.txt",sep = "\t",row.names = FALSE, quote=FALSE)

meta2d(infile="wol2_rna_pipeline_rpob/filtered_rna_NA.txt", cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir="wol2_rna_pipeline_rpob/NA_metacycle",timepoints=rep(seq(1, 21, by=4),each=3),minper=20,maxper=24)


####################################################
#rarefied norm to rpob
sampledata<- fread("~/scratch/TRF_multiomics/metatranscript/woltka2_m_results/metaT_metadata_ztcat_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  arrange(zt_time)

rna <-fread("~/scratch/TRF_multiomics/metatranscript/woltka2_m_results/g-diversity-core-metrics114K/rarefied_table/pfam-144k_clean_noNT_normRPOB.txt")

filter_rna_FA<-rna%>%dplyr::select(1,2:13)
write.table(filter_rna_FA,"wol2_rna_pipeline_144k_rpob/filtered_rna_FA.txt",sep = "\t",row.names = FALSE, quote=FALSE)

meta2d(infile="wol2_rna_pipeline_144k_rpob/filtered_rna_FA.txt", cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir="wol2_rna_pipeline_144k_rpob/FA_metacycle",timepoints=rep(seq(1, 21, by=4),each=2),minper=20,maxper=24)

filter_rna_FT<-rna%>%dplyr::select(1,14:25)
write.table(filter_rna_FT,"wol2_rna_pipeline_144k_rpob/filtered_rna_FT.txt",sep = "\t",row.names = FALSE, quote=FALSE)

meta2d(infile="wol2_rna_pipeline_144k_rpob/filtered_rna_FT.txt", cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir="wol2_rna_pipeline_144k_rpob/FT_metacycle",timepoints=rep(seq(1, 21, by=4),each=2),minper=20,maxper=24)

filter_rna_NA<-rna%>%dplyr::select(1,26:43)
write.table(filter_rna_NA,"wol2_rna_pipeline_144k_rpob/filtered_rna_NA.txt",sep = "\t",row.names = FALSE, quote=FALSE)

meta2d(infile="wol2_rna_pipeline_144k_rpob/filtered_rna_NA.txt", cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir="wol2_rna_pipeline_144k_rpob/NA_metacycle",timepoints=rep(seq(1, 21, by=4),each=3),minper=20,maxper=24)

