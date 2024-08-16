#metagenomic cycling 

setwd("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metagenomic/woltka2_results/filtered_metaG/cyclic_analysis/")
library(MetaCycle)
library(tidyverse)
library(data.table)
####################################################
#not norm to rpob

sampledata<- fread("~/scratch/TRF_multiomics/metagenomic/metaG_metadata_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  arrange(zt_time)

dna <-fread("~/scratch/TRF_multiomics/metagenomic/woltka2_results/filtered_metaG/pfam_TPM/pfam_clean_noNT_TPM.txt")%>%
  dplyr::rename(FeatureID=`#OTU ID`)

filter_dna_FA<-dna%>%dplyr::select(1,8:13,2:7)
write.table(filter_dna_FA,"wol2_dna_pipeline/filtered_dna_FA.txt",sep = "\t",row.names = FALSE, quote=FALSE)

meta2d(infile="wol2_dna_pipeline/filtered_dna_FA.txt", cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir="wol2_dna_pipeline/FA_metacycle",timepoints=rep(seq(1, 21, by=4),each=2),minper=20,maxper=24)

filter_dna_FT<-dna%>%dplyr::select(1,20:25,14:19)
write.table(filter_dna_FT,"wol2_dna_pipeline/filtered_dna_FT.txt",sep = "\t",row.names = FALSE, quote=FALSE)

meta2d(infile="wol2_dna_pipeline/filtered_dna_FT.txt", cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir="wol2_dna_pipeline/FT_metacycle",timepoints=rep(seq(1, 21, by=4),each=2),minper=20,maxper=24)

filter_dna_NA<-dna%>%dplyr::select(1,35:43,26:34)
write.table(filter_dna_NA,"wol2_dna_pipeline/filtered_dna_NA.txt",sep = "\t",row.names = FALSE, quote=FALSE)

meta2d(infile="wol2_dna_pipeline/filtered_dna_NA.txt", cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir="wol2_dna_pipeline/NA_metacycle",timepoints=rep(seq(1, 21, by=4),each=3),minper=20,maxper=24)

#rpob is not cycling
####################################################

#norm to rpob
sampledata<- fread("~/scratch/TRF_multiomics/metagenomic/metaG_metadata_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  arrange(zt_time)

dna <-fread("~/scratch/TRF_multiomics/metagenomic/woltka2_results/filtered_metaG/pfam_TPM/pfam_clean_noNT_TPM_normRPOB.txt")

filter_dna_FA<-dna%>%dplyr::select(1,8:13,2:7)
write.table(filter_dna_FA,"wol2_dna_pipeline_rpob/filtered_dna_FA.txt",sep = "\t",row.names = FALSE, quote=FALSE)

meta2d(infile="wol2_dna_pipeline_rpob/filtered_dna_FA.txt", cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir="wol2_dna_pipeline_rpob/FA_metacycle",timepoints=rep(seq(1, 21, by=4),each=2),minper=20,maxper=24)

filter_dna_FT<-dna%>%dplyr::select(1,20:25,14:19)
write.table(filter_dna_FT,"wol2_dna_pipeline_rpob/filtered_dna_FT.txt",sep = "\t",row.names = FALSE, quote=FALSE)

meta2d(infile="wol2_dna_pipeline_rpob/filtered_dna_FT.txt", cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir="wol2_dna_pipeline_rpob/FT_metacycle",timepoints=rep(seq(1, 21, by=4),each=2),minper=20,maxper=24)

filter_dna_NA<-dna%>%dplyr::select(1,35:43,26:34)
write.table(filter_dna_NA,"wol2_dna_pipeline_rpob/filtered_dna_NA.txt",sep = "\t",row.names = FALSE, quote=FALSE)

meta2d(infile="wol2_dna_pipeline_rpob/filtered_dna_NA.txt", cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir="wol2_dna_pipeline_rpob/NA_metacycle",timepoints=rep(seq(1, 21, by=4),each=3),minper=20,maxper=24)


####################################################

#run it on the genus -not norm
sampledata<- fread("~/scratch/TRF_multiomics/metagenomic/metaG_metadata_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  arrange(zt_time)

dna <-fread("~/scratch/TRF_multiomics/metagenomic/woltka2_results/filtered_metaG/genus_clean_noNT.txt")

filter_dna_FA<-dna%>%dplyr::select(1,8:13,2:7)
write.table(filter_dna_FA,"wol2_dna_pipeline_genus/filtered_dna_FA.txt",sep = "\t",row.names = FALSE, quote=FALSE)

meta2d(infile="wol2_dna_pipeline_genus/filtered_dna_FA.txt", cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir="wol2_dna_pipeline_genus/FA_metacycle",timepoints=rep(seq(1, 21, by=4),each=2),minper=20,maxper=24)

filter_dna_FT<-dna%>%dplyr::select(1,20:25,14:19)
write.table(filter_dna_FT,"wol2_dna_pipeline_genus/filtered_dna_FT.txt",sep = "\t",row.names = FALSE, quote=FALSE)

meta2d(infile="wol2_dna_pipeline_genus/filtered_dna_FT.txt", cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir="wol2_dna_pipeline_genus/FT_metacycle",timepoints=rep(seq(1, 21, by=4),each=2),minper=20,maxper=24)

filter_dna_NA<-dna%>%dplyr::select(1,35:43,26:34)
write.table(filter_dna_NA,"wol2_dna_pipeline_genus/filtered_dna_NA.txt",sep = "\t",row.names = FALSE, quote=FALSE)

meta2d(infile="wol2_dna_pipeline_genus/filtered_dna_NA.txt", cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir="wol2_dna_pipeline_genus/NA_metacycle",timepoints=rep(seq(1, 21, by=4),each=3),minper=20,maxper=24)

#run it on the genus -rclr norm
sampledata<- fread("~/scratch/TRF_multiomics/metagenomic/metaG_metadata_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  arrange(zt_time)

dna <-fread("~/scratch/TRF_multiomics/metagenomic/woltka2_results/filtered_metaG/genus_clean_noNT_rclr.txt")

filter_dna_FA<-dna%>%dplyr::select(1,8:13,2:7)
write.table(filter_dna_FA,"wol2_dna_pipeline_genus-rclr/filtered_dna_FA.txt",sep = "\t",row.names = FALSE, quote=FALSE)

meta2d(infile="wol2_dna_pipeline_genus-rclr/filtered_dna_FA.txt", cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir="wol2_dna_pipeline_genus-rclr/FA_metacycle",timepoints=rep(seq(1, 21, by=4),each=2),minper=20,maxper=24)

filter_dna_FT<-dna%>%dplyr::select(1,20:25,14:19)
write.table(filter_dna_FT,"wol2_dna_pipeline_genus-rclr/filtered_dna_FT.txt",sep = "\t",row.names = FALSE, quote=FALSE)

meta2d(infile="wol2_dna_pipeline_genus-rclr/filtered_dna_FT.txt", cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir="wol2_dna_pipeline_genus-rclr/FT_metacycle",timepoints=rep(seq(1, 21, by=4),each=2),minper=20,maxper=24)

filter_dna_NA<-dna%>%dplyr::select(1,35:43,26:34)
write.table(filter_dna_NA,"wol2_dna_pipeline_genus-rclr/filtered_dna_NA.txt",sep = "\t",row.names = FALSE, quote=FALSE)

meta2d(infile="wol2_dna_pipeline_genus-rclr/filtered_dna_NA.txt", cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir="wol2_dna_pipeline_genus-rclr/NA_metacycle",timepoints=rep(seq(1, 21, by=4),each=3),minper=20,maxper=24)

#run it on the genus -prop norm
sampledata<- fread("~/scratch/TRF_multiomics/metagenomic/metaG_metadata_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  arrange(zt_time)

dna <-fread("~/scratch/TRF_multiomics/metagenomic/woltka2_results/filtered_metaG/genus_clean_noNT_prop.txt")

filter_dna_FA<-dna%>%dplyr::select(1,8:13,2:7)
write.table(filter_dna_FA,"wol2_dna_pipeline_genus-prop/filtered_dna_FA.txt",sep = "\t",row.names = FALSE, quote=FALSE)

meta2d(infile="wol2_dna_pipeline_genus-prop/filtered_dna_FA.txt", cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir="wol2_dna_pipeline_genus-prop/FA_metacycle",timepoints=rep(seq(1, 21, by=4),each=2),minper=20,maxper=24)

filter_dna_FT<-dna%>%dplyr::select(1,20:25,14:19)
write.table(filter_dna_FT,"wol2_dna_pipeline_genus-prop/filtered_dna_FT.txt",sep = "\t",row.names = FALSE, quote=FALSE)

meta2d(infile="wol2_dna_pipeline_genus-prop/filtered_dna_FT.txt", cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir="wol2_dna_pipeline_genus-prop/FT_metacycle",timepoints=rep(seq(1, 21, by=4),each=2),minper=20,maxper=24)

filter_dna_NA<-dna%>%dplyr::select(1,35:43,26:34)
write.table(filter_dna_NA,"wol2_dna_pipeline_genus-prop/filtered_dna_NA.txt",sep = "\t",row.names = FALSE, quote=FALSE)

meta2d(infile="wol2_dna_pipeline_genus-prop/filtered_dna_NA.txt", cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir="wol2_dna_pipeline_genus-prop/NA_metacycle",timepoints=rep(seq(1, 21, by=4),each=3),minper=20,maxper=24)
####################################################

#run it on the species-not norm
sampledata<- fread("~/scratch/TRF_multiomics/metagenomic/metaG_metadata_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  arrange(zt_time)

dna <-fread("~/scratch/TRF_multiomics/metagenomic/woltka2_results/filtered_metaG/species_clean_noNT.txt")

filter_dna_FA<-dna%>%dplyr::select(1,8:13,2:7)
write.table(filter_dna_FA,"wol2_dna_pipeline_species/filtered_dna_FA.txt",sep = "\t",row.names = FALSE, quote=FALSE)

meta2d(infile="wol2_dna_pipeline_species/filtered_dna_FA.txt", cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir="wol2_dna_pipeline_species/FA_metacycle",timepoints=rep(seq(1, 21, by=4),each=2),minper=20,maxper=24)

filter_dna_FT<-dna%>%dplyr::select(1,20:25,14:19)
write.table(filter_dna_FT,"wol2_dna_pipeline_species/filtered_dna_FT.txt",sep = "\t",row.names = FALSE, quote=FALSE)

meta2d(infile="wol2_dna_pipeline_species/filtered_dna_FT.txt", cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir="wol2_dna_pipeline_species/FT_metacycle",timepoints=rep(seq(1, 21, by=4),each=2),minper=20,maxper=24)

filter_dna_NA<-dna%>%dplyr::select(1,35:43,26:34)
write.table(filter_dna_NA,"wol2_dna_pipeline_species/filtered_dna_NA.txt",sep = "\t",row.names = FALSE, quote=FALSE)

meta2d(infile="wol2_dna_pipeline_species/filtered_dna_NA.txt", cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir="wol2_dna_pipeline_species/NA_metacycle",timepoints=rep(seq(1, 21, by=4),each=3),minper=20,maxper=24)

#run it on the species-norm to rclr
sampledata<- fread("~/scratch/TRF_multiomics/metagenomic/metaG_metadata_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  arrange(zt_time)

dna <-fread("~/scratch/TRF_multiomics/metagenomic/woltka2_results/filtered_metaG/species_clean_noNT_rclr.txt")

filter_dna_FA<-dna%>%dplyr::select(1,8:13,2:7)
write.table(filter_dna_FA,"wol2_dna_pipeline_species-rclr/filtered_dna_FA.txt",sep = "\t",row.names = FALSE, quote=FALSE)

meta2d(infile="wol2_dna_pipeline_species-rclr/filtered_dna_FA.txt", cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir="wol2_dna_pipeline_species-rclr/FA_metacycle",timepoints=rep(seq(1, 21, by=4),each=2),minper=20,maxper=24)

filter_dna_FT<-dna%>%dplyr::select(1,20:25,14:19)
write.table(filter_dna_FT,"wol2_dna_pipeline_species-rclr/filtered_dna_FT.txt",sep = "\t",row.names = FALSE, quote=FALSE)

meta2d(infile="wol2_dna_pipeline_species-rclr/filtered_dna_FT.txt", cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir="wol2_dna_pipeline_species-rclr/FT_metacycle",timepoints=rep(seq(1, 21, by=4),each=2),minper=20,maxper=24)

filter_dna_NA<-dna%>%dplyr::select(1,35:43,26:34)
write.table(filter_dna_NA,"wol2_dna_pipeline_species-rclr/filtered_dna_NA.txt",sep = "\t",row.names = FALSE, quote=FALSE)

meta2d(infile="wol2_dna_pipeline_species-rclr/filtered_dna_NA.txt", cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir="wol2_dna_pipeline_species-rclr/NA_metacycle",timepoints=rep(seq(1, 21, by=4),each=3),minper=20,maxper=24)

#run it on the species-norm to prop
sampledata<- fread("~/scratch/TRF_multiomics/metagenomic/metaG_metadata_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  arrange(zt_time)

dna <-fread("~/scratch/TRF_multiomics/metagenomic/woltka2_results/filtered_metaG/species_clean_noNT_prop.txt")

filter_dna_FA<-dna%>%dplyr::select(1,8:13,2:7)
write.table(filter_dna_FA,"wol2_dna_pipeline_species-prop/filtered_dna_FA.txt",sep = "\t",row.names = FALSE, quote=FALSE)

meta2d(infile="wol2_dna_pipeline_species-prop/filtered_dna_FA.txt", cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir="wol2_dna_pipeline_species-prop/FA_metacycle",timepoints=rep(seq(1, 21, by=4),each=2),minper=20,maxper=24)

filter_dna_FT<-dna%>%dplyr::select(1,20:25,14:19)
write.table(filter_dna_FT,"wol2_dna_pipeline_species-prop/filtered_dna_FT.txt",sep = "\t",row.names = FALSE, quote=FALSE)

meta2d(infile="wol2_dna_pipeline_species-prop/filtered_dna_FT.txt", cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir="wol2_dna_pipeline_species-prop/FT_metacycle",timepoints=rep(seq(1, 21, by=4),each=2),minper=20,maxper=24)

filter_dna_NA<-dna%>%dplyr::select(1,35:43,26:34)
write.table(filter_dna_NA,"wol2_dna_pipeline_species-prop/filtered_dna_NA.txt",sep = "\t",row.names = FALSE, quote=FALSE)

meta2d(infile="wol2_dna_pipeline_species-prop/filtered_dna_NA.txt", cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir="wol2_dna_pipeline_species-prop/NA_metacycle",timepoints=rep(seq(1, 21, by=4),each=3),minper=20,maxper=24)

####################################################

#run it at the OGU level

#rclr

sampledata<- fread("~/scratch/TRF_multiomics/metagenomic/metaG_metadata_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  arrange(zt_time)

dna <-fread("~/scratch/TRF_multiomics/metagenomic/woltka2_results/filtered_metaG/none_clean_noNT_rclr.txt")

filter_dna_FA<-dna%>%dplyr::select(1,8:13,2:7)
write.table(filter_dna_FA,"wol2_dna_pipeline_ogu-rclr/filtered_dna_FA.txt",sep = "\t",row.names = FALSE, quote=FALSE)

meta2d(infile="wol2_dna_pipeline_ogu-rclr/filtered_dna_FA.txt", cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir="wol2_dna_pipeline_ogu-rclr/FA_metacycle",timepoints=rep(seq(1, 21, by=4),each=2),minper=20,maxper=24)

filter_dna_FT<-dna%>%dplyr::select(1,20:25,14:19)
write.table(filter_dna_FT,"wol2_dna_pipeline_ogu-rclr/filtered_dna_FT.txt",sep = "\t",row.names = FALSE, quote=FALSE)

meta2d(infile="wol2_dna_pipeline_ogu-rclr/filtered_dna_FT.txt", cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir="wol2_dna_pipeline_ogu-rclr/FT_metacycle",timepoints=rep(seq(1, 21, by=4),each=2),minper=20,maxper=24)

filter_dna_NA<-dna%>%dplyr::select(1,35:43,26:34)
write.table(filter_dna_NA,"wol2_dna_pipeline_ogu-rclr/filtered_dna_NA.txt",sep = "\t",row.names = FALSE, quote=FALSE)

meta2d(infile="wol2_dna_pipeline_ogu-rclr/filtered_dna_NA.txt", cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir="wol2_dna_pipeline_ogu-rclr/NA_metacycle",timepoints=rep(seq(1, 21, by=4),each=3),minper=20,maxper=24)

#prop

sampledata<- fread("~/scratch/TRF_multiomics/metagenomic/metaG_metadata_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  arrange(zt_time)

dna <-fread("~/scratch/TRF_multiomics/metagenomic/woltka2_results/filtered_metaG/none_clean_noNT_prop.txt")

filter_dna_FA<-dna%>%dplyr::select(1,8:13,2:7)
write.table(filter_dna_FA,"wol2_dna_pipeline_ogu-prop/filtered_dna_FA.txt",sep = "\t",row.names = FALSE, quote=FALSE)

meta2d(infile="wol2_dna_pipeline_ogu-prop/filtered_dna_FA.txt", cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir="wol2_dna_pipeline_ogu-prop/FA_metacycle",timepoints=rep(seq(1, 21, by=4),each=2),minper=20,maxper=24)

filter_dna_FT<-dna%>%dplyr::select(1,20:25,14:19)
write.table(filter_dna_FT,"wol2_dna_pipeline_ogu-prop/filtered_dna_FT.txt",sep = "\t",row.names = FALSE, quote=FALSE)

meta2d(infile="wol2_dna_pipeline_ogu-prop/filtered_dna_FT.txt", cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir="wol2_dna_pipeline_ogu-prop/FT_metacycle",timepoints=rep(seq(1, 21, by=4),each=2),minper=20,maxper=24)

filter_dna_NA<-dna%>%dplyr::select(1,35:43,26:34)
write.table(filter_dna_NA,"wol2_dna_pipeline_ogu-prop/filtered_dna_NA.txt",sep = "\t",row.names = FALSE, quote=FALSE)

meta2d(infile="wol2_dna_pipeline_ogu-prop/filtered_dna_NA.txt", cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir="wol2_dna_pipeline_ogu-prop/NA_metacycle",timepoints=rep(seq(1, 21, by=4),each=3),minper=20,maxper=24)

