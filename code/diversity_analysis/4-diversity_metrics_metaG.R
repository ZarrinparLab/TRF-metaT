setwd("~/Notebooks/sfloresr/TRF-metaT/")

library(tidyverse)
library(data.table)
library("qiime2R")
library(ggpubr)

###########################################################
#paths
dat_metadata_ori<-"data/metaG_metadata.txt"
dat_metadata<-"data/metaG_metadata_noNT.txt"
notnorm_data<-"data/pfam_metaG/pfam.tsv"
dat_ord<-"data/diversity_analysis/rpca_results_metaG/ordination.qza"
dat_path<-"data/pfam_metaG/"
res_path<-"data/diversity_analysis/rpca_results_metaG/"
fig_path<-"figures/diversity_analysis/"
###########################################################
#functions

subset_dat<-function(dt){
  samps_list<-c("#FeatureID","FA1b_S36","FA1c_S37","FA2b_S38","FA2c_S39","FA3b_S40","FA3c_S41","FA4b_S42","FA4c_S43",  
                "FA5b_S44","FA5c_S45","FA6b_S46","FA6c_S47","FT1b_S48","FT1c_S49","FT2b_S50","FT2c_S51","FT3b_S52",
                "FT3c_S53","FT4b_S54","FT4c_S55","FT5b_S56","FT5c_S57","FT6b_S58","FT6c_S59","NA1a_S1","NA1b_S2",  
                "NA1c_S3","NA2a_S4","NA2b_S5","NA2c_S6","NA3a_S7","NA3b_S8","NA3c_S9","NA4a_S10","NA4b_S11",  
                "NA4c_S12","NA5a_S13","NA5b_S14","NA5c_S15","NA6a_S16","NA6b_S17","NA6c_S18")
  dat<-fread(dt)%>%dplyr::select(all_of(samps_list))
  names(dat) <- c("FeatureID","FA1b","FA1c","FA2b","FA2c","FA3b","FA3c","FA4b","FA4c","FA5b","FA5c",
                  "FA6b","FA6c","FT1b","FT1c","FT2b","FT2c","FT3b","FT3c","FT4b","FT4c","FT5b",
                  "FT5c","FT6b","FT6c","NA1a","NA1b","NA1c","NA2a","NA2b","NA2c","NA3a","NA3b","NA3c",
                  "NA4a","NA4b","NA4c","NA5a","NA5b","NA5c","NA6a","NA6b","NA6c")
  return(dat)
}

rpcacond_plt<-function(dat,cond){
  
  rpca_sub<-dat%>%
    filter(condition==cond)%>%
    mutate(zt_time=factor(zt_time,levels=c("1","5","9","13","17","21")))
  
  p<-rpca_sub %>%
    ggplot(aes(x=PC1, y=PC2, fill=zt_time)) +
    geom_point(alpha=1.0,size=3,shape=21) + 
    theme_pubr() +
    scale_fill_manual(values=c("#67001f","#d6604d","#ffba92","#8ccff3","#4393c3","#053061"))+
    labs(color="ZT time",
         x =paste("PC1 (",round(ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
         y =paste("PC2 (",round(ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""))+ggtitle(paste("MGX",cond,sep=" "))+ theme(plot.title = element_text(face = "bold"))
  
  return(p)
}

###########################################################

#clean metadata
md<-fread(dat_metadata_ori)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  filter(condition!="NT")%>%
  mutate(cond_phase=paste(condition,lightdark,sep="_"),
         cond_zt=paste(condition,zt_time,sep="_"),
         fasted=ifelse(zt_time %in% c(5,9,13), "fasted","fed"),
         cond_fasted=paste(condition,fasted,sep="_"))
write.table(md,dat_metadata,sep = "\t",row.names = FALSE, quote=FALSE)

#clean annotations
annotations<-fread(notnorm_data) %>%dplyr::select(c(1,61))%>%
  dplyr::rename(FeatureID="#FeatureID")
write.csv(annotations,paste0(dat_path,"pfam_annotationkey.csv"), row.names = FALSE)

#subset data
dat_nn<-subset_dat(notnorm_data)
write.table(dat_nn,paste0(dat_path,"pfam_clean_noNT.tsv"),sep = "\t",row.names = FALSE, quote=FALSE)

#####################################################################################
#used emperor to plot RPCA--Figure 1C
#permamova values are in condition-significance.qzv files 
#need to run rna_rpca_metaG.sh script first

###########################################################
#plot the ZT time for each condition--Figure 1E-G (middle)
#need to run rna_rpca_metaG.sh script first

ord <- read_qza(dat_ord)

samp_ord<-ord$data$Vectors
write.table(samp_ord,paste0(res_path,"sample_ordination.txt"),sep = "\t",row.names = FALSE,quote=FALSE)
feat_ord<-ord$data$Species
write.table(feat_ord,paste0(res_path,"feature_ordination.txt"),sep = "\t",row.names = FALSE,quote=FALSE)

md<-fread(dat_metadata)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))

rpca<-ord$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2,PC3)%>%
  dplyr::rename(sample_name=SampleID)%>%
  left_join(md,by="sample_name")%>%
  mutate(condition=factor(condition,levels=c("FT","FA","NA")),
         phase=factor(lightdark,levels=c("light","dark")))

p<-rpcacond_plt(rpca,"NA")
ggsave(paste0(fig_path,"SFR24_0412_mgx_NA_RPCA.pdf"), plot=p,height=3.5, width=3.5)

p<-rpcacond_plt(rpca,"FA")
ggsave(paste0(fig_path,"SFR24_0412_mgx_FA_RPCA.pdf"), plot=p,height=3.5, width=3.5)

p<-rpcacond_plt(rpca,"FT")
ggsave(paste0(fig_path,"SFR24_0412_mgx_FT_RPCA.pdf"), plot=p,height=3.5, width=3.5)

#permamova values are in cond_fasted-significance.qzv and condphase-significance.qzv files 