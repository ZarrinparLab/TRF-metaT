setwd("~/Notebooks/sfloresr/TRF-metaT/")

library(tidyverse)
library(data.table)
library("qiime2R")
library(ggpubr)

###########################################################
#paths
dat_metadata_ori<-"data/trf_cecal_metadata_noNT.tsv"
dat_metadata<-"data/metaT_metadata_ztcat_noNT.txt"
notnorm_data<-"data/pfam_metaT/pfam.tsv"
dat_ord<-"data/diversity_analysis/rpca_results_metaT/ordination.qza"
dat_path<-"data/pfam_metaT/"
res_path<-"data/diversity_analysis/rpca_results_metaT/"
fig_path<-"figures/diversity_analysis/"
###########################################################
#functions
subset_dat<-function(dt){
  samps_list<-c("#FeatureID","cFA01a_S39","cFA01b_S45","cFA05a_S40","cFA05b_S46","cFA09a_S41",
                "cFA09b_S47","cFA13a_S36","cFA13b_S42","cFA17a_S37","cFA17b_S43",
                "cFA21a_S38","cFA21b_S44","cFT01a_S51","cFT01b_S57","cFT05a_S52",
                "cFT05b_S58","cFT09a_S53","cFT09b_S59","cFT13a_S48","cFT13b_S54",
                "cFT17a_S49","cFT17b_S55","cFT21a_S50","cFT21b_S56","cNA01a_S4",
                "cNA01b_S10","cNA01c_S16","cNA05a_S5","cNA05b_S11","cNA05c_S17",
                "cNA09a_S6","cNA09b_S12","cNA09c_S18","cNA13a_S1","cNA13b_S7",
                "cNA13c_S13","cNA17a_S2","cNA17b_S8","cNA17c_S14","cNA21a_S3",
                "cNA21b_S9","cNA21c_S15")
  dat<-fread(dt)%>%dplyr::select(all_of(samps_list))
  names(dat) <- c("FeatureID","cFA01a","cFA01b","cFA05a","cFA05b","cFA09a","cFA09b","cFA13a","cFA13b","cFA17a","cFA17b",
                  "cFA21a","cFA21b","cFT01a","cFT01b","cFT05a","cFT05b","cFT09a","cFT09b","cFT13a","cFT13b",
                  "cFT17a","cFT17b","cFT21a","cFT21b","cNA01a","cNA01b","cNA01c","cNA05a","cNA05b","cNA05c",
                  "cNA09a","cNA09b","cNA09c","cNA13a","cNA13b","cNA13c","cNA17a","cNA17b","cNA17c","cNA21a",
                  "cNA21b","cNA21c")
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
         y =paste("PC2 (",round(ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""))+ggtitle(paste("MTX",cond,sep=" "))+ theme(plot.title = element_text(face = "bold"))
  
  return(p)
}
###########################################################
#clean metadata
md<-fread(dat_metadata_ori)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition),
         phase=ifelse(zt_time<13,"light","dark"),
         zt_time_cat=paste("ZT",zt_time, sep=""),
         cond_phase=paste(condition,phase,sep="_"),
         cond_zt=paste(condition,zt_time,sep="_"),
         fasted=ifelse(zt_time_cat %in% c("ZT5","ZT9","ZT13"), "fasted","fed"),
         cond_fasted=paste(condition,fasted,sep="_"))
write.table(md, dat_metadata, sep = "\t",row.names = FALSE, quote=FALSE)

#clean annotations
annotations<-fread(notnorm_data) %>%dplyr::select(c(1,61))%>%
  dplyr::rename(FeatureID="#FeatureID")
write.csv(annotations,paste0(dat_path,"pfam_annotationkey.csv"), row.names = FALSE)

#subset data
dat_nn<-subset_dat(notnorm_data)
write.table(dat_nn,paste0(dat_path,"pfam_clean_noNT.tsv"),sep = "\t",row.names = FALSE, quote=FALSE)

#####################################################################################
#used emperor to plot RPCA--Figure 1D
#permamova values are in condition-significance.qzv files 
#need to run rna_rpca_metaT.sh script first

##########################################################
#plot the ZT time for each condition--Figure 1E-G (right)
#need to run rna_rpca_metaT.sh script first

ord <- read_qza(dat_ord)

samp_ord<-ord$data$Vectors
write.table(samp_ord,paste0(res_path,"sample_ordination.txt"),sep = "\t",row.names = FALSE,quote=FALSE)
feat_ord<-ord$data$Species
write.table(feat_ord,paste0(res_path,"feature_ordination.txt"),sep = "\t",row.names = FALSE,quote=FALSE)

md<-fread(dat_metadata)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))

rpca<-ord$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2, PC3)%>%
  dplyr::rename(sample_name=SampleID)%>%
  left_join(md,by="sample_name")%>%
  mutate(condition=factor(condition,levels=c("FT","FA","NA")),
         phase=factor(phase,levels=c("light","dark")))

p<-rpcacond_plt(rpca,"NA")
ggsave(paste0(fig_path,"SFR24_0412_mtx_NA_RPCA.pdf"), plot=p,height=3.5, width=3.5)

p<-rpcacond_plt(rpca,"FA")
ggsave(paste0(fig_path,"SFR24_0412_mtx_FA_RPCA.pdf"), plot=p,height=3.5, width=3.5)

p<-rpcacond_plt(rpca,"FT")
ggsave(paste0(fig_path,"SFR24_0412_mtx_FT_RPCA.pdf"), plot=p,height=3.5, width=3.5)

#permamova values are in cond_fasted-significance.qzv and condphase-significance.qzv files 