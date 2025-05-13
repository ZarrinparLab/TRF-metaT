setwd("~/Notebooks/sfloresr/TRF-metaT/")

library(tidyverse)
library(data.table)
library("qiime2R")
library(ggpubr)

#####################################################################################
#paths
dat_metadata_ori<-"data/metadata.TRF_combined.tab"
dat_metadata<-"data/metadata.TRF_combined_wLD.tab"
dat_ord<-"data/diversity_analysis/rpca_results_m16s/ordination.qza"
res_path<-"data/diversity_analysis/rpca_results_m16s/"
fig_path<-"figures/diversity_analysis/"
#####################################################################################
#functions
rpcacond_plt<-function(dat,cond){
  
  rpca_sub<-dat%>%
    filter(condition==cond)%>%
    mutate(zt=factor(zt,levels=c("1","5","9","13","17","21")))
  
  p<-rpca_sub %>%
    ggplot(aes(x=PC1, y=PC2, fill=zt)) +
    geom_point(alpha=1.0,size=3,shape=21) + 
    theme_pubr() +
    scale_fill_manual(values=c("#67001f","#d6604d","#ffba92","#8ccff3","#4393c3","#053061"))+
    labs(color="ZT time",
         x =paste("PC1 (",round(ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
         y =paste("PC2 (",round(ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""))+ggtitle(paste("16S",cond,sep=" "))+ theme(plot.title = element_text(face = "bold"))
  
  return(p)
}
#####################################################################################

#clean metadata
md<-fread(dat_metadata_ori)%>%
  mutate(phase=ifelse(zt<13,"light","dark"),
         cond_phase=paste(condition,phase,sep="_"),
         fasted=ifelse(zt %in% c(5,9,13), "fasted","fed"),
         cond_fasted=paste(condition,fasted,sep="_"))
write.table(md,dat_metadata,sep = "\t",row.names = FALSE, quote=FALSE)

#####################################################################################
#used emperor to plot RPCA--Figure 1B
#permamova values are in condition-significance.qzv files 
#need to run rna_rpca_m16s.sh script first

#####################################################################################
#plot the ZT time for each condition--Figure 1E-G (left)
#need to run rna_rpca_m16s.sh script first

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
  filter(condition!="NT")%>%
  mutate(condition=factor(condition,levels=c("FT","FA","NA")),
         phase=factor(phase,levels=c("light","dark")))

p<-rpcacond_plt(rpca,"NA")
ggsave(paste0(fig_path,"SFR24_0412_m16s_NA_RPCA.pdf"), plot=p,height=3.5, width=3.5)

p<-rpcacond_plt(rpca,"FA")
ggsave(paste0(fig_path,"SFR24_0412_m16s_FA_RPCA.pdf"), plot=p,height=3.5, width=3.5)

p<-rpcacond_plt(rpca,"FT")
ggsave(paste0(fig_path,"SFR24_0412_m16s_FT_RPCA.pdf"), plot=p,height=3.5, width=3.5)

#permamova values are in cond_fasted-significance.qzv and condphase-significance.qzv files 