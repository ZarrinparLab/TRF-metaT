setwd("~/Notebooks/sfloresr/TRF-metaT/")

library(MetaCycle)
library(tidyverse)
library(data.table)
library(ggpubr)
library(RColorBrewer)

##########################################################
#paths
dat_annot<-"data/pfam_metaG/pfam_annotationkey.csv"
dat_mgx_FTcyc<-"data/cycling_analysis/mgx_cyclic_analysis/wol2_dna_pipeline_rpob/FT_metacycle/meta2d_filtered_dna_FT.txt"
dat_mgx_FAcyc<-"data/cycling_analysis/mgx_cyclic_analysis/wol2_dna_pipeline_rpob/FA_metacycle/meta2d_filtered_dna_FA.txt"
dat_mgx_NAcyc<-"data/cycling_analysis/mgx_cyclic_analysis/wol2_dna_pipeline_rpob/NA_metacycle/meta2d_filtered_dna_NA.txt"
fig_path<-"figures/cycling_analysis/"
##########################################################
#functions

load_cyc_dat<-function(dat){
  cyc_dat<-fread(dat)%>%
    dplyr::rename(FeatureID=CycID) %>%
    left_join(.,pfam_annot, by="FeatureID")%>%
    arrange(JTK_adjphase)%>%
    mutate(label_name=paste(FeatureID, Name, sep=" "))
  return(cyc_dat)
}

get_summ_dat<-function(d1,d2,d3,d4){
  cyc_summ<-data.frame(condition=c("NA","FA","FT"),
                       cycling_n=c(nrow(d1),nrow(d2),nrow(d3)),
                       not_cycling_n=c(nrow(d4)-nrow(d1),nrow(d4)-nrow(d2),nrow(d4)-nrow(d3)))%>%
    mutate(cycling_perc=(cycling_n/nrow(d4))*100,
           not_cycling_perc=(not_cycling_n/nrow(d4))*100)
  return(cyc_summ)
}

get_summ_plt<-function(dat){
  dat_plot<-dat%>%gather(cycling_grp, value, -condition)%>%
    filter(cycling_grp=="cycling_perc"|cycling_grp=="not_cycling_perc")%>%
    mutate(condition=factor(condition,levels=c("NA","FA","FT")))%>%
    filter(cycling_grp=="cycling_perc")
  
  p<-ggplot(data=dat_plot, aes(x=condition, y=value, fill=condition)) +
    geom_bar(stat="identity") + theme_classic() +
    theme(legend.position = "top")+
    scale_y_continuous(expand=c(0,0), limits=c(0,10)) +
    labs(y="transcripts (%)") +
    scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))
  return(p)
}

get_chi<-function(dat,compar){
  if(compar=="NAFA"){
    mtx_chi<-dat%>%dplyr::select(1:3)%>%
      filter(condition!="FT")%>%
      column_to_rownames("condition")%>%as.matrix()
  }
  else if (compar=="FAFT"){
    mtx_chi<-dat%>%dplyr::select(1:3)%>%
      filter(condition!="NA")%>%
      column_to_rownames("condition")%>%as.matrix()
  }
  else{
    mtx_chi<-dat%>%dplyr::select(1:3)%>%
      filter(condition!="FA")%>%
      column_to_rownames("condition")%>%as.matrix()
  }
  res<-chisq.test(mtx_chi,simulate.p.value=TRUE, B=2000) 
  return(res)
}
##########################################################
#load MGX cycling files

#annotations
pfam_annot<-fread(dat_annot)

FT_metacyc<-load_cyc_dat(dat_mgx_FTcyc)
FA_metacyc<-load_cyc_dat(dat_mgx_FAcyc)
NA_metacyc<-load_cyc_dat(dat_mgx_NAcyc)

##########################################################
#make a bar plot summarizing the number of cycling and non-cycling hits (removing DUF)
# for mgx-- Figure 3A (left)

sigFT<-FT_metacyc%>%filter(JTK_pvalue<0.05,!grepl("DUF",Name)) #45
sigFA<-FA_metacyc%>%filter(JTK_pvalue<0.05,!grepl("DUF",Name)) #42
sigNA<-NA_metacyc%>%filter(JTK_pvalue<0.05,!grepl("DUF",Name)) #308
noduf<-FT_metacyc%>%filter(!grepl("DUF",Name)) #3473

cyc_summ<-get_summ_dat(sigNA,sigFA,sigFT,noduf)
p<-get_summ_plt(cyc_summ)
ggsave(paste0(fig_path,"SFR23_0619_MG_cycling_summ_barplot.pdf"), plot=p,height=3, width=2.5)

get_chi(cyc_summ,"NAFA")
#X-squared = 212.89, df = NA, p-value = 0.0004998 NA FA
get_chi(cyc_summ,"FAFT")
#X-squared = 0.10476, df = NA, p-value = 0.8281 no FA vs FT diff
get_chi(cyc_summ,"NAFT")
#X-squared = 206.44, df = NA, p-value = 0.0004998 NA FT

##########################################################
#make venn diagram showing overlap of whats cycling for mgx --Figure 3B (left)

list_venn <- list(NA_ = sigNA$FeatureID,
                  FA = sigFA$FeatureID,
                  FT = sigFT$FeatureID)

draw.venn(sigNA$FeatureID,sigFA$FeatureID,sigFT$FeatureID,
          title="",subtitle="",
          xtitle="NA",ytitle="FA",ztitle="FT",
          xt_s=3,yt_s=3,zt_s=3, nr_s=3,
          x_c="#0072B2",y_c="#D55E00",z_c="#009E73",
          output="pdf",
          filename=paste0(fig_path,"SFR23_0619_MG_venn_overlap.pdf"))
