setwd("~/Notebooks/sfloresr/TRF-metaT/")

library(MetaCycle)
library(tidyverse)
library(data.table)
library(ggpubr)
library(gplots)
library(BioVenn)
library(RColorBrewer)

##########################################################
#paths
dat_annot<-"data/pfam_metaT/pfam_annotationkey.csv"
dat_metadata<-"data/metaT_metadata_ztcat_noNT.txt"
dat_go<-"data/pfam_metaT/go_name.txt"
dat_pfamtogo<-"data/pfam_metaT/pfam-to-go-process.map"
dat_mtx_FTcyc_144k<-"data/cycling_analysis/mtx_cyclic_analysis/wol2_rna_pipeline_144k_rpob/FT_metacycle/meta2d_filtered_rna_FT.txt"
dat_mtx_FAcyc_144k<-"data/cycling_analysis/mtx_cyclic_analysis/wol2_rna_pipeline_144k_rpob/FA_metacycle/meta2d_filtered_rna_FA.txt"
dat_mtx_NAcyc_144k<-"data/cycling_analysis/mtx_cyclic_analysis/wol2_rna_pipeline_144k_rpob/NA_metacycle/meta2d_filtered_rna_NA.txt"
dat_mtx_FTcyc<-"data/cycling_analysis/mtx_cyclic_analysis/wol2_rna_pipeline_rpob/FT_metacycle/meta2d_filtered_rna_FT.txt"
dat_mtx_FAcyc<-"data/cycling_analysis/mtx_cyclic_analysis/wol2_rna_pipeline_rpob/FA_metacycle/meta2d_filtered_rna_FA.txt"
dat_mtx_NAcyc<-"data/cycling_analysis/mtx_cyclic_analysis/wol2_rna_pipeline_rpob/NA_metacycle/meta2d_filtered_rna_NA.txt"
dat_mtx_nrpob<-"data/pfam_metaT/pfam-TPM_clean_noNT_normRPOB.txt"
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

get_pfamtogo<-function(datpfamgo){
  pfamGOp<-read.table(datpfamgo,header = FALSE, sep = "\t",
                      col.names = paste0("V",seq_len(4)), fill = TRUE)%>%
    gather(column,GO_Term,-V1)%>%
    dplyr::select(1,3)%>%
    dplyr::rename(FeatureID=V1)%>%
    left_join(.,gonames,by="GO_Term")%>%
    filter(!is.na(name))
  return(pfamGOp)
}

stderror <- function(x) sd(x)/sqrt(length(x))

ZT_dist <- function(data,id) {
  
  data_sub<-data%>%filter(FeatureID==id)
  p<-data_sub%>%
    ggplot(aes(x=zt_time, y=log10(mn_TPM+1), color=condition)) +
    geom_point(alpha=1.0) + geom_line() +
    theme_pubr() +
    scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+
    scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
    scale_x_continuous(breaks=c(1,5,9,13,17,21))+
    geom_ribbon(aes(ymin = log10(mn_TPM-sem+1), ymax = log10(mn_TPM+sem+1), fill=condition),alpha=0.3,colour = NA)+
    labs(color="condition",
         y ="log10(avgTPM)",
         x ="ZT time",
         title=id)+ theme(plot.title = element_text(face = "bold"))
  return(p)
}

##########################################################
#load MTX cycling files

#annotations
pfam_annot<-fread(dat_annot)

#rarefied
FT_metacyc_144k<-load_cyc_dat(dat_mtx_FTcyc_144k)
FA_metacyc_144k<-load_cyc_dat(dat_mtx_FAcyc_144k)
NA_metacyc_144k<-load_cyc_dat(dat_mtx_NAcyc_144k)

#not rarefied
FT_metacyc<-load_cyc_dat(dat_mtx_FTcyc)
FA_metacyc<-load_cyc_dat(dat_mtx_FAcyc)
NA_metacyc<-load_cyc_dat(dat_mtx_NAcyc)

##########################################################
#make a bar plot summarizing the number of cycling and 
#non-cycling hits (removing DUF) for rarefied mtx--Figure S3A

sigFT_144k<-FT_metacyc_144k%>%filter(JTK_pvalue<0.05,!grepl("DUF",Name)) #122
sigFA_144k<-FA_metacyc_144k%>%filter(JTK_pvalue<0.05,!grepl("DUF",Name)) #56
sigNA_144k<-NA_metacyc_144k%>%filter(JTK_pvalue<0.05,!grepl("DUF",Name)) #292
noduf_144k<-FT_metacyc_144k%>%filter(!grepl("DUF",Name)) #4501

cyc_summ<-get_summ_dat(sigNA_144k,sigFA_144k,sigFT_144k,noduf_144k)
p<-get_summ_plt(cyc_summ)
ggsave(paste0(fig_path,"SFR23_0620_MT_cycling_summ_144k_barplot.pdf"), plot=p,height=3, width=2.5)

get_chi(cyc_summ,"NAFA")
#X-squared = 166.48, df = NA, p-value = 0.0004998 FA NA
get_chi(cyc_summ,"FAFT")
#X-squared = 24.966, df = NA, p-value = 0.0004998 FA FT
get_chi(cyc_summ,"NAFT")
#X-squared = 73.172, df = NA, p-value = 0.0004998 FT NA

##########################################################

#make venn diagram showing overlap of whats cycling for rarefied mtx--Figure S3B

list_venn <- list(NA_ = sigNA_144k$FeatureID,
                  FA = sigFA_144k$FeatureID,
                  FT = sigFT_144k$FeatureID)

draw.venn(sigNA_144k$FeatureID,sigFA_144k$FeatureID,sigFT_144k$FeatureID,
          title="",subtitle="",
          xtitle="NA",ytitle="FA",ztitle="FT",
          xt_s=3,yt_s=3,zt_s=3, nr_s=3,
          x_c="#0072B2",y_c="#D55E00",z_c="#009E73",
          output="pdf",
          filename=paste0(fig_path,"SFR23_0620_MT_144k_venn_overlap.pdf"))

##########################################################
#make a bar plot summarizing the number of cycling and non-cycling hits (removing DUF)
# for non-rarefied mtx-- Figure 3A (right)

sigFT<-FT_metacyc%>%filter(JTK_pvalue<0.05,!grepl("DUF",Name)) #662
sigFA<-FA_metacyc%>%filter(JTK_pvalue<0.05,!grepl("DUF",Name)) #117
sigNA<-NA_metacyc%>%filter(JTK_pvalue<0.05,!grepl("DUF",Name)) #571
noduf<-FT_metacyc%>%filter(!grepl("DUF",Name)) #7113

cyc_summ<-get_summ_dat(sigNA,sigFA,sigFT,noduf)
p<-get_summ_plt(cyc_summ)
ggsave(paste0(fig_path,"SFR23_0619_MT_cycling_summ_barplot.pdf"), plot=p,height=3, width=2.5)

get_chi(cyc_summ,"NAFA")
#X-squared = 314.81, df = NA, p-value = 0.0004998 FA NA
get_chi(cyc_summ,"FAFT")
#X-squared = 403.38, df = NA, p-value = 0.0004998 FA FT
get_chi(cyc_summ,"NAFT")
#X-squared = 7.3535, df = NA, p-value = 0.008996 NA FT

##########################################################
#make venn diagram showing overlap of whats cycling for non-rarefied mtx --Figure 3B (right)

list_venn <- list(NA_ = sigNA$FeatureID,
                  FA = sigFA$FeatureID,
                  FT = sigFT$FeatureID)

draw.venn(sigNA$FeatureID,sigFA$FeatureID,sigFT$FeatureID,
          title="",subtitle="",
          xtitle="NA",ytitle="FA",ztitle="FT",
          xt_s=3,yt_s=3,zt_s=3, nr_s=3,
          x_c="#0072B2",y_c="#D55E00",z_c="#009E73",
          output="pdf",
          filename=paste0(fig_path,"SFR23_0619_MT_venn_overlap.pdf"))

##########################################################
##see cycling_comb.R script for Figure 3C and S3C
##########################################################
#create heatmap describing shared patterns--Figure 3D

ItemsList <- venn(list_venn, show.plot = FALSE)
all<-attributes(ItemsList)$intersections

shared_cyc<-plyr::ldply (all, data.frame)%>%
  dplyr::rename(pattern=".id", FeatureID="X..i..")

tworhythmic<-shared_cyc%>%filter(pattern %in% c("NA_:FA:FT","FA:FT","NA_:FA","NA_:FT"))

md<-fread(dat_metadata)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))
go_annot<-fread(dat_go)%>%
  dplyr::rename(FeatureID=GO_Term)

pfamTPM_all<-fread(dat_mtx_nrpob)%>%
  filter(FeatureID %in% tworhythmic$FeatureID)%>%
  gather(sample_name,TPM_counts,-FeatureID) %>%
  left_join(.,pfam_annot, by ="FeatureID") %>%
  mutate(label_name=paste(FeatureID, Name, sep=" "))%>%
  left_join(.,md,by="sample_name") %>%
  left_join(.,tworhythmic,by="FeatureID")%>%
  group_by(label_name,condition,zt_time, pattern)%>%dplyr::summarise(mn_TPM=mean(TPM_counts))%>%
  group_by(label_name,condition)%>%mutate(Zscore=(mn_TPM - mean(mn_TPM))/sd(mn_TPM))%>%
   mutate(zt_time=factor(zt_time, levels = c("1","5","9","13","17","21")),
          condition=factor(condition,levels=c("NA","FA","FT")),
          pattern=factor(pattern,levels=c("NA_:FA:FT","FA:FT","NA_:FA","NA_:FT")))

#get custom ordering of JTK_adjphase

NA_metacyc_woFTFA<-NA_metacyc%>%filter(!(FeatureID %in% all$`FA:FT`))
FT_metacyc_justFTFA<-FT_metacyc%>%filter(FeatureID %in% all$`FA:FT`)
orderlist<-rbind(FT_metacyc_justFTFA,NA_metacyc_woFTFA)

pfamTPM_all$label_name <- factor(pfamTPM_all$label_name,levels = orderlist$label_name)

plt<-ggplot(pfamTPM_all,aes(x=zt_time, y=label_name)) +theme_classic()+
  geom_tile(aes(fill=Zscore))+
  scale_x_discrete(expand = c(0, 0))+
  facet_grid(pattern~condition,scales="free",space="free")+
  theme(axis.ticks.y=element_blank(),panel.spacing.x=unit(1, "lines"),axis.text.y = element_text(size = 4),
        panel.spacing.y=unit(0.3, "lines"),strip.text.y = element_text(angle = 0))+
  scale_fill_viridis(option="inferno") + 
  xlab("ZT Time") + ylab("transcripts")+labs(fill='TPM Z-score')+ggtitle("")

ggsave(paste0(fig_path,"SFR23_0619_MT_heatmap_byphase_RRN_zbycond.pdf"), plot=plt,height=4, width=9)

##########################################################
#find what cycling is shared by FT FA and NA and FT--Figure 3E-F

ItemsList <- venn(list_venn, show.plot = FALSE)
all<-attributes(ItemsList)$intersections

FTNA_shared_cyc<-pfam_annot%>%filter(FeatureID %in% all$`NA_:FT`) #86
FTFA_shared_cyc<-pfam_annot%>%filter(FeatureID %in% all$`FA:FT`) #38
NA_cyc<-pfam_annot%>%filter(FeatureID %in% all$`NA_`) #470
FA_cyc<-pfam_annot%>%filter(FeatureID %in% all$`FA`) #64
FT_cyc<-pfam_annot%>%filter(FeatureID %in% all$`FT`) #534

#FT FA
FT_38sharedcyc<-sigFT%>%filter(FeatureID %in% all$`FA:FT`)%>%
  dplyr::select(FeatureID,JTK_adjphase)%>%
  mutate(condition="FT")

FA_38sharedcyc<-sigFA%>%filter(FeatureID %in% all$`FA:FT`)%>%
  dplyr::select(FeatureID,JTK_adjphase)%>%
  mutate(condition="FA")

comb_38sharedcyc<-rbind(FT_38sharedcyc,FA_38sharedcyc)%>%
  mutate(condition=factor(condition,levels=c("FA","FT")))

ggplot(comb_38sharedcyc, aes(x=JTK_adjphase, color=condition,fill=condition)) +
  geom_density(alpha=0.3)+ theme_pubr()+
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_color_manual(values=c("#D55E00","#009E73"))+
  scale_fill_manual(values=c("#D55E00","#009E73")) 

ggsave(paste0(fig_path,"SFR23_0619_MT_densityPhasenoFA_38cycFTFA.pdf"), width = 4.5, height = 4) #Fig 3E

#FT NA
FT_86sharedcyc<-sigFT%>%filter(FeatureID %in% all$`NA_:FT`)%>%
  dplyr::select(FeatureID,JTK_adjphase)%>%
  mutate(condition="FT")

NA_86sharedcyc<-sigNA%>%filter(FeatureID %in% all$`NA_:FT`)%>%
  dplyr::select(FeatureID,JTK_adjphase)%>%
  mutate(condition="NA")

comb_86sharedcyc<-rbind(FT_86sharedcyc,NA_86sharedcyc)%>%
  mutate(condition=factor(condition,levels=c("NA","FT")))

ggplot(comb_86sharedcyc, aes(x=JTK_adjphase, color=condition,fill=condition)) +
  geom_density(alpha=0.3)+ theme_pubr()+
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_color_manual(values=c("#0072B2","#009E73"))+
  scale_fill_manual(values=c("#0072B2","#009E73"))

ggsave(paste0(fig_path,"SFR23_0619_MT_densityPhasenoFA_86cycFTNA.pdf"), width = 4.5, height = 4) #Fig 3F

##########################################################
#look at what type of cycling hits were just in NA--Figure 3G

gonames<-fread(dat_go)

NA_cyc<-pfam_annot%>%filter(FeatureID %in% all$NA_) #470

pfamGOp_NA<-get_pfamtogo(dat_pfamtogo)%>%
  filter(FeatureID %in% NA_cyc$FeatureID)

NAcyc_summ<-pfamGOp_NA%>%
  group_by(name)%>%summarise(n=n())%>%
  mutate(condition="NA")%>%
  filter(n>1)%>%
  arrange(n)

NAcyc_summ$name <- factor(NAcyc_summ$name,levels = NAcyc_summ$name )

ggplot(data=NAcyc_summ, aes(x=name, y=n, fill=condition)) +
  geom_bar(stat="identity", position=position_dodge()) + coord_flip() +
  scale_fill_manual(values=c("#0072B2")) +theme_pubr() +
  scale_y_continuous(expand=c(0,0), limits=c(0,20))

ggsave(paste0(fig_path,"SFR23_0619_MT_cycNA_GOterms.pdf"),height=3.5, width=6)
##########################################################
#look at what type of cycling hits were just in FT--Figure 3I

FT_cyc<-pfam_annot%>%filter(FeatureID %in% all$FT) #534

pfamGOp_FT<-get_pfamtogo(dat_pfamtogo)%>%
  filter(FeatureID %in% FT_cyc$FeatureID)

FTcyc_summ<-pfamGOp_FT%>%
  group_by(name)%>%summarise(n=n())%>%
  mutate(condition="FT")%>%
  filter(n>2)%>%
  arrange(n)

FTcyc_summ$name <- factor(FTcyc_summ$name,levels = FTcyc_summ$name )

ggplot(data=FTcyc_summ, aes(x=name, y=n, fill=condition)) +
  geom_bar(stat="identity", position=position_dodge())  +coord_flip() +
  scale_fill_manual(values=c("#009E73")) +theme_pubr() +
  scale_y_continuous(expand=c(0,0), limits=c(0,20))

ggsave(paste0(fig_path,"SFR23_0619_MT_cycFT_GOterms.pdf"),height=3.5, width=6)
##########################################################
#look at what type of cycling hits were just in FA--Figure S3F

FA_cyc<-pfam_annot%>%filter(FeatureID %in% all$FA) #64

pfamGOp_FA<-get_pfamtogo(dat_pfamtogo)%>%
  filter(FeatureID %in% FA_cyc$FeatureID)

FAcyc_summ<-pfamGOp_FA%>%
  group_by(name)%>%summarise(n=n())%>%
  mutate(condition="FA")%>%
  arrange(n)

FAcyc_summ$name <- factor(FAcyc_summ$name,levels = FAcyc_summ$name )

ggplot(data=FAcyc_summ, aes(x=name, y=n, fill=condition)) +
  geom_bar(stat="identity", position=position_dodge()) + coord_flip() +
  scale_fill_manual(values=c("#D55E00")) +theme_pubr() +
  scale_y_continuous(expand=c(0,0), limits=c(0,20))

ggsave(paste0(fig_path,"SFR23_0619_MT_cycFA_GOterms.pdf"),height=3, width=8)
##########################################################
#plot specific examples over ZT time--Figure 3H-J, S3D-E,G

pfamZT<-fread(dat_mtx_nrpob)%>%
  gather(sample_name,TPM_counts,-FeatureID) %>%
  left_join(.,pfam_annot, by ="FeatureID") %>%
  mutate(label_name=paste(FeatureID, Name, sep=" "))%>%
  left_join(.,md,by="sample_name") %>%
  group_by(FeatureID,condition,zt_time)%>%summarise(mn_TPM=mean(TPM_counts),sem=stderror(TPM_counts))%>%
  mutate(condition=factor(condition, levels = c("NA","FA","FT")))

#functions that cycled in both FA and FT--Fig S3D
p<-ZT_dist(pfamZT,"PF00251.23")
ggsave(paste0(fig_path,"FTFA_shared_cyc_MT_zt/PF00251.23.pdf"), p, height=3, width=3)
p<-ZT_dist(pfamZT,"PF04463.15")
ggsave(paste0(fig_path,"FTFA_shared_cyc_MT_zt/PF04463.15.pdf"), p, height=3, width=3)
p<-ZT_dist(pfamZT,"PF08244.15")
ggsave(paste0(fig_path,"FTFA_shared_cyc_MT_zt/PF08244.15.pdf"), p, height=3, width=3)

#functions that cycled in both NA and FT--Fig S3E
p<-ZT_dist(pfamZT,"PF02449.18")
ggsave(paste0(fig_path,"FTNA_shared_cyc_MT_zt/PF02449.18.pdf"), p, height=3, width=3)
p<-ZT_dist(pfamZT,"PF03808.16")
ggsave(paste0(fig_path,"FTNA_shared_cyc_MT_zt/PF03808.16.pdf"), p, height=3, width=3)
p<-ZT_dist(pfamZT,"PF10437.12")
ggsave(paste0(fig_path,"FTNA_shared_cyc_MT_zt/PF10437.12.pdf"), p, height=3, width=3)

#functions that cycled just in NA--Fig 3H
p<-ZT_dist(pfamZT,"PF02770.22")
ggsave(paste0(fig_path,"NA_cyc_MT_zt/PF02770.22.pdf"), p, height=3, width=3)
p<-ZT_dist(pfamZT,"PF19277.2")
ggsave(paste0(fig_path,"NA_cyc_MT_zt/PF19277.2.pdf"), p, height=3, width=3)

#functions that cycled just in FT--Fig 3I
p<-ZT_dist(pfamZT,"PF01219.22")
ggsave(paste0(fig_path,"FT_cyc_MT_zt/PF01219.22.pdf"), p, height=3, width=3)
p<-ZT_dist(pfamZT,"PF13561.9")
ggsave(paste0(fig_path,"FT_cyc_MT_zt/PF13561.9.pdf"), p, height=3, width=3)

#functions that cycled just in FA--Fig S3G
p<-ZT_dist(pfamZT,"PF05139.17")
ggsave(paste0(fig_path,"FA_cyc_MT_zt/PF05139.17.pdf"), p, height=3, width=3)
p<-ZT_dist(pfamZT,"PF09704.13")
ggsave(paste0(fig_path,"FA_cyc_MT_zt/PF09704.13.pdf"), p, height=3, width=3)
