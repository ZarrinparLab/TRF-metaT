setwd("~/scratch/TRF_multiomics/metagenomic/woltka2_results/filtered_metaG/")

library(MetaCycle)
library(tidyverse)
library(data.table)
library(ggpubr)
library(ggvenn)
library(gplots)
library(viridis)
library(RColorBrewer)
library(BioVenn) 
library("qiime2R")

##########################################################
#MG-species

#annotations
pfam_annot<-fread("pfam_notnorm/pfam_annotationkey.csv")

#read in FT metacycle hits
#FT_metacyc<-fread("cyclic_analysis/wol2_dna_pipeline_ogu-rclr/FT_metacycle/meta2d_filtered_dna_FT.txt")%>%
FT_metacyc<-fread("cyclic_analysis/wol2_dna_pipeline_ogu-prop/FT_metacycle/meta2d_filtered_dna_FT.txt")%>%
#FT_metacyc<-fread("cyclic_analysis/wol2_dna_pipeline_species-rclr/FT_metacycle/meta2d_filtered_dna_FT.txt")%>%
#FT_metacyc<-fread("cyclic_analysis/wol2_dna_pipeline_species-prop/FT_metacycle/meta2d_filtered_dna_FT.txt")%>%
#FT_metacyc<-fread("cyclic_analysis/wol2_dna_pipeline_genus-rclr/FT_metacycle/meta2d_filtered_dna_FT.txt")%>%
#FT_metacyc<-fread("cyclic_analysis/wol2_dna_pipeline_genus-rclr/FT_metacycle/meta2d_filtered_dna_FT.txt")%>%
  dplyr::rename(FeatureID=CycID) %>%
  left_join(.,pfam_annot, by="FeatureID")%>%
  arrange(JTK_adjphase)%>%
  mutate(label_name=paste(FeatureID, Name, sep=" "))

#read in FA metacycle hits
#FA_metacyc<-fread("cyclic_analysis/wol2_dna_pipeline_ogu-rclr/FA_metacycle/meta2d_filtered_dna_FA.txt")%>%
FA_metacyc<-fread("cyclic_analysis/wol2_dna_pipeline_ogu-prop/FA_metacycle/meta2d_filtered_dna_FA.txt")%>%
#FA_metacyc<-fread("cyclic_analysis/wol2_dna_pipeline_species-rclr/FA_metacycle/meta2d_filtered_dna_FA.txt")%>%
#FA_metacyc<-fread("cyclic_analysis/wol2_dna_pipeline_species-prop/FA_metacycle/meta2d_filtered_dna_FA.txt")%>%
#FA_metacyc<-fread("cyclic_analysis/wol2_dna_pipeline_genus-rclr/FA_metacycle/meta2d_filtered_dna_FA.txt")%>%
#FA_metacyc<-fread("cyclic_analysis/wol2_dna_pipeline_genus-prop/FA_metacycle/meta2d_filtered_dna_FA.txt")%>%
  dplyr::rename(FeatureID=CycID) %>%
  left_join(.,pfam_annot, by="FeatureID")%>%
  arrange(JTK_adjphase)%>%
  mutate(label_name=paste(FeatureID, Name, sep=" "))

#read in NA metacycle hits
#NA_metacyc<-fread("cyclic_analysis/wol2_dna_pipeline_ogu-rclr/NA_metacycle/meta2d_filtered_dna_NA.txt")%>%
NA_metacyc<-fread("cyclic_analysis/wol2_dna_pipeline_ogu-prop/NA_metacycle/meta2d_filtered_dna_NA.txt")%>%
#NA_metacyc<-fread("cyclic_analysis/wol2_dna_pipeline_species-rclr/NA_metacycle/meta2d_filtered_dna_NA.txt")%>%
#NA_metacyc<-fread("cyclic_analysis/wol2_dna_pipeline_species-prop/NA_metacycle/meta2d_filtered_dna_NA.txt")%>%
#NA_metacyc<-fread("cyclic_analysis/wol2_dna_pipeline_genus-rclr/NA_metacycle/meta2d_filtered_dna_NA.txt")%>%
#NA_metacyc<-fread("cyclic_analysis/wol2_dna_pipeline_genus-prop/NA_metacycle/meta2d_filtered_dna_NA.txt")%>%
  dplyr::rename(FeatureID=CycID) %>%
  left_join(.,pfam_annot, by="FeatureID")%>%
  arrange(JTK_adjphase)%>%
  mutate(label_name=paste(FeatureID, Name, sep=" "))

#make a bar plot summarizing the number of cycling and non-cycling hits (removing DUF)

sigFT<-FT_metacyc%>%filter(JTK_pvalue<0.1) #9 8 4 | 10 3 5
sigFA<-FA_metacyc%>%filter(JTK_pvalue<0.1) #18 16 10 | 11 15 12
sigNA<-NA_metacyc%>%filter(JTK_pvalue<0.1) #35 33 22 | 43 41 25

# cyc_summ<-data.frame(condition=c("NA","FA","FT"),
#                      cycling_n=c(35,18,9),
#                      not_cycling_n=c(570-35,570-18,570-9))%>%
#   mutate(cycling_perc=(cycling_n/570)*100,
#          not_cycling_perc=(not_cycling_n/570)*100)
# write.table(cyc_summ,"cyclic_analysis/wol2_dna_pipeline_ogu-rclr/SFR23_0823_MG_cycling_summary_table.txt",sep = "\t",row.names = FALSE, quote=FALSE)

cyc_summ<-data.frame(condition=c("NA","FA","FT"),
                     cycling_n=c(43,11,10),
                     not_cycling_n=c(570-43,570-11,570-10))%>%
  mutate(cycling_perc=(cycling_n/570)*100,
         not_cycling_perc=(not_cycling_n/570)*100)
write.table(cyc_summ,"cyclic_analysis/wol2_dna_pipeline_ogu-prop/SFR23_0828_MG_cycling_summary_table.txt",sep = "\t",row.names = FALSE, quote=FALSE)


# cyc_summ<-data.frame(condition=c("NA","FA","FT"),
#                      cycling_n=c(33,16,8),
#                      not_cycling_n=c(473-33,473-16,473-8))%>%
#   mutate(cycling_perc=(cycling_n/473)*100,
#          not_cycling_perc=(not_cycling_n/473)*100)
# write.table(cyc_summ,"cyclic_analysis/wol2_dna_pipeline_species-rclr/SFR23_0823_MG_cycling_summary_table.txt",sep = "\t",row.names = FALSE, quote=FALSE)

# cyc_summ<-data.frame(condition=c("NA","FA","FT"),
#                      cycling_n=c(41,15,3),
#                      not_cycling_n=c(473-41,473-15,473-3))%>%
#   mutate(cycling_perc=(cycling_n/473)*100,
#          not_cycling_perc=(not_cycling_n/473)*100)
# write.table(cyc_summ,"cyclic_analysis/wol2_dna_pipeline_species-prop/SFR23_0828_MG_cycling_summary_table.txt",sep = "\t",row.names = FALSE, quote=FALSE)

# cyc_summ<-data.frame(condition=c("NA","FA","FT"),
#                      cycling_n=c(22,10,4),
#                      not_cycling_n=c(262-22,262-10,262-4))%>%
#   mutate(cycling_perc=(cycling_n/262)*100,
#          not_cycling_perc=(not_cycling_n/262)*100)
# write.table(cyc_summ,"cyclic_analysis/wol2_dna_pipeline_genus-rclr/SFR23_0823_MG_cycling_summary_table.txt",sep = "\t",row.names = FALSE, quote=FALSE)

# cyc_summ<-data.frame(condition=c("NA","FA","FT"),
#                      cycling_n=c(25,12,5),
#                      not_cycling_n=c(262-25,262-12,262-5))%>%
#   mutate(cycling_perc=(cycling_n/262)*100,
#          not_cycling_perc=(not_cycling_n/262)*100)
# write.table(cyc_summ,"cyclic_analysis/wol2_dna_pipeline_genus-prop/SFR23_0828_MG_cycling_summary_table.txt",sep = "\t",row.names = FALSE, quote=FALSE)

mtx_chi<-cyc_summ%>%dplyr::select(1:3)%>%
  filter(condition!="NA")%>%
  column_to_rownames("condition")%>%as.matrix()

chisq.test(mtx_chi,simulate.p.value=TRUE, B=2000) 

#OGU -rclr
#X-squared = 17.506, df = NA, p-value = 0.0004998
#X-squared = 3.0728, df = NA, p-value = 0.1264 no FA vs FT diff
#X-squared = 15.98, df = NA, p-value = 0.0004998 NA vs FT diff
#X-squared = 5.7187, df = NA, p-value = 0.02199 NA vs FA diff

#OGU -prop
#X-squared = 34.316, df = NA, p-value = 0.0004998
#X-squared = 0.048513, df = NA, p-value = 1 no FA vs FT diff
#X-squared = 21.549, df = NA, p-value = 0.0004998 NA vs FT diff
#X-squared = 19.906, df = NA, p-value = 0.0004998 NA vs FA diff

#species -rclr
#X-squared = 17.876, df = NA, p-value = 0.0004998
#X-squared = 2.7361, df = NA, p-value = 0.1569 no FA vs FT diff
#X-squared = 15.935, df = NA, p-value = 0.0004998 NA vs FT diff
#X-squared = 6.2201, df = NA, p-value = 0.02299 NA vs FA diff

#species -prop
#X-squared = 40.038, df = NA, p-value = 0.0004998
#X-squared = 8.1552, df = NA, p-value = 0.01099 FA vs FT diff
#X-squared = 34.419, df = NA, p-value = 0.0004998 NA vs FT diff
#X-squared = 12.831, df = NA, p-value = 0.0004998 NA vs FA diff

#genus -rclr
#X-squared = 15.315, df = NA, p-value = 0.0004998
#X-squared = 2.642, df = NA, p-value = 0.1789 no FA vs FT diff
#X-squared = 13.662, df = NA, p-value = 0.0004998 NA vs FT diff 
#X-squared = 4.7927, df = NA, p-value = 0.04598 NA vs FA diff

#genus -prop
#X-squared = 15.545, df = NA, p-value = 0.0004998
#X-squared = 2.979, df = NA, p-value = 0.1349 no FA vs FT diff
#X-squared = 14.143, df = NA, p-value = 0.0009995 NA vs FT diff 
#X-squared = 4.9146, df = NA, p-value = 0.03148 NA vs FA diff

cyc_summ_plot<-cyc_summ%>%gather(cycling_grp, value, -condition)%>%
  filter(cycling_grp=="cycling_perc"|cycling_grp=="not_cycling_perc")%>%
  mutate(condition=factor(condition,levels=c("NA","FA","FT")))%>%
  filter(cycling_grp=="cycling_perc")

p<-ggplot(data=cyc_summ_plot, aes(x=condition, y=value, fill=condition)) +
  geom_bar(stat="identity") + theme_classic() +
  theme(legend.position = "top")+
  scale_y_continuous(expand=c(0,0), limits=c(0,10)) +
  labs(title="OGU",y="genes (%)") +
  #labs(title="Species",y="genes (%)") +
  #labs(title="Genus",y="genes (%)") +
  scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))

#ggsave("cyclic_analysis/wol2_dna_pipeline_ogu-rclr/SFR23_0823_MG_cycling_summ_barplot.pdf", plot=p,height=3, width=2.5)
ggsave("cyclic_analysis/wol2_dna_pipeline_ogu-prop/SFR23_0828_MG_cycling_summ_barplot.pdf", plot=p,height=3, width=2.5)
#ggsave("cyclic_analysis/wol2_dna_pipeline_species-rclr/SFR23_0823_MG_cycling_summ_barplot.pdf", plot=p,height=3, width=2.5)
#ggsave("cyclic_analysis/wol2_dna_pipeline_species-prop/SFR23_0828_MG_cycling_summ_barplot.pdf", plot=p,height=3, width=2.5)
#ggsave("cyclic_analysis/wol2_dna_pipeline_genus-rclr/SFR23_0823_MG_cycling_summ_barplot.pdf", plot=p,height=3, width=2.5)
#ggsave("cyclic_analysis/wol2_dna_pipeline_genus-prop/SFR23_0828_MG_cycling_summ_barplot.pdf", plot=p,height=3, width=2.5)

#make venn diagram showing overlap of whats cycling

list_venn <- list(NA_ = sigNA$FeatureID,
                  FA = sigFA$FeatureID,
                  FT = sigFT$FeatureID)

draw.venn(sigNA$FeatureID,sigFA$FeatureID,sigFT$FeatureID,
          title="",subtitle="",
          xtitle="NA",ytitle="FA",ztitle="FT",
          xt_s=3,yt_s=3,zt_s=3, nr_s=3,
          x_c="#0072B2",y_c="#D55E00",z_c="#009E73",
          output="pdf",
          #filename="cyclic_analysis/wol2_dna_pipeline_ogu-rclr/SFR23_0823_MG_venn_overlap.pdf")
          filename="cyclic_analysis/wol2_dna_pipeline_ogu-prop/SFR23_0828_MG_venn_overlap.pdf")
          #filename="cyclic_analysis/wol2_dna_pipeline_species-rclr/SFR23_0823_MG_venn_overlap.pdf")
          #filename="cyclic_analysis/wol2_dna_pipeline_species-prop/SFR23_0828_MG_venn_overlap.pdf")
          #filename="cyclic_analysis/wol2_dna_pipeline_genus-rclr/SFR23_0823_MG_venn_overlap.pdf")
          #filename="cyclic_analysis/wol2_dna_pipeline_genus-prop/SFR23_0828_MG_venn_overlap.pdf")

##########################################################
#create heatmap describing shared patterns

ItemsList <- venn(list_venn, show.plot = FALSE)
all<-attributes(ItemsList)$intersections

shared_cyc<-plyr::ldply (all, data.frame)%>%
  dplyr::rename(pattern=".id", FeatureID="X..i..")

tworhythmic<-shared_cyc%>%filter(pattern %in% c("NA_:FA:FT","FA:FT","NA_:FA","NA_:FT"))

md<-fread("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metagenomic/metaG_metadata_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))
go_annot<-fread("~/scratch/TRF_multiomics/metatranscript/woltka2_results/go_name.txt")%>%
  dplyr::rename(FeatureID=GO_Term)

pfamTPM_all<-fread("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metagenomic/woltka2_results/filtered_metaG/pfam_TPM/pfam_clean_noNT_TPM_normRPOB.txt")%>%
  filter(FeatureID %in% tworhythmic$FeatureID)%>%
  gather(sample_name,TPM_counts,-FeatureID) %>%
  left_join(.,pfam_annot, by ="FeatureID") %>%
  mutate(label_name=paste(FeatureID, Name, sep=" "))%>%
  left_join(.,md,by="sample_name") %>%
  left_join(.,tworhythmic,by="FeatureID")%>%
  group_by(label_name,condition,zt_time, pattern)%>%dplyr::summarise(mn_TPM=mean(TPM_counts))%>%
  group_by(label_name,condition)%>%mutate(Zscore=(mn_TPM - mean(mn_TPM))/sd(mn_TPM))%>%
   mutate(zt_time=factor(zt_time, levels = c("1", "5","9","13","17","21")),
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

ggsave("cyclic_analysis/wol2_dna_pipeline_rpob/SFR23_0619_heatmap_byphase_RRN_zbycond.pdf", plot=plt,height=3, width=9)

ItemsList <- venn(list_venn, show.plot = FALSE)
all<-attributes(ItemsList)$intersections

#pfam

FTNA_shared_cyc<-pfam_annot%>%filter(FeatureID %in% all$`NA_:FT`) #3
FTFA_shared_cyc<-pfam_annot%>%filter(FeatureID %in% all$`FA:FT`) #1
NA_cyc<-pfam_annot%>%filter(FeatureID %in% all$`NA_`) #302
FA_cyc<-pfam_annot%>%filter(FeatureID %in% all$`FA`) #38
FT_cyc<-pfam_annot%>%filter(FeatureID %in% all$`FT`) #41

#plot specific examples over ZT time
stderror <- function(x) sd(x)/sqrt(length(x))

pfamZT<-fread("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metagenomic/woltka2_results/filtered_metaG/pfam_TPM/pfam_clean_noNT_TPM_normRPOB.txt")%>%
  filter(FeatureID %in% FTNA_shared_cyc$FeatureID)%>% 
  gather(sample_name,TPM_counts,-FeatureID) %>%
  left_join(.,pfam_annot, by ="FeatureID") %>%
  mutate(label_name=paste(FeatureID, Name, sep=" "))%>%
  left_join(.,md,by="sample_name") %>%
  group_by(FeatureID,condition,zt_time)%>%summarise(mn_TPM=mean(TPM_counts),sem=stderror(TPM_counts))%>%
  mutate(condition=factor(condition, levels = c("NA","FA","FT")))

ZT_dist <- function(data,FeatureID) {
  
  data%>%
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
         title=FeatureID)+ theme(plot.title = element_text(face = "bold"))
  ggsave(paste("cyclic_analysis/wol2_dna_pipeline_rpob/FTNA_shared_cyc_zt/",FeatureID,".pdf",sep=""), height=3, width=3)
}


pfamZT_nested <- pfamZT %>% 
  group_by(FeatureID) %>% 
  nest()

pfamZT_plots <- 
  pfamZT_nested %>% 
  mutate(plot = map2(data, FeatureID,  ~ ZT_dist(.x,.y)))

#find what cycling is shared by FT FA
FT_1sharedcyc<-sigFT%>%filter(FeatureID %in% all$`FA:FT`)%>%
  dplyr::select(FeatureID,JTK_adjphase)%>%
  mutate(condition="FT")

FA_1sharedcyc<-sigFA%>%filter(FeatureID %in% all$`FA:FT`)%>%
  dplyr::select(FeatureID,JTK_adjphase)%>%
  mutate(condition="FA")

comb_1sharedcyc<-rbind(FT_1sharedcyc,FA_1sharedcyc)%>%
  mutate(condition=factor(condition,levels=c("FA","FT")))

write.table(comb_1sharedcyc,"cyclic_analysis/wol2_dna_pipeline_rpob/SFR23_0619_MG_1cycFTFA.txt",sep = "\t",row.names = FALSE, quote=FALSE)


ggplot(comb_1sharedcyc, aes(x=JTK_adjphase, color=condition,fill=condition)) +
  geom_density(alpha=0.3)+ theme_pubr()+
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_color_manual(values=c("#D55E00","#009E73"))+
  scale_fill_manual(values=c("#D55E00","#009E73")) #cant run with just one feature


#find what cycling is shared by FT NA
FT_3sharedcyc<-sigFT%>%filter(FeatureID %in% all$`NA_:FT`)%>%
  dplyr::select(FeatureID,JTK_adjphase)%>%
  mutate(condition="FT")

NA_3sharedcyc<-sigNA%>%filter(FeatureID %in% all$`NA_:FT`)%>%
  dplyr::select(FeatureID,JTK_adjphase)%>%
  mutate(condition="NA")

comb_3sharedcyc<-rbind(FT_3sharedcyc,NA_3sharedcyc)%>%
  mutate(condition=factor(condition,levels=c("NA","FT")))

write.table(comb_3sharedcyc,"cyclic_analysis/wol2_dna_pipeline_rpob/SFR23_0619_MG_3cycFTNA.txt",sep = "\t",row.names = FALSE, quote=FALSE)


ggplot(comb_3sharedcyc, aes(x=JTK_adjphase, color=condition,fill=condition)) +
  geom_density(alpha=0.3)+ theme_pubr()+
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_color_manual(values=c("#0072B2","#009E73"))+
  scale_fill_manual(values=c("#0072B2","#009E73"))

ggsave("cyclic_analysis/wol2_dna_pipeline_rpob/SFR23_0619_MG_densityPhasenoFA_3cycFTNA.pdf", width = 4.5, height = 4)


#find what cycling is shared by FA NA
FA_3sharedcyc<-sigFA%>%filter(FeatureID %in% all$`NA_:FA`)%>%
  dplyr::select(FeatureID,JTK_adjphase)%>%
  mutate(condition="FA")

NA_3sharedcyc<-sigNA%>%filter(FeatureID %in% all$`NA_:FA`)%>%
  dplyr::select(FeatureID,JTK_adjphase)%>%
  mutate(condition="NA")

comb_3sharedcyc<-rbind(FA_3sharedcyc,NA_3sharedcyc)%>%
  mutate(condition=factor(condition,levels=c("NA","FA")))

write.table(comb_3sharedcyc,"cyclic_analysis/wol2_dna_pipeline_rpob/SFR23_0619_MG_3cycFANA.txt",sep = "\t",row.names = FALSE, quote=FALSE)


ggplot(comb_3sharedcyc, aes(x=JTK_adjphase, color=condition,fill=condition)) +
  geom_density(alpha=0.3)+ theme_pubr()+
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_color_manual(values=c("#0072B2","#D55E00"))+
  scale_fill_manual(values=c("#0072B2","#D55E00")) 

ggsave("cyclic_analysis/wol2_dna_pipeline_rpob/SFR23_0619_MG_densityPhasenoFA_3cycFANA.pdf", width = 4.5, height = 4)


#find what cycling is shared by FT FA NA -there is nothing cycling in all 3
FT_0sharedcyc<-sigFT%>%filter(FeatureID %in% all$`NA_:FA:FT`)%>%
  dplyr::select(FeatureID,JTK_adjphase)%>% 
  mutate(condition="FT")

FA_0sharedcyc<-sigFA%>%filter(FeatureID %in% all$`NA_:FA:FT`)%>%
  dplyr::select(FeatureID,JTK_adjphase)%>%
  mutate(condition="FA")

NA_0sharedcyc<-sigNA%>%filter(FeatureID %in% all$`NA_:FA:FT`)%>%
  dplyr::select(FeatureID,JTK_adjphase)%>%
  mutate(condition="NA")

#look at what type of cycling hits were just in FT

#pfam 
FT_cyc<-pfam_annot%>%filter(FeatureID %in% all$FT) #41

gonames<-fread("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_results/go_name.txt")
pfamGOp<-read.table("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_results/pfam-to-go-process.map",header = FALSE, sep = "\t",
                      col.names = paste0("V",seq_len(4)), fill = TRUE)%>%
  gather(column,GO_Term,-V1)%>%
  dplyr::select(1,3)%>%
  dplyr::rename(FeatureID=V1)%>%
  filter(FeatureID %in% FT_cyc$FeatureID)%>%
  left_join(.,gonames,by="GO_Term")%>%
  filter(!is.na(name))

FTcyc_summ<-pfamGOp%>%
  group_by(name)%>%summarise(n=n())%>%
  mutate(condition="FT")%>%
  arrange(n)

FTcyc_summ$name <- factor(FTcyc_summ$name,levels = FTcyc_summ$name )

ggplot(data=FTcyc_summ, aes(x=name, y=n, fill=condition)) +
  geom_bar(stat="identity", position=position_dodge())  +coord_flip() +
  scale_fill_manual(values=c("#009E73")) +theme_pubr() +
  scale_y_continuous(expand=c(0,0), limits=c(0,10))

ggsave("cyclic_analysis/wol2_dna_pipeline_rpob/SFR23_0619_cycFT_GOterms.pdf",height=3.5, width=8)

FA_cyc<-pfam_annot%>%filter(FeatureID %in% all$FA) #38

pfamGOp2<-read.table("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_results/pfam-to-go-process.map",header = FALSE, sep = "\t",
                    col.names = paste0("V",seq_len(4)), fill = TRUE)%>%
  gather(column,GO_Term,-V1)%>%
  dplyr::select(1,3)%>%
  dplyr::rename(FeatureID=V1)%>%
  filter(FeatureID %in% FA_cyc$FeatureID)%>%
  left_join(.,gonames,by="GO_Term")%>%
  filter(!is.na(name))

FAcyc_summ<-pfamGOp2%>%
  group_by(name)%>%summarise(n=n())%>%
  mutate(condition="FA")%>%
  arrange(n)

FAcyc_summ$name <- factor(FAcyc_summ$name,levels = FAcyc_summ$name )

ggplot(data=FAcyc_summ, aes(x=name, y=n, fill=condition)) +
  geom_bar(stat="identity", position=position_dodge()) + coord_flip() +
  scale_fill_manual(values=c("#D55E00")) +theme_pubr() +
  scale_y_continuous(expand=c(0,0), limits=c(0,10))

ggsave("cyclic_analysis/wol2_dna_pipeline_rpob/SFR23_0619_cycFA_GOterms.pdf",height=3, width=6)

NA_cyc<-pfam_annot%>%filter(FeatureID %in% all$NA_) #302

pfamGOp3<-read.table("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_results/pfam-to-go-process.map",header = FALSE, sep = "\t",
                     col.names = paste0("V",seq_len(4)), fill = TRUE)%>%
  gather(column,GO_Term,-V1)%>%
  dplyr::select(1,3)%>%
  dplyr::rename(FeatureID=V1)%>%
  filter(FeatureID %in% NA_cyc$FeatureID)%>%
  left_join(.,gonames,by="GO_Term")%>%
  filter(!is.na(name))

NAcyc_summ<-pfamGOp3%>%
  group_by(name)%>%summarise(n=n())%>%
  mutate(condition="NA")%>%
  filter(n>1)%>%
   arrange(n)

NAcyc_summ$name <- factor(NAcyc_summ$name,levels = NAcyc_summ$name )

ggplot(data=NAcyc_summ, aes(x=name, y=n, fill=condition)) +
  geom_bar(stat="identity", position=position_dodge()) + coord_flip() +
  scale_fill_manual(values=c("#0072B2")) +theme_pubr() +
  scale_y_continuous(expand=c(0,0), limits=c(0,20))

ggsave("cyclic_analysis/wol2_dna_pipeline_rpob/SFR23_0619_cycNA_GOterms.pdf",height=3.5, width=6)


FTFA_cyc<-pfam_annot%>%filter(FeatureID %in% all$`FA:FT`) #1

pfamGOp4<-read.table("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_results/pfam-to-go-process.map",header = FALSE, sep = "\t",
                     col.names = paste0("V",seq_len(4)), fill = TRUE)%>%
  gather(column,GO_Term,-V1)%>%
  dplyr::select(1,3)%>%
  dplyr::rename(FeatureID=V1)%>%
  filter(FeatureID %in% FTFA_cyc$FeatureID)%>%
  left_join(.,gonames,by="GO_Term")%>%
  filter(!is.na(name)) #no hits

FTNA_cyc<-pfam_annot%>%filter(FeatureID %in% all$`NA_:FT`) #3

pfamGOp5<-read.table("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_results/pfam-to-go-process.map",header = FALSE, sep = "\t",
                     col.names = paste0("V",seq_len(4)), fill = TRUE)%>%
  gather(column,GO_Term,-V1)%>%
  dplyr::select(1,3)%>%
  dplyr::rename(FeatureID=V1)%>%
  filter(FeatureID %in% FTNA_cyc$FeatureID)%>%
  left_join(.,gonames,by="GO_Term")%>%
  filter(!is.na(name)) #no hits
