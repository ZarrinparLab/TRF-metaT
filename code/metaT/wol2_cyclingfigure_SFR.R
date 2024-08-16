setwd("~/scratch/TRF_multiomics/metatranscript/woltka2_m_results/pfam/")

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
#MT

#annotations
pfam_annot<-fread("pfam_annotationkey.csv")

#read in FT metacycle hits
FT_metacyc<-fread("cyclic_analysis/wol2_rna_pipeline/FT_metacycle/meta2d_filtered_rna_FT.txt")%>%
  dplyr::rename(FeatureID=CycID) %>%
  left_join(.,pfam_annot, by="FeatureID")%>%
  arrange(JTK_adjphase)%>%
  mutate(label_name=paste(FeatureID, Name, sep=" "))

#read in FA metacycle hits
FA_metacyc<-fread("cyclic_analysis/wol2_rna_pipeline/FA_metacycle/meta2d_filtered_rna_FA.txt")%>%
  dplyr::rename(FeatureID=CycID) %>%
  left_join(.,pfam_annot, by="FeatureID")%>%
  arrange(JTK_adjphase)%>%
  mutate(label_name=paste(FeatureID, Name, sep=" "))

#read in NA metacycle hits
NA_metacyc<-fread("cyclic_analysis/wol2_rna_pipeline/NA_metacycle/meta2d_filtered_rna_NA.txt")%>%
  dplyr::rename(FeatureID=CycID) %>%
  left_join(.,pfam_annot, by="FeatureID")%>%
  arrange(JTK_adjphase)%>%
  mutate(label_name=paste(FeatureID, Name, sep=" "))

#rpob is not cycling in any of the conditions

##########################################################
#MT

#annotations
pfam_annot<-fread("pfam_annotationkey.csv")

#read in FT metacycle hits
FT_metacyc<-fread("cyclic_analysis/wol2_rna_pipeline_144k_rpob/FT_metacycle/meta2d_filtered_rna_FT.txt")%>%
  dplyr::rename(FeatureID=CycID) %>%
  left_join(.,pfam_annot, by="FeatureID")%>%
  arrange(JTK_adjphase)%>%
  mutate(label_name=paste(FeatureID, Name, sep=" "))

#read in FA metacycle hits
FA_metacyc<-fread("cyclic_analysis/wol2_rna_pipeline_144k_rpob/FA_metacycle/meta2d_filtered_rna_FA.txt")%>%
  dplyr::rename(FeatureID=CycID) %>%
  left_join(.,pfam_annot, by="FeatureID")%>%
  arrange(JTK_adjphase)%>%
  mutate(label_name=paste(FeatureID, Name, sep=" "))

#read in NA metacycle hits
NA_metacyc<-fread("cyclic_analysis/wol2_rna_pipeline_144k_rpob/NA_metacycle/meta2d_filtered_rna_NA.txt")%>%
  dplyr::rename(FeatureID=CycID) %>%
  left_join(.,pfam_annot, by="FeatureID")%>%
  arrange(JTK_adjphase)%>%
  mutate(label_name=paste(FeatureID, Name, sep=" "))

#make a bar plot summarizing the number of cycling and non-cycling hits (removing DUF)

sigFT<-FT_metacyc%>%filter(JTK_pvalue<0.05,!grepl("DUF",Name)) #122
sigFA<-FA_metacyc%>%filter(JTK_pvalue<0.05,!grepl("DUF",Name)) #56
sigNA<-NA_metacyc%>%filter(JTK_pvalue<0.05,!grepl("DUF",Name)) #292
noduf<-FT_metacyc%>%filter(!grepl("DUF",Name)) #4501

cyc_summ<-data.frame(condition=c("NA","FA","FT"),
                     cycling_n=c(292,56,122),
                     not_cycling_n=c(4501-292,4501-56,4501-122))%>%
  mutate(cycling_perc=(cycling_n/4501)*100,
         not_cycling_perc=(not_cycling_n/4501)*100)
write.table(cyc_summ,"cyclic_analysis/wol2_rna_pipeline_144k_rpob/SFR23_0620_MT_cycling_summary_table.txt",sep = "\t",row.names = FALSE, quote=FALSE)

mtx_chi<-cyc_summ%>%dplyr::select(1:3)%>%
  #filter(condition!="NA")%>%
  #filter(condition!="FT")%>%
  filter(condition!="FA")%>%
  column_to_rownames("condition")%>%as.matrix()

chisq.test(mtx_chi,simulate.p.value=TRUE, B=2000) #X-squared = 196.08, df = NA, p-value = 0.0004998
#X-squared = 24.966, df = NA, p-value = 0.0004998 just FA FT
#X-squared = 166.48, df = NA, p-value = 0.0004998 FA NA
#X-squared = 73.172, df = NA, p-value = 0.0004998 FT NA

mtx_chi<-cyc_summ%>%dplyr::select(1:3)%>%
  filter(condition!="FA")%>%
  column_to_rownames("condition")%>%as.matrix()

chisq.test(mtx_chi,simulate.p.value=TRUE, B=2000) #X-squared = 73.172, df = NA, p-value = 0.0004998 just NA FT

mtx_chi<-cyc_summ%>%dplyr::select(1:3)%>%
  filter(condition!="FT")%>%
  column_to_rownames("condition")%>%as.matrix()

chisq.test(mtx_chi,simulate.p.value=TRUE, B=2000) #X-squared = 166.48, df = NA, p-value = 0.0004998

cyc_summ_plot<-cyc_summ%>%gather(cycling_grp, value, -condition)%>%
  filter(cycling_grp=="cycling_perc"|cycling_grp=="not_cycling_perc")%>%
  mutate(condition=factor(condition,levels=c("NA","FA","FT")))%>%
  filter(cycling_grp=="cycling_perc")

p<-ggplot(data=cyc_summ_plot, aes(x=condition, y=value, fill=condition)) +
  geom_bar(stat="identity") + theme_classic() +
  theme(legend.position = "top")+
  scale_y_continuous(expand=c(0,0), limits=c(0,10)) +
  labs(y="transcripts (%)") +
  scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))

ggsave("cyclic_analysis/wol2_rna_pipeline_144k_rpob/SFR23_0620_MT_cycling_summ_barplot.pdf", plot=p,height=3, width=2.5)

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
          filename="cyclic_analysis/wol2_rna_pipeline_144k_rpob/SFR23_0620_MT_venn_overlap.pdf")
##########################################################
#MT

#annotations
pfam_annot<-fread("pfam_annotationkey.csv")

#read in FT metacycle hits
FT_metacyc<-fread("cyclic_analysis/wol2_rna_pipeline_rpob/FT_metacycle/meta2d_filtered_rna_FT.txt")%>%
  dplyr::rename(FeatureID=CycID) %>%
  left_join(.,pfam_annot, by="FeatureID")%>%
  arrange(JTK_adjphase)%>%
  mutate(label_name=paste(FeatureID, Name, sep=" "))

#read in FA metacycle hits
FA_metacyc<-fread("cyclic_analysis/wol2_rna_pipeline_rpob/FA_metacycle/meta2d_filtered_rna_FA.txt")%>%
  dplyr::rename(FeatureID=CycID) %>%
  left_join(.,pfam_annot, by="FeatureID")%>%
  arrange(JTK_adjphase)%>%
  mutate(label_name=paste(FeatureID, Name, sep=" "))

#read in NA metacycle hits
NA_metacyc<-fread("cyclic_analysis/wol2_rna_pipeline_rpob/NA_metacycle/meta2d_filtered_rna_NA.txt")%>%
  dplyr::rename(FeatureID=CycID) %>%
  left_join(.,pfam_annot, by="FeatureID")%>%
  arrange(JTK_adjphase)%>%
  mutate(label_name=paste(FeatureID, Name, sep=" "))

#make a bar plot summarizing the number of cycling and non-cycling hits (removing DUF)

sigFT<-FT_metacyc%>%filter(JTK_pvalue<0.05,!grepl("DUF",Name)) #662
sigFA<-FA_metacyc%>%filter(JTK_pvalue<0.05,!grepl("DUF",Name)) #117
sigNA<-NA_metacyc%>%filter(JTK_pvalue<0.05,!grepl("DUF",Name)) #571
noduf<-FT_metacyc%>%filter(!grepl("DUF",Name)) #7113

cyc_summ<-data.frame(condition=c("NA","FA","FT"),
                     cycling_n=c(571,117,662),
                     not_cycling_n=c(7113-571,7113-117,7113-662))%>%
  mutate(cycling_perc=(cycling_n/7113)*100,
         not_cycling_perc=(not_cycling_n/7113)*100)
write.table(cyc_summ,"cyclic_analysis/wol2_rna_pipeline_rpob/SFR23_0619_MT_cycling_summary_table.txt",sep = "\t",row.names = FALSE, quote=FALSE)

mtx_chi<-cyc_summ%>%dplyr::select(1:3)%>%
  #filter(condition!="NA")%>%
  #filter(condition!="FT")%>%
  filter(condition!="FA")%>%
  column_to_rownames("condition")%>%as.matrix()

chisq.test(mtx_chi,simulate.p.value=TRUE, B=2000) #X-squared = 404.42, df = NA, p-value = 0.0004998
#X-squared = 403.38, df = NA, p-value = 0.0004998 just FA FT
#X-squared = 314.81, df = NA, p-value = 0.0004998 just FA NA
#X-squared = 7.3535, df = NA, p-value = 0.008996 just NA FT

cyc_summ_plot<-cyc_summ%>%gather(cycling_grp, value, -condition)%>%
  filter(cycling_grp=="cycling_perc"|cycling_grp=="not_cycling_perc")%>%
  mutate(condition=factor(condition,levels=c("NA","FA","FT")))%>%
  filter(cycling_grp=="cycling_perc")

p<-ggplot(data=cyc_summ_plot, aes(x=condition, y=value, fill=condition)) +
  geom_bar(stat="identity") + theme_classic() +
  theme(legend.position = "top")+
  scale_y_continuous(expand=c(0,0), limits=c(0,10)) +
  labs(y="transcripts (%)") +
  scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))

ggsave("cyclic_analysis/wol2_rna_pipeline_rpob/SFR23_0619_MT_cycling_summ_barplot.pdf", plot=p,height=3, width=2.5)

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
          filename="cyclic_analysis/wol2_rna_pipeline_rpob/SFR23_0619_MT_venn_overlap.pdf")



#create heatmap describing shared patterns

ItemsList <- venn(list_venn, show.plot = FALSE)
all<-attributes(ItemsList)$intersections

shared_cyc<-plyr::ldply (all, data.frame)%>%
  dplyr::rename(pattern=".id", FeatureID="X..i..")

tworhythmic<-shared_cyc%>%filter(pattern %in% c("NA_:FA:FT","FA:FT","NA_:FA","NA_:FT"))

md<-fread("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_m_results/metaT_metadata_ztcat_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))
go_annot<-fread("~/scratch/TRF_multiomics/metatranscript/woltka2_results/go_name.txt")%>%
  dplyr::rename(FeatureID=GO_Term)

pfamTPM_all<-fread("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_results/pfam_TPM/pfam_clean_noNT_normRPOB.txt")%>%
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

ggsave("cyclic_analysis/wol2_rna_pipeline_rpob/SFR23_0619_heatmap_byphase_RRN_zbycond.pdf", plot=plt,height=4, width=9)

ItemsList <- venn(list_venn, show.plot = FALSE)
all<-attributes(ItemsList)$intersections

#pfam

FTNA_shared_cyc<-pfam_annot%>%filter(FeatureID %in% all$`NA_:FT`) #86
FTFA_shared_cyc<-pfam_annot%>%filter(FeatureID %in% all$`FA:FT`) #38
NA_cyc<-pfam_annot%>%filter(FeatureID %in% all$`NA_`) #470
FA_cyc<-pfam_annot%>%filter(FeatureID %in% all$`FA`) #64
FT_cyc<-pfam_annot%>%filter(FeatureID %in% all$`FT`) #534

#plot specific examples over ZT time
stderror <- function(x) sd(x)/sqrt(length(x))

pfamZT<-fread("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_m_results/pfam/pfam-TPM_clean_noNT_normRPOB.txt")%>%
  filter(FeatureID %in% FT_cyc$FeatureID)%>% 
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
  ggsave(paste("cyclic_analysis/wol2_rna_pipeline_rpob/FT_cyc_zt/",FeatureID,".pdf",sep=""), height=3, width=3)
}


pfamZT_nested <- pfamZT %>% 
  group_by(FeatureID) %>% 
  nest()

pfamZT_plots <- 
  pfamZT_nested %>% 
  mutate(plot = map2(data, FeatureID,  ~ ZT_dist(.x,.y)))

#find what cycling is shared by FT FA
FT_38sharedcyc<-sigFT%>%filter(FeatureID %in% all$`FA:FT`)%>%
  dplyr::select(FeatureID,JTK_adjphase)%>%
  mutate(condition="FT")

FA_38sharedcyc<-sigFA%>%filter(FeatureID %in% all$`FA:FT`)%>%
  dplyr::select(FeatureID,JTK_adjphase)%>%
  mutate(condition="FA")

comb_38sharedcyc<-rbind(FT_38sharedcyc,FA_38sharedcyc)%>%
  mutate(condition=factor(condition,levels=c("FA","FT")))

write.table(comb_38sharedcyc,"cyclic_analysis/wol2_rna_pipeline_rpob/SFR23_0619_MT_38cycFTFA.txt",sep = "\t",row.names = FALSE, quote=FALSE)


ggplot(comb_38sharedcyc, aes(x=JTK_adjphase, color=condition,fill=condition)) +
  geom_density(alpha=0.3)+ theme_pubr()+
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_color_manual(values=c("#D55E00","#009E73"))+
  scale_fill_manual(values=c("#D55E00","#009E73")) 

ggsave("cyclic_analysis/wol2_rna_pipeline_rpob/SFR23_0619_MT_densityPhasenoFA_38cycFTFA.pdf", width = 4.5, height = 4)

#find what cycling is shared by FT NA
FT_86sharedcyc<-sigFT%>%filter(FeatureID %in% all$`NA_:FT`)%>%
  dplyr::select(FeatureID,JTK_adjphase)%>%
  mutate(condition="FT")

NA_86sharedcyc<-sigNA%>%filter(FeatureID %in% all$`NA_:FT`)%>%
  dplyr::select(FeatureID,JTK_adjphase)%>%
  mutate(condition="NA")

comb_86sharedcyc<-rbind(FT_86sharedcyc,NA_86sharedcyc)%>%
  mutate(condition=factor(condition,levels=c("NA","FT")))

write.table(comb_86sharedcyc,"cyclic_analysis/wol2_rna_pipeline_rpob/SFR23_0619_MT_86cycFTNA.txt",sep = "\t",row.names = FALSE, quote=FALSE)


ggplot(comb_86sharedcyc, aes(x=JTK_adjphase, color=condition,fill=condition)) +
  geom_density(alpha=0.3)+ theme_pubr()+
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_color_manual(values=c("#0072B2","#009E73"))+
  scale_fill_manual(values=c("#0072B2","#009E73"))

ggsave("cyclic_analysis/wol2_rna_pipeline_rpob/SFR23_0619_MT_densityPhasenoFA_86cycFTNA.pdf", width = 4.5, height = 4)


#find what cycling is shared by FA NA
FA_11sharedcyc<-sigFA%>%filter(FeatureID %in% all$`NA_:FA`)%>%
  dplyr::select(FeatureID,JTK_adjphase)%>%
  mutate(condition="FA")

NA_11sharedcyc<-sigNA%>%filter(FeatureID %in% all$`NA_:FA`)%>%
  dplyr::select(FeatureID,JTK_adjphase)%>%
  mutate(condition="NA")

comb_11sharedcyc<-rbind(FA_11sharedcyc,NA_11sharedcyc)%>%
  mutate(condition=factor(condition,levels=c("NA","FA")))

write.table(comb_11sharedcyc,"cyclic_analysis/wol2_rna_pipeline_rpob/SFR23_0619_MT_11cycFANA.txt",sep = "\t",row.names = FALSE, quote=FALSE)


ggplot(comb_11sharedcyc, aes(x=JTK_adjphase, color=condition,fill=condition)) +
  geom_density(alpha=0.3)+ theme_pubr()+
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_color_manual(values=c("#0072B2","#D55E00"))+
  scale_fill_manual(values=c("#0072B2","#D55E00")) 

ggsave("cyclic_analysis/wol2_rna_pipeline_rpob/SFR23_0619_MT_densityPhasenoFA_11cycFANA.pdf", width = 4.5, height = 4)


#find what cycling is shared by FT FA NA
FT_4sharedcyc<-sigFT%>%filter(FeatureID %in% all$`NA_:FA:FT`)%>%
  dplyr::select(FeatureID,JTK_adjphase)%>%
  mutate(condition="FT")

FA_4sharedcyc<-sigFA%>%filter(FeatureID %in% all$`NA_:FA:FT`)%>%
  dplyr::select(FeatureID,JTK_adjphase)%>%
  mutate(condition="FA")

NA_4sharedcyc<-sigNA%>%filter(FeatureID %in% all$`NA_:FA:FT`)%>%
  dplyr::select(FeatureID,JTK_adjphase)%>%
  mutate(condition="NA")

comb_4sharedcyc<-rbind(FT_4sharedcyc,FA_4sharedcyc,NA_4sharedcyc)%>%
  mutate(condition=factor(condition,levels=c("NA","FA","FT")))

write.table(comb_4sharedcyc,"cyclic_analysis/wol2_rna_pipeline_rpob/SFR23_0619_MT_4cycFTFANA.txt",sep = "\t",row.names = FALSE, quote=FALSE)


ggplot(comb_4sharedcyc, aes(x=JTK_adjphase, color=condition,fill=condition)) +
  geom_density(alpha=0.3)+ theme_pubr()+
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+
  scale_fill_manual(values=c("#0072B2","#D55E00","#009E73")) 

ggsave("cyclic_analysis/wol2_rna_pipeline_rpob/SFR23_0619_MT_densityPhasenoFA_3cycFTFANA.pdf", width = 4.5, height = 4)

#look at what type of cycling hits were just in FT

#pfam 
FT_cyc<-pfam_annot%>%filter(FeatureID %in% all$FT) #534

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
  filter(n>2)%>%
  arrange(n)

FTcyc_summ$name <- factor(FTcyc_summ$name,levels = FTcyc_summ$name )

ggplot(data=FTcyc_summ, aes(x=name, y=n, fill=condition)) +
  geom_bar(stat="identity", position=position_dodge())  +coord_flip() +
  scale_fill_manual(values=c("#009E73")) +theme_pubr() +
  scale_y_continuous(expand=c(0,0), limits=c(0,20))

ggsave("cyclic_analysis/wol2_rna_pipeline_rpob/SFR23_0619_cycFT_GOterms.pdf",height=3.5, width=6)

FA_cyc<-pfam_annot%>%filter(FeatureID %in% all$FA) #64

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
  scale_y_continuous(expand=c(0,0), limits=c(0,20))

ggsave("cyclic_analysis/wol2_rna_pipeline_rpob/SFR23_0619_cycFA_GOterms.pdf",height=3, width=8)

NA_cyc<-pfam_annot%>%filter(FeatureID %in% all$NA_) #470

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

ggsave("cyclic_analysis/wol2_rna_pipeline_rpob/SFR23_0619_cycNA_GOterms.pdf",height=3.5, width=6)


FTFA_cyc<-pfam_annot%>%filter(FeatureID %in% all$`FA:FT`) #38

pfamGOp4<-read.table("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_results/pfam-to-go-process.map",header = FALSE, sep = "\t",
                     col.names = paste0("V",seq_len(4)), fill = TRUE)%>%
  gather(column,GO_Term,-V1)%>%
  dplyr::select(1,3)%>%
  dplyr::rename(FeatureID=V1)%>%
  filter(FeatureID %in% FTFA_cyc$FeatureID)%>%
  left_join(.,gonames,by="GO_Term")%>%
  filter(!is.na(name))

FTNA_cyc<-pfam_annot%>%filter(FeatureID %in% all$`NA_:FT`) #86

pfamGOp5<-read.table("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_results/pfam-to-go-process.map",header = FALSE, sep = "\t",
                     col.names = paste0("V",seq_len(4)), fill = TRUE)%>%
  gather(column,GO_Term,-V1)%>%
  dplyr::select(1,3)%>%
  dplyr::rename(FeatureID=V1)%>%
  filter(FeatureID %in% FTNA_cyc$FeatureID)%>%
  left_join(.,gonames,by="GO_Term")%>%
  filter(!is.na(name))

#create feature metadata file with Go terms for qurro

pfam2GO<-read.table("/mnt/zarrinpar/scratch/sfloresr/metatranscript/metatranscript/woltka2_results/pfam-to-go-process.map",header = FALSE, sep = "\t",
                       col.names = paste0("V",seq_len(4)), fill = TRUE)%>%
  gather(column,GO_Term,-V1)%>%
  dplyr::select(1,3)%>%
  dplyr::rename(FeatureID=V1)%>%
  left_join(.,gonames,by="GO_Term")%>%
  filter(!is.na(name))%>%
  group_by(FeatureID)%>%
  slice(1)%>%
  select(-GO_Term)%>%
  rename(GO_Term=name)
  #summarise(n=n())
  
feat_annot<-fread("/mnt/zarrinpar/scratch/sfloresr/metatranscript/metatranscript/woltka2_results/pfam/pfam_annotationkey.csv")%>%
  left_join(.,pfam2GO,by="FeatureID")

write.table(feat_annot,"/mnt/zarrinpar/scratch/sfloresr/metatranscript/metatranscript/woltka2_results/pfam/pfam_annotationkey_wGO.txt",sep = "\t",row.names = FALSE, quote=FALSE)

natlog_lipid<-fread("cycling_diversity/ordination_deicode_TRFmetaT_cyc_noNT/qurro_plot_data_cyc_lipid.tsv")%>%
  select(1:3)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  mutate(condition=factor(condition,levels=c("FT","FA","NA")))

p<-ggplot(natlog_lipid, aes(x=condition, y=Current_Natural_Log_Ratio,fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  theme_minimal()+ scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
  labs(x="condition",y="Natural Log Ratio",title="Lipid Metabolism")+
  theme(legend.position = "none")

pairwise.wilcox.test(natlog_lipid$Current_Natural_Log_Ratio, natlog_lipid$condition,
                     p.adjust.method="fdr")
# FT    FA   
# FA 0.031 -    
#   NA 0.002 0.095

ggsave("cycling_diversity/ordination_deicode_TRFmetaT_cyc_noNT/SFR22_1102_natlog_lipidvstransl_cyc.pdf", plot=p,height=3.5, width=3.5)

natlog_carbo<-fread("cycling_diversity/ordination_deicode_TRFmetaT_cyc_noNT/qurro_plot_data_cyc_carbo.tsv")%>%
  select(1:3)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  mutate(condition=factor(condition,levels=c("FT","FA","NA")))

p<-ggplot(natlog_carbo, aes(x=condition, y=Current_Natural_Log_Ratio,fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  theme_minimal()+ scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
  labs(x="condition",y="Natural Log Ratio",title="Carbohydrate Metabolic Process")+
  theme(legend.position = "none")

pairwise.wilcox.test(natlog_carbo$Current_Natural_Log_Ratio, natlog_carbo$condition,
                     p.adjust.method="fdr")
# FT   FA  
# FA 0.62 -   
#   NA 0.66 0.62

ggsave("cycling_diversity/ordination_deicode_TRFmetaT_cyc_noNT/SFR22_1102_natlog_carbovstransl_cyc.pdf", plot=p,height=3.5, width=3.5)

natlog_protlys<-fread("cycling_diversity/ordination_deicode_TRFmetaT_cyc_noNT/qurro_plot_data_cyc_proteolysis.tsv")%>%
  select(1:3)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  mutate(condition=factor(condition,levels=c("FT","FA","NA")))

p<-ggplot(natlog_protlys, aes(x=condition, y=Current_Natural_Log_Ratio,fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  theme_minimal()+ scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
  labs(x="condition",y="Natural Log Ratio",title="Proteolysis")+
  theme(legend.position = "none")

pairwise.wilcox.test(natlog_protlys$Current_Natural_Log_Ratio, natlog_protlys$condition,
                     p.adjust.method="fdr")
# FT      FA     
# FA 0.03628 -      
#   NA 0.00014 0.11337

ggsave("cycling_diversity/ordination_deicode_TRFmetaT_cyc_noNT/SFR22_1102_natlog_protlysvstransl_cyc.pdf", plot=p,height=3.5, width=3.5)