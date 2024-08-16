setwd("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_m_results/")

library(tidyverse)
library(data.table)
library("qiime2R")
library("Biostrings")
library(ggrepel)
library(ggpubr)
library(ggbreak)
library(scatterplot3d)

###########################################################
md<-fread("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/trf_cecal_metadata_noNT.tsv")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition),
         phase=ifelse(zt_time<13,"light","dark"),
         zt_time_cat=paste("ZT",zt_time, sep=""),
         cond_phase=paste(condition,phase,sep="_"),
         cond_zt=paste(condition,zt_time,sep="_"),
         fasted=ifelse(zt_time_cat %in% c("ZT5","ZT9","ZT13"), "fasted","fed"),
         cond_fasted=paste(condition,fasted,sep="_"))
#write.table(md, "metaT_metadata_ztcat_noNT.txt",sep = "\t",row.names = FALSE, quote=FALSE)
write.table(md, "metaT_metadata_ztcat_noNT_wfasting.txt",sep = "\t",row.names = FALSE, quote=FALSE)

samps_list<-c("#FeatureID","cFA01a_S39","cFA01b_S45","cFA05a_S40","cFA05b_S46","cFA09a_S41",
              "cFA09b_S47","cFA13a_S36","cFA13b_S42","cFA17a_S37","cFA17b_S43",
              "cFA21a_S38","cFA21b_S44","cFT01a_S51","cFT01b_S57","cFT05a_S52",
              "cFT05b_S58","cFT09a_S53","cFT09b_S59","cFT13a_S48","cFT13b_S54",
              "cFT17a_S49","cFT17b_S55","cFT21a_S50","cFT21b_S56","cNA01a_S4",
              "cNA01b_S10","cNA01c_S16","cNA05a_S5","cNA05b_S11","cNA05c_S17",
              "cNA09a_S6","cNA09b_S12","cNA09c_S18","cNA13a_S1","cNA13b_S7",
              "cNA13c_S13","cNA17a_S2","cNA17b_S8","cNA17c_S14","cNA21a_S3",
              "cNA21b_S9","cNA21c_S15")

#dat<-fread("pfam/pfam.tsv") %>%dplyr::select(all_of(samps_list))
#dat<-fread("pfam/pfam-TPM.tsv") %>%dplyr::select(all_of(samps_list))
dat<-fread("genome.tsv") %>%dplyr::select(all_of(samps_list))
names(dat) <- c("FeatureID","cFA01a","cFA01b","cFA05a","cFA05b","cFA09a","cFA09b","cFA13a","cFA13b","cFA17a","cFA17b",
                "cFA21a","cFA21b","cFT01a","cFT01b","cFT05a","cFT05b","cFT09a","cFT09b","cFT13a","cFT13b",
                "cFT17a","cFT17b","cFT21a","cFT21b","cNA01a","cNA01b","cNA01c","cNA05a","cNA05b","cNA05c",
                "cNA09a","cNA09b","cNA09c","cNA13a","cNA13b","cNA13c","cNA17a","cNA17b","cNA17c","cNA21a",
                "cNA21b","cNA21c")
#write.table(dat,"pfam/pfam_clean_noNT.tsv",sep = "\t",row.names = FALSE, quote=FALSE)
#write.table(dat,"pfam/pfam-TPM_clean_noNT.tsv",sep = "\t",row.names = FALSE, quote=FALSE)
write.table(dat,"genome_clean_noNT.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

annotations<-fread("pfam/pfam.tsv") %>%dplyr::select(c(1,61))%>%
  dplyr::rename(FeatureID="#FeatureID")
write.csv(annotations,"pfam/pfam_annotationkey.csv", row.names = FALSE)


#remove potential outlier FA09b

dat<-fread("pfam/pfam_clean_noNT.tsv")%>%
  dplyr::select(-cFA09b)
write.table(dat,"pfam/pfam_clean_noNT_rmFA09b.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

###########################################################
###shannon-rarefied
md<-fread("metaT_metadata_ztcat_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))

alpha <- read_qza("pfam/diversity-core-metrics2.3M/shannon_vector.qza")$data %>%
  rownames_to_column("sample_name") %>%
  left_join(.,md,by="sample_name") %>%
  mutate(condition=factor(condition,levels=c("NA","FA","FT")),
         phase=factor(phase,levels=c("light","dark")))

p<-ggplot(alpha, aes(x=condition, y=shannon_entropy, fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
  theme_classic()+
  labs(x="condition",y="shannon distance", title="metaT alpha diversity")+
  theme(legend.position = "none")

ggsave("pfam/diversity-core-metrics2.3M/shannon_metaT.pdf", plot=p,height=2.5, width=3.5)

pairwise.wilcox.test(alpha$shannon, alpha$condition,p.adjust.method="fdr")

# NA      FA  
# FA 3.4e-06 -   
#   FT 6.9e-08 0.22

p<-ggplot(alpha, aes(x=condition, y=shannon_entropy, fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
  facet_wrap(~phase)+
  theme_classic()+
  labs(x="condition",y="shannon distance", title="metaT Alpha Diversity")+
  theme(legend.position = "none", plot.title = element_text(face = "bold"))
ggsave("pfam/diversity-core-metrics2.3M/shannon_metaT_lightdark.pdf", plot=p,height=2.5, width=3.5)

alphaL<-alpha%>%filter(phase=="light")
pairwise.wilcox.test(alphaL$shannon, alphaL$condition,
                     p.adjust.method="fdr")

# NA     FA    
# FA 0.0012 -     
#   FT 0.0012 0.8182


alphaD<-alpha%>%filter(phase=="dark")
pairwise.wilcox.test(alphaD$shannon, alphaD$condition,
                     p.adjust.method="fdr")

# NA     FA    
# FA 0.0072 -     
#   FT 0.0012 0.0931

alphaFT<-alpha%>%filter(condition=="FT")
pairwise.wilcox.test(alphaFT$shannon, alphaFT$phase,
                     p.adjust.method="fdr")
# light
# dark 0.48  

alphaFA<-alpha%>%filter(condition=="FA")
pairwise.wilcox.test(alphaFA$shannon, alphaFA$phase,
                     p.adjust.method="fdr")
# light
# dark 0.7 

alphaNA<-alpha%>%filter(condition=="NA")
pairwise.wilcox.test(alphaNA$shannon, alphaNA$phase,
                     p.adjust.method="fdr")
# light
# dark 0.04

#shannon-rarefied but looking at fasted vs. non-fasted

md<-fread("metaT_metadata_ztcat_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))

alpha <- read_qza("pfam/diversity-core-metrics2.3M/shannon_vector.qza")$data %>%
  rownames_to_column("sample_name") %>%
  left_join(.,md,by="sample_name") %>%
  mutate(fast=ifelse(zt_time<17,"fasted","not_fasted"))%>%
  mutate(condition=factor(condition,levels=c("NA","FA","FT")),
         fast=factor(fast,levels=c("fasted","not_fasted")))

p<-ggplot(alpha, aes(x=condition, y=shannon_entropy, fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  scale_fill_manual(values=c("#009E73","#D55E00","#0072B2"))+
  facet_wrap(~fast)+
  theme_classic()+
  labs(x="condition",y="shannon distance", title="metaT Alpha Diversity")+
  theme(legend.position = "none", plot.title = element_text(face = "bold"))
ggsave("pfam/diversity-core-metrics2.3M/shannon_metaT_fastnotfast.pdf", plot=p,height=2.5, width=3.5)

alphaF<-alpha%>%filter(fast=="fasted")
pairwise.wilcox.test(alphaF$shannon, alphaF$condition,
                     p.adjust.method="fdr")
# NA      FA  
# FA 4.8e-05 -   
#   FT 4.8e-05 0.65

alphaNF<-alpha%>%filter(fast=="not_fasted")
pairwise.wilcox.test(alphaNF$shannon, alphaNF$condition,
                     p.adjust.method="fdr")

# NA    FA   
# FA 0.100 -    
#   FT 0.029 0.114

alphaFT<-alpha%>%filter(condition=="FT")
pairwise.wilcox.test(alphaFT$shannon, alphaFT$fast,
                     p.adjust.method="fdr")
# fasted
# not_fasted 0.81  

alphaFA<-alpha%>%filter(condition=="FA")
pairwise.wilcox.test(alphaFA$shannon, alphaFA$fast,
                     p.adjust.method="fdr")
# fasted
# not_fasted 0.46

alphaNA<-alpha%>%filter(condition=="NA")
pairwise.wilcox.test(alphaNA$shannon, alphaNA$fast,
                     p.adjust.method="fdr")
# fasted
# not_fasted 0.12  
##########################################################
#rpca results

ord <- read_qza("pfam/rpca_results/ordination.qza")

samp_ord<-ord$data$Vectors
write.table(samp_ord,"pfam/rpca_results/sample_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)
feat_ord<-ord$data$Species
write.table(feat_ord,"pfam/rpca_results/feature_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)

md<-fread("metaT_metadata_ztcat_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))

rpca<-ord$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2, PC3)%>%
  dplyr::rename(sample_name=SampleID)%>%
  left_join(md,by="sample_name")%>%
  mutate(condition=factor(condition,levels=c("FT","FA","NA")),
         phase=factor(phase,levels=c("light","dark")))

pdf("pfam/rpca_results/SFR23_0619_metaT_RPCA.pdf",width = 5, height = 5)

# Create a 3D PCA plot
colors <- c("#009E73","#D55E00","#0072B2")
colors <- colors[as.numeric(rpca$condition)]

scatterplot3d(rpca[,2:4],angle = 230,pch = 16,color=colors,grid=FALSE,
              xlab = "PC1", ylab = "PC2", zlab = "PC3",
              main = "MTX")
legend("bottom", legend = levels(rpca$condition),
       col =  c("#009E73","#D55E00","#0072B2"), pch = 16, 
       inset = -0.25, xpd = TRUE, horiz = TRUE)

dev.off()

#plot the ZT time per condition
NA_rpca<-rpca%>%
  filter(condition=="NA")%>%
  mutate(zt_time=factor(zt_time,levels=c("1","5","9","13","17","21")))

p<-NA_rpca %>%
  ggplot(aes(x=PC1, y=PC2, fill=zt_time)) +
  geom_point(alpha=1.0,size=3,shape=21) + 
  theme_pubr() +
  scale_fill_manual(values=c("#67001f","#d6604d","#ffba92","#8ccff3","#4393c3","#053061"))+
  labs(color="ZT time",
       x =paste("PC1 (",round(ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""))+ggtitle("MTX NA")+ theme(plot.title = element_text(face = "bold"))

ggsave("pfam/rpca_results/SFR24_0412_mtx_NA_RPCA.pdf", plot=p,height=3.5, width=3.5)

FA_rpca<-rpca%>%
  filter(condition=="FA")%>%
  mutate(zt_time=factor(zt_time,levels=c("1","5","9","13","17","21")))

p<-FA_rpca %>%
  ggplot(aes(x=PC1, y=PC2, fill=zt_time)) +
  geom_point(alpha=1.0,size=3,shape=21) + 
  theme_pubr() +
  scale_fill_manual(values=c("#67001f","#d6604d","#ffba92","#8ccff3","#4393c3","#053061"))+
  labs(color="ZT time",
       x =paste("PC1 (",round(ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""))+ggtitle("MTX FA")+ theme(plot.title = element_text(face = "bold"))

ggsave("pfam/rpca_results/SFR24_0412_mtx_FA_RPCA.pdf", plot=p,height=3.5, width=3.5)

FT_rpca<-rpca%>%
  filter(condition=="FT")%>%
  mutate(zt_time=factor(zt_time,levels=c("1","5","9","13","17","21")))

p<-FT_rpca %>%
  ggplot(aes(x=PC1, y=PC2, fill=zt_time)) +
  geom_point(alpha=1.0,size=3,shape=21) + 
  theme_pubr() +
  scale_fill_manual(values=c("#67001f","#d6604d","#ffba92","#8ccff3","#4393c3","#053061"))+
  labs(color="ZT time",
       x =paste("PC1 (",round(ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""))+ggtitle("MTX FT")+ theme(plot.title = element_text(face = "bold"))

ggsave("pfam/rpca_results/SFR24_0412_mtx_FT_RPCA.pdf", plot=p,height=3.5, width=3.5)

###########################################################

#make nat log plots

md<-fread("metaT_metadata_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%dplyr::rename(`Sample ID`=sample_name)


natlog_lipid<-fread("ordination_deicode_TRFmetaT_noNT/qurro_sample_plot_data_lipid_biometab_transcrpt.tsv")%>%
  dplyr::select(1:2)%>%
  left_join(.,md,by="Sample ID")%>%
  mutate(fast=ifelse(zt_time<17,"fasted","not_fasted"))%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  mutate(condition=factor(condition,levels=c("NA","FA","FT")),
         lightdark=factor(lightdark,levels=c("light","dark")),
         fast=factor(fast,levels=c("fasted","not_fasted")))

p<-ggplot(natlog_lipid, aes(x=condition, y=Current_Natural_Log_Ratio,fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  theme_minimal()+ scale_color_manual(values=c("#009E73","#D55E00","#0072B2"))+scale_fill_manual(values=c("#009E73","#D55E00","#0072B2"))+
  labs(x="condition",y="Natural Log Ratio",title="Lipid Biosynthetic & Metabolic Process")+
  theme(legend.position = "none")

pairwise.wilcox.test(natlog_lipid$Current_Natural_Log_Ratio, natlog_lipid$condition,
                     p.adjust.method="fdr")
# NA     FA    
# FA 0.2486 -     
#   FT 0.0049 0.0674

ggsave("ordination_deicode_TRFmetaT_noNT/SFR23_0405_natlog_lipidvstranscrpt.pdf", plot=p,height=3.5, width=3.5)

p<-ggplot(natlog_lipid, aes(x=condition, y=Current_Natural_Log_Ratio,fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  facet_wrap(~lightdark)+
  theme_classic()+ scale_color_manual(values=c("#009E73","#D55E00","#0072B2"))+scale_fill_manual(values=c("#009E73","#D55E00","#0072B2"))+
  labs(x="condition",y="Natural Log Ratio",title="Lipid Biosynthetic & Metabolic Process")+
  theme(legend.position = "none")
ggsave("ordination_deicode_TRFmetaT_noNT/SFR23_0405_natlog_lipidvstranscrpt_lightdark.pdf", plot=p,height=3.5, width=5)

natlog_lipidL<-natlog_lipid%>%filter(lightdark=="light")
pairwise.wilcox.test(natlog_lipidL$Current_Natural_Log_Ratio, natlog_lipidL$condition,
                     p.adjust.method="fdr")

# NA   FA  
# FA 0.78 -   
#   FT 0.20 0.73

natlog_lipidD<-natlog_lipid%>%filter(lightdark=="dark")
pairwise.wilcox.test(natlog_lipidD$Current_Natural_Log_Ratio, natlog_lipidD$condition,
                     p.adjust.method="fdr")

# NA    FA   
# FA 0.388 -    
#   FT 0.054 0.054

natlog_lipidFT<-natlog_lipid%>%filter(condition=="FT")
pairwise.wilcox.test(natlog_lipidFT$Current_Natural_Log_Ratio, natlog_lipidFT$lightdark,
                     p.adjust.method="fdr")

# light
# dark 0.18

natlog_lipidFA<-natlog_lipid%>%filter(condition=="FA")
pairwise.wilcox.test(natlog_lipidFA$Current_Natural_Log_Ratio, natlog_lipidFA$lightdark,
                     p.adjust.method="fdr")

# light
# dark 0.7

natlog_lipidNA<-natlog_lipid%>%filter(condition=="NA")
pairwise.wilcox.test(natlog_lipidNA$Current_Natural_Log_Ratio, natlog_lipidNA$lightdark,
                     p.adjust.method="fdr")

# light
# dark 0.67 

p<-ggplot(natlog_lipid, aes(x=condition, y=Current_Natural_Log_Ratio,fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  facet_wrap(~fast)+
  theme_classic()+ scale_color_manual(values=c("#009E73","#D55E00","#0072B2"))+scale_fill_manual(values=c("#009E73","#D55E00","#0072B2"))+
  labs(x="condition",y="Natural Log Ratio",title="Lipid Biosynthetic & Metabolic Process")+
  theme(legend.position = "none")
ggsave("ordination_deicode_TRFmetaT_noNT/SFR23_0405_natlog_lipidvstranscrpt_fastnotfast.pdf", plot=p,height=3.5, width=5)

natlog_lipidF<-natlog_lipid%>%filter(fast=="fasted")
pairwise.wilcox.test(natlog_lipidF$Current_Natural_Log_Ratio, natlog_lipidF$condition,
                     p.adjust.method="fdr")

# NA    FA   
# FA 0.234 -    
#   FT 0.009 0.234

natlog_lipidNF<-natlog_lipid%>%filter(fast=="not_fasted")
pairwise.wilcox.test(natlog_lipidNF$Current_Natural_Log_Ratio, natlog_lipidNF$condition,
                     p.adjust.method="fdr")

# NA   FA  
# FA 0.91 -   
#   FT 0.39 0.34

natlog_lipidFT<-natlog_lipid%>%filter(condition=="FT")
pairwise.wilcox.test(natlog_lipidFT$Current_Natural_Log_Ratio, natlog_lipidFT$fast,
                     p.adjust.method="fdr")

# fasted
# not_fasted 0.68

natlog_lipidFA<-natlog_lipid%>%filter(condition=="FA")
pairwise.wilcox.test(natlog_lipidFA$Current_Natural_Log_Ratio, natlog_lipidFA$fast,
                     p.adjust.method="fdr")

# fasted
# not_fasted 0.93 

natlog_lipidNA<-natlog_lipid%>%filter(condition=="NA")
pairwise.wilcox.test(natlog_lipidNA$Current_Natural_Log_Ratio, natlog_lipidNA$fast,
                     p.adjust.method="fdr")

# fasted
# not_fasted 0.55 

natlog_carbo<-fread("ordination_deicode_TRFmetaT_noNT/qurro_sample_plot_data_carb_metabproc_transcrpt.tsv")%>%
  dplyr::select(1:2)%>%
  left_join(.,md,by="Sample ID")%>%
  mutate(fast=ifelse(zt_time<17,"fasted","not_fasted"))%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  mutate(condition=factor(condition,levels=c("NA","FA","FT")),
         lightdark=factor(lightdark,levels=c("light","dark")),
         fast=factor(fast,levels=c("fasted","not_fasted")))

p<-ggplot(natlog_carbo, aes(x=condition, y=Current_Natural_Log_Ratio,fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  theme_minimal()+ scale_color_manual(values=c("#009E73","#D55E00","#0072B2"))+scale_fill_manual(values=c("#009E73","#D55E00","#0072B2"))+
  labs(x="condition",y="Natural Log Ratio",title="Carbohydrate Metabolic Process")+
  theme(legend.position = "none")

ggsave("ordination_deicode_TRFmetaT_noNT/SFR23_0405_natlog_carbovstranscrpt.pdf", plot=p,height=3.5, width=3.5)

pairwise.wilcox.test(natlog_carbo$Current_Natural_Log_Ratio, natlog_carbo$condition,
                     p.adjust.method="fdr")
# NA      FA     
# FA 0.09741 -      
#   FT 0.00018 0.11350

p<-ggplot(natlog_carbo, aes(x=condition, y=Current_Natural_Log_Ratio,fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  facet_wrap(~lightdark)+
  theme_classic()+ scale_color_manual(values=c("#009E73","#D55E00","#0072B2"))+scale_fill_manual(values=c("#009E73","#D55E00","#0072B2"))+
  labs(x="condition",y="Natural Log Ratio",title="Carbohydrate Metabolic Process")+
  theme(legend.position = "none")

ggsave("ordination_deicode_TRFmetaT_noNT/SFR23_0405_natlog_carbovstranscrpt_lightdark.pdf", plot=p,height=3.5, width=5)

natlog_carboL<-natlog_carbo%>%filter(lightdark=="light")
pairwise.wilcox.test(natlog_carboL$Current_Natural_Log_Ratio, natlog_carboL$condition,
                     p.adjust.method="fdr")

# NA     FA    
# FA 0.3095 -     
#   FT 0.0012 0.3095

natlog_carboD<-natlog_carbo%>%filter(lightdark=="dark")
pairwise.wilcox.test(natlog_carboD$Current_Natural_Log_Ratio, natlog_carboD$condition,
                     p.adjust.method="fdr")

# NA   FA  
# FA 0.33 -   
#   FT 0.15 0.33

natlog_carboFT<-natlog_carbo%>%filter(condition=="FT")
pairwise.wilcox.test(natlog_carboFT$Current_Natural_Log_Ratio, natlog_carboFT$lightdark,
                     p.adjust.method="fdr")

# light
# dark 0.82

natlog_carboFA<-natlog_carbo%>%filter(condition=="FA")
pairwise.wilcox.test(natlog_carboFA$Current_Natural_Log_Ratio, natlog_carboFA$lightdark,
                     p.adjust.method="fdr")

# light
# dark 0.82

natlog_carboNA<-natlog_carbo%>%filter(condition=="NA")
pairwise.wilcox.test(natlog_carboNA$Current_Natural_Log_Ratio, natlog_carboNA$lightdark,
                     p.adjust.method="fdr")

# light
# dark 0.024

p<-ggplot(natlog_carbo, aes(x=condition, y=Current_Natural_Log_Ratio,fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  facet_wrap(~fast)+
  theme_classic()+ scale_color_manual(values=c("#009E73","#D55E00","#0072B2"))+scale_fill_manual(values=c("#009E73","#D55E00","#0072B2"))+
  labs(x="condition",y="Natural Log Ratio",title="Carbohydrate Metabolic Process")+
  theme(legend.position = "none")

ggsave("ordination_deicode_TRFmetaT_noNT/SFR23_0405_natlog_carbovstranscrpt_fastnotfast.pdf", plot=p,height=3.5, width=5)

natlog_carboF<-natlog_carbo%>%filter(fast=="fasted")
pairwise.wilcox.test(natlog_carboF$Current_Natural_Log_Ratio, natlog_carboF$condition,
                     p.adjust.method="fdr")

# NA     FA    
# FA 0.3054 -     
#   FT 0.0032 0.1245

natlog_carboF<-natlog_carbo%>%filter(fast=="not_fasted")
pairwise.wilcox.test(natlog_carboF$Current_Natural_Log_Ratio, natlog_carboF$condition,
                     p.adjust.method="fdr")

# NA   FA  
# FA 0.10 -   
#   FT 0.10 0.89

natlog_carboFT<-natlog_carbo%>%filter(condition=="FT")
pairwise.wilcox.test(natlog_carboFT$Current_Natural_Log_Ratio, natlog_carboFT$fast,
                     p.adjust.method="fdr")

# fasted
# not_fasted 1 

natlog_carboFA<-natlog_carbo%>%filter(condition=="FA")
pairwise.wilcox.test(natlog_carboFA$Current_Natural_Log_Ratio, natlog_carboFA$fast,
                     p.adjust.method="fdr")

# fasted
# not_fasted 0.37

natlog_carboNA<-natlog_carbo%>%filter(condition=="NA")
pairwise.wilcox.test(natlog_carboNA$Current_Natural_Log_Ratio, natlog_carboNA$fast,
                     p.adjust.method="fdr")

# fasted
# not_fasted 0.44

natlog_protlys<-fread("ordination_deicode_TRFmetaT_noNT/qurro_sample_plot_data_proteolysis_transcrpt.tsv")%>%
  dplyr::select(1:2)%>%
  left_join(.,md,by="Sample ID")%>%
  mutate(fast=ifelse(zt_time<17,"fasted","not_fasted"))%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  mutate(condition=factor(condition,levels=c("NA","FA","FT")),
         lightdark=factor(lightdark,levels=c("light","dark")),
         fast=factor(fast,levels=c("fasted","not_fasted")))

p<-ggplot(natlog_protlys, aes(x=condition, y=Current_Natural_Log_Ratio,fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  theme_minimal()+ scale_color_manual(values=c("#009E73","#D55E00","#0072B2"))+scale_fill_manual(values=c("#009E73","#D55E00","#0072B2"))+
  labs(x="condition",y="Natural Log Ratio",title="Proteolysis")+
  theme(legend.position = "none")

pairwise.wilcox.test(natlog_protlys$Current_Natural_Log_Ratio, natlog_protlys$condition,
                     p.adjust.method="fdr")
# NA     FA    
# FA 0.9834 -     
#   FT 0.0068 0.0259

ggsave("ordination_deicode_TRFmetaT_noNT/SFR23_0405_natlog_protlysvstranscrpt.pdf", plot=p,height=3.5, width=3.5)

p<-ggplot(natlog_protlys, aes(x=condition, y=Current_Natural_Log_Ratio,fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  facet_wrap(~lightdark)+
  theme_classic()+ scale_color_manual(values=c("#009E73","#D55E00","#0072B2"))+scale_fill_manual(values=c("#009E73","#D55E00","#0072B2"))+
  labs(x="condition",y="Natural Log Ratio",title="Proteolysis")+
  theme(legend.position = "none")

ggsave("ordination_deicode_TRFmetaT_noNT/SFR23_0405_natlog_protlysvstranscrpt_lightdark.pdf", plot=p,height=3.5, width=5)

natlog_protlysL<-natlog_protlys%>%filter(lightdark=="light")
pairwise.wilcox.test(natlog_protlysL$Current_Natural_Log_Ratio, natlog_protlysL$condition,
                     p.adjust.method="fdr")

# NA   FA  
# FA 0.61 -   
#   FT 0.59 0.59

natlog_protlysD<-natlog_protlys%>%filter(lightdark=="dark")
pairwise.wilcox.test(natlog_protlysD$Current_Natural_Log_Ratio, natlog_protlysD$condition,
                     p.adjust.method="fdr")

# NA     FA    
# FA 0.6889 -     
#   FT 0.0048 0.0617

natlog_protlysFT<-natlog_protlys%>%filter(condition=="FT")
pairwise.wilcox.test(natlog_protlysFT$Current_Natural_Log_Ratio, natlog_protlysFT$lightdark,
                     p.adjust.method="fdr")

# light
# dark 0.48

natlog_protlysFA<-natlog_protlys%>%filter(condition=="FA")
pairwise.wilcox.test(natlog_protlysFA$Current_Natural_Log_Ratio, natlog_protlysFA$lightdark,
                     p.adjust.method="fdr")

# light
# dark 0.24

natlog_protlysNA<-natlog_protlys%>%filter(condition=="NA")
pairwise.wilcox.test(natlog_protlysNA$Current_Natural_Log_Ratio, natlog_protlysNA$lightdark,
                     p.adjust.method="fdr")

# light
# dark 0.22

p<-ggplot(natlog_protlys, aes(x=condition, y=Current_Natural_Log_Ratio,fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  facet_wrap(~fast)+
  theme_classic()+ scale_color_manual(values=c("#009E73","#D55E00","#0072B2"))+scale_fill_manual(values=c("#009E73","#D55E00","#0072B2"))+
  labs(x="condition",y="Natural Log Ratio",title="Proteolysis")+
  theme(legend.position = "none")

ggsave("ordination_deicode_TRFmetaT_noNT/SFR23_0405_natlog_protlysvstranscrpt_fastnotfast.pdf", plot=p,height=3.5, width=5)

natlog_protlysF<-natlog_protlys%>%filter(fast=="fasted")
pairwise.wilcox.test(natlog_protlysF$Current_Natural_Log_Ratio, natlog_protlysF$condition,
                     p.adjust.method="fdr")

# NA   FA  
# FA 0.31 -   
#   FT 0.12 0.33

natlog_protlysNF<-natlog_protlys%>%filter(fast=="not_fasted")
pairwise.wilcox.test(natlog_protlysNF$Current_Natural_Log_Ratio, natlog_protlysNF$condition,
                     p.adjust.method="fdr")

# NA    FA   
# FA 0.067 -    
#   FT 0.057 0.057

natlog_protlysFT<-natlog_protlys%>%filter(condition=="FT")
pairwise.wilcox.test(natlog_protlysFT$Current_Natural_Log_Ratio, natlog_protlysFT$fast,
                     p.adjust.method="fdr")

# fasted
# not_fasted 0.81

natlog_protlysFA<-natlog_protlys%>%filter(condition=="FA")
pairwise.wilcox.test(natlog_protlysFA$Current_Natural_Log_Ratio, natlog_protlysFA$fast,
                     p.adjust.method="fdr")

# fasted
# not_fasted 0.048

natlog_protlysNA<-natlog_protlys%>%filter(condition=="NA")
pairwise.wilcox.test(natlog_protlysNA$Current_Natural_Log_Ratio, natlog_protlysNA$fast,
                     p.adjust.method="fdr")

# fasted
# not_fasted 0.62 

#now plot BSH over transcript genes

natlog_bsh<-fread("ordination_deicode_TRFmetaT_noNT/qurro_sample_plot_data_BSH_transcrpt.tsv")%>%
#natlog_bsh<-fread("ordination_deicode_TRFmetaT_noNT/qurro_sample_plot_data_BSH_RPOB.tsv")%>%
  dplyr::select(1:2)%>%
  left_join(.,md,by="Sample ID")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  mutate(condition=factor(condition,levels=c("FT","FA","NA")),
         lightdark=factor(lightdark,levels=c("light","dark")),
         Current_Natural_Log_Ratio=as.numeric(Current_Natural_Log_Ratio))

p<-ggplot(natlog_bsh, aes(x=condition, y=Current_Natural_Log_Ratio,fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  theme_minimal()+ scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
  labs(x="condition",y="Natural Log Ratio",title="Bile Salt Hydrolase (BSH)")+
  theme(legend.position = "none")

pairwise.wilcox.test(natlog_bsh$Current_Natural_Log_Ratio, natlog_bsh$condition,
                     p.adjust.method="fdr")
# FT      FA     
# FA 0.79    -      
#   NA 1.7e-07 2.4e-07

ggsave("ordination_deicode_TRFmetaT_noNT/SFR23_0220_natlog_BSHvstranscrpt.pdf", plot=p,height=3.5, width=3.5)

p<-ggplot(natlog_bsh, aes(x=condition, y=Current_Natural_Log_Ratio,fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  facet_wrap(~lightdark)+
  theme_classic()+ scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
  labs(x="condition",y="Natural Log Ratio",title="Bile Salt Hydrolase (BSH)")+
  theme(legend.position = "none")

ggsave("ordination_deicode_TRFmetaT_noNT/SFR23_0220_natlog_BSHvstranscrpt_lightdark.pdf", plot=p,height=3.5, width=5)

natlog_bshL<-natlog_bsh%>%filter(lightdark=="light")
pairwise.wilcox.test(natlog_bshL$Current_Natural_Log_Ratio, natlog_bshL$condition,
                     p.adjust.method="fdr")

# FT     FA    
# FA 0.3290 -     
#   NA 0.0015 0.0015

natlog_bshD<-natlog_bsh%>%filter(lightdark=="dark")
pairwise.wilcox.test(natlog_bshD$Current_Natural_Log_Ratio, natlog_bshD$condition,
                     p.adjust.method="fdr")

# FT     FA    
# FA 0.3939 -     
#   NA 0.0006 0.0006

natlog_bshFT<-natlog_bsh%>%filter(condition=="FT")
pairwise.wilcox.test(natlog_bshFT$Current_Natural_Log_Ratio, natlog_bshFT$lightdark,
                     p.adjust.method="fdr")

# light
# dark 0.33

natlog_bshFA<-natlog_bsh%>%filter(condition=="FA")
pairwise.wilcox.test(natlog_bshFA$Current_Natural_Log_Ratio, natlog_bshFA$lightdark,
                     p.adjust.method="fdr")

# light
# dark 0.24

natlog_bshNA<-natlog_bsh%>%filter(condition=="NA")
pairwise.wilcox.test(natlog_bshNA$Current_Natural_Log_Ratio, natlog_bshNA$lightdark,
                     p.adjust.method="fdr")

# light
# dark 0.34
