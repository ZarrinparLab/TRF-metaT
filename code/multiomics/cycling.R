setwd("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/multiomics")

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

#load MTX cycling file

#annotations
pfam_annot<-fread("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_m_results/pfam/pfam_annotationkey.csv")

#read in FT metacycle hits
#FT_metacyc_MT<-fread("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_m_results/pfam/cyclic_analysis/wol2_rna_pipeline_rpob/FT_metacycle/meta2d_filtered_rna_FT.txt")%>%
FT_metacyc_MT<-fread("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_m_results/pfam/cyclic_analysis/wol2_rna_pipeline_144k_rpob/FT_metacycle/meta2d_filtered_rna_FT.txt")%>%
  dplyr::rename(FeatureID=CycID) %>%
  left_join(.,pfam_annot, by="FeatureID")%>%
  arrange(JTK_adjphase)%>%
  mutate(label_name=paste(FeatureID, Name, sep=" "))
sigFTMT<-FT_metacyc_MT%>%filter(JTK_pvalue<0.05,!grepl("DUF",Name)) #662 #122

#FA_metacyc_MT<-fread("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_m_results/pfam/cyclic_analysis/wol2_rna_pipeline_rpob/FA_metacycle/meta2d_filtered_rna_FA.txt")%>%
FA_metacyc_MT<-fread("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_m_results/pfam/cyclic_analysis/wol2_rna_pipeline_144k_rpob/FA_metacycle/meta2d_filtered_rna_FA.txt")%>%
  dplyr::rename(FeatureID=CycID) %>%
  left_join(.,pfam_annot, by="FeatureID")%>%
  arrange(JTK_adjphase)%>%
  mutate(label_name=paste(FeatureID, Name, sep=" "))
sigFAMT<-FA_metacyc_MT%>%filter(JTK_pvalue<0.05,!grepl("DUF",Name)) #117 #56

#NA_metacyc_MT<-fread("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_m_results/pfam/cyclic_analysis/wol2_rna_pipeline_rpob/NA_metacycle/meta2d_filtered_rna_NA.txt")%>%
NA_metacyc_MT<-fread("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_m_results/pfam/cyclic_analysis/wol2_rna_pipeline_144k_rpob/NA_metacycle/meta2d_filtered_rna_NA.txt")%>%
  dplyr::rename(FeatureID=CycID) %>%
  left_join(.,pfam_annot, by="FeatureID")%>%
  arrange(JTK_adjphase)%>%
  mutate(label_name=paste(FeatureID, Name, sep=" "))
sigNAMT<-NA_metacyc_MT%>%filter(JTK_pvalue<0.05,!grepl("DUF",Name)) #571 #292

#load MGX cycling figure

FT_metacyc_MG<-fread("~/scratch/TRF_multiomics/metatranscript/woltka2_results/cyclic_analysis/wol2_genome_pipeline_rpob_SFR/FT_metacycle/meta2d_filtered_genome_FT.txt")%>%
  dplyr::rename(FeatureID=CycID) %>%
  left_join(.,pfam_annot, by="FeatureID")%>%
  arrange(JTK_adjphase)%>%
  mutate(label_name=paste(FeatureID, Name, sep=" "))
sigFTMG<-FT_metacyc_MG%>%filter(JTK_pvalue<0.05,!grepl("DUF",Name)) #45

FA_metacyc_MG<-fread("~/scratch/TRF_multiomics/metatranscript/woltka2_results/cyclic_analysis/wol2_genome_pipeline_rpob_SFR/FA_metacycle/meta2d_filtered_genome_FA.txt")%>%
  dplyr::rename(FeatureID=CycID) %>%
  left_join(.,pfam_annot, by="FeatureID")%>%
  arrange(JTK_adjphase)%>%
  mutate(label_name=paste(FeatureID, Name, sep=" "))
sigFAMG<-FA_metacyc_MG%>%filter(JTK_pvalue<0.05,!grepl("DUF",Name)) #42

NA_metacyc_MG<-fread("~/scratch/TRF_multiomics/metatranscript/woltka2_results/cyclic_analysis/wol2_genome_pipeline_rpob_SFR/NA_metacycle/meta2d_filtered_genome_NA.txt")%>%
  dplyr::rename(FeatureID=CycID) %>%
  left_join(.,pfam_annot, by="FeatureID")%>%
  arrange(JTK_adjphase)%>%
  mutate(label_name=paste(FeatureID, Name, sep=" "))
sigNAMG<-NA_metacyc_MG%>%filter(JTK_pvalue<0.05,!grepl("DUF",Name)) #308

jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

j_val<-c(jaccard(sigNAMG$FeatureID,sigNAMT$FeatureID),jaccard(sigNAMG$FeatureID,sigFAMT$FeatureID),jaccard(sigNAMG$FeatureID,sigFTMT$FeatureID),
         jaccard(sigFAMG$FeatureID,sigFAMT$FeatureID),jaccard(sigFAMG$FeatureID,sigFTMT$FeatureID),
         jaccard(sigFTMG$FeatureID,sigFTMT$FeatureID))


var1<-c(rep("NA-MGX",3),rep("FA-MGX",2),"FT-MGX")
var2<-c("NA-MTX","FA-MTX","FT-MTX",
        "FA-MTX","FT-MTX",
        "FT-MTX")
jdf<-data.frame(var1,var2,j_val)%>%mutate(var1=factor(var1,levels=c("NA-MGX","FA-MGX","FT-MGX")),
                                          var2=factor(var2,levels=c("FT-MTX","FA-MTX","NA-MTX")))

p<-ggplot(jdf, aes(x = var1, y = var2, fill = j_val)) +
  geom_tile() +
  geom_text(aes(label = signif(j_val,digits=3))) +
  scale_fill_distiller(palette = "Blues",direction=1) +
  theme_minimal() +
  theme(panel.grid = element_blank(),legend.position = "top")+labs(title="Jaccard Index",x="",y="")
#ggsave("SFR23_0619_jaccard_metacy_MTrMG_summ_small.pdf",plot=p, width = 3, height = 3.5)
ggsave("SFR23_0620_jaccard_metacy_MT144krMG_summ_small.pdf",plot=p, width = 3, height = 3.5)


##############################################

#show MTX and MGX summ of unique cycling hits

#for MTX
sigFTMT<-FT_metacyc_MT%>%filter(JTK_pvalue<0.05,!grepl("DUF",Name)) #369 #573
sigFAMT<-FA_metacyc_MT%>%filter(JTK_pvalue<0.05,!grepl("DUF",Name)) #105 #134
sigNAMT<-NA_metacyc_MT%>%filter(JTK_pvalue<0.05,!grepl("DUF",Name)) #513 #518

list_venn <- list(NA_ = sigNAMT$FeatureID,
                  FA = sigFAMT$FeatureID,
                  FT = sigFTMT$FeatureID)

ItemsList <- venn(list_venn, show.plot = FALSE)
all<-attributes(ItemsList)$intersections


#for MGX
sigFTMG<-FT_metacyc_MG%>%filter(JTK_pvalue<0.05,!grepl("DUF",Name)) #59 #45
sigFAMG<-FA_metacyc_MG%>%filter(JTK_pvalue<0.05,!grepl("DUF",Name)) #68 #42
sigNAMG<-NA_metacyc_MG%>%filter(JTK_pvalue<0.05,!grepl("DUF",Name)) #331 #398

list_vennG <- list(NA_ = sigNAMG$FeatureID,
                  FA = sigFAMG$FeatureID,
                  FT = sigFTMG$FeatureID)

ItemsListG <- venn(list_vennG, show.plot = FALSE)
allG<-attributes(ItemsListG)$intersections

#FT
FT_cyc<-pfam_annot%>%filter(FeatureID %in% all$FT) #293 #474

gonames<-fread("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_results/go_name.txt")
pfamGOp<-read.table("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_results/pfam-to-go-process.map",header = FALSE, sep = "\t",
                    col.names = paste0("V",seq_len(4)), fill = TRUE)%>%
  gather(column,GO_Term,-V1)%>%
  dplyr::select(1,3)%>%
  dplyr::rename(FeatureID=V1)%>%
  filter(FeatureID %in% FT_cyc$FeatureID)%>% #we dont get any GO Term matches
  left_join(.,gonames,by="GO_Term")%>%
  filter(!is.na(name))

FTcyc_summ<-pfamGOp%>%
  group_by(name)%>%summarise(n=n())%>%
  mutate(condition="FT")%>%
  arrange(n) %>%
  mutate(method="MTX")

FT_cycG<-pfam_annot%>%filter(FeatureID %in% allG$FT) #51 #41

gonames<-fread("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_results/go_name.txt")
pfamGOp2<-read.table("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_results/pfam-to-go-process.map",header = FALSE, sep = "\t",
                    col.names = paste0("V",seq_len(4)), fill = TRUE)%>%
  gather(column,GO_Term,-V1)%>%
  dplyr::select(1,3)%>%
  dplyr::rename(FeatureID=V1)%>%
  filter(FeatureID %in% FT_cycG$FeatureID)%>% 
  left_join(.,gonames,by="GO_Term")%>%
  filter(!is.na(name))

FTcyc_summ2<-pfamGOp2%>%
  group_by(name)%>%summarise(n=n())%>%
  mutate(condition="FT")%>%
  arrange(n) %>%
  mutate(method="MGX")

FTcyc_summComb<-rbind(FTcyc_summ,FTcyc_summ2)%>%
  group_by(name)%>%mutate(total=sum(n))%>%
  arrange(total)#%>%
  #filter(n>1)

FTcyc_summComb$name <- factor(FTcyc_summComb$name,levels = unique(FTcyc_summComb$name))

ggplot(data=FTcyc_summComb, aes(x=name, y=n, fill=method)) +
  geom_bar(stat="identity", position="stack") + coord_flip() +
  scale_fill_manual(values=c("red","blue")) +theme_pubr() +
  theme(axis.text = element_text(size = 7))  +
  scale_y_continuous(expand=c(0,0), limits=c(0,10))

ggsave("SFR23_0601_cycFTMTMG_GOterms.pdf",height=7, width=6)
#ggsave("SFR23_0522_cycFTMTMG_GOterms_nmor1.pdf",height=4, width=8)


#FA
FA_cyc<-pfam_annot%>%filter(FeatureID %in% all$FA) #71 #89

gonames<-fread("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_results/go_name.txt")
pfamGOp<-read.table("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_results/pfam-to-go-process.map",header = FALSE, sep = "\t",
                    col.names = paste0("V",seq_len(4)), fill = TRUE)%>%
  gather(column,GO_Term,-V1)%>%
  dplyr::select(1,3)%>%
  dplyr::rename(FeatureID=V1)%>%
  filter(FeatureID %in% FA_cyc$FeatureID)%>% 
  left_join(.,gonames,by="GO_Term")%>%
  filter(!is.na(name))

FAcyc_summ<-pfamGOp%>%
  group_by(name)%>%summarise(n=n())%>%
  mutate(condition="FA")%>%
  arrange(n) %>%
  mutate(method="MTX")

FA_cycG<-pfam_annot%>%filter(FeatureID %in% allG$FA) #60 #38

gonames<-fread("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_results/go_name.txt")
pfamGOp2<-read.table("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_results/pfam-to-go-process.map",header = FALSE, sep = "\t",
                     col.names = paste0("V",seq_len(4)), fill = TRUE)%>%
  gather(column,GO_Term,-V1)%>%
  dplyr::select(1,3)%>%
  dplyr::rename(FeatureID=V1)%>%
  filter(FeatureID %in% FA_cycG$FeatureID)%>% 
  left_join(.,gonames,by="GO_Term")%>%
  filter(!is.na(name))

FAcyc_summ2<-pfamGOp2%>%
  group_by(name)%>%summarise(n=n())%>%
  mutate(condition="FA")%>%
  arrange(n) %>%
  mutate(method="MGX")

FAcyc_summComb<-rbind(FAcyc_summ,FAcyc_summ2)%>%
  group_by(name)%>%mutate(total=sum(n))%>%
  arrange(total)

FAcyc_summComb$name <- factor(FAcyc_summComb$name,levels = unique(FAcyc_summComb$name))

ggplot(data=FAcyc_summComb, aes(x=name, y=n, fill=method)) +
  geom_bar(stat="identity", position="stack") + coord_flip() +
  scale_fill_manual(values=c("red","blue")) +theme_pubr() +
  scale_y_continuous(expand=c(0,0), limits=c(0,10))

ggsave("SFR23_0601_cycFAMTMG_GOterms.pdf",height=5, width=8)

#NA
NA_cyc<-pfam_annot%>%filter(FeatureID %in% all$NA_) #459 #433

gonames<-fread("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_results/go_name.txt")
pfamGOp<-read.table("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_results/pfam-to-go-process.map",header = FALSE, sep = "\t",
                    col.names = paste0("V",seq_len(4)), fill = TRUE)%>%
  gather(column,GO_Term,-V1)%>%
  dplyr::select(1,3)%>%
  dplyr::rename(FeatureID=V1)%>%
  filter(FeatureID %in% NA_cyc$FeatureID)%>% 
  left_join(.,gonames,by="GO_Term")%>%
  filter(!is.na(name))

NAcyc_summ<-pfamGOp%>%
  group_by(name)%>%summarise(n=n())%>%
  mutate(condition="NA")%>%
  arrange(n) %>%
  mutate(method="MTX")

NA_cycG<-pfam_annot%>%filter(FeatureID %in% allG$NA_) #315

gonames<-fread("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_results/go_name.txt")
pfamGOp2<-read.table("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_results/pfam-to-go-process.map",header = FALSE, sep = "\t",
                     col.names = paste0("V",seq_len(4)), fill = TRUE)%>%
  gather(column,GO_Term,-V1)%>%
  dplyr::select(1,3)%>%
  dplyr::rename(FeatureID=V1)%>%
  filter(FeatureID %in% NA_cycG$FeatureID)%>% 
  left_join(.,gonames,by="GO_Term")%>%
  filter(!is.na(name))

NAcyc_summ2<-pfamGOp2%>%
  group_by(name)%>%summarise(n=n())%>%
  mutate(condition="NA")%>%
  arrange(n) %>%
  mutate(method="MGX")

NAcyc_summComb<-rbind(NAcyc_summ,NAcyc_summ2)%>%
  group_by(name)%>%mutate(total=sum(n))%>%
  arrange(total)#%>%
  #filter(n>1)

NAcyc_summComb$name <- factor(NAcyc_summComb$name,levels = unique(NAcyc_summComb$name))

ggplot(data=NAcyc_summComb, aes(x=name, y=n, fill=method)) +
  geom_bar(stat="identity", position="identity") + coord_flip() +
  scale_fill_manual(values=c("red","blue")) +theme_pubr() +
  scale_y_continuous(expand=c(0,0),limits=c(0,25))

ggsave("SFR23_0601_cycNAMTMG_GOterms.pdf",height=14, width=8)

