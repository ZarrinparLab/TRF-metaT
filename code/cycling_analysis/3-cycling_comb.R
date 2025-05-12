setwd("~/Notebooks/sfloresr/TRF-metaT/")

library(tidyverse)
library(data.table)

##########################################################
#paths
dat_annot<-"data/pfam_metaT/pfam_annotationkey.csv"
dat_go<-"data/pfam_metaT/go_name.txt"
dat_pfamtogo<-"data/pfam_metaT/pfam-to-go-process.map"
dat_mtx_FTcyc<-"data/cycling_analysis/mtx_cyclic_analysis/wol2_rna_pipeline_rpob/FT_metacycle/meta2d_filtered_rna_FT.txt"
dat_mtx_FAcyc<-"data/cycling_analysis/mtx_cyclic_analysis/wol2_rna_pipeline_rpob/FA_metacycle/meta2d_filtered_rna_FA.txt"
dat_mtx_NAcyc<-"data/cycling_analysis/mtx_cyclic_analysis/wol2_rna_pipeline_rpob/NA_metacycle/meta2d_filtered_rna_NA.txt"
dat_mgx_FTcyc<-"data/cycling_analysis/mgx_cyclic_analysis/wol2_dna_pipeline_rpob/FT_metacycle/meta2d_filtered_dna_FT.txt"
dat_mgx_FAcyc<-"data/cycling_analysis/mgx_cyclic_analysis/wol2_dna_pipeline_rpob/FA_metacycle/meta2d_filtered_dna_FA.txt"
dat_mgx_NAcyc<-"data/cycling_analysis/mgx_cyclic_analysis/wol2_dna_pipeline_rpob/NA_metacycle/meta2d_filtered_dna_NA.txt"
dat_mtx_FTcyc_144k<-"data/cycling_analysis/mtx_cyclic_analysis/wol2_rna_pipeline_144k_rpob/FT_metacycle/meta2d_filtered_rna_FT.txt"
dat_mtx_FAcyc_144k<-"data/cycling_analysis/mtx_cyclic_analysis/wol2_rna_pipeline_144k_rpob/FA_metacycle/meta2d_filtered_rna_FA.txt"
dat_mtx_NAcyc_144k<-"data/cycling_analysis/mtx_cyclic_analysis/wol2_rna_pipeline_144k_rpob/NA_metacycle/meta2d_filtered_rna_NA.txt"
path_rnac<-"data/cycling_analysis/mtx_cyclic_analysis/wol2_rna_pipeline_rpob/"
path_dnac<-"data/cycling_analysis/mgx_cyclic_analysis/wol2_dna_pipeline_rpob/"
dat_path<-"data/cycling_analysis/"
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

jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}
##########################################################
#load MTX cycling files--not rarefied

#annotations
pfam_annot<-fread(dat_annot)

FT_metacyc_MT<-load_cyc_dat(dat_mtx_FTcyc)
sigFTMT<-FT_metacyc_MT%>%filter(JTK_pvalue<0.05,!grepl("DUF",Name)) #662

FA_metacyc_MT<-load_cyc_dat(dat_mtx_FAcyc)
sigFAMT<-FA_metacyc_MT%>%filter(JTK_pvalue<0.05,!grepl("DUF",Name)) #117

NA_metacyc_MT<-load_cyc_dat(dat_mtx_NAcyc)
sigNAMT<-NA_metacyc_MT%>%filter(JTK_pvalue<0.05,!grepl("DUF",Name)) #571
##########################################################
#load MGX cycling files
FT_metacyc_MG<-load_cyc_dat(dat_mgx_FTcyc)
sigFTMG<-FT_metacyc_MG%>%filter(JTK_pvalue<0.05,!grepl("DUF",Name)) #45

FA_metacyc_MG<-load_cyc_dat(dat_mgx_FAcyc)
sigFAMG<-FA_metacyc_MG%>%filter(JTK_pvalue<0.05,!grepl("DUF",Name)) #42

NA_metacyc_MG<-load_cyc_dat(dat_mgx_NAcyc)
sigNAMG<-NA_metacyc_MG%>%filter(JTK_pvalue<0.05,!grepl("DUF",Name)) #308
##########################################################
#load MTX cycling files--rarefied

#annotations
pfam_annot<-fread(dat_annot)

FT_metacyc_144k_MT<-load_cyc_dat(dat_mtx_FTcyc_144k)
sigFTMT_144k<-FT_metacyc_144k_MT%>%filter(JTK_pvalue<0.05,!grepl("DUF",Name)) #122

FA_metacyc_144k_MT<-load_cyc_dat(dat_mtx_FAcyc_144k)
sigFAMT_144k<-FA_metacyc_144k_MT%>%filter(JTK_pvalue<0.05,!grepl("DUF",Name)) #56

NA_metacyc_144k_MT<-load_cyc_dat(dat_mtx_NAcyc_144k)
sigNAMT_144k<-NA_metacyc_144k_MT%>%filter(JTK_pvalue<0.05,!grepl("DUF",Name)) #292

##########################################################
#get jaccard distance MGX vs. MTX (non-rarefied)--Figure 3C
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
ggsave(paste0(fig_path,"SFR23_0619_jaccard_metacy_MTrMG_summ_small.pdf"),plot=p, width = 3, height = 3.5)

##########################################################
#get jaccard distance MGX vs. MTX (rarefied)--Figure S3C
j_val<-c(jaccard(sigNAMG$FeatureID,sigNAMT_144k$FeatureID),jaccard(sigNAMG$FeatureID,sigFAMT_144k$FeatureID),jaccard(sigNAMG$FeatureID,sigFTMT_144k$FeatureID),
         jaccard(sigFAMG$FeatureID,sigFAMT_144k$FeatureID),jaccard(sigFAMG$FeatureID,sigFTMT_144k$FeatureID),
         jaccard(sigFTMG$FeatureID,sigFTMT_144k$FeatureID))


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

ggsave(paste0(fig_path,"SFR23_0620_jaccard_metacy_MT144krMG_summ_small.pdf"),plot=p, width = 3, height = 3.5)

##########################################################
#combine results in one table --Table S3

gonames<-fread(dat_go)
pfam2GO<-read.table(dat_pfamtogo,header = FALSE, sep = "\t",
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

feat_annot<-fread(dat_annot)%>%
  left_join(.,pfam2GO,by="FeatureID")

cNAt<-fread(paste0(path_rnac,"NA_metacycle/meta2d_filtered_rna_NA.txt"))%>%
  mutate(condition="NA",datatype="metaT")%>%
  arrange(JTK_pvalue)
cFAt<-fread(paste0(path_rnac,"FA_metacycle/meta2d_filtered_rna_FA.txt"))%>%
  mutate(condition="FA",datatype="metaT")%>%
  arrange(JTK_pvalue)
cFTt<-fread(paste0(path_rnac,"FT_metacycle/meta2d_filtered_rna_FT.txt"))%>%
  mutate(condition="FT",datatype="metaT")%>%
  arrange(JTK_pvalue)

cNAm<-fread(paste0(path_dnac,"NA_metacycle/meta2d_filtered_dna_NA.txt"))%>%
  mutate(condition="NA",datatype="metaG")%>%
  arrange(JTK_pvalue)
cFAm<-fread(paste0(path_dnac,"FA_metacycle/meta2d_filtered_dna_FA.txt"))%>%
  mutate(condition="FA",datatype="metaG")%>%
  arrange(JTK_pvalue)
cFTm<-fread(paste0(path_dnac,"FT_metacycle/meta2d_filtered_dna_FT.txt"))%>%
  mutate(condition="FT",datatype="metaG")%>%
  arrange(JTK_pvalue)

comcyc<-rbind(cNAt,cFAt,cFTt,cNAm,cFAm,cFTm)%>%
  dplyr::rename(FeatureID=CycID)%>%
  left_join(.,feat_annot,by="FeatureID")%>%
  dplyr::select(FeatureID,Name,GO_Term,condition,datatype,everything())%>%
  arrange(JTK_pvalue)
write.table(comcyc, file = paste0(dat_path,"metacycle_results_combined.txt"), #Table S3
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
