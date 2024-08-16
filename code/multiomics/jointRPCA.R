setwd("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/multiomics/")

library(tidyverse)
library(data.table)
library("qiime2R")
library("Biostrings")
library(ggcorrplot)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(scatterplot3d)
library(VennDiagram)
####################################################################

#clean up mtx file
mtx_md<-fread("~/scratch/TRF_multiomics/metatranscript/woltka2_m_results/metaT_metadata_ztcat_noNT.txt")%>%
  filter(zt_time==5 |zt_time==17 | zt_time==9 | zt_time==21)%>%
  filter(sample_name!="cNA09b"& sample_name!="cNA09c"& sample_name!="cNA21b"& sample_name!="cNA21c")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))

test_list<-c("FA-Day-1","FA-Night-2","FT-Day-2","FT-Night-1","NA-Day-3","NA-Night-2")
newids<-c(1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4)

mtx_md_cln<-mtx_md%>%
  cbind(newids)%>%
  mutate(DayNight=ifelse(phase=="light","Day","Night"),
         new_samp_id=paste(condition,DayNight,newids,sep="-"),
         traintest=ifelse(new_samp_id %in% test_list, "test","train"))
write.table(mtx_md_cln, "jrpca/jrpca_mtx_key.txt",sep = "\t",row.names = FALSE, quote=FALSE)

jrpca_md<-mtx_md_cln%>%
  dplyr::select(new_samp_id,organ, diet,schedule,condition,phase,traintest)%>%
  dplyr::rename(sample_name=new_samp_id)

write.table(jrpca_md, "jrpca/jrpca_metadata.txt",sep = "\t",row.names = FALSE, quote=FALSE)

mtx<-fread("~/scratch/TRF_multiomics/metatranscript/woltka2_m_results/pfam/pfam_clean_noNT.tsv")%>%
  dplyr::select(FeatureID,all_of(mtx_md$sample_name))

names(mtx) <- c("FeatureID","FA-Day-1","FA-Day-2","FA-Day-3","FA-Day-4","FA-Night-1","FA-Night-2","FA-Night-3","FA-Night-4","FT-Day-1","FT-Day-2", 
                "FT-Day-3","FT-Day-4","FT-Night-1","FT-Night-2","FT-Night-3","FT-Night-4","NA-Day-1","NA-Day-2","NA-Day-3","NA-Day-4","NA-Night-1",
                "NA-Night-2","NA-Night-3","NA-Night-4")
write.table(mtx, "jrpca/jrpca_mtx_feattab.txt",sep = "\t",row.names = FALSE, quote=FALSE)

#make mtx feature key
pfam_annot<-fread("~/scratch/TRF_multiomics/metatranscript/woltka2_m_results/pfam/pfam_annotationkey.csv")%>%
  dplyr::rename(pfam_name=Name)
gonames<-fread("~/scratch/TRF_multiomics/metatranscript/woltka2_results/go_name.txt")%>%
  dplyr::rename(GO_name=name)
pfam2go<-read.table("~/scratch/TRF_multiomics/metatranscript/woltka2_results/pfam-to-go-process.map",header = FALSE, sep = "\t",
                      col.names = paste0("V",seq_len(4)), fill = TRUE)%>%
  gather(column,GO_Term,-V1)%>%
  dplyr::select(1,3)%>%
  dplyr::rename(FeatureID=V1)%>%
  distinct(FeatureID,.keep_all = TRUE)

mtx_annot<-pfam_annot%>%
  left_join(.,pfam2go, by="FeatureID")%>%
  left_join(.,gonames,by="GO_Term")
write.table(mtx_annot, "jrpca/jrpca_mtx_annotation_key.txt",sep = "\t",row.names = FALSE, quote=FALSE)

#clean up mtg file
mtg_md<-fread("~/scratch/TRF_multiomics/metagenomic/metaG_metadata_noNT.txt")%>%
  filter(zt_time==5 |zt_time==17 | zt_time==9 | zt_time==21)%>%
  filter(sample_name!="NA6b"& sample_name!="NA6c"& sample_name!="NA3b"& sample_name!="NA3c")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  dplyr::rename(phase=lightdark)%>%
  arrange(zt_time)%>%
  arrange(condition)

test_list<-c("FA-Day-1","FA-Night-2","FT-Day-2","FT-Night-1","NA-Day-3","NA-Night-2")
newids<-c(1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4)

mtg_md_cln<-mtg_md%>%
  cbind(newids)%>%
  mutate(DayNight=ifelse(phase=="light","Day","Night"),
         new_samp_id=paste(condition,DayNight,newids,sep="-"),
         traintest=ifelse(new_samp_id %in% test_list, "test","train"))
write.table(mtg_md_cln, "jrpca_mtg/jrpca_mtg_key.txt",sep = "\t",row.names = FALSE, quote=FALSE)

jrpca_md<-mtg_md_cln%>%
  dplyr::select(new_samp_id, diet,schedule,condition,phase,traintest)%>%
  dplyr::rename(sample_name=new_samp_id)

write.table(jrpca_md, "jrpca_mtg/jrpca_metadata.txt",sep = "\t",row.names = FALSE, quote=FALSE)

mtg<-fread("~/scratch/TRF_multiomics/metagenomic/woltka2_results/filtered_metaG/pfam_notnorm/pfam_clean_noNT.txt")%>%
  dplyr::select(FeatureID,all_of(mtg_md$sample_name))

names(mtg) <- c("FeatureID","FA-Day-1","FA-Day-2","FA-Day-3","FA-Day-4","FA-Night-1","FA-Night-2","FA-Night-3","FA-Night-4","FT-Day-1","FT-Day-2", 
                "FT-Day-3","FT-Day-4","FT-Night-1","FT-Night-2","FT-Night-3","FT-Night-4","NA-Day-1","NA-Day-2","NA-Day-3","NA-Day-4","NA-Night-1",
                "NA-Night-2","NA-Night-3","NA-Night-4")
write.table(mtg, "jrpca_mtg/jrpca_mtg_feattab.txt",sep = "\t",row.names = FALSE, quote=FALSE)

#make mtg feature key
  pfam_annot<-fread("~/scratch/TRF_multiomics/metagenomic/woltka2_results/filtered_metaG/pfam_notnorm/pfam_annotationkey.csv")%>%
  dplyr::rename(pfam_name=Name)
gonames<-fread("~/scratch/TRF_multiomics/metatranscript/woltka2_results/go_name.txt")%>%
  dplyr::rename(GO_name=name)
pfam2go<-read.table("~/scratch/TRF_multiomics/metatranscript/woltka2_results/pfam-to-go-process.map",header = FALSE, sep = "\t",
                    col.names = paste0("V",seq_len(4)), fill = TRUE)%>%
  gather(column,GO_Term,-V1)%>%
  dplyr::select(1,3)%>%
  dplyr::rename(FeatureID=V1)%>%
  distinct(FeatureID,.keep_all = TRUE)

mtg_annot<-pfam_annot%>%
  left_join(.,pfam2go, by="FeatureID")%>%
  left_join(.,gonames,by="GO_Term")
write.table(mtg_annot, "jrpca_mtg/jrpca_mtg_annotation_key.txt",sep = "\t",row.names = FALSE, quote=FALSE)

#clean up mtb file
mtb_md<-fread("~/scratch/TRF_multiomics/metabolom/metabolom_metadata_clean.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  mutate(new_samp_id=paste(GROUP,replicate,sep="-"),
         traintest=ifelse(new_samp_id %in% test_list, "test","train"))%>%
  arrange(condition)%>%
  filter(replicate!=5)
#write.table(mtb_md, "jrpca/jrpca_mtb_key.txt",sep = "\t",row.names = FALSE, quote=FALSE)
write.table(mtb_md, "jrpca_mtg/jrpca_mtb_key.txt",sep = "\t",row.names = FALSE, quote=FALSE)

mtb<-fread("~/scratch/TRF_multiomics/metabolom/NightDay_metabolom_cln.txt")%>%
  dplyr::select(FeatureID,all_of(mtb_md$sample_name))

names(mtb) <- c("FeatureID","FA-Day-1","FA-Day-2","FA-Day-3","FA-Day-4","FA-Night-1","FA-Night-2","FA-Night-3","FA-Night-4","FT-Day-1","FT-Day-2", 
                "FT-Day-3","FT-Day-4","FT-Night-1","FT-Night-2","FT-Night-3","FT-Night-4","NA-Day-1","NA-Day-2","NA-Day-3","NA-Day-4","NA-Night-1",
                "NA-Night-2","NA-Night-3","NA-Night-4")
#write.table(mtb, "jrpca/jrpca_mtb_feattab.txt",sep = "\t",row.names = FALSE, quote=FALSE)
write.table(mtb, "jrpca_mtg/jrpca_mtb_feattab.txt",sep = "\t",row.names = FALSE, quote=FALSE)

#make mtb feature key
mtb_key<-fread("~/scratch/TRF_multiomics/metabolom/UCSD-03-13VW_OrigScale_mod.csv")%>%
  dplyr::select(c(5,1:4,6:14))%>% dplyr::rename(FeatureID=COMP_ID)
#write.table(mtb_key, "jrpca/jrpca_mtb_annotation_key.txt",sep = "\t",row.names = FALSE, quote=FALSE)
write.table(mtb_key, "jrpca_mtg/jrpca_mtb_annotation_key.txt",sep = "\t",row.names = FALSE, quote=FALSE)
####################################################################

#create a table with the selected features 5% mtx and carb/scfa mtb hits

mtx_feat<-fread("jrpca/qurro_mtx/selected_features_mtx5per.tsv")

mtx<-fread("jrpca/jrpca_mtx_feattab.txt")%>%
  filter(FeatureID %in% mtx_feat$`Feature ID`)
write.table(mtx, "jrpca/jrpca_mtx5per_feattab.txt",sep = "\t",row.names = FALSE, quote=FALSE)

mtb_feat<-fread("jrpca/qurro_mtb/selected_features_mtbcarb_LCFA.tsv")

mtb<-fread("jrpca/jrpca_mtb_feattab.txt")%>%
  filter(FeatureID %in% mtb_feat$`Feature ID`)
write.table(mtb, "jrpca/jrpca_mtbcarbSCFA_feattab.txt",sep = "\t",row.names = FALSE, quote=FALSE)

corr<-fread("jrpca/correlation_qurro_table/Correlation.tsv")%>%
  dplyr::select(c(1:48))
corr<-corr[c(48:707),]#PF00144.27

write.table(corr, "jrpca/correlation_qurro_table/Correlation_sub.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

#create a table with the selected features from aldex2  from top/bottom5%

mtx_feat1<-fread("jrpca/qurro_mtx/selected_features_mtx5per.tsv")
mtx_feat2<-fread("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_m_results/pfam/aldex/SFR23_0620_FANA_450uniquehits_BHless0.1.txt")
mtx<-fread("jrpca/jrpca_mtx_feattab.txt")%>%
  filter(FeatureID %in% mtx_feat1$`Feature ID`)%>%
  filter(FeatureID %in% mtx_feat2$FeatureID)
write.table(mtx, "jrpca/jrpca_mtxFANA_450aldex_feattab.txt",sep = "\t",row.names = FALSE, quote=FALSE)

corr<-fread("jrpca/correlation_qurro4_table/Correlation.tsv")%>%
  dplyr::select(c(1:48))
corr<-corr[c(48:99),]#PF01926.26

write.table(corr, "jrpca/correlation_qurro4_table/Correlation_sub.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

#create a table with the selected features from aldex2 

mtx_feat<-fread("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_m_results/pfam/aldex/SFR23_0620_FAFTD_73uniquehits_BHless0.1.txt")

mtx<-fread("jrpca/jrpca_mtx_feattab.txt")%>%
  filter(FeatureID %in% mtx_feat$FeatureID)
write.table(mtx, "jrpca/jrpca_mtxFAFTD73aldex_feattab.txt",sep = "\t",row.names = FALSE, quote=FALSE)

mtx_feat1<-fread("jrpca/qurro_mtx/selected_features_mtx5per.tsv")
mtx_feat2<-fread("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_m_results/pfam/aldex/SFR23_0620_FAFTD_73uniquehits_BHless0.1.txt")
mtx<-fread("jrpca/jrpca_mtx_feattab.txt")%>%
  filter(FeatureID %in% mtx_feat1$`Feature ID`)%>%
  filter(FeatureID %in% mtx_feat2$FeatureID)
write.table(mtg, "jrpca/jrpca_mtxFAFTD_73aldex_feattab_5peroverlap.txt",sep = "\t",row.names = FALSE, quote=FALSE)
#only 4 feat overlap between the 5perc and FAFTD

corr<-fread("jrpca/correlation_qurro3_table/Correlation.tsv")%>%
  dplyr::select(c(1:48))
corr<-corr[c(48:120),]#PF13561.9

write.table(corr, "jrpca/correlation_qurro3_table/Correlation_sub.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

#create a table with the selected features carb/fattyacid mtx and carb/scfa mtb hits

mtx_feat<-fread("jrpca/qurro_mtx/selected_features_mtxcarb_fattyacid.tsv")

mtx<-fread("jrpca/jrpca_mtx_feattab.txt")%>%
  filter(FeatureID %in% mtx_feat$`Feature ID`)
write.table(mtx, "jrpca/jrpca_mtxcarb_fattyacid_feattab.txt",sep = "\t",row.names = FALSE, quote=FALSE)

corr<-fread("jrpca/correlation_qurro2_table/Correlation.tsv")%>%
  dplyr::select(c(1:48))
corr<-corr[c(48:122),]#PF03033.23

write.table(corr, "jrpca/correlation_qurro2_table/Correlation_sub.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

mtx_key<-fread("jrpca/jrpca_mtx_annotation_key.txt")%>%
  dplyr::select(FeatureID,pfam_name)%>%
  dplyr::rename(name=pfam_name)%>%
  mutate(omics="MTX")
mtb_key<-fread("jrpca/jrpca_mtb_annotation_key.txt")%>%
  dplyr::select(FeatureID,BIOCHEMICAL)%>%
  dplyr::rename(name=BIOCHEMICAL)%>%
  mutate(omics="metab")

key<-rbind(mtx_key,mtb_key)%>%
  mutate(annot=paste(FeatureID,name,sep=" "))%>%
  dplyr::select(-name)

corr_m<-fread("jrpca/correlation_qurro2_table/Correlation.tsv")%>%
  dplyr::rename(FeatureID=featureid)%>%
  left_join(.,key,by="FeatureID")%>%
  dplyr::select(c(125,2:123))%>%
  column_to_rownames("annot")%>%t()%>%as.data.frame()%>%
  rownames_to_column("FeatureID")%>%
  left_join(.,key,by="FeatureID")%>%
  dplyr::select(c(125,2:123))%>%column_to_rownames("annot")
 
key_m<-key%>%dplyr::select(annot,omics)%>%
  mutate(omics=factor(omics,levels=c("MTX","metab")))%>%
           column_to_rownames("annot")

my_colour = list(omics = c(MTX = "blue", metab = "green4"))

p<-pheatmap(corr_m, annotation_row = key_m, annotation_col = key_m, annotation_colors = my_colour,color= inferno(10))
ggsave("jrpca/correlation_qurro2_table/SFR23_0628_heatmap_wannot.pdf", p,width = 30, height = 30)

#corr_n<-fread("jrpca/correlation_qurro2_table/Correlation_sub.tsv",header=TRUE)%>%
#corr_n<-fread("jrpca/correlation_qurro3_table/Correlation_sub.tsv",header=TRUE)%>%
corr_n<-fread("jrpca/correlation_qurro4_table/Correlation_sub.tsv",header=TRUE)%>%
  dplyr::rename(FeatureID=featureid)%>%
  left_join(.,key,by="FeatureID")%>%
  dplyr::select(c(50,2:48))%>%
  column_to_rownames("annot")%>%t()%>%as.data.frame()%>%
  rownames_to_column("FeatureID")%>%
  left_join(.,key,by="FeatureID")%>%
  #dplyr::select(c(78,2:76))%>%
  #dplyr::select(c(76,2:74))%>%
  dplyr::select(c(55,2:53))%>%
  column_to_rownames("annot")

#selfet_mtx<-fread("jrpca/qurro_mtx/selected_features_mtxcarbmtbproc_fattyacid.tsv")
#selfet_mtx<-fread("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_m_results/pfam/aldex/SFR23_0620_FAFTD_73uniquehits_BHless0.1.txt")
selfet_mtx<-fread("jrpca/jrpca_mtxFANA_450aldex_feattab.txt")
mtx_key<-fread("jrpca/jrpca_mtx_annotation_key.txt")%>%
  #filter(FeatureID %in% selfet_mtx$`Feature ID`)%>%
  filter(FeatureID %in% selfet_mtx$FeatureID)%>%
  mutate(annot=paste(FeatureID,pfam_name,sep=" "))%>%
  dplyr::rename(GO_process=GO_name)%>%
  dplyr::select(annot,GO_process)%>%
  column_to_rownames("annot")%>%
  mutate(GO_process=ifelse(is.na(GO_process),"unknown",GO_process))

selfet_mtb<-fread("jrpca/qurro_mtb/selected_features_mtbcarb_LCFA.tsv")
mtb_key<-fread("jrpca/jrpca_mtb_annotation_key.txt")%>%
  filter(FeatureID %in% selfet_mtb$`Feature ID`)%>%
  mutate(annot=paste(FeatureID,BIOCHEMICAL,sep=" "))%>%
  dplyr::rename(pathway=SUPER_PATHWAY)%>%
  dplyr::select(annot,pathway)%>%
  column_to_rownames("annot")

# my_colour = list(pathway = c(MTX = "blue", metab = "green4"),
#                  GO_process=c(Lipid=,))annotation_colors = my_colour

draw_colnames_45 <- function (coln, gaps, ...) {
  coord <- pheatmap:::find_coordinates(length(coln), gaps)
  x     <- coord$coord - 0.5 * coord$size
  res   <- grid::textGrob(
    coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"),
    vjust = 0.75, hjust = 1, rot = 90, gp = grid::gpar(...)
  )
  return(res)
}

assignInNamespace(
  x = "draw_colnames",
  value = "draw_colnames_45",
  ns = asNamespace("pheatmap")
)

p<-pheatmap(corr_n, border_color= NA,annotation_row = mtb_key, annotation_col = mtx_key,color= inferno(10))
#ggsave("jrpca/correlation_qurro2_table/SFR23_0628_heatmap_wbio_annot.pdf", p,width = 20, height = 13)
#ggsave("jrpca/correlation_qurro3_table/SFR23_0703_heatmap_FAFTDmtx_wbio_annot.pdf", p,width = 20, height = 13)
ggsave("jrpca/correlation_qurro4_table/SFR23_0703_heatmap_FANAmtx_wbio_annot.pdf", p,width = 16, height = 13)



# corr<-fread("jrpca/correlation_table/Correlation.tsv")%>%column_to_rownames("featureid")
# ggcorrplot(corr, hc.order = TRUE)
# 
# corr<-fread("jrpca/correlation_mtb_table/Correlation.tsv")%>%column_to_rownames("featureid")
# ggcorrplot(corr, hc.order = TRUE)
# 
# corr<-fread("jrpca/correlation_mtx_table/Correlation.tsv")%>%column_to_rownames("featureid")
# ggcorrplot(corr, hc.order = TRUE)

#######################################################################################

#jrpca-mtx and mtb
ord <- read_qza("jrpca/joint_biplot.qza")

samp_ord<-ord$data$Vectors
write.table(samp_ord,"jrpca/joint_biplot_sample_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)
feat_ord<-ord$data$Species
write.table(feat_ord,"jrpca/joint_biplot_feature_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)

md<-fread("jrpca/jrpca_metadata.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))

rpca<-ord$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2, PC3)%>%
  dplyr::rename(sample_name=SampleID)%>%
  left_join(md,by="sample_name")%>%
  mutate(condition=factor(condition,levels=c("NA","FA","FT")),
         phase=factor(phase,levels=c("light","dark")))

pdf("jrpca/SFR23_0628_jrpcabiplot.pdf",width = 5, height = 5)

# Create a 3D PCA plot
colors <- c("#0072B2","#D55E00","#009E73")
colors <- colors[as.numeric(rpca$condition)]

scatterplot3d(rpca[,2:4],angle = 300,pch = 16,color=colors,grid=FALSE,
              xlab = "PC1", ylab = "PC2", zlab = "PC3",
              main = "MTX + Targeted Metabolites")
legend("bottom", legend = levels(rpca$condition),
       col =  c("#0072B2","#D55E00","#009E73"), pch = 16, 
       inset = -0.25, xpd = TRUE, horiz = TRUE)

dev.off()

#jrpca-mtg and mtb
ord <- read_qza("jrpca_mtg/joint_biplot.qza")

samp_ord<-ord$data$Vectors
write.table(samp_ord,"jrpca_mtg/joint_biplot_sample_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)
feat_ord<-ord$data$Species
write.table(feat_ord,"jrpca_mtg/joint_biplot_feature_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)

md<-fread("jrpca_mtg/jrpca_metadata.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))

rpca<-ord$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2, PC3)%>%
  dplyr::rename(sample_name=SampleID)%>%
  left_join(md,by="sample_name")%>%
  mutate(condition=factor(condition,levels=c("NA","FA","FT")),
         phase=factor(phase,levels=c("light","dark")))

pdf("jrpca_mtg/SFR23_0725_jrpcabiplot.pdf",width = 5, height = 5)

# Create a 3D PCA plot
colors <- c("#0072B2","#D55E00","#009E73")
colors <- colors[as.numeric(rpca$condition)]

scatterplot3d(rpca[,2:4],angle = 300,pch = 16,color=colors,grid=FALSE,
              xlab = "PC1", ylab = "PC2", zlab = "PC3",
              main = "MTG + Targeted Metabolites")
legend("bottom", legend = levels(rpca$condition),
       col =  c("#0072B2","#D55E00","#009E73"), pch = 16, 
       inset = -0.25, xpd = TRUE, horiz = TRUE)

dev.off()

#######################################################################################

#plot qurro nat log ratios-jrcpa mtx and mtb comparison

md<-fread("jrpca/jrpca_metadata.txt")%>%mutate(condition=ifelse(is.na(condition),"NA",condition))

#mtb

#mtb carbohydrate vs long chain fatty acid
#qr1<-fread("jrpca/qurro_mtb/sample_plot_data_mtbcarb_LCFA.tsv")%>% #same p-values for mtx-mtb and mtg-mtb
qr1<-fread("jrpca_mtg/qurro_mtb/sample_plot_data_carb_LCFA.tsv")%>%
  dplyr::select(c(1:2))%>%dplyr::rename(sample_name=`Sample ID`)%>%
  left_join(.,md,by="sample_name")%>%
  mutate(condition=factor(condition,levels=c("NA","FA","FT")),
         phase=factor(phase,levels=c("light","dark")))


p<-ggplot(qr1, aes(x=condition, y=Current_Natural_Log_Ratio,fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  theme_minimal()+ scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
  labs(x="condition",y="Natural Log Ratio (Carbohydrate/LCFA)",title="targeted metabolomics")+
  theme(legend.position = "none")

ggsave("jrpca/qurro_mtb/natlog_mtbcarb_LCFA.pdf", width = 3, height = 3)

pairwise.wilcox.test(qr1$Current_Natural_Log_Ratio, qr1$condition,
                     p.adjust.method="fdr")

# NA      FA     
# FA 0.00023 -      
#   FT 0.00023 0.03792

p<-ggplot(qr1, aes(x=condition, y=Current_Natural_Log_Ratio,fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  facet_wrap(~phase)+
  theme_classic()+ scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
  labs(x="condition",y="Natural Log Ratio (Carbohydrate/LCFA)",title="")+
  theme(legend.position = "none")
ggsave("jrpca/qurro_mtb/natlog_mtbcarb_LCFA_lightdark.pdf", plot=p,height=3, width=4)

qr1L<-qr1%>%filter(phase=="light")
pairwise.wilcox.test(qr1L$Current_Natural_Log_Ratio, qr1L$condition,
                     p.adjust.method="fdr")

# NA    FA   
# FA 0.043 -    
#   FT 0.043 0.343

qr1D<-qr1%>%filter(phase=="dark")
pairwise.wilcox.test(qr1D$Current_Natural_Log_Ratio, qr1D$condition,
                     p.adjust.method="fdr")

# NA    FA   
# FA 0.043 -    
#   FT 0.043 0.114

qr1FT<-qr1%>%filter(condition=="FT")
pairwise.wilcox.test(qr1FT$Current_Natural_Log_Ratio, qr1FT$phase,
                     p.adjust.method="fdr")

# light
# dark 0.49

qr1FA<-qr1%>%filter(condition=="FA")
pairwise.wilcox.test(qr1FA$Current_Natural_Log_Ratio, qr1FA$phase,
                     p.adjust.method="fdr")

# light
# dark 0.11

qr1NA<-qr1%>%filter(condition=="NA")
pairwise.wilcox.test(qr1NA$Current_Natural_Log_Ratio, qr1NA$phase,
                     p.adjust.method="fdr")
# light
# dark 0.69 

#plot 2BA to 1BA

#qr2<-fread("jrpca/qurro_mtb/sample_plot_data_mtb2BA_1BA.tsv")%>% #same p-values for mtx-mtb and mtg-mtb
qr2<-fread("jrpca_mtg/qurro_mtb/sample_plot_data_2BA_1BA.tsv")%>%
  dplyr::select(c(1:2))%>%dplyr::rename(sample_name=`Sample ID`)%>%
  left_join(.,md,by="sample_name")%>%
  mutate(condition=factor(condition,levels=c("NA","FA","FT")),
         phase=factor(phase,levels=c("light","dark")))


p<-ggplot(qr2, aes(x=condition, y=Current_Natural_Log_Ratio,fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  theme_minimal()+ scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
  labs(x="condition",y="Natural Log Ratio (2BA/1BA)",title="targeted metabolomics")+
  theme(legend.position = "none")

ggsave("jrpca/qurro_mtb/natlog_mtb2BA_1BA.pdf", width = 3, height = 3)

pairwise.wilcox.test(qr2$Current_Natural_Log_Ratio, qr2$condition,
                     p.adjust.method="fdr")

# NA     FA    
# FA 0.0748 -     
#   FT 0.0089 0.3282

p<-ggplot(qr2, aes(x=condition, y=Current_Natural_Log_Ratio,fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  facet_wrap(~phase)+
  theme_classic()+ scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
  labs(x="condition",y="Natural Log Ratio (Carbohydrate/LCFA)",title="")+
  theme(legend.position = "none")
ggsave("jrpca/qurro_mtb/natlog_mtb2BA_1BA_lightdark.pdf", plot=p,height=3, width=4)

qr2L<-qr2%>%filter(phase=="light")
pairwise.wilcox.test(qr2L$Current_Natural_Log_Ratio, qr2L$condition,
                     p.adjust.method="fdr")

# NA    FA   
# FA 0.343 -    
#   FT 0.086 0.171

qr2D<-qr2%>%filter(phase=="dark")
pairwise.wilcox.test(qr2D$Current_Natural_Log_Ratio, qr2D$condition,
                     p.adjust.method="fdr")

# NA   FA  
# FA 0.17 -   
#   FT 0.17 1.00

qr2FT<-qr2%>%filter(condition=="FT")
pairwise.wilcox.test(qr2FT$Current_Natural_Log_Ratio, qr2FT$phase,
                     p.adjust.method="fdr")

# light
# dark 0.34

qr2FA<-qr2%>%filter(condition=="FA")
pairwise.wilcox.test(qr2FA$Current_Natural_Log_Ratio, qr2FA$phase,
                     p.adjust.method="fdr")

# light
# dark 0.029

qr2NA<-qr2%>%filter(condition=="NA")
pairwise.wilcox.test(qr2NA$Current_Natural_Log_Ratio, qr2NA$phase,
                     p.adjust.method="fdr")
# light
# dark 0.89

#######################################################################################

#mtx qurro plot

qr3<-fread("jrpca/qurro_mtx/sample_plot_data_mtxcarbmtbproc_fattyacid.tsv")%>%
  dplyr::select(c(1:2))%>%dplyr::rename(sample_name=`Sample ID`)%>%
  left_join(.,md,by="sample_name")%>%
  mutate(condition=factor(condition,levels=c("NA","FA","FT")),
         phase=factor(phase,levels=c("light","dark")))


p<-ggplot(qr3, aes(x=condition, y=Current_Natural_Log_Ratio,fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  theme_minimal()+ scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
  labs(x="condition",y="Natural Log Ratio\n(carbohydrate metab proc/fatty acid metab)",title="MTX")+
  theme(legend.position = "none")

ggsave("jrpca/qurro_mtx/natlog_mtxcarbmtbproc_fattyacid.pdf", width = 3, height = 3)

pairwise.wilcox.test(qr3$Current_Natural_Log_Ratio, qr3$condition,
                     p.adjust.method="fdr")
# NA   FA  
# FA 0.88 -   
#   FT 0.49 0.49

p<-ggplot(qr3, aes(x=condition, y=Current_Natural_Log_Ratio,fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  facet_wrap(~phase)+
  theme_classic()+ scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
  labs(x="condition",y="Natural Log Ratio\n(carbohydrate metab proc/fatty acid metab)",title="")+
  theme(legend.position = "none")
ggsave("jrpca/qurro_mtx/natlog_mtxcarbmtbproc_fattyacid_lightdark.pdf", plot=p,height=3, width=4)

qr3L<-qr3%>%filter(phase=="light")
pairwise.wilcox.test(qr3L$Current_Natural_Log_Ratio, qr3L$condition,
                     p.adjust.method="fdr")

# NA   FA  
# FA 0.34 -   
#   FT 0.51 1.00

qr3D<-qr3%>%filter(phase=="dark")
pairwise.wilcox.test(qr3D$Current_Natural_Log_Ratio, qr3D$condition,
                     p.adjust.method="fdr")

# NA   FA  
# FA 0.51 -   
#   FT 0.69 0.34

qr3FT<-qr3%>%filter(condition=="FT")
pairwise.wilcox.test(qr3FT$Current_Natural_Log_Ratio, qr3FT$phase,
                     p.adjust.method="fdr")

# light
# dark 0.11

qr3FA<-qr3%>%filter(condition=="FA")
pairwise.wilcox.test(qr3FA$Current_Natural_Log_Ratio, qr3FA$phase,
                     p.adjust.method="fdr")

# light
# dark 1

qr3NA<-qr3%>%filter(condition=="NA")
pairwise.wilcox.test(qr3NA$Current_Natural_Log_Ratio, qr3NA$phase,
                     p.adjust.method="fdr")
# light
# dark 0.029

qr4<-fread("jrpca/qurro_mtx/sample_plot_data_mtxfattyacid_transcript.tsv")%>%
  dplyr::select(c(1:2))%>%dplyr::rename(sample_name=`Sample ID`)%>%
  left_join(.,md,by="sample_name")%>%
  mutate(condition=factor(condition,levels=c("NA","FA","FT")),
         phase=factor(phase,levels=c("light","dark")))

p<-ggplot(qr4, aes(x=condition, y=Current_Natural_Log_Ratio,fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  theme_minimal()+ scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
  labs(x="condition",y="Natural Log Ratio\n(fatty acid metab/transcription)",title="MTX")+
  theme(legend.position = "none")

ggsave("jrpca/qurro_mtx/natlog_mtxcarbmtbproc_transcript.pdf", width = 3, height = 3)

pairwise.wilcox.test(qr4$Current_Natural_Log_Ratio, qr4$condition,
                     p.adjust.method="fdr")
# NA   FA  
# FA 0.88 -   
#   FT 0.16 0.16

qr5<-fread("jrpca/qurro_mtx/sample_plot_data_mtxcarbmtbproc_transcript.tsv")%>%
  dplyr::select(c(1:2))%>%dplyr::rename(sample_name=`Sample ID`)%>%
  left_join(.,md,by="sample_name")%>%
  mutate(condition=factor(condition,levels=c("NA","FA","FT")),
         phase=factor(phase,levels=c("light","dark")))

p<-ggplot(qr5, aes(x=condition, y=Current_Natural_Log_Ratio,fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  theme_minimal()+ scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
  labs(x="condition",y="Natural Log Ratio\n(carbohydrate metab proc/transcription)",title="MTX")+
  theme(legend.position = "none")

ggsave("jrpca/qurro_mtx/natlog_mtxcarbfattyacid_transcript.pdf", width = 3, height = 3)

pairwise.wilcox.test(qr5$Current_Natural_Log_Ratio, qr5$condition,
                     p.adjust.method="fdr")

# NA    FA   
# FA 0.798 -    
#   FT 0.057 0.057


qr6<-fread("jrpca/qurro_mtx/sample_plot_data_mtxprotelysis_transcript.tsv")%>%
  dplyr::select(c(1:2))%>%dplyr::rename(sample_name=`Sample ID`)%>%
  left_join(.,md,by="sample_name")%>%
  mutate(condition=factor(condition,levels=c("NA","FA","FT")),
         phase=factor(phase,levels=c("light","dark")))

p<-ggplot(qr6, aes(x=condition, y=Current_Natural_Log_Ratio,fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  theme_minimal()+ scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
  labs(x="condition",y="Natural Log Ratio (proteolysis/transcription)",title="MTX")+
  theme(legend.position = "none")

ggsave("jrpca/qurro_mtx/natlog_mtxproteolysis_transcript.pdf", width = 3, height = 3)

pairwise.wilcox.test(qr6$Current_Natural_Log_Ratio, qr6$condition,
                     p.adjust.method="fdr")


qr7<-fread("jrpca/qurro_mtx/sample_plot_data_mtx5per.tsv")%>%
  dplyr::select(c(1:2))%>%dplyr::rename(sample_name=`Sample ID`)%>%
  left_join(.,md,by="sample_name")%>%
  mutate(condition=factor(condition,levels=c("NA","FA","FT")),
         phase=factor(phase,levels=c("light","dark")))

p<-ggplot(qr7, aes(x=condition, y=Current_Natural_Log_Ratio,fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  theme_minimal()+ scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
  labs(x="condition",y="Natural Log Ratio (top5/bottom5)",title="MTX")+
  theme(legend.position = "none")

ggsave("jrpca/qurro_mtx/natlog_mtx5per.pdf", width = 3, height = 3)

pairwise.wilcox.test(qr7$Current_Natural_Log_Ratio, qr7$condition,
                     p.adjust.method="fdr")

# NA      FA     
# FA 0.00093 -      
#   FT 0.00047 0.27863

p<-ggplot(qr7, aes(x=condition, y=Current_Natural_Log_Ratio,fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  facet_wrap(~phase)+
  theme_classic()+ scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
  labs(x="condition",y="Natural Log Ratio (top5/bottom5)",title="")+
  theme(legend.position = "none")

ggsave("jrpca/qurro_mtx/natlog_mtx5per_lightdark.pdf", width = 4, height = 3)

qr7L<-qr7%>%filter(phase=="light")
pairwise.wilcox.test(qr7L$Current_Natural_Log_Ratio, qr7L$condition,
                     p.adjust.method="fdr")

# NA    FA   
# FA 0.043 -    
#   FT 0.043 0.686

qr7D<-qr7%>%filter(phase=="dark")
pairwise.wilcox.test(qr7D$Current_Natural_Log_Ratio, qr7D$condition,
                     p.adjust.method="fdr")

# NA    FA   
# FA 0.029 -    
#   FT 0.029 0.029

#######################################################################################

#mtg qurro plot

qr3<-fread("jrpca_mtg/qurro_mtg/sample_plot_data_mtgtop5perc.tsv")%>%
  dplyr::select(c(1:2))%>%dplyr::rename(sample_name=`Sample ID`)%>%
  left_join(.,md,by="sample_name")%>%
  mutate(condition=factor(condition,levels=c("NA","FA","FT")),
         phase=factor(phase,levels=c("light","dark")))


p<-ggplot(qr3, aes(x=condition, y=Current_Natural_Log_Ratio,fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  theme_minimal()+ scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
  labs(x="condition",y="Natural Log Ratio(topbottom5perc)",title="MTG")+
  theme(legend.position = "none")

ggsave("jrpca_mtg/qurro_mtg/natlog_mtgtopbottom5perc.pdf", width = 3, height = 3)

pairwise.wilcox.test(qr3$Current_Natural_Log_Ratio, qr3$condition,
                     p.adjust.method="fdr")
# NA      FA     
# FA 0.00023 -      
#   FT 0.00023 0.72090

p<-ggplot(qr3, aes(x=condition, y=Current_Natural_Log_Ratio,fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  facet_wrap(~phase)+
  theme_classic()+ scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
  labs(x="condition",y="Natural Log Ratio(topbottom5perc)",title="MTG")+
  theme(legend.position = "none")
ggsave("jrpca_mtg/qurro_mtg/natlog_mtgtopbottom5perc_lightdark.pdf", plot=p,height=3, width=4)

qr3L<-qr3%>%filter(phase=="light")
pairwise.wilcox.test(qr3L$Current_Natural_Log_Ratio, qr3L$condition,
                     p.adjust.method="fdr")

# NA    FA   
# FA 0.043 -    
#   FT 0.043 0.114

qr3D<-qr3%>%filter(phase=="dark")
pairwise.wilcox.test(qr3D$Current_Natural_Log_Ratio, qr3D$condition,
                     p.adjust.method="fdr")

# NA    FA   
# FA 0.043 -    
#   FT 0.043 0.486

qr3FT<-qr3%>%filter(condition=="FT")
pairwise.wilcox.test(qr3FT$Current_Natural_Log_Ratio, qr3FT$phase,
                     p.adjust.method="fdr")

# light
# dark 0.2  

qr3FA<-qr3%>%filter(condition=="FA")
pairwise.wilcox.test(qr3FA$Current_Natural_Log_Ratio, qr3FA$phase,
                     p.adjust.method="fdr")

# light
# dark 0.34 

qr3NA<-qr3%>%filter(condition=="NA")
pairwise.wilcox.test(qr3NA$Current_Natural_Log_Ratio, qr3NA$phase,
                     p.adjust.method="fdr")
# light
# dark 0.029

qr4<-fread("jrpca_mtg/qurro_mtg/sample_plot_data_mtgtop10perc.tsv")%>%
  dplyr::select(c(1:2))%>%dplyr::rename(sample_name=`Sample ID`)%>%
  left_join(.,md,by="sample_name")%>%
  mutate(condition=factor(condition,levels=c("NA","FA","FT")),
         phase=factor(phase,levels=c("light","dark")))

p<-ggplot(qr4, aes(x=condition, y=Current_Natural_Log_Ratio,fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  theme_minimal()+ scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
  labs(x="condition",y="Natural Log Ratio(topbottom10perc)",title="MTG")+
  theme(legend.position = "none")

ggsave("jrpca_mtg/qurro_mtg/natlog_mtgtopbottom10perc.pdf", width = 3, height = 3)

pairwise.wilcox.test(qr4$Current_Natural_Log_Ratio, qr4$condition,
                     p.adjust.method="fdr")
# NA      FA     
# FA 0.00047 -      
#   FT 0.00047 0.79845

p<-ggplot(qr4, aes(x=condition, y=Current_Natural_Log_Ratio,fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  facet_wrap(~phase)+
  theme_classic()+ scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
  labs(x="condition",y="Natural Log Ratio(topbottom10perc)",title="MTG")+
  theme(legend.position = "none")
ggsave("jrpca_mtg/qurro_mtg/natlog_mtgtopbottom10perc_lightdark.pdf", plot=p,height=3, width=4)

qr4L<-qr4%>%filter(phase=="light")
pairwise.wilcox.test(qr4L$Current_Natural_Log_Ratio, qr4L$condition,
                     p.adjust.method="fdr")

# NA    FA   
# FA 0.043 -    
#   FT 0.043 0.057

qr4D<-qr4%>%filter(phase=="dark")
pairwise.wilcox.test(qr4D$Current_Natural_Log_Ratio, qr4D$condition,
                     p.adjust.method="fdr")

# NA    FA   
# FA 0.086 -    
#   FT 0.086 0.343

qr4FT<-qr4%>%filter(condition=="FT")
pairwise.wilcox.test(qr4FT$Current_Natural_Log_Ratio, qr4FT$phase,
                     p.adjust.method="fdr")

# light
# dark 0.11  

qr4FA<-qr4%>%filter(condition=="FA")
pairwise.wilcox.test(qr4FA$Current_Natural_Log_Ratio, qr4FA$phase,
                     p.adjust.method="fdr")

# light
# dark 0.2

qr4NA<-qr4%>%filter(condition=="NA")
pairwise.wilcox.test(qr4NA$Current_Natural_Log_Ratio, qr4NA$phase,
                     p.adjust.method="fdr")
# light
# dark 0.057

#######################################################################################

#correlational plots of mtg and mtb

mtg_feat<-fread("jrpca_mtg/qurro_mtg/selected_features_mtgtop5perc.tsv")

mtg<-fread("jrpca_mtg/jrpca_mtg_feattab.txt")%>%
  filter(FeatureID %in% mtg_feat$`Feature ID`)
write.table(mtg, "jrpca_mtg/jrpca_mtg5per_feattab.txt",sep = "\t",row.names = FALSE, quote=FALSE)

corr<-fread("jrpca_mtg/correlation_qurro_table/Correlation.tsv")%>%
  dplyr::select(c(1:48))
corr<-corr[c(48:389),]#PF00012.23

write.table(corr, "jrpca_mtg/correlation_qurro_table/Correlation_sub.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

mtg_feat1<-fread("jrpca_mtg/qurro_mtg/selected_features_mtgtop5perc.tsv")
mtg_feat2<-fread("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_m_results/pfam/aldex/SFR23_0620_FANA_450uniquehits_BHless0.1.txt")
mtg<-fread("jrpca_mtg/jrpca_mtg_feattab.txt")%>%
  filter(FeatureID %in% mtg_feat1$`Feature ID`)%>%
  filter(FeatureID %in% mtg_feat2$FeatureID)
write.table(mtg, "jrpca_mtg/jrpca_mtgFANA_450aldex_feattab.txt",sep = "\t",row.names = FALSE, quote=FALSE)

corr<-fread("jrpca_mtg/correlation_qurro2_table/Correlation.tsv")%>%
  dplyr::select(c(1:48))
corr<-corr[c(48:81),]#PF08245.15

write.table(corr, "jrpca_mtg/correlation_qurro2_table/Correlation_sub.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

corr_n<-fread("jrpca_mtg/correlation_qurro2_table/Correlation_sub.tsv",header=TRUE)%>%
  dplyr::rename(FeatureID=featureid)%>%
  left_join(.,key,by="FeatureID")%>%
  dplyr::select(c(50,2:48))%>%
  column_to_rownames("annot")%>%t()%>%as.data.frame()%>%
  rownames_to_column("FeatureID")%>%
  left_join(.,key,by="FeatureID")%>%
  dplyr::select(c(37,2:35))%>%
  column_to_rownames("annot")

selfet_mtg<-fread("jrpca_mtg/jrpca_mtgFANA_450aldex_feattab.txt")
mtg_key<-fread("jrpca_mtg/jrpca_mtg_annotation_key.txt")%>%
  filter(FeatureID %in% selfet_mtg$FeatureID)%>%
  mutate(annot=paste(FeatureID,pfam_name,sep=" "))%>%
  dplyr::rename(GO_process=GO_name)%>%
  dplyr::select(annot,GO_process)%>%
  column_to_rownames("annot")%>%
  mutate(GO_process=ifelse(is.na(GO_process),"unknown",GO_process))

selfet_mtb<-fread("jrpca/qurro_mtb/selected_features_mtbcarb_LCFA.tsv")
mtb_key<-fread("jrpca/jrpca_mtb_annotation_key.txt")%>%
  filter(FeatureID %in% selfet_mtb$`Feature ID`)%>%
  mutate(annot=paste(FeatureID,BIOCHEMICAL,sep=" "))%>%
  dplyr::rename(pathway=SUPER_PATHWAY)%>%
  dplyr::select(annot,pathway)%>%
  column_to_rownames("annot")

# my_colour = list(pathway = c(MTX = "blue", metab = "green4"),
#                  GO_process=c(Lipid=,))annotation_colors = my_colour

draw_colnames_45 <- function (coln, gaps, ...) {
  coord <- pheatmap:::find_coordinates(length(coln), gaps)
  x     <- coord$coord - 0.5 * coord$size
  res   <- grid::textGrob(
    coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"),
    vjust = 0.75, hjust = 1, rot = 90, gp = grid::gpar(...)
  )
  return(res)
}

assignInNamespace(
  x = "draw_colnames",
  value = "draw_colnames_45",
  ns = asNamespace("pheatmap")
)

p<-pheatmap(corr_n, border_color= NA,annotation_row = mtb_key, annotation_col = mtg_key,color= inferno(10))
ggsave("jrpca_mtg/correlation_qurro2_table/SFR23_0728_heatmap_FANAmtg_wbio_annot.pdf", p,width = 16, height = 13)

####################################################################

#run jrpca for mtx vs mtg

mtx<-fread("jrpca_mtx_mtg/jrpca_mtx_feattab.txt")%>%
  mutate(FeatureID=paste(FeatureID,"mtx",sep="-"))

write.table(mtx, "jrpca_mtx_mtg/jrpca_mtx_annot_feattab.txt",sep = "\t",row.names = FALSE, quote=FALSE)

mtg<-fread("jrpca_mtx_mtg/jrpca_mtg_feattab.txt")%>%
  mutate(FeatureID=paste(FeatureID,"mtg",sep="-"))

write.table(mtx, "jrpca_mtx_mtg/jrpca_mtg_annot_feattab.txt",sep = "\t",row.names = FALSE, quote=FALSE)

#run jrpca for mtx vs mtg just 5%

mtx<-fread("jrpca_mtx_mtg/jrpca_mtx5per_feattab.txt")%>%
  mutate(FeatureID=paste(FeatureID,"mtx",sep="-"))

write.table(mtx, "jrpca_mtx_mtg/jrpca_mtx5per_annot_feattab.txt",sep = "\t",row.names = FALSE, quote=FALSE)

mtg<-fread("jrpca_mtx_mtg/jrpca_mtg5per_feattab.txt")%>%
  mutate(FeatureID=paste(FeatureID,"mtg",sep="-"))

write.table(mtg, "jrpca_mtx_mtg/jrpca_mtg5per_annot_feattab.txt",sep = "\t",row.names = FALSE, quote=FALSE)

corr<-fread("jrpca_mtx_mtg//correlation_5per_table/Correlation.tsv")%>%
  dplyr::select(c(1:661))
corr<-corr[c(661:1002),]#PF00012.23-mtg

write.table(corr, "jrpca_mtx_mtg/correlation_5per_table/Correlation_sub.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

list_venn <- list(mtx = mtx$FeatureID,
                  mtg = mtg$FeatureID)

p<-venn.diagram(list_venn, fill = c("white", "black"),height = 10,
                width = 10,alpha = c(0.5, 0.5), lwd =0,filename = NULL)

grid.draw(p)
pdf(file="jrpca_mtx_mtg/correlation_5per_table/venn5per.pdf")
grid.draw(p)
dev.off()
