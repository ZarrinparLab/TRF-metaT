setwd("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/multiomics")

library(tidyverse)
library(data.table)
library(ggpubr)
library(ggvenn)
library(gplots)
library(viridis)
library(RColorBrewer)
library("qiime2R")
library(ggpattern)
##########################################################

#load MTX aldex summary file
mtx_sum<-fread("~/scratch/TRF_multiomics/metatranscript/woltka2_m_results/pfam/aldex/SFR23_0620_summaldexhits_BHless0.1.txt")%>%
  mutate(method="MTX")
#load MGX aldex summary file
mgx_sum<-fread("~/scratch/TRF_multiomics/metagenomic/woltka2_results/filtered_metaG/pfam_notnorm/aldex/SFR23_0606_summaldexhits_BHless0.1.txt")%>%
  mutate(method="MGX")

sum_file<-rbind(mtx_sum,mgx_sum)%>%
  mutate(comparison=factor(comparison, levels = c("NAFA","NAFT","FAFT")))

ggplot(sum_file, aes(x=comparison, y=perc,fill=method)) + 
  geom_bar(stat="identity", position="dodge") +
  scale_fill_manual(values=c("red","blue"))+
  theme_minimal()+scale_y_continuous(limits=c(0,30))+
  theme(legend.position = "top")

ggsave("SFR23_0620_summaldexhits_MTMG_BHless0.1.pdf", width = 2.2, height = 3)

set.seed(1234)
#FANA
chi<-sum_file%>%
  filter(comparison=="NAFA")%>%
  column_to_rownames("method")%>%
  dplyr::select(3,4)%>%
  as.matrix()

chisq.test(chi,simulate.p.value=TRUE, B=2000) #X-squared = 4.6056, df = NA, p-value = 0.03148

#FTNA
chi<-sum_file%>%
  filter(comparison=="NAFT")%>%
  column_to_rownames("method")%>%
  dplyr::select(3,4)%>%
  as.matrix()

chisq.test(chi,simulate.p.value=TRUE, B=2000) #X-squared = 197.81, df = NA, p-value = 0.0004998

#FAFT
chi<-sum_file%>%
  filter(comparison=="FAFT")%>%
  column_to_rownames("method")%>%
  dplyr::select(3,4)%>%
  as.matrix()

chisq.test(chi, simulate.p.value=TRUE, B=2000) #X-squared = 8.6248, df = NA, p-value = 0.003998

##########################################################

#load MTX aldex summary file
mtx_sum<-fread("~/scratch/TRF_multiomics/metatranscript/woltka2_m_results/pfam/aldex/SFR23_0620_summaldexhitsLD_BHless0.1.txt")%>%
  mutate(method="MTX")
#load MGX aldex summary file
mgx_sum<-fread("~/scratch/TRF_multiomics/metagenomic/woltka2_results/filtered_metaG/pfam_notnorm/aldex/SFR23_0606_summaldexhitsLD_BHless0.1.txt")%>%
  mutate(method="MGX")

sum_file<-rbind(mtx_sum,mgx_sum)%>%
  mutate(comparison=factor(comparison, levels = c("NAFA","NAFT","FAFT")),
         phase=factor(phase, levels = c("light","dark")),
         perc=abs(perc))

ggplot(sum_file, aes(x=comparison, y=perc,fill=method)) + 
  geom_bar(stat="identity", position="dodge") +
  facet_wrap(vars(phase), nrow = 1)+
  scale_fill_manual(values=c("red","blue"))+
  theme_bw()+scale_y_continuous(expand=c(0,0),limits=c(0,25))+
  theme(legend.position = "top")

ggsave("SFR24_0528_summaldexhitsLD_MTMG_BHless0.1_facet.pdf", width = 4, height = 3)

# ggplot(sum_file, aes(comparison, perc)) +
#   geom_col_pattern(
#     aes(pattern = phase, fill = method, group=method), 
#     colour                   = 'black',
#     pattern_density          = 0.1, 
#     pattern_key_scale_factor = 1.3,
#     pattern_fill = 'black',
#     position = 'dodge') +
#   theme_minimal() +
#   labs(y="percent (%)") +
#   scale_fill_manual(values = c(MTX='blue', MGX='red')) + 
#   scale_pattern_manual(values = c("none", "stripe")) +  
#   theme(legend.position = 'right')
# 
# ggsave("SFR24_0307_summaldexhitsLD_MTMG_BHless0.1_hashed.pdf", width = 4, height = 3)

# sum_file %>% 
#   group_by(method, comparison) %>% 
#   mutate(comparison=factor(comparison, levels = c("NAFA","NAFT","FAFT")),
#          phase=factor(phase, levels = c("light","dark")))%>%
#   ggplot(aes(comparison, perc, fill =method)) + 
#   geom_col(data = . %>% filter( phase=="light"), position = position_dodge(width = 0.9), alpha = 0.4) +
#   geom_col(data = . %>% filter( phase=="dark"), position = position_dodge(width = 0.9), alpha = 1) +
#   geom_tile(aes(y=NA_integer_, alpha = factor(phase))) + 
#   scale_alpha_manual(values = c(0.4,1)) +
#   scale_fill_manual(values=c("red","blue"))+theme_minimal()

#ggsave("SFR23_0620_summaldexhitsLD_MTMG_BHless0.1.pdf", width = 4, height = 3)

#FANAL
chi<-sum_file%>%
  filter(comparison=="NAFA" & phase=="light")%>%
  column_to_rownames("method")%>%
  dplyr::select(4:5)%>%
  as.matrix()

chisq.test(chi, simulate.p.value=TRUE, B=2000) #X-squared = 17.052, df = NA, p-value = 0.0004998

#FTNAL
chi<-sum_file%>%
  filter(comparison=="NAFT" & phase=="light")%>%
  column_to_rownames("method")%>%
  dplyr::select(4:5)%>%
  as.matrix()

chisq.test(chi, simulate.p.value=TRUE, B=2000) #X-squared = 440.09, df = NA, p-value = 0.0004998

#FAFTL
chi<-sum_file%>%
  filter(comparison=="FAFT" & phase=="light")%>%
  column_to_rownames("method")%>%
  dplyr::select(4:5)%>%
  as.matrix()

chisq.test(chi,simulate.p.value=TRUE, B=2000) #X-squared = NaN, df = 1, p-value = NA

#FANAD
chi<-sum_file%>%
  filter(comparison=="NAFA" & phase=="dark")%>%
  column_to_rownames("method")%>%
  dplyr::select(4:5)%>%
  as.matrix()

chisq.test(chi,simulate.p.value=TRUE, B=2000) #X-squared = 70.183, df = NA, p-value = 0.0004998

#FTNAD
chi<-sum_file%>%
  filter(comparison=="NAFT" & phase=="dark")%>%
  column_to_rownames("method")%>%
  dplyr::select(4:5)%>%
  as.matrix()

chisq.test(chi,simulate.p.value=TRUE, B=2000) #X-squared = 564.07, df = NA, p-value = 0.0004998

#FAFTD
chi<-sum_file%>%
  filter(comparison=="FAFT" & phase=="dark")%>%
  column_to_rownames("method")%>%
  dplyr::select(4:5)%>%
  as.matrix()

chisq.test(chi,simulate.p.value=TRUE, B=2000) #X-squared = 199.89, df = NA, p-value = 0.0004998
##########################################################
annot<-fread("~/scratch/TRF_multiomics/metatranscript/woltka2_m_results/pfam/pfam_annotationkey.csv")

#read in light aldex hits-MTX

#FAFTL
mt_faftl<-fread("~/scratch/TRF_multiomics/metatranscript/woltka2_m_results/pfam/aldex/SFR23_0620_FAFTL_ald_effectwpval.txt")%>%
  left_join(.,annot, by ="FeatureID") %>%
  filter(!grepl("DUF",Name))%>%
  filter(wi.eBH<0.1) #0

#FTNAL
mt_ftnal<-fread("~/scratch/TRF_multiomics/metatranscript/woltka2_m_results/pfam/aldex/SFR23_0620_FTNAL_ald_effectwpval.txt")%>%
  left_join(.,annot, by ="FeatureID") %>%
  filter(!grepl("DUF",Name))%>%
  filter(wi.eBH<0.1) #1342

#FANAL
mt_fanal<-fread("~/scratch/TRF_multiomics/metatranscript/woltka2_m_results/pfam/aldex/SFR23_0620_FANAL_ald_effectwpval.txt")%>%
  left_join(.,annot, by ="FeatureID") %>%
  filter(!grepl("DUF",Name))%>%
  filter(wi.eBH<0.1) #1119

#read in light aldex hits-MGX

#FAFTL
mg_faftl<-fread("~/scratch/TRF_multiomics/metagenomic/woltka2_results/filtered_metaG/pfam_notnorm/aldex/SFR23_0606_FAFTL_ald_effectwpval.txt")%>%
  left_join(.,annot, by ="FeatureID") %>%
  filter(!grepl("DUF",Name))%>%
  filter(wi.eBH<0.1) #0

#FTNAL
mg_ftnal<-fread("~/scratch/TRF_multiomics/metagenomic/woltka2_results/filtered_metaG/pfam_notnorm/aldex/SFR23_0606_FTNAL_ald_effectwpval.txt")%>%
  left_join(.,annot, by ="FeatureID") %>%
  filter(!grepl("DUF",Name))%>%
  filter(wi.eBH<0.1) #162

#FANAL
mg_fanal<-fread("~/scratch/TRF_multiomics/metagenomic/woltka2_results/filtered_metaG/pfam_notnorm/aldex/SFR23_0606_FANAL_ald_effectwpval.txt")%>%
  left_join(.,annot, by ="FeatureID") %>%
  filter(!grepl("DUF",Name))%>%
  filter(wi.eBH<0.1) #472


jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

j_val<-c(jaccard(mg_fanal$FeatureID,mt_fanal$FeatureID),jaccard(mg_fanal$FeatureID,mt_ftnal$FeatureID),jaccard(mg_fanal$FeatureID,mt_faftl$FeatureID),
         jaccard(mg_ftnal$FeatureID,mt_ftnal$FeatureID),jaccard(mg_ftnal$FeatureID,mt_faftl$FeatureID),
         jaccard(mg_faftl$FeatureID,mg_faftl$FeatureID))


var1<-c(rep("NAFA-MGX",3),rep("NAFT-MGX",2),"FAFT-MGX")
var2<-c("NAFA-MTX","NAFT-MTX","FAFT-MTX",
        "NAFT-MTX","FAFT-MTX",
        "FAFT-MTX")
jdf<-data.frame(var1,var2,j_val)%>%mutate(var1=factor(var1,levels=c("NAFA-MGX","NAFT-MGX","FAFT-MGX")),
                                          var2=factor(var2,levels=c("FAFT-MTX","NAFT-MTX","NAFA-MTX")))

p<-ggplot(jdf, aes(x = var1, y = var2, fill = j_val)) +
  geom_tile() +
  geom_text(aes(label = signif(j_val,digits=3))) +
  scale_fill_distiller(palette = "Oranges",direction=1) +
  theme_minimal() +
  theme(panel.grid = element_blank(),legend.position = "right")+labs(title="Jaccard Index",x="",y="")
ggsave("SFR23_0620_jaccard_aldex_MTLrMGL_summ_small.pdf",plot=p, width = 4, height = 3)


#read in dark aldex hits-MTX

#FAFTL
mt_faftd<-fread("~/scratch/TRF_multiomics/metatranscript/woltka2_m_results/pfam/aldex/SFR23_0620_FAFTD_ald_effectwpval.txt")%>%
  left_join(.,annot, by ="FeatureID") %>%
  filter(!grepl("DUF",Name))%>%
  filter(wi.eBH<0.1) #311

#FTNAL
mt_ftnad<-fread("~/scratch/TRF_multiomics/metatranscript/woltka2_m_results/pfam/aldex/SFR23_0620_FTNAD_ald_effectwpval.txt")%>%
  left_join(.,annot, by ="FeatureID") %>%
  filter(!grepl("DUF",Name))%>%
  filter(wi.eBH<0.1) #1085

#FANAL
mt_fanad<-fread("~/scratch/TRF_multiomics/metatranscript/woltka2_m_results/pfam/aldex/SFR23_0620_FANAD_ald_effectwpval.txt")%>%
  left_join(.,annot, by ="FeatureID") %>%
  filter(!grepl("DUF",Name))%>%
  filter(wi.eBH<0.1) #1133

#read in aldex hits-MGX

#FAFTL
mg_faftd<-fread("~/scratch/TRF_multiomics/metagenomic/woltka2_results/filtered_metaG/pfam_notnorm/aldex/SFR23_0606_FAFTD_ald_effectwpval.txt")%>%
  left_join(.,annot, by ="FeatureID") %>%
  filter(!grepl("DUF",Name))%>%
  filter(wi.eBH<0.1) #0

#FTNAL
mg_ftnad<-fread("~/scratch/TRF_multiomics/metagenomic/woltka2_results/filtered_metaG/pfam_notnorm/aldex/SFR23_0606_FTNAD_ald_effectwpval.txt")%>%
  left_join(.,annot, by ="FeatureID") %>%
  filter(!grepl("DUF",Name))%>%
  filter(wi.eBH<0.1) #53

#FANAL
mg_fanad<-fread("~/scratch/TRF_multiomics/metagenomic/woltka2_results/filtered_metaG/pfam_notnorm/aldex/SFR23_0606_FANAD_ald_effectwpval.txt")%>%
  left_join(.,annot, by ="FeatureID") %>%
  filter(!grepl("DUF",Name))%>%
  filter(wi.eBH<0.1) #419


jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

j_val<-c(jaccard(mg_fanad$FeatureID,mt_fanad$FeatureID),jaccard(mg_fanad$FeatureID,mt_ftnad$FeatureID),jaccard(mg_fanad$FeatureID,mt_faftd$FeatureID),
         jaccard(mg_ftnad$FeatureID,mt_ftnad$FeatureID),jaccard(mg_ftnad$FeatureID,mt_faftd$FeatureID),
         jaccard(mg_faftd$FeatureID,0))


var1<-c(rep("NAFA-MGX",3),rep("NAFT-MGX",2),"FAFT-MGX")
var2<-c("NAFA-MTX","NAFT-MTX","FAFT-MTX",
        "NAFT-MTX","FAFT-MTX",
        "FAFT-MTX")
jdf<-data.frame(var1,var2,j_val)%>%mutate(var1=factor(var1,levels=c("NAFA-MGX","NAFT-MGX","FAFT-MGX")),
                                          var2=factor(var2,levels=c("FAFT-MTX","NAFT-MTX","NAFA-MTX")))

p<-ggplot(jdf, aes(x = var1, y = var2, fill = j_val)) +
  geom_tile() +
  geom_text(aes(label = signif(j_val,digits=3))) +
  scale_fill_distiller(palette = "Purples",direction=1) +
  theme_minimal() +
  theme(panel.grid = element_blank(),legend.position = "right")+labs(title="Jaccard Index",x="",y="")
ggsave("SFR23_0620_jaccard_aldex_MTDrMGD_summ_small.pdf",plot=p, width = 4, height = 3)

##########################################################

#get the jaccard distance for just mtx light vs. dark
j_val<-c(jaccard(mt_fanal$FeatureID,mt_fanad$FeatureID),jaccard(mt_fanal$FeatureID,mt_ftnad$FeatureID),jaccard(mt_fanal$FeatureID,mt_faftd$FeatureID),
         jaccard(mt_ftnal$FeatureID,mt_ftnad$FeatureID),jaccard(mt_ftnal$FeatureID,mt_faftd$FeatureID),
         jaccard(mt_faftl$FeatureID,mt_faftd$FeatureID))


var1<-c(rep("NAFA-L",3),rep("NAFT-L",2),"FAFT-L")
var2<-c("NAFA-D","NAFT-D","FAFT-D",
        "NAFT-D","FAFT-D",
        "FAFT-D")
jdf<-data.frame(var1,var2,j_val)%>%mutate(var1=factor(var1,levels=c("NAFA-L","NAFT-L","FAFT-L")),
                                          var2=factor(var2,levels=c("FAFT-D","NAFT-D","NAFA-D")))

p<-ggplot(jdf, aes(x = var1, y = var2, fill = j_val)) +
  geom_tile() +
  geom_text(aes(label = signif(j_val,digits=3))) +
  scale_fill_distiller(palette = "Blues",direction=1) +
  theme_minimal() +
  theme(panel.grid = element_blank(),legend.position = "top")+labs(title="Jaccard Index",x="",y="")
ggsave("SFR23_0607_jaccard_aldex_MTLD_summ_small.pdf",plot=p, width = 3.5, height = 4)

##########################################################

#load MTX aldex summary file for LD comparison
mtx_sum<-fread("~/scratch/TRF_multiomics/metatranscript/woltka2_m_results/pfam/aldex/SFR23_0623_summaldexhits_justLD_BHless0.1.txt")%>%
  mutate(method="MTX",
         comparison=ifelse(is.na(comparison),"NA",comparison))
#load MGX aldex summary file for LD comparison
mgx_sum<-fread("~/scratch/TRF_multiomics/metagenomic/woltka2_results/filtered_metaG/pfam_notnorm/aldex/SFR23_0623_summaldexhits_justLD_BHless0.1.txt")%>%
  mutate(method="MGX",
         comparison=ifelse(is.na(comparison),"NA",comparison))

sum_file<-rbind(mtx_sum,mgx_sum)%>%
  mutate(comparison=factor(comparison, levels = c("NA","FA","FT")))

ggplot(sum_file, aes(x=comparison, y=perc,fill=method)) + 
  geom_bar(stat="identity", position="dodge") +
  scale_fill_manual(values=c("red","blue"))+
  theme_minimal()+scale_y_continuous(limits=c(0,10))+
  theme(legend.position = "top")

ggsave("SFR23_0623_summaldexhits_MTMG_justLD_BHless0.1.pdf", width = 2.2, height = 2.5)

set.seed(1234)

#NA (light vs. dark)
chi<-sum_file%>%
  filter(comparison=="NA")%>%
  column_to_rownames("method")%>%
  dplyr::select(3,4)%>%
  as.matrix()

chisq.test(chi,simulate.p.value=TRUE, B=2000) #X-squared = 134.71, df = NA, p-value = 0.0004998
