setwd("~/Notebooks/sfloresr/TRF-metaT/")

library(tidyverse)
library(data.table)
library(ggpubr)
library(RColorBrewer)

##########################################################
#paths
dat_mtx_s<-"data/DE_analysis/aldex_metaT/SFR23_0620_summaldexhits_BHless0.1.txt"
dat_mgx_s<-"data/DE_analysis/aldex_metaG/SFR23_0606_summaldexhits_BHless0.1.txt"
dat_mtx_lds<-"data/DE_analysis/aldex_metaT/SFR23_0620_summaldexhitsLD_BHless0.1.txt"
dat_mgx_lds<-"data/DE_analysis/aldex_metaG/SFR23_0606_summaldexhitsLD_BHless0.1.txt"
res_path_mtx<-"data/DE_analysis/aldex_metaT/SFR23_0620_"
res_path_mgx<-"data/DE_analysis/aldex_metaG/SFR23_0606_"
dat_annot<-"data/pfam_metaT/pfam_annotationkey.csv"
fig_path<-"figures/DE_analysis/"
##########################################################
#functions

get_chi<-function(dat,compar){
  chi<-dat%>%
    filter(comparison==compar)%>%
    column_to_rownames("method")%>%
    dplyr::select(3,4)%>%
    as.matrix()
  res<-chisq.test(chi,simulate.p.value=TRUE, B=2000) 
  return(res)
}

get_chi_LD<-function(dat,compar,ld){
  chi<-dat%>%
    filter(comparison==compar & phase==ld)%>%
    column_to_rownames("method")%>%
    dplyr::select(4:5)%>%
    as.matrix()
  res<-chisq.test(chi,simulate.p.value=TRUE, B=2000) 
  return(res)
}

clean_aldex_res<-function(dat){
  sub_dat<-fread(dat)%>%
    filter(!grepl("DUF",Name))%>%
    filter(wi.eBH<0.1)
  return(sub_dat)
}

jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

##########################################################
#load aldex summary files

mtx_sum<-fread(dat_mtx_s)%>%mutate(method="MTX")
mgx_sum<-fread(dat_mgx_s)%>%mutate(method="MGX")
mtx_sum_ld<-fread(dat_mtx_lds)%>%mutate(method="MTX")
mgx_sum_ld<-fread(dat_mgx_lds)%>%mutate(method="MGX")

#########################################################
#plot mtx and mgx aldex summary files together --Figure 2A
sum_file<-rbind(mtx_sum,mgx_sum)%>%
  mutate(comparison=factor(comparison, levels = c("NAFA","NAFT","FAFT")))

ggplot(sum_file, aes(x=comparison, y=perc,fill=method)) + 
  geom_bar(stat="identity", position="dodge") +
  scale_fill_manual(values=c("red","blue"))+
  theme_minimal()+scale_y_continuous(limits=c(0,30))+
  theme(legend.position = "top")

ggsave(paste0(fig_path,"SFR23_0620_summaldexhits_MTMG_BHless0.1.pdf"), width = 2.2, height = 3)

set.seed(1234)
get_chi(sum_file,"NAFA") #X-squared = 4.6056, df = NA, p-value = 0.03148
get_chi(sum_file,"NAFT") #X-squared = 197.81, df = NA, p-value = 0.0004998
get_chi(sum_file,"FAFT") #X-squared = 8.6248, df = NA, p-value = 0.003998

##########################################################
#plot mtx and mgx aldex summary LD files together --Figure 2B

sum_file<-rbind(mtx_sum_ld,mgx_sum_ld)%>%
  mutate(comparison=factor(comparison, levels = c("NAFA","NAFT","FAFT")),
         phase=factor(phase, levels = c("light","dark")),
         perc=abs(perc))

ggplot(sum_file, aes(x=comparison, y=perc,fill=method)) + 
  geom_bar(stat="identity", position="dodge") +
  facet_wrap(vars(phase), nrow = 1)+
  scale_fill_manual(values=c("red","blue"))+
  theme_bw()+scale_y_continuous(expand=c(0,0),limits=c(0,25))+
  theme(legend.position = "top")

ggsave(paste0(fig_path,"SFR24_0528_summaldexhitsLD_MTMG_BHless0.1_facet.pdf"), width = 4, height = 3)

get_chi_LD(sum_file,"NAFA","light") #X-squared = 17.052, df = NA, p-value = 0.0004998
get_chi_LD(sum_file,"NAFT","light") #X-squared = 440.09, df = NA, p-value = 0.0004998
get_chi_LD(sum_file,"FAFT","light") #X-squared = NaN, df = 1, p-value = NA
get_chi_LD(sum_file,"NAFA","dark") #X-squared = 70.183, df = NA, p-value = 0.0004998
get_chi_LD(sum_file,"NAFT","dark") #X-squared = 564.07, df = NA, p-value = 0.0004998
get_chi_LD(sum_file,"FAFT","dark") #X-squared = 199.89, df = NA, p-value = 0.0004998

##########################################################
#jaccard distance between mgx and mtx for the light phase--Figure 2C (top)

#read in light aldex hits-MTX
mt_faftl<-clean_aldex_res(paste0(res_path_mtx,"FAFTL_ald_effectwpval_wannot.txt")) #0
mt_ftnal<-clean_aldex_res(paste0(res_path_mtx,"FTNAL_ald_effectwpval_wannot.txt")) #1342
mt_fanal<-clean_aldex_res(paste0(res_path_mtx,"FANAL_ald_effectwpval_wannot.txt")) #1119

#read in light aldex hits-MGX
mg_faftl<-clean_aldex_res(paste0(res_path_mgx,"FAFTL_ald_effectwpval_wannot.txt"))#0
mg_ftnal<-clean_aldex_res(paste0(res_path_mgx,"FTNAL_ald_effectwpval_wannot.txt"))#162
mg_fanal<-clean_aldex_res(paste0(res_path_mgx,"FANAL_ald_effectwpval_wannot.txt"))#472

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
ggsave(paste0(fig_path,"SFR23_0620_jaccard_aldex_MTLrMGL_summ_small.pdf"),plot=p, width = 4, height = 3)

##########################################################
#jaccard distance between mgx and mtx for the light phase--Figure 2C (bottom)

#read in dark aldex hits-MTX
mt_faftd<-clean_aldex_res(paste0(res_path_mtx,"FAFTD_ald_effectwpval_wannot.txt")) #311
mt_ftnad<-clean_aldex_res(paste0(res_path_mtx,"FTNAD_ald_effectwpval_wannot.txt")) #1085
mt_fanad<-clean_aldex_res(paste0(res_path_mtx,"FANAD_ald_effectwpval_wannot.txt")) #1133

#read in dark aldex hits-MGX
mg_faftd<-clean_aldex_res(paste0(res_path_mgx,"FAFTD_ald_effectwpval_wannot.txt"))#0
mg_ftnad<-clean_aldex_res(paste0(res_path_mgx,"FTNAD_ald_effectwpval_wannot.txt"))#53
mg_fanad<-clean_aldex_res(paste0(res_path_mgx,"FANAD_ald_effectwpval_wannot.txt"))#419

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
ggsave(paste0(fig_path,"SFR23_0620_jaccard_aldex_MTDrMGD_summ_small.pdf"),plot=p, width = 4, height = 3)

##########################################################
#jaccard distance or just mtx light vs. dark--Figure S2B

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
ggsave(paste0(fig_path,"SFR23_0607_jaccard_aldex_MTLD_summ_small.pdf"),plot=p, width = 3.5, height = 4)
