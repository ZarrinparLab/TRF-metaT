setwd("~/Notebooks/sfloresr/TRF-metaT/")

library(tidyverse)
library(data.table)
library(ALDEx2)
library(ggpubr)

###########################################################
#inputs
dat_nn<-"g-diversity-core-metrics114K/rarefied_table/pfam-144k.tsv"
dat_metadata<-"data/metaT_metadata_ztcat_noNT.txt"
dat_annot<-"data/pfam_metaT/pfam_annotationkey.csv"
dat_path<-"data/pfam_metaT/"
res_path<-"data/DE_analysis/aldex_metaT/SFR24_0529_"
fig_path<-"figures/DE_analysis/"
###########################################################
#functions
run_aldex<-function(dat,condsAB){
  AB.ald<-aldex.clr(round(dat),condsAB, mc.samples=500, denom="all", verbose=F)
  AB.ttest<-aldex.ttest(AB.ald,paired.test = FALSE, hist.plot=FALSE)%>%rownames_to_column("FeatureID")
  AB.effect<-aldex.effect(AB.ald)%>%rownames_to_column("FeatureID")%>%left_join(.,AB.ttest,by="FeatureID")
  return(AB.effect)
}

get_chi_LD<-function(dat,compar){
  chi<-dat%>%
    filter(comparison==compar)%>%
    column_to_rownames("phase")%>%
    dplyr::select(2:3)%>%
    as.matrix()
  res<-chisq.test(chi,simulate.p.value=TRUE, B=2000) 
  return(res)
}
###########################################################
#create pairwise files

#FA vs FT
pfamFAFT<-fread(dat_nn)%>%dplyr::select(c(1:25))
write.table(pfamFAFT,"g-diversity-core-metrics114K/rarefied_table/feature-table-FAFT.txt",sep = "\t",row.names = FALSE, quote=FALSE)

#just light
pfamFAFTL<-fread(dat_nn)%>%dplyr::select(c(1:7,14:19))
write.table(pfamFAFTL,"g-diversity-core-metrics114K/rarefied_table/feature-table-FAFTL.txt",sep = "\t",row.names = FALSE, quote=FALSE)

#just dark
pfamFAFTD<-fread(dat_nn)%>%dplyr::select(c(1,8:13,20:25))
write.table(pfamFAFTD,"g-diversity-core-metrics114K/rarefied_table/feature-table-FAFTD.txt",sep = "\t",row.names = FALSE, quote=FALSE)

#FT vs NA
pfamFTNA<-fread(dat_nn)%>%dplyr::select(c(1,14:43))
write.table(pfamFTNA,"g-diversity-core-metrics114K/rarefied_table/feature-table-FTNA.txt",sep = "\t",row.names = FALSE, quote=FALSE)

#light
pfamFTNAL<-fread(dat_nn)%>%dplyr::select(c(1,14:19,26:34))
write.table(pfamFTNAL,"g-diversity-core-metrics114K/rarefied_table/feature-table-FTNAL.txt",sep = "\t",row.names = FALSE, quote=FALSE)

#dark
pfamFTNAD<-fread(dat_nn)%>%dplyr::select(c(1,20:25,35:43))
write.table(pfamFTNAD,"g-diversity-core-metrics114K/rarefied_table/feature-table-FTNAD.txt",sep = "\t",row.names = FALSE, quote=FALSE)

#FA vs NA
pfamFANA<-fread(dat_nn)%>%dplyr::select(c(1:13,26:43))
write.table(pfamFANA,"g-diversity-core-metrics114K/rarefied_table/feature-table-FANA.txt",sep = "\t",row.names = FALSE, quote=FALSE)

#light
pfamFANAL<-fread(dat_nn)%>%dplyr::select(c(1:7,26:34))
write.table(pfamFANAL,"g-diversity-core-metrics114K/rarefied_table/feature-table-FANAL.txt",sep = "\t",row.names = FALSE, quote=FALSE)

#dark
pfamFANAD<-fread(dat_nn)%>%dplyr::select(c(1,8:13,35:43))
write.table(pfamFANAD,"g-diversity-core-metrics114K/rarefied_table/feature-table-FANAD.txt",sep = "\t",row.names = FALSE, quote=FALSE)

###########################################################
#run aldex
annot<-fread(dat_annot)
md<-fread(dat_metadata)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))
conds <- md$condition

##FA vs. FT
pfamFAFT<-fread("g-diversity-core-metrics114K/rarefied_table/feature-table-FAFT.txt")%>%column_to_rownames("FeatureID")
condsFAvsFT<-conds[1:24]
FAFT.effect<-run_aldex(pfamFAFT,condsFAvsFT)
write.table(FAFT.effect.annot,paste0(res_path,"FAFT_ald_effectwpval_144k.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

FAFT.effect<-fread("pfam/aldex/SFR24_0529_FAFT_ald_effectwpval_144K.txt")

#FA vs FT (light)
md<-fread(dat_metadata)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition)) %>%
  filter(zt_time<13)
conds <- md$condition

condsFAvsFTL<-conds[1:12]
pfamFAFTL<-fread("g-diversity-core-metrics114K/rarefied_table/feature-table-FAFTL.txt")%>%column_to_rownames("FeatureID")
FAFTL.effect<-run_aldex(pfamFAFTL,condsFAvsFTL)
write.table(FAFTL.effect,paste0(res_path,"FAFTL_ald_effectwpval_144K.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

#FA vs FT (dark)
md<-fread(dat_metadata)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition)) %>%
  filter(zt_time>9)
conds <- md$condition

condsFAvsFTD<-conds[1:12]
pfamFAFTD<-fread("g-diversity-core-metrics114K/rarefied_table/feature-table-FAFTD.txt")%>%column_to_rownames("FeatureID")
FAFTD.effect<-run_aldex(pfamFAFTD,condsFAvsFTD)
write.table(FAFTD.effect,paste0(res_path,"FAFTD_ald_effectwpval_144K.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

##NA vs FT
md<-fread(dat_metadata)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))
conds <- md$condition

pfamFTNA<-fread("g-diversity-core-metrics114K/rarefied_table/feature-table-FTNA.txt")%>%column_to_rownames("FeatureID")
condsFTvsNA<-conds[c(13:42)]
FTNA.effect<-run_aldex(pfamFTNA,condsFTvNA)
write.table(FTNA.effect,paste0(res_path,"FTNA_ald_effectwpval_144K.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

##NA vs FT (light)
md<-fread(dat_metadata)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition)) %>%
  filter(zt_time<13)
conds <- md$condition

condsFTvsNAL<-conds[c(7:21)]
pfamFTNAL<-fread("g-diversity-core-metrics114K/rarefied_table/feature-table-FTNAL.txt")%>%column_to_rownames("FeatureID")
FTNAL.effect<-run_aldex(pfamFTNAL,condsFTvNAL)
write.table(FTNAL.effect,paste0(res_path,"FTNAL_ald_effectwpval_144K.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

#NA vs FT (dark)
md<-fread(dat_metadata)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition)) %>%
  filter(zt_time>9)
conds <- md$condition

condsFTvsNAD<-conds[c(7:21)]
pfamFTNAD<-fread("g-diversity-core-metrics114K/rarefied_table/feature-table-FTNAD.txt")%>%column_to_rownames("FeatureID")
FTNALD.effect<-run_aldex(pfamFTNAD,condsFTvNAD)
write.table(FTNAD.effect,paste0(res_path,"FTNAD_ald_effectwpval_144K.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

##NA vs FA
md<-fread(dat_metadata)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))
conds <- md$condition

pfamFANA<-fread("g-diversity-core-metrics114K/rarefied_table/feature-table-FANA.txt")%>%column_to_rownames("FeatureID")
condsFAvsNA<-conds[c(1:12,25:42)]
FANA.effect<-run_aldex(pfamFANA,condsFAvsNA)
write.table(FANA.effect,paste0(res_path,"FANA_ald_effectwpval_144K.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

##NA vs FA (light)
md<-fread(dat_metadata)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition)) %>%
  filter(zt_time<13)
conds <- md$condition

condsFAvsNAL<-conds[c(1:6,13:21)]
pfamFANAL<-fread("g-diversity-core-metrics114K/rarefied_table/feature-table-FANAL.txt")%>%column_to_rownames("FeatureID")
FANAL.effect<-run_aldex(pfamFANAL,condsFAvsNAL)
write.table(FANAL.effect,paste0(res_path,"FANAL_ald_effectwpval_144K.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

#NA vs FA (dark)
md<-fread(dat_metadata)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition)) %>%
  filter(zt_time>9)
conds <- md$condition

condsFAvsNAD<-conds[c(1:6,13:21)]
pfamFANAD<-fread("g-diversity-core-metrics114K/rarefied_table/feature-table-FANAD.txt")%>%column_to_rownames("FeatureID")
FANAD.effect<-run_aldex(pfamFANAD,condsFAvsNAD)
write.table(FANAD.effect,paste0(res_path,"FANAD_ald_effectwpval_144K.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

###########################################################
#summary of aldex hits for rarefied mtx--Figure S2A(left)
a<-(nrow(FANA.effect.annot%>%filter(wi.eBH<0.1))/nrow(FANA.effect.annot))*100
b<-(nrow(FTNA.effect.annot%>%filter(wi.eBH<0.1))/nrow(FTNA.effect.annot))*100
c<-(nrow(FAFT.effect.annot%>%filter(wi.eBH<0.1))/nrow(FAFT.effect.annot))*100

aldex_summhits <- data.frame (comparison  = c("NAFA","NAFT","FAFT"),
                              perc= c(a,b,c)) %>%
  mutate(comparison=factor(comparison, levels = c("NAFA","NAFT","FAFT")))

ggplot(aldex_summhits, aes(x=comparison, y=perc)) + 
  geom_bar(stat="identity", position="identity",fill="gray30") +
  theme_minimal() 

ggsave(paste0(fig_path,"SFR24_0529_summaldexhits_pless0.05_144K.pdf"), width = 2.2, height = 3)

###########################################################
#summary of aldex hits for rarefied mtx by light and dark--Figure S2A(right)

al<-(nrow(FANAL.effect.annot%>%filter(wi.eBH<0.1))/nrow(FANAL.effect.annot))*100
bl<-(nrow(FTNAL.effect.annot%>%filter(wi.eBH<0.1))/nrow(FTNAL.effect.annot))*100
cl<-(nrow(FAFTL.effect.annot%>%filter(wi.eBH<0.1))/nrow(FAFTL.effect.annot))*100

ad<-(nrow(FANAD.effect.annot%>%filter(wi.eBH<0.1))/nrow(FANAD.effect.annot))*100
bd<-(nrow(FTNAD.effect.annot%>%filter(wi.eBH<0.1))/nrow(FTNAD.effect.annot))*100
cd<-(nrow(FAFTD.effect.annot%>%filter(wi.eBH<0.1))/nrow(FAFTD.effect.annot))*100

aldex_summhitsLD <- data.frame (comparison  = c("NAFA","NAFT",,"FAFT","NAFA","NAFT","FAFT"),
                              phase= c("light","light","light","dark","dark","dark"),
                              perc= c(al,bl,cl,ad,bd,cd)) %>%
  mutate(comparison=factor(comparison, levels = c("NAFA","NAFT","FAFT")),
         phase=factor(phase, levels = c("light","dark")))

ggplot(aldex_summhitsLD, aes(x=comparison, y=perc,fill=phase)) + 
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_manual(values=c("gray70","gray10")) +
  theme_minimal() +
  theme(legend.position = "top")

ggsave(paste0(fig_path,"SFR24_0529_summaldexhitsLD_pless0.05_144K.pdf"), width = 2.2, height = 3)

#to calc chi-squared

acl<-nrow(FANAL.effect.annot%>%filter(wi.eBH<0.1))
bcl<-nrow(FTNAL.effect.annot%>%filter(wi.eBH<0.1))
ccl<-nrow(FAFTL.effect.annot%>%filter(wi.eBH<0.1))

acd<-nrow(FANAD.effect.annot%>%filter(wi.eBH<0.1))
bcd<-nrow(FTNAD.effect.annot%>%filter(wi.eBH<0.1))
ccd<-nrow(FAFTD.effect.annot%>%filter(wi.eBH<0.1))

aldex_summhitsLD_stat <- data.frame (comparison  = c("NAFA","NAFT","FAFT","NAFA","NAFT","FAFT"),
                                phase= c("light","light","light","dark","dark","dark"),
                                sig= c(nrow(FANAL.effect.annot%>%filter(wi.eBH<0.1)),
                                       nrow(FTNAL.effect.annot%>%filter(wi.eBH<0.1)),
                                       nrow(FAFTL.effect.annot%>%filter(wi.eBH<0.1)),
                                       nrow(FANAD.effect.annot%>%filter(wi.eBH<0.1)),
                                       nrow(FTNAD.effect.annot%>%filter(wi.eBH<0.1)),
                                       nrow(FAFTD.effect.annot%>%filter(wi.eBH<0.1))),
                                not_sig=c(nrow(FANAL.effect.annot)-nrow(FANAL.effect.annot%>%filter(wi.eBH<0.1)),
                                          nrow(FTNAL.effect.annot)-nrow(FTNAL.effect.annot%>%filter(wi.eBH<0.1)),
                                          nrow(FAFTL.effect.annot)-nrow(FAFTL.effect.annot%>%filter(wi.eBH<0.1)),
                                          nrow(FANAD.effect.annot)-nrow(FANAD.effect.annot%>%filter(wi.eBH<0.1)),
                                          nrow(FTNAD.effect.annot)-nrow(FTNAD.effect.annot%>%filter(wi.eBH<0.1)),
                                          nrow(FAFTD.effect.annot)-nrow(FAFTD.effect.annot%>%filter(wi.eBH<0.1))))

get_chi_LD(aldex_summhitsLD_stat,"NAFA") #X-squared = 4.4864, df = NA, p-value = 0.03898
get_chi_LD(aldex_summhitsLD_stat,"NAFT") #X-squared = 41.471, df = NA, p-value = 0.0004998
get_chi_LD(aldex_summhitsLD_stat,"FAFT") #X-squared = 47.982, df = NA, p-value = 0.0004998
