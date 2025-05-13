setwd("~/Notebooks/sfloresr/TRF-metaT/")

library(tidyverse)
library(data.table)
library(ALDEx2)

###########################################################
#paths
dat_nn<-"data/pfam_metaG/pfam_clean_noNT.txt"
dat_metadata<-"data/metaG_metadata_noNT.txt"
dat_annot<-"data/pfam_metaG/pfam_annotationkey.csv"
dat_path<-"data/pfam_metaG/"
res_path<-"data/DE_analysis/aldex_metaG/SFR23_0606_"
###########################################################
#functions

run_aldex<-function(dat,condsAB){
  AB.ald<-aldex.clr(round(dat),condsAB, mc.samples=500, denom="all", verbose=F)
  AB.ttest<-aldex.ttest(AB.ald,paired.test = FALSE, hist.plot=FALSE)%>%rownames_to_column("FeatureID")
  AB.effect<-aldex.effect(AB.ald)%>%rownames_to_column("FeatureID")%>%left_join(.,AB.ttest,by="FeatureID")
  
  AB.effect.annot<-AB.effect %>% left_join(.,annot, by ="FeatureID") %>%
    filter(!grepl("DUF",Name))%>%
    mutate(diffexpr=ifelse(wi.eBH<0.1 & diff.btw>0, "up",
                           ifelse(wi.eBH<0.1 & diff.btw< 0,"down","none")))%>%
    dplyr::select(FeatureID,Name, everything())
  return(AB.effect.annot)
}
###########################################################
#create pairwise files

#FA vs FT
pfamFAFT<-fread(dat_nn)%>%dplyr::select(c(1:25))
write.table(pfamFAFT,paste0(dat_path,"pfam_FAFT.txt"),sep = "\t", row.names = FALSE, quote=FALSE)

#just light
pfamFAFTL<-fread(dat_nn)%>%dplyr::select(c(1,8:13,20:25))
write.table(pfamFAFTL,paste0(dat_path,"pfam_FAFTL.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

#just dark
pfamFAFTD<-fread(dat_nn)%>%dplyr::select(c(1:7,14:19))
write.table(pfamFAFTD,paste0(dat_path,"pfam_FAFTD.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

#FT vs NA
pfamFTNA<-fread(dat_nn)%>%dplyr::select(c(1,14:43))
write.table(pfamFTNA,paste0(dat_path,"pfam_FTNA.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

#light
pfamFTNAL<-fread(dat_nn)%>%dplyr::select(c(1,20:25,35:43))
write.table(pfamFTNAL,paste0(dat_path,"pfam_FTNAL.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

#dark
pfamFTNAD<-fread(dat_nn)%>%dplyr::select(c(1,14:19,26:34))
write.table(pfamFTNAD,paste0(dat_path,"pfam_FTNAD.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

#FA vs NA
pfamFANA<-fread(dat_nn)%>%dplyr::select(c(1:13,26:43))
write.table(pfamFANA,paste0(dat_path,"pfam_FANA.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

#light
pfamFANAL<-fread(dat_nn)%>%dplyr::select(c(1,8:13,35:43))
write.table(pfamFANAL,paste0(dat_path,"pfam_FANAL.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

#dark
pfamFANAD<-fread(dat_nn)%>%dplyr::select(c(1:7,26:34))
write.table(pfamFANAD,paste0(dat_path,"pfam_FANAD.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

###########################################################
#run aldex

annot<-fread(dat_annot)
md<-fread(dat_metadata)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  arrange(condition,zt_time)
conds <- md$condition

##FA vs. FT
pfamFAFT<-fread(paste0(dat_path,"pfam_FAFT.txt"))%>%column_to_rownames("FeatureID")
condsFAvsFT<-conds[1:24]
FAFT.effect.annot<-run_aldex(pfamFAFT,condsFAvsFT)
write.table(FAFT.effect.annot,paste0(res_path,"FAFT_ald_effectwpval_wannot.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

FAFT.effect.annot %>% group_by(diffexpr) %>% tally(sort = TRUE)

# diffexpr     n
# <chr>    <int>
#   1 none      3096

#FA vs FT (light)
md<-fread(md_file)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition)) %>%
  filter(zt_time<13)%>%
  arrange(condition,zt_time)
conds <- md$condition

condsFAvsFTL<-conds[1:12]
pfamFAFTL<-fread(paste0(dat_path,"pfam_FAFTL.txt"))%>%column_to_rownames("FeatureID")
FAFTL.effect.annot<-run_aldex(pfamFAFTL,condsFAvsFTL)
write.table(FAFTL.effect.annot,paste0(res_path,"FAFTL_ald_effectwpval_wannot.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

FAFTL.effect.annot %>% group_by(diffexpr) %>% tally(sort = TRUE)

# diffexpr     n
# <chr>    <int>
#   1 none      2858

#FA vs FT (dark)
md<-fread(md_file)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition)) %>%
  filter(zt_time>9)%>%
  arrange(condition,zt_time)
conds <- md$condition

condsFAvsFTD<-conds[1:12]
pfamFAFTD<-fread(paste0(dat_path,"pfam_FAFTD.txt"))%>%column_to_rownames("FeatureID")
FAFTD.effect.annot<-run_aldex(pfamFAFTD,condsFAvsFTD)
write.table(FAFTD.effect.annot,paste0(res_path,"FAFTD_ald_effectwpval_wannot.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

FAFTD.effect.annot %>% group_by(diffexpr) %>% tally(sort = TRUE)

# diffexpr     n
# <chr>    <int>
#   1 none      2885

##NA vs FT
md<-fread(md_file)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  arrange(condition,zt_time)
conds <- md$condition

pfamFTNA<-fread(paste0(dat_path,"pfam_FTNA.txt"))%>%column_to_rownames("FeatureID")
condsFTvsNA<-conds[c(13:42)]
FTNA.effect.annot<-run_aldex(pfamFTNA,condsFTvsNA)
write.table(FTNA.effect.annot,paste0(res_path,"FTNA_ald_effectwpval_wannot.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

FTNA.effect.annot %>% group_by(diffexpr) %>% tally(sort = TRUE)

# diffexpr     n
# <chr>    <int>
#   1 none      2922
# 2 up         369
# 3 down       137

##NA vs FT (light)
md<-fread(md_file)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition)) %>%
  filter(zt_time<13)%>%
  arrange(condition,zt_time)
conds <- md$condition

condsFTvsNAL<-conds[c(7:21)]
pfamFTNAL<-fread(paste0(dat_path,"pfam_FTNAL.txt"))%>%column_to_rownames("FeatureID")
FTNAL.effect.annot<-run_aldex(pfamFTNAL,condsFTvsNAL)
write.table(FTNAL.effect.annot,paste0(res_path,"FTNAL_ald_effectwpval_wannot.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

FTNAL.effect.annot %>% group_by(diffexpr) %>% tally(sort = TRUE)

# diffexpr     n
# <chr>    <int>
#   1 none      3148
# 2 up         139
# 3 down        23

#NA vs FT (dark)
md<-fread(md_file)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition)) %>%
  filter(zt_time>9)%>%
  arrange(condition,zt_time)
conds <- md$condition

condsFTvsNAD<-conds[c(7:21)]
pfamFTNAD<-fread(paste0(dat_path,"pfam_FTNAD.txt"))%>%column_to_rownames("FeatureID")
FTNAD.effect.annot<-run_aldex(pfamFTNAD,condsFTvsNAD)
write.table(FTNAD.effect.annot,paste0(res_path,"FTNAD_ald_effectwpval_wannot.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

FTNAD.effect.annot %>% group_by(diffexpr) %>% tally(sort = TRUE)

# diffexpr     n
# <chr>    <int>
#   1 none      3059
# 2 up          38
# 3 down        15

##NA vs FA
md<-fread(md_file)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  arrange(condition,zt_time)
conds <- md$condition

pfamFANA<-fread(paste0(dat_path,"pfam_FANA.txt"))%>%column_to_rownames("FeatureID")
condsFAvsNA<-conds[c(1:12,25:42)]
FANA.effect.annot<-run_aldex(pfamFANA,condsFAvsNA)
write.table(FANA.effect.annot,paste0(res_path,"FANA_ald_effectwpval_wannot.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

FANA.effect.annot %>% group_by(diffexpr) %>% tally(sort = TRUE)

# diffexpr     n
# <chr>    <int>
#   1 none      2376
# 2 up         746
# 3 down       258

##NA vs FA (light)
md<-fread(md_file)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition)) %>%
  filter(zt_time<13)%>%
  arrange(condition,zt_time)
conds <- md$condition

condsFAvsNAL<-conds[c(1:6,13:21)]
pfamFANAL<-fread(paste0(dat_path,"pfam_FANAL.txt"))%>%column_to_rownames("FeatureID")
FANAL.effect.annot<-run_aldex(pfamFANAL,condsFAvsNAL)
write.table(FANAL.effect.annot,paste0(res_path,"FANAL_ald_effectwpval_wannot.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

FANAL.effect.annot %>% group_by(diffexpr) %>% tally(sort = TRUE)

# diffexpr     n
# <chr>    <int>
#   1 none      2802
# 2 up         387
# 3 down        85

#NA vs FA (dark)
md<-fread(md_file)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition)) %>%
  filter(zt_time>9)%>%
  arrange(condition,zt_time)
conds <- md$condition

condsFAvsNAD<-conds[c(1:6,13:21)]
pfamFANAD<-fread(paste0(dat_path,"pfam_FANAD.txt"))%>%column_to_rownames("FeatureID")
FANAD.effect.annot<-run_aldex(pfamFANAD,condsFAvsNAD)
write.table(FANAD.effect.annot,paste0(res_path,"FANAD_ald_effectwpval_wannot.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

FANAD.effect.annot %>% group_by(diffexpr) %>% tally(sort = TRUE)

# diffexpr     n
# <chr>    <int>
#   1 none      2680
# 2 up         337
# 3 down        82
###########################################################
#create summary of hits

#FANA
a_dw<-(nrow(FANA.effect.annot[FANA.effect.annot$diffexpr=="down"])/nrow(FANA.effect.annot))*100
a_up<-(nrow(FANA.effect.annot[FANA.effect.annot$diffexpr=="up"])/nrow(FANA.effect.annot))*100

#FTNA
b_dw<-(nrow(FTNA.effect.annot[FTNA.effect.annot$diffexpr=="down"])/nrow(FTNA.effect.annot))*100
b_up<-(nrow(FTNA.effect.annot[FTNA.effect.annot$diffexpr=="up"])/nrow(FTNA.effect.annot))*100

#FAFT
c_dw<-(nrow(FAFT.effect.annot[FAFT.effect.annot$diffexpr=="down"])/nrow(FAFT.effect.annot))*100
c_up<-(nrow(FAFT.effect.annot[FAFT.effect.annot$diffexpr=="up"])/nrow(FAFT.effect.annot))*100

aldex_summhits <- data.frame (comparison  = c("NAFA","NAFA","NAFT","NAFT","FAFT","FAFT"),
                              direction = c("down","up","down","up","down","up"),
                              perc= c(a_dw,a_up,b_dw,b_up,c_dw,c_up)) %>%
  group_by(comparison)%>%
  summarise(perc=sum(perc))%>%
  mutate(hits=c(nrow(FAFT.effect.annot[FAFT.effect.annot$diffexpr!="none"]),
                nrow(FANA.effect.annot[FANA.effect.annot$diffexpr!="none"]),
                nrow(FTNA.effect.annot[FTNA.effect.annot$diffexpr!="none"])),
         not_hits=c(nrow(FAFT.effect.annot[FAFT.effect.annot$diffexpr=="none"]),
                    nrow(FANA.effect.annot[FANA.effect.annot$diffexpr=="none"]),
                    nrow(FTNA.effect.annot[FTNA.effect.annot$diffexpr=="none"])),
         comparison=factor(comparison, levels = c("NAFA","NAFT","FAFT")))

write.table(aldex_summhits,paste0(res_path,"summaldexhits_BHless0.1.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

#Summarise hits for all three pairwise comparisons based on light and dark

#FANAL
al_dw<-(nrow(FANAL.effect.annot[FANAL.effect.annot$diffexpr=="down"])/nrow(FANAL.effect.annot))*100
al_up<-(nrow(FANAL.effect.annot[FANAL.effect.annot$diffexpr=="up"])/nrow(FANAL.effect.annot))*100

#FANAD
ad_dw<-(nrow(FANAD.effect.annot[FANAD.effect.annot$diffexpr=="down"])/nrow(FANAD.effect.annot))*100
ad_up<-(nrow(FANAD.effect.annot[FANAD.effect.annot$diffexpr=="up"])/nrow(FANAD.effect.annot))*100

#FTNAL
bl_dw<-(nrow(FTNAL.effect.annot[FTNAL.effect.annot$diffexpr=="down"])/nrow(FTNAL.effect.annot))*100
bl_up<-(nrow(FTNAL.effect.annot[FTNAL.effect.annot$diffexpr=="up"])/nrow(FTNAL.effect.annot))*100

#FTNAD
bd_dw<-(nrow(FTNAD.effect.annot[FTNAD.effect.annot$diffexpr=="down"])/nrow(FTNAD.effect.annot))*100
bd_up<-(nrow(FTNAD.effect.annot[FTNAD.effect.annot$diffexpr=="up"])/nrow(FTNAD.effect.annot))*100


#FAFTL
cl_dw<-(nrow(FAFTL.effect.annot[FAFTL.effect.annot$diffexpr=="down"])/nrow(FAFTL.effect.annot))*100
cl_up<-(nrow(FAFTL.effect.annot[FAFTL.effect.annot$diffexpr=="up"])/nrow(FAFTL.effect.annot))*100

#FAFTD
cd_dw<-(nrow(FAFTD.effect.annot[FAFTD.effect.annot$diffexpr=="down"])/nrow(FAFTD.effect.annot))*100
cd_up<-(nrow(FAFTD.effect.annot[FAFTD.effect.annot$diffexpr=="up"])/nrow(FAFTD.effect.annot))*100


aldex_summhitsLD <- data.frame (comparison  = c("NAFA","NAFA","NAFT","NAFT","FAFT","FAFT","NAFA","NAFA","NAFT","NAFT","FAFT","FAFT"),
                              direction = c("down","up","down","up","down","up","down","up","down","up","down","up"),
                              phase= c("light","light","light","light","light","light","dark","dark","dark","dark","dark","dark"),
                              perc= c(al_dw,al_up,bl_dw,bl_up,cl_dw,cl_up,-ad_dw,-ad_up,-bd_dw,-bd_up,-cd_dw,-cd_up)) %>%
  group_by(comparison,phase)%>%
  summarise(perc=sum(perc))%>%
  ungroup()%>%
  mutate(hits=c(nrow(FAFTD.effect.annot[FAFTD.effect.annot$diffexpr!="none"]),
                nrow(FAFTL.effect.annot[FAFTL.effect.annot$diffexpr!="none"]),
                nrow(FANAD.effect.annot[FANAD.effect.annot$diffexpr!="none"]),
                nrow(FANAL.effect.annot[FANAL.effect.annot$diffexpr!="none"]),
                nrow(FTNAD.effect.annot[FTNAD.effect.annot$diffexpr!="none"]),
                nrow(FTNAL.effect.annot[FTNAL.effect.annot$diffexpr!="none"])),
         not_hits=c(nrow(FAFTD.effect.annot[FAFTD.effect.annot$diffexpr=="none"]),
                    nrow(FAFTL.effect.annot[FAFTL.effect.annot$diffexpr=="none"]),
                    nrow(FANAD.effect.annot[FANAD.effect.annot$diffexpr=="none"]),
                    nrow(FANAL.effect.annot[FANAL.effect.annot$diffexpr=="none"]),
                    nrow(FTNAD.effect.annot[FTNAD.effect.annot$diffexpr=="none"]),
                    nrow(FTNAL.effect.annot[FTNAL.effect.annot$diffexpr=="none"])),
         comparison=factor(comparison, levels = c("NAFA","NAFT","FAFT")),
         phase=factor(phase, levels = c("light","dark")))

write.table(aldex_summhitsLD,paste0(res_path,"summaldexhitsLD_BHless0.1.txt"),sep = "\t",row.names = FALSE, quote=FALSE)
