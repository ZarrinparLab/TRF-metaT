setwd("~/scratch/TRF_multiomics/metagenomic/woltka2_results/filtered_metaG/")

library(tidyverse)
library(data.table)
library("Biostrings")
library(ggrepel)
library(ggpubr)
library(ggbreak)
library(ALDEx2)
library(viridis)
library(RColorBrewer)
library(ggvenn)
library(gplots)
library(DESeq2)
library(qiime2R)
library(VennDiagram)
###########################################################
#inputs

pfam_mtb<-"pfam_notnorm/pfam_clean_noNT.txt"
md_file<-"~/scratch/TRF_multiomics/metagenomic/metaG_metadata_noNT.txt"
dir<-"pfam_notnorm"

###########################################################
#pairwise files

#FA vs FT
pfamFAFT<-fread(pfam_mtb)%>%dplyr::select(c(1:25))
write.table(pfamFAFT,paste0(dir,"/pfam_FAFT.txt"),sep = "\t", row.names = FALSE, quote=FALSE)

#just light
pfamFAFTL<-fread(pfam_mtb)%>%dplyr::select(c(1,8:13,20:25))
write.table(pfamFAFTL,paste0(dir,"/pfam_FAFTL.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

#just dark
pfamFAFTD<-fread(pfam_mtb)%>%dplyr::select(c(1:7,14:19))
write.table(pfamFAFTD,paste0(dir,"/pfam_FAFTD.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

#FT vs NA
pfamFTNA<-fread(pfam_mtb)%>%dplyr::select(c(1,14:43))
write.table(pfamFTNA,paste0(dir,"/pfam_FTNA.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

#light
pfamFTNAL<-fread(pfam_mtb)%>%dplyr::select(c(1,20:25,35:43))
write.table(pfamFTNAL,paste0(dir,"/pfam_FTNAL.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

#dark
pfamFTNAD<-fread(pfam_mtb)%>%dplyr::select(c(1,14:19,26:34))
write.table(pfamFTNAD,paste0(dir,"/pfam_FTNAD.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

#FA vs NA
pfamFANA<-fread(pfam_mtb)%>%dplyr::select(c(1:13,26:43))
write.table(pfamFANA,paste0(dir,"/pfam_FANA.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

#light
pfamFANAL<-fread(pfam_mtb)%>%dplyr::select(c(1,8:13,35:43))
write.table(pfamFANAL,paste0(dir,"/pfam_FANAL.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

#dark
pfamFANAD<-fread(pfam_mtb)%>%dplyr::select(c(1:7,26:34))
write.table(pfamFANAD,paste0(dir,"/pfam_FANAD.txt"),sep = "\t",row.names = FALSE, quote=FALSE)
###########################################################
#aldex
annot<-fread(paste0(dir,"/pfam_annotationkey.csv"))
md<-fread(md_file)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  arrange(condition,zt_time)
conds <- md$condition

dir2<-"pfam_notnorm/aldex/SFR23_0606_"

##FA vs. FT
pfamFAFT<-fread(paste0(dir,"/pfam_FAFT.txt"))%>%column_to_rownames("FeatureID")
condsFAvsFT<-conds[1:24]
FAFT.ald<-aldex.clr(round(pfamFAFT),condsFAvsFT, mc.samples=500, denom="all", verbose=F)
FAFT.ttest<-aldex.ttest(FAFT.ald,paired.test = FALSE, hist.plot=FALSE)%>%rownames_to_column("FeatureID")
FAFT.effect<-aldex.effect(FAFT.ald)%>%rownames_to_column("FeatureID")%>%left_join(.,FAFT.ttest,by="FeatureID")
write.table(FAFT.effect,paste0(dir2,"FAFT_ald_effectwpval.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

FAFT.effect<-fread(paste0(dir2,"FAFT_ald_effectwpval.txt"))

FAFT.effect.annot<-FAFT.effect %>% left_join(.,annot, by ="FeatureID") %>%
  filter(!grepl("DUF",Name))%>%
  mutate(diffexpr=ifelse(wi.eBH<0.1 & diff.btw>0, "up",
                         ifelse(wi.eBH<0.1 & diff.btw< -0,"down","none")))

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

condsFAvsFT<-conds[1:12]
pfamFAFT<-fread(paste0(dir,"/pfam_FAFTL.txt"))%>%column_to_rownames("FeatureID")
FAFTL.ald<-aldex.clr(round(pfamFAFT),condsFAvsFT, mc.samples=500, denom="all", verbose=F)

FAFTL.ttest<-aldex.ttest(FAFTL.ald,paired.test = FALSE, hist.plot=FALSE)%>%rownames_to_column("FeatureID")
FAFTL.effect<-aldex.effect(FAFTL.ald)%>%rownames_to_column("FeatureID")%>%left_join(.,FAFTL.ttest,by="FeatureID")
write.table(FAFTL.effect,paste0(dir2,"FAFTL_ald_effectwpval.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

FAFTL.effect<-fread(paste0(dir2,"FAFTL_ald_effectwpval.txt"))

FAFTL.effect.annot<-FAFTL.effect %>% left_join(.,annot, by ="FeatureID") %>%
  filter(!grepl("DUF",Name))%>%
  mutate(diffexpr=ifelse(wi.eBH<0.1 & diff.btw>0, "up",
                         ifelse(wi.eBH<0.1 & diff.btw< -0,"down","none")))

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

condsFAvsFT<-conds[1:12]
pfamFAFT<-fread(paste0(dir,"/pfam_FAFTD.txt"))%>%column_to_rownames("FeatureID")
FAFTD.ald<-aldex.clr(round(pfamFAFT),condsFAvsFT, mc.samples=500, denom="all", verbose=F)

FAFTD.ttest<-aldex.ttest(FAFTD.ald,paired.test = FALSE, hist.plot=FALSE)%>%rownames_to_column("FeatureID")
FAFTD.effect<-aldex.effect(FAFTD.ald)%>%rownames_to_column("FeatureID")%>%left_join(.,FAFTD.ttest,by="FeatureID")
write.table(FAFTD.effect,paste0(dir2,"FAFTD_ald_effectwpval.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

FAFTD.effect<-fread(paste0(dir2,"FAFTD_ald_effectwpval.txt"))
FAFTD.effect.annot<-FAFTD.effect %>% left_join(.,annot, by ="FeatureID") %>%
  filter(!grepl("DUF",Name))%>%
  mutate(diffexpr=ifelse(wi.eBH<0.1 & diff.btw>0, "up",
                          ifelse(wi.eBH<0.1 & diff.btw< -0,"down","none")))

FAFTD.effect.annot %>% group_by(diffexpr) %>% tally(sort = TRUE)

# diffexpr     n
# <chr>    <int>
#   1 none      2885

##NA vs FT
md<-fread(md_file)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  arrange(condition,zt_time)
conds <- md$condition

pfamFTNA<-fread(paste0(dir,"/pfam_FTNA.txt"))%>%column_to_rownames("FeatureID")
condsFTvsNA<-conds[c(13:42)]
FTNA.ald<-aldex.clr(round(pfamFTNA),condsFTvsNA, mc.samples=500, denom="all", verbose=F)
FTNA.ttest<-aldex.ttest(FTNA.ald,paired.test = FALSE, hist.plot=FALSE)%>%rownames_to_column("FeatureID")
FTNA.effect<-aldex.effect(FTNA.ald)%>%rownames_to_column("FeatureID")%>%left_join(.,FTNA.ttest,by="FeatureID")%>%
  mutate(diff.btw=-(diff.btw))
write.table(FTNA.effect,paste0(dir2,"FTNA_ald_effectwpval.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

FTNA.effect<-fread(paste0(dir2,"FTNA_ald_effectwpval.txt"))

FTNA.effect.annot<-FTNA.effect %>% left_join(.,annot, by ="FeatureID") %>%
  filter(!grepl("DUF",Name))%>%
  mutate(diffexpr=ifelse(wi.eBH<0.1 & diff.btw>0, "up",
                         ifelse(wi.eBH<0.1 & diff.btw< -0,"down","none")))

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

condsFTvsNA<-conds[c(7:21)]
pfamFTNA<-fread(paste0(dir,"/pfam_FTNAL.txt"))%>%column_to_rownames("FeatureID")
FTNAL.ald<-aldex.clr(round(pfamFTNA),condsFTvsNA, mc.samples=500, denom="all", verbose=F)

FTNAL.ttest<-aldex.ttest(FTNAL.ald,paired.test = FALSE, hist.plot=FALSE)%>%rownames_to_column("FeatureID")
FTNAL.effect<-aldex.effect(FTNAL.ald)%>%rownames_to_column("FeatureID")%>%left_join(.,FTNAL.ttest,by="FeatureID")%>%
  mutate(diff.btw=-(diff.btw))
write.table(FTNAL.effect,paste0(dir2,"FTNAL_ald_effectwpval.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

FTNAL.effect<-fread(paste0(dir2,"FTNAL_ald_effectwpval.txt"))

FTNAL.effect.annot<-FTNAL.effect %>% left_join(.,annot, by ="FeatureID") %>%
  filter(!grepl("DUF",Name))%>%
  mutate(diffexpr=ifelse(wi.eBH<0.1 & diff.btw>0, "up",
                         ifelse(wi.eBH<0.1 & diff.btw< -0,"down","none")))

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

condsFTvsNA<-conds[c(7:21)]
pfamFTNA<-fread(paste0(dir,"/pfam_FTNAD.txt"))%>%column_to_rownames("FeatureID")
FTNAD.ald<-aldex.clr(round(pfamFTNA),condsFTvsNA, mc.samples=500, denom="all", verbose=F)

FTNAD.ttest<-aldex.ttest(FTNAD.ald,paired.test = FALSE, hist.plot=FALSE)%>%rownames_to_column("FeatureID")
FTNAD.effect<-aldex.effect(FTNAD.ald)%>%rownames_to_column("FeatureID")%>%left_join(.,FTNAD.ttest,by="FeatureID")%>%
  mutate(diff.btw=-(diff.btw))
write.table(FTNAD.effect,paste0(dir2,"FTNAD_ald_effectwpval.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

FTNAD.effect<-fread(paste0(dir2,"FTNAD_ald_effectwpval.txt"))

FTNAD.effect.annot<-FTNAD.effect %>% left_join(.,annot, by ="FeatureID") %>%
  filter(!grepl("DUF",Name))%>%
  mutate(diffexpr=ifelse(wi.eBH<0.1 & diff.btw>0, "up",
                         ifelse(wi.eBH<0.1 & diff.btw< -0,"down","none")))

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

pfamFANA<-fread(paste0(dir,"/pfam_FANA.txt"))%>%column_to_rownames("FeatureID")
condsFAvsNA<-conds[c(1:12,25:42)]
FANA.ald<-aldex.clr(round(pfamFANA),condsFAvsNA, mc.samples=500, denom="all", verbose=F)
FANA.ttest<-aldex.ttest(FANA.ald,paired.test = FALSE, hist.plot=FALSE)%>%rownames_to_column("FeatureID")
FANA.effect<-aldex.effect(FANA.ald)%>%rownames_to_column("FeatureID")%>%left_join(.,FANA.ttest,by="FeatureID")%>%
  mutate(diff.btw=-(diff.btw))
write.table(FANA.effect,paste0(dir2,"FANA_ald_effectwpval.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

FANA.effect<-fread(paste0(dir2,"FANA_ald_effectwpval.txt"))

FANA.effect.annot<-FANA.effect %>% left_join(.,annot, by ="FeatureID") %>%
  filter(!grepl("DUF",Name))%>%
  mutate(diffexpr=ifelse(wi.eBH<0.1 & diff.btw>0, "up",
                         ifelse(wi.eBH<0.1 & diff.btw< -0,"down","none")))

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

condsFAvsNA<-conds[c(1:6,13:21)]
pfamFANA<-fread(paste0(dir,"/pfam_FANAL.txt"))%>%column_to_rownames("FeatureID")
FANAL.ald<-aldex.clr(round(pfamFANA),condsFAvsNA, mc.samples=500, denom="all", verbose=F)

FANAL.ttest<-aldex.ttest(FANAL.ald,paired.test = FALSE, hist.plot=FALSE)%>%rownames_to_column("FeatureID")
FANAL.effect<-aldex.effect(FANAL.ald)%>%rownames_to_column("FeatureID")%>%left_join(.,FANAL.ttest,by="FeatureID")%>%
  mutate(diff.btw=-(diff.btw))
write.table(FANAL.effect,paste0(dir2,"FANAL_ald_effectwpval.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

FANAL.effect<-fread(paste0(dir2,"FANAL_ald_effectwpval.txt"))

FANAL.effect.annot<-FANAL.effect %>% left_join(.,annot, by ="FeatureID") %>%
  filter(!grepl("DUF",Name))%>%
  mutate(diffexpr=ifelse(wi.eBH<0.1 & diff.btw>0, "up",
                         ifelse(wi.eBH<0.1 & diff.btw< -0,"down","none")))

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

condsFAvsNA<-conds[c(1:6,13:21)]
pfamFANA<-fread(paste0(dir,"/pfam_FANAD.txt"))%>%column_to_rownames("FeatureID")
FANAD.ald<-aldex.clr(round(pfamFANA),condsFAvsNA, mc.samples=500, denom="all", verbose=F)

FANAD.ttest<-aldex.ttest(FANAD.ald,paired.test = FALSE, hist.plot=FALSE)%>%rownames_to_column("FeatureID")
FANAD.effect<-aldex.effect(FANAD.ald)%>%rownames_to_column("FeatureID")%>%left_join(.,FANAD.ttest,by="FeatureID")%>%
  mutate(diff.btw=-(diff.btw))
write.table(FANAD.effect,paste0(dir2,"FANAD_ald_effectwpval.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

FANAD.effect<-fread(paste0(dir2,"FANAD_ald_effectwpval.txt"))

FANAD.effect.annot<-FANAD.effect %>% left_join(.,annot, by ="FeatureID") %>%
  filter(!grepl("DUF",Name))%>%
  mutate(diffexpr=ifelse(wi.eBH<0.1 & diff.btw>0, "up",
                         ifelse(wi.eBH<0.1 & diff.btw< -0,"down","none")))

FANAD.effect.annot %>% group_by(diffexpr) %>% tally(sort = TRUE)

# diffexpr     n
# <chr>    <int>
#   1 none      2680
# 2 up         337
# 3 down        82
###########################################################

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

write.table(aldex_summhits,paste0(dir2,"summaldexhits_BHless0.1.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

ggplot(aldex_summhits, aes(x=comparison, y=perc)) + 
  geom_bar(stat="identity", position="identity",fill="gray30") +
  theme_minimal() +scale_y_continuous(limits=c(0,30))

ggsave(paste0(dir2,"summaldexhits_BHless0.1.pdf"), width = 2.2, height = 3)

#Summarise hits for all three pairwise comparisons based on light and dark (all-not just unique)

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

write.table(aldex_summhitsLD,paste0(dir2,"summaldexhitsLD_BHless0.1.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

ggplot(aldex_summhitsLD, aes(x=comparison, y=perc,fill=phase)) + 
  geom_bar(stat="identity", position="identity") +
  scale_fill_manual(values=c("gray70","gray10")) +
  theme_minimal() + scale_y_continuous(limits=c(-22,22))+
  theme(legend.position = "top")

ggsave(paste0(dir2,"summaldexhitsLD_BHless0.1.pdf"), width = 2.2, height = 3)

###########################################################
#plot summary of the sig hits, just unique, no overlap

list_venn <- list(FAFT = (FAFT.effect.annot%>%filter(wi.eBH<0.1))$FeatureID,
                  FANA = (FANA.effect.annot%>%filter(wi.eBH<0.1))$FeatureID,
                  FTNA = (FTNA.effect.annot%>%filter(wi.eBH<0.1))$FeatureID)

ItemsList <- venn(list_venn, show.plot = FALSE)
all<-attributes(ItemsList)$intersections

#FANA
a<-(length(all$FANA)/nrow(FANA.effect.annot))*100

#FTNA
b<-(length(all$FTNA)/nrow(FTNA.effect.annot))*100

#FAFT
c<-(length(all$FAFT)/nrow(FAFT.effect.annot))*100

aldex_summhits_unique <- data.frame (comparison  = c("NAFA","NAFT","FAFT"),
                              perc= c(a,b,c))%>%
  mutate(comparison=factor(comparison, levels = c("NAFA","NAFT","FAFT")))

write.table(aldex_summhits_unique,paste0(dir2,"summaldexhitsUQ_BHless0.1.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

ggplot(aldex_summhits_unique, aes(x=comparison, y=perc)) + 
  geom_bar(stat="identity", position="identity",fill="gray30") +
  theme_minimal()+scale_y_continuous(limits=c(0,25))

ggsave(paste0(dir2,"summaldexhitsUQ_BHless0.1.pdf"), width = 2.2, height = 3)

list_venn <- list(FAFT = (FAFTL.effect.annot%>%filter(wi.eBH<0.1))$FeatureID,
                  FANA = (FANAL.effect.annot%>%filter(wi.eBH<0.1))$FeatureID,
                  FTNA = (FTNAL.effect.annot%>%filter(wi.eBH<0.1))$FeatureID)

ItemsList <- venn(list_venn, show.plot = FALSE)
allL<-attributes(ItemsList)$intersections

#FANAL
al<-(length(allL$FANA)/nrow(FANAL.effect.annot))*100

#FTNAL
bl<-(length(allL$FTNA)/nrow(FTNAL.effect.annot))*100

#FAFTL
cl<-(length(allL$FAFT)/nrow(FAFTL.effect.annot))*100

list_venn <- list(FAFT = (FAFTD.effect.annot%>%filter(wi.eBH<0.1))$FeatureID,
                  FANA = (FANAD.effect.annot%>%filter(wi.eBH<0.1))$FeatureID,
                  FTNA = (FTNAD.effect.annot%>%filter(wi.eBH<0.1))$FeatureID)

ItemsList <- venn(list_venn, show.plot = FALSE)
allD<-attributes(ItemsList)$intersections

#FANAD
ad<-(length(allD$FANA)/nrow(FANAD.effect.annot))*100

#FTNAD
bd<-(length(allD$FTNA)/nrow(FTNAD.effect.annot))*100

#FAFTD
cd<-(length(allD$FAFT)/nrow(FAFTD.effect.annot))*100

aldex_summhitsLD <- data.frame (comparison  = c("NAFA","NAFT","FAFT","NAFA","NAFT","FAFT"),
                                phase= c("light","light","light","dark","dark","dark"),
                                perc= c(al,bl,cl,
                                        -ad,-bd,-cd))%>%
  mutate(comparison=factor(comparison, levels = c("NAFA","NAFT","FAFT")),
         phase=factor(phase, levels = c("light","dark")))

write.table(aldex_summhitsLD,paste0(dir2,"summaldexhitsLD_unique_BHless0.1.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

ggplot(aldex_summhitsLD, aes(x=comparison, y=perc,fill=phase)) + 
  geom_bar(stat="identity", position="identity") +
  scale_fill_manual(values=c("gray70","gray10")) +
  theme_minimal() +scale_y_continuous(limits=c(-12,12))+
  theme(legend.position = "top")

ggsave(paste0(dir2,"summaldexhitsLD_unique_BHless0.1.pdf"), width = 2.2, height = 3)

###########################################################
#unique hits/venn diagram
list_venn <- list(FAFT = (FAFT.effect.annot%>%filter(wi.eBH<0.1))$FeatureID,
                  FANA = (FANA.effect.annot%>%filter(wi.eBH<0.1))$FeatureID,
                  FTNA = (FTNA.effect.annot%>%filter(wi.eBH<0.1))$FeatureID)

ItemsList <- venn(list_venn, show.plot = FALSE)
all<-attributes(ItemsList)$intersections
unique_hits_FAFT<-FAFT.effect.annot%>%filter(FeatureID %in% all$FAFT)
unique_hits_FTNA<-FTNA.effect.annot%>%filter(FeatureID %in% all$FTNA)
unique_hits_FANA<-FANA.effect.annot%>%filter(FeatureID %in% all$FANA)
write.table(unique_hits_FAFT,paste0(dir2,"FAFT_",length(all$FAFT),"uniquehits_BHless0.1.txt"),sep = "\t",row.names = FALSE, quote=FALSE)
write.table(unique_hits_FTNA,paste0(dir2,"FTNA_",length(all$FTNA),"uniquehits_BHless0.1.txt"),sep = "\t",row.names = FALSE, quote=FALSE)
write.table(unique_hits_FANA,paste0(dir2,"FANA_",length(all$FANA),"uniquehits_BHless0.1.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

list_venn <- list(FAFT = (FAFTL.effect.annot%>%filter(wi.eBH<0.1))$FeatureID,
                  FANA = (FANAL.effect.annot%>%filter(wi.eBH<0.1))$FeatureID,
                  FTNA = (FTNAL.effect.annot%>%filter(wi.eBH<0.1))$FeatureID)

ItemsList <- venn(list_venn, show.plot = FALSE)
light<-attributes(ItemsList)$intersections
unique_hits_light<-FAFTL.effect.annot%>%filter(FeatureID %in% light$FAFT)
unique_hits_lightFANA<-FANAL.effect.annot%>%filter(FeatureID %in% light$FANA)
unique_hits_lightFTNA<-FTNAL.effect.annot%>%filter(FeatureID %in% light$FTNA)
write.table(unique_hits_light,paste0(dir2,"FAFTL_",length(light$FAFT),"uniquehits_BHless0.1.txt"),sep = "\t",row.names = FALSE, quote=FALSE)
write.table(unique_hits_lightFTNA,paste0(dir2,"FTNAL_",length(light$FTNA),"uniquehits_BHless0.1.txt"),sep = "\t",row.names = FALSE, quote=FALSE)
write.table(unique_hits_lightFANA,paste0(dir2,"FANAL_",length(light$FANA),"uniquehits_BHless0.1.txt"),sep = "\t",row.names = FALSE, quote=FALSE)


list_venn <- list(FAFT = (FAFTD.effect.annot%>%filter(wi.eBH<0.1))$FeatureID,
                  FANA = (FANAD.effect.annot%>%filter(wi.eBH<0.1))$FeatureID,
                  FTNA = (FTNAD.effect.annot%>%filter(wi.eBH<0.1))$FeatureID)

ItemsList <- venn(list_venn, show.plot = FALSE)
dark<-attributes(ItemsList)$intersections
unique_hits_dark<-FAFTD.effect.annot%>%filter(FeatureID %in% dark$FAFT)
unique_hits_darkFANA<-FANAD.effect.annot%>%filter(FeatureID %in% dark$FANA)
unique_hits_darkFTNA<-FTNAD.effect.annot%>%filter(FeatureID %in% dark$FTNA)
write.table(unique_hits_dark,paste0(dir2,"FAFTD_",length(dark$FAFT),"uniquehits_BHless0.1.txt"),sep = "\t",row.names = FALSE, quote=FALSE)
write.table(unique_hits_darkFTNA,paste0(dir2,"FTNAD_",length(dark$FTNA),"uniquehits_BHless0.1.txt"),sep = "\t",row.names = FALSE, quote=FALSE)
write.table(unique_hits_darkFANA,paste0(dir2,"FANAD_",length(dark$FANA),"uniquehits_BHless0.1.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

# list_venn <- list(light = unique_hits_light$FeatureID,
#                   dark = unique_hits_dark$FeatureID)
# 
# p<-venn.diagram(list_venn, fill = c("white", "black"),height = 10,
#                  width = 10,alpha = c(0.5, 0.5), lwd =0,filename = NULL)
# 
# grid.draw(p)
# pdf(file=paste0(dir2,"FAFT_LDintersect_BHless0.05.pdf"))
# grid.draw(p)
# dev.off()

list_venn <- list(light = unique_hits_lightFTNA$FeatureID,
                  dark = unique_hits_darkFTNA$FeatureID)

p<-venn.diagram(list_venn, fill = c("white", "black"),height = 10,
                width = 10,alpha = c(0.5, 0.5), lwd =0,filename = NULL)

grid.draw(p)
pdf(file=paste0(dir2,"FTNA_LDintersect_BHless0.05.pdf"))
grid.draw(p)
dev.off()

list_venn <- list(light = unique_hits_lightFANA$FeatureID,
                  dark = unique_hits_darkFANA$FeatureID)

p<-venn.diagram(list_venn, fill = c("white", "black"),height = 10,
                width = 10,alpha = c(0.5, 0.5), lwd =0,filename = NULL)

grid.draw(p)
pdf(file=paste0(dir2,"FANA_LDintersect_BHless0.05.pdf"))
grid.draw(p)
dev.off()
