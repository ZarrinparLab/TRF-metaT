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
#files for light vs dark comparison

#FA
pfamFA<-fread(pfam_mtb)%>%dplyr::select(c(1,8:13,2:7))
write.table(pfamFA,paste0(dir,"/pfam_FA.txt"),sep = "\t", row.names = FALSE, quote=FALSE)

#FT vs NA
pfamFT<-fread(pfam_mtb)%>%dplyr::select(c(1,20:25,14:19))
write.table(pfamFT,paste0(dir,"/pfam_FT.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

#FA vs NA
pfamNA<-fread(pfam_mtb)%>%dplyr::select(c(1,35:43,26:34))
write.table(pfamNA,paste0(dir,"/pfam_NA.txt"),sep = "\t",row.names = FALSE, quote=FALSE)
###########################################################

#aldex
annot<-fread(paste0(dir,"/pfam_annotationkey.csv"))
md<-fread(md_file)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  arrange(condition,zt_time)
phase <- md$lightdark

dir2<-"pfam_notnorm/aldex/SFR23_0623_"

##FA light vs. dark
pfamFA<-fread(paste0(dir,"/pfam_FA.txt"))%>%column_to_rownames("FeatureID")
phaseFA<-phase[1:12]
FA.ald<-aldex.clr(round(pfamFA),phaseFA, mc.samples=500, denom="all", verbose=F)
FA.ttest<-aldex.ttest(FA.ald,paired.test = FALSE, hist.plot=FALSE)%>%rownames_to_column("FeatureID")
FA.effect<-aldex.effect(FA.ald)%>%rownames_to_column("FeatureID")%>%left_join(.,FA.ttest,by="FeatureID")%>%
  mutate(diff.btw=-(diff.btw))
write.table(FA.effect,paste0(dir2,"FA_LD_ald_effectwpval.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

FA.effect<-fread(paste0(dir2,"FA_LD_ald_effectwpval.txt"))

FA.effect.annot<-FA.effect %>% left_join(.,annot, by ="FeatureID") %>%
  filter(!grepl("DUF",Name))%>%
  mutate(diffexpr=ifelse(wi.eBH<0.1 & diff.btw>0, "up",
                         ifelse(wi.eBH<0.1 & diff.btw< -0,"down","none")))

FA.effect.annot %>% group_by(diffexpr) %>% tally(sort = TRUE)

# diffexpr     n
# <chr>    <int>
#   1 none      2791

##FT light dark
md<-fread(md_file)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  arrange(condition,zt_time)
phase <- md$lightdark

pfamFT<-fread(paste0(dir,"/pfam_FT.txt"))%>%column_to_rownames("FeatureID")
phaseFT<-phase[c(13:24)]
FT.ald<-aldex.clr(round(pfamFT),phaseFT, mc.samples=500, denom="all", verbose=F)
FT.ttest<-aldex.ttest(FT.ald,paired.test = FALSE, hist.plot=FALSE)%>%rownames_to_column("FeatureID")
FT.effect<-aldex.effect(FT.ald)%>%rownames_to_column("FeatureID")%>%left_join(.,FT.ttest,by="FeatureID")%>%
  mutate(diff.btw=-(diff.btw))
write.table(FT.effect,paste0(dir2,"FT_LD_ald_effectwpval.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

FT.effect<-fread(paste0(dir2,"FT_LD_ald_effectwpval.txt"))

FT.effect.annot<-FT.effect %>% left_join(.,annot, by ="FeatureID") %>%
  filter(!grepl("DUF",Name))%>%
  mutate(diffexpr=ifelse(wi.eBH<0.1 & diff.btw>0, "up",
                         ifelse(wi.eBH<0.1 & diff.btw< -0,"down","none")))

FT.effect.annot %>% group_by(diffexpr) %>% tally(sort = TRUE)

# diffexpr     n
# <chr>    <int>
#   1 none      2886

##NA light vs. dark
md<-fread(md_file)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  arrange(condition,zt_time)
phase <- md$lightdark

pfamNA<-fread(paste0(dir,"/pfam_NA.txt"))%>%column_to_rownames("FeatureID")
phaseNA<-phase[c(25:42)]
NA.ald<-aldex.clr(round(pfamNA),phaseNA, mc.samples=500, denom="all", verbose=F)
NA.ttest<-aldex.ttest(NA.ald,paired.test = FALSE, hist.plot=FALSE)%>%rownames_to_column("FeatureID")
NA.effect<-aldex.effect(NA.ald)%>%rownames_to_column("FeatureID")%>%left_join(.,NA.ttest,by="FeatureID")%>%
  mutate(diff.btw=-(diff.btw))
write.table(NA.effect,paste0(dir2,"NA_LD_ald_effectwpval.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

NA.effect<-fread(paste0(dir2,"NA_LD_ald_effectwpval.txt"))

NA.effect.annot<-NA.effect %>% left_join(.,annot, by ="FeatureID") %>%
  filter(!grepl("DUF",Name))%>%
  mutate(diffexpr=ifelse(wi.eBH<0.1 & diff.btw>0, "up",
                         ifelse(wi.eBH<0.1 & diff.btw< -0,"down","none")))

NA.effect.annot %>% group_by(diffexpr) %>% tally(sort = TRUE)

# diffexpr     n
# <chr>    <int>
#   1 none      3257
###########################################################

#NA
a_dw<-(nrow(NA.effect.annot[NA.effect.annot$diffexpr=="down"])/nrow(NA.effect.annot))*100
a_up<-(nrow(NA.effect.annot[NA.effect.annot$diffexpr=="up"])/nrow(NA.effect.annot))*100

#FA
c_dw<-(nrow(FA.effect.annot[FA.effect.annot$diffexpr=="down"])/nrow(FA.effect.annot))*100
c_up<-(nrow(FA.effect.annot[FA.effect.annot$diffexpr=="up"])/nrow(FA.effect.annot))*100

#FT
b_dw<-(nrow(FT.effect.annot[FT.effect.annot$diffexpr=="down"])/nrow(FT.effect.annot))*100
b_up<-(nrow(FT.effect.annot[FT.effect.annot$diffexpr=="up"])/nrow(FT.effect.annot))*100

aldex_summhits <- data.frame (comparison  = c("NA","NA","FA","FA","FT","FT"),
                              direction = c("light","dark","light","dark","light","dark"),
                              perc= c(a_dw,a_up,c_dw,c_up,b_dw,b_up)) %>%
  group_by(comparison)%>%
  summarise(perc=sum(perc))%>%
  mutate(hits=c(nrow(FA.effect.annot[FA.effect.annot$diffexpr!="none"]),
                nrow(FT.effect.annot[FT.effect.annot$diffexpr!="none"]),
                nrow(NA.effect.annot[NA.effect.annot$diffexpr!="none"])),
         not_hits=c(nrow(FA.effect.annot[FA.effect.annot$diffexpr=="none"]),
                    nrow(FT.effect.annot[FT.effect.annot$diffexpr=="none"]),
                    nrow(NA.effect.annot[NA.effect.annot$diffexpr=="none"])),
         comparison=factor(comparison, levels = c("NA","FA","FT")))

write.table(aldex_summhits,paste0(dir2,"summaldexhits_justLD_BHless0.1.txt"),sep = "\t",row.names = FALSE, quote=FALSE)
