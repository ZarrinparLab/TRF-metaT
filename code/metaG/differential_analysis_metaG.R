setwd("~/scratch/metatranscript/metagenomic/woltka2_results")

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
#get the total read amt for PF04563.18 (rpob)
mtb<-fread("filtered_metaG/pfam/pfam_clean_noNT.txt")%>%
  filter(FeatureID=="PF04563.18")%>%
  column_to_rownames("FeatureID")%>%
  rowSums() #662 is the depth of RPOB

###########################################################
#FA vs FT
pfamFAFT<-fread("filtered_metaG/pfam/pfam_clean_noNT.txt")%>%dplyr::select(c(1:25))
write.table(pfamFAFT,"filtered_metaG/pfam/pfam_FAFT.txt",sep = "\t",row.names = FALSE, quote=FALSE)

mdFAFT<-fread("~/scratch/metatranscript/metagenomic/metaG_metadata_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  filter(condition!="NA")
write.table(mdFAFT,"metaG_FAFT_metadata.txt",sep = "\t",row.names = FALSE, quote=FALSE)

#just light
pfamFAFTL<-fread("filtered_metaG/pfam/pfam_clean_noNT.txt")%>%dplyr::select(c(1,8:13,20:25))
write.table(pfamFAFTL,"filtered_metaG/pfam/pfam_FAFTL.txt",sep = "\t",row.names = FALSE, quote=FALSE)

mdFAFTL<-fread("~/scratch/metatranscript/metagenomic/metaG_metadata_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  filter(condition!="NA" & lightdark=="light")
write.table(mdFAFTL,"metaG_FAFTL_metadata.txt",sep = "\t",row.names = FALSE, quote=FALSE)

#just dark
pfamFAFTD<-fread("filtered_metaG/pfam/pfam_clean_noNT.txt")%>%dplyr::select(c(1:7,14:19))
write.table(pfamFAFTD,"filtered_metaG/pfam/pfam_FAFTD.txt",sep = "\t",row.names = FALSE, quote=FALSE)

mdFAFTD<-fread("~/scratch/metatranscript/metagenomic/metaG_metadata_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  filter(condition!="NA" & lightdark=="dark")
write.table(mdFAFTD,"metaG_FAFTD_metadata.txt",sep = "\t",row.names = FALSE, quote=FALSE)

#FT vs NA
pfamFTNA<-fread("filtered_metaG/pfam/pfam_clean_noNT.txt")%>%dplyr::select(c(1,14:43))
write.table(pfamFTNA,"filtered_metaG/pfam/pfam_FTNA.txt",sep = "\t",row.names = FALSE, quote=FALSE)

mdFTNA<-fread("~/scratch/metatranscript/metagenomic/metaG_metadata_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  filter(condition!="FA")
write.table(mdFTNA,"metaG_FTNA_metadata.txt",sep = "\t",row.names = FALSE, quote=FALSE)

#light
pfamFTNAL<-fread("filtered_metaG/pfam/pfam_clean_noNT.txt")%>%dplyr::select(c(1,20:25,35:43))
write.table(pfamFTNAL,"filtered_metaG/pfam/pfam_FTNAL.txt",sep = "\t",row.names = FALSE, quote=FALSE)

mdFTNAL<-fread("~/scratch/metatranscript/metagenomic/metaG_metadata_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  filter(condition!="FA"& lightdark=="light")
write.table(mdFTNAL,"metaG_FTNAL_metadata.txt",sep = "\t",row.names = FALSE, quote=FALSE)

#dark
pfamFTNAD<-fread("filtered_metaG/pfam/pfam_clean_noNT.txt")%>%dplyr::select(c(1,14:19,26:34))
write.table(pfamFTNAD,"filtered_metaG/pfam/pfam_FTNAD.txt",sep = "\t",row.names = FALSE, quote=FALSE)

mdFTNAD<-fread("~/scratch/metatranscript/metagenomic/metaG_metadata_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  filter(condition!="FA"& lightdark=="dark")
write.table(mdFTNAD,"metaG_FTNAD_metadata.txt",sep = "\t",row.names = FALSE, quote=FALSE)

#FA vs NA
pfamFANA<-fread("filtered_metaG/pfam/pfam_clean_noNT.txt")%>%dplyr::select(c(1:13,26:43))
write.table(pfamFANA,"filtered_metaG/pfam/pfam_FANA.txt",sep = "\t",row.names = FALSE, quote=FALSE)

mdFANA<-fread("~/scratch/metatranscript/metagenomic/metaG_metadata_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  filter(condition!="FT")
write.table(mdFANA,"metaG_FANA_metadata.txt",sep = "\t",row.names = FALSE, quote=FALSE)

#light
pfamFANAL<-fread("filtered_metaG/pfam/pfam_clean_noNT.txt")%>%dplyr::select(c(1,8:13,35:43))
write.table(pfamFANAL,"filtered_metaG/pfam/pfam_FANAL.txt",sep = "\t",row.names = FALSE, quote=FALSE)

mdFANAL<-fread("~/scratch/metatranscript/metagenomic/metaG_metadata_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  filter(condition!="FT"& lightdark=="light")
write.table(mdFANAL,"metaG_FANAL_metadata.txt",sep = "\t",row.names = FALSE, quote=FALSE)

#dark
pfamFANAD<-fread("filtered_metaG/pfam/pfam_clean_noNT.txt")%>%dplyr::select(c(1:7,26:34))
write.table(pfamFANAD,"filtered_metaG/pfam/pfam_FANAD.txt",sep = "\t",row.names = FALSE, quote=FALSE)

mdFANAD<-fread("~/scratch/metatranscript/metagenomic/metaG_metadata_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  filter(condition!="FT"& lightdark=="dark")
write.table(mdFANAD,"metaG_FANAD_metadata.txt",sep = "\t",row.names = FALSE, quote=FALSE)

###########################################################
#aldex
annot<-fread("~/scratch/metatranscript/metatranscript/woltka2_results/pfam/pfam_annotationkey.csv")
md<-fread("~/scratch/metatranscript/metagenomic/metaG_metadata_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  arrange(condition,zt_time)
conds <- md$condition

##FA vs. FT
pfamFAFT<-fread("filtered_metaG/pfam/pfam_FAFT.txt")%>%column_to_rownames("FeatureID")
condsFAvsFT<-conds[1:24]
FAFT.ald<-aldex.clr(round(pfamFAFT),condsFAvsFT, mc.samples=500, denom="all", verbose=F)
FAFT.ttest<-aldex.ttest(FAFT.ald,paired.test = FALSE, hist.plot=FALSE)%>%rownames_to_column("FeatureID")
FAFT.effect<-aldex.effect(FAFT.ald)%>%rownames_to_column("FeatureID")%>%left_join(.,FAFT.ttest,by="FeatureID")
write.table(FAFT.effect,"filtered_metaG/aldex/SFR23_0523_FAFT_ald_effectwpval.txt",sep = "\t",row.names = FALSE, quote=FALSE)

FAFT.effect<-fread("filtered_metaG/aldex/SFR23_0523_FAFT_ald_effectwpval.txt")

FAFT.effect.annot<-FAFT.effect %>% left_join(.,annot, by ="FeatureID") %>%
  filter(!grepl("DUF",Name))%>%
  mutate(diffexpr=ifelse(wi.ep<0.05 & diff.btw>0, "up",
                         ifelse(wi.ep<0.05 & diff.btw< -0,"down","none"))) #%>%
  #filter(!str_detect(Name,"unknown"))

FAFT.effect.annot %>% group_by(diffexpr) %>% tally(sort = TRUE)

# diffexpr     n
# <chr>    <int>
#   1 none      3041
# 2 down        23
# 3 up           3

##FA vs FT (light)
md<-fread("~/scratch/metatranscript/metagenomic/metaG_metadata_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition)) %>%
  filter(zt_time<13)%>%
  arrange(condition,zt_time)
conds <- md$condition

condsFAvsFT<-conds[1:12]
pfamFAFT<-fread("filtered_metaG/pfam/pfam_FAFTL.txt")%>%column_to_rownames("FeatureID")
FAFTL.ald<-aldex.clr(round(pfamFAFT),condsFAvsFT, mc.samples=500, denom="all", verbose=F)

FAFTL.ttest<-aldex.ttest(FAFTL.ald,paired.test = FALSE, hist.plot=FALSE)%>%rownames_to_column("FeatureID")
FAFTL.effect<-aldex.effect(FAFTL.ald)%>%rownames_to_column("FeatureID")%>%left_join(.,FAFTL.ttest,by="FeatureID")
write.table(FAFTL.effect,"filtered_metaG/aldex/SFR23_0523_FAFTL_ald_effectwpval.txt",sep = "\t",row.names = FALSE, quote=FALSE)

FAFTL.effect<-fread("filtered_metaG/aldex/SFR23_0523_FAFTL_ald_effectwpval.txt")

FAFTL.effect.annot<-FAFTL.effect %>% left_join(.,annot, by ="FeatureID") %>%
  filter(!grepl("DUF",Name))%>%
  mutate(diffexpr=ifelse(wi.ep<0.05 & diff.btw>0, "up",
                         ifelse(wi.ep<0.05 & diff.btw< -0,"down","none"))) #%>%
  # filter(!str_detect(Name,"unknown"))

FAFTL.effect.annot %>% group_by(diffexpr) %>% tally(sort = TRUE)

# diffexpr     n
# <chr>    <int>
#   1 none      2816
# 2 down        34
# 3 up           9

#FA vs FT (dark)
md<-fread("~/scratch/metatranscript/metagenomic/metaG_metadata_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition)) %>%
  filter(zt_time>9)%>%
  arrange(condition,zt_time)
conds <- md$condition

condsFAvsFT<-conds[1:12]
pfamFAFT<-fread("filtered_metaG/pfam/pfam_FAFTD.txt")%>%column_to_rownames("FeatureID")
FAFTD.ald<-aldex.clr(round(pfamFAFT),condsFAvsFT, mc.samples=500, denom="all", verbose=F)

FAFTD.ttest<-aldex.ttest(FAFTD.ald,paired.test = FALSE, hist.plot=FALSE)%>%rownames_to_column("FeatureID")
FAFTD.effect<-aldex.effect(FAFTD.ald)%>%rownames_to_column("FeatureID")%>%left_join(.,FAFTD.ttest,by="FeatureID")
write.table(FAFTD.effect,"filtered_metaG/aldex/SFR23_0523_FAFTD_ald_effectwpval.txt",sep = "\t",row.names = FALSE, quote=FALSE)

FAFTD.effect<-fread("filtered_metaG/aldex/SFR23_0523_FAFTD_ald_effectwpval.txt")

FAFTD.effect.annot<-FAFTD.effect %>% left_join(.,annot, by ="FeatureID") %>%
  filter(!grepl("DUF",Name))%>%
  mutate(diffexpr=ifelse(wi.ep<0.05 & diff.btw>0, "up",
                          ifelse(wi.ep<0.05 & diff.btw< -0,"down","none"))) #%>%
  # filter(!str_detect(Name,"unknown"))

FAFTD.effect.annot %>% group_by(diffexpr) %>% tally(sort = TRUE)

# diffexpr     n
# <chr>    <int>
#   1 none      2857
# 2 down         3

##NA vs FT
md<-fread("~/scratch/metatranscript/metagenomic/metaG_metadata_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  arrange(condition,zt_time)
conds <- md$condition

pfamFTNA<-fread("filtered_metaG/pfam/pfam_FTNA.txt")%>%column_to_rownames("FeatureID")
condsFTvsNA<-conds[c(13:42)]
FTNA.ald<-aldex.clr(round(pfamFTNA),condsFTvsNA, mc.samples=500, denom="all", verbose=F)
FTNA.ttest<-aldex.ttest(FTNA.ald,paired.test = FALSE, hist.plot=FALSE)%>%rownames_to_column("FeatureID")
FTNA.effect<-aldex.effect(FTNA.ald)%>%rownames_to_column("FeatureID")%>%left_join(.,FTNA.ttest,by="FeatureID")%>%
  mutate(diff.btw=-(diff.btw))

write.table(FTNA.effect,"filtered_metaG/aldex/SFR23_0523_FTNA_ald_effectwpval.txt",sep = "\t",row.names = FALSE, quote=FALSE)

FTNA.effect<-fread("filtered_metaG/aldex/SFR23_0523_FTNA_ald_effectwpval.txt")

FTNA.effect.annot<-FTNA.effect %>% left_join(.,annot, by ="FeatureID") %>%
  filter(!grepl("DUF",Name))%>%
  mutate(diffexpr=ifelse(wi.ep<0.05 & diff.btw>0, "up",
                         ifelse(wi.ep<0.05 & diff.btw< -0,"down","none"))) #%>%
  #filter(!str_detect(Name,"unknown"))

FTNA.effect.annot %>% group_by(diffexpr) %>% tally(sort = TRUE)

# diffexpr     n
# <chr>    <int>
#   1 none      2847
# 2 up         427
# 3 down       153

##NA vs FT (light)
md<-fread("~/scratch/metatranscript/metagenomic/metaG_metadata_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition)) %>%
  filter(zt_time<13)%>%
  arrange(condition,zt_time)
conds <- md$condition

condsFTvsNA<-conds[c(7:21)]
pfamFTNA<-fread("filtered_metaG/pfam/pfam_FTNAL.txt")%>%column_to_rownames("FeatureID")
FTNAL.ald<-aldex.clr(round(pfamFTNA),condsFTvsNA, mc.samples=500, denom="all", verbose=F)

FTNAL.ttest<-aldex.ttest(FTNAL.ald,paired.test = FALSE, hist.plot=FALSE)%>%rownames_to_column("FeatureID")
FTNAL.effect<-aldex.effect(FTNAL.ald)%>%rownames_to_column("FeatureID")%>%left_join(.,FTNAL.ttest,by="FeatureID")%>%
  mutate(diff.btw=-(diff.btw))
write.table(FTNAL.effect,"filtered_metaG/aldex/SFR23_0524_FTNAL_ald_effectwpval.txt",sep = "\t",row.names = FALSE, quote=FALSE)

FTNAL.effect<-fread("filtered_metaG/aldex/SFR23_0524_FTNAL_ald_effectwpval.txt")

FTNAL.effect.annot<-FTNAL.effect %>% left_join(.,annot, by ="FeatureID") %>%
  filter(!grepl("DUF",Name))%>%
  mutate(diffexpr=ifelse(wi.ep<0.05 & diff.btw>0, "up",
                         ifelse(wi.ep<0.05 & diff.btw< -0,"down","none"))) #%>%
  #filter(!str_detect(Name,"unknown"))

FTNAL.effect.annot %>% group_by(diffexpr) %>% tally(sort = TRUE)

# diffexpr     n
# <chr>    <int>
#   1 none      3006
# 2 up         241
# 3 down        62

#NA vs FT (dark)
md<-fread("~/scratch/metatranscript/metagenomic/metaG_metadata_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition)) %>%
  filter(zt_time>9)%>%
  arrange(condition,zt_time)
conds <- md$condition

condsFTvsNA<-conds[c(7:21)]
pfamFTNA<-fread("filtered_metaG/pfam/pfam_FTNAD.txt")%>%column_to_rownames("FeatureID")
FTNAD.ald<-aldex.clr(round(pfamFTNA),condsFTvsNA, mc.samples=500, denom="all", verbose=F)

FTNAD.ttest<-aldex.ttest(FTNAD.ald,paired.test = FALSE, hist.plot=FALSE)%>%rownames_to_column("FeatureID")
FTNAD.effect<-aldex.effect(FTNAD.ald)%>%rownames_to_column("FeatureID")%>%left_join(.,FTNAD.ttest,by="FeatureID")%>%
  mutate(diff.btw=-(diff.btw))
write.table(FTNAD.effect,"filtered_metaG/aldex/SFR23_0524_FTNAD_ald_effectwpval.txt",sep = "\t",row.names = FALSE, quote=FALSE)

FTNAD.effect<-fread("filtered_metaG/aldex/SFR23_0524_FTNAD_ald_effectwpval.txt")

FTNAD.effect.annot<-FTNAD.effect %>% left_join(.,annot, by ="FeatureID") %>%
  filter(!grepl("DUF",Name))%>%
  mutate(diffexpr=ifelse(wi.ep<0.05 & diff.btw>0, "up",
                         ifelse(wi.ep<0.05 & diff.btw< -0,"down","none"))) #%>%
  #filter(!str_detect(Name,"unknown"))

FTNAD.effect.annot %>% group_by(diffexpr) %>% tally(sort = TRUE)

# diffexpr     n
# <chr>    <int>
#   1 none      2878
# 2 up         159
# 3 down        65

#NA vs FA
md<-fread("~/scratch/metatranscript/metagenomic/metaG_metadata_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  arrange(condition,zt_time)
conds <- md$condition

pfamFANA<-fread("filtered_metaG/pfam/pfam_FANA.txt")%>%column_to_rownames("FeatureID")
condsFAvsNA<-conds[c(1:12,25:42)]
FANA.ald<-aldex.clr(round(pfamFANA),condsFAvsNA, mc.samples=500, denom="all", verbose=F)
FANA.ttest<-aldex.ttest(FANA.ald,paired.test = FALSE, hist.plot=FALSE)%>%rownames_to_column("FeatureID")
FANA.effect<-aldex.effect(FANA.ald)%>%rownames_to_column("FeatureID")%>%left_join(.,FANA.ttest,by="FeatureID")%>%
  mutate(diff.btw=-(diff.btw))
write.table(FANA.effect,"filtered_metaG/aldex/SFR23_0524_FANA_ald_effectwpval.txt",sep = "\t",row.names = FALSE, quote=FALSE)

FANA.effect<-fread("filtered_metaG/aldex/SFR23_0524_FANA_ald_effectwpval.txt")

FANA.effect.annot<-FANA.effect %>% left_join(.,annot, by ="FeatureID") %>%
  filter(!grepl("DUF",Name))%>%
  mutate(diffexpr=ifelse(wi.ep<0.05 & diff.btw>0, "up",
                         ifelse(wi.ep<0.05 & diff.btw< -0,"down","none"))) #%>%
  #filter(!str_detect(Name,"unknown"))

FANA.effect.annot %>% group_by(diffexpr) %>% tally(sort = TRUE)

# diffexpr     n
# <chr>    <int>
#   1 none      2406
# 2 up         729
# 3 down       235

##NA vs FA (light)
md<-fread("~/scratch/metatranscript/metagenomic/metaG_metadata_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition)) %>%
  filter(zt_time<13)%>%
  arrange(condition,zt_time)
conds <- md$condition

condsFAvsNA<-conds[c(1:6,13:21)]
pfamFANA<-fread("filtered_metaG/pfam/pfam_FANAL.txt")%>%column_to_rownames("FeatureID")
FANAL.ald<-aldex.clr(round(pfamFANA),condsFAvsNA, mc.samples=500, denom="all", verbose=F)

FANAL.ttest<-aldex.ttest(FANAL.ald,paired.test = FALSE, hist.plot=FALSE)%>%rownames_to_column("FeatureID")
FANAL.effect<-aldex.effect(FANAL.ald)%>%rownames_to_column("FeatureID")%>%left_join(.,FANAL.ttest,by="FeatureID")%>%
  mutate(diff.btw=-(diff.btw))
write.table(FANAL.effect,"filtered_metaG/aldex/SFR23_0524_FANAL_ald_effectwpval.txt",sep = "\t",row.names = FALSE, quote=FALSE)

FANAL.effect<-fread("filtered_metaG/aldex/SFR23_0524_FANAL_ald_effectwpval.txt")

FANAL.effect.annot<-FANAL.effect %>% left_join(.,annot, by ="FeatureID") %>%
  filter(!grepl("DUF",Name))%>%
  mutate(diffexpr=ifelse(wi.ep<0.05 & diff.btw>0, "up",
                         ifelse(wi.ep<0.05 & diff.btw< -0,"down","none"))) #%>%
  #filter(!str_detect(Name,"unknown"))

FANAL.effect.annot %>% group_by(diffexpr) %>% tally(sort = TRUE)

# diffexpr     n
# <chr>    <int>
#   1 none      2718
# 2 up         441
# 3 down       107

#NA vs FA (dark)
md<-fread("~/scratch/metatranscript/metagenomic/metaG_metadata_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition)) %>%
  filter(zt_time>9)%>%
  arrange(condition,zt_time)
conds <- md$condition

# condsFAvsNA<-conds[c(1:6,13:21)]
pfamFANA<-fread("filtered_metaG/pfam/pfam_FANAD.txt")%>%column_to_rownames("FeatureID")
FANAD.ald<-aldex.clr(round(pfamFANA),condsFAvsNA, mc.samples=500, denom="all", verbose=F)

FANAD.ttest<-aldex.ttest(FANAD.ald,paired.test = FALSE, hist.plot=FALSE)%>%rownames_to_column("FeatureID")
FANAD.effect<-aldex.effect(FANAD.ald)%>%rownames_to_column("FeatureID")%>%left_join(.,FANAD.ttest,by="FeatureID")%>%
  mutate(diff.btw=-(diff.btw))
write.table(FANAD.effect,"filtered_metaG/aldex/SFR23_0524_FANAD_ald_effectwpval.txt",sep = "\t",row.names = FALSE, quote=FALSE)

FANAD.effect<-fread("filtered_metaG/aldex/SFR23_0524_FANAD_ald_effectwpval.txt")

FANAD.effect.annot<-FANAD.effect %>% left_join(.,annot, by ="FeatureID") %>%
  filter(!grepl("DUF",Name))%>%
  mutate(diffexpr=ifelse(wi.ep<0.05 & diff.btw>0, "up",
                         ifelse(wi.ep<0.05 & diff.btw< -0,"down","none"))) #%>%
  #filter(!str_detect(Name,"unknown"))

FANAD.effect.annot %>% group_by(diffexpr) %>% tally(sort = TRUE)

# diffexpr     n
# <chr>    <int>
#   1 none      2607
# 2 up         379
# 3 down        98
##########################################################

#FANA
235/3370
729/3370

#FTNA
153/3427
427/3427

#FAFT
23/3067
3/3067

aldex_summhits <- data.frame (comparison  = c("NAFA","NAFA","NAFT","NAFT","FAFT","FAFT"),
                              direction = c("down","up","down","up","down","up"),
                              perc= c(6.97,21.6,4.46,12.5,0.750,0.0978)) %>%
  mutate(comparison=factor(comparison, levels = c("NAFA","NAFT","FAFT")))

ggplot(aldex_summhits, aes(x=comparison, y=perc)) + 
  geom_bar(stat="identity", position="identity",fill="gray30") +
  theme_minimal() 

ggsave("filtered_metaG/aldex/SFR23_0524_summaldexhits_pless0.05.pdf", width = 2.2, height = 3)

#Summarise hits for all three pairwise comparisons based on light and dark (all-not just unique)

#FANAL
107/3266
441/3266

#FANAD
98/3084
379/3084

#FTNAL
62/3309
241/3309

#FTNAD
65/3102
159/3102

#FAFTL
34/2859
9/2859

#FAFTD
3/2860
0/2860

aldex_summhitsLD <- data.frame (comparison  = c("NAFA","NAFA","NAFT","NAFT","FAFT","FAFT","NAFA","NAFA","NAFT","NAFT","FAFT","FAFT"),
                              direction = c("down","up","down","up","down","up","down","up","down","up","down","up"),
                              phase= c("light","light","light","light","light","light","dark","dark","dark","dark","dark","dark"),
                              perc= c(3.28,13.5,1.87,7.28,1.19,0.315,-3.18,-12.3,-2.10,-5.13,-0.105,0)) %>%
  mutate(comparison=factor(comparison, levels = c("NAFA","NAFT","FAFT")),
         phase=factor(phase, levels=c("light","dark")))

ggplot(aldex_summhitsLD, aes(x=comparison, y=perc,fill=phase)) + 
  geom_bar(stat="identity", position="identity") +
  scale_fill_manual(values=c("gray70","gray10")) +
  theme_minimal() +
  theme(legend.position = "top")

ggsave("filtered_metaG/aldex/SFR23_0524_summaldexhitsLD_pless0.05.pdf", width = 2.2, height = 3)

###########################################################
#plot summary of the sig hits, just unique, no overlap

list_venn <- list(FAFT = (FAFT.effect.annot%>%filter(wi.ep<0.05))$FeatureID,
                  FANA = (FANA.effect.annot%>%filter(wi.ep<0.05))$FeatureID,
                  FTNA = (FTNA.effect.annot%>%filter(wi.ep<0.05))$FeatureID)

ggvenn(list_venn, c("FAFT", "FANA", "FTNA"),
       fill_color = c("white","black","gray"))
ggsave("filtered_metaG/aldex/SFR23_0524_venn_pless0.05.pdf", width = 3.5, height = 3)

424/3370
49/3427
10/3067

aldex_summhits_unique <- data.frame (comparison  = c("NAFA","NAFT","FAFT"),
                              perc= c(12.6,1.43,0.326))%>%
  mutate(comparison=factor(comparison, levels = c("NAFA","NAFT","FAFT")))

ggplot(aldex_summhits_unique, aes(x=comparison, y=perc)) + 
  geom_bar(stat="identity", position="identity",fill="gray30") +
  theme_minimal()

ggsave("filtered_metaG/aldex/SFR23_0524_summaldexhitsUQ_pless0.05.pdf", width = 2.2, height = 3)

#FANAL
278/3266

#FTNAL
41/3309

#FAFTL
33/2859

#FANAD
288/3084

#FTNAD
36/3102

#FAFTD
2/2860

aldex_summhitsLD <- data.frame (comparison  = c("NAFA","NAFT","FAFT","NAFA","NAFT","FAFT"),
                                phase= c("light","light","light","dark","dark","dark"),
                                perc= c(8.51,1.24,1.15,
                                        -9.34,-1.16,-0.0699))%>%
  mutate(comparison=factor(comparison, levels = c("NAFA","NAFT","FAFT")),
         phase=factor(phase, levels=c("light","dark")))

ggplot(aldex_summhitsLD, aes(x=comparison, y=perc,fill=phase)) + 
  geom_bar(stat="identity", position="identity") +
  scale_fill_manual(values=c("gray70","gray10")) +
  theme_minimal() +
  theme(legend.position = "top")

ggsave("filtered_metaG/aldex/SFR23_0524_summaldexhitsLD_unique_pless0.05.pdf", width = 2.2, height = 3)

###########################################################
#unique hits/venn diagram
list_venn <- list(FAFT = (FAFT.effect.annot%>%filter(wi.ep<0.05))$FeatureID,
                  FANA = (FANA.effect.annot%>%filter(wi.ep<0.05))$FeatureID,
                  FTNA = (FTNA.effect.annot%>%filter(wi.ep<0.05))$FeatureID)

ggvenn(list_venn, c("FAFT", "FANA", "FTNA"),
       fill_color = c("white","black","gray"))

ItemsList <- venn(list_venn, show.plot = FALSE)
all<-attributes(ItemsList)$intersections
#unique_hits_all<-FAFT.effect.annot%>%filter(FeatureID %in% all$FAFT)
#unique_hits_all<-FTNA.effect.annot%>%filter(FeatureID %in% all$FTNA)
unique_hits_all<-FANA.effect.annot%>%filter(FeatureID %in% all$FANA)
#write.table(unique_hits_all,"filtered_metaG/aldex/SFR23_0524_FAFT_10uniquehits_pless0.05.txt",sep = "\t",row.names = FALSE, quote=FALSE)
#write.table(unique_hits_all,"filtered_metaG/aldex/SFR23_0524_FTNA_49uniquehits_pless0.05.txt",sep = "\t",row.names = FALSE, quote=FALSE)
write.table(unique_hits_all,"filtered_metaG/aldex/SFR23_0524_FANA_424uniquehits_pless0.05.txt",sep = "\t",row.names = FALSE, quote=FALSE)

list_venn <- list(FAFT = (FAFTL.effect.annot%>%filter(wi.ep<0.05))$FeatureID,
                  FANA = (FANAL.effect.annot%>%filter(wi.ep<0.05))$FeatureID,
                  FTNA = (FTNAL.effect.annot%>%filter(wi.ep<0.05))$FeatureID)

ggvenn(list_venn, c("FAFT", "FANA", "FTNA"),
       fill_color = c("white","black","gray"))
ggsave("filtered_metaG/aldex/SFR23_0524_Light_venn_pless0.05.pdf", width = 3.5, height = 3)


ItemsList <- venn(list_venn, show.plot = FALSE)
light<-attributes(ItemsList)$intersections
#unique_hits_light<-FAFTL.effect.annot%>%filter(FeatureID %in% light$FAFT)
#unique_hits_lightFANA<-FANAL.effect.annot%>%filter(FeatureID %in% light$FANA)
unique_hits_lightFTNA<-FTNAL.effect.annot%>%filter(FeatureID %in% light$FTNA)
#write.table(unique_hits_light,"filtered_metaG/aldex/SFR23_0524_FAFTL_33uniquehits_pless0.05.txt",sep = "\t",row.names = FALSE, quote=FALSE)
#write.table(unique_hits_lightFANA,"filtered_metaG/aldex/SFR23_0524_FANAL_278uniquehits_pless0.05.txt",sep = "\t",row.names = FALSE, quote=FALSE)
write.table(unique_hits_lightFTNA,"filtered_metaG/aldex/SFR23_0524_FTNAL_41uniquehits_pless0.05.txt",sep = "\t",row.names = FALSE, quote=FALSE)

list_venn <- list(FAFT = (FAFTD.effect.annot%>%filter(wi.ep<0.05))$FeatureID,
                  FANA = (FANAD.effect.annot%>%filter(wi.ep<0.05))$FeatureID,
                  FTNA = (FTNAD.effect.annot%>%filter(wi.ep<0.05))$FeatureID)

ggvenn(list_venn, c("FAFT", "FANA", "FTNA"),
       fill_color = c("white","black","gray"))
ggsave("filtered_metaG/aldex/SFR23_0524_Dark_venn_pless0.05.pdf", width = 3.5, height = 3)

ItemsList <- venn(list_venn, show.plot = FALSE)
dark<-attributes(ItemsList)$intersections
#unique_hits_dark<-FAFTD.effect.annot%>%filter(FeatureID %in% dark$FAFT)
#unique_hits_darkFANA<-FANAD.effect.annot%>%filter(FeatureID %in% dark$FANA)
unique_hits_darkFTNA<-FTNAD.effect.annot%>%filter(FeatureID %in% dark$FTNA)
#write.table(unique_hits_dark,"filtered_metaG/aldex/SFR22_0524_FAFTD_2uniquehits_pless0.05.txt",sep = "\t",row.names = FALSE, quote=FALSE)
#write.table(unique_hits_darkFANA,"filtered_metaG/aldex/SFR23_0524_FANAD_288uniquehits_pless0.05.txt",sep = "\t",row.names = FALSE, quote=FALSE)
write.table(unique_hits_darkFTNA,"filtered_metaG/aldex/SFR23_0524_FTNAD_36uniquehits_pless0.05.txt",sep = "\t",row.names = FALSE, quote=FALSE)

list_venn <- list(light = unique_hits_light$FeatureID,
                  dark = unique_hits_dark$FeatureID)

p<-venn.diagram(list_venn, fill = c("white", "black"),height = 10,
                 width = 10,alpha = c(0.5, 0.5), lwd =0,filename = NULL)

grid.draw(p)
pdf(file="filtered_metaG/aldex/SFR23_0524_FAFT_LDintersect_pless0.05.pdf")
grid.draw(p)
dev.off()

list_venn <- list(light = unique_hits_lightFTNA$FeatureID,
                  dark = unique_hits_darkFTNA$FeatureID)

p<-venn.diagram(list_venn, fill = c("white", "black"),height = 10,
                width = 10,alpha = c(0.5, 0.5), lwd =0,filename = NULL)

grid.draw(p)
pdf(file="filtered_metaG/aldex/SFR23_0524_FTNA_LDintersect_pless0.05.pdf")
grid.draw(p)
dev.off()

list_venn <- list(light = unique_hits_lightFANA$FeatureID,
                  dark = unique_hits_darkFANA$FeatureID)

p<-venn.diagram(list_venn, fill = c("white", "black"),height = 10,
                width = 10,alpha = c(0.5, 0.5), lwd =0,filename = NULL)

grid.draw(p)
pdf(file="filtered_metaG/aldex/SFR23_0524_FANA_LDintersect_pless0.05.pdf")
grid.draw(p)
dev.off()
###########################################################

#volcano
gonames<-fread("/mnt/zarrinpar/scratch/sfloresr/metatranscript/metatranscript/woltka2_results/go_name.txt")
pfamFAFTD<-read.table("/mnt/zarrinpar/scratch/sfloresr/metatranscript/metatranscript/woltka2_results/pfam-to-go-process.map",header = FALSE, sep = "\t",
                    col.names = paste0("V",seq_len(4)), fill = TRUE)%>%
  gather(column,GO_Term,-V1)%>%
  dplyr::select(1,3)%>%
  dplyr::rename(FeatureID=V1)%>%
  filter(FeatureID %in% unique_hits_dark$FeatureID)%>%
  left_join(.,gonames,by="GO_Term")%>%
  filter(!is.na(name))%>%
  mutate(phase="dark")

pfamFAFTL<-read.table("/mnt/zarrinpar/scratch/sfloresr/metatranscript/metatranscript/woltka2_results/pfam-to-go-process.map",header = FALSE, sep = "\t",
                      col.names = paste0("V",seq_len(4)), fill = TRUE)%>%
  gather(column,GO_Term,-V1)%>%
  dplyr::select(1,3)%>%
  dplyr::rename(FeatureID=V1)%>%
  filter(FeatureID %in% unique_hits_light$FeatureID)%>%
  left_join(.,gonames,by="GO_Term")%>%
  filter(!is.na(name))%>%
  mutate(phase="light")

pfamFAFT<-rbind(pfamFAFTD,pfamFAFTL)


#pfamFAFT_summ<-pfamFAFTD%>%
#pfamFAFT_summ<-pfamFAFTL%>%
pfamFAFT_summ<-pfamFAFT%>%
  group_by(name,phase)%>%summarise(n=n())%>%
  group_by(name)%>%mutate(sum=sum(n))%>%
  mutate(n=ifelse(phase=="light", n*-1,n))%>%
  arrange(sum)

pfamFAFT_summ$name <- factor(pfamFAFT_summ$name,levels = unique(pfamFAFT_summ$name))

ggplot(data=pfamFAFT_summ, aes(x=name, y=n, fill=phase)) +
  geom_bar(stat="identity", position="identity") +
  scale_fill_manual(values=c("gray70", "gray10")) +coord_flip() +theme_pubr() +
  #ggtitle("FAFTD")+
  #ggtitle("FAFTL")+
  ggtitle("FAFT")+
  scale_y_continuous(expand=c(0,0))

ggsave("filtered_metaG/aldex/SFR23_0525_FAFT_GOterms.pdf",height=6.2, width=6.5)

#volcano plot for FAvFT

# pfamFAFTL_top5<-pfamFAFTL%>%filter(name %in% c("translation","carbohydrate metabolic process",
#                                                "transmembrane transport","regulation of transcription, DNA-templated",
#                                                "proteolysis"))
# FAFTL.effect.annot<-FAFTL.effect%>%left_join(.,pfamFAFTL_top5,by="FeatureID")%>%
#   mutate(name=ifelse(is.na(name) & (FeatureID %in% unique_hits_light$FeatureID), "sig",name))%>%
#   mutate(name=factor(name, levels=c("translation","carbohydrate metabolic process",
#                                     "transmembrane transport","regulation of transcription, DNA-templated",
#                                     "proteolysis", "sig")))
# 
# p<-ggplot(data=FAFTL.effect.annot, aes(x=diff.btw, y=-log10(wi.ep),color=name)) +
#   ggtitle("FA vs. FT (Light)")+
#   geom_point(size = 0.1) +
#   xlim(-8,8) +
#   theme_pubr() +
#   theme(legend.position = "bottom") +
#   geom_hline(yintercept=-log10(0.05), linetype = "dashed",col="grey10")+
#   scale_color_manual(values=c("orange","black"))
# 
# ggsave("aldex/SFR23_0409_FAFTL_noDUF_volcanoplot_wgoannot_pless0.05.pdf", width = 3.5, height = 3)
# 
# pfamFAFTD_top5<-pfamFAFTD%>%filter(name %in% c("translation","carbohydrate metabolic process",
#                                                "transmembrane transport","regulation of transcription, DNA-templated",
#                                                "proteolysis"))
# FAFTD.effect.annot<-FAFTD.effect%>%left_join(.,pfamFAFTD_top5,by="FeatureID")%>%
#   mutate(name=ifelse(is.na(name) & (FeatureID %in% unique_hits_dark$FeatureID), "sig",name))%>%
#   mutate(name=factor(name, levels=c("translation","carbohydrate metabolic process",
#                                     "transmembrane transport","regulation of transcription, DNA-templated",
#                                     "proteolysis", "sig")))
# 
# p<-ggplot(data=FAFTD.effect.annot, aes(x=diff.btw, y=-log10(wi.ep),color=name)) +
#   ggtitle("FA vs. FT (Dark)")+
#   geom_point(size = 0.1) +  #ylim(0,8) + 
#   xlim(-8,8) +
#   theme_pubr() +
#   theme(legend.position = "bottom") +
#   geom_hline(yintercept=-log10(0.05), linetype = "dashed",col="grey10")+
#   scale_color_manual(values=c("red","chartreuse3","orange","purple","black"))
# 
# ggsave("aldex/SFR23_0409_FAFTD_noDUF_volcanoplot_wgoannot_pless0.05.pdf", width = 3.5, height = 3)

pfamFANAD<-read.table("/mnt/zarrinpar/scratch/sfloresr/metatranscript/metatranscript/woltka2_results/pfam-to-go-process.map",header = FALSE, sep = "\t",
                      col.names = paste0("V",seq_len(4)), fill = TRUE)%>%
  gather(column,GO_Term,-V1)%>%
  dplyr::select(1,3)%>%
  dplyr::rename(FeatureID=V1)%>%
  filter(FeatureID %in% unique_hits_darkFANA$FeatureID)%>%
  left_join(.,gonames,by="GO_Term")%>%
  filter(!is.na(name))%>%
  mutate(phase="dark")

pfamFANAL<-read.table("/mnt/zarrinpar/scratch/sfloresr/metatranscript/metatranscript/woltka2_results/pfam-to-go-process.map",header = FALSE, sep = "\t",
                      col.names = paste0("V",seq_len(4)), fill = TRUE)%>%
  gather(column,GO_Term,-V1)%>%
  dplyr::select(1,3)%>%
  dplyr::rename(FeatureID=V1)%>%
  filter(FeatureID %in% unique_hits_lightFANA$FeatureID)%>%
  left_join(.,gonames,by="GO_Term")%>%
  filter(!is.na(name))%>%
  mutate(phase="light")

pfamFANA<-rbind(pfamFANAD,pfamFANAL)

#pfamFANA_summ<-pfamFANAD%>%
#pfamFANA_summ<-pfamFANAL%>%
pfamFANA_summ<-pfamFANA%>%
  group_by(name,phase)%>%summarise(n=n())%>%
  #group_by(name)%>%summarise(n=n())%>%
  group_by(name)%>%mutate(sum=sum(n))%>%
  mutate(n=ifelse(phase=="light", n*-1,n),
         phase=factor(phase, levels=c("light","dark")))%>%
  arrange(sum)%>%
  filter(sum>1)

pfamFANA_summ$name <- factor(pfamFANA_summ$name,levels = unique(pfamFANA_summ$name))
#pfamFANA_summ$name <- factor(pfamFANA_summ$name,levels = pfamFANA_summ$name)

ggplot(data=pfamFANA_summ, aes(x=name, y=n, fill=phase)) +
  geom_bar(stat="identity", position="identity") +
  scale_fill_manual(values=c("gray70", "gray10")) +
  coord_flip() +theme_pubr() +
  #ggtitle("FANAD")+
  #ggtitle("FANAL")+
  ggtitle("FANA")+
  scale_y_continuous(expand=c(0,0))

ggsave("filtered_metaG/aldex/SFR23_0524_FANA_GOterms.pdf",height=6.5, width=6.5)

#make a volcano plot for the light and dark FANA hits showing the GO process as colors

# pfamFANAL_top5<-pfamFANAL%>%filter(name %in% c("translation","carbohydrate metabolic process",
#                                                "transmembrane transport","regulation of transcription, DNA-templated",
#                                                "proteolysis"))
# FANAL.effect.annot<-FANAL.effect%>%left_join(.,pfamFANAL_top5,by="FeatureID")%>%
#   mutate(name=ifelse(is.na(name) & (FeatureID %in% unique_hits_lightFANA$FeatureID), "sig",name))%>%
#   mutate(name=factor(name, levels=c("translation","carbohydrate metabolic process",
#                                     "transmembrane transport","regulation of transcription, DNA-templated",
#                                     "proteolysis", "sig")))
# 
# p<-ggplot(data=FANAL.effect.annot, aes(x=diff.btw, y=-log10(wi.ep),color=name)) +
#   ggtitle("NA vs. FA (Light)")+
#   geom_point(size = 0.1) +  #ylim(0,8) + 
#   xlim(-8,8) +
#   theme_pubr() +
#   theme(legend.position = "bottom") +
#   geom_hline(yintercept=-log10(0.05), linetype = "dashed",col="grey10")+
#   scale_color_manual(values=c("blue", "red","chartreuse3","orange","purple","black"))
# 
# ggsave("aldex/SFR23_0409_FANAL_noDUF_volcanoplot_wgoannot_pless0.05.pdf", width = 3.5, height = 3)
# 
# pfamFANAD_top5<-pfamFANAD%>%filter(name %in% c("translation","carbohydrate metabolic process",
#                                                "transmembrane transport","regulation of transcription, DNA-templated",
#                                                "proteolysis"))
# FANAD.effect.annot<-FANAD.effect%>%left_join(.,pfamFANAD_top5,by="FeatureID")%>%
#   mutate(name=ifelse(is.na(name) & (FeatureID %in% unique_hits_darkFANA$FeatureID), "sig",name))%>%
#   mutate(name=factor(name, levels=c("translation","carbohydrate metabolic process",
#                                     "transmembrane transport","regulation of transcription, DNA-templated",
#                                     "proteolysis", "sig")))
# 
# p<-ggplot(data=FANAD.effect.annot, aes(x=diff.btw, y=-log10(wi.ep),color=name)) +
#   ggtitle("NA vs. FA (Dark)")+
#   geom_point(size = 0.1) +  #ylim(0,8) + 
#   xlim(-8,8) +
#   theme_pubr() +
#   theme(legend.position = "bottom") +
#   geom_hline(yintercept=-log10(0.05), linetype = "dashed",col="grey10")+
#   scale_color_manual(values=c("blue", "red","chartreuse3","orange","purple","black"))
# 
# ggsave("aldex/SFR23_0409_FANAD_noDUF_volcanoplot_wgoannot_pless0.05.pdf", width = 3.5, height = 3)

pfamFTNAD<-read.table("/mnt/zarrinpar/scratch/sfloresr/metatranscript/metatranscript/woltka2_results/pfam-to-go-process.map",header = FALSE, sep = "\t",
                      col.names = paste0("V",seq_len(4)), fill = TRUE)%>%
  gather(column,GO_Term,-V1)%>%
  dplyr::select(1,3)%>%
  dplyr::rename(FeatureID=V1)%>%
  filter(FeatureID %in% unique_hits_darkFTNA$FeatureID)%>%
  left_join(.,gonames,by="GO_Term")%>%
  filter(!is.na(name))%>%
  mutate(phase="dark")

pfamFTNAL<-read.table("/mnt/zarrinpar/scratch/sfloresr/metatranscript/metatranscript/woltka2_results/pfam-to-go-process.map",header = FALSE, sep = "\t",
                      col.names = paste0("V",seq_len(4)), fill = TRUE)%>%
  gather(column,GO_Term,-V1)%>%
  dplyr::select(1,3)%>%
  dplyr::rename(FeatureID=V1)%>%
  filter(FeatureID %in% unique_hits_lightFTNA$FeatureID)%>%
  left_join(.,gonames,by="GO_Term")%>%
  filter(!is.na(name))%>%
  mutate(phase="light")


pfamFTNA<-rbind(pfamFTNAD,pfamFTNAL)

#pfamFTNA_summ<-pfamFTNAD%>%
#pfamFTNA_summ<-pfamFTNAL%>%
pfamFTNA_summ<-pfamFTNA%>%
  group_by(name,phase)%>%summarise(n=n())%>%
  #group_by(name)%>%summarise(n=n())%>%
  group_by(name)%>%mutate(sum=sum(n))%>%
  mutate(n=ifelse(phase=="light", n*-1,n),
         phase=factor(phase, levels=c("light","dark")))%>%
  arrange(sum)#%>%
  #arrange(n)
  #filter(sum>1)

pfamFTNA_summ$name <- factor(pfamFTNA_summ$name,levels = unique(pfamFTNA_summ$name))
#pfamFTNA_summ$name <- factor(pfamFTNA_summ$name,levels = pfamFTNA_summ$name)

ggplot(data=pfamFTNA_summ, aes(x=name, y=n, fill=phase)) +
  geom_bar(stat="identity", position="identity") +
  scale_fill_manual(values=c("gray70", "gray10")) +
  coord_flip() +theme_pubr() +
  #ggtitle("FTNAD")+
  #ggtitle("FTNAL")+
  ggtitle("FTNA")+
  scale_y_continuous(expand=c(0,0))

ggsave("filtered_metaG/aldex/SFR23_0524_FTNA_GOterms.pdf",height=6, width=8)

#make a volcano plot for the light and dark FTNA hits showing the GO process as colors

pfamFTNAL_top5<-pfamFTNAL%>%filter(name %in% c("translation","carbohydrate metabolic process",
                                               "transmembrane transport","regulation of transcription, DNA-templated",
                                               "proteolysis"))
FTNAL.effect.annot<-FTNAL.effect%>%left_join(.,pfamFTNAL_top5,by="FeatureID")%>%
  mutate(name=ifelse(is.na(name) & (FeatureID %in% unique_hits_lightFTNA$FeatureID), "sig",name))%>%
  mutate(name=factor(name, levels=c("translation","carbohydrate metabolic process",
                                    "transmembrane transport","regulation of transcription, DNA-templated",
                                    "proteolysis", "sig")))

p<-ggplot(data=FTNAL.effect.annot, aes(x=diff.btw, y=-log10(wi.ep),color=name)) +
  ggtitle("NA vs. FT (Light)")+
  geom_point(size = 0.1) +  #ylim(0,8) + 
  xlim(-8,8) +
  theme_pubr() +
  theme(legend.position = "bottom") +
  geom_hline(yintercept=-log10(0.05), linetype = "dashed",col="grey10")+
  scale_color_manual(values=c("blue", "red","chartreuse3","orange","purple","black"))

ggsave("aldex/SFR23_0409_FTNAL_noDUF_volcanoplot_wgoannot_pless0.05.pdf", width = 3.5, height = 3)

pfamFTNAD_top5<-pfamFTNAD%>%filter(name %in% c("translation","carbohydrate metabolic process",
                                               "transmembrane transport","regulation of transcription, DNA-templated",
                                               "proteolysis"))

FTNAD.effect.annot<-FTNAD.effect%>%left_join(.,pfamFTNAD_top5,by="FeatureID")%>%
  mutate(name=ifelse(is.na(name) & (FeatureID %in% unique_hits_darkFTNA$FeatureID), "sig",name))%>%
  mutate(name=factor(name, levels=c("translation","carbohydrate metabolic process",
                                    "transmembrane transport","regulation of transcription, DNA-templated",
                                    "proteolysis", "sig")))

p<-ggplot(data=FTNAD.effect.annot, aes(x=diff.btw, y=-log10(wi.ep),color=name)) +
  ggtitle("NA vs. FT (Dark)")+
  geom_point(size = 0.1) +  #ylim(0,8) + 
  xlim(-8,8) +
  theme_pubr() +
  theme(legend.position = "bottom") +
  geom_hline(yintercept=-log10(0.05), linetype = "dashed",col="grey10")+
  scale_color_manual(values=c("blue", "red","chartreuse3","orange","purple","black"))

ggsave("aldex/SFR23_0409_FTNAD_noDUF_volcanoplot_wgoannot_pless0.05.pdf", width = 3.5, height = 3)

###########################################################

unique_hits_light<-fread("aldex/SFR23_0409_FAFTL_25uniquehits_pless0.05.txt")
unique_hits_dark<-fread("aldex/SFR22_0409_FAFTD_122uniquehits_pless0.05.txt")

list_venn <- list(light = unique_hits_light$FeatureID,
                  dark = unique_hits_dark$FeatureID)

# unique_hits_lightFANA<-fread("aldex/SFR23_0111_FANAL_188uniquehits_pless0.05.txt")
# unique_hits_lightFTNA<-fread("aldex/SFR23_0409_FTNAL_461uniquehits_pless0.05.txt")
# unique_hits_darkFANA<-fread("aldex/SFR23_0409_FANAD_293uniquehits_pless0.05.txt")
# unique_hits_darkFTNA<-fread("aldex/SFR23_0409_FTNAD_259uniquehits_pless0.05.txt")

ItemsList <- venn(list_venn, show.plot = FALSE)
FAFTonly<-attributes(ItemsList)$intersections
lightuniquehits<-FAFTL.effect.annot%>%filter((FeatureID %in% FAFTonly$light))
write.table(lightuniquehits,"aldex/SFR23_0410_FAFTL_20uniquehitsFAFTonly_pless0.05.txt",sep = "\t",row.names = FALSE, quote=FALSE)
darkuniquehits<-FAFTD.effect.annot%>%filter((FeatureID %in% FAFTonly$dark))
write.table(darkuniquehits,"aldex/SFR23_0410_FAFTD_117uniquehitsFAFTonly_pless0.05.txt",sep = "\t",row.names = FALSE, quote=FALSE)
bothuniquehits<-FAFTD.effect.annot%>%filter((FeatureID %in% FAFTonly$`light:dark`))
write.table(bothuniquehits,"aldex/SFR23_0410_FAFTLD_5uniquehitsFAFTonly_pless0.05.txt",sep = "\t",row.names = FALSE, quote=FALSE)

#pfam to go process
lightuniquehits<-fread("aldex/SFR23_0410_FAFTL_20uniquehitsFAFTonly_pless0.05.txt")%>%
  #mutate(FeatureID=gsub("\\..*","",FeatureID),sig_level="light_sig")%>%
  mutate(sig_level="light_sig")%>%
  dplyr::select(FeatureID,sig_level)

darkuniquehits<-fread("aldex/SFR23_0410_FAFTD_117uniquehitsFAFTonly_pless0.05.txt")%>%
  #mutate(FeatureID=gsub("\\..*","",FeatureID),sig_level="dark_sig")%>%
  mutate(sig_level="dark_sig")%>%
  dplyr::select(FeatureID,sig_level)

bothuniquehits<-fread("aldex/SFR23_0410_FAFTLD_5uniquehitsFAFTonly_pless0.05.txt")%>%
  #mutate(FeatureID=gsub("\\..*","",FeatureID),sig_level="both_sig")%>%
  mutate(sig_level="both_sig")%>%
  dplyr::select(FeatureID,sig_level)

all_uniquehits<-rbind(darkuniquehits, lightuniquehits, bothuniquehits)

gonames<-fread("go_name.txt")
#max(count.fields("pfam-to-go-process.map", sep = '\t'))

#both
pfamGOp<-read.table("pfam-to-go-process.map",header = FALSE, sep = "\t", 
                    col.names = paste0("V",seq_len(4)), fill = TRUE)%>%
  gather(column,GO_Term,-V1)%>%
  dplyr::select(1,3)%>%
  dplyr::rename(FeatureID=V1)%>%
  filter(FeatureID %in% all_uniquehits$FeatureID)%>%
  left_join(.,gonames,by="GO_Term")%>%
  left_join(.,all_uniquehits, by="FeatureID")%>%
  filter(!is.na(name))

# pitabB<-pfamGOp%>%group_by(name,sig_level)%>%summarise(n=n())%>%
#   mutate(perc = n / sum(n)) %>% 
#   arrange(perc) %>%
#   mutate(labels = scales::percent(perc))

# p<-ggplot(data=pitabB, aes(x=reorder(name, n), y=n,fill=sig_level)) +
#   geom_bar(stat="identity") + coord_flip()

#light
pfamGOp<-read.table("pfam-to-go-process-2.map", header=TRUE)%>%
  dplyr::rename(FeatureID=pfam, GO_Term=go_term)%>%
  filter(FeatureID %in% lightuniquehits$FeatureID)%>%
    #filter(FeatureID %in% darkuniquehits$FeatureID)
  left_join(.,gonames,by="GO_Term")

# pitabL<-pfamGOp%>%group_by(name)%>%summarise(n=n())%>%
#   mutate(perc = n / sum(n)) %>% 
#   arrange(perc) %>%
#   mutate(labels = scales::percent(perc))

#dark
pfamGOp2<-read.table("pfam-to-go-process-2.map", header=TRUE)%>%
  dplyr::rename(FeatureID=pfam, GO_Term=go_term)%>%
  filter(FeatureID %in% darkuniquehits$FeatureID)%>%
  left_join(.,gonames,by="GO_Term")
  
# pitabD<-pfamGOp2%>%group_by(name)%>%summarise(n=n())%>%
#   mutate(perc = n / sum(n)) %>% 
#   arrange(perc) %>%
#   mutate(labels = scales::percent(perc))

#plot specific functional examples
mdT<-fread("metaT_metadata_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))
#PF13561.9 ACP, fatty acid elongation, PF02036.20 lipid transport
metaT_dat<-fread("pfam/pfam_clean_noNT_TPM.txt") %>%
  dplyr::filter(grepl("PF02253.18",FeatureID))%>%
  #dplyr::filter(grepl("PF01219.22",FeatureID))%>%
  #dplyr::filter(grepl("PF13561.9",FeatureID))%>%
  #dplyr::filter(grepl("PF02036.20",FeatureID))%>%
  gather(sample_name,TPM_counts,-FeatureID)%>%
  left_join(.,mdT,by="sample_name")%>%
  mutate(condition=factor(condition, levels = c("NA","FA","FT")),
         lightdark=factor(lightdark, levels = c("light","dark")))

p<-ggplot(metaT_dat,aes(x=condition, y=log10(TPM_counts+1), fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  scale_fill_manual(values=c("#009E73","#D55E00","#0072B2"))+
  facet_grid(~lightdark,scales = "free_y")+
  theme_bw()+
  labs(x="condition",y="log10(TPM)",title="PF02253.18: Phospholipase A1")+
  #labs(x="condition",y="log10(TPM)",title="PF01219.22: Prokaryotic diacylglycerol kinase")+
  #labs(x="condition",y="log10(TPM)",title="PF13561.9: Enoyl-(Acyl carrier protein) reductase")+
  #labs(x="condition",y="log10(TPM)",title="PF02036.20: SCP-2 sterol transfer family")+
  theme(legend.position = "none")

#ggsave("aldex/specific_hits/SFR23_0410_PF02036.20_boxplotLD.pdf", width = 3.5, height = 3)
#ggsave("aldex/specific_hits/SFR23_0410_PF13561.9_boxplotLD.pdf", width = 3.5, height = 3)
#ggsave("aldex/specific_hits/SFR23_0410_PF01219.22_boxplotLD.pdf", width = 3.5, height = 3)
ggsave("aldex/specific_hits/SFR23_0410_PF02253.18_boxplotLD.pdf", width = 3.5, height = 3)

# Basic piechart

colourCount = length(pitabL$name)
getPalette = colorRampPalette(brewer.pal(6, "Dark2"))

ggplot(pitabL, aes(x="", y=n, fill=fct_reorder(name, n)))+
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + guides(fill=guide_legend(title="GO process")) +
  geom_text(aes(x=1.7, label = labels),
            position = position_stack(vjust = 0.5))+
  theme(legend.position="right",
        legend.title = element_text(size = 5, face="bold"), 
        legend.text = element_text(size = 5),
        legend.key.size = unit(0.5,"line")) +
  scale_fill_manual(values = getPalette(colourCount))

#ggsave("aldex/SFR22_0906_FAFTD_118uniquehitsFAFTonly_GOp.pdf", width = 6, height = 3)
ggsave("aldex/SFR22_0906_FAFTL_20uniquehitsFAFTonly_GOp.pdf", width = 3.5, height = 3)

#heatmap of all sig results
# md<-fread("metaT_metadata_noNT.txt")%>%
#   mutate(condition=ifelse(is.na(condition),"NA",condition))
# annot<-fread("pfam/pfam_annotationkey.csv")
# 
# lightuniquehits<-fread("aldex/SFR22_0829_FAFTL_20uniquehitsFAFTonly_pless0.05.txt")%>%
#   mutate(sig_level="light_sig")%>%dplyr::select(FeatureID,sig_level)
# darkuniquehits<-fread("aldex/SFR22_0829_FAFTD_118uniquehitsFAFTonly_pless0.05.txt")%>%
#   mutate(sig_level="dark_sig")%>%dplyr::select(FeatureID,sig_level)
# bothuniquehits<-fread("aldex/SFR22_1026_FAFTLD_4uniquehitsFAFTonly_pless0.05.txt")%>%
#   mutate(sig_level="both_sig")%>%dplyr::select(FeatureID,sig_level)
# 
# all_uniquehits<-rbind(darkuniquehits, lightuniquehits, bothuniquehits)
# 
# pfamTPM_all<-fread("pfam/pfam_clean_noNT_TPM.txt")%>%
#   filter(FeatureID %in% all_uniquehits$FeatureID)%>% 
#   gather(sample_name,TPM_counts,-FeatureID) %>%
#   left_join(.,annot, by ="FeatureID") %>%
#   mutate(label_name=paste(FeatureID, Name, sep=" "))%>%
#   left_join(.,md,by="sample_name") %>%
#   left_join(.,all_uniquehits,by="FeatureID") %>%
#   group_by(label_name,condition,lightdark,sig_level)%>%summarise(mn_TPM=mean(TPM_counts))%>%
#   group_by(label_name)%>%mutate(Zscore=scale(mn_TPM))%>%
#   mutate(condition=factor(condition, levels = c("FT", "FA","NA")),
#          lightdark=factor(lightdark, levels = c("light","dark")),
#          sig_level=factor(sig_level, levels = c("light_sig","dark_sig","both_sig")))
# 
# sigFAFTLD_mtx<-pfamTPM_all%>%filter(lightdark=="dark")%>%
#   dplyr::select(label_name,condition,Zscore)%>%spread(condition,Zscore)
# sigFAFTLD_mtx_mtx <- sigFAFTLD_mtx%>%column_to_rownames("label_name")%>%as.matrix()
# sigFAFTLD_mtx.dendro <- as.dendrogram(hclust(d = dist(x = sigFAFTLD_mtx_mtx)))
# sigFAFTLD_mtx_order <- rownames(sigFAFTLD_mtx_mtx)[order.dendrogram(sigFAFTLD_mtx.dendro)]
# pfamTPM_all$label_name <- factor(pfamTPM_all$label_name,levels = sigFAFTLD_mtx_order )
# 
# plt<-ggplot(pfamTPM_all,aes(x=condition, y=label_name)) +theme_classic()+
#   geom_tile(aes(fill=Zscore))+
#   scale_x_discrete(expand = c(0, 0))+ facet_grid(sig_level~lightdark,scales="free",space="free")+
#   theme(axis.ticks.y=element_blank(),panel.spacing.x=unit(0, "lines"),axis.text.y = element_text(size = 5))+
#   scale_fill_viridis(option="inferno") + 
#   xlab("condition") + ylab("transcripts")+labs(fill='TPM Z-score')+ggtitle("")
# ggsave("aldex/SFR22_1026_heatmap_lightHits_Zscore_pless0.05.pdf", plot=plt,height=6, width=6)


#heatmap of light results
# md<-fread("metaT_metadata_noNT.txt")%>%
#   mutate(condition=ifelse(is.na(condition),"NA",condition))
# annot<-fread("pfam/pfam_annotationkey.csv")
# 
# pfamTPM_all<-fread("pfam/pfam_clean_noNT_TPM.txt")%>%
#   filter(FeatureID %in% lightuniquehits$FeatureID)%>% 
#   gather(sample_name,TPM_counts,-FeatureID) %>%
#   left_join(.,annot, by ="FeatureID") %>%
#   mutate(label_name=paste(FeatureID, Name, sep=" "))%>%
#   left_join(.,md,by="sample_name") %>%
#   group_by(label_name,condition,lightdark)%>%summarise(mn_TPM=mean(TPM_counts))%>%
#   group_by(label_name)%>%mutate(Zscore=scale(mn_TPM))%>%
#   mutate(condition=factor(condition, levels = c("FT", "FA","NA")),
#          lightdark=factor(lightdark, levels = c("light","dark")))
# 
# 
# sigFAFTLD_mtx<-pfamTPM_all%>%filter(lightdark=="light")%>%
#   dplyr::select(label_name,condition,Zscore)%>%spread(condition,Zscore)
# sigFAFTLD_mtx_mtx <- sigFAFTLD_mtx%>%column_to_rownames("label_name")%>%as.matrix()
# sigFAFTLD_mtx.dendro <- as.dendrogram(hclust(d = dist(x = sigFAFTLD_mtx_mtx)))
# sigFAFTLD_mtx_order <- rownames(sigFAFTLD_mtx_mtx)[order.dendrogram(sigFAFTLD_mtx.dendro)]
# pfamTPM_all$label_name <- factor(pfamTPM_all$label_name,levels = sigFAFTLD_mtx_order )
# 
# plt<-ggplot(pfamTPM_all,aes(x=condition, y=label_name)) +theme_classic()+
#   geom_tile(aes(fill=Zscore))+
#   scale_x_discrete(expand = c(0, 0))+ facet_wrap(~lightdark,scales="free_x")+
#   theme(axis.ticks.y=element_blank(),panel.spacing.x=unit(0, "lines"),axis.text.y = element_text(size = 5))+
#   scale_fill_viridis(option="inferno") + 
#   xlab("condition") + ylab("transcripts")+labs(fill='TPM Z-score')+ggtitle("")
# ggsave("aldex/SFR22_0829_heatmap_lightHits_Zscore_pless0.05.pdf", plot=plt,height=4, width=6)

#heatmap of dark results
# md<-fread("metaT_metadata_noNT.txt")%>%
#   mutate(condition=ifelse(is.na(condition),"NA",condition))
# annot<-fread("pfam/pfam_annotationkey.csv")
# 
# pfamTPM_all<-fread("pfam/pfam_clean_noNT_TPM.txt")%>%
#   filter(FeatureID %in% darkuniquehits$FeatureID)%>% 
#   gather(sample_name,TPM_counts,-FeatureID) %>%
#   left_join(.,annot, by ="FeatureID") %>%
#   mutate(label_name=paste(FeatureID, Name, sep=" "))%>%
#   left_join(.,md,by="sample_name") %>%
#   group_by(label_name,condition,lightdark)%>%summarise(mn_TPM=mean(TPM_counts))%>%
#   group_by(label_name)%>%mutate(Zscore=scale(mn_TPM))%>%
#   mutate(condition=factor(condition, levels = c("FT", "FA","NA")),
#          lightdark=factor(lightdark, levels = c("light","dark")))
# 
# 
# sigFAFTLD_mtx<-pfamTPM_all%>%filter(lightdark=="dark")%>%
#   dplyr::select(label_name,condition,Zscore)%>%spread(condition,Zscore)
# sigFAFTLD_mtx_mtx <- sigFAFTLD_mtx%>%column_to_rownames("label_name")%>%as.matrix()
# sigFAFTLD_mtx.dendro <- as.dendrogram(hclust(d = dist(x = sigFAFTLD_mtx_mtx)))
# sigFAFTLD_mtx_order <- rownames(sigFAFTLD_mtx_mtx)[order.dendrogram(sigFAFTLD_mtx.dendro)]
# pfamTPM_all$label_name <- factor(pfamTPM_all$label_name,levels = sigFAFTLD_mtx_order )
# 
# plt<-ggplot(pfamTPM_all,aes(x=condition, y=label_name)) +theme_classic()+
#   geom_tile(aes(fill=Zscore))+
#   scale_x_discrete(expand = c(0, 0))+ facet_wrap(~lightdark,scales="free_x")+
#   theme(axis.ticks.y=element_blank(),panel.spacing.x=unit(0, "lines"),axis.text.y = element_text(size = 5))+
#   scale_fill_viridis(option="inferno") + 
#   xlab("condition") + ylab("transcripts")+labs(fill='TPM Z-score')+ggtitle("")
# ggsave("aldex/SFR22_0829_heatmap_darkHits_Zscore_pless0.05.pdf", plot=plt,height=10, width=6)

#Is RPOB PF04563.18 cyclical?
mdT<-fread("metaT_metadata_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))
annot<-fread("pfam/pfam_annotationkey.csv")

stderror <- function(x) sd(x)/sqrt(length(x))

rpob_MT<-fread("pfam/pfam_clean_noNT_TPM.txt")%>%
  filter(FeatureID=="PF04563.18")%>%
  gather(sample_name,TPM,-FeatureID)%>%
  left_join(.,mdT,by="sample_name") %>%
  group_by(condition,zt_time)%>%mutate(meanAbun=mean(TPM),sem=stderror(TPM))%>%
  mutate(condition=factor(condition,levels=c("FT","FA","NA")))

p<-rpob_MT%>%
  ggplot(aes(x=zt_time, y=meanAbun, color=condition)) +
  geom_point(alpha=1.0) + geom_line() +
  theme_pubr() +
  scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+
  scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
  geom_ribbon(aes(ymin = meanAbun-sem, ymax = meanAbun+sem, fill=condition),alpha=0.5,colour = NA)+
  labs(color="condition",
       y ="MTX TPM",
       x ="ZT time")+ theme(plot.title = element_text(face = "bold"))
ggsave("MTvsMG/SFR22_0912_rpob_MTvsZT_err.pdf", plot=p,height=4, width=5)

mdG<-fread("/mnt/zarrinpar/scratch/sfloresr/metatranscript/anvi_metagenomic/woltka2_results/metaG_metadata_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))
annot<-fread("pfam/pfam_annotationkey.csv")

rpob_MG<-fread("/mnt/zarrinpar/scratch/sfloresr/metatranscript/anvi_metagenomic/woltka2_results/pfam/pfam_clean_noNT_TPM.txt")%>%
  filter(FeatureID=="PF04563.18")%>%
  gather(sample_name,TPM,-FeatureID)%>%
  left_join(.,mdG,by="sample_name") %>%
  group_by(condition,zt_time)%>%mutate(meanAbun=mean(TPM),sem=stderror(TPM))%>%
  mutate(condition=factor(condition,levels=c("FT","FA","NA")))

p<-rpob_MG%>%
  ggplot(aes(x=zt_time, y=meanAbun, color=condition)) +
  geom_point(alpha=1.0) + geom_line() +
  theme_pubr() +
  scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+
  scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
  geom_ribbon(aes(ymin = meanAbun-sem, ymax = meanAbun+sem, fill=condition),alpha=0.5,colour = NA)+
  labs(color="condition",
       y ="MGX TPM",
       x ="ZT time")+ theme(plot.title = element_text(face = "bold"))
ggsave("MTvsMG/SFR22_0911_rpob_MGvsZT_err.pdf", plot=p,height=4, width=5)

MTXcounts<-fread("pfam/pfam_clean_noNT.txt")%>%column_to_rownames("FeatureID")%>%round()%>%
  rownames_to_column("FeatureID")
MTX_deseq <- DESeqDataSetFromMatrix(countData=MTXcounts, 
                              colData=mdT, 
                              design=~condition, tidy = TRUE)
dds <- DESeq(MTX_deseq)
isControl <- rownames(dds) %in% "PF04563.18"
dds<- estimateSizeFactors(dds, controlGenes=isControl)
MTX_norm<-counts(dds, normalized=TRUE)%>%as.data.frame()%>%
  rownames_to_column("FeatureID")
write.table(MTX_norm,"pfam/pfam_normrpob_noNT.txt",sep = "\t",row.names = FALSE, quote=FALSE)

MGXcounts<-fread("/mnt/zarrinpar/scratch/sfloresr/metatranscript/anvi_metagenomic/woltka2_results/pfam/pfam_clean_noNT.txt")%>%column_to_rownames("FeatureID")%>%round()%>%
  rownames_to_column("FeatureID")
MGX_deseq <- DESeqDataSetFromMatrix(countData=MGXcounts, 
                                    colData=mdG, 
                                    design=~condition, tidy = TRUE)
ddsG <- DESeq(MGX_deseq)
isControl <- rownames(ddsG) %in% "PF04563.18"
ddsG<- estimateSizeFactors(ddsG, controlGenes=isControl)
MGX_norm<-counts(ddsG, normalized=TRUE)%>%as.data.frame()%>%
  rownames_to_column("FeatureID")
write.table(MGX_norm,"/mnt/zarrinpar/scratch/sfloresr/metatranscript/anvi_metagenomic/woltka2_results/pfam/pfam_normrpob_noNT.txt",sep = "\t",row.names = FALSE, quote=FALSE)
ggsave("MTvsMG/SFR22_0911_rpob_MTvsZT.pdf", plot=p,height=4, width=5)

#plot normalize abundance of sig hits mgx vs mtx (light)
mdT<-fread("metaT_metadata_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))
annot<-fread("pfam/pfam_annotationkey.csv")
lightuniquehits<-fread("aldex/SFR22_0829_FAFTL_20uniquehitsFAFTonly_pless0.05.txt")

TpfamTPM_sigL_zt<-fread("pfam/pfam_normrpob_noNT.txt")%>%
  filter(FeatureID %in% lightuniquehits$FeatureID)%>%
  gather(sample_name,norm_counts,-FeatureID) %>%
  left_join(.,mdT,by="sample_name") %>%
  filter(condition=="FT" & lightdark=="light")%>%
  group_by(FeatureID, condition, zt_time)%>%
  mutate(meanAbun=mean(norm_counts),str=sd(norm_counts)/sqrt(length(norm_counts)),
         t.score = qt(p=0.05/2, df=length(norm_counts)-1,lower.tail=F),
         margin.error=t.score*str)%>%
  mutate(condition=factor(condition,levels=c("FT","FA","NA")))
TpfamTPM_sigL_zt$FeatureID<-fct_reorder(TpfamTPM_sigL_zt$FeatureID, TpfamTPM_sigL_zt$meanAbun)

p<-ggplot(data=TpfamTPM_sigL_zt, aes(x=FeatureID, y=log10(meanAbun+1), fill=as.factor(zt_time)))+ 
  geom_bar(stat="identity") + theme_minimal()+
  ylim(0,25) +
  coord_flip() +ggtitle("TRF unique (light) MTX expression") #position=position_dodge()
ggsave("MTvsMG/SFR22_0911_MTX_FTL_barplots_zt.pdf", plot=p,height=4, width=5)

TpfamTPM_sigL<-fread("pfam/pfam_normrpob_noNT.txt")%>%
  filter(FeatureID %in% lightuniquehits$FeatureID)%>%
  gather(sample_name,norm_counts,-FeatureID) %>%
  left_join(.,mdT,by="sample_name") %>%
  filter(condition=="FT" & lightdark=="light")%>%
  group_by(FeatureID, condition)%>%
  mutate(meanAbun=mean(norm_counts),str=sd(norm_counts)/sqrt(length(norm_counts)),
         t.score = qt(p=0.05/2, df=length(norm_counts)-1,lower.tail=F),
         margin.error=t.score*str)%>%
  mutate(condition=factor(condition,levels=c("FT","FA","NA")))
TpfamTPM_sigL$FeatureID<-fct_reorder(TpfamTPM_sigL$FeatureID, TpfamTPM_sigL$meanAbun)

p<-ggplot(data=TpfamTPM_sigL, aes(x=reorder(FeatureID, log10(meanAbun+1)), y=log10(meanAbun+1)))+ 
  geom_bar(stat="identity",fill="grey60") + theme_minimal()+
  ylim(0,25) +
  coord_flip() +ggtitle("TRF unique (light) MTX expression") #position=position_dodge()
ggsave("MTvsMG/SFR22_0911_MTX_FTL_barplots.pdf", plot=p,height=4, width=5)


mdG<-fread("/mnt/zarrinpar/scratch/sfloresr/metatranscript/anvi_metagenomic/woltka2_results/metaG_metadata_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))
annot<-fread("pfam/pfam_annotationkey.csv")
lightuniquehits<-fread("aldex/SFR22_0829_FAFTL_20uniquehitsFAFTonly_pless0.05.txt")

MpfamTPM_sigL_zt<-fread("/mnt/zarrinpar/scratch/sfloresr/metatranscript/anvi_metagenomic/woltka2_results/pfam/pfam_normrpob_noNT.txt")%>%
  filter(FeatureID %in% lightuniquehits$FeatureID)%>%
  gather(sample_name,norm_counts,-FeatureID) %>%
  left_join(.,mdG,by="sample_name") %>%
  filter(condition=="FT" & lightdark=="light")%>%
  group_by(FeatureID, condition, zt_time)%>%
  mutate(meanAbun=mean(norm_counts),str=sd(norm_counts)/sqrt(length(norm_counts)),
         t.score = qt(p=0.05/2, df=length(norm_counts)-1,lower.tail=F),
         margin.error=t.score*str)%>%
  mutate(condition=factor(condition,levels=c("FT","FA","NA")),
         FeatureID=factor(FeatureID,levels=levels(TpfamTPM_sigL_zt$FeatureID)))

p<-ggplot(data=MpfamTPM_sigL_zt, aes(x=FeatureID, y=log10(meanAbun+1), fill=as.factor(zt_time)))+ 
  geom_bar(stat="identity") + theme_minimal()+
  ylim(0,25) +
  coord_flip() +ggtitle("TRF unique (light) MGX expression") #position=position_dodge()
ggsave("MTvsMG/SFR22_0911_MGX_FTL_barplots_zt.pdf", plot=p,height=4, width=5)

MpfamTPM_sigL<-fread("/mnt/zarrinpar/scratch/sfloresr/metatranscript/anvi_metagenomic/woltka2_results/pfam/pfam_normrpob_noNT.txt")%>%
  filter(FeatureID %in% lightuniquehits$FeatureID)%>%
  gather(sample_name,norm_counts,-FeatureID) %>%
  left_join(.,mdG,by="sample_name") %>%
  filter(condition=="FT" & lightdark=="light")%>%
  group_by(FeatureID, condition)%>%
  mutate(meanAbun=mean(norm_counts),str=sd(norm_counts)/sqrt(length(norm_counts)),
         t.score = qt(p=0.05/2, df=length(norm_counts)-1,lower.tail=F),
         margin.error=t.score*str)%>%
  mutate(condition=factor(condition,levels=c("FT","FA","NA")),
         FeatureID=factor(FeatureID,levels=levels(TpfamTPM_sigL$FeatureID)))

p<-ggplot(data=MpfamTPM_sigL, aes(x=FeatureID, y=log10(meanAbun+1)))+ 
  geom_bar(stat="identity",fill="grey60") + theme_minimal()+
  ylim(0,25) +
  coord_flip() +ggtitle("TRF unique (light) MGX expression") #position=position_dodge()
ggsave("MTvsMG/SFR22_0911_MGX_FTL_barplots.pdf", plot=p,height=4, width=5)

#plot normalize abundance of sig hits mgx vs mtx (dark)
mdT<-fread("metaT_metadata_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))
annot<-fread("pfam/pfam_annotationkey.csv")
darkuniquehits<-fread("aldex/SFR22_0829_FAFTD_118uniquehitsFAFTonly_pless0.05.txt")

TpfamTPM_sigD_zt<-fread("pfam/pfam_normrpob_noNT.txt")%>%
  filter(FeatureID %in% darkuniquehits$FeatureID)%>%
  gather(sample_name,norm_counts,-FeatureID) %>%
  left_join(.,mdT,by="sample_name") %>%
  filter(condition=="FT" & lightdark=="light")%>%
  group_by(FeatureID, condition, zt_time)%>%
  mutate(meanAbun=mean(norm_counts),str=sd(norm_counts)/sqrt(length(norm_counts)),
         t.score = qt(p=0.05/2, df=length(norm_counts)-1,lower.tail=F),
         margin.error=t.score*str)%>%
  mutate(condition=factor(condition,levels=c("FT","FA","NA")))
TpfamTPM_sigD_zt$FeatureID<-fct_reorder(TpfamTPM_sigD_zt$FeatureID, TpfamTPM_sigD_zt$meanAbun)

p<-ggplot(data=TpfamTPM_sigD_zt, aes(x=FeatureID, y=log10(meanAbun+1), fill=as.factor(zt_time)))+ 
  geom_bar(stat="identity") + theme_minimal()+ theme(axis.text.y = element_text(size = 5)) +
  ylim(0,25) +
  coord_flip() +ggtitle("TRF unique (dark) MTX expression") #position=position_dodge()
ggsave("MTvsMG/SFR22_0911_MTX_FTD_barplots_zt.pdf", plot=p,height=10, width=5)

TpfamTPM_sigD<-fread("pfam/pfam_normrpob_noNT.txt")%>%
  filter(FeatureID %in% darkuniquehits$FeatureID)%>%
  gather(sample_name,norm_counts,-FeatureID) %>%
  left_join(.,mdT,by="sample_name") %>%
  filter(condition=="FT" & lightdark=="light")%>%
  group_by(FeatureID, condition)%>%
  mutate(meanAbun=mean(norm_counts),str=sd(norm_counts)/sqrt(length(norm_counts)),
         t.score = qt(p=0.05/2, df=length(norm_counts)-1,lower.tail=F),
         margin.error=t.score*str)%>%
  mutate(condition=factor(condition,levels=c("FT","FA","NA")))
TpfamTPM_sigD$FeatureID<-fct_reorder(TpfamTPM_sigD$FeatureID, TpfamTPM_sigD$meanAbun)

p<-ggplot(data=TpfamTPM_sigD, aes(x=reorder(FeatureID, log10(meanAbun+1)), y=log10(meanAbun+1)))+ 
  geom_bar(stat="identity",fill="grey10") + theme_minimal()+ theme(axis.text.y = element_text(size = 5)) +
  ylim(0,25) +
  coord_flip() +ggtitle("TRF unique (dark) MTX expression") #position=position_dodge()
ggsave("MTvsMG/SFR22_0911_MTX_FTD_barplots.pdf", plot=p,height=10, width=5)


mdG<-fread("/mnt/zarrinpar/scratch/sfloresr/metatranscript/anvi_metagenomic/woltka2_results/metaG_metadata_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))
annot<-fread("pfam/pfam_annotationkey.csv")
darkuniquehits<-fread("aldex/SFR22_0829_FAFTD_118uniquehitsFAFTonly_pless0.05.txt")

MpfamTPM_sigD_zt<-fread("/mnt/zarrinpar/scratch/sfloresr/metatranscript/anvi_metagenomic/woltka2_results/pfam/pfam_normrpob_noNT.txt")%>%
  filter(FeatureID %in% darkuniquehits$FeatureID)%>%
  gather(sample_name,norm_counts,-FeatureID) %>%
  left_join(.,mdG,by="sample_name") %>%
  filter(condition=="FT" & lightdark=="light")%>%
  group_by(FeatureID, condition, zt_time)%>%
  mutate(meanAbun=mean(norm_counts),str=sd(norm_counts)/sqrt(length(norm_counts)),
         t.score = qt(p=0.05/2, df=length(norm_counts)-1,lower.tail=F),
         margin.error=t.score*str)%>%
  mutate(condition=factor(condition,levels=c("FT","FA","NA")),
         FeatureID=factor(FeatureID,levels=levels(TpfamTPM_sigD_zt$FeatureID)))

p<-ggplot(data=MpfamTPM_sigD_zt, aes(x=FeatureID, y=log10(meanAbun+1), fill=as.factor(zt_time)))+ 
  geom_bar(stat="identity") + theme_minimal()+theme(axis.text.y = element_text(size = 5)) +
  ylim(0,25) +
  coord_flip() +ggtitle("TRF unique (dark) MGX expression") #position=position_dodge()
ggsave("MTvsMG/SFR22_0911_MGX_FTD_barplots_zt.pdf", plot=p,height=10, width=5)

MpfamTPM_sigD<-fread("/mnt/zarrinpar/scratch/sfloresr/metatranscript/anvi_metagenomic/woltka2_results/pfam/pfam_normrpob_noNT.txt")%>%
  filter(FeatureID %in% darkuniquehits$FeatureID)%>%
  gather(sample_name,norm_counts,-FeatureID) %>%
  left_join(.,mdG,by="sample_name") %>%
  filter(condition=="FT" & lightdark=="light")%>%
  group_by(FeatureID, condition)%>%
  mutate(meanAbun=mean(norm_counts),str=sd(norm_counts)/sqrt(length(norm_counts)),
         t.score = qt(p=0.05/2, df=length(norm_counts)-1,lower.tail=F),
         margin.error=t.score*str)%>%
  mutate(condition=factor(condition,levels=c("FT","FA","NA")),
         FeatureID=factor(FeatureID,levels=levels(TpfamTPM_sigD$FeatureID)))

p<-ggplot(data=MpfamTPM_sigD, aes(x=FeatureID, y=log10(meanAbun+1)))+ 
  geom_bar(stat="identity",fill="grey10") + theme_minimal()+ theme(axis.text.y = element_text(size = 5)) +
  ylim(0,25) +
  coord_flip() +ggtitle("TRF unique (light) MGX expression") #position=position_dodge()
ggsave("MTvsMG/SFR22_0911_MGX_FTD_barplots.pdf", plot=p,height=10, width=5)

#metaT vs metaG for sig hits (light)
mdT<-fread("metaT_metadata_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))
annot<-fread("pfam/pfam_annotationkey.csv")

TpfamTPM_sigL<-fread("pfam/pfam_normrpob_noNT.txt")%>%
  filter(FeatureID %in% lightuniquehits$FeatureID)%>%
  gather(sample_name,TPM_counts,-FeatureID) %>%
  left_join(.,mdT,by="sample_name") %>%
  filter(condition=="FT" & lightdark=="light")%>%
  dplyr::select(FeatureID,sample_name,TPM_counts)%>%
  dplyr::rename(MT_TPM_counts=TPM_counts)

mdM<-fread("/mnt/zarrinpar/scratch/sfloresr/metatranscript/anvi_metagenomic/woltka2_results/metaG_metadata_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))
annot<-fread("pfam/pfam_annotationkey.csv")

MpfamTPM_sigL<-fread("/mnt/zarrinpar/scratch/sfloresr/metatranscript/anvi_metagenomic/woltka2_results/pfam/pfam_normrpob_noNT.txt")%>%
  filter(FeatureID %in% lightuniquehits$FeatureID)%>%
  gather(sample_name,TPM_counts,-FeatureID) %>%
  left_join(.,mdM,by="sample_name") %>%
  mutate(sample_name=ifelse(zt_time<13,
                     paste("c",condition,"0",zt_time,replicate, sep=""),
                     paste("c",condition,zt_time,replicate, sep=""))) %>%
  filter(condition=="FT" & lightdark=="light")%>%
  dplyr::select(FeatureID,sample_name,TPM_counts)%>%
  dplyr::rename(MG_TPM_counts=TPM_counts)

#only have matching samples for replicate b
TpMpTPM<-merge(TpfamTPM_sigL,MpfamTPM_sigL,by=c("FeatureID","sample_name"))%>%
  left_join(.,mdT,by="sample_name") %>%
  filter(MG_TPM_counts>0)


p<-ggplot(TpMpTPM, aes(x=log10(MG_TPM_counts+1), y=log10(MT_TPM_counts+1),
                       colour=as.factor(zt_time))) + geom_point() +
  geom_abline(intercept = 0, slope = 1, color="black",
              linetype="dashed")+
  geom_text(label=TpMpTPM$FeatureID,hjust = 0, nudge_x = 0.05)+
  theme_light()+ggtitle("FAFTL")
ggsave("aldex/SFR22_0911_xyplot_lightHits_MTvsMG_pless0.05_rpob.pdf", plot=p,height=4, width=6)

p<-ggplot(TpMpTPM, aes(x=log10(MG_TPM_counts+1), y=log10(MT_TPM_counts+1),
                       colour=FeatureID)) + geom_point() +
  geom_abline(intercept = 0, slope = 1, color="black",
              linetype="dashed")+
  geom_text(label=TpMpTPM$zt_time,hjust = 0, nudge_x = 0.05)+
  theme_light()+ggtitle("FAFTL") + theme(legend.key.size = unit(0, 'lines'))
ggsave("aldex/SFR22_0911_xyplot_lightHits_MTvsMG_samplN_pless0.05_rpob.pdf", plot=p,height=4, width=6)

#metaT vs metaG for sig hits (dark)
mdT<-fread("metaT_metadata_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))
annot<-fread("pfam/pfam_annotationkey.csv")

TpfamTPM_sigD<-fread("pfam/pfam_normrpob_noNT.txt")%>%
  filter(FeatureID %in% darkuniquehits$FeatureID)%>%
  gather(sample_name,TPM_counts,-FeatureID) %>%
  left_join(.,mdT,by="sample_name") %>%
  filter(condition=="FT" & lightdark=="dark")%>%
  dplyr::select(FeatureID,sample_name,TPM_counts)%>%
  dplyr::rename(MT_TPM_counts=TPM_counts)

mdM<-fread("/mnt/zarrinpar/scratch/sfloresr/metatranscript/anvi_metagenomic/woltka2_results/metaG_metadata_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))
annot<-fread("pfam/pfam_annotationkey.csv")

MpfamTPM_sigD<-fread("/mnt/zarrinpar/scratch/sfloresr/metatranscript/anvi_metagenomic/woltka2_results/pfam/pfam_normrpob_noNT.txt")%>%
  filter(FeatureID %in% darkuniquehits$FeatureID)%>%
  gather(sample_name,TPM_counts,-FeatureID) %>%
  left_join(.,mdM,by="sample_name") %>%
  mutate(sample_name=ifelse(zt_time<13,
                            paste("c",condition,"0",zt_time,replicate, sep=""),
                            paste("c",condition,zt_time,replicate, sep=""))) %>%
  filter(condition=="FT" & lightdark=="dark")%>%
  dplyr::select(FeatureID,sample_name,TPM_counts)%>%
  dplyr::rename(MG_TPM_counts=TPM_counts)

#only have matching samples for replicate b
TpMpTPM<-merge(TpfamTPM_sigD,MpfamTPM_sigD,by=c("FeatureID","sample_name"))%>%
  left_join(.,mdT,by="sample_name") %>%
  filter(MG_TPM_counts>0)

p<-ggplot(TpMpTPM, aes(x=log10(MG_TPM_counts+1), y=log10(MT_TPM_counts+1),
                       colour=as.factor(zt_time))) + geom_point() +
  geom_abline(intercept = 0, slope = 1, color="black",
              linetype="dashed")+
  geom_text(data=subset(TpMpTPM, (log10(MG_TPM_counts+1) < 1.5 & log10(MT_TPM_counts+1) < 1)| (log10(MG_TPM_counts+1) < 2 & log10(MT_TPM_counts+1) > 3)),
            aes(log10(MG_TPM_counts+1),log10(MT_TPM_counts+1),label=FeatureID))+
  theme_light()+ggtitle("FAFTD")
ggsave("aldex/SFR22_0911_xyplot_darkHits_MTvsMG_pless0.05_rpob.pdf", plot=p,height=4, width=6)

p<-ggplot(TpMpTPM, aes(x=log10(MG_TPM_counts+1), y=log10(MT_TPM_counts+1),
                       colour=FeatureID)) + geom_point() +
  geom_abline(intercept = 0, slope = 1, color="black",
              linetype="dashed")+
  geom_text(label=TpMpTPM$zt_time,hjust = 0, nudge_x = 0.05)+
  theme_light()+ggtitle("FAFTD") + theme(legend.key.size = unit(0, 'lines'))
ggsave("aldex/SFR22_0911_xyplot_darkHits_MTvsMG_samplN_pless0.05_rpob.pdf", plot=p,height=6, width=13)

#just BSH

BSHpfamTPM<-fread("pfam/pfam_normrpob_noNT.txt")%>%
  filter(FeatureID=="PF02275.21")%>% 
  gather(sample_name,TPM_counts,-FeatureID) %>%
  left_join(.,mdT,by="sample_name") %>%
  #filter(condition=="FT" & lightdark=="dark")%>%
  dplyr::select(FeatureID,sample_name,TPM_counts)%>%
  dplyr::rename(MT_TPM_counts=TPM_counts)

BSH_mgpfamTPM<-fread("/mnt/zarrinpar/scratch/sfloresr/metatranscript/anvi_metagenomic/woltka2_results/pfam/pfam_normrpob_noNT.txt")%>%
  filter(FeatureID=="PF02275.21")%>% 
  gather(sample_name,TPM_counts,-FeatureID) %>%
  left_join(.,mdM,by="sample_name") %>%
  mutate(sample_name=ifelse(zt_time<13,
                            paste("c",condition,"0",zt_time,replicate, sep=""),
                            paste("c",condition,zt_time,replicate, sep=""))) %>%
  #filter(condition=="FT" & lightdark=="dark")%>%
  dplyr::select(FeatureID,sample_name,TPM_counts)%>%
  dplyr::rename(MG_TPM_counts=TPM_counts)

BSH_TpMpTPM<-merge(BSHpfamTPM,BSH_mgpfamTPM,by=c("FeatureID","sample_name"))%>%
  left_join(.,mdT,by="sample_name")%>%
  mutate(condition=factor(condition, levels = c("FT", "FA","NA")),
         lightdark=factor(lightdark, levels = c("light","dark")))

p<-ggplot(BSH_TpMpTPM, aes(x=log10(MG_TPM_counts+1), y=log10(MT_TPM_counts+1), 
                       colour=condition,shape=lightdark)) + geom_point() +
  scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+
  geom_abline(intercept = 0, slope = 1, color="black", 
              linetype="dashed")+
  scale_shape_manual(values=c(3,16)) +
  geom_text(label=BSH_TpMpTPM$zt_time,hjust = 0, nudge_x = 0.05)+
  theme_light()+ggtitle("BSH")
ggsave("BSH/SFR22_0911_xyplot_BSH_MTvsMG.pdf", plot=p,height=4, width=6)


p<-ggplot(BSH_TpMpTPM, aes(x=condition, y=log10(MT_TPM_counts+1), fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
  theme_classic()+
  labs(x="condition",y="log10(relative abundance)", title="PF02275.21:BSH MT")+
  theme(legend.position = "none")
ggsave("BSH/SFR22_0911_BSHMT_boxplot.pdf", plot=p,height=2.5, width=3.5)

pairwise.wilcox.test(BSH_TpMpTPM$MT_TPM_counts, BSH_TpMpTPM$condition,p.adjust.method="fdr")

# FT      FA     
# FA 0.39    -      
#   NA 2.2e-05 2.2e-05

p<-ggplot(BSH_TpMpTPM, aes(x=condition, y=log10(MT_TPM_counts+1), fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
  facet_wrap(~lightdark)+
  theme_classic()+
  labs(x="condition",y="log10(relative abundance)", title="PF02275.21:BSH MT")+
  theme(legend.position = "none")
ggsave("BSH/SFR22_1012_BSHMT_boxplot.pdf", plot=p,height=2.5, width=3.5)

light<-BSH_TpMpTPM%>%filter(lightdark=="light")
pairwise.wilcox.test(light$MT_TPM_counts, light$condition,p.adjust.method="fdr")

# FT    FA   
# FA 0.200 -    
#   NA 0.014 0.014

dark<-BSH_TpMpTPM%>%filter(lightdark=="dark")
pairwise.wilcox.test(dark$MT_TPM_counts, dark$condition,p.adjust.method="fdr")

# FT    FA   
# FA 1.000 -    
#   NA 0.014 0.014

FT<-BSH_TpMpTPM%>%filter(condition=="FT")
pairwise.wilcox.test(FT$MT_TPM_counts, FT$lightdark,p.adjust.method="fdr")

# light
# dark 0.7

FA<-BSH_TpMpTPM%>%filter(condition=="FA")
pairwise.wilcox.test(FA$MT_TPM_counts, FA$lightdark,p.adjust.method="fdr")

# light
# dark 0.1

pNA<-BSH_TpMpTPM%>%filter(condition=="NA")
pairwise.wilcox.test(pNA$MT_TPM_counts, pNA$lightdark,p.adjust.method="fdr")

# light
# dark 0.094

p<-ggplot(BSH_TpMpTPM, aes(x=condition, y=log10(MG_TPM_counts+1), fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
  theme_classic()+
  labs(x="condition",y="log10(relative abundance)", title="PF02275.21:BSH MG")+
  theme(legend.position = "none")
ggsave("BSH/SFR22_0911_BSHMG_boxplot.pdf", plot=p,height=2.5, width=3.5)

pairwise.wilcox.test(BSH_TpMpTPM$MG_TPM_counts, BSH_TpMpTPM$condition,p.adjust.method="fdr")

# FT   FA  
# FA 0.63 -   
#   NA 0.26 0.26

#TPM of all sig hits over ZT time

salmonZT<-fread("qiime/rmdoubletons_noNT_metaT/salmonTPM_normcounts_rmdoubletons_noNT.txt")%>%
  filter( Name %in% sigFAFTLD$`OTU ID`)%>% 
  gather(sample_name,TPM_counts,-Name) %>%
  left_join(.,annot, by ="Name") %>%
  mutate(label_name=paste(Name, annotation, sep=" "))%>%
  left_join(.,md,by="sample_name") %>%
  group_by(Name,condition,zt_time)%>%summarise(mn_TPM=mean(TPM_counts))%>%
  mutate(condition=factor(condition, levels = c("FT", "FA","NA")))

#example<-salmonZT%>%filter(Name=="EIFNIBFL_171339")

ZT_dist <- function(data,Name) {
  
  data%>%
    ggplot(aes(x=zt_time, y=log10(mn_TPM+1), color=condition)) +
    geom_point(alpha=1.0) + geom_line() +
    theme_pubr() +
    scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+
    scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
    scale_x_continuous(breaks=c(1,5,9,13,17,21))+
    #geom_ribbon(aes(ymin = meanAbun-margin.error, ymax = meanAbun+margin.error, fill=condition),alpha=0.5,colour = NA)+
    labs(color="condition",
         y ="avg log10TPM",
         x ="ZT time",
         title=Name)+ theme(plot.title = element_text(face = "bold"))
  ggsave(paste("aldex/sighits_vsZT/",Name,".pdf",sep=""), height=3, width=3)
}

salmonZT_nested <- salmonZT %>% 
  group_by(Name) %>% 
  nest()

salmonZT_plots <- 
  salmonZT_nested %>% 
  mutate(plot = map2(data, Name,  ~ ZT_dist(.x,.y)))

#FA vs NA
condsFAvsNA<-conds[c(1:12,25:42)]
salmonFANA<-fread("qiime/rmdoubletons_noNT_metaT/songbird/FANA/FANA_salmonTPM_normcounts_rmdoubletons_noNT.txt")%>%
  column_to_rownames("Name")
FANA.ald<-aldex.clr(round(salmonFANA),condsFAvsNA, mc.samples=500, denom="all", verbose=F)

FANA.ttest<-aldex.ttest(FANA.ald,paired.test = FALSE, hist.plot=FALSE)%>%rownames_to_column("OTU ID")
FANA.effect<-aldex.effect(FANA.ald)%>%rownames_to_column("OTU ID")%>%left_join(.,FANA.ttest,by="OTU ID")
write.table(FANA.effect,"aldex/SFR22_0520_FANA_ald_effectwpval.txt",sep = "\t",row.names = FALSE, quote=FALSE)

aldex.plot(FANA.effect, type="MA", test="wilcox",cutoff.pval=0.1)
aldex.plot(FANA.effect, type="MW", test="wilcox",cutoff.pval=0.1)

#FA vs NA (light)
md<-fread("/mnt/zarrinpar/scratch/sfloresr/metatranscript/metatranscript/trf_cecal_metadata_noNT.tsv")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition)) %>%
  filter(zt_time<13)
conds <- md$condition

condsFAvsNA<-conds[c(1:6,13:21)]
salmonFANA<-fread("qiime/rmdoubletons_noNT_metaT/songbird/FANA/FANA_salmonTPM_normcounts_rmdoubletons_noNT.txt")%>%
  column_to_rownames("Name") %>%dplyr::select(c(1:6,13:21))
FANAL.ald<-aldex.clr(round(salmonFANA),condsFAvsNA, mc.samples=500, denom="all", verbose=F)

FANAL.ttest<-aldex.ttest(FANAL.ald,paired.test = FALSE, hist.plot=FALSE)%>%rownames_to_column("OTU ID")
FANAL.effect<-aldex.effect(FANAL.ald)%>%rownames_to_column("OTU ID")%>%left_join(.,FANAL.ttest,by="OTU ID")
write.table(FANAL.effect,"aldex/SFR22_0628_FANAL_ald_effectwpval.txt",sep = "\t",row.names = FALSE, quote=FALSE)

FANAL.effect<-fread("aldex/SFR22_0622_FANAL_ald_effectwpval.txt")

aldex.plot(FANAL.effect, type="MA", test="wilcox",cutoff.pval=0.1)
aldex.plot(FANAL.effect, type="MW", test="wilcox",cutoff.pval=0.1)

#FA vs NA (dark)
md<-fread("/mnt/zarrinpar/scratch/sfloresr/metatranscript/metatranscript/trf_cecal_metadata_noNT.tsv")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition)) %>%
  filter(zt_time>9)
conds <- md$condition

condsFAvsNA<-conds[c(1:6,13:21)]
salmonFANA<-fread("qiime/rmdoubletons_noNT_metaT/songbird/FANA/FANA_salmonTPM_normcounts_rmdoubletons_noNT.txt")%>%
  column_to_rownames("Name") %>%dplyr::select(c(7:12,22:30))
FANAD.ald<-aldex.clr(round(salmonFANA),condsFAvsNA, mc.samples=500, denom="all", verbose=F)

FANAD.ttest<-aldex.ttest(FANAD.ald,paired.test = FALSE, hist.plot=FALSE)%>%rownames_to_column("OTU ID")
FANAD.effect<-aldex.effect(FANAD.ald)%>%rownames_to_column("OTU ID")%>%left_join(.,FANAD.ttest,by="OTU ID")
write.table(FANAD.effect,"aldex/SFR22_0628_FANAD_ald_effectwpval.txt",sep = "\t",row.names = FALSE, quote=FALSE)

aldex.plot(FANAD.effect, type="MA", test="wilcox",cutoff.pval=0.1)
aldex.plot(FANAD.effect, type="MW", test="wilcox",cutoff.pval=0.1)

#FT vs NA
condsFTvsNA<-conds[c(13:42)]
salmonFTNA<-fread("qiime/rmdoubletons_noNT_metaT/songbird/FTNA/FTNA_salmonTPM_normcounts_rmdoubletons_noNT.txt")%>%
  column_to_rownames("Name")
FTNA.ald<-aldex.clr(round(salmonFTNA),condsFTvsNA, mc.samples=500, denom="all", verbose=F)

FTNA.ttest<-aldex.ttest(FTNA.ald,paired.test = FALSE, hist.plot=FALSE)%>%rownames_to_column("OTU ID")
FTNA.effect<-aldex.effect(FTNA.ald)%>%rownames_to_column("OTU ID")%>%left_join(.,FTNA.ttest,by="OTU ID")
write.table(FTNA.effect,"aldex/SFR22_0528_FANAD_ald_effectwpval.txt",sep = "\t",row.names = FALSE, quote=FALSE)

aldex.plot(FTNA.effect, type="MA", test="wilcox",cutoff.pval=0.1)
aldex.plot(FTNA.effect, type="MW", test="wilcox",cutoff.pval=0.1)

#FT vs NA (light)
md<-fread("/mnt/zarrinpar/scratch/sfloresr/metatranscript/metatranscript/trf_cecal_metadata_noNT.tsv")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition)) %>%
  filter(zt_time<13)
conds <- md$condition

condsFTvsNA<-conds[c(7:21)]
salmonFTNA<-fread("qiime/rmdoubletons_noNT_metaT/songbird/FTNA/FTNA_salmonTPM_normcounts_rmdoubletons_noNT.txt")%>%
  column_to_rownames("Name") %>%dplyr::select(c(1:6,13:21))
FTNAL.ald<-aldex.clr(round(salmonFTNA),condsFTvsNA, mc.samples=500, denom="all", verbose=F)

FTNAL.ttest<-aldex.ttest(FTNAL.ald,paired.test = FALSE, hist.plot=FALSE)%>%rownames_to_column("OTU ID")
FTNAL.effect<-aldex.effect(FTNAL.ald)%>%rownames_to_column("OTU ID")%>%left_join(.,FTNAL.ttest,by="OTU ID")
write.table(FTNAL.effect,"aldex/SFR22_0628_FTNAL_ald_effectwpval.txt",sep = "\t",row.names = FALSE, quote=FALSE)

FANAL.effect<-fread("aldex/SFR22_0628_FTNAL_ald_effectwpval.txt")

aldex.plot(FANAL.effect, type="MA", test="wilcox",cutoff.pval=0.1)
aldex.plot(FANAL.effect, type="MW", test="wilcox",cutoff.pval=0.1)

#light vs dark FT
md<-fread("/mnt/zarrinpar/scratch/sfloresr/metatranscript/metatranscript/trf_cecal_metadata_noNT.tsv")%>%
  filter(condition=="FT") %>%
  #mutate(condition=ifelse(is.na(condition),"NA",condition)) %>%
  mutate(lightdark=ifelse(zt_time<13,"light","dark")) 
conds <- md$lightdark

salmonFT<-fread("qiime/rmdoubletons_noNT_metaT/songbird/FAFT/FAFT_salmonTPM_normcounts_rmdoubletons_noNT.txt")%>%
  column_to_rownames("Name") %>%dplyr::select(c(13:24))
FTLD.ald<-aldex.clr(round(salmonFT),conds, mc.samples=500, denom="all", verbose=F)

FTLD.ttest<-aldex.ttest(FTLD.ald,paired.test = FALSE, hist.plot=FALSE)%>%rownames_to_column("OTU ID")
FTLD.effect<-aldex.effect(FTLD.ald)%>%rownames_to_column("OTU ID")%>%left_join(.,FTLD.ttest,by="OTU ID")
write.table(FTLD.effect,"aldex/SFR22_0524_FTlighdark_ald_effectwpval.txt",sep = "\t",row.names = FALSE, quote=FALSE)

FTLD.effect<-fread("aldex/SFR22_0524_FTlighdark_ald_effectwpval.txt")

annot<-fread("annotationKey.csv") %>%
  dplyr::rename(`OTU ID`=CycID)
mult_annot<-fread("prot_search/PROKKA_09232021_microAnnot/annotation_results/PROKKA_09232021.faa.annot")%>%
  dplyr::rename(`OTU ID`=query_id)
# GO_annot<-fread("/mnt/zarrinpar/scratch/sfloresr/metatranscript/metatranscript/woltka_results/species_process/name.txt")%>%
#   dplyr::rename(process_go="GO term")

FTLD.effect.annot<-FTLD.effect %>% left_join(.,annot, by ="OTU ID") %>%
  left_join(.,mult_annot, by ="OTU ID") %>%
  mutate(diffexpr=ifelse(wi.ep<0.05 & diff.btw>0, "up",
                         ifelse(wi.ep<0.05 & diff.btw< -0,"down","none"))) %>%
  filter(!str_detect(annotation,"hypothetical|ribosomal|Uncharac"))

p<-ggplot(data=FTLD.effect.annot, aes(x=diff.btw, y=-log10(wi.ep),color=diffexpr)) +
  ggtitle("FT (dark vs. light)")+
  geom_point(size = 0.5) + 
  theme_minimal() +
  #geom_text_repel(data = subset(deseq.faft.mt, diffexpr!="none" & annotation!="hypothetical protein"), 
  #                aes(label = annotation), size=1)+
  theme(legend.position = "none") +
  scale_color_manual(values=c("blue", "black", "red"))# +
#geom_vline(xintercept=c(-3.5, 3.5), linetype = "dashed", col="grey36") +
#geom_hline(yintercept=-log10(0.05), linetype = "dashed",col="grey36")

ggsave("aldex/SFR22_0524_FTLD_volcanoplot.png", width = 3.5, height = 3)

sigFTLD<-FTLD.effect.annot%>% filter(wi.ep<0.05)
#11 up (light)
#147 down (dark)


# Basic piechart

colourCount = length(unique(sigFAFTL_goClean$annotation.y))
getPalette = colorRampPalette(brewer.pal(9, "Dark2"))

ggplot(sigFAFTL_goClean, aes(x="", y=n, fill=fct_reorder(annotation.y, n)))+
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + guides(fill=guide_legend(title="GO process")) +
  theme(legend.position="right",
        legend.title = element_text(size = 5, face="bold"), 
        legend.text = element_text(size = 5),
        legend.key.size = unit(0.5,"line")) +
  scale_fill_manual(values = getPalette(colourCount))


ggsave("aldex/SFR22_0602_FAFTD_GOpless0.01.pdf", width = 3, height = 2)

# p<-ggplot(data=FAFT.effect.annot, aes(x=diff.btw, y=-log10(wi.ep),color=diffexpr)) +
#   ggtitle("FA vs. FT")+
#   geom_point(size = 0.5) +
#   ylim(0,4) + xlim(-6,8) +
#   theme_minimal() +
#   theme(legend.position = "none") +
#   scale_color_manual(values=c("blue", "black", "red")) +
#geom_vline(xintercept=c(-3.5, 3.5), linetype = "dashed", col="grey36") +
#geom_hline(yintercept=-log10(0.05), linetype = "dashed",col="grey36")

#ggsave("aldex/SFR23_0110_FAFT_noDUF_volcanoplot_pless0.05.pdf", width = 3.5, height = 3)

# p<-ggplot(data=FAFTL.effect.annot, aes(x=diff.btw, y=-log10(wi.ep),color=diffexpr)) +
#   ggtitle("FA vs. FT (Light)")+
#   geom_point(size = 0.5) + 
#   ylim(0,4) + xlim(-6,8) +
#   theme_minimal() +
#   theme(legend.position = "none") +
#   scale_color_manual(values=c("blue", "black", "red")) +
#geom_vline(xintercept=c(-3.5, 3.5), linetype = "dashed", col="grey36") +
#geom_hline(yintercept=-log10(0.05), linetype = "dashed",col="grey36")

#ggsave("aldex/SFR23_0110_FAFTL_noDUF_volcanoplot_pless0.05.pdf", width = 3.5, height = 3)

# p<-ggplot(data=FAFTD.effect.annot, aes(x=diff.btw, y=-log10(wi.ep),color=diffexpr)) +
#   ggtitle("FA vs. FT (Dark)")+
#   geom_point(size = 0.5) + 
#   ylim(0,4) + xlim(-6,8) +
#   theme_minimal() +
#   theme(legend.position = "none") +
#   scale_color_manual(values=c("blue", "black", "red")) +
#geom_vline(xintercept=c(-3.5, 3.5), linetype = "dashed", col="grey36") +
#geom_hline(yintercept=-log10(0.05), linetype = "dashed",col="grey36")

#ggsave("aldex/SFR23_0110_FAFTD_noDUF_volcanoplot_pless0.05.pdf", width = 3.5, height = 3)

# list_venn <- list(light = (FAFTL.effect.annot%>%filter(wi.ep<0.05))$FeatureID,
#                   dark = (FAFTD.effect.annot%>%filter(wi.ep<0.05))$FeatureID,
#                   all = (FAFT.effect.annot%>%filter(wi.ep<0.05))$FeatureID)
# 
# ggvenn(list_venn, c("light", "dark", "all"),
#        fill_color = c("white", "black", "gray"))
# ggsave("aldex/SFR22_0824_FAFTLDvenn_pless0.05.pdf", width = 3.5, height = 3)

# p<-ggplot(data=FTNA.effect.annot, aes(x=diff.btw, y=-log10(wi.ep),color=diffexpr)) +
#   ggtitle("NA vs. FT")+
#   geom_point(size = 0.5) +
#   ylim(0,8) + xlim(-8,8) +
#   theme_minimal() +
#   theme(legend.position = "none") +
#   scale_color_manual(values=c("blue", "black", "red")) +
#geom_vline(xintercept=c(-3.5, 3.5), linetype = "dashed", col="grey36") +
#geom_hline(yintercept=-log10(0.05), linetype = "dashed",col="grey36")

#ggsave("aldex/SFR23_0110_FTNA_noDUF_volcanoplot_pless0.05.pdf", width = 3.5, height = 3)

# p<-ggplot(data=FTNAL.effect.annot, aes(x=diff.btw, y=-log10(wi.ep),color=diffexpr)) +
#   ggtitle("NA vs. FT (Light)")+
#   geom_point(size = 0.5) +
#   ylim(0,8) + xlim(-8,8) +
#   theme_minimal() +
#   theme(legend.position = "none") +
#   scale_color_manual(values=c("blue", "black", "red")) +
#geom_vline(xintercept=c(-3.5, 3.5), linetype = "dashed", col="grey36") +
#geom_hline(yintercept=-log10(0.05), linetype = "dashed",col="grey36")

#ggsave("aldex/SFR23_0110_FTNAL_noDUF_volcanoplot_pless0.05.pdf", width = 3.5, height = 3)

# p<-ggplot(data=FTNAD.effect.annot, aes(x=diff.btw, y=-log10(wi.ep),color=diffexpr)) +
#   ggtitle("NA vs. FT (Dark)")+
#   geom_point(size = 0.5) +
#   ylim(0,8) + xlim(-8,8) +
#   theme_minimal() +
#   theme(legend.position = "none") +
#   scale_color_manual(values=c("blue", "black", "red")) +
#geom_vline(xintercept=c(-3.5, 3.5), linetype = "dashed", col="grey36") +
#geom_hline(yintercept=-log10(0.05), linetype = "dashed",col="grey36")

#ggsave("aldex/SFR23_0110_FTNAD_noDUF_volcanoplot_pless0.05.pdf", width = 3.5, height = 3)

# list_venn <- list(light = (FTNAL.effect.annot%>%filter(wi.ep<0.05))$FeatureID,
#                   dark = (FTNAD.effect.annot%>%filter(wi.ep<0.05))$FeatureID,
#                   all = (FTNA.effect.annot%>%filter(wi.ep<0.05))$FeatureID)
# 
# ggvenn(list_venn, c("light", "dark", "all"),
#        fill_color = c("white", "black", "gray"))
# ggsave("aldex/SFR22_0824_FTNALDvenn_pless0.05.pdf", width = 3.5, height = 3)

# p<-ggplot(data=FANA.effect.annot, aes(x=diff.btw, y=-log10(wi.ep),color=diffexpr)) +
#   ggtitle("NA vs. FA")+
#   geom_point(size = 0.5) +
#   ylim(0,8) + xlim(-8,8) +
#   theme_minimal() +
#   theme(legend.position = "none") +
#   scale_color_manual(values=c("blue", "black", "red")) +
#geom_vline(xintercept=c(-3.5, 3.5), linetype = "dashed", col="grey36") +
#geom_hline(yintercept=-log10(0.05), linetype = "dashed",col="grey36")

#ggsave("aldex/SFR23_0110_FANA_noDUF_volcanoplot_pless0.05.pdf", width = 3.5, height = 3)

# p<-ggplot(data=FANAL.effect.annot, aes(x=diff.btw, y=-log10(wi.ep),color=diffexpr)) +
#   ggtitle("NA vs. FA (Light)")+
#   geom_point(size = 0.5) +  #ylim(0,8) + 
#   xlim(-8,8) +
#   theme_minimal() +
#   theme(legend.position = "none") +
#   scale_color_manual(values=c("blue", "black", "red")) +
#geom_vline(xintercept=c(-3.5, 3.5), linetype = "dashed", col="grey36") +
#geom_hline(yintercept=-log10(0.05), linetype = "dashed",col="grey36")

#ggsave("aldex/SFR23_0110_FANAL_noDUF_volcanoplot_pless0.05.pdf", width = 3.5, height = 3)

# p<-ggplot(data=FANAD.effect.annot, aes(x=diff.btw, y=-log10(wi.ep),color=diffexpr)) +
#   ggtitle("NA vs. FA (Dark)")+
#   geom_point(size = 0.5) + #ylim(0,8) 
#   xlim(-8,8) +
#   theme_minimal() +
#   theme(legend.position = "none") +
#   scale_color_manual(values=c("blue", "black", "red")) +
#geom_vline(xintercept=c(-3.5, 3.5), linetype = "dashed", col="grey36") +
#geom_hline(yintercept=-log10(0.05), linetype = "dashed",col="grey36")

#ggsave("aldex/SFR23_0110_FANAD_noDUF_volcanoplot_pless0.05.pdf", width = 3.5, height = 3)

#Summarise hits for all three pairwise comparisons (all-not just unique)

# list_venn <- list(light = (FANAL.effect.annot%>%filter(wi.ep<0.05))$FeatureID,
#                   dark = (FANAD.effect.annot%>%filter(wi.ep<0.05))$FeatureID,
#                   all = (FANA.effect.annot%>%filter(wi.ep<0.05))$FeatureID)
# 
# ggvenn(list_venn, c("light", "dark", "all"),
#        fill_color = c("white", "black", "gray"))
# ggsave("aldex/SFR22_0824_FANALDvenn_pless0.05.pdf", width = 3.5, height = 3)

# ggvenn(list_venn, c("light", "dark"),
#        fill_color = c("white", "black"))
# ggsave("aldex/SFR23_0111_FAFT_LDintersect_pless0.05.pdf", width = 3.5, height = 3)

# list_venn <- list(light = unique_hits_light$FeatureID,
#                   dark = unique_hits_dark$FeatureID,
#                   all = unique_hits_all$FeatureID)
# 
# ggvenn(list_venn, c("light", "dark", "all"),
#        fill_color = c("white", "black", "gray"))
# ggsave("aldex/SFR22_0829_FAFT_LDallintersect_pless0.05.pdf", width = 3.5, height = 3)
