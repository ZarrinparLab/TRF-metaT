setwd("~/Notebooks/sfloresr/TRF-metaT/")

library(tidyverse)
library(data.table)
library(ALDEx2)
library(ggpubr)
library(scales)
library(ggvenn)

###########################################################
#inputs
dat_nn<-"data/pfam_metaT/pfam_clean_noNT.txt"
dat_metadata<-"data/metaT_metadata_ztcat_noNT.txt"
dat_go<-"data/pfam_metaT/go_name.txt"
dat_pfamtogo<-"data/pfam_metaT/pfam-to-go-process.map"
dat_annot<-"data/pfam_metaT/pfam_annotationkey.csv"
dat_mtx_nrpob<-"data/pfam_metaT/pfam-TPM_clean_noNT_normRPOB.txt"
dat_path<-"data/pfam_metaT/"
res_path<-"data/DE_analysis/aldex_metaT/SFR23_0620_"
fig_path<-"figures/DE_analysis/"
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

get_pfamtogo<-function(datpfamgo){
  pfamGOp<-read.table(datpfamgo,header = FALSE, sep = "\t",
                      col.names = paste0("V",seq_len(4)), fill = TRUE)%>%
    gather(column,GO_Term,-V1)%>%
    dplyr::select(1,3)%>%
    dplyr::rename(FeatureID=V1)%>%
    left_join(.,gonames,by="GO_Term")%>%
    filter(!is.na(name))
  return(pfamGOp)
}

ZT_dist <- function(data,id) {
  
  data_sub<-data%>%filter(FeatureID==id)
  p<-data_sub%>%
    ggplot(aes(x=condition, y=log10(TPM_counts+1),fill=condition))+
    geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                           position=position_dodge(1)) +
    theme_minimal()+ scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+
    scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
    theme_pubr() +
    facet_grid(~phase)+
    labs(color="condition",
         y ="log10(Transcript/RPOB TPM)",
         x ="ZT time",
         title=id)+ theme(plot.title = element_text(face = "bold"),
                                 legend.position = "right")
  return(p)
}
###########################################################
#create pairwise files

#FA vs FT
pfamFAFT<-fread(dat_nn)%>%dplyr::select(c(1:25))
write.table(pfamFAFT,paste0(dat_path,"pfam_FAFT.txt"),sep = "\t", row.names = FALSE, quote=FALSE)

#just light
pfamFAFTL<-fread(dat_nn)%>%dplyr::select(c(1:7,14:19))
write.table(pfamFAFTL,paste0(dat_path,"pfam_FAFTL.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

#just dark
pfamFAFTD<-fread(dat_nn)%>%dplyr::select(c(1,8:13,20:25))
write.table(pfamFAFTD,paste0(dat_path,"pfam_FAFTD.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

#FT vs NA
pfamFTNA<-fread(dat_nn)%>%dplyr::select(c(1,14:43))
write.table(pfamFTNA,paste0(dat_path,"pfam_FTNA.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

#light
pfamFTNAL<-fread(dat_nn)%>%dplyr::select(c(1,14:19,26:34))
write.table(pfamFTNAL,paste0(dat_path,"pfam_FTNAL.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

#dark
pfamFTNAD<-fread(dat_nn)%>%dplyr::select(c(1,20:25,35:43))
write.table(pfamFTNAD,paste0(dat_path,"pfam_FTNAD.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

#FA vs NA
pfamFANA<-fread(dat_nn)%>%dplyr::select(c(1:13,26:43))
write.table(pfamFANA,paste0(dat_path,"pfam_FANA.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

#light
pfamFANAL<-fread(dat_nn)%>%dplyr::select(c(1:7,26:34))
write.table(pfamFANAL,paste0(dat_path,"pfam_FANAL.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

#dark
pfamFANAD<-fread(dat_nn)%>%dplyr::select(c(1,8:13,35:43))
write.table(pfamFANAD,paste0(dat_path,"pfam_FANAD.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

###########################################################
#run aldex

annot<-fread(dat_annot)
md<-fread(dat_metadata)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))
conds <- md$condition

##FA vs. FT
pfamFAFT<-fread(paste0(dat_path,"pfam_FAFT.txt"))%>%column_to_rownames("FeatureID")
condsFAvsFT<-conds[1:24]
FAFT.effect.annot<-run_aldex(pfamFAFT,condsFAvsFT)
write.table(FAFT.effect.annot,paste0(res_path,"FAFT_ald_effectwpval_wannot.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

FAFT.effect.annot %>% group_by(diffexpr) %>% tally(sort = TRUE)

# diffexpr     n
# <chr>    <int>
#   1 none      5379
# 2 down         9
# 3 up           6

#FA vs FT (light)
md<-fread(dat_metadata)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition)) %>%
  filter(zt_time<13)
conds <- md$condition

condsFAvsFTL<-conds[1:12]
pfamFAFTL<-fread(paste0(dat_path,"pfam_FAFTL.txt"))%>%column_to_rownames("FeatureID")
FAFTL.effect.annot<-run_aldex(pfamFAFTL,condsFAvsFTL)
write.table(FAFTL.effect.annot,paste0(res_path,"FAFTL_ald_effectwpval_wannot.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

FAFTL.effect.annot %>% group_by(diffexpr) %>% tally(sort = TRUE)

#   diffexpr     n
#   <chr>    <int>
# 1 none      5118
  
#FA vs FT (dark)
md<-fread(dat_metadata)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition)) %>%
  filter(zt_time>9)
conds <- md$condition

condsFAvsFTD<-conds[1:12]
pfamFAFTD<-fread(paste0(dat_path,"pfam_FAFTD.txt"))%>%column_to_rownames("FeatureID")
FAFTD.effect.annot<-run_aldex(pfamFAFTD,condsFAvsFTD)
write.table(FAFTD.effect.annot,paste0(res_path,"FAFTD_ald_effectwpval_wannot.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

FAFTD.effect.annot %>% group_by(diffexpr) %>% tally(sort = TRUE)

# A tibble: 3 × 2
# diffexpr     n
# <chr>    <int>
#   1 none      4370
# 2 up         162
# 3 down       149

##NA vs FT
md<-fread(dat_metadata)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))
conds <- md$condition

pfamFTNA<-fread(paste0(dat_path,"pfam_FTNA.txt"))%>%column_to_rownames("FeatureID")
condsFTvsNA<-conds[c(13:42)]
FTNA.effect.annot<-run_aldex(pfamFTNA,condsFTvsNA)
write.table(FTNA.effect.annot,paste0(res_path,"FTNA_ald_effectwpval_wannot.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

FTNA.effect.annot %>% group_by(diffexpr) %>% tally(sort = TRUE)

# diffexpr     n
# <chr>    <int>
#   1 none      4739
# 2 up        1013
# 3 down       761

##NA vs FT (light)
md<-fread(dat_metadata)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition)) %>%
  filter(zt_time<13)
conds <- md$condition

condsFTvsNAL<-conds[c(7:21)]
pfamFTNAL<-fread(paste0(dat_path,"pfam_FTNAL.txt"))%>%column_to_rownames("FeatureID")
FTNAL.effect.annot<-run_aldex(pfamFTNAL,condsFTvsNAL)
write.table(FTNAL.effect.annot,paste0(res_path,"FTNAL_ald_effectwpval_wannot.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

FTNAL.effect.annot %>% group_by(diffexpr) %>% tally(sort = TRUE)

# diffexpr     n
# <chr>    <int>
#   1 none      4978
# 2 up         849
# 3 down       493

#NA vs FT (dark)
md<-fread(dat_metadata)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition)) %>%
  filter(zt_time>9)
conds <- md$condition

condsFTvsNAD<-conds[c(7:21)]
pfamFTNAD<-fread(paste0(dat_path,"pfam_FTNAD.txt"))%>%column_to_rownames("FeatureID")
FTNAD.effect.annot<-run_aldex(pfamFTNAD,condsFTvsNAD)
write.table(FTNAD.effect.annot,paste0(res_path,"FTNAD_ald_effectwpval_wannot.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

FTNAD.effect.annot %>% group_by(diffexpr) %>% tally(sort = TRUE)

# diffexpr     n
# <chr>    <int>
#   1 none      4407
# 2 up         657
# 3 down       428

##NA vs FA
md<-fread(dat_metadata)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))
conds <- md$condition

pfamFANA<-fread(paste0(dat_path,"pfam_FANA.txt"))%>%column_to_rownames("FeatureID")
condsFAvsNA<-conds[c(1:12,25:42)]
FANA.effect.annot<-run_aldex(pfamFANA,condsFAvsNA)
write.table(FANA.effect.annot,paste0(res_path,"FANA_ald_effectwpval_wannot.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

FANA.effect.annot %>% group_by(diffexpr) %>% tally(sort = TRUE)

# diffexpr     n
# <chr>    <int>
#   1 none      4704
# 2 up        1099
# 3 down       69

##NA vs FA (light)
md<-fread(dat_metadata)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition)) %>%
  filter(zt_time<13)
conds <- md$condition

condsFAvsNAL<-conds[c(1:6,13:21)]
pfamFANAL<-fread(paste0(dat_path,"pfam_FANAL.txt"))%>%column_to_rownames("FeatureID")
FANA.effect.annot<-run_aldex(pfamFANAL,condsFAvsNAL)
write.table(FANAL.effect.annot,paste0(res_path,"FANAL_ald_effectwpval_wannot.txt"),sep = "\t",row.names = FALSE, quote=FALSE)


FANAL.effect.annot %>% group_by(diffexpr) %>% tally(sort = TRUE)

# diffexpr     n
# <chr>    <int>
#   1 none      5194
# 2 up         748
# 3 down       371

#NA vs FA (dark)
md<-fread(dat_metadata)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition)) %>%
  filter(zt_time>9)
conds <- md$condition

condsFAvsNAD<-conds[c(1:6,13:21)]
pfamFANAD<-fread(paste0(dat_path,"pfam_FANAD.txt"))%>%column_to_rownames("FeatureID")
FANAD.effect.annot<-run_aldex(pfamFANAD,condsFAvsNAD)
write.table(FANAD.effect.annot,paste0(res_path,"FANAD_ald_effectwpval_wannot.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

FANAD.effect.annot %>% group_by(diffexpr) %>% tally(sort = TRUE)

# diffexpr     n
# <chr>    <int>
#   1 none      4318
# 2 up         676
# 3 down       457
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

###########################################################
#unique hits per pairwise comparison

list_venn <- list(FAFT = (FAFTL.effect.annot%>%filter(wi.eBH<0.1))$FeatureID,
                  FANA = (FANAL.effect.annot%>%filter(wi.eBH<0.1))$FeatureID,
                  FTNA = (FTNAL.effect.annot%>%filter(wi.eBH<0.1))$FeatureID)

ItemsList <- venn(list_venn, show.plot = FALSE)
light<-attributes(ItemsList)$intersections

unique_hits_light<-FAFTL.effect.annot%>%filter(FeatureID %in% light$FAFT)
unique_hits_lightFANA<-FANAL.effect.annot%>%filter(FeatureID %in% light$FANA)
unique_hits_lightFTNA<-FTNAL.effect.annot%>%filter(FeatureID %in% light$FTNA)

list_venn <- list(FAFT = (FAFTD.effect.annot%>%filter(wi.eBH<0.1))$FeatureID,
                  FANA = (FANAD.effect.annot%>%filter(wi.eBH<0.1))$FeatureID,
                  FTNA = (FTNAD.effect.annot%>%filter(wi.eBH<0.1))$FeatureID)

ItemsList <- venn(list_venn, show.plot = FALSE)
dark<-attributes(ItemsList)$intersections

unique_hits_dark<-FAFTD.effect.annot%>%filter(FeatureID %in% dark$FAFT)
unique_hits_darkFANA<-FANAD.effect.annot%>%filter(FeatureID %in% dark$FANA)
unique_hits_darkFTNA<-FTNAD.effect.annot%>%filter(FeatureID %in% dark$FTNA)
###########################################################
#create GO summary plot for NAFA--Figure 2D

gonames<-fread(dat_go)

pfamFANAD<-get_pfamtogo(dat_pfamtogo)%>%
  filter(FeatureID %in% unique_hits_darkFANA$FeatureID)%>%
  mutate(phase="dark")

pfamFANAL<-get_pfamtogo(dat_pfamtogo)%>%
  filter(FeatureID %in% unique_hits_lightFANA$FeatureID)%>%
  mutate(phase="light")

pfamFANA<-rbind(pfamFANAD,pfamFANAL)

pfamFANA_summ<-pfamFANA%>%
  group_by(name,phase)%>%summarise(n=n())%>%
  group_by(name)%>%mutate(sum=sum(n))%>%
  mutate(n=ifelse(phase=="light", n*-1,n),
         phase=factor(phase, levels = c("light","dark")))%>%
  mutate(n=as.integer(n))%>%
  arrange(sum)%>%
  filter(sum>1)

pfamFANA_summ$name <- factor(pfamFANA_summ$name,levels = unique(pfamFANA_summ$name))

ggplot(data=pfamFANA_summ, aes(x=name, y=n, fill=phase)) +
  geom_bar(stat="identity", position="identity") +
  scale_fill_manual(values=c("gray70","gray10")) +
  scale_y_continuous(expand=c(0,0),breaks= pretty_breaks())+
  coord_flip() +theme_pubr() +
  ggtitle("FANA")

ggsave(paste0(fig_path,"FANA_GOterms.pdf"),height=6.5, width=8)
###########################################################
#create GO summary plot for NAFA--Figure S2C

pfamFTNAD<-get_pfamtogo(dat_pfamtogo)%>%
  filter(FeatureID %in% unique_hits_darkFTNA$FeatureID)%>%
  mutate(phase="dark")

pfamFTNAL<-get_pfamtogo(dat_pfamtogo)%>%
  filter(FeatureID %in% unique_hits_lightFTNA$FeatureID)%>%
  mutate(phase="light")

pfamFTNA<-rbind(pfamFTNAD,pfamFTNAL)

pfamFTNA_summ<-pfamFTNA%>%
  group_by(name,phase)%>%summarise(n=n())%>%
  group_by(name)%>%mutate(sum=sum(n))%>%
  mutate(n=ifelse(phase=="light", n*-1,n),
         phase=factor(phase, levels = c("light","dark")))%>%
  arrange(sum)%>%
  filter(sum>1)

pfamFTNA_summ$name <- factor(pfamFTNA_summ$name,levels = unique(pfamFTNA_summ$name))

ggplot(data=pfamFTNA_summ, aes(x=name, y=n, fill=phase)) +
  geom_bar(stat="identity", position="identity") +
  scale_fill_manual(values=c("gray70","gray10")) +
  coord_flip() +theme_pubr() +
  ggtitle("FTNA")+
  scale_y_continuous(expand=c(0,0))

ggsave(paste0(fig_path,"FTNA_GOterms.pdf"),height=8, width=8)
###########################################################
#create GO summary plot for FAFT--Figure S2D

pfamFAFTD<-get_pfamtogo(dat_pfamtogo)%>%
  filter(FeatureID %in% unique_hits_dark$FeatureID)%>%
  mutate(phase="dark")

pfamFAFTL<-get_pfamtogo(dat_pfamtogo)%>%
  filter(FeatureID %in% unique_hits_light$FeatureID)%>%
  mutate(phase="light")

pfamFAFT<-rbind(pfamFAFTD,pfamFAFTL)

pfamFAFT_summ<-pfamFAFT%>%
  group_by(name,phase)%>%summarise(n=n())%>%
  group_by(name)%>%mutate(sum=sum(n))%>%
  mutate(n=ifelse(phase=="light", n*-1,n),
         phase=factor(phase, levels = c("light","dark")))%>%
  arrange(sum)

pfamFAFT_summ$name <- factor(pfamFAFT_summ$name,levels = unique(pfamFAFT_summ$name))

ggplot(data=pfamFAFT_summ, aes(x=name, y=n, fill=phase)) +
  geom_bar(stat="identity", position="identity") +
  scale_fill_manual(values=c("gray10")) +coord_flip() +theme_pubr() +
  ggtitle("FAFT")+
  scale_y_continuous(expand=c(0,0))

ggsave(paste0(fig_path,"FAFT_GOterms_wsingleton.pdf"),height=6, width=7)

##################################################################
#plot specific examples but light and dark phase--Figure 2E-F, S2E-F

pfamZT<-fread(dat_mtx_nrpob)%>%
  gather(sample_name,TPM_counts,-FeatureID) %>%
  left_join(.,annot, by ="FeatureID") %>%
  mutate(label_name=paste(FeatureID, Name, sep=" "))%>%
  left_join(.,md,by="sample_name") %>%
  mutate(condition=factor(condition, levels = c("NA","FA","FT")),
         phase=factor(phase, levels = c("light","dark")))

#functions that cycled in both FA and FT--Fig S3D
p<-ZT_dist(pfamZT,"PF01219.22")
ggsave(paste0(fig_path,"FAFTD_zt/PF01219.22.pdf"), p, height=3, width=5) #Fig 2E
p<-ZT_dist(pfamZT,"PF01704.21")
ggsave(paste0(fig_path,"FAFTD_zt/PF01704.21.pdf"), p, height=3, width=5) #Fig 2F

p<-ZT_dist(pfamZT,"PF13561.9")
ggsave(paste0(fig_path,"FAFTD_zt/PF13561.9.pdf"), p, height=3, width=5) #Fig S2E
p<-ZT_dist(pfamZT,"PF03705.18")
ggsave(paste0(fig_path,"FAFTD_zt/PF03705.18.pdf"), p, height=3, width=5) #Fig S2F