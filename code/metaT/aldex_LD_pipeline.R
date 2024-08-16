setwd("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_m_results/")

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
library(scales)
###########################################################
#inputs

pfam_mtb<-"pfam/pfam_clean_noNT.tsv"
md_file<-"metaT_metadata_ztcat_noNT.txt"
dir<-"pfam"

###########################################################
#files for light vs dark comparison

#FA
pfamFA<-fread(pfam_mtb)%>%dplyr::select(c(1:13))
write.table(pfamFA,paste0(dir,"/pfam_FA.txt"),sep = "\t", row.names = FALSE, quote=FALSE)

#FT vs NA
pfamFT<-fread(pfam_mtb)%>%dplyr::select(c(1,14:25))
write.table(pfamFT,paste0(dir,"/pfam_FT.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

#FA vs NA
pfamNA<-fread(pfam_mtb)%>%dplyr::select(c(1,26:43))
write.table(pfamNA,paste0(dir,"/pfam_NA.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

###########################################################
#aldex
annot<-fread(paste0(dir,"/pfam_annotationkey.csv"))
md<-fread(md_file)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))
phase <- md$phase

dir2<-"pfam/aldex/SFR23_0623_"

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
#   1 none      4878

##FT light dark
md<-fread(md_file)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))
phase <- md$phase

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

# A tibble: 1 × 2
# diffexpr     n
# <chr>    <int>
#   1 none      4958

##NA light vs. dark
md<-fread(md_file)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))
phase <- md$phase

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
#   1 none      6127
# 2 up         211
# 3 down        46

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

###########################################################
#GO summary plot
gonames<-fread("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_results/go_name.txt")
hitsL<-NA.effect.annot%>%filter(wi.eBH<0.1 & diff.btw< -0)
hitsD<-NA.effect.annot%>%filter(wi.eBH<0.1 & diff.btw> 0)

#NA light vs. dark
pfamNAL<-read.table("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_results/pfam-to-go-process.map",header = FALSE, sep = "\t",
                      col.names = paste0("V",seq_len(4)), fill = TRUE)%>%
  gather(column,GO_Term,-V1)%>%
  dplyr::select(1,3)%>%
  dplyr::rename(FeatureID=V1)%>%
  filter(FeatureID %in% hitsL$FeatureID)%>%
  left_join(.,gonames,by="GO_Term")%>%
  filter(!is.na(name))%>%
  mutate(phase="light")

pfamNAD<-read.table("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_results/pfam-to-go-process.map",header = FALSE, sep = "\t",
                    col.names = paste0("V",seq_len(4)), fill = TRUE)%>%
  gather(column,GO_Term,-V1)%>%
  dplyr::select(1,3)%>%
  dplyr::rename(FeatureID=V1)%>%
  filter(FeatureID %in% hitsD$FeatureID)%>%
  left_join(.,gonames,by="GO_Term")%>%
  filter(!is.na(name))%>%
  mutate(phase="dark")

pfamNA<-rbind(pfamNAL,pfamNAD)

pfamNA_summ<-pfamNA%>%
  group_by(name,phase)%>%summarise(n=n())%>%
  group_by(name)%>%mutate(sum=sum(n))%>%
  mutate(n=ifelse(phase=="light", n*-1,n),
         phase=factor(phase, levels = c("light","dark")))%>%
  mutate(n=as.integer(n))%>%
  arrange(sum)%>%
  filter(sum>1)

pfamNA_summ$name <- factor(pfamNA_summ$name,levels = unique(pfamNA_summ$name))

ggplot(data=pfamNA_summ, aes(x=name, y=n, fill=phase)) +
  geom_bar(stat="identity", position="identity") +
  scale_fill_manual(values=c("gray70","gray10")) +
  scale_y_continuous(expand=c(0,0),breaks= pretty_breaks())+
  coord_flip() +theme_pubr() +
  ggtitle("NA (light vs. dark)")

ggsave(paste0(dir2,"NALD_GOterms.pdf"),height=3.4, width=6)

###########################################################

#plot the individual hits
pfamZT<-fread("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_m_results/pfam/pfam-TPM_clean_noNT_normRPOB.txt")%>%
  #filter(FeatureID %in% pfamNAL$FeatureID)%>% 
  filter(FeatureID %in% pfamNAD$FeatureID)%>% 
  gather(sample_name,TPM_counts,-FeatureID) %>%
  left_join(.,annot, by ="FeatureID") %>%
  mutate(label_name=paste(FeatureID, Name, sep=" "))%>%
  left_join(.,md,by="sample_name") %>%
  mutate(condition=factor(condition, levels = c("NA","FA","FT")),
         phase=factor(phase, levels = c("light","dark")))

ZT_dist <- function(data,FeatureID) {
  
  data%>%
    ggplot(aes(x=condition, y=log10(TPM_counts+1),fill=interaction(condition,phase), dodge=condition))+
    geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                             position=position_dodge(1)) +
    theme_minimal()+ scale_color_manual(values=c("#0072B2","#D55E00","#009E73","#004166","#893c00","#00523b"))+
    scale_fill_manual(values=c("#0072B2","#D55E00","#009E73","#004166","#893c00","#00523b"))+
    theme_pubr() +
    labs(color="condition",
         y ="log10(Transcript/RPOB TPM)",
         x ="ZT time",
         title=FeatureID)+ theme(plot.title = element_text(face = "bold"),
                                 legend.position = "right")

  ggsave(paste("pfam/aldex/NA_dark_zt/",FeatureID,".pdf",sep=""), height=3, width=5)
}


pfamZT_nested <- pfamZT %>% 
  group_by(FeatureID) %>% 
  nest()

pfamZT_plots <- 
  pfamZT_nested %>% 
  mutate(plot = map2(data, FeatureID,  ~ ZT_dist(.x,.y)))
