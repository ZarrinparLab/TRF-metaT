setwd("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_m_results/")

library(tidyverse)
library(data.table)
library("qiime2R")
library("Biostrings")
library(ggrepel)
library(ggpubr)
library(gplots)
library(ggvenn)
library(viridis)
library(pheatmap)
library(ggpubfigs)
library(rstatix)
library(UpSetR)
library(ComplexHeatmap)
library(RColorBrewer)

library(nlme)
library(multcomp)
library(emmeans)
library(EnhancedVolcano)
###########################################################

metaTf_lme <- function(mtb) {
  tryCatch({
    m1 <- lme(log_10_peakabun ~ ENB * Phase,
              random = ~1 | mouseid,
              data = mtb)
    aov <- anova(m1)
    ph <- emmeans(m1, pairwise ~ ENB | Phase)
    ph_est <- as.data.frame(ph$contrasts)$estimate
    ph_pval <- as.data.frame(ph$contrasts)$p.value
    df <- data.frame(Intercept_pval = aov$`p-value`[[1]],
                     ENB_pval = aov$`p-value`[[2]],
                     Phase_pval = aov$`p-value`[[3]],
                     ENB_Phase_pval = aov$`p-value`[[4]],
                     light_AZ51vAZ52=ph_est[[1]],
                     light_AZ51vDnBSH1=ph_est[[2]],
                     light_AZ51vLgBSH=ph_est[[3]],
                     light_AZ52vDnBSH1=ph_est[[4]],
                     light_AZ52vLgBSH=ph_est[[5]],
                     light_DnBSH1vLgBSH=ph_est[[6]],
                     dark_AZ51vAZ52=ph_est[[7]],
                     dark_AZ51vDnBSH1=ph_est[[8]],
                     dark_AZ51vLgBSH=ph_est[[9]],
                     dark_AZ52vDnBSH1=ph_est[[10]],
                     dark_AZ52vLgBSH=ph_est[[11]],
                     dark_DnBSH1vLgBSH=ph_est[[12]],
                     light_AZ51vAZ52_pval=ph_pval[[1]],
                     light_AZ51vDnBSH1_pval=ph_pval[[2]],
                     light_AZ51vLgBSH_pval=ph_pval[[3]],
                     light_AZ52vDnBSH1_pval=ph_pval[[4]],
                     light_AZ52vLgBSH_pval=ph_pval[[5]],
                     light_DnBSH1vLgBSH_pval=ph_pval[[6]],
                     dark_AZ51vAZ52_pval=ph_pval[[7]],
                     dark_AZ51vDnBSH1_pval=ph_pval[[8]],
                     dark_AZ51vLgBSH_pval=ph_pval[[9]],
                     dark_AZ52vDnBSH1_pval=ph_pval[[10]],
                     dark_AZ52vLgBSH_pval=ph_pval[[11]],
                     dark_DnBSH1vLgBSH_pval=ph_pval[[12]])
    return(df)
  }, error = function(e) {
    # Handle the case of overfitting (or any other error)
    warning("Model fitting failed: ", conditionMessage(e))
    return(data.frame(Intercept_pval = NA,ENB_pval = NA,Phase_pval = NA,ENB_Phase_pval = NA,light_AZ51vAZ52=NA,
                      light_AZ51vDnBSH1=NA,light_AZ51vLgBSH=NA,light_AZ52vDnBSH1=NA,light_AZ52vLgBSH=NA,
                      light_DnBSH1vLgBSH=NA,dark_AZ51vAZ52=NA,dark_AZ51vDnBSH1=NA,dark_AZ51vLgBSH=NA,
                      dark_AZ52vDnBSH1=NA,dark_AZ52vLgBSH=NA,dark_DnBSH1vLgBSH=NA,light_AZ51vAZ52_pval=NA,
                      light_AZ51vDnBSH1_pval=NA,light_AZ51vLgBSH_pval=NA,light_AZ52vDnBSH1_pval=NA,light_AZ52vLgBSH_pval=NA,
                      light_DnBSH1vLgBSH_pval=NA,dark_AZ51vAZ52_pval=NA,dark_AZ51vDnBSH1_pval=NA,dark_AZ51vLgBSH_pval=NA,
                      dark_AZ52vDnBSH1_pval=NA,dark_AZ52vLgBSH_pval=NA,dark_DnBSH1vLgBSH_pval=NA))
  })
}

metaTf_lmecomb <- function(mtb) {
  tryCatch({
    m1 <- lme(log_10_peakabun ~ ENB,
              random = ~1 | mouseid,
              data = mtb)
    aov <- anova(m1)
    ph <- emmeans(m1, pairwise ~ ENB)
    ph_est <- as.data.frame(ph$contrasts)$estimate
    ph_pval <- as.data.frame(ph$contrasts)$p.value
    df <- data.frame(Intercept_pval = aov$`p-value`[[1]],
                     ENB_pval = aov$`p-value`[[2]],
                     AZ51vAZ52=ph_est[[1]],
                     AZ51vDnBSH1=ph_est[[2]],
                     AZ51vLgBSH=ph_est[[3]],
                     AZ52vDnBSH1=ph_est[[4]],
                     AZ52vLgBSH=ph_est[[5]],
                     DnBSH1vLgBSH=ph_est[[6]],
                     AZ51vAZ52_pval=ph_pval[[1]],
                     AZ51vDnBSH1_pval=ph_pval[[2]],
                     AZ51vLgBSH_pval=ph_pval[[3]],
                     AZ52vDnBSH1_pval=ph_pval[[4]],
                     AZ52vLgBSH_pval=ph_pval[[5]],
                     DnBSH1vLgBSH_pval=ph_pval[[6]])
    return(df)
  }, error = function(e) {
    # Handle the case of overfitting (or any other error)
    warning("Model fitting failed: ", conditionMessage(e))
    return(data.frame(Intercept_pval = NA,ENB_pval = NA,AZ51vAZ52=NA,AZ51vDnBSH1=NA,AZ51vLgBSH=NA,AZ52vDnBSH1=NA,AZ52vLgBSH=NA,DnBSH1vLgBSH=NA,
                      AZ51vAZ52_pval=NA,AZ51vDnBSH1_pval=NA,AZ51vLgBSH_pval=NA,AZ52vDnBSH1_pval=NA,AZ52vLgBSH_pval=NA,DnBSH1vLgBSH_pval=NA))
  })
}
###########################################################

#annot<-fread("BSH/BSH_ENB/GNPS_fecal_mtb/nf_output/library/merged_results_with_gnps.tsv")%>%
annot<-fread("BSH/BSH_ENB/GNPS_fecal_mtb/nf_output/library/merged_results_with_gnps_collapse.tsv")%>%
  mutate(Compound_Name = gsub("\\|.*", "", Compound_Name))%>%
  #mutate(Compound_Name = gsub("\\(delta mass.*", "", Compound_Name))%>%
  # mutate(conjugation=sub("^\\|", "", conjugation))%>%
  # mutate(conjugation=sub("(^[^|]+).*", "\\1", conjugation))%>%
  mutate(FeatureID=as.character(FeatureID))

#mtb<-fread("BSH/BSH_ENB/GNPS_fecal_mtb/nf_output/clustering/featuretable_reformated_justBA_summcollapse.csv")
mtb<-fread("BSH/BSH_ENB/GNPS_fecal_mtb/nf_output/clustering/featuretable_reformated_summcollapse.csv")

md<-fread("BSH/BSH_ENB/GNPS_fecal_mtb/metadata_filename/BSH_ENB_invivo_metadata.txt")
###########################################################
#get the standards
mtb_std<-mtb%>%
  dplyr::select("row.ID",contains("std") | contains("aminoacids"))%>%
  dplyr::rename(FeatureID=`row.ID`)%>%
  left_join(.,annot,by="FeatureID")%>%
  dplyr::select(c(1,37,2:22))

md_sub<-md%>%
  filter(collection_time!="")

mtb_sub<-mtb%>%
  dplyr::rename(FeatureID=`row.ID`)%>%
  dplyr::select(FeatureID, all_of(md_sub$filename))%>%
  gather(filename,peak_abun,-FeatureID)%>%
  left_join(.,annot,by="FeatureID")%>%
  left_join(.,md,by="filename")%>%
  mutate(log_10_peakabun=log10(peak_abun+1),
         Phase=factor(Phase,levels=c("Light","Dark")))

mtb_subL<-mtb_sub%>%
  filter(Phase=="Light")

mtb_subD<-mtb_sub%>%
  filter(Phase=="Dark")

TRF_unpr_ttest<-function(mtb){
  x<-pairwise.t.test(mtb$log_10_peakabun,mtb$ENB,p.adjust.method = "fdr")
  df<-data.frame(AZ51vAZ52h=x$p.value[1],
                 AZ51vDnBSH1=x$p.value[2],
                 AZ51vLgBSH=x$p.value[3],
                 AZ52vDnBSH1=x$p.value[5],
                 AZ52vLgBSH=x$p.value[6],
                 DnBSH1vLgBSH=x$p.value[9])
  return(df)
}

#stat.test<-mtb_sub%>%
#stat.testL<-mtb_subL%>%
stat.testD<-mtb_subD%>%
  group_by(FeatureID,Compound_Name)%>%
  nest()%>%
  mutate(pvals = map(data, ~ TRF_unpr_ttest(.x)))%>%
  dplyr::select(-data)%>%
  unnest()
write.table(stat.test,"BSH/BSH_ENB/GNPS_fecal_mtb/results_allmtb/ttest_unpr_ttest_enb_pvals.txt",sep = "\t",row.names = FALSE, quote=FALSE)
write.table(stat.test,"BSH/BSH_ENB/GNPS_fecal_mtb/results_allmtb/ttest_unpr_ttest_enbL_pvals.txt",sep = "\t",row.names = FALSE, quote=FALSE)
write.table(stat.test,"BSH/BSH_ENB/GNPS_fecal_mtb/results_allmtb/ttest_unpr_ttest_enbD_pvals.txt",sep = "\t",row.names = FALSE, quote=FALSE)

ttest_res<-fread("BSH/BSH_ENB/GNPS_fecal_mtb/results_allmtb/ttest_unpr_ttest_enb_pvals.txt")
dnbsh1sig<-ttest_res%>%filter(AZ51vDnBSH1<0.05 | AZ51vAZ52h<0.05| AZ51vLgBSH<0.05)%>%
  mutate(phase="combined")%>%
  filter(!is.na(Compound_Name))
ttest_resL<-fread("BSH/BSH_ENB/GNPS_fecal_mtb/results_allmtb/ttest_unpr_ttest_enbL_pvals.txt")
dnbsh1sigL<-ttest_resL%>%filter(AZ51vDnBSH1<0.05| AZ51vAZ52h<0.05| AZ51vLgBSH<0.05)%>%
  mutate(phase="light")%>%
  filter(!is.na(Compound_Name))
ttest_resD<-fread("BSH/BSH_ENB/GNPS_fecal_mtb/results_allmtb/ttest_unpr_ttest_enbD_pvals.txt")
dnbsh1sigD<-ttest_resD%>%filter(AZ51vDnBSH1<0.05| AZ51vAZ52h<0.05| AZ51vLgBSH<0.05)%>%
  mutate(phase="dark")%>%
  filter(!is.na(Compound_Name))

mtb_int<-rbind(dnbsh1sig,dnbsh1sigL,dnbsh1sigD)
write.table(mtb_int,"BSH/BSH_ENB/GNPS_fecal_mtb/results_allmtb/ttest_unpr_ttest_enb_pvals_pless0.05_combLD.txt",sep = "\t",row.names = FALSE, quote=FALSE)

mtb_intLD<-mtb_int%>%
  filter(phase!="combined")

mtb_int<-mtb_sub%>%
  filter(FeatureID %in% dnbsh1sig$FeatureID)%>%
  mutate(name=paste("mtb",FeatureID,Compound_Name,sep=" "))%>%
  group_by(FeatureID,name,ENB)%>%
  summarise(mn_peak_abun=mean(peak_abun))%>%
  mutate(log_10_peakabun=log10(mn_peak_abun+1))%>%
  group_by(FeatureID)%>%mutate(Zscore=(log_10_peakabun - mean(log_10_peakabun))/sd(log_10_peakabun))

mtb_int<-mtb_sub%>%
  filter(FeatureID %in% mtb_intLD$FeatureID)%>%
  filter(!is.na(Compound_Name))%>%
  mutate(name=paste("mtb",FeatureID,Compound_Name,sep=" "))%>%
  group_by(FeatureID,name,ENB,collection_time)%>%
  summarise(mn_peak_abun=mean(peak_abun))%>%
  mutate(log_10_peakabun=log10(mn_peak_abun+1))%>%
  group_by(FeatureID,collection_time)%>%mutate(Zscore=(log_10_peakabun - mean(log_10_peakabun))/sd(log_10_peakabun))%>%
  mutate(collection_time=factor(collection_time, levels=c("ZT3","ZT15")))

plt<-ggplot(mtb_int,aes(x=name, y=fct_rev(ENB))) +theme_classic()+
  #geom_tile(aes(fill=log_10_peakabun))+
  geom_tile(aes(fill=Zscore))+
  scale_x_discrete(expand = c(0, 0))+
  facet_grid(collection_time~.,scales="free",space="free")+
  theme(axis.ticks.y=element_blank(),panel.spacing.x=unit(0.3, "lines"),axis.text.y = element_text(size = 8),
        panel.spacing.y=unit(0.3, "lines"),axis.text.x = element_text(angle=90,hjust = 1))+
  scale_fill_distiller(palette = "Spectral", direction = -1)+
  labs(x="Bile acids",y="BSH strains")

ggsave("BSH/BSH_ENB/GNPS_fecal_mtb/results_allmtb/SFR25_0310_sighitsENBLDphase_heatmap_zscore.pdf", plt,width = 6.5, height = 16)

#plot these hits individually
indv_plot<-function(mtb,mtb_id){
  featmtb<-mtb%>%filter(FeatureID==mtb_id)
  p <- ggplot(featmtb, aes(x=ENB,y=log_10_peakabun,colour=ENB,fill=ENB)) +
    geom_boxplot(alpha=0.3)+
    geom_point()+
    facet_wrap(~Phase)+
    scale_fill_manual(values = c("gray70","#CC79A7","#009E73","#D55E00"))+
    scale_colour_manual(values = c("gray70","#CC79A7","#009E73","#D55E00"))+
    theme_pubr()+labs(title=str_wrap(paste(unique(featmtb$Compound_Name), " (", mtb_id, ")", sep = ""),width = 40))+
    theme(legend.position = "none")
  return(p)
}


for(i in dnbsh1sigD$FeatureID){
  p<-indv_plot(mtb_sub,i)
  #ggsave(paste("BSH/BSH_ENB/GNPS_fecal_mtb/results_allmtb/indv_plt/sigDnBSH10.05/",i,".pdf",sep=""),p,width = 5, height = 3)
  #ggsave(paste("BSH/BSH_ENB/GNPS_fecal_mtb/results_allmtb/indv_plt/sigDnBSH1L0.05/",i,".pdf",sep=""),p,width = 5, height = 3)
  ggsave(paste("BSH/BSH_ENB/GNPS_fecal_mtb/results_allmtb/indv_plt/sigDnBSH1D0.05/",i,".pdf",sep=""),p,width = 5, height = 4)
}

stat.test_lme<-mtb_sub%>%
  group_by(FeatureID,Compound_Name)%>%
  nest()%>%
  mutate(pvals = map(data, ~ metaTf_lme(.x)))%>%
  dplyr::select(-data)%>%
  unnest()
write.table(stat.test_lme,"BSH/BSH_ENB/GNPS_fecal_mtb/results_allmtb/lme_enb_pvals.txt",sep = "\t",row.names = FALSE, quote=FALSE)

#pval vs. fold change

#light AZ51 vs. DnBSH1
summ_AZ51vDnBSH1<-fread("BSH/BSH_ENB/GNPS_fecal_mtb/results_allmtb/lme_enb_pvals.txt")%>%
  dplyr::select(c(1,2,8,20))%>%
  mutate(light_AZ51vDnBSH1=-1*light_AZ51vDnBSH1)
write.table(summ_AZ51vDnBSH1,"BSH/BSH_ENB/GNPS_fecal_mtb/results_allmtb/lme_enb_pvals_lightAZ51vDnBSH1_ONLY.txt",sep = "\t",row.names = FALSE,quote=FALSE)

summ_AZ51vDnBSH1<-fread("BSH/BSH_ENB/GNPS_fecal_mtb/lme_enb_pvals_lightAZ51vDnBSH1_ONLY.txt")

#Dark AZ51 vs. DnBSH1
summ_AZ51vDnBSH1<-fread("BSH/BSH_ENB/GNPS_fecal_mtb/results_allmtb/lme_enb_pvals.txt")%>%
  dplyr::select(c(1,2,14,26))%>%
  mutate(dark_AZ51vDnBSH1=-1*dark_AZ51vDnBSH1)
write.table(summ_AZ51vDnBSH1,"BSH/BSH_ENB/GNPS_fecal_mtb/results_allmtb/lme_enb_pvals_darkAZ51vDnBSH1_ONLY.txt",sep = "\t",row.names = FALSE,quote=FALSE)

#light AZ51 vs. AZ52
summ_AZ51vAZ52<-fread("BSH/BSH_ENB/GNPS_fecal_mtb/results_allmtb/lme_enb_pvals.txt")%>%
  dplyr::select(c(1,2,7,19))%>%
  mutate(light_AZ51vAZ52=-1*light_AZ51vAZ52)
write.table(summ_AZ51vAZ52,"BSH/BSH_ENB/GNPS_fecal_mtb/results_allmtb/lme_enb_pvals_lightAZ51vAZ52_ONLY.txt",sep = "\t",row.names = FALSE,quote=FALSE)

#dark AZ51 vs. AZ52
summ_AZ51vAZ52<-fread("BSH/BSH_ENB/GNPS_fecal_mtb/results_allmtb/lme_enb_pvals.txt")%>%
  dplyr::select(c(1,2,13,25))%>%
  mutate(dark_AZ51vAZ52=-1*dark_AZ51vAZ52)
write.table(summ_AZ51vAZ52,"BSH/BSH_ENB/GNPS_fecal_mtb/results_allmtb/lme_enb_pvals_darkAZ51vAZ52_ONLY.txt",sep = "\t",row.names = FALSE,quote=FALSE)

#light AZ51 vs. LgBSH
summ_AZ51vLgBSH<-fread("BSH/BSH_ENB/GNPS_fecal_mtb/results_allmtb/lme_enb_pvals.txt")%>%
  dplyr::select(c(1,2,9,21))%>%
  mutate(light_AZ51vLgBSH=-1*light_AZ51vLgBSH)
write.table(summ_AZ51vLgBSH,"BSH/BSH_ENB/GNPS_fecal_mtb/results_allmtb/lme_enb_pvals_lightAZ51vLgBSH_ONLY.txt",sep = "\t",row.names = FALSE,quote=FALSE)

#dark AZ51 vs. LgBSH
summ_AZ51vLgBSH<-fread("BSH/BSH_ENB/GNPS_fecal_mtb/results_allmtb/lme_enb_pvals.txt")%>%
  dplyr::select(c(1,2,15,27))%>%
  mutate(dark_AZ51vLgBSH=-1*dark_AZ51vLgBSH)
write.table(summ_AZ51vLgBSH,"BSH/BSH_ENB/GNPS_fecal_mtb/results_allmtb/lme_enb_pvals_darkAZ51vLgBSH_ONLY.txt",sep = "\t",row.names = FALSE,quote=FALSE)

p<-EnhancedVolcano(
  # summ_AZ51vLgBSH,
  # lab = summ_AZ51vLgBSH$Compound_Name,
  # summ_AZ51vAZ52,
  # lab = summ_AZ51vAZ52$Compound_Name,
  summ_AZ51vDnBSH1,
  lab = summ_AZ51vDnBSH1$Compound_Name,
  xlim=c(-2,2),
  ylim=c(0,4),
  title = 'Light AZ-51 vs. DnBSH1',
  # title = 'Light AZ-51 vs. AZ-52',
  # title = 'Light AZ-51 vs. LgBSH',
  subtitle=NA,
  caption=NA,
  x = "light_AZ51vDnBSH1",
  y = "light_AZ51vDnBSH1_pval",
  # x = "light_AZ51vAZ52",
  # y = "light_AZ51vAZ52_pval",
  # x = "light_AZ51vLgBSH",
  # y = "light_AZ51vLgBSH_pval",
  labSize = 1,
  pointSize = 2,
  colAlpha=1,
  col=c('gray37', 'gray37','blue','red'),
  pCutoff = 0.05,
  FCcutoff = 0,
  legendPosition = "none",
  border="full",
  gridlines.major=FALSE,
  gridlines.minor=FALSE
) +theme(text=element_text(size=7))

ggsave("BSH/BSH_ENB/GNPS_fecal_mtb/results_allmtb/light_AZ51vDnBSH1_volplot.pdf",height=5, width=5)
ggsave("BSH/BSH_ENB/GNPS_fecal_mtb/results_allmtb/dark_AZ51vDnBSH1_volplot.pdf",height=5, width=5)
ggsave("BSH/BSH_ENB/GNPS_fecal_mtb/results_allmtb/light_AZ51vAZ52_volplot.pdf",height=5, width=5)
ggsave("BSH/BSH_ENB/GNPS_fecal_mtb/results_allmtb/dark_AZ51vAZ52_volplot.pdf",height=5, width=5)
ggsave("BSH/BSH_ENB/GNPS_fecal_mtb/results_allmtb/light_AZ51vLgBSH_volplot.pdf",height=5, width=5)
ggsave("BSH/BSH_ENB/GNPS_fecal_mtb/results_allmtb/dark_AZ51vLgBSH_volplot.pdf",height=5, width=5)

#combining light and dark

stat.test_lme_comb<-mtb_sub%>%
  group_by(FeatureID,Compound_Name)%>%
  nest()%>%
  mutate(pvals = map(data, ~ metaTf_lmecomb(.x)))%>%
  dplyr::select(-data)%>%
  unnest()
write.table(stat.test_lme_comb,"BSH/BSH_ENB/GNPS_fecal_mtb/results_allmtb/lme_enb_comb_pvals.txt",sep = "\t",row.names = FALSE, quote=FALSE)

#pval vs. fold change

#AZ51 vs. DnBSH1
summ_AZ51vDnBSH1<-fread("BSH/BSH_ENB/GNPS_fecal_mtb/results_allmtb/lme_enb_comb_pvals.txt")%>%
  dplyr::select(c(1,2,6,12))%>%
  mutate(AZ51vDnBSH1=-1*AZ51vDnBSH1)
write.table(summ_AZ51vDnBSH1,"BSH/BSH_ENB/GNPS_fecal_mtb/results_allmtb/lme_enb_pvals_AZ51vDnBSH1_ONLY.txt",sep = "\t",row.names = FALSE,quote=FALSE)

#AZ51 vs. AZ52
summ_AZ51vAZ52<-fread("BSH/BSH_ENB/GNPS_fecal_mtb/results_allmtb/lme_enb_comb_pvals.txt")%>%
  dplyr::select(c(1,2,5,11))%>%
  mutate(AZ51vAZ52=-1*AZ51vAZ52)
write.table(summ_AZ51vAZ52,"BSH/BSH_ENB/GNPS_fecal_mtb/results_allmtb/lme_enb_pvals_AZ51vAZ52_ONLY.txt",sep = "\t",row.names = FALSE,quote=FALSE)

#AZ51 vs. LgBSH
summ_AZ51vLgBSH<-fread("BSH/BSH_ENB/GNPS_fecal_mtb/results_allmtb/lme_enb_comb_pvals.txt")%>%
  dplyr::select(c(1,2,7,13))%>%
  mutate(AZ51vLgBSH=-1*AZ51vLgBSH)
write.table(summ_AZ51vLgBSH,"BSH/BSH_ENB/GNPS_fecal_mtb/results_allmtb/lme_enb_pvals_AZ51vLgBSH_ONLY.txt",sep = "\t",row.names = FALSE,quote=FALSE)

p<-EnhancedVolcano(
  summ_AZ51vLgBSH,
  lab = summ_AZ51vLgBSH$Compound_Name,
  # summ_AZ51vAZ52,
  # lab = summ_AZ51vAZ52$Compound_Name,
  # summ_AZ51vDnBSH1,
  # lab = summ_AZ51vDnBSH1$Compound_Name,
  xlim=c(-2,2),
  ylim=c(0,4),
  # title = 'AZ-51 vs. DnBSH1',
  # title = 'AZ-51 vs. AZ-52',
  title = 'AZ-51 vs. LgBSH',
  subtitle=NA,
  caption=NA,
  # x = "AZ51vDnBSH1",
  # y = "AZ51vDnBSH1_pval",
  # x = "AZ51vAZ52",
  # y = "AZ51vAZ52_pval",
  x = "AZ51vLgBSH",
  y = "AZ51vLgBSH_pval",
  labSize = 1,
  pointSize = 2,
  colAlpha=1,
  col=c('gray37', 'gray37','blue','red'),
  pCutoff = 0.05,
  FCcutoff = 0,
  legendPosition = "none",
  border="full",
  gridlines.major=FALSE,
  gridlines.minor=FALSE
) +theme(text=element_text(size=7))

ggsave("BSH/BSH_ENB/GNPS_fecal_mtb/results_allmtb/AZ51vDnBSH1_volplot.pdf",height=5, width=5)
ggsave("BSH/BSH_ENB/GNPS_fecal_mtb/results_allmtb/AZ51vAZ52_volplot.pdf",height=5, width=5)
ggsave("BSH/BSH_ENB/GNPS_fecal_mtb/results_allmtb/AZ51vLgBSH_volplot.pdf",height=5, width=5)
