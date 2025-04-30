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

###########################################################

annot<-fread("BSH/BSH_ENB/GNPS_fecal_mtb/nf_output/library/merged_results_with_gnps_justBA_curatedconjannot_collapse.tsv")%>%
  mutate(Compound_Name = gsub("\\|.*", "", Compound_Name))%>%
  #mutate(Compound_Name = gsub("\\(delta mass.*", "", Compound_Name))%>%
  mutate(conjugation=sub("^\\|", "", conjugation))%>%
  mutate(conjugation=sub("(^[^|]+).*", "\\1", conjugation))%>%
  mutate(FeatureID=as.character(FeatureID))

mtb<-fread("BSH/BSH_ENB/GNPS_fecal_mtb/nf_output/clustering/featuretable_reformated_justBA_summcollapse.csv")

md<-fread("BSH/BSH_ENB/GNPS_fecal_mtb/metadata_filename/BSH_ENB_invivo_metadata.txt")
###########################################################
#get the standards
mtb_std<-mtb%>%
  select("row.ID",contains("std") | contains("aminoacids"))%>%
  dplyr::rename(FeatureID=`row.ID`)%>%
  left_join(.,annot,by="FeatureID")%>%
  dplyr::select(c(1,37,2:22))

#heatmap of all hits w/replicates 
md_sub<-md%>%
  filter(collection_time!="")

mtb_sub<-mtb%>%
  dplyr::rename(FeatureID=`row.ID`)%>%
  dplyr::select(FeatureID, all_of(md_sub$filename))%>%
  gather(filename,peak_abun,-FeatureID)%>%
  left_join(.,annot,by="FeatureID")%>%
  left_join(.,md,by="filename")%>%
  group_by(mouseid,ENB,conjugation)%>%
  summarise(mn_peak_abun=mean(peak_abun))%>%
  mutate(log_10_peakabun=log10(mn_peak_abun+1),
         ENB=factor(ENB,levels=c("AZ-51","AZ-52","DnBSH1","LgBSH")))%>%
  #filter(!(conjugation %in% c("candidate")))%>%
  group_by(mouseid) %>%
  mutate(total_peak_abun = sum(mn_peak_abun, na.rm = TRUE),
    conjugation_pct = (mn_peak_abun / total_peak_abun) * 100)

ordered_mice <- mtb_sub %>%
  filter(conjugation == "unconjugated") %>%
  arrange(conjugation_pct) %>%
  pull(mouseid)

mtb_sub <- mtb_sub %>%
  mutate(mouseid_ordered = factor(mouseid, levels = unique(ordered_mice)))

plt <- ggplot(mtb_sub, aes(x = mouseid_ordered, y = conjugation_pct, fill = conjugation)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(~ENB, scales = "free")+
  scale_x_discrete(expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0))+
  theme_pubr()+theme(axis.text.x = element_blank())

ggsave("BSH/BSH_ENB/GNPS_fecal_mtb/SFR25_0121_summconjugation_bysample.pdf.pdf", plt,width = 10, height = 4)
ggsave("BSH/BSH_ENB/GNPS_fecal_mtb/SFR25_0121_summconjugation_bysamplev2.pdf", plt,width = 10, height = 4)
ggsave("BSH/BSH_ENB/GNPS_fecal_mtb/SFR25_0121_summconjugation_bymicev2.pdf", plt,width = 10, height = 4)
ggsave("BSH/BSH_ENB/GNPS_fecal_mtb/SFR25_0121_summconjugation_bymice.pdf", plt,width = 10, height = 4)

#show abun/prevalnece of BBAA
prev_bbaa<-mtb%>%
  dplyr::rename(FeatureID=`row.ID`)%>%
  dplyr::select(FeatureID, all_of(md_sub$filename))%>%
  gather(filename,peak_abun,-FeatureID)%>%
  left_join(.,annot,by="FeatureID")%>%
  left_join(.,md,by="filename")%>%
  filter(conjugation=="conjugated")%>%
  mutate(log_10_peakabun=log10(peak_abun+1),
         ENB=factor(ENB,levels=c("AZ-51","AZ-52","DnBSH1","LgBSH")))%>%
  group_by(FeatureID,Compound_Name,mouseid,ENB)%>%
  summarise(mn_peak_abun=mean(peak_abun))%>%
  filter(mn_peak_abun>0)%>%
  group_by(FeatureID,Compound_Name,ENB)%>%
  summarise(n=n())

mtb_inAZ51<-(prev_bbaa%>%filter(ENB=="AZ-51"& n>0))$FeatureID #all that are in AZ-51

bbaa_noAZ51<-mtb%>%
  dplyr::rename(FeatureID=`row.ID`)%>%
  dplyr::select(FeatureID, all_of(md_sub$filename))%>%
  gather(filename,peak_abun,-FeatureID)%>%
  left_join(.,annot,by="FeatureID")%>%
  left_join(.,md,by="filename")%>%
  filter(conjugation=="conjugated")%>%
  mutate(log_10_peakabun=log10(peak_abun+1),
         ENB=factor(ENB,levels=c("AZ-51","AZ-52","DnBSH1","LgBSH")))%>%
  group_by(mouseid,ENB,FeatureID)%>%
  summarise(mn_peak_abun=mean(peak_abun))%>%
  filter(!(FeatureID %in% mtb_inAZ51))%>%
  filter(mn_peak_abun>0)

unique_bbaa<-unique(bbaa_noAZ51$FeatureID) #feautures no in AZ51

bbaa_highestabun<-mtb%>%
  dplyr::rename(FeatureID=`row.ID`)%>%
  dplyr::select(FeatureID, all_of(md_sub$filename))%>%
  gather(filename,peak_abun,-FeatureID)%>%
  left_join(.,annot,by="FeatureID")%>%
  left_join(.,md,by="filename")%>%
  filter(conjugation=="conjugated")%>%
  mutate(log_10_peakabun=log10(peak_abun+1),
         ENB=factor(ENB,levels=c("AZ-51","AZ-52","DnBSH1","LgBSH")))%>%
  group_by(FeatureID)%>%
  summarise(mn_peak_abun=mean(peak_abun))%>%
  arrange(desc(mn_peak_abun))%>%
  top_n(7)

top_bbaa<-unique(bbaa_highestabun$FeatureID)

bbaa<-mtb%>%
  dplyr::rename(FeatureID=`row.ID`)%>%
  dplyr::select(FeatureID, all_of(md_sub$filename))%>%
  gather(filename,peak_abun,-FeatureID)%>%
  left_join(.,annot,by="FeatureID")%>%
  left_join(.,md,by="filename")%>%
  filter(conjugation=="conjugated")%>%
  #mutate(conjid=ifelse(FeatureID %in% unique_bbaa,paste(Compound_Name,FeatureID,sep=" "),"z_other"))%>%
  mutate(conjid=ifelse(FeatureID %in% top_bbaa,paste(Compound_Name,FeatureID,sep=" "),"z_other"))%>%
  group_by(conjid,mouseid,ENB)%>%
  summarise(mn_peak_abun=mean(peak_abun))%>%
  mutate(ENB=factor(ENB,levels=c("AZ-51","AZ-52","DnBSH1","LgBSH")))

plt <- ggplot(bbaa, aes(x = mouseid, y = mn_peak_abun, fill=conjid)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(~ENB, scales = "free")+
  scale_x_discrete(expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0))+
  theme_pubr()+theme(legend.position="right",axis.text.x = element_blank())
#ggsave("BSH/BSH_ENB/GNPS_fecal_mtb/SFR25_0127_bbaanotAZ51_bymice.pdf", plt,width = 20, height = 4)
ggsave("BSH/BSH_ENB/GNPS_fecal_mtb/SFR25_0127_bbaatopabun_bymice.pdf", plt,width = 20, height = 4)

# plt<-ggplot(mtb_sub,aes(x=FeatureID, y=ENB)) +theme_classic()+
#   geom_tile(aes(fill=log_10_peakabun))+
#   scale_x_discrete(expand = c(0, 0))+
#   facet_grid(~Phase,scales="free",space="free")+
#   theme(axis.ticks.y=element_blank(),panel.spacing.x=unit(0.3, "lines"),axis.text.y = element_text(size = 8),
#         panel.spacing.y=unit(0.3, "lines"),strip.text.y = element_text(angle = 0),axis.text.x = element_blank())+
#   scale_fill_distiller(palette = "Spectral", direction = -1)
# 
# ggsave("BSH/BSH_ENB/GNPS_fecal_mtb/SFR25_017_heatmap_LDENB.pdf", plt,width = 20, height = 4)

mtb_sub<-mtb%>%
  dplyr::rename(FeatureID=`row.ID`)%>%
  dplyr::select(FeatureID, all_of(md_sub$filename))%>%
  gather(filename,peak_abun,-FeatureID)%>%
  left_join(.,annot,by="FeatureID")%>%
  left_join(.,md,by="filename")%>%
  mutate(log_10_peakabun=log10(peak_abun+1),
         Phase=factor(Phase,levels=c("Light","Dark")))

#36372 UDCA
#29648 TUDCA
mtb_tudcacondg<-mtb_sub%>%
  filter(FeatureID==36372 |FeatureID==29648)%>%
  dplyr::select(FeatureID,52:58,peak_abun)%>%
  spread(FeatureID,peak_abun)%>%
  dplyr::rename(UDCA=`36372`,TUDCA=`29648`)%>%
  mutate(uncongcongrat=UDCA/TUDCA)%>%
  #mutate(uncongcongrat=ifelse(uncongcongrat==Inf,2,uncongcongrat))%>%
  mutate(log_ratio=log10(uncongcongrat+1))

p <- ggplot(mtb_tudcacondg, aes(x=ENB,y=uncongcongrat,colour=ENB,fill=ENB)) +
  geom_boxplot(alpha=0.3)+
  geom_point()+
  facet_wrap(~Phase)+
  scale_fill_manual(values = c("gray70","#CC79A7","#009E73","#D55E00"))+
  scale_colour_manual(values = c("gray70","#CC79A7","#009E73","#D55E00"))+
  theme_pubr()+labs(y="UDCA/TUDCA")+
  theme(legend.position = "none")


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

stat.test<-mtb_sub%>%
#stat.testL<-mtb_subL%>%
#stat.testD<-mtb_subD%>%
  group_by(FeatureID,Compound_Name)%>%
  nest()%>%
  mutate(pvals = map(data, ~ TRF_unpr_ttest(.x)))%>%
  dplyr::select(-data)%>%
  unnest()
#write.table(stat.test,"BSH/BSH_ENB/GNPS_fecal_mtb/ttest_unpr_ttest_enb_pvals.txt",sep = "\t",row.names = FALSE, quote=FALSE)
#write.table(stat.test,"BSH/BSH_ENB/GNPS_fecal_mtb/ttest_unpr_ttest_enbL_pvals.txt",sep = "\t",row.names = FALSE, quote=FALSE)
write.table(stat.test,"BSH/BSH_ENB/GNPS_fecal_mtb/ttest_unpr_ttest_enbD_pvals.txt",sep = "\t",row.names = FALSE, quote=FALSE)
write.table(stat.test,"BSH/BSH_ENB/GNPS_fecal_mtb/ttest_unpr_ttest_enb_pvals_nofdr.txt",sep = "\t",row.names = FALSE, quote=FALSE)

ttest_res<-fread("BSH/BSH_ENB/GNPS_fecal_mtb/ttest_unpr_ttest_enb_pvals.txt")
dnbsh1sig<-ttest_res%>%filter(AZ51vDnBSH1<0.05 | AZ51vAZ52h<0.05| AZ51vLgBSH<0.05)%>%
  mutate(phase="combined")
ttest_resL<-fread("BSH/BSH_ENB/GNPS_fecal_mtb/ttest_unpr_ttest_enbL_pvals.txt")
dnbsh1sigL<-ttest_resL%>%filter(AZ51vDnBSH1<0.05| AZ51vAZ52h<0.05| AZ51vLgBSH<0.05)%>%
  mutate(phase="light")
ttest_resD<-fread("BSH/BSH_ENB/GNPS_fecal_mtb/ttest_unpr_ttest_enbD_pvals.txt")
dnbsh1sigD<-ttest_resD%>%filter(AZ51vDnBSH1<0.05| AZ51vAZ52h<0.05| AZ51vLgBSH<0.05)%>%
  mutate(phase="dark")

bbaa_int<-rbind(dnbsh1sig,dnbsh1sigL,dnbsh1sigD)
write.table(bbaa_int,"BSH/BSH_ENB/GNPS_fecal_mtb/ttest_unpr_ttest_enb_pvals_pless0.05_combLD.txt",sep = "\t",row.names = FALSE, quote=FALSE)

bbaa_intLD<-bbaa_int%>%
  filter(phase!="combined")

mtb_int<-mtb_sub%>%
  filter(FeatureID %in% dnbsh1sig$FeatureID)%>%
  mutate(name=paste("mtb",FeatureID,Compound_Name,sep=" "))%>%
  group_by(FeatureID,name,ENB)%>%
  summarise(mn_peak_abun=mean(peak_abun))%>%
  mutate(log_10_peakabun=log10(mn_peak_abun+1))%>%
  group_by(FeatureID)%>%mutate(Zscore=(log_10_peakabun - mean(log_10_peakabun))/sd(log_10_peakabun))

mtb_int<-mtb_sub%>%
  filter(FeatureID %in% bbaa_intLD$FeatureID)%>%
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

#ggsave("BSH/BSH_ENB/GNPS_fecal_mtb/SFR25_0210_sighitsBBAAENBcombphase_heatmap_log10.pdf", plt,width = 10, height = 15)
#ggsave("BSH/BSH_ENB/GNPS_fecal_mtb/SFR25_0210_sighitsBBAAENBcombphase_heatmap_zscore.pdf", plt,width = 10, height = 15)
#ggsave("BSH/BSH_ENB/GNPS_fecal_mtb/SFR25_0210_sighitsBBAAENBLDphase_heatmap_log10.pdf", plt,width = 10, height = 15)
ggsave("BSH/BSH_ENB/GNPS_fecal_mtb/SFR25_0210_sighitsBBAAENBLDphase_heatmap_zscore.pdf", plt,width = 6, height = 16)



# bbaa_results <- stat.test %>%
#   filter(!str_detect(Compound_Name, "Candidate"))

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


for(i in dnbsh1sig$FeatureID){
  p<-indv_plot(mtb_sub,i)
  #ggsave(paste("BSH/BSH_ENB/GNPS_fecal_mtb/indv_plt/sigDnBSH10.05/",i,".pdf",sep=""),p,width = 5, height = 3)
  #ggsave(paste("BSH/BSH_ENB/GNPS_fecal_mtb/indv_plt/sigDnBSH1L0.1/",i,".pdf",sep=""),p,width = 5, height = 3)
  ggsave(paste("BSH/BSH_ENB/GNPS_fecal_mtb/indv_plt/sigDnBSH1D0.1/",i,".pdf",sep=""),p,width = 5, height = 4)
}

p441_i<-indv_plot(mtb_sub,"441_i")
ggsave("BSH/BSH_ENB/GNPS_fecal_mtb/indv_plt/mtb_alsoinvitro/441_i.pdf",p441_i,width = 4, height = 3)
ggsave("BSH/BSH_ENB/GNPS_fecal_mtb/indv_plt/mtb_alsoinvitro/441_iLD.pdf",p441_i,width = 5, height = 3)

p37410<-indv_plot(mtb_sub,"37410")
ggsave("BSH/BSH_ENB/GNPS_fecal_mtb/indv_plt/mtb_alsoinvitro/37410.pdf",p37410,width = 4, height = 3)
ggsave("BSH/BSH_ENB/GNPS_fecal_mtb/indv_plt/mtb_alsoinvitro/37410LD.pdf",p37410,width = 5, height = 3)

p31640<-indv_plot(mtb_sub,"31640")
ggsave("BSH/BSH_ENB/GNPS_fecal_mtb/indv_plt/mtb_alsoinvitro/31640.pdf",p31640,width = 4, height = 3)
ggsave("BSH/BSH_ENB/GNPS_fecal_mtb/indv_plt/mtb_alsoinvitro/31640LD.pdf",p31640,width = 5, height = 3)

p23888<-indv_plot(mtb_sub,"23888")
ggsave("BSH/BSH_ENB/GNPS_fecal_mtb/indv_plt/mtb_alsoinvitro/23888.pdf",p23888,width = 4, height = 3)
ggsave("BSH/BSH_ENB/GNPS_fecal_mtb/indv_plt/mtb_alsoinvitro/23888LD.pdf",p23888,width = 5, height = 3)

p31501<-indv_plot(mtb_sub,"31501")
ggsave("BSH/BSH_ENB/GNPS_fecal_mtb/indv_plt/mtb_alsoinvitro/31501.pdf",p31501,width = 4, height = 3)
ggsave("BSH/BSH_ENB/GNPS_fecal_mtb/indv_plt/mtb_alsoinvitro/31501LD.pdf",p31501,width = 5, height = 3)


stat.test_lme<-mtb_sub%>%
  group_by(FeatureID,Compound_Name)%>%
  nest()%>%
  mutate(pvals = map(data, ~ metaTf_lme(.x)))%>%
  dplyr::select(-data)%>%
  unnest()
write.table(stat.test_lme,"BSH/BSH_ENB/GNPS_fecal_mtb/lme_enb_pvals.txt",sep = "\t",row.names = FALSE, quote=FALSE)

#pval vs. fold change

#light AZ51 vs. DnBSH1
summ_AZ51vDnBSH1<-fread("BSH/BSH_ENB/GNPS_fecal_mtb/lme_enb_pvals.txt")%>%
  dplyr::select(c(1,2,8,20))%>%
  mutate(light_AZ51vDnBSH1=-1*light_AZ51vDnBSH1)
write.table(summ_AZ51vDnBSH1,"BSH/BSH_ENB/GNPS_fecal_mtb/lme_enb_pvals_lightAZ51vDnBSH1_ONLY.txt",sep = "\t",row.names = FALSE,quote=FALSE)

summ_AZ51vDnBSH1<-fread("BSH/BSH_ENB/GNPS_fecal_mtb/lme_enb_pvals_lightAZ51vDnBSH1_ONLY.txt")

#Dark AZ51 vs. DnBSH1
summ_AZ51vDnBSH1<-fread("BSH/BSH_ENB/GNPS_fecal_mtb/lme_enb_pvals.txt")%>%
  dplyr::select(c(1,2,14,26))%>%
  mutate(dark_AZ51vDnBSH1=-1*dark_AZ51vDnBSH1)
write.table(summ_AZ51vDnBSH1,"BSH/BSH_ENB/GNPS_fecal_mtb/lme_enb_pvals_darkAZ51vDnBSH1_ONLY.txt",sep = "\t",row.names = FALSE,quote=FALSE)

#light AZ51 vs. AZ52
summ_AZ51vAZ52<-fread("BSH/BSH_ENB/GNPS_fecal_mtb/lme_enb_pvals.txt")%>%
  dplyr::select(c(1,2,7,19))%>%
  mutate(light_AZ51vAZ52=-1*light_AZ51vAZ52)
write.table(summ_AZ51vAZ52,"BSH/BSH_ENB/GNPS_fecal_mtb/lme_enb_pvals_lightAZ51vAZ52_ONLY.txt",sep = "\t",row.names = FALSE,quote=FALSE)

#dark AZ51 vs. AZ52
summ_AZ51vAZ52<-fread("BSH/BSH_ENB/GNPS_fecal_mtb/lme_enb_pvals.txt")%>%
  dplyr::select(c(1,2,13,25))%>%
  mutate(dark_AZ51vAZ52=-1*dark_AZ51vAZ52)
write.table(summ_AZ51vAZ52,"BSH/BSH_ENB/GNPS_fecal_mtb/lme_enb_pvals_darkAZ51vAZ52_ONLY.txt",sep = "\t",row.names = FALSE,quote=FALSE)

#light AZ51 vs. LgBSH
summ_AZ51vLgBSH<-fread("BSH/BSH_ENB/GNPS_fecal_mtb/lme_enb_pvals.txt")%>%
  dplyr::select(c(1,2,9,21))%>%
  mutate(light_AZ51vLgBSH=-1*light_AZ51vLgBSH)
write.table(summ_AZ51vLgBSH,"BSH/BSH_ENB/GNPS_fecal_mtb/lme_enb_pvals_lightAZ51vLgBSH_ONLY.txt",sep = "\t",row.names = FALSE,quote=FALSE)

#dark AZ51 vs. LgBSH
summ_AZ51vLgBSH<-fread("BSH/BSH_ENB/GNPS_fecal_mtb/lme_enb_pvals.txt")%>%
  dplyr::select(c(1,2,15,27))%>%
  mutate(dark_AZ51vLgBSH=-1*dark_AZ51vLgBSH)
write.table(summ_AZ51vLgBSH,"BSH/BSH_ENB/GNPS_fecal_mtb/lme_enb_pvals_darkAZ51vLgBSH_ONLY.txt",sep = "\t",row.names = FALSE,quote=FALSE)

p<-EnhancedVolcano(
  summ_AZ51vLgBSH,
  lab = summ_AZ51vLgBSH$Compound_Name,
  # summ_AZ51vAZ52,
  # lab = summ_AZ51vAZ52$Compound_Name,
  # summ_AZ51vDnBSH1,
  # lab = summ_AZ51vDnBSH1$Compound_Name,
  xlim=c(-2,2),
  ylim=c(0,7),
  # title = 'Dark AZ-51 vs. DnBSH1',
  # title = 'Dark AZ-51 vs. AZ-52',
  title = 'Dark AZ-51 vs. LgBSH',
  subtitle=NA,
  caption=NA,
  # x = "dark_AZ51vDnBSH1",
  # y = "dark_AZ51vDnBSH1_pval",
  # x = "dark_AZ51vAZ52",
  # y = "dark_AZ51vAZ52_pval",
  x = "dark_AZ51vLgBSH",
  y = "dark_AZ51vLgBSH_pval",
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

ggsave("BSH/BSH_ENB/GNPS_fecal_mtb/light_AZ51vDnBSH1_volplot.pdf",height=5, width=5)
ggsave("BSH/BSH_ENB/GNPS_fecal_mtb/dark_AZ51vDnBSH1_volplot.pdf",height=5, width=5)
ggsave("BSH/BSH_ENB/GNPS_fecal_mtb/light_AZ51vAZ52_volplot.pdf",height=5, width=5)
ggsave("BSH/BSH_ENB/GNPS_fecal_mtb/dark_AZ51vAZ52_volplot.pdf",height=5, width=5)
ggsave("BSH/BSH_ENB/GNPS_fecal_mtb/light_AZ51vLgBSH_volplot.pdf",height=5, width=5)
ggsave("BSH/BSH_ENB/GNPS_fecal_mtb/dark_AZ51vLgBSH_volplot.pdf",height=5, width=5)












BA_suppl_subset<-BA_t0ht48h_rmctrl%>%filter(BA_suppl=="GCA")
mtb_sub<-mtb%>%
  mutate(condition=case_when(BSH_strain=="Dny-BSH1"|BSH_strain=="Dny-BSH2"~"FT",
                             BSH_strain=="LCAG-95-BSH1208"|BSH_strain=="Ep-BSH101"~"NA",
                             BSH_strain=="Lgasseri-BSH"~"FA",
                              .default = "none"))%>%
  mutate(baseline_chng=ifelse(`0h`+`48h`==0, NA,baseline_chng))%>%
  mutate(cmpd_name_new=paste(FeatureID,Compound_Name,sep="_"))%>%
  filter(FeatureID %in% BA_suppl_subset$FeatureID)%>%
  mutate(BSH_strain=factor(BSH_strain,levels=c("EcAZ-1-cat","AZ-52","Dny-BSH1","Dny-BSH2","Lgasseri-BSH",
                                        "LCAG-95-BSH1208","Ep-BSH101","25MeOH","50MeOH","ctrl")),
         condition=factor(condition,levels=c("NA","FA","FT")))
  

indv_plot <- mtb_sub %>%
  group_by(FeatureID) %>%
  nest() %>%
  mutate(pval = map2(data, FeatureID,  ~ show_deconj(.x,.y)))


#plot three examples of conjugated bile acids

mtb_sub<-mtb%>%
  dplyr::rename(FeatureID=`row.ID`)%>%
  dplyr::select(FeatureID, all_of(md_sub$sample_name))%>%
  gather(sample_name,peak_abun,-FeatureID)%>%
  left_join(.,md,by="sample_name")%>%
  mutate(unique_id=paste(FeatureID,BA_suppl,sep="_"))%>%
  filter(unique_id %in% unique(BA_t0ht48h_rmctrl$unique_id))%>%
  left_join(.,annot,by="FeatureID")%>%
  mutate(log_10_peakabun=log10(peak_abun+1),
         timepoint=factor(timepoint, levels = c("0h","24h","48h")),
         BSH_strain=factor(BSH_strain,levels=c("EcAZ-1-cat","AZ-52","Dny-BSH1","Dny-BSH2","Lgasseri-BSH",
                                               "LCAG-95-BSH1208","Ep-BSH101","25MeOH","50MeOH","ctrl")),
         BA_suppl=factor(BA_suppl,levels=c("BHI-50DMSO","BHI","GCA","GCDCA","GDCA","GLCA","GUDCA","TCA",
                                           "TCDCA","TDCA","TLCA","TUDCA")))%>%
  filter(!(BA_suppl=="BHI-50DMSO"|BA_suppl=="BHI"))%>%
  filter(timepoint!="24h")


plt_sig_mtbs<-function(mtb,mtb_id){
  cmp_name<-unique(mtb$Compound_Name)
  p<-ggplot(data=mtb, aes(x=timepoint, y=log_10_peakabun, colour=BA_suppl, fill=BA_suppl)) +
    geom_boxplot(color = "black",alpha=0.3)+
    geom_point(shape=21, color = "black", position=position_jitterdodge())+
    facet_wrap(~BSH_strain, nrow=1)+
    theme_classic()+
    labs(title=paste(cmp_name," (",mtb_id,")",sep=""),
         x="timepoint", y="log10(peak_area+1)")+
    theme(legend.position = "right", plot.title = element_text(size = 12),axis.title.x = element_text(size = 10))
  ggsave(paste("BSH/BSH_ENB/GNPS/ttest_unpr_results/indv_mtb_plt/SFR24_0501_mtb",mtb_id,"_notpaired.pdf",sep=""), p,width = 8, height = 2)
  
}

indv_plot<-mtb_sub%>%
  group_by(FeatureID)%>%
  nest()%>%
  mutate(pval = map2(data, FeatureID,  ~ plt_sig_mtbs(.x,.y)))

###########################################################
#get change from baseline
annot<-fread("BSH/BSH_ENB/GNPS/nf_output/networking/library-results-merged_results_with_gnps_justBA_collapse.tsv")%>%
  mutate(Compound_Name = gsub("\\|.*", "", Compound_Name))%>%
  mutate(Compound_Name = gsub("\\(delta mass.*", "", Compound_Name))%>%
  dplyr::select(FeatureID,Compound_Name)

md<-fread("BSH/BSH_ENB/GNPS/metadata_filename/ZT_BSH_metadata_justBA.tsv")%>%
  mutate(cult_plt=case_when(plate=="P4"|plate=="P5"|plate=="P6"~"P2",
                            plate=="P7"|plate=="P8"|plate=="P9"~"P3",
                            .default = "P1"))%>%
  filter(timepoint!="24h")
  
mtb<-fread("BSH/BSH_ENB/GNPS/nf_output/clustering/featuretable_reformated_justBA_summcollapse_cln.csv")%>%
  dplyr::rename(FeatureID=`row.ID`)%>%
  dplyr::select(FeatureID,any_of(md$sample_name))%>%
  gather(sample_name,peak_area,-FeatureID)%>%
  left_join(.,md,by="sample_name")%>%
  mutate(name_new=paste(BA_suppl,BSH_strain,sep="-"))%>%
  dplyr::select(FeatureID,name_new,BA_suppl,BSH_strain,cult_plt,plate_location,timepoint,peak_area)%>%
  filter(timepoint!=""& !(BA_suppl%in% c("BHI","BHI-50DMSO")))%>%
  pivot_wider(names_from = timepoint, values_from = peak_area)%>%
  mutate(baseline_chng=`48h`-`0h`)%>%
  #mutate(baseline_chng=`48h`/(`0h`+1))%>%
  right_join(.,annot,by="FeatureID")%>%
  filter(`0h`+`48h`>0)#%>%
  # filter(baseline_chng>0)

p <- ggplot(mtb, aes(x=log10(baseline_chng+0.000001), colour=BSH_strain)) + 
  geom_density()+ facet_wrap(~BA_suppl, nrow=2)+
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")+
  theme_pubr()+labs(title="all BAs")

#ggsave("BSH/BSH_ENB/GNPS/change_baseline/allBA_density.pdf", p,width = 10, height = 4)
ggsave("BSH/BSH_ENB/GNPS/change_baseline/allBA_density_48-0h.pdf", p,width = 10, height = 4)

BA_amine<-annot%>% filter(str_detect(Compound_Name, "[^- ]-[^- ]"))

BA_candi<-mtb%>%filter(!(FeatureID %in% BA_amine$FeatureID))%>%
  mutate(condition=case_when(BSH_strain=="Dny-BSH1"|BSH_strain=="Dny-BSH2"~"FT",
                                      BSH_strain=="LCAG-95-BSH1208"|BSH_strain=="Ep-BSH101"~"NA",
                                      BSH_strain=="Lgasseri-BSH"~"FA",
                                      .default = "none"),
         BSH_strain=factor(BSH_strain,levels=c("EcAZ-1-cat","AZ-52","Dny-BSH1","Dny-BSH2","Lgasseri-BSH",
                                               "LCAG-95-BSH1208","Ep-BSH101","25MeOH","50MeOH","ctrl")))%>%
    mutate(condition=factor(condition,levels=c("NA","FA","FT")))

p <- ggplot(BA_candi, aes(x=log10(baseline_chng+0.000001), colour=BA_suppl)) + 
  geom_density()+ facet_wrap(~BSH_strain, nrow=2)+
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")+
  theme_pubr()+labs(title="candidate BA")
#ggsave("BSH/BSH_ENB/GNPS/change_baseline/candBA_density.pdf", p,width = 8, height = 4)
ggsave("BSH/BSH_ENB/GNPS/change_baseline/candBA_density_48h-0h.pdf", p,width = 8, height = 4)


show_deconj<-function(mtb, BA,mtb_id){
  BA_candi_x<-mtb%>%filter(BA_suppl==BA)%>%
    mutate(cmpd_name_new=paste(FeatureID,Compound_Name,sep="_"))%>%
    filter(FeatureID==mtb_id)
    #filter(`0h`>0)
  p <- ggplot(BA_candi_x, aes(x=fct_rev(BSH_strain),y=log10(baseline_chng+0.000001),colour=condition)) +
    geom_boxplot()+
    #geom_point()+
    facet_wrap(~BA_suppl)+ coord_flip()+
    scale_colour_manual(values=c("#0072B2","#D55E00","#009E73"))+
    geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
    theme_pubr()+labs(title=unique(BA_candi_x$cmpd_name_new))+theme(legend.position = "none")
  
  return(p)
}
BA_candi_GCA<-show_deconj(BA_candi,"GCA","238_i")
ggsave("BSH/BSH_ENB/GNPS/change_baseline/GCA_238_i_deconj.pdf", BA_candi_GCA,width =4, height = 3)
BA_candi_GCA<-show_deconj(BA_candi,"GCA","185_i")
ggsave("BSH/BSH_ENB/GNPS/change_baseline/GCA_185_i_deconj.pdf", BA_candi_GCA,width =4, height = 3)
BA_candi_GCA<-show_deconj(BA_candi,"GCA","43992")
ggsave("BSH/BSH_ENB/GNPS/change_baseline/GCA_43992_deconj.pdf", BA_candi_GCA,width =4, height = 3)

BA_candi_TLCA<-show_deconj(BA_candi,"TLCA","240_i")
ggsave("BSH/BSH_ENB/GNPS/change_baseline/TLCA_240_i_deconj.pdf", BA_candi_TLCA,width =4, height = 3)
BA_candi_TDCA<-show_deconj(BA_candi,"TDCA","88_i")
ggsave("BSH/BSH_ENB/GNPS/change_baseline/TDCA_88_i_deconj.pdf", BA_candi_TDCA,width =4, height = 3)
BA_candi_TCA<-show_deconj(BA_candi,"TCA","201_i")
ggsave("BSH/BSH_ENB/GNPS/change_baseline/TCA_201_i_deconj.pdf", BA_candi_TCA,width =4, height = 3)
# BA_candi_TCA<-show_deconj(BA_candi,"TCA","50529")




BA_amine<-mtb%>%filter(FeatureID %in% BA_amine$FeatureID)%>%
  mutate(BSH_strain=factor(BSH_strain,levels=c("EcAZ-1-cat","AZ-52","Dny-BSH1","Dny-BSH2","Lgasseri-BSH",
                                        "LCAG-95-BSH1208","Ep-BSH101","25MeOH","50MeOH","ctrl")),
         BA_suppl=factor(BA_suppl,levels=c("GCA","TCA","GCDCA","TCDCA","GDCA","TDCA","GLCA","TLCA","GUDCA","TUDCA")))%>%
  group_by(FeatureID,name_new,BSH_strain,BA_suppl)%>%
  summarise(mn_chg_baseline=mean(baseline_chng))%>%
  left_join(annot,by="FeatureID")%>%
  mutate(cmpd_name_new=paste(FeatureID,Compound_Name))

BA_amine <- BA_amine %>%
  mutate(cmpd_name_new = factor(cmpd_name_new, levels = unique(BA_amine$cmpd_name_new[order(BA_amine$Compound_Name)])))

p <- ggplot(BA_amine, aes(x=fct_rev(BSH_strain), y=log10(mn_chg_baseline), fill=cmpd_name_new)) + 
  geom_bar(stat="identity")+
  facet_wrap(~BA_suppl, nrow=5)+
  scale_fill_manual(values = c("#B4A7D6","#582C83","#A6CEE3","#3789AD","#87CA6A","#3DA533","#238B45","#979C62",
                               "#F98F8E","#E94330","#FDA440","#FF7F00","#783F04")) +
  #scale_y_continuous(expand=c(0,0), breaks = seq(-10, 20,by = 5))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
  theme_bw()+labs(title="amine BAs")+coord_flip()+
  theme(legend.position = "right")

#ggsave(paste("BSH/BSH_ENB/GNPS/change_baseline/amine_conj_summ_tall.pdf",sep=""), p,width = 6, height = 7)
ggsave(paste("BSH/BSH_ENB/GNPS/change_baseline/amine_conj_summ_tall_48h-0h.pdf",sep=""), p,width = 6, height = 7)

###########################################################

#load the counts table from Ipsita to create heatmap

annot<-fread("BSH/BSH_ENB/BSH_Zarrinpar_libhit.csv")%>%
  dplyr::rename(FeatureID=Scan)%>%
  mutate(FeatureID=as.character(FeatureID))%>%
  dplyr::select(FeatureID,Compound_Name)

order<-c("TCA_EcAZ-1-cat_0h","TCA_EcAZ-1-cat_48h",
         "TCA_AZ-52_0h","TCA_AZ-52_48h",
         "TCA_Dny_BSH1_0h","TCA_Dny_BSH1_48h",
         "TCA_Dny_BSH2_0h","TCA_Dny_BSH2_48h",
         "TCA_Lgasseri_0h","TCA_Lgasseri_48h", 
         "TCA_Ep_BSH101_0h","TCA_Ep_BSH101_48h",
         "TCA_LCAG-95_0h","TCA_LCAG-95_48h",
         "TDCA_EcAZ-1-cat_0h","TDCA_EcAZ-1-cat_48h",
         "TDCA_AZ-52_0h", "TDCA_AZ-52_48h",
         "TDCA_Dny_BSH1_0h","TDCA_Dny_BSH1_48h",
         "TDCA_Dny_BSH2_0h","TDCA_Dny_BSH2_48h",
         "TDCA_Lgasseri_0h","TDCA_Lgasseri_48h",
         "TDCA_Ep_BSH101_0h","TDCA_Ep_BSH101_48h",
         "TDCA_LCAG-95_0h","TDCA_LCAG-95_48h",
         "Dny_BSH1_blank_0h","Dny_BSH1_blank_48h",
         "Dny_BSH2_blank_0h","Dny_BSH2_blank_48h",
         "Lgasseri_blank_0h", "Lgasseri_blank_48h",
         "Ep_BSH101_blank_0h","Ep_BSH101_blank_48h",
         "LCAG-95_blank_0h","LCAG-95_blank_48h",
         "Media_TCA_blank_nobacteria_0h","Media_TCA_blank_nobacteria_48h",
         "Media_TDCA_blank_nobacteria_0h","Media_TDCA_blank_nobacteria_48h",
         "Media_blank_nobacteria_0h","Media_blank_nobacteria_48h",
         "EcAZ-1-cat_blank_0h","EcAZ-1-cat_blank_48h",
         "AZ-52_blank_0h","AZ-52_blank_48h",
         "Blank")

mtb<-read_csv("BSH/BSH_ENB/BSH_feature_with_metadat_sub.csv")%>%
  group_by(ATTRIBUTE_name) %>%
  summarize_all(mean)

mtb <- mtb[order(match(mtb$ATTRIBUTE_name, order)), ]

mtb<-mtb%>%column_to_rownames("ATTRIBUTE_name")%>% as.matrix()%>%
  log10()

mtb[ mtb== -Inf] <- 0

# mtb<-read_csv("BSH/BSH_ENB/BSH_feature_with_metadat_sub.csv")%>%
#   gather(FeatureID,peak_area,-ATTRIBUTE_name)%>%
#   left_join(.,annot, by="FeatureID")

p<-pheatmap(mtb,cluster_rows=FALSE,row_order = order)
ggsave("BSH/BSH_ENB/SFR24_0307_heatmapTCATDCA.pdf", p,width = 10, height = 8)


mtb_annot<-mtb%>%as.data.frame()%>%
  rownames_to_column("ATTRIBUTE_name")%>%
  gather(FeatureID,log_peak,-ATTRIBUTE_name)%>%
  left_join(.,annot,by="FeatureID")%>%
  mutate(name=paste(FeatureID,Compound_Name, sep=" "))%>%
  dplyr::select(name,ATTRIBUTE_name,log_peak)%>%
  spread(key=name, value=log_peak)
  
mtb_annot <- mtb_annot[order(match(mtb_annot$ATTRIBUTE_name, order)), ]
mtb_annot<-mtb_annot%>%
  rownames_to_column()%>%
  dplyr::select(-rowname)%>%
  column_to_rownames("ATTRIBUTE_name")%>% as.matrix()

p<-pheatmap(mtb_annot,cluster_rows=FALSE,row_order = order)
ggsave("BSH/BSH_ENB/SFR24_0307_heatmapTCATDCA_annot.pdf", p,width = 10, height = 30)

#subset heatmap 

order_new<-c("TCA_EcAZ-1-cat_0h",
         "TCA_AZ-52_0h",
         "TCA_Dny_BSH1_0h",
         "TCA_Dny_BSH2_0h",
         "TCA_Lgasseri_0h",
         "TCA_Ep_BSH101_0h",
         "TCA_LCAG-95_0h",
         "TDCA_EcAZ-1-cat_0h",
         "TDCA_AZ-52_0h", 
         "TDCA_Dny_BSH1_0h",
         "TDCA_Dny_BSH2_0h",
         "TDCA_Lgasseri_0h",
         "TDCA_Ep_BSH101_0h",
         "TDCA_LCAG-95_0h",
         "TCA_EcAZ-1-cat_48h",
         "TCA_AZ-52_48h",
         "TCA_Dny_BSH1_48h",
         "TCA_Dny_BSH2_48h",
         "TCA_Lgasseri_48h",
         "TCA_Ep_BSH101_48h",
         "TCA_LCAG-95_48h",
         "TDCA_EcAZ-1-cat_48h",
         "TDCA_AZ-52_48h",
         "TDCA_Dny_BSH1_48h",
         "TDCA_Dny_BSH2_48h",
         "TDCA_Lgasseri_48h",
         "TDCA_Ep_BSH101_48h",
         "TDCA_LCAG-95_48h")

mtb_noblanks<-mtb[1:28,]

mtb_noblanks<-mtb_noblanks%>%as.data.frame()%>%
  rownames_to_column("ATTRIBUTE_name")

mtb_noblanks<- mtb_noblanks[order(match(mtb_noblanks$ATTRIBUTE_name, order_new)), ]

mtb_noblanks<-mtb_noblanks%>%
  rownames_to_column()%>%
  dplyr::select(-rowname)%>%
  column_to_rownames("ATTRIBUTE_name")%>% as.matrix()

p<-pheatmap(mtb_noblanks,cluster_rows=FALSE)
ggsave("BSH/BSH_ENB/SFR24_0307_heatmapTCATDCA_noblanks_048h.pdf", p,width = 10, height = 4)

mtb_noblanks_annot<-mtb_noblanks%>%as.data.frame()%>%
  rownames_to_column("ATTRIBUTE_name")%>%
  gather(FeatureID,log_peak,-ATTRIBUTE_name)%>%
  left_join(.,annot,by="FeatureID")%>%
  mutate(name=paste(FeatureID,Compound_Name, sep=" "))%>%
  dplyr::select(name,ATTRIBUTE_name,log_peak)%>%
  spread(key=name, value=log_peak)

mtb_noblanks_annot <- mtb_noblanks_annot[order(match(mtb_noblanks_annot$ATTRIBUTE_name, order_new)), ]
mtb_noblanks_annot<-mtb_noblanks_annot%>%
  rownames_to_column()%>%
  dplyr::select(-rowname)%>%
  column_to_rownames("ATTRIBUTE_name")%>% as.matrix()

p<-pheatmap(mtb_noblanks_annot,cluster_rows=FALSE,row_order = order)
ggsave("BSH/BSH_ENB/SFR24_0307_heatmapTCATDCA_noblanks_048h_annot.pdf", p,width = 10, height = 25)

#TCA and TDCA

TCA<-c("12997","12865","9195","12880","12814","12897","12891","12901","12878","12861","12922","12902",
       "12905","12875","11900","9161","9273")
TDCA<-c("16285","16162","16167","16231","16330","16226","16235")

mtb_forplt<-read_csv("BSH/BSH_ENB/BSH_feature_with_metadat_sub.csv")%>%
  gather(FeatureID,peak_area,-ATTRIBUTE_name)%>%
  filter(FeatureID %in% TCA)%>%
  #filter(FeatureID %in% TDCA)%>%
  filter(!grepl("blank", ATTRIBUTE_name,ignore.case = TRUE))%>%
  separate(ATTRIBUTE_name,c("BA","BSH_origin","BSH_name","timepoint"),sep="_",remove=FALSE)%>%
  mutate(timepoint=ifelse(is.na(timepoint),BSH_name,timepoint),
         BSH_name=ifelse(BSH_name=="0h"|BSH_name=="48h",NA,BSH_name),
         name=sub("_[^_]+$", "", ATTRIBUTE_name),
         peak_area=log10(peak_area+1))%>%
  filter(BA=="TCA")%>%
  #filter(BA=="TDCA")%>%
  mutate(name=factor(name,levels=c("TCA_EcAZ-1-cat","TCA_AZ-52","TCA_Dny_BSH1","TCA_Dny_BSH2",
                                   "TCA_Lgasseri","TCA_Ep_BSH101","TCA_LCAG-95")))%>%
  # mutate(name=factor(name,levels=c("TDCA_EcAZ-1-cat","TDCA_AZ-52","TDCA_Dny_BSH1","TDCA_Dny_BSH2",
  #                                  "TDCA_Lgasseri","TDCA_Ep_BSH101","TDCA_LCAG-95")))%>%
  dplyr::select(-ATTRIBUTE_name)%>%
  arrange(FeatureID)

mtb_forplt<-mtb_forplt%>%spread(timepoint,peak_area)

stat.test<-mtb_forplt %>%
  group_by(name) %>%
  wilcox_test(data =., peak_area ~ timepoint) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")

stat.test

p<-ggpaired(mtb_forplt, cond1="0h",cond2="48h",fill="name", id="FeatureID")+
  facet_wrap(~name, nrow=1)+
  labs(title="TDCA", x="timepoint", y="log10(peak_area+1)")+
  theme(legend.position = "none")
ggsave("BSH/BSH_ENB/SFR24_0403_TCA048h_quant_paired.pdf", p,width = 12, height = 3)
ggsave("BSH/BSH_ENB/SFR24_0403_TDCA048h_quant_paired.pdf", p,width = 12, height = 3)

p<-ggplot(data=mtb_forplt, aes(x=timepoint, y=peak_area, colour=name, fill=name)) +
  geom_boxplot(color = "black",alpha=0.3)+
  geom_point(shape=21, color = "black", position=position_jitterdodge())+
  facet_wrap(~name, nrow=1)+
  theme_classic()+
  labs(title="TCA", x="timepoint", y="log10(peak_area+1)")+
  theme(legend.position = "none", plot.title = element_text(size = 12),axis.title.x = element_text(size = 10))
ggsave("BSH/BSH_ENB/SFR24_0403_TCA048h_quant_notpaired.pdf", p,width = 8, height = 2)
ggsave("BSH/BSH_ENB/SFR24_0403_TDCA048h_quant_notpaired.pdf", p,width = 8, height = 2)
