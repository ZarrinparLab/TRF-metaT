setwd("~/Notebooks/sfloresr/TRF-metaT/")

library(tidyverse)
library(data.table)
library(EnhancedVolcano)
library(emmeans)
library(nlme)
library(multcomp)

###########################################################
#paths
dat_feattab<-"data/ENB_invivo_metab/GNPS_fecal_mtb/GNPS_outputs/clustering/featuretable_reformated_justBA_summcollapse.csv"
dat_metadata<-"data/ENB_invivo_metab/GNPS_fecal_mtb/BSH_ENB_invivo_metadata.txt"
dat_library<-"data/ENB_invivo_metab/GNPS_fecal_mtb/GNPS_outputs/library/merged_results_with_gnps_justBA_collapse.tsv"
dat_path<-"data/ENB_invivo_metab/GNPS_fecal_mtb/"
fig_path<-"figures/ENB_invivo_metab/"
###########################################################
#functions

metaTf_lme <- function(mtb) {
  tryCatch({
    m1 <- lme(log_10_peakabun ~ ENB * Phase,
              random = ~1 | mouseid,
              data = mtb)
    aov <- anova(m1)
    ph <- emmeans(m1, pairwise ~ ENB | Phase)
    ph_est <- as.data.frame(ph$contrasts)$estimate
    ph_pval <- as.data.frame(ph$contrasts)$p.value
    df <- data.frame(Intercept_pval = aov$p-value[[1]],
                     ENB_pval = aov$p-value[[2]],
                     Phase_pval = aov$p-value[[3]],
                     ENB_Phase_pval = aov$p-value[[4]],
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

volcano_plt<-function(dat,title_n,x_n,y_n,max1,max2){
  p<-EnhancedVolcano(
    dat,
    lab = dat$Compound_Name,
    xlim=c(-1*max1,max1),
    ylim=c(0,max2),
    title = title_n,
    subtitle=NA,
    caption=NA,
    x = x_n,
    y = y_n,
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
  
  return(p)
}

###########################################################
#cleaning files

annot<-fread(dat_library)%>%
  mutate(Compound_Name = gsub("\\|.*", "", Compound_Name))%>%
  mutate(FeatureID=as.character(FeatureID))

mtb<-fread(dat_feattab)
md<-fread(dat_metadata)

###########################################################
#run lme
md_sub<-md%>%
  filter(collection_time!="")

mtb_sub<-mtb%>%
  dplyr::rename(FeatureID=row.ID)%>%
  dplyr::select(FeatureID, all_of(md_sub$filename))%>%
  gather(filename,peak_abun,-FeatureID)%>%
  left_join(.,annot,by="FeatureID")%>%
  left_join(.,md,by="filename")%>%
  mutate(log_10_peakabun=log10(peak_abun+1),
         Phase=factor(Phase,levels=c("Light","Dark")))

stat.test_lme<-mtb_sub%>%
  group_by(FeatureID,Compound_Name)%>%
  nest()%>%
  mutate(pvals = map(data, ~ metaTf_lme(.x)))%>%
  dplyr::select(-data)%>%
  unnest()
write.table(stat.test_lme,paste0(dat_path,"lme_enb_pvals.txt"), #Table S6
            sep = "\t",row.names = FALSE, quote=FALSE)

###########################################################
#plot volcano of pval vs. fold change--Figure S6B

#light AZ51 vs. AZ52
summ_AZ51vAZ52<-fread(paste0(dat_path,"lme_enb_pvals.txt"))%>%
  dplyr::select(c(1,2,7,19))%>%
  mutate(light_AZ51vAZ52=-1*light_AZ51vAZ52)

p<-volcano_plt(summ_AZ51vAZ52,'Light AZ-51 vs. AZ-52',"light_AZ51vAZ52","light_AZ51vAZ52_pval",1.5,4)
ggsave(paste0(fig_path,"light_AZ51vAZ52_volplot.pdf"),height=5, width=5)

#dark AZ51 vs. AZ52
summ_AZ51vAZ52<-fread(paste0(dat_path,"lme_enb_pvals.txt"))%>%
  dplyr::select(c(1,2,13,25))%>%
  mutate(dark_AZ51vAZ52=-1*dark_AZ51vAZ52)

p<-volcano_plt(summ_AZ51vAZ52,'Dark AZ-51 vs. AZ-52',"dark_AZ51vAZ52","dark_AZ51vAZ52_pval",1.5,6)
ggsave(paste0(fig_path,"dark_AZ51vAZ52_volplot.pdf"),height=5, width=5)

#light AZ51 vs. DnBSH1
summ_AZ51vDnBSH1<-fread(paste0(dat_path,"lme_enb_pvals.txt"))%>%
  dplyr::select(c(1,2,8,20))%>%
  mutate(light_AZ51vDnBSH1=-1*light_AZ51vDnBSH1)

p<-volcano_plt(summ_AZ51vDnBSH1,'Light AZ-51 vs. DnBSH1',"light_AZ51vDnBSH1","light_AZ51vDnBSH1_pval",1.5,4)
ggsave(paste0(fig_path,"light_AZ51vDnBSH1_volplot.pdf"),height=5, width=5)

#Dark AZ51 vs. DnBSH1
summ_AZ51vDnBSH1<-fread(paste0(dat_path,"lme_enb_pvals.txt"))%>%
  dplyr::select(c(1,2,14,26))%>%
  mutate(dark_AZ51vDnBSH1=-1*dark_AZ51vDnBSH1)

p<-volcano_plt(summ_AZ51vDnBSH1,'Dark AZ-51 vs. DnBSH1',"dark_AZ51vDnBSH1","dark_AZ51vDnBSH1_pval",1.5,4)
ggsave(paste0(fig_path,"dark_AZ51vDnBSH1_volplot.pdf"),height=5, width=5)

#light AZ51 vs. LgBSH
summ_AZ51vLgBSH<-fread(paste0(dat_path,"lme_enb_pvals.txt"))%>%
  dplyr::select(c(1,2,9,21))%>%
  mutate(light_AZ51vLgBSH=-1*light_AZ51vLgBSH)

p<-volcano_plt(summ_AZ51vLgBSH,'Light AZ-51 vs. LgBSH',"light_AZ51vLgBSH","light_AZ51vLgBSH_pval",1.5,4)
ggsave(paste0(fig_path,"light_AZ51vLgBSH_volplot.pdf"),height=5, width=5)

#dark AZ51 vs. LgBSH
summ_AZ51vLgBSH<-fread(paste0(dat_path,"lme_enb_pvals.txt"))%>%
  dplyr::select(c(1,2,15,27))%>%
  mutate(dark_AZ51vLgBSH=-1*dark_AZ51vLgBSH)

p<-volcano_plt(summ_AZ51vLgBSH,'Dark AZ-51 vs. LgBSH',"dark_AZ51vLgBSH","dark_AZ51vLgBSH_pval",2,7)
ggsave(paste0(fig_path,"dark_AZ51vLgBSH_volplot.pdf"),height=5, width=5)