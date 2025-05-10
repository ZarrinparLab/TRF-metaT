setwd("~/Notebooks/sfloresr/TRF-metaT/")

library(tidyverse)
library(data.table)
library(ggpubr)
library(viridis)

###########################################################
#paths
notnorm_data<-"data/bsh_analysis/species_pfam_metaT/species_pfam.tsv"
norm_data<-"data/bsh_analysis/species_pfam_metaT/species_pfam-TPM.tsv"
bdm_NAFA_light<-"data/bsh_analysis/species_pfam_metaT/birdman_outputs/species_pfam_BSHonly_clean_rmzero_noNT_lightNA.beta_var.tsv"
bdm_NAFA_dark<-"data/bsh_analysis/species_pfam_metaT/birdman_outputs/species_pfam_BSHonly_clean_rmzero_noNT_darkNA.beta_var.tsv"
bdm_FAFT_light<-"data/bsh_analysis/species_pfam_metaT/birdman_outputs/species_pfam_BSHonly_clean_rmzero_noNT_light.beta_var.tsv"
bdm_FAFT_dark<-"data/bsh_analysis/species_pfam_metaT/birdman_outputs/species_pfam_BSHonly_clean_rmzero_noNT_dark.beta_var.tsv"
dat_metadata<-"data/bsh_analysis/metaT_metadata_ztcat_noNT.txt"
dat_path<-"data/bsh_analysis/"
fig_path<-"figures/bsh_analysis/"
###########################################################
#functions

subset_dat<-function(dt){
  samps_list<-c("#FeatureID","cFA01a_S39.","cFA01b_S45.","cFA05a_S40.","cFA05b_S46.","cFA09a_S41.",
                "cFA09b_S47.","cFA13a_S36.","cFA13b_S42.","cFA17a_S37.","cFA17b_S43.",
                "cFA21a_S38.","cFA21b_S44.","cFT01a_S51.","cFT01b_S57.","cFT05a_S52.",
                "cFT05b_S58.","cFT09a_S53.","cFT09b_S59.","cFT13a_S48.","cFT13b_S54.",
                "cFT17a_S49.","cFT17b_S55.","cFT21a_S50.","cFT21b_S56.","cNA01a_S4.",
                "cNA01b_S10.","cNA01c_S16.","cNA05a_S5.","cNA05b_S11.","cNA05c_S17.",
                "cNA09a_S6.","cNA09b_S12.","cNA09c_S18.","cNA13a_S1.","cNA13b_S7.",
                "cNA13c_S13.","cNA17a_S2.","cNA17b_S8.","cNA17c_S14.","cNA21a_S3.",
                "cNA21b_S9.","cNA21c_S15.")
  dat<-fread(dt)%>%dplyr::select(all_of(samps_list))
  names(dat) <- c("FeatureID","cFA01a","cFA01b","cFA05a","cFA05b","cFA09a","cFA09b","cFA13a","cFA13b","cFA17a","cFA17b",
                  "cFA21a","cFA21b","cFT01a","cFT01b","cFT05a","cFT05b","cFT09a","cFT09b","cFT13a","cFT13b",
                  "cFT17a","cFT17b","cFT21a","cFT21b","cNA01a","cNA01b","cNA01c","cNA05a","cNA05b","cNA05c",
                  "cNA09a","cNA09b","cNA09c","cNA13a","cNA13b","cNA13c","cNA17a","cNA17b","cNA17c","cNA21a",
                  "cNA21b","cNA21c")
  return(dat)
}

subset_dat_bshLD<-function(dt,phase){
  if(phase=="light"){
    samps_list<-c("FeatureID","cFA01a","cFA01b","cFA05a","cFA05b","cFA09a","cFA09b",
                  "cFT01a","cFT01b","cFT05a","cFT05b","cFT09a","cFT09b",
                  "cNA01a","cNA01b","cNA01c","cNA05a","cNA05b","cNA05c",
                  "cNA09a","cNA09b","cNA09c")
  }
  else{
    samps_list<-c("FeatureID","cFA13a","cFA13b","cFA17a","cFA17b",
                  "cFA21a","cFA21b","cFT13a","cFT13b",
                  "cFT17a","cFT17b","cFT21a","cFT21b",
                  "cNA13a","cNA13b","cNA13c","cNA17a","cNA17b","cNA17c","cNA21a",
                  "cNA21b","cNA21c")
  }
  dat<-dt%>%dplyr::select(all_of(samps_list))
  return(dat)
}

bsh_dist_dat<-function(dat,phase){
  bsh_dat<-fread(dat)%>%
    filter(FeatureID %in% selfeat$`Feature ID`)%>%
    mutate(FeatureID=gsub("\\|PF02275.21","",FeatureID))%>%
    gather(sample_name,TPM_counts,-FeatureID)%>%
    left_join(.,mdT,by="sample_name")%>%
    group_by(FeatureID,condition,phase)%>%summarise(mn_TPM=mean(TPM_counts), log_mn_TPM=log10(mean(TPM_counts)+1))%>%
    mutate(condition=factor(condition, levels = c("FT", "FA","NA")),
           phase=factor(phase, levels = c("light","dark")))
  
  if(phase=="light"){
    bsh_dat<-bsh_dat%>%filter(phase=="light")
  }
  else{
    bsh_dat<-bsh_dat%>%filter(phase=="dark")
  }
  bsh_dat<-bsh_dat%>%
    group_by(FeatureID)%>%mutate(sum_TPM=sum(mn_TPM), sum_log_TPM=sum(log_mn_TPM))%>%
    filter(mn_TPM > 0, !grepl('[[:digit:]]+', FeatureID))
  
  return(bsh_dat)
}

bsh_dist_plt<-function(dat,phase){
  p<-ggplot(data=dat, aes(x=reorder(FeatureID, sum_log_TPM), y=log_mn_TPM, fill=condition)) +
    geom_bar(stat="identity") + coord_flip() + theme_classic()+
    scale_fill_manual(values=c("#009E73","#D55E00","#0072B2"))+
    labs(y="log10(mean_TPM)",x="Species",title=phase)+
    theme(legend.position = "none")+
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0))
  return(p)
}

bdm_nafa_dat<-function(dat,phase){
  bsh_bdmdat<-fread(dat)%>%
    dplyr::rename(ratio=`C(condition, Treatment('NA'))[T.FA]_mean`,
                  FeatureID=Feature)%>%
    mutate(`C(condition, Treatment('NA'))[T.FA]_hdi`=gsub("[(]|[)]","",`C(condition, Treatment('NA'))[T.FA]_hdi`))%>%
    separate(`C(condition, Treatment('NA'))[T.FA]_hdi`,c("min","max"), sep=",")%>%
    mutate(min=as.numeric(min),
           max=as.numeric(max),
           diff=max-min)%>%
    dplyr::select(FeatureID, ratio, min, max, diff)
  if(phase=="light"){
    bsh_bdmdat<-bsh_bdmdat%>%mutate(phase="light")
  }
  else{
    bsh_bdmdat<-bsh_bdmdat%>%mutate(phase="dark")
  }
  return(bsh_bdmdat)
}

bdm_faft_dat<-function(dat,phase){
  bsh_bdmdat<-fread(dat)%>%
    dplyr::rename(ratio=`C(condition, Treatment('FA'))[T.FT]_mean`,
                  FeatureID=Feature)%>%
    mutate(`C(condition, Treatment('FA'))[T.FT]_hdi`=gsub("[(]|[)]","",`C(condition, Treatment('FA'))[T.FT]_hdi`))%>%
    separate(`C(condition, Treatment('FA'))[T.FT]_hdi`,c("min","max"), sep=",")%>%
    mutate(min=as.numeric(min),
           max=as.numeric(max),
           diff=max-min)%>%
    dplyr::select(FeatureID, ratio, min, max,diff)

  if(phase=="light"){
    bsh_bdmdat<-bsh_bdmdat%>%mutate(phase="light")
  }
  else{
    bsh_bdmdat<-bsh_bdmdat%>%mutate(phase="dark")
  }
  return(bsh_bdmdat)
}

bdm_plt<-function(dat){
  p<-ggplot(dat, aes(x =Feature , y = ratio, ymin = min, ymax = max, color=phase, group=phase)) + 
    geom_linerange(position = position_dodge(width = 0.8)) + theme_pubr()+
    geom_pointrange(position = position_dodge(width = 0.8)) + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    coord_flip() +scale_color_manual(values=c("gray70","gray10"))
  return(p)
}
###########################################################
#subset data
dat_nn<-subset_dat(notnorm_data)
write.table(dat_nn,paste0(dat_path,"species_pfam_metaT/species_pfam_clean_noNT.tsv"),sep = "\t",row.names = FALSE, quote=FALSE)

dat_n<-subset_dat(norm_data)
write.table(dat_n,paste0(dat_path,"species_pfam_metaT/species_pfam_clean_noNT-TPM.tsv"),sep = "\t",row.names = FALSE, quote=FALSE)
###########################################################
#subset the data to just have BSH by light and dark-->to run birdman

bsh<-fread(paste0(dat_path,"species_pfam_metaT/species_pfam_clean_noNT.tsv"))%>%
  dplyr::filter(grepl("PF02275.21",FeatureID))%>%
  separate(FeatureID,c("FeatureID",NA), sep="\\|PF", extra="drop")%>%
  mutate(FeatureID=gsub(" ","_",FeatureID))

bsh<-bsh%>%column_to_rownames("FeatureID")
bsh_rmz <- bsh[!(rowSums(bsh != 0) ==0), ]
bsh_rmz<-bsh_rmz%>%rownames_to_column("FeatureID")

#just light
dat<-subset_dat_bshLD(bsh_rmz,"light")
write.table(dat,paste0(dat_path,"species_pfam_metaT/species_pfam_BSHonly_clean_rmzero_noNT_light.tsv"),sep = "\t",row.names = FALSE, quote=FALSE)
#just dark
dat<-subset_dat_bshLD(bsh_rmz,"dark")
write.table(dat,paste0(dat_path,"species_pfam_metaT/species_pfam_BSHonly_clean_rmzero_noNT_dark.tsv"),sep = "\t",row.names = FALSE, quote=FALSE)
###########################################################
#subset the data to just have BSH and RPOB-->to run qurro

bsh<-fread(paste0(dat_path,"species_pfam_metaT/species_pfam_clean_noNT.tsv"))%>%
  dplyr::filter(grepl("PF02275.21",FeatureID)|grepl("PF04563.18",FeatureID))

bsh<-bsh%>%column_to_rownames("FeatureID")
bsh_rmz <- bsh[!(rowSums(bsh != 0) ==0), ]
bsh_rmz<-bsh_rmz%>%rownames_to_column("FeatureID") #67 BSH 1085 rpob
write.table(bsh_rmz,paste0(dat_path,"species_pfam_metaT/species_pfam_BSH_RPOB_clean_rmzero_noNT.tsv"),sep = "\t",row.names = FALSE, quote=FALSE)
###########################################################
#natural log of BSH vs RPOB--Figure 4B

md<-fread(dat_metadata)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  dplyr::rename(`Sample ID`=sample_name)

natlog_bsh<-fread(paste0(dat_path,"species_pfam_metaT/rpca_results_BSH_RPOB_rmzero/sample_plot_data_speciespfam_BSH_RPOB.tsv"))%>%
  dplyr::select(1:2)%>%
  left_join(.,md,by="Sample ID")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  mutate(condition=factor(condition,levels=c("NA","FA","FT")),
         phase=factor(phase,levels=c("light","dark")),
         Current_Natural_Log_Ratio=as.numeric(Current_Natural_Log_Ratio))%>%
  filter(phase=="light")

p<-ggplot(natlog_bsh, aes(x=condition, y=Current_Natural_Log_Ratio,fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  facet_wrap(~phase)+
  theme_classic()+ scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
  labs(x="condition",y="Natural Log Ratio (BSH/RPOB)",title="Bile Salt Hydrolase (BSH)")+
  theme(legend.position = "none")

ggsave(paste0(fig_path,"SFR23_0705_natlog_speciesBSHvsRPOB_lightdark.pdf"), plot=p,height=3, width=4)

natlog_bshL<-natlog_bsh%>%filter(phase=="light")
pairwise.wilcox.test(natlog_bshL$Current_Natural_Log_Ratio, natlog_bshL$condition,
                     p.adjust.method="fdr")
#species BSH_RPOB
# NA     FA    
# FA 0.0012 -     
#   FT 0.0012 0.5887

natlog_bshD<-natlog_bsh%>%filter(phase=="dark")
pairwise.wilcox.test(natlog_bshD$Current_Natural_Log_Ratio, natlog_bshD$condition,
                     p.adjust.method="fdr")
#species BSH_RPOB
# NA     FA    
# FA 0.0006 -     
#   FT 0.0006 0.0649

###########################################################
#plot the distribution of BSH by species--Figure S4B

mdT<-fread(dat_metadata)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))

selfeat<-fread(paste0(dat_path,"species_pfam_metaT/rpca_results_BSH_RPOB_rmzero/selected_features_speciespfam_BSH_RPOB.tsv"))%>%
  dplyr::filter(grepl("PF02275.21",`Feature ID`))

bsh_datL<-bsh_dist_dat(paste0(dat_path,"species_pfam_metaT/species_pfam_clean_noNT-TPM.tsv"),"light")
lst_keep<-((bsh_datL%>%arrange(-sum_log_TPM))$FeatureID%>%unique())[1:15]
bsh_datL <- bsh_datL %>%filter(FeatureID %in% lst_keep)

p<-bsh_dist_plt(bsh_datL,"Light")
ggsave(paste0(fig_path,"SFR23_0719_BSHlight_logTPM.pdf"), p, height=3, width=3.5)

bsh_datD<-bsh_dist_dat(paste0(dat_path,"species_pfam_metaT/species_pfam_clean_noNT-TPM.tsv"),"dark")
lst_keep<-((bsh_datD%>%arrange(-sum_log_TPM))$FeatureID%>%unique())[1:15]
bsh_datD <- bsh_datD %>%filter(FeatureID %in% lst_keep)

p<-bsh_dist_plt(bsh_datD,"Dark")
ggsave(paste0(fig_path,"SFR23_0719_BSHdark_logTPM.pdf"), p, height=3, width=3.5)
###########################################################
#birdman results (NA vs. FA)--Figure S4C-D

#light arranged hdi<10
BSHlight<-bdm_nafa_dat(bdm_NAFA_light,"light")%>%
  filter(diff<10)%>%
  arrange(ratio)

BSHdark<-bdm_nafa_dat(bdm_NAFA_dark,"dark")%>%
  filter(FeatureID %in% BSHlight$Feature)

BSH<-rbind(BSHlight,BSHdark)%>%mutate(phase=factor(phase,levels=c("light","dark")))
BSH$Feature <- factor(BSH$Feature,levels = BSHlight$Feature)

p<-bdm_plt(BSH)
ggsave(paste0(fig_path,"FANALD_smhdi.pdf"),p,height=5, width=6)

#dark arranged hdi<10
BSHdark<-bdm_nafa_dat(bdm_NAFA_dark,"dark")%>%
  filter(diff<10)%>%
  arrange(ratio)

BSHlight<-bdm_nafa_dat(bdm_NAFA_light,"light")%>%
  filter(FeatureID %in% BSHdark$Feature)

BSH<-rbind(BSHlight,BSHdark)%>%mutate(phase=factor(phase,levels=c("light","dark")))
BSH$Feature <- factor(BSH$Feature,levels = BSHdark$Feature)

p<-bdm_plt(BSH)
ggsave(paste0(fig_path,"FANALD_darr_smhdi.pdf"),p,height=5, width=6)

###########################################################
#birdman results (FA vs FT)--Figure 4D

BSHlight<-bdm_faft_dat(bdm_FAFT_light,"light")%>%
  filter(diff<10)

BSHdark<-bdm_faft_dat(bdm_FAFT_dark,"dark")%>%
  filter(FeatureID %in% BSHlight$FeatureID)%>%
  arrange(ratio)

BSH<-rbind(BSHlight,BSHdark)%>%mutate(phase=factor(phase,levels=c("light","dark")))

BSH$Feature <- factor(BSH$Feature,levels = BSHdark$Feature)

p<-bdm_plt(BSH)
ggsave(paste0(fig_path,"FAFTLD_darr_smhdiL.pdf"),p,height=5, width=6)

###########################################################
#summarized Birdman results table--Table S4
#need to run the targetted bsh search script first

BSHlight<-bdm_faft_dat(bdm_FAFT_light,"light")
BSHdark<-bdm_faft_dat(bdm_FAFT_dark,"dark")
BSH<-rbind(BSHlight,BSHdark)%>%mutate(phase=factor(phase,levels=c("light","dark")))%>%dplyr::select(-diff)
uFAFT<-BSH%>%mutate(comparison= "FAvFT",search_method="untargetted")%>%
  dplyr::rename(BSHspecies=FeatureID)

BSHlight<-bdm_nafa_dat(bdm_NAFA_light,"light")
BSHdark<-bdm_nafa_dat(bdm_NAFA_dark,"dark")
BSH<-rbind(BSHlight,BSHdark)%>%mutate(phase=factor(phase,levels=c("light","dark")))%>%dplyr::select(-diff)
uNAFA<-BSH%>%mutate(comparison= "NAvFA",search_method="untargetted")%>%
  dplyr::rename(BSHspecies=FeatureID)

tFAFT<-fread(paste0(dat_path,"BSH_proteindb/birdman_outputs/genomeLD_noNT_rmdbton.beta_var.tsv"))%>%
  mutate(comparison= "FAvFT",search_method="targetted")%>%
  dplyr::rename(BSHspecies=name)

bsh_birdman<-rbind(uFAFT,uNAFA,tFAFT)%>%
  arrange(ratio)

write.table(bsh_birdman, file = paste0(dat_path,"birdman_results_combined.txt"),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
