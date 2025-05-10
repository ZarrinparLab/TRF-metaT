setwd("~/Notebooks/sfloresr/TRF-metaT/")

library(tidyverse)
library(data.table)
library(ggpubr)

###########################################################
#paths
notnorm_data<-"data/bsh_analysis/species_pfam_metaG/species_pfam.tsv"
dat_metadata<-"data/bsh_analysis/metaG_metadata_noNT.txt"
dat_path<-"data/bsh_analysis/"
fig_path<-"figures/bsh_analysis/"
###########################################################
#functions
subset_dat<-function(dt){
  samps_list<-c("#FeatureID","FA1b_S36","FA1c_S37","FA2b_S38","FA2c_S39","FA3b_S40","FA3c_S41","FA4b_S42","FA4c_S43",  
                "FA5b_S44","FA5c_S45","FA6b_S46","FA6c_S47","FT1b_S48","FT1c_S49","FT2b_S50","FT2c_S51","FT3b_S52",
                "FT3c_S53","FT4b_S54","FT4c_S55","FT5b_S56","FT5c_S57","FT6b_S58","FT6c_S59","NA1a_S1","NA1b_S2",  
                "NA1c_S3","NA2a_S4","NA2b_S5","NA2c_S6","NA3a_S7","NA3b_S8","NA3c_S9","NA4a_S10","NA4b_S11",  
                "NA4c_S12","NA5a_S13","NA5b_S14","NA5c_S15","NA6a_S16","NA6b_S17","NA6c_S18")
  dat<-fread(dt)%>%dplyr::select(all_of(samps_list))
  names(dat) <- c("FeatureID","FA1b","FA1c","FA2b","FA2c","FA3b","FA3c","FA4b","FA4c",  
                  "FA5b","FA5c","FA6b","FA6c","FT1b","FT1c","FT2b","FT2c","FT3b",
                  "FT3c","FT4b","FT4c","FT5b","FT5c","FT6b","FT6c","NA1a","NA1b",  
                  "NA1c","NA2a","NA2b","NA2c","NA3a","NA3b","NA3c","NA4a","NA4b",  
                  "NA4c","NA5a","NA5b","NA5c","NA6a","NA6b","NA6c")
  return(dat)
}

subset_dat_bshLD<-function(dt,phase){
  if(phase=="light"){
    samps_list<-c("FeatureID","FA4b","FA4c","FA5b","FA5c","FA6b","FA6c",
                  "FT4b","FT4c","FT5b","FT5c","FT6b","FT6c","NA4a","NA4b",  
                  "NA4c","NA5a","NA5b","NA5c","NA6a","NA6b","NA6c")
  }
  else{
    samps_list<-c("FeatureID","FA1b","FA1c","FA2b","FA2c","FA3b","FA3c",
                  "FT1b","FT1c","FT2b","FT2c","FT3b","FT3c","NA1a","NA1b",  
                  "NA1c","NA2a","NA2b","NA2c","NA3a","NA3b","NA3c")
  }
  dat<-dt%>%dplyr::select(all_of(samps_list))
  return(dat)
}
###########################################################
#subset data
dat_nn<-subset_dat(notnorm_data)
write.table(dat_nn,paste0(dat_path,"species_pfam_metaG/species_pfam_clean_noNT.tsv"),sep = "\t",row.names = FALSE, quote=FALSE)

###########################################################
#subset the data to just have BSH by light and dark-->to run birdman

bsh<-fread(paste0(dat_path,"species_pfam_metaG/species_pfam_clean_noNT.tsv"))%>%
  dplyr::filter(grepl("PF02275.21",FeatureID))%>%
  separate(FeatureID,c("FeatureID",NA), sep="\\|PF", extra="drop")%>%
  mutate(FeatureID=gsub(" ","_",FeatureID))

bsh<-bsh%>%column_to_rownames("FeatureID")
bsh_rmz <- bsh[!(rowSums(bsh != 0) ==0), ]
bsh_rmz<-bsh_rmz%>%rownames_to_column("FeatureID")

#just light
dat<-subset_dat_bshLD(bsh_rmz,"light")
write.table(dat,paste0(dat_path,"species_pfam_metaG/species_pfam_BSHonly_clean_rmzero_noNT_light.tsv"),sep = "\t",row.names = FALSE, quote=FALSE)

#just dark
dat<-subset_dat_bshLD(bsh_rmz,"dark")
write.table(dat,paste0(dat_path,"species_pfam_metaG/species_pfam_BSHonly_clean_rmzero_noNT_dark.tsv"),sep = "\t",row.names = FALSE, quote=FALSE)
###########################################################
#subset the data to just have BSH and RPOB-->to run qurro

bsh<-fread(paste0(dat_path,"species_pfam_metaG/species_pfam_clean_noNT.tsv"))%>%
  dplyr::filter(grepl("PF02275.21",FeatureID)|grepl("PF04563.18",FeatureID))

bsh<-bsh%>%column_to_rownames("FeatureID")
bsh_rmz <- bsh[!(rowSums(bsh != 0) ==0), ]
bsh_rmz<-bsh_rmz%>%rownames_to_column("FeatureID")
write.table(bsh_rmz,paste0(dat_path,"species_pfam_metaG/species_pfam_BSH_RPOB_clean_rmzero_noNT.tsv"),sep = "\t",row.names = FALSE, quote=FALSE)
###########################################################
#natural log of BSH vs RPOB--Figure S4E

md<-fread(dat_metadata)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  dplyr::rename(phase=lightdark,`Sample ID`=sample_name)

natlog_bsh<-fread(paste0(dat_path,"species_pfam_metaG/rpca_results_BSH_RPOB_rmzero/sample_plot_data_speciespfam_BSH_RPOB.tsv"))%>%
  dplyr::select(1:2)%>%
  left_join(.,md,by="Sample ID")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  mutate(condition=factor(condition,levels=c("NA","FA","FT")),
         phase=factor(phase,levels=c("light","dark")),
         Current_Natural_Log_Ratio=as.numeric(Current_Natural_Log_Ratio))

p<-ggplot(natlog_bsh, aes(x=condition, y=Current_Natural_Log_Ratio,fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  facet_wrap(~phase)+
  theme_classic()+ scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
  labs(x="condition",y="Natural Log Ratio (BSH/RPOB)",title="Bile Salt Hydrolase (BSH)")+
  theme(legend.position = "none")

ggsave(paste0(fig_path,"SFR25_0108_natlog_speciesBSHvsRPOB_lightdark_metaG.pdf"), plot=p,height=3, width=4)
natlog_bshL<-natlog_bsh%>%filter(phase=="light")
pairwise.wilcox.test(natlog_bshL$Current_Natural_Log_Ratio, natlog_bshL$condition,
                     p.adjust.method="fdr")

#species BSH_RPOB
# NA   FA  
# FA 0.37 -   
#   FT 0.64 0.37

natlog_bshD<-natlog_bsh%>%filter(phase=="dark")
pairwise.wilcox.test(natlog_bshD$Current_Natural_Log_Ratio, natlog_bshD$condition,
                     p.adjust.method="fdr")
#species BSH_RPOB
# NA   FA  
# FA 0.95 -   
#   FT 0.95 0.95

###########################################################
 
#note birdman run kept failing bc of all the zeroes for NA vs FA and FA vs FT overall and light/dark
