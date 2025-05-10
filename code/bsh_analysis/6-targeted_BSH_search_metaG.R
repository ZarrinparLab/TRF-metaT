setwd("~/Notebooks/sfloresr/TRF-metaT/")

library(tidyverse)
library(data.table)
library(ggpubr)

##############################################################################
#paths
notnorm_data<-"data/bsh_analysis/BSH_proteindb_metaG/genome.tsv"
norm_data<-"data/bsh_analysis/BSH_proteindb_metaG/genome-TPM.tsv"
bdm_FAFT_light<-"data/bsh_analysis/BSH_proteindb_metaG/birdman_outputs/genomeL_noNT_rmdbton.beta_var.tsv"
bdm_FAFT_dark<-"data/bsh_analysis/BSH_proteindb_metaG/birdman_outputs/genomeD_noNT_rmdbton.beta_var.tsv"
dat_metadata<-"data/bsh_analysis/metaG_metadata_noNT.txt"
prot_metadata<-"data/bsh_analysis/BSH_db_metadata.txt"
dat_path<-"data/bsh_analysis/"
fig_path<-"figures/bsh_analysis/"
##############################################################################
#functions
subset_dat<-function(dt){
  samps_list<-c("#FeatureID","FA1b_S36","FA1c_S37","FA2b_S38","FA2c_S39","FA3b_S40","FA3c_S41","FA4b_S42","FA4c_S43",  
                "FA5b_S44","FA5c_S45","FA6b_S46","FA6c_S47","FT1b_S48","FT1c_S49","FT2b_S50","FT2c_S51","FT3b_S52",
                "FT3c_S53","FT4b_S54","FT4c_S55","FT5b_S56","FT5c_S57","FT6b_S58","FT6c_S59","NA1a_S1","NA1b_S2",  
                "NA1c_S3","NA2a_S4","NA2b_S5","NA2c_S6","NA3a_S7","NA3b_S8","NA3c_S9","NA4a_S10","NA4b_S11",  
                "NA4c_S12","NA5a_S13","NA5b_S14","NA5c_S15","NA6a_S16","NA6b_S17","NA6c_S18")
  dat<-fread(dt)%>%dplyr::select(all_of(samps_list))%>%
    mutate(`#FeatureID`=gsub("\\|","_",`#FeatureID`))
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
  dat<-dt%>%dplyr::select(all_of(samps_list))%>%column_to_rownames("FeatureID")
  dat<- dat[!(rowSums(dat != 0) <2), ]
  dat<-dat%>%rownames_to_column("FeatureID")
  
  return(dat)
}

bdm_faft_dat<-function(dat,phase){
  bsh_bdmdat<-fread(dat)%>%
    dplyr::rename(ratio=`C(condition, Treatment('FA'))[T.FT]_mean`,
                  FeatureID=Feature)%>%
    mutate(`C(condition, Treatment('FA'))[T.FT]_hdi`=gsub("[(]|[)]","",`C(condition, Treatment('FA'))[T.FT]_hdi`))%>%
    separate(`C(condition, Treatment('FA'))[T.FT]_hdi`,c("min","max"), sep=",")%>%
    mutate(min=as.numeric(min),
           max=as.numeric(max))%>%
    left_join(.,pmd,by="FeatureID")%>%
    dplyr::select(name, ratio, min, max)
  
  if(phase=="light"){
    bsh_bdmdat<-bsh_bdmdat%>%mutate(phase="light")
  }
  else{
    bsh_bdmdat<-bsh_bdmdat%>%mutate(phase="dark")
  }
  return(bsh_bdmdat)
}
##############################################################################
#subset data-->to run qurro
dat_nn<-subset_dat(notnorm_data)
write.table(dat_nn,paste0(dat_path,"BSH_proteindb_metaG/genome_noNT.tsv"),sep = "\t",row.names = FALSE, quote=FALSE)

dat<-dat_nn%>%column_to_rownames("FeatureID")
dat_rmzs <- dat[!(rowSums(dat != 0) <2), ]
dat_rmzs<-dat_rmzs%>%rownames_to_column("FeatureID")
write.table(dat_rmzs,paste0(dat_path,"BSH_proteindb_metaG/genome_noNT_rmdbton.tsv"),sep = "\t",row.names = FALSE, quote=FALSE)

##############################################################################
#subset by light and dark-->to run birdman

#just light
dat_rmzs_l<-subset_dat_bshLD(dat_rmzs,"light")
write.table(dat_rmzs_l,paste0(dat_path,"BSH_proteindb_metaG/genomeL_noNT_rmdbton.tsv"),sep = "\t",row.names = FALSE, quote=FALSE)

#just dark
dat_rmzs_d<-subset_dat_bshLD(dat_rmzs,"dark")
write.table(dat_rmzs_d,paste0(dat_path,"BSH_proteindb_metaG/genomeD_noNT_rmdbton.tsv"),sep = "\t",row.names = FALSE, quote=FALSE)

##############################################################################
#natural log of top and bottom 20% of bsh--Figure S4F

md<-fread(dat_metadata)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  dplyr::rename(phase=lightdark)

natlog_bsh<-fread(paste0(dat_path,"BSH_proteindb_metaG/rpca_results_genome/sample_plot_data_20tp.tsv"))%>%
  dplyr::select(1:2)%>%
  dplyr::rename(sample_name=`Sample ID`)%>%
  left_join(.,md,by="sample_name")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  mutate(condition=factor(condition,levels=c("NA","FA","FT")),
         phase=factor(phase,levels=c("light","dark")),
         Current_Natural_Log_Ratio=as.numeric(Current_Natural_Log_Ratio))

p<-ggplot(natlog_bsh, aes(x=condition, y=Current_Natural_Log_Ratio,fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  facet_wrap(~phase)+
  theme_classic()+ scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
  labs(x="condition",y="Natural Log Ratio (top20/bottom20 BSH)",title="Bile Salt Hydrolase (BSH)")+
  theme(legend.position = "none")

ggsave(paste0(fig_path,"nat_log_20tp_metaG.pdf"),height=3, width=4)

natlog_bshL<-natlog_bsh%>%filter(phase=="light")
pairwise.wilcox.test(natlog_bshL$Current_Natural_Log_Ratio, natlog_bshL$condition,
                     p.adjust.method="fdr")

# NA  FA 
# FA 0.2 -  
#   FT 0.4 0.2


natlog_bshD<-natlog_bsh%>%filter(phase=="dark")
pairwise.wilcox.test(natlog_bshD$Current_Natural_Log_Ratio, natlog_bshD$condition,
                     p.adjust.method="fdr")

# NA   FA  
# FA 0.61 -   
#   FT 0.57 0.57
##############################################################################
#birdman results (FA vs FT)--Figure S4G

pmd<-fread(prot_metadata)%>%
  mutate(FeatureID=paste("tr_",Entry,"_",`Entry Name`,sep=""),
         name=paste(Organism,Entry,sep="-"))%>%
  dplyr::select(FeatureID, everything())

BSHlight<-bdm_faft_dat(bdm_FAFT_light,"light")
BSHdark<-bdm_faft_dat(bdm_FAFT_dark,"dark")%>%
  filter(name %in% BSHlight$name)%>%
  arrange(ratio)
BSHlight<-bdm_faft_dat(bdm_FAFT_light,"light")%>%
  filter(name %in% BSHdark$name)

BSH<-rbind(BSHlight,BSHdark)%>%
  mutate(phase=factor(phase,levels=c("light","dark")))

BSH$name <- factor(BSH$name,levels = BSHdark$name )

ggplot(BSH, aes(x =name , y = ratio, ymin = min, ymax = max, color=phase, group=phase)) + 
  geom_linerange(position = position_dodge(width = 0.8)) + theme_pubr()+
  geom_pointrange(position = position_dodge(width = 0.8)) + 
  scale_y_continuous(limits=c(-12.5,10))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  coord_flip() +scale_color_manual(values=c("gray70","gray10"))

ggsave(paste0(fig_path,"FAFTLD_darkarrange_metaG.pdf"),height=7, width=14)
