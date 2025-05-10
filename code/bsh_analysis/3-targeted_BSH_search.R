setwd("~/Notebooks/sfloresr/TRF-metaT/")

library(tidyverse)
library(data.table)
library(ggpubr)

##############################################################################
#paths
notnorm_data<-"data/bsh_analysis/BSH_proteindb/genome.tsv"
norm_data<-"data/bsh_analysis/BSH_proteindb/genome-TPM.tsv"
bdm_FAFT_light<-"data/bsh_analysis/BSH_proteindb/birdman_outputs/genomeL_noNT_rmdbton.beta_var.tsv"
bdm_FAFT_dark<-"data/bsh_analysis/BSH_proteindb/birdman_outputs/genomeD_noNT_rmdbton.beta_var.tsv"
dat_metadata<-"data/bsh_analysis/metaT_metadata_ztcat_noNT.txt"
prot_metadata<-"data/bsh_analysis/BSH_db_metadata.txt"
dat_path<-"data/bsh_analysis/"
fig_path<-"figures/bsh_analysis/"
##############################################################################
#functions
subset_dat<-function(dt){
  samps_list<-c("#FeatureID","cFA01a_S39","cFA01b_S45","cFA05a_S40","cFA05b_S46","cFA09a_S41",
                "cFA09b_S47","cFA13a_S36","cFA13b_S42","cFA17a_S37","cFA17b_S43",
                "cFA21a_S38","cFA21b_S44","cFT01a_S51","cFT01b_S57","cFT05a_S52",
                "cFT05b_S58","cFT09a_S53","cFT09b_S59","cFT13a_S48","cFT13b_S54",
                "cFT17a_S49","cFT17b_S55","cFT21a_S50","cFT21b_S56","cNA01a_S4",
                "cNA01b_S10","cNA01c_S16","cNA05a_S5","cNA05b_S11","cNA05c_S17",
                "cNA09a_S6","cNA09b_S12","cNA09c_S18","cNA13a_S1","cNA13b_S7",
                "cNA13c_S13","cNA17a_S2","cNA17b_S8","cNA17c_S14","cNA21a_S3",
                "cNA21b_S9","cNA21c_S15")
  dat<-fread(dt)%>%dplyr::select(all_of(samps_list))%>%
    mutate(`#FeatureID`=gsub("\\|","_",`#FeatureID`))
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
write.table(dat_nn,paste0(dat_path,"BSH_proteindb/genome_noNT.tsv"),sep = "\t",row.names = FALSE, quote=FALSE)

dat<-dat_nn%>%column_to_rownames("FeatureID")
dat_rmzs <- dat[!(rowSums(dat != 0) <2), ]
dat_rmzs<-dat_rmzs%>%rownames_to_column("FeatureID")
write.table(dat_rmzs,paste0(dat_path,"BSH_proteindb/genome_noNT_rmdbton.tsv"),sep = "\t",row.names = FALSE, quote=FALSE)

dat_n<-subset_dat(norm_data)
write.table(dat_n,paste0(dat_path,"BSH_proteindb/genome-TPM_noNT.tsv"),sep = "\t",row.names = FALSE, quote=FALSE)
###########################################################
#subset by light and dark-->to run birdman

#just light
dat_rmzs_l<-subset_dat_bshLD(dat_rmzs,"light")
write.table(dat_rmzs_l,paste0("BSH_proteindb/genomeL_noNT_rmdbton.tsv"),sep = "\t",row.names = FALSE, quote=FALSE)

#just dark
dat_rmzs_d<-subset_dat_bshLD(dat_rmzs,"dark")
write.table(dat_rmzs_d,"BSH_proteindb/genomeD_noNT_rmdbton.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

##############################################################################
#natural log of top and bottom 10% of bsh--Figure 4C

md<-fread(dat_metadata)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))

natlog_bsh<-fread(paste0(dat_path,"BSH_proteindb/rpca_results_genome/sample_plot_data_10tp.tsv"))%>%
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
  labs(x="condition",y="Natural Log Ratio (top10/bottom10 BSH)",title="Bile Salt Hydrolase (BSH)")+
  theme(legend.position = "none")

ggsave(paste0(fig_path,"nat_log_10tp.pdf"),height=3, width=4)

natlog_bshL<-natlog_bsh%>%filter(phase=="light")
pairwise.wilcox.test(natlog_bshL$Current_Natural_Log_Ratio, natlog_bshL$condition,
                     p.adjust.method="fdr")

# NA     FA    
# FA 0.0006 -     
#   FT 0.0006 0.5887

natlog_bshD<-natlog_bsh%>%filter(phase=="dark")
pairwise.wilcox.test(natlog_bshD$Current_Natural_Log_Ratio, natlog_bshD$condition,
                     p.adjust.method="fdr")
# NA     FA    
# FA 0.0012 -     
#   FT 0.0012 0.1797

##############################################################################
#birdman results (FA vs FT)--Figure 4E

pmd<-fread(prot_metadata)%>%
  mutate(FeatureID=paste("tr_",Entry,"_",`Entry Name`,sep=""),
         name=paste(Organism,Entry,sep="-"))%>%
  dplyr::select(FeatureID, everything())

BSHdark<-bdm_faft_dat(bdm_FAFT_dark,"dark")
BSHlight<-bdm_faft_dat(bdm_FAFT_light,"light")
BSH<-rbind(BSHlight,BSHdark)%>%mutate(phase=factor(phase,levels=c("light","dark")))
write.table(BSH,paste0(dat_path,"BSH_proteindb/birdman_outputs/genomeLD_noNT_rmdbton.beta_var.tsv"),sep = "\t",row.names = FALSE, quote=FALSE)

BSHdark<-bdm_faft_dat(bdm_FAFT_dark,"dark")%>%
  mutate(credible=ifelse(min>0|max<0,"yes","no"))%>%
  filter(credible=="yes")%>%
  dplyr::select(name, ratio, min, max, phase)%>%
  arrange(ratio)
                
BSHlight<-bdm_faft_dat(bdm_FAFT_light,"light")%>%
  filter(name %in% BSHdark$name)

selfeat<-rbind(BSHlight,BSHdark)%>%select(1:2,5)%>%
  spread(phase,ratio)%>%
  mutate(diff_m=abs(light-dark))%>%
  filter(diff_m>3)%>%
  arrange(desc(diff_m))

BSH<-rbind(BSHlight,BSHdark)%>%
  mutate(phase=factor(phase,levels=c("light","dark")))%>%
  filter(name %in% selfeat$name)

BSH$name <- factor(BSH$name,levels = BSHdark$name )

ggplot(BSH, aes(x =name , y = ratio, ymin = min, ymax = max, color=phase, group=phase)) + 
  geom_linerange(position = position_dodge(width = 0.8)) + theme_pubr()+
  geom_pointrange(position = position_dodge(width = 0.8)) + 
  scale_y_continuous(limits=c(-12.5,10))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  coord_flip() +scale_color_manual(values=c("gray70","gray10"))

ggsave(paste0(fig_path,"FAFTLD_Dcredbgdiff_smhdi_darkarranged.pdf"),height=7, width=14.5)
