setwd("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metagenomic/woltka2_results/filtered_metaG/")

library(tidyverse)
library(data.table)
library("qiime2R")
library("Biostrings")
library(ggrepel)
library(ggpubr)
library(gplots)
library(ggvenn)
library(ggbreak)
library(viridis)
###########################################################
samps_list<-c("#FeatureID","FA1b_S36","FA1c_S37","FA2b_S38","FA2c_S39","FA3b_S40","FA3c_S41","FA4b_S42","FA4c_S43",  
              "FA5b_S44","FA5c_S45","FA6b_S46","FA6c_S47","FT1b_S48","FT1c_S49","FT2b_S50","FT2c_S51","FT3b_S52",
              "FT3c_S53","FT4b_S54","FT4c_S55","FT5b_S56","FT5c_S57","FT6b_S58","FT6c_S59","NA1a_S1","NA1b_S2",  
              "NA1c_S3","NA2a_S4","NA2b_S5","NA2c_S6","NA3a_S7","NA3b_S8","NA3c_S9","NA4a_S10","NA4b_S11",  
              "NA4c_S12","NA5a_S13","NA5b_S14","NA5c_S15","NA6a_S16","NA6b_S17","NA6c_S18")

dat<-fread("species_pfam/species_pfam.tsv")%>%dplyr::select(all_of(samps_list))
names(dat) <- c("FeatureID","FA1b","FA1c","FA2b","FA2c","FA3b","FA3c","FA4b","FA4c",  
                "FA5b","FA5c","FA6b","FA6c","FT1b","FT1c","FT2b","FT2c","FT3b",
                "FT3c","FT4b","FT4c","FT5b","FT5c","FT6b","FT6c","NA1a","NA1b",  
                "NA1c","NA2a","NA2b","NA2c","NA3a","NA3b","NA3c","NA4a","NA4b",  
                "NA4c","NA5a","NA5b","NA5c","NA6a","NA6b","NA6c")
write.table(dat,"species_pfam/species_pfam_clean_noNT.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

#remove singletons and doubletons
dat<-dat%>%column_to_rownames("FeatureID")

dat_rmz <- dat[!(rowSums(dat != 0) == 0), ]
dat_rmz<-dat_rmz%>%rownames_to_column("FeatureID")
write.table(dat_rmz,"species_pfam/species_pfam_clean_rmzero_noNT.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

dat_rmzs <- dat[!(rowSums(dat != 0) <2), ]
dat_rmzs<-dat_rmzs%>%rownames_to_column("FeatureID")
write.table(dat_rmzs,"species_pfam/species_pfam_clean_rmsingletons_noNT.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

dat_rmzd <- dat[!(rowSums(dat != 0) < 3), ]
dat_rmzd<-dat_rmzd%>%rownames_to_column("FeatureID")
write.table(dat_rmzd,"species_pfam/species_pfam_clean_rmdoubletons_noNT.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

#just light
samps_list<-c("FeatureID","FA4b","FA4c","FA5b","FA5c","FA6b","FA6c",
              "FT4b","FT4c","FT5b","FT5c","FT6b","FT6c","NA4a","NA4b",  
              "NA4c","NA5a","NA5b","NA5c","NA6a","NA6b","NA6c")

dat<-fread("species_pfam/species_pfam_clean_rmdoubletons_noNT.tsv")%>%dplyr::select(all_of(samps_list))
write.table(dat,"species_pfam/species_pfam_clean_rmdoubletons_noNT_light.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

#just dark
samps_list<-c("FeatureID","FA1b","FA1c","FA2b","FA2c","FA3b","FA3c",
              "FT1b","FT1c","FT2b","FT2c","FT3b","FT3c","NA1a","NA1b",  
              "NA1c","NA2a","NA2b","NA2c","NA3a","NA3b","NA3c")

dat<-fread("species_pfam/species_pfam_clean_rmdoubletons_noNT.tsv")%>%dplyr::select(all_of(samps_list))
write.table(dat,"species_pfam/species_pfam_clean_rmdoubletons_noNT_dark.tsv",sep = "\t",row.names = FALSE, quote=FALSE)
###########################################################

#subset the data to just have BSH

bsh<-fread("species_pfam/species_pfam_clean_noNT.tsv")%>%
  dplyr::filter(grepl("PF02275.21",FeatureID))%>%
  separate(FeatureID,c("FeatureID",NA), sep="\\|PF", extra="drop")%>%
  mutate(FeatureID=gsub(" ","_",FeatureID))

bsh<-bsh%>%column_to_rownames("FeatureID")
bsh_rmz <- bsh[!(rowSums(bsh != 0) ==0), ]
bsh_rmz<-bsh_rmz%>%rownames_to_column("FeatureID")
write.table(bsh_rmz,"species_pfam/species_pfam_BSHonly_clean_rmzero_noNT.tsv",sep = "\t",row.names = FALSE, quote=FALSE)
#just light
samps_list<-c("FeatureID","FA4b","FA4c","FA5b","FA5c","FA6b","FA6c",
              "FT4b","FT4c","FT5b","FT5c","FT6b","FT6c","NA4a","NA4b",  
              "NA4c","NA5a","NA5b","NA5c","NA6a","NA6b","NA6c")

dat<-fread("species_pfam/species_pfam_BSHonly_clean_rmzero_noNT.tsv")%>%dplyr::select(all_of(samps_list))
dat_rmz <- dat[!(rowSums(dat != 0) <2), ]
write.table(dat,"species_pfam/species_pfam_BSHonly_clean_rmzero_noNT_light.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

#just dark
samps_list<-c("FeatureID","FA1b","FA1c","FA2b","FA2c","FA3b","FA3c",
              "FT1b","FT1c","FT2b","FT2c","FT3b","FT3c","NA1a","NA1b",  
              "NA1c","NA2a","NA2b","NA2c","NA3a","NA3b","NA3c")

dat<-fread("species_pfam/species_pfam_BSHonly_clean_rmzero_noNT.tsv")%>%dplyr::select(all_of(samps_list))
dat_rmz <- dat[!(rowSums(dat != 0) <2), ]
write.table(dat,"species_pfam/species_pfam_BSHonly_clean_rmzero_noNT_dark.tsv",sep = "\t",row.names = FALSE, quote=FALSE)
###########################################################

#subset the data to just have BSH and RPOB
# bsh<-fread("genus_pfam/genus_pfam_clean_noNT.tsv")%>%
#   dplyr::filter(grepl("PF02275.21",FeatureID)|grepl("PF04563.18",FeatureID))

bsh<-fread("species_pfam/species_pfam_clean_noNT.tsv")%>%
  dplyr::filter(grepl("PF02275.21",FeatureID)|grepl("PF04563.18",FeatureID))

bsh<-bsh%>%column_to_rownames("FeatureID")
bsh_rmz <- bsh[!(rowSums(bsh != 0) ==0), ]
bsh_rmz<-bsh_rmz%>%rownames_to_column("FeatureID") #67 BSH 1085 rpob
write.table(bsh_rmz,"species_pfam/species_pfam_BSH_RPOB_clean_rmzero_noNT.tsv",sep = "\t",row.names = FALSE, quote=FALSE)
###########################################################

##look at RPCA

md<-fread("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metagenomic/metaG_metadata_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  dplyr::rename(phase=lightdark)

ord <- read_qza("species_pfam/rpca_results_rmdoubletons/ordination.qza")

samp_ord<-ord$data$Vectors
write.table(samp_ord,"species_pfam/rpca_results_rmdoubletons/sample_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)
feat_ord<-ord$data$Species
write.table(feat_ord,"species_pfam/rpca_results_rmdoubletons/feature_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)

rpca<-ord$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2)%>%
  dplyr::rename(sample_name=SampleID)%>%
  left_join(md,by="sample_name")%>%
  mutate(condition=factor(condition,levels=c("NA","FA","FT")),
         phase=factor(phase,levels=c("light","dark")))

p<-rpca %>%
  ggplot(aes(x=PC1, y=PC2, color=condition, shape=phase)) +
  geom_point(alpha=1.0) + 
  theme_pubr() +#stat_ellipse(type = "t", linetype = 2,aes(group = condition))+
  scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+
  scale_shape_manual(values=c(3,16)) +
  labs(color="condition",
       x =paste("PC1 (",round(ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""))+ggtitle("genus|pfam RPCA")+ theme(plot.title = element_text(face = "bold"))
ggsave("species_pfam/rpca_results_rmdoubletons/SFR25_0108_speciespfam_RPCA.pdf", plot=p,height=3.5, width=3.5)

###########################################################

#natlog of BSH vs RPOB
md<-fread("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metagenomic/metaG_metadata_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  dplyr::rename(phase=lightdark,`Sample ID`=sample_name)

natlog_bsh<-fread("species_pfam/rpca_results_BSH_RPOB_rmzero/sample_plot_data_speciespfam_BSH_RPOB.tsv")%>%
  dplyr::select(1:2)%>%
  left_join(.,md,by="Sample ID")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  mutate(condition=factor(condition,levels=c("NA","FA","FT")),
         phase=factor(phase,levels=c("light","dark")),
         Current_Natural_Log_Ratio=as.numeric(Current_Natural_Log_Ratio))

p<-ggplot(natlog_bsh, aes(x=condition, y=Current_Natural_Log_Ratio,fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  theme_minimal()+ scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
  labs(x="condition",y="Natural Log Ratio (BSH/RPOB)",title="Bile Salt Hydrolase (BSH)")+
  theme(legend.position = "none")

pairwise.wilcox.test(natlog_bsh$Current_Natural_Log_Ratio, natlog_bsh$condition,
                     p.adjust.method="none")
#species BSH_RPOB
# NA   FA  
# FA 0.46 -   
#   FT 0.61 0.65

ggsave("species_pfam/rpca_results_BSH_RPOB_rmzero/SFR25_0108_natlog_speciesBSHvsRPOB.pdf", plot=p,height=3.5, width=3.5)
       
p<-ggplot(natlog_bsh, aes(x=condition, y=Current_Natural_Log_Ratio,fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  facet_wrap(~phase)+
  theme_classic()+ scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
  labs(x="condition",y="Natural Log Ratio (BSH/RPOB)",title="Bile Salt Hydrolase (BSH)")+
  theme(legend.position = "none")

ggsave("species_pfam/rpca_results_BSH_RPOB_rmzero/SFR25_0108_natlog_speciesBSHvsRPOB_lightdark.pdf", plot=p,height=3, width=4)
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

natlog_bshFT<-natlog_bsh%>%filter(condition=="FT")
pairwise.wilcox.test(natlog_bshFT$Current_Natural_Log_Ratio, natlog_bshFT$phase,
                     p.adjust.method="none")

# light
# dark 0.54  

natlog_bshFA<-natlog_bsh%>%filter(condition=="FA")
pairwise.wilcox.test(natlog_bshFA$Current_Natural_Log_Ratio, natlog_bshFA$phase,
                     p.adjust.method="none")

# light
# dark 0.24 

natlog_bshNA<-natlog_bsh%>%filter(condition=="NA")
pairwise.wilcox.test(natlog_bshNA$Current_Natural_Log_Ratio, natlog_bshNA$phase,
                     p.adjust.method="none")

# light
# dark 0.26

###########################################################

#plot the distribution of BSH by species (Light/Dark)

mdG<-fread("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metagenomic/metaG_metadata_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  dplyr::rename(phase=lightdark)

selfeat<-fread("species_pfam/rpca_results_BSH_RPOB_rmzero/selected_features_speciespfam_BSH_RPOB.tsv")%>%
  dplyr::filter(grepl("PF02275.21",`Feature ID`))

bsh_dat<-fread("species_pfam/species_pfam_clean_noNT_TPM.tsv")%>%
  filter(FeatureID %in% selfeat$`Feature ID`)%>%
  mutate(FeatureID=gsub("\\|PF02275.21","",FeatureID))%>%
  gather(sample_name,TPM_counts,-FeatureID)%>%
  left_join(.,mdG,by="sample_name")%>%
  group_by(FeatureID,condition,phase)%>%summarise(mn_TPM=mean(TPM_counts), log_mn_TPM=log10(mean(TPM_counts)+1))%>%
  mutate(condition=factor(condition, levels = c("FT", "FA","NA")),
         phase=factor(phase, levels = c("light","dark")))%>%
  #filter(phase=="light")%>%
  filter(phase=="dark")%>%
  group_by(FeatureID)%>%mutate(sum_TPM=sum(mn_TPM), sum_log_TPM=sum(log_mn_TPM))

p<-ggplot(data=bsh_dat, aes(x=reorder(FeatureID, sum_log_TPM), y=log_mn_TPM, fill=condition)) +
  geom_bar(stat="identity") + coord_flip() + theme_classic()+
  scale_fill_manual(values=c("#009E73","#D55E00","#0072B2"))+
  #labs(y="log10(mean_TPM)",x="Species",title="Light")+
  labs(y="log10(mean_TPM)",x="Species",title="Dark")+
  theme(legend.position = "none")+
  scale_x_discrete(expand = c(0, 0)) +
  #scale_y_continuous(limits=c(0,70),expand = c(0, 0))
  scale_y_continuous(expand = c(0, 0))

ggsave("species_pfam/SFR25_0108_BSHlight_logTPM.pdf", height=3, width=3.5)
ggsave("species_pfam/SFR25_0108_BSHdark_logTPM.pdf", height=3, width=3.5)

###########################################################
#plot each of the 13 individual

bsh_dat<-fread("species_pfam/species_pfam_clean_noNT_TPM.tsv")%>%
  filter(FeatureID %in% selfeat$`Feature ID`)%>%
  mutate(FeatureID=gsub("\\|PF02275.21","",FeatureID))%>%
  gather(sample_name,TPM_counts,-FeatureID)%>%
  left_join(.,mdG,by="sample_name")%>%
  mutate(log_TPM=log10(TPM_counts+1),
         condition=factor(condition, levels = c("NA", "FA","FT")),
         phase=factor(phase, levels = c("light","dark")))

p<-ggplot(bsh_dat, aes(x=condition, y=log_TPM,fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  facet_wrap(~FeatureID+phase, ncol=6)+
  theme_classic()+ scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
  labs(x="condition",y="log10 TPM")+
  theme(legend.position = "none")+ theme(plot.title = element_text(face = "bold",hjust = 0.5))
ggsave("species_pfam/indv_abun_plot_faceted.pdf", height=6, width=6)

indv_plot<-function(data,FeatureID){
  
  data%>%
    ggplot(aes(x=condition, y=log_TPM,fill=condition)) +
    geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                           position=position_dodge(1)) +
    facet_wrap(~phase)+
    theme_classic()+ scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
    labs(x="condition",y="log10 TPM",title=stringr::str_wrap(FeatureID, width=30))+
    theme(legend.position = "none")+ theme(plot.title = element_text(face = "bold",hjust = 0.5))
  ggsave(paste("species_pfam/indv_abun_plot/",FeatureID,".pdf",sep=""), height=3, width=3.5)
}

BSHld_pval<-function(df,phase){
  if(phase=="light"){
    subset_df<-df%>%filter(phase=="light")
    x<-pairwise.wilcox.test(subset_df$log_TPM, subset_df$condition, p.adjust.method = "fdr")
    df_pvals <- data.frame(
      NtgALFvTgALF = x[["p.value"]][1],
      NtgALFvTgTRF = x[["p.value"]][2],
      TgALFvTgTRF = x[["p.value"]][4])}
  else{
    subset_df<-df%>%filter(phase=="dark")
    x<-pairwise.wilcox.test(subset_df$log_TPM, subset_df$condition, p.adjust.method = "fdr")
    df_pvals <- data.frame(
      NtgALFvTgALF = x[["p.value"]][1],
      NtgALFvTgTRF = x[["p.value"]][2],
      TgALFvTgTRF = x[["p.value"]][4])}
    return(df_pvals)
  }

mtgBSH_nested <- bsh_dat %>%
  group_by(FeatureID) %>% 
  nest()

mtgBSH_plots <- 
  mtgBSH_nested %>% 
  mutate(plot = map2(data, FeatureID,  ~ indv_plot(.x,.y)))

mtbbdm_pvals_L <- 
  mtgBSH_nested %>% 
  mutate(pvals = map2(data,"light", ~BSHld_pval(.x,.y)))%>%
  dplyr::select(-data)%>%
  unnest()%>%
  mutate(phase="light")

mtbbdm_pvals_D <- 
  mtgBSH_nested %>% 
  mutate(pvals = map2(data,"dark", ~BSHld_pval(.x,.y)))%>%
  dplyr::select(-data)%>%
  unnest()%>%
  mutate(phase="dark")

mtbbdm_pvals<-rbind(mtbbdm_pvals_L,mtbbdm_pvals_D)

#this shows no significance for TRF condition and no light/dark differences

write.table(mtbbdm_pvals,"species_pfam/indv_abun_plot/mtg-BSH-conditionLD-wilcox-pvals.txt",sep = "\t",row.names = FALSE,quote=FALSE)
###########################################################
 
#note birdman run kept failing bc of all the zeroes for NA vs FA and FA vs FT overall and light/dark
