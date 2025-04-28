setwd("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metagenomic/woltka2_results/filtered_metaG/")

library(tidyverse)
library(data.table)
library(qiime2R)
library(ggpubr)

##############################################################################
samps_list<-c("#FeatureID","FA1b_S36","FA1c_S37","FA2b_S38","FA2c_S39","FA3b_S40","FA3c_S41","FA4b_S42","FA4c_S43",  
              "FA5b_S44","FA5c_S45","FA6b_S46","FA6c_S47","FT1b_S48","FT1c_S49","FT2b_S50","FT2c_S51","FT3b_S52",
              "FT3c_S53","FT4b_S54","FT4c_S55","FT5b_S56","FT5c_S57","FT6b_S58","FT6c_S59","NA1a_S1","NA1b_S2",  
              "NA1c_S3","NA2a_S4","NA2b_S5","NA2c_S6","NA3a_S7","NA3b_S8","NA3c_S9","NA4a_S10","NA4b_S11",  
              "NA4c_S12","NA5a_S13","NA5b_S14","NA5c_S15","NA6a_S16","NA6b_S17","NA6c_S18")

dat<-fread("BSH/genome.tsv")%>%dplyr::select(all_of(samps_list))%>%
  mutate(`#FeatureID`=gsub("\\|","_",`#FeatureID`))
names(dat) <- c("FeatureID","FA1b","FA1c","FA2b","FA2c","FA3b","FA3c","FA4b","FA4c",  
                "FA5b","FA5c","FA6b","FA6c","FT1b","FT1c","FT2b","FT2c","FT3b",
                "FT3c","FT4b","FT4c","FT5b","FT5c","FT6b","FT6c","NA1a","NA1b",  
                "NA1c","NA2a","NA2b","NA2c","NA3a","NA3b","NA3c","NA4a","NA4b",  
                "NA4c","NA5a","NA5b","NA5c","NA6a","NA6b","NA6c")

write.table(dat,"BSH/genome_noNT.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

dat<-dat%>%column_to_rownames("FeatureID")
dat_rmzs <- dat[!(rowSums(dat != 0) <2), ]
dat_rmzs<-dat_rmzs%>%rownames_to_column("FeatureID")
write.table(dat_rmzs,"BSH/genome_noNT_rmdbton.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

#just light
samps_list<-c("FeatureID","FA4b","FA4c","FA5b","FA5c","FA6b","FA6c",
              "FT4b","FT4c","FT5b","FT5c","FT6b","FT6c","NA4a","NA4b",  
              "NA4c","NA5a","NA5b","NA5c","NA6a","NA6b","NA6c")

dat_rmzs_l<-dat_rmzs%>%dplyr::select(all_of(samps_list))%>%column_to_rownames("FeatureID")
dat_rmzs_l <- dat_rmzs_l[!(rowSums(dat_rmzs_l != 0) <2), ]
dat_rmzs_l<-dat_rmzs_l%>%rownames_to_column("FeatureID")
write.table(dat_rmzs_l,"BSH/genomeL_noNT_rmdbton.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

#just dark
samps_list<-c("FeatureID","FA1b","FA1c","FA2b","FA2c","FA3b","FA3c",
              "FT1b","FT1c","FT2b","FT2c","FT3b","FT3c","NA1a","NA1b",  
              "NA1c","NA2a","NA2b","NA2c","NA3a","NA3b","NA3c")

dat_rmzs_d<-dat_rmzs%>%dplyr::select(all_of(samps_list))%>%column_to_rownames("FeatureID")
dat_rmzs_d <- dat_rmzs_d[!(rowSums(dat_rmzs_d != 0) <2), ]
dat_rmzs_d<-dat_rmzs_d%>%rownames_to_column("FeatureID")
write.table(dat_rmzs_d,"BSH/genomeD_noNT_rmdbton.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

md<-fread("BSH/BSH_db_metadata_cln.txt")
##############################################################################

#load the birdman results
BSHlight<-fread("BSH/birdman/genomeL_noNT_rmdbton.beta_var.tsv")%>%
  dplyr::rename(ratio=`C(condition, Treatment('FA'))[T.FT]_mean`,
                FeatureID=Feature)%>%
  mutate(`C(condition, Treatment('FA'))[T.FT]_hdi`=gsub("[(]|[)]","",`C(condition, Treatment('FA'))[T.FT]_hdi`))%>%
  separate(`C(condition, Treatment('FA'))[T.FT]_hdi`,c("min","max"), sep=",")%>%
  mutate(min=as.numeric(min),
         max=as.numeric(max),
         credible=ifelse(min>0|max<0,"yes","no"),
         phase="light")%>%
  left_join(.,md,by="FeatureID")%>%
  #filter(credible=="yes")%>%
  select(name, ratio, min, max, phase)%>%
  arrange(ratio)%>%
  filter(name %in% BSHdark$name)


BSHdark<-fread("BSH/birdman/genomeD_noNT_rmdbton.beta_var.tsv")%>%
  dplyr::rename(ratio=`C(condition, Treatment('FA'))[T.FT]_mean`,
                FeatureID=Feature)%>%
  mutate(`C(condition, Treatment('FA'))[T.FT]_hdi`=gsub("[(]|[)]","",`C(condition, Treatment('FA'))[T.FT]_hdi`))%>%
  separate(`C(condition, Treatment('FA'))[T.FT]_hdi`,c("min","max"), sep=",")%>%
  mutate(min=as.numeric(min),
         max=as.numeric(max),
         credible=ifelse(min>0|max<0,"yes","no"),
         phase="dark")%>%
  left_join(.,md,by="FeatureID")%>%
  #filter(credible=="yes")%>%
  select(name, ratio, min, max, phase)%>%
  filter(name %in% BSHlight$name)%>%
  arrange(ratio)

BSH<-rbind(BSHlight,BSHdark)%>%
  mutate(phase=factor(phase,levels=c("light","dark")))%>%
  arrange(ratio)

#BSH$name <- factor(BSH$name,levels = BSHlight$name )
BSH$name <- factor(BSH$name,levels = BSHdark$name )

ggplot(BSH, aes(x =name , y = ratio, ymin = min, ymax = max, color=phase, group=phase)) + 
  geom_linerange(position = position_dodge(width = 0.8)) + theme_pubr()+
  geom_pointrange(position = position_dodge(width = 0.8)) + 
  scale_y_continuous(limits=c(-12.5,10))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  coord_flip() +scale_color_manual(values=c("gray70","gray10"))

ggsave("BSH/birdman/FAFTLD_darkarrange.pdf",height=7, width=14)

##############################################################################
samps_list<-c("#FeatureID","FA1b_S36","FA1c_S37","FA2b_S38","FA2c_S39","FA3b_S40","FA3c_S41","FA4b_S42","FA4c_S43",  
              "FA5b_S44","FA5c_S45","FA6b_S46","FA6c_S47","FT1b_S48","FT1c_S49","FT2b_S50","FT2c_S51","FT3b_S52",
              "FT3c_S53","FT4b_S54","FT4c_S55","FT5b_S56","FT5c_S57","FT6b_S58","FT6c_S59","NA1a_S1","NA1b_S2",  
              "NA1c_S3","NA2a_S4","NA2b_S5","NA2c_S6","NA3a_S7","NA3b_S8","NA3c_S9","NA4a_S10","NA4b_S11",  
              "NA4c_S12","NA5a_S13","NA5b_S14","NA5c_S15","NA6a_S16","NA6b_S17","NA6c_S18")

dat<-fread("BSH/genome_TPM.tsv")%>%dplyr::select(all_of(samps_list))%>%
  mutate(`#FeatureID`=gsub("\\|","_",`#FeatureID`))
names(dat) <- c("FeatureID","FA1b","FA1c","FA2b","FA2c","FA3b","FA3c","FA4b","FA4c",  
                "FA5b","FA5c","FA6b","FA6c","FT1b","FT1c","FT2b","FT2c","FT3b",
                "FT3c","FT4b","FT4c","FT5b","FT5c","FT6b","FT6c","NA1a","NA1b",  
                "NA1c","NA2a","NA2b","NA2c","NA3a","NA3b","NA3c","NA4a","NA4b",  
                "NA4c","NA5a","NA5b","NA5c","NA6a","NA6b","NA6c")

write.table(dat,"BSH/genome-TPM_noNT.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

dat<-dat%>%column_to_rownames("FeatureID")
dat_rmzs <- dat[!(rowSums(dat != 0) <2), ]
dat_rmzs<-dat_rmzs%>%rownames_to_column("FeatureID")
write.table(dat_rmzs,"BSH/genome-TPM_noNT_rmdbton.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

md<-fread("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metagenomic/metaG_metadata_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  dplyr::rename(phase=lightdark)

bsh_TPM<-dat%>%rownames_to_column("FeatureID")%>%
  gather(sample_name,TPM,-FeatureID)%>%
  left_join(.,md,by="sample_name")

p<-ggplot(bsh_TPM, aes(x=condition, y=log10(TPM+1),fill=condition)) +
  geom_boxplot(alpha=0.3) + #geom_dotplot(binaxis='y', stackdir='center',
                                        #position=position_dodge(1)) +
  theme_minimal()+ scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
  labs(x="condition",y="BSH TPM",title="Bile Salt Hydrolase (BSH)")+
  theme(legend.position = "none")

p<-ggplot(bsh_TPM, aes(x=condition, y=log10(TPM+1),fill=condition)) +
  geom_boxplot(alpha=0.3) +
  facet_wrap(~phase)+
  theme_classic()+ scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
  labs(x="condition",y="BSH TPM",title="Bile Salt Hydrolase (BSH)")+
  theme(legend.position = "none")

#natlog_bsh<-fread("BSH/rpca_results_genome/sample_plot_data_10tp.tsv")%>%
natlog_bsh<-fread("BSH/rpca_results_genome/sample_plot_data_20tp.tsv")%>%
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
  theme_minimal()+ scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
  #labs(x="condition",y="Natural Log Ratio (top bottom 10 percent)",title="Bile Salt Hydrolase (BSH)")+
  labs(x="condition",y="Natural Log Ratio (top bottom 20 percent)",title="Bile Salt Hydrolase (BSH)")+
  theme(legend.position = "none")

pairwise.wilcox.test(natlog_bsh$Current_Natural_Log_Ratio, natlog_bsh$condition,
                     p.adjust.method="none")
#10
# NA   FA  
# FA 0.22 -   
#   FT 0.76 0.33

#20
# NA    FA   
# FA 0.562 -    
#   FT 0.485 0.093

p<-ggplot(natlog_bsh, aes(x=condition, y=Current_Natural_Log_Ratio,fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  facet_wrap(~phase)+
  theme_classic()+ scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
  #labs(x="condition",y="Natural Log Ratio (top10/bottom10 BSH)",title="Bile Salt Hydrolase (BSH)")+
  labs(x="condition",y="Natural Log Ratio (top20/bottom20 BSH)",title="Bile Salt Hydrolase (BSH)")+
  theme(legend.position = "none")

#ggsave("BSH/rpca_results_genome/nat_log_10tp.pdf",height=3, width=4)
ggsave("BSH/rpca_results_genome/nat_log_20tp.pdf",height=3, width=4)

natlog_bshL<-natlog_bsh%>%filter(phase=="light")
pairwise.wilcox.test(natlog_bshL$Current_Natural_Log_Ratio, natlog_bshL$condition,
                     p.adjust.method="fdr")

#10
# NA  FA 
# FA 0.4 -  
#   FT 0.4 0.4

#20
# NA  FA 
# FA 0.2 -  
#   FT 0.4 0.2


natlog_bshD<-natlog_bsh%>%filter(phase=="dark")
pairwise.wilcox.test(natlog_bshD$Current_Natural_Log_Ratio, natlog_bshD$condition,
                     p.adjust.method="fdr")
#10
# NA FA
# FA 1  - 
#   FT 1  1 

#20
# NA   FA  
# FA 0.61 -   
#   FT 0.57 0.57


#plot the distribution of BSH by genus (Light)

md<-fread("BSH/BSH_db_metadata_cln.txt")

mdG<-fread("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metagenomic/metaG_metadata_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  dplyr::rename(phase=lightdark)


bsh_dat<-fread("BSH/genome-TPM_noNT_rmdbton.tsv") %>%
  gather(sample_name,TPM_counts,-FeatureID)%>%
  mutate(FeatureID=gsub("sp_","tr_",FeatureID))%>%
  left_join(.,md,by="FeatureID")%>%
  left_join(.,mdG,by="sample_name")%>%
  group_by(name,condition,phase)%>%summarise(mn_TPM=mean(TPM_counts), log_mn_TPM=log10(mean(TPM_counts)+1))%>%
  mutate(condition=factor(condition, levels = c("FT", "FA","NA")),
         phase=factor(phase, levels = c("light","dark")))%>%
  filter(phase=="light")%>%
  #ilter(phase=="dark")%>%
  group_by(name)%>%mutate(sum_TPM=sum(mn_TPM), sum_log_TPM=sum(log_mn_TPM))%>%
  filter(mn_TPM>0,
         name %in% lst_keep)

lst_keep<-((bsh_dat%>%arrange(-sum_log_TPM))$name%>%unique())[1:35]

p<-ggplot(data=bsh_dat, aes(x=reorder(name, sum_log_TPM), y=log_mn_TPM, fill=condition)) +
  geom_bar(stat="identity") + coord_flip() + theme_classic()+
  scale_fill_manual(values=c("#009E73","#D55E00","#0072B2"))+
  labs(y="log10(mean_TPM)",x="Species",title="Light")+
  #labs(y="log10(mean_TPM)",x="Species",title="Dark")+
  theme(legend.position = "none")+
  scale_x_discrete(expand = c(0, 0)) +
  #scale_y_continuous(limits=c(0,70),expand = c(0, 0))
  scale_y_continuous(expand = c(0, 0))

ggsave("BSH/SFR25_0109_BSHlight_logTPM.pdf", height=5, width=10)
ggsave("BSH/SFR25_0109_BSHdark_logTPM.pdf", height=5, width=10)
