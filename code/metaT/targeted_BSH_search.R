setwd("~/scratch/TRF_multiomics/metatranscript/woltka2_m_results")

library(tidyverse)
library(data.table)
library(qiime2R)
library(ggpubr)

##############################################################################
samps_list<-c("#FeatureID","cFA01a_S39","cFA01b_S45","cFA05a_S40","cFA05b_S46","cFA09a_S41",
              "cFA09b_S47","cFA13a_S36","cFA13b_S42","cFA17a_S37","cFA17b_S43",
              "cFA21a_S38","cFA21b_S44","cFT01a_S51","cFT01b_S57","cFT05a_S52",
              "cFT05b_S58","cFT09a_S53","cFT09b_S59","cFT13a_S48","cFT13b_S54",
              "cFT17a_S49","cFT17b_S55","cFT21a_S50","cFT21b_S56","cNA01a_S4",
              "cNA01b_S10","cNA01c_S16","cNA05a_S5","cNA05b_S11","cNA05c_S17",
              "cNA09a_S6","cNA09b_S12","cNA09c_S18","cNA13a_S1","cNA13b_S7",
              "cNA13c_S13","cNA17a_S2","cNA17b_S8","cNA17c_S14","cNA21a_S3",
              "cNA21b_S9","cNA21c_S15")

dat<-fread("BSH/genome.tsv")%>%dplyr::select(all_of(samps_list))%>%
  mutate(`#FeatureID`=gsub("\\|","_",`#FeatureID`))
names(dat) <- c("FeatureID","cFA01a","cFA01b","cFA05a","cFA05b","cFA09a","cFA09b","cFA13a","cFA13b","cFA17a","cFA17b",
                "cFA21a","cFA21b","cFT01a","cFT01b","cFT05a","cFT05b","cFT09a","cFT09b","cFT13a","cFT13b",
                "cFT17a","cFT17b","cFT21a","cFT21b","cNA01a","cNA01b","cNA01c","cNA05a","cNA05b","cNA05c",
                "cNA09a","cNA09b","cNA09c","cNA13a","cNA13b","cNA13c","cNA17a","cNA17b","cNA17c","cNA21a",
                "cNA21b","cNA21c")

write.table(dat,"BSH/genome_noNT.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

dat<-dat%>%column_to_rownames("FeatureID")
dat_rmzs <- dat[!(rowSums(dat != 0) <2), ]
dat_rmzs<-dat_rmzs%>%rownames_to_column("FeatureID")
write.table(dat_rmzs,"BSH/genome_noNT_rmdbton.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

#just light
samps_list<-c("FeatureID","cFA01a","cFA01b","cFA05a","cFA05b","cFA09a","cFA09b",
              "cFT01a","cFT01b","cFT05a","cFT05b","cFT09a","cFT09b",
              "cNA01a","cNA01b","cNA01c","cNA05a","cNA05b","cNA05c",
              "cNA09a","cNA09b","cNA09c")

dat_rmzs_l<-dat_rmzs%>%dplyr::select(all_of(samps_list))%>%column_to_rownames("FeatureID")
dat_rmzs_l <- dat_rmzs_l[!(rowSums(dat_rmzs_l != 0) <2), ]
dat_rmzs_l<-dat_rmzs_l%>%rownames_to_column("FeatureID")
write.table(dat_rmzs_l,"BSH/genomeL_noNT_rmdbton.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

#just dark
samps_list<-c("FeatureID","cFA13a","cFA13b","cFA17a","cFA17b",
              "cFA21a","cFA21b","cFT13a","cFT13b",
              "cFT17a","cFT17b","cFT21a","cFT21b",
              "cNA13a","cNA13b","cNA13c","cNA17a","cNA17b","cNA17c","cNA21a",
              "cNA21b","cNA21c")

dat_rmzs_d<-dat_rmzs%>%dplyr::select(all_of(samps_list))%>%column_to_rownames("FeatureID")
dat_rmzs_d <- dat_rmzs_d[!(rowSums(dat_rmzs_d != 0) <2), ]
dat_rmzs_d<-dat_rmzs_d%>%rownames_to_column("FeatureID")
write.table(dat_rmzs_d,"BSH/genomeD_noNT_rmdbton.tsv",sep = "\t",row.names = FALSE, quote=FALSE)


md<-fread("BSH/BSH_db_metadata.txt")%>%
  mutate(FeatureID=paste("tr_",Entry,"_",`Entry Name`,sep=""),
         name=paste(Organism,Entry,sep="-"))%>%
  dplyr::select(FeatureID, everything())
  #filter(FeatureID %in% dat_rmzs$FeatureID)

write.table(md,"BSH/BSH_db_metadata_cln.txt",sep = "\t",row.names = FALSE, quote=FALSE)
##############################################################################

#songbird results -FAvFT

FAFTL_diff <- read_qza("BSH/songbird_genomeL_FA/differentials.qza")$data%>%
  dplyr::rename(feature_ranking=`C(condition, Treatment('FA'))[T.FT]`)%>%
  mutate(condition=ifelse(feature_ranking<0,"FA","FT")) %>%
  mutate(condition=factor(condition,levels=c("FA","FT")),
         phase="light")%>%
  dplyr::rename(FeatureID=Feature.ID)%>%
  mutate(FeatureID=gsub("sp_","tr_",FeatureID))%>%
  left_join(.,md, by="FeatureID")%>%
  arrange(feature_ranking)

FAFTL_diff_min<-FAFTL_diff%>%slice_min(feature_ranking, n=10)
FAFTL_diff_max<-FAFTL_diff%>%slice_max(feature_ranking, n=10)

combFAFTL<-rbind(FAFTL_diff_min,FAFTL_diff_max)%>%arrange(feature_ranking)

FAFTD_diff <- read_qza("BSH/songbird_genomeD_FA/differentials.qza")$data%>%
  dplyr::rename(feature_ranking=`C(condition, Treatment('FA'))[T.FT]`)%>%
  mutate(condition=ifelse(feature_ranking<0,"FA","FT")) %>%
  mutate(condition=factor(condition,levels=c("FA","FT")),
         phase="dark")%>%
  dplyr::rename(FeatureID=Feature.ID)%>%
  mutate(FeatureID=gsub("sp_","tr_",FeatureID))%>%
  left_join(.,md, by="FeatureID")%>%
  filter(FeatureID %in% combFAFTL$FeatureID)

combFAFT<-rbind(combFAFTL,FAFTD_diff)%>%arrange(feature_ranking)%>%
  mutate(phase=factor(phase,levels=c("light","dark")))
  
#plot the differentials

combFAFT$name <- factor(combFAFT$name,levels = combFAFTL$name)

ggplot(data=combFAFT, aes(x=name, y=feature_ranking, fill=phase)) +
  geom_bar(stat="identity", position=position_dodge()) +coord_flip() +theme_pubr() +
  scale_fill_manual(values=c("gray70","gray10")) +
  scale_y_continuous(expand=c(0,0))

ggsave("BSH/songbird_genomeL_FA/SFR23_0717_TgBSHsbdiffFTFA_tb15.pdf",height=5, width=9)

#songbird results -FAvNA

FANAL_diff <- read_qza("BSH/songbird_genomeL/differentials.qza")$data%>%
  dplyr::rename(feature_ranking=`C(condition, Treatment('NA'))[T.FA]`)%>%
  mutate(condition=ifelse(feature_ranking<0,"NA","FA")) %>%
  mutate(condition=factor(condition,levels=c("NA","FA")),
         phase="light")%>%
  dplyr::rename(FeatureID=Feature.ID)%>%
  mutate(FeatureID=gsub("sp_","tr_",FeatureID))%>%
  left_join(.,md, by="FeatureID")%>%
  arrange(feature_ranking)

FANAL_diff_min<-FANAL_diff%>%slice_min(feature_ranking, n=10)
FANAL_diff_max<-FANAL_diff%>%slice_max(feature_ranking, n=10)

combFANAL<-rbind(FANAL_diff_min,FANAL_diff_max)%>%arrange(feature_ranking)

FANAD_diff <- read_qza("BSH/songbird_genomeD/differentials.qza")$data%>%
  dplyr::rename(feature_ranking=`C(condition, Treatment('NA'))[T.FA]`)%>%
  mutate(condition=ifelse(feature_ranking<0,"NA","FA")) %>%
  mutate(condition=factor(condition,levels=c("NA","FA")),
         phase="dark")%>%
  dplyr::rename(FeatureID=Feature.ID)%>%
  mutate(FeatureID=gsub("sp_","tr_",FeatureID))%>%
  left_join(.,md, by="FeatureID")%>%
  filter(FeatureID %in% combFANAL$FeatureID)

combFANA<-rbind(combFANAL,FANAD_diff)%>%arrange(feature_ranking)%>%
  mutate(phase=factor(phase,levels=c("light","dark")))

#plot the differentials

combFANA$name <- factor(combFANA$name,levels = combFANAL$name)

ggplot(data=combFANA, aes(x=name, y=feature_ranking, fill=phase)) +
  geom_bar(stat="identity", position=position_dodge()) +coord_flip() +theme_pubr() +
  scale_fill_manual(values=c("gray70","gray10")) +
  scale_y_continuous(expand=c(0,0))

ggsave("BSH/songbird_genomeL/SFR23_0717_TgBSHsbdiffFANA_tb15.pdf",height=5, width=9)

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
  dplyr::select(name, ratio, min, max, phase)#%>%
  #arrange(ratio)%>%
  #filter(name %in% BSHdark$name)


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
  dplyr::select(name, ratio, min, max, phase)#%>%
  #filter(name %in% BSHlight$name)
  #arrange(ratio)

selfeat<-rbind(BSHlight,BSHdark)%>%select(1:2,5)%>%
  spread(phase,ratio)%>%
  mutate(diff_m=abs(light-dark))%>%
  filter(diff_m>3)%>%
  arrange(desc(diff_m))

BSH<-rbind(BSHlight,BSHdark)%>%
  mutate(phase=factor(phase,levels=c("light","dark")))#%>%
  #filter(name %in% selfeat$name)
write.table(BSH,"BSH/birdman/genomeLD_noNT_rmdbton.beta_var.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

BSH$name <- factor(BSH$name,levels = selfeat$name)
BSH$name <- factor(BSH$name,levels = BSHlight$name )
#BSH$name <- factor(BSH$name,levels = BSHdark$name )

ggplot(BSH, aes(x =name , y = ratio, ymin = min, ymax = max, color=phase, group=phase)) + 
  geom_linerange(position = position_dodge(width = 0.8)) + theme_pubr()+
  geom_pointrange(position = position_dodge(width = 0.8)) + 
  scale_y_continuous(limits=c(-12.5,10))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  coord_flip() +scale_color_manual(values=c("gray70","gray10"))

#ggsave("BSH/birdman/FAFTLD_Dcredbgdiff_smhdi.pdf",height=7, width=13)
ggsave("BSH/birdman/FAFTLD_Dcredbgdiff_smhdi_v2.pdf",height=7, width=14.5)
ggsave("BSH/birdman/FAFTLD_Dcredbgdiff_smhdi_darkarranged.pdf",height=7, width=14.5)
#ggsave("BSH/birdman/FAFTLD_Lcred_smhdi.pdf",height=7, width=7)
#ggsave("BSH/birdman/FAFTLD_Dcred_smhdi.pdf",height=15, width=13)


##############################################################################
samps_list<-c("#FeatureID","cFA01a_S39","cFA01b_S45","cFA05a_S40","cFA05b_S46","cFA09a_S41",
              "cFA09b_S47","cFA13a_S36","cFA13b_S42","cFA17a_S37","cFA17b_S43",
              "cFA21a_S38","cFA21b_S44","cFT01a_S51","cFT01b_S57","cFT05a_S52",
              "cFT05b_S58","cFT09a_S53","cFT09b_S59","cFT13a_S48","cFT13b_S54",
              "cFT17a_S49","cFT17b_S55","cFT21a_S50","cFT21b_S56","cNA01a_S4",
              "cNA01b_S10","cNA01c_S16","cNA05a_S5","cNA05b_S11","cNA05c_S17",
              "cNA09a_S6","cNA09b_S12","cNA09c_S18","cNA13a_S1","cNA13b_S7",
              "cNA13c_S13","cNA17a_S2","cNA17b_S8","cNA17c_S14","cNA21a_S3",
              "cNA21b_S9","cNA21c_S15")

dat<-fread("BSH/genome-TPM.tsv")%>%dplyr::select(all_of(samps_list))%>%
  mutate(`#FeatureID`=gsub("\\|","_",`#FeatureID`))
names(dat) <- c("FeatureID","cFA01a","cFA01b","cFA05a","cFA05b","cFA09a","cFA09b","cFA13a","cFA13b","cFA17a","cFA17b",
                "cFA21a","cFA21b","cFT01a","cFT01b","cFT05a","cFT05b","cFT09a","cFT09b","cFT13a","cFT13b",
                "cFT17a","cFT17b","cFT21a","cFT21b","cNA01a","cNA01b","cNA01c","cNA05a","cNA05b","cNA05c",
                "cNA09a","cNA09b","cNA09c","cNA13a","cNA13b","cNA13c","cNA17a","cNA17b","cNA17c","cNA21a",
                "cNA21b","cNA21c")

write.table(dat,"BSH/genome-TPM_noNT.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

dat<-dat%>%column_to_rownames("FeatureID")
dat_rmzs <- dat[!(rowSums(dat != 0) <2), ]
dat_rmzs<-dat_rmzs%>%rownames_to_column("FeatureID")
write.table(dat_rmzs,"BSH/genome-TPM_noNT_rmdbton.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

md<-fread("metaT_metadata_ztcat_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))

bsh_TPM<-dat%>%rownames_to_column("FeatureID")%>%
  gather(sample_name,TPM,-FeatureID)%>%
  left_join(.,md,by="sample_name")

p<-ggplot(bsh_TPM, aes(x=condition, y=log10(TPM+1),fill=condition)) +
  geom_boxplot(alpha=0.3) + #geom_dotplot(binaxis='y', stackdir='center',
                                        # position=position_dodge(1)) +
  theme_minimal()+ scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
  labs(x="condition",y="BSH TPM",title="Bile Salt Hydrolase (BSH)")+
  theme(legend.position = "none")

p<-ggplot(bsh_TPM, aes(x=condition, y=log10(TPM+1),fill=condition)) +
  geom_boxplot(alpha=0.3) +
  facet_wrap(~phase)+
  theme_classic()+ scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
  labs(x="condition",y="BSH TPM",title="Bile Salt Hydrolase (BSH)")+
  theme(legend.position = "none")

#natlog_bsh<-fread("BSH/rpca_results_genome/sample_plot_data_suboclos.tsv")%>%
natlog_bsh<-fread("BSH/rpca_results_genome/sample_plot_data_10tp.tsv")%>%
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
  labs(x="condition",y="Natural Log Ratio (BSH/RPOB)",title="Bile Salt Hydrolase (BSH)")+
  theme(legend.position = "none")

pairwise.wilcox.test(natlog_bsh$Current_Natural_Log_Ratio, natlog_bsh$condition,
                     p.adjust.method="none")

# NA      FA  
# FA 2.3e-08 -   
#   FT 1.6e-07 0.38

p<-ggplot(natlog_bsh, aes(x=condition, y=Current_Natural_Log_Ratio,fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  facet_wrap(~phase)+
  theme_classic()+ scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
  labs(x="condition",y="Natural Log Ratio (top10/bottom10 BSH)",title="Bile Salt Hydrolase (BSH)")+
  theme(legend.position = "none")

ggsave("BSH/rpca_results_genome/nat_log_10tp.pdf",height=3, width=4)

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


#plot the distribution of BSH by genus (Light)

md<-fread("BSH/BSH_db_metadata.txt")%>%
  mutate(FeatureID=paste("tr_",Entry,"_",`Entry Name`,sep=""),
         name=paste(Organism,Entry,sep="-"))%>%
  select(FeatureID, everything())

mdT<-fread("metaT_metadata_ztcat_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))


bsh_dat<-fread("BSH/genome-TPM_noNT_rmdbton.tsv") %>%
  gather(sample_name,TPM_counts,-FeatureID)%>%
  mutate(FeatureID=gsub("sp_","tr_",FeatureID))%>%
  left_join(.,md,by="FeatureID")%>%
  left_join(.,mdT,by="sample_name")%>%
  group_by(name,condition,phase)%>%summarise(mn_TPM=mean(TPM_counts), log_mn_TPM=log10(mean(TPM_counts)+1))%>%
  mutate(condition=factor(condition, levels = c("FT", "FA","NA")),
         phase=factor(phase, levels = c("light","dark")))%>%
  #filter(phase=="light")%>%
  filter(phase=="dark")%>%
  group_by(name)%>%mutate(sum_TPM=sum(mn_TPM), sum_log_TPM=sum(log_mn_TPM))%>%
  filter(mn_TPM>0,
         name %in% lst_keep)

lst_keep<-((bsh_dat%>%arrange(-sum_log_TPM))$name%>%unique())[1:35]

p<-ggplot(data=bsh_dat, aes(x=reorder(name, sum_log_TPM), y=log_mn_TPM, fill=condition)) +
  geom_bar(stat="identity") + coord_flip() + theme_classic()+
  scale_fill_manual(values=c("#009E73","#D55E00","#0072B2"))+
  #labs(y="log10(mean_TPM)",x="Species",title="Light")+
  labs(y="log10(mean_TPM)",x="Species",title="Dark")+
  theme(legend.position = "none")+
  scale_x_discrete(expand = c(0, 0)) +
  #scale_y_continuous(limits=c(0,70),expand = c(0, 0))
  scale_y_continuous(expand = c(0, 0))

ggsave("BSH/SFR23_0719_BSHlight_logTPM.pdf", height=5, width=7)
ggsave("BSH/SFR23_0719_BSHdark_logTPM.pdf", height=5, width=7)

##############################################################################

#just dubosiella newyorkensis BSH abundance-TPM

bsh_dat<-fread("BSH/genome-TPM_noNT_rmdbton.tsv") %>%
  gather(sample_name,TPM_counts,-FeatureID)%>%
  mutate(FeatureID=gsub("sp_","tr_",FeatureID))%>%
  left_join(.,md,by="FeatureID")%>%
  left_join(.,mdT,by="sample_name")%>%
  mutate(log_TPM=log2(TPM_counts+1),
         condition=factor(condition, levels = c("NA", "FA","FT")),
         phase=factor(phase, levels = c("light","dark")))%>%
  filter(Organism=="Dubosiella newyorkensis")

bsh1<-bsh_dat%>%filter(Entry=="A0A1U7NKD7")
bsh2<-bsh_dat%>%filter(Entry=="A0A1U7NP31")

p<-ggplot(bsh1, aes(x=condition, y=TPM_counts,fill=condition)) +
#<-ggplot(bsh1, aes(x=condition, y=log_TPM,fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  facet_wrap(~phase)+
  theme_classic()+ scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
  labs(x="condition",y="log2(TPM)",title="D. newyorkensis BSH A0A1U7NKD7")+
  theme(legend.position = "none")

ggsave("BSH/SFR25_0103_DnewyorkensisBSH1_TPM.pdf", height=3, width=4)
ggsave("BSH/SFR25_0103_DnewyorkensisBSH1_log2TPM.pdf", height=3, width=4)

bsh1L<-bsh1%>%filter(phase=="light")
pairwise.wilcox.test(bsh1L$TPM_counts, bsh1L$condition,p.adjust.method="fdr")
pairwise.wilcox.test(bsh1L$log_TPM, bsh1L$condition,p.adjust.method="fdr")
# NA   FA  
# FA 0.73 -   
#   FT 0.73 0.73

bsh1D<-bsh1%>%filter(phase=="dark")
pairwise.wilcox.test(bsh1D$TPM_counts, bsh1D$condition,p.adjust.method="fdr")
pairwise.wilcox.test(bsh1D$log_TPM, bsh1D$condition,p.adjust.method="fdr")
# NA    FA   
# FA 0.077 -    
#   FT 0.082 0.029

p<-ggplot(bsh2, aes(x=condition, y=TPM_counts,fill=condition)) +
#p<-ggplot(bsh2, aes(x=condition, y=log_TPM,fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  facet_wrap(~phase)+
  theme_classic()+ scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
  labs(x="condition",y="log2(TPM)",title="D. newyorkensis BSH A0A1U7NP31")+
  theme(legend.position = "none")

ggsave("BSH/SFR25_0103_DnewyorkensisBSH2_TPM.pdf", height=3, width=4)
ggsave("BSH/SFR25_0103_DnewyorkensisBSH2_log2TPM.pdf", height=3, width=4)

bsh2L<-bsh2%>%filter(phase=="light")
pairwise.wilcox.test(bsh2L$TPM_counts, bsh2L$condition,p.adjust.method="fdr")
pairwise.wilcox.test(bsh2L$log_TPM, bsh2L$condition,p.adjust.method="fdr")
# NA   FA  
# FA 0.66 -   
#   FT 0.66 1.00

bsh2D<-bsh2%>%filter(phase=="dark")
pairwise.wilcox.test(bsh2D$TPM_counts, bsh2D$condition,p.adjust.method="fdr")
pairwise.wilcox.test(bsh2D$log_TPM, bsh2D$condition,p.adjust.method="fdr")
# NA    FA   
# FA 0.192 -    
#   FT 0.192 0.029

##############################################################################

#just dubosiella newyorkensis BSH abundance-norm

#natlog_nkd7<-fread("BSH/rpca_results_genome/sample_plot_data_NKD7b10.tsv")%>%
natlog_np31<-fread("BSH/rpca_results_genome/sample_plot_data_NP31b10.tsv")%>%
  dplyr::select(1:2)%>%
  dplyr::rename(sample_name=`Sample ID`)%>%
  left_join(.,mdT,by="sample_name")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  mutate(condition=factor(condition,levels=c("NA","FA","FT")),
         phase=factor(phase,levels=c("light","dark")),
         Current_Natural_Log_Ratio=as.numeric(Current_Natural_Log_Ratio))

p<-ggplot(natlog_np31, aes(x=condition, y=Current_Natural_Log_Ratio,fill=condition)) +
#p<-ggplot(natlog_nkd7, aes(x=condition, y=Current_Natural_Log_Ratio,fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  theme_minimal()+ scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
  labs(x="condition",y="Natural Log Ratio (NP31/b10%)",title="Bile Salt Hydrolase (BSH)")+
  theme(legend.position = "none")

pairwise.wilcox.test(natlog_nkd7$Current_Natural_Log_Ratio, natlog_nkd7$condition,
                     p.adjust.method="fdr")

# NA   FA  
# FA 0.67 -   
#   FT 1.00 1.00

pairwise.wilcox.test(natlog_np31$Current_Natural_Log_Ratio, natlog_np31$condition,
                     p.adjust.method="fdr")

# NA  FA 
# FA 0.8 -  
#   FT 1.0 0.8


p<-ggplot(natlog_np31, aes(x=condition, y=Current_Natural_Log_Ratio,fill=condition)) +
#p<-ggplot(natlog_nkd7, aes(x=condition, y=Current_Natural_Log_Ratio,fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  facet_wrap(~phase)+
  theme_classic()+ scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
  labs(x="condition",y="Natural Log Ratio (NP31/bottom10 BSH)",title="Bile Salt Hydrolase (BSH)")+
  theme(legend.position = "none")

ggsave("BSH/rpca_results_genome/nat_log_NKD710b.pdf",height=3, width=4)
ggsave("BSH/rpca_results_genome/nat_log_NP3110b.pdf",height=3, width=4)

natlog_nkd7L<-natlog_nkd7%>%filter(phase=="light")
pairwise.wilcox.test(natlog_nkd7L$Current_Natural_Log_Ratio, natlog_nkd7L$condition,
                     p.adjust.method="fdr")

# NA FA
# FA 1  - 
#   FT 1  1 

natlog_nkd7D<-natlog_nkd7%>%filter(phase=="dark")
pairwise.wilcox.test(natlog_nkd7D$Current_Natural_Log_Ratio, natlog_nkd7D$condition,
                     p.adjust.method="fdr")

# Error in wilcox.test.default(xi, xj, paired = paired, ...) : 
#   not enough (non-missing) 'x' observations

natlog_np31L<-natlog_np31%>%filter(phase=="light")
pairwise.wilcox.test(natlog_np31L$Current_Natural_Log_Ratio, natlog_np31L$condition,
                     p.adjust.method="fdr")
# NA   FA  
# FA 0.86 -   
#   FT 0.86 0.86

natlog_np31D<-natlog_np31%>%filter(phase=="dark")
pairwise.wilcox.test(natlog_np31D$Current_Natural_Log_Ratio, natlog_np31D$condition,
                     p.adjust.method="fdr")

# Error in wilcox.test.default(xi, xj, paired = paired, ...) : 
#   not enough (non-missing) 'x' observations

##############################################################################

#plot ANI for BSH of interest

ani<-fread("BSH/BSH_ENB/prot_seqn_analysis_foley/BSH_ani_subset.txt")%>%
  dplyr::select(-V1)

# names(ani) <- c("Group1","G000364225_101","prot_1","prot_4","prot_5","prot_6","tr|A0A1U7NKD7|A0A1U7NKD7_9FIRM",
#                 "tr|A0A1U7NP31|A0A1U7NP31_9FIRM","EEU29635","KFOAOCCA_00205","LGAS_RS04710","KFOAOCCA_00001",
#                 "WP_003547395.1","KRL87916","KRN10896.1","ABQ82994","LBA0892","WP_005718943.1","WP_150399422.1",                
#                 "JIJPODMP_01796","AZ52_04447","SFE73109.1","G009917455_1208","WP_056959219.1","EFK28582.1",                    
#                 "prot_3","prot_2","KRM51566.1","KFOAOCCA_01485","LGAS_RS00260","G000014425_51")

names(ani)<-c("Group1","G000364225_101","prot_1","prot_4","prot_2","KRM51566.1","KFOAOCCA_01485","LGAS_RS00260",
              "G000014425_51","tr|A0A1U7NKD7|A0A1U7NKD7_9FIRM","prot_5","prot_6","tr|A0A1U7NP31|A0A1U7NP31_9FIRM",
              "JIJPODMP_01796","AZ52_04447","SFE73109.1","G009917455_1208","WP_056959219.1","EFK28582.1","prot_3" )

ord_list<-c("G000364225_101","prot_1","prot_4","prot_2","KRM51566.1","KFOAOCCA_01485","LGAS_RS00260",
           "G000014425_51","tr|A0A1U7NKD7|A0A1U7NKD7_9FIRM","prot_5","prot_6","tr|A0A1U7NP31|A0A1U7NP31_9FIRM",
           "JIJPODMP_01796","AZ52_04447","SFE73109.1","G009917455_1208","WP_056959219.1","EFK28582.1","prot_3")

ani<-ani%>%
  gather(Group2,ani_values,-Group1)%>%
  mutate(Group1=factor(Group1,levels=ord_list),
         Group2=factor(Group2,levels=ord_list))


p<-ggplot(ani, aes(x = Group2, y = Group1, fill = ani_values)) +
  geom_tile() +
  geom_text(aes(label = signif(ani_values,digits=3))) +
  scale_fill_distiller(palette = "Reds",direction=1) +
  theme_minimal() +
  theme(panel.grid = element_blank(),legend.position = "top")+labs(title="Percent Indentity Matrix",x="",y="")
ggsave("BSH/BSH_ENB/prot_seqn_analysis_foley/SFR25_0207_anisub.pdf",plot=p, width = 10, height = 10)





