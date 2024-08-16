setwd("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_m_results/")

library(tidyverse)
library(MetaCycle)
library(data.table)
library("qiime2R")
library("Biostrings")
library(ggrepel)
library(ggpubr)
library(gplots)
library(ggvenn)
library(ggbreak)
library(viridis)
library(scatterplot3d)
###########################################################
samps_list<-c("#FeatureID","cFA01a_S39.","cFA01b_S45.","cFA05a_S40.","cFA05b_S46.","cFA09a_S41.",
              "cFA09b_S47.","cFA13a_S36.","cFA13b_S42.","cFA17a_S37.","cFA17b_S43.",
              "cFA21a_S38.","cFA21b_S44.","cFT01a_S51.","cFT01b_S57.","cFT05a_S52.",
              "cFT05b_S58.","cFT09a_S53.","cFT09b_S59.","cFT13a_S48.","cFT13b_S54.",
              "cFT17a_S49.","cFT17b_S55.","cFT21a_S50.","cFT21b_S56.","cNA01a_S4.",
              "cNA01b_S10.","cNA01c_S16.","cNA05a_S5.","cNA05b_S11.","cNA05c_S17.",
              "cNA09a_S6.","cNA09b_S12.","cNA09c_S18.","cNA13a_S1.","cNA13b_S7.",
              "cNA13c_S13.","cNA17a_S2.","cNA17b_S8.","cNA17c_S14.","cNA21a_S3.",
              "cNA21b_S9.","cNA21c_S15.")

#dat<-fread("species_pfam/species_pfam.tsv")%>%dplyr::select(all_of(samps_list))
dat<-fread("species_pfam/species_pfam-TPM.tsv")%>%dplyr::select(all_of(samps_list))
names(dat) <- c("FeatureID","cFA01a","cFA01b","cFA05a","cFA05b","cFA09a","cFA09b","cFA13a","cFA13b","cFA17a","cFA17b",
                "cFA21a","cFA21b","cFT01a","cFT01b","cFT05a","cFT05b","cFT09a","cFT09b","cFT13a","cFT13b",
                "cFT17a","cFT17b","cFT21a","cFT21b","cNA01a","cNA01b","cNA01c","cNA05a","cNA05b","cNA05c",
                "cNA09a","cNA09b","cNA09c","cNA13a","cNA13b","cNA13c","cNA17a","cNA17b","cNA17c","cNA21a",
                "cNA21b","cNA21c")
#write.table(dat,"species_pfam/species_pfam_clean_noNT.tsv",sep = "\t",row.names = FALSE, quote=FALSE)
write.table(dat,"species_pfam/species_pfam_clean_noNT-TPM.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

#remove singletons and doubletons
dat<-dat%>%column_to_rownames("FeatureID")

dat_rmz <- dat[!(rowSums(dat != 0) == 0), ]
dat_rmz<-dat_rmz%>%rownames_to_column("FeatureID")
#write.table(dat_rmz,"species_pfam/species_pfam_clean_rmzero_noNT.tsv",sep = "\t",row.names = FALSE, quote=FALSE)
write.table(dat_rmz,"species_pfam/species_pfam_clean_rmzero_noNT-TPM.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

dat_rmzs <- dat[!(rowSums(dat != 0) <2), ]
dat_rmzs<-dat_rmzs%>%rownames_to_column("FeatureID")
#write.table(dat_rmzs,"species_pfam/species_pfam_clean_rmsingletons_noNT.tsv",sep = "\t",row.names = FALSE, quote=FALSE)
write.table(dat_rmzs,"species_pfam/species_pfam_clean_rmsingletons_noNT-TPM.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

dat_rmzd <- dat[!(rowSums(dat != 0) < 3), ]
dat_rmzd<-dat_rmzd%>%rownames_to_column("FeatureID")
#write.table(dat_rmzd,"species_pfam/species_pfam_clean_rmdoubletons_noNT.tsv",sep = "\t",row.names = FALSE, quote=FALSE)
write.table(dat_rmzd,"species_pfam/species_pfam_clean_rmdoubletons_noNT-TPM.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

#just light
samps_list<-c("FeatureID","cFA01a","cFA01b","cFA05a","cFA05b","cFA09a","cFA09b",
              "cFT01a","cFT01b","cFT05a","cFT05b","cFT09a","cFT09b",
              "cNA01a","cNA01b","cNA01c","cNA05a","cNA05b","cNA05c",
              "cNA09a","cNA09b","cNA09c")

#dat<-fread("genus_pfam/genus_pfam_clean_rmdoubletons_noNT.tsv")%>%dplyr::select(all_of(samps_list))
dat<-fread("species_pfam/species_pfam_clean_rmdoubletons_noNT.tsv")%>%dplyr::select(all_of(samps_list))
#write.table(dat,"genus_pfam/genus_pfam_clean_rmdoubletons_noNT_light.tsv",sep = "\t",row.names = FALSE, quote=FALSE)
write.table(dat,"species_pfam/species_pfam_clean_rmdoubletons_noNT_light.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

#just dark
samps_list<-c("FeatureID","cFA13a","cFA13b","cFA17a","cFA17b",
              "cFA21a","cFA21b","cFT13a","cFT13b",
              "cFT17a","cFT17b","cFT21a","cFT21b",
              "cNA13a","cNA13b","cNA13c","cNA17a","cNA17b","cNA17c","cNA21a",
              "cNA21b","cNA21c")

#dat<-fread("genus_pfam/genus_pfam_clean_rmdoubletons_noNT.tsv")%>%dplyr::select(all_of(samps_list))
dat<-fread("species_pfam/species_pfam_clean_rmdoubletons_noNT.tsv")%>%dplyr::select(all_of(samps_list))
#write.table(dat,"genus_pfam/genus_pfam_clean_rmdoubletons_noNT_dark.tsv",sep = "\t",row.names = FALSE, quote=FALSE)
write.table(dat,"species_pfam/species_pfam_clean_rmdoubletons_noNT_dark.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

dat_rmz5 <- dat[!(rowSums(dat != 0) < 5), ]
dat_rmz5<-dat_rmz5%>%rownames_to_column("FeatureID")
#write.table(dat_rmz5,"genus_pfam/genus_pfam_clean_rminless5samps_noNT.tsv",sep = "\t",row.names = FALSE, quote=FALSE)
write.table(dat_rmz5,"species_pfam/species_pfam_clean_rminless5samps_noNT.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

dat_rmz8 <- dat[!(rowSums(dat != 0) < 8), ] #20%
dat_rmz8<-dat_rmz8%>%rownames_to_column("FeatureID")
#write.table(dat_rmz8,"genus_pfam/genus_pfam_clean_rminless8samps_noNT.tsv",sep = "\t",row.names = FALSE, quote=FALSE)
write.table(dat_rmz8,"species_pfam/species_pfam_clean_rminless8samps_noNT.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

###########################################################

#subset the data to just have BSH
# bsh<-fread("genus_pfam/genus_pfam_clean_noNT.tsv")%>%
#   dplyr::filter(grepl("PF02275.21",FeatureID))

bsh<-fread("species_pfam/species_pfam_clean_noNT.tsv")%>%
  dplyr::filter(grepl("PF02275.21",FeatureID))%>%
  separate(FeatureID,c("FeatureID",NA), sep="\\|PF", extra="drop")%>%
  mutate(FeatureID=gsub(" ","_",FeatureID))

bsh<-bsh%>%column_to_rownames("FeatureID")
bsh_rmz <- bsh[!(rowSums(bsh != 0) ==0), ]
bsh_rmz<-bsh_rmz%>%rownames_to_column("FeatureID")
#write.table(bsh_rmz,"genus_pfam/genus_pfam_BSHonly_clean_rmzero_noNT.tsv",sep = "\t",row.names = FALSE, quote=FALSE)
write.table(bsh_rmz,"species_pfam/species_pfam_BSHonly_clean_rmzero_noNT.tsv",sep = "\t",row.names = FALSE, quote=FALSE)
#just light
samps_list<-c("FeatureID","cFA01a","cFA01b","cFA05a","cFA05b","cFA09a","cFA09b",
              "cFT01a","cFT01b","cFT05a","cFT05b","cFT09a","cFT09b",
              "cNA01a","cNA01b","cNA01c","cNA05a","cNA05b","cNA05c",
              "cNA09a","cNA09b","cNA09c")

#dat<-fread("genus_pfam/genus_pfam_BSHonly_clean_rmzero_noNT.tsv")%>%dplyr::select(all_of(samps_list))
#write.table(dat,"genus_pfam/genus_pfam_BSHonly_clean_rmzero_noNT_light.tsv",sep = "\t",row.names = FALSE, quote=FALSE)
dat<-fread("species_pfam/species_pfam_BSHonly_clean_rmzero_noNT.tsv")%>%dplyr::select(all_of(samps_list))
dat_rmz <- dat[!(rowSums(dat != 0) <2), ]
write.table(dat,"species_pfam/species_pfam_BSHonly_clean_rmzero_noNT_light.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

#just dark
samps_list<-c("FeatureID","cFA13a","cFA13b","cFA17a","cFA17b",
              "cFA21a","cFA21b","cFT13a","cFT13b",
              "cFT17a","cFT17b","cFT21a","cFT21b",
              "cNA13a","cNA13b","cNA13c","cNA17a","cNA17b","cNA17c","cNA21a",
              "cNA21b","cNA21c")

#dat<-fread("genus_pfam/genus_pfam_BSHonly_clean_rmzero_noNT.tsv")%>%dplyr::select(all_of(samps_list))
#write.table(dat,"genus_pfam/genus_pfam_BSHonly_clean_rmzero_noNT_dark.tsv",sep = "\t",row.names = FALSE, quote=FALSE)
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
#write.table(bsh_rmz,"genus_pfam/genus_pfam_BSH_RPOB_clean_rmzero_noNT.tsv",sep = "\t",row.names = FALSE, quote=FALSE)
write.table(bsh_rmz,"species_pfam/species_pfam_BSH_RPOB_clean_rmzero_noNT.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

#just light
samps_list<-c("FeatureID","cFA01a","cFA01b","cFA05a","cFA05b","cFA09a","cFA09b",
              "cFT01a","cFT01b","cFT05a","cFT05b","cFT09a","cFT09b",
              "cNA01a","cNA01b","cNA01c","cNA05a","cNA05b","cNA05c",
              "cNA09a","cNA09b","cNA09c")

#dat<-fread("genus_pfam/genus_pfam_BSHonly_clean_rmzero_noNT.tsv")%>%dplyr::select(all_of(samps_list))
#write.table(dat,"genus_pfam/genus_pfam_BSHonly_clean_rmzero_noNT_light.tsv",sep = "\t",row.names = FALSE, quote=FALSE)
dat<-fread("species_pfam/species_pfam_BSHonly_clean_rmzero_noNT.tsv")%>%dplyr::select(all_of(samps_list))
write.table(dat,"species_pfam/species_pfam_BSHonly_clean_rmzero_noNT_light.tsv",sep = "\t",row.names = FALSE, quote=FALSE)


#just dark
samps_list<-c("FeatureID","cFA13a","cFA13b","cFA17a","cFA17b",
              "cFA21a","cFA21b","cFT13a","cFT13b",
              "cFT17a","cFT17b","cFT21a","cFT21b",
              "cNA13a","cNA13b","cNA13c","cNA17a","cNA17b","cNA17c","cNA21a",
              "cNA21b","cNA21c")

#dat<-fread("genus_pfam/genus_pfam_BSHonly_clean_rmzero_noNT.tsv")%>%dplyr::select(all_of(samps_list))
#write.table(dat,"genus_pfam/genus_pfam_BSHonly_clean_rmzero_noNT_dark.tsv",sep = "\t",row.names = FALSE, quote=FALSE)
dat<-fread("species_pfam/species_pfam_BSHonly_clean_rmzero_noNT.tsv")%>%dplyr::select(all_of(samps_list))
write.table(dat,"species_pfam/species_pfam_BSHonly_clean_rmzero_noNT_dark.tsv",sep = "\t",row.names = FALSE, quote=FALSE)
###########################################################

##look at RPCA

md<-fread("metaT_metadata_ztcat_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))

ord <- read_qza("genus_pfam/rpca_results_rmdoubletons/ordination.qza")

samp_ord<-ord$data$Vectors
write.table(samp_ord,"genus_pfam/rpca_results_rmdoubletons/sample_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)
feat_ord<-ord$data$Species
write.table(feat_ord,"genus_pfam/rpca_results_rmdoubletons/feature_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)

rpca<-ord$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2)%>%
  dplyr::rename(sample_name=SampleID)%>%
  left_join(md,by="sample_name")%>%
  mutate(condition=factor(condition,levels=c("NA","FA","FT")),
         phase=factor(phase,levels=c("light","dark")))

p<-rpca %>%
  ggplot(aes(x=PC1, y=PC2, color=condition, shape=phase)) +
  geom_point(alpha=1.0) + 
  theme_pubr() +stat_ellipse(type = "t", linetype = 2,aes(group = condition))+
  scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+
  scale_shape_manual(values=c(3,16)) +
  labs(color="condition",
       x =paste("PC1 (",round(ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""))+ggtitle("genus|pfam RPCA")+ theme(plot.title = element_text(face = "bold"))
ggsave("genus_pfam/rpca_results_rmdoubletons/SFR23_0630_genupfam_RPCA.pdf", plot=p,height=3.5, width=3.5)


pdf("genus_pfam/rpca_results_rmdoubletons/SFR23_0630_metaT_genupfam_RPCA.pdf",width = 5, height = 5)
# Create a 3D PCA plot
colors <- c("#0072B2","#D55E00","#009E73")
colors <- colors[as.numeric(rpca$condition)]

scatterplot3d(rpca[,2:4],angle = 320,pch = 16,color=colors,grid=FALSE,
              xlab = "PC1", ylab = "PC2", zlab = "PC3",
              main = "genus|pfam MTX")
legend("bottom", legend = levels(rpca$condition),
       col =  c("#0072B2","#D55E00","#009E73"), pch = 16, 
       inset = -0.25, xpd = TRUE, horiz = TRUE)

dev.off()

###shannon-rarefied
md<-fread("metaT_metadata_ztcat_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))

alpha <- read_qza("genus_pfam/diversity-core-metrics1.6M/shannon_vector.qza")$data %>%
  rownames_to_column("sample_name") %>%
  left_join(.,md,by="sample_name") %>%
  mutate(condition=factor(condition,levels=c("NA","FA","FT")),
         phase=factor(phase,levels=c("light","dark")))

p<-ggplot(alpha, aes(x=condition, y=shannon_entropy, fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
  theme_classic()+
  labs(x="condition",y="shannon distance", title="genus|pfam metaT alpha diversity")+
  theme(legend.position = "none")

ggsave("genus_pfam/diversity-core-metrics1.6M/shannon_metaT.pdf", plot=p,height=2.5, width=3.5)

pairwise.wilcox.test(alpha$shannon, alpha$condition,p.adjust.method="fdr")

# NA      FA  
# FA 1.8e-05 -   
#   FT 2.7e-05 0.85

p<-ggplot(alpha, aes(x=condition, y=shannon_entropy, fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
  facet_wrap(~phase)+
  theme_classic()+
  labs(x="condition",y="shannon distance", title="metaT Alpha Diversity")+
  theme(legend.position = "none", plot.title = element_text(face = "bold"))
ggsave("genus_pfam/diversity-core-metrics1.6M/shannon_metaT_lightdark.pdf", plot=p,height=2.5, width=3.5)

alphaL<-alpha%>%filter(phase=="light")
pairwise.wilcox.test(alphaL$shannon, alphaL$condition,
                     p.adjust.method="fdr")

# NA     FA    
# FA 0.0065 -     
#   FT 0.0065 0.4206


alphaD<-alpha%>%filter(phase=="dark")
pairwise.wilcox.test(alphaD$shannon, alphaD$condition,
                     p.adjust.method="fdr")

# NA     FA    
# FA 0.0130 -     
#   FT 0.0065 0.2403
###########################################################

#natlog of BSH vs RPOB
md<-fread("metaT_metadata_ztcat_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  dplyr::rename(`Sample ID`=sample_name)

natlog_bsh<-fread("species_pfam/rpca_results_BSH_RPOB_rmzero/sample_plot_data_speciespfam_BSH_RPOB.tsv")%>%
#natlog_bsh<-fread("genus_pfam/rpca_results_BSH_RPOB_rmzero/sample_plot_data_BSH_RPOB.tsv")%>%
#natlog_bsh<-fread("genus_pfam/rpca_results_rmdoubletons/sample_plot_data_genusBSH_RPOB.tsv")%>%
#natlog_bsh<-fread("species_pfam/rpca_results_rminless5samps/sample_plot_data.tsv")%>%
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
#NA      FA  
# FA 6.9e-07 -   
#   FT 4.6e-08 0.63

#ggsave("genus_pfam/rpca_results_rmdoubletons/SFR23_0630_natlog_genusBSHvsRPOB.pdf", plot=p,height=3.5, width=3.5)
#ggsave("species_pfam/rpca_results_rminless5samps/SFR23_0703_natlog_speciesBSHvsRPOB.pdf", plot=p,height=3.5, width=3.5)
ggsave("species_pfam/rpca_results_BSH_RPOB_rmzero/SFR23_0705_natlog_speciesBSHvsRPOB.pdf", plot=p,height=3.5, width=3.5)
       
p<-ggplot(natlog_bsh, aes(x=condition, y=Current_Natural_Log_Ratio,fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  facet_wrap(~phase)+
  theme_classic()+ scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
  labs(x="condition",y="Natural Log Ratio (BSH/RPOB)",title="Bile Salt Hydrolase (BSH)")+
  theme(legend.position = "none")

#ggsave("genus_pfam/rpca_results_rmdoubletons/SFR23_0630_natlog_genusBSHvsRPOB_lightdark.pdf", plot=p,height=3.5, width=5)
#ggsave("species_pfam/rpca_results_rminless5samps/SFR23_0703_natlog_speciesBSHvsRPOB_lightdark.pdf", plot=p,height=3.5, width=5)
ggsave("species_pfam/rpca_results_BSH_RPOB_rmzero/SFR23_0705_natlog_speciesBSHvsRPOB_lightdark.pdf", plot=p,height=3, width=4)
natlog_bshL<-natlog_bsh%>%filter(phase=="light")
pairwise.wilcox.test(natlog_bshL$Current_Natural_Log_Ratio, natlog_bshL$condition,
                     p.adjust.method="fdr")

#species BSH_RPOB
# NA     FA    
# FA 0.0012 -     
#   FT 0.0012 0.5887

# pairwise.t.test(natlog_bshL$Current_Natural_Log_Ratio, natlog_bshL$condition, p.adjust.method="none",paired=FALSE)
# 
# # NA      FA 
# # FA 2.0e-05 -  
# #   FT 1.1e-05 0.8


natlog_bshD<-natlog_bsh%>%filter(phase=="dark")
pairwise.wilcox.test(natlog_bshD$Current_Natural_Log_Ratio, natlog_bshD$condition,
                     p.adjust.method="fdr")
#species BSH_RPOB
# NA     FA    
# FA 0.0006 -     
#   FT 0.0006 0.0649

# pairwise.t.test(natlog_bshD$Current_Natural_Log_Ratio, natlog_bshD$condition, p.adjust.method="none",paired=FALSE)
# 
# # NA      FA  
# # FA 8.9e-05 -   
# #   FT 3.5e-06 0.17

natlog_bshFT<-natlog_bsh%>%filter(condition=="FT")
pairwise.wilcox.test(natlog_bshFT$Current_Natural_Log_Ratio, natlog_bshFT$phase,
                     p.adjust.method="none")

# light
# dark 0.31 

natlog_bshFA<-natlog_bsh%>%filter(condition=="FA")
pairwise.wilcox.test(natlog_bshFA$Current_Natural_Log_Ratio, natlog_bshFA$phase,
                     p.adjust.method="none")

# light
# dark 0.24 

natlog_bshNA<-natlog_bsh%>%filter(condition=="NA")
pairwise.wilcox.test(natlog_bshNA$Current_Natural_Log_Ratio, natlog_bshNA$phase,
                     p.adjust.method="none")

# light
# dark 0.19 

###########################################################
#natlog but using the songbird results


# diff<-read_qza("genus_pfam/songbird_rmdoubletons_FA/differentials.qza")$data
# diff<-read_qza("genus_pfam/songbird_rmdoubletons/differentials.qza")$data
# diff<-read_qza("genus_pfam/songbird_rmdl_light_FA/differentials.qza")$data
# diff<-read_qza("genus_pfam/songbird_rmdl_dark_FA/differentials.qza")$data
# diff<-read_qza("genus_pfam/songbird_rmdl_light/differentials.qza")$data
# diff<-read_qza("genus_pfam/songbird_rmdl_dark/differentials.qza")$data

diff<-read_qza("genus_pfam/songbird_BSH_rmzero/differentials.qza")$data
diff<-read_qza("genus_pfam/songbird_BSH_rmzero_FA/differentials.qza")$data
diff<-read_qza("genus_pfam/songbird_BSH_rmzero_light_FA/differentials.qza")$data
diff<-read_qza("genus_pfam/songbird_BSH_rmzero_dark_FA/differentials.qza")$data
diff<-read_qza("genus_pfam/songbird_BSH_rmzero_dark/differentials.qza")$data
diff<-read_qza("genus_pfam/songbird_BSH_rmzero_light/differentials.qza")$data

#natlog of BSH vs RPOB
md<-fread("metaT_metadata_ztcat_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  dplyr::rename(`Sample ID`=sample_name)

#natlog_bsh<-fread("genus_pfam/songbird_rmdoubletons/sample_plot_data_sb_BSH_RPOB_NA.tsv")%>%
natlog_bsh<-fread("genus_pfam/songbird_rmdoubletons/sample_plot_data_sb_BSH_RPOB_FA.tsv")%>%
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
                     p.adjust.method="fdr")
# NA      FA  
# FA 6.9e-08 -   
#   FT 6.9e-08 0.76

ggsave("genus_pfam/songbird_rmdoubletons/SFR23_0703_natlog_genusBSHvsRPOB.pdf", plot=p,height=3.5, width=3.5)

p<-ggplot(natlog_bsh, aes(x=condition, y=Current_Natural_Log_Ratio,fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  facet_wrap(~phase)+
  theme_classic()+ scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
  labs(x="condition",y="Natural Log Ratio (BSH/RPOB)",title="Bile Salt Hydrolase (BSH)")+
  theme(legend.position = "none")

ggsave("genus_pfam/songbird_rmdoubletons/SFR23_0703_natlog_genusBSHvsRPOB_lightdark.pdf", plot=p,height=3.5, width=5)

natlog_bshL<-natlog_bsh%>%filter(phase=="light")
pairwise.wilcox.test(natlog_bshL$Current_Natural_Log_Ratio, natlog_bshL$condition,
                     p.adjust.method="fdr")

# NA     FA    
# FA 0.0006 -     
#   FT 0.0006 0.2403

natlog_bshD<-natlog_bsh%>%filter(phase=="dark")
pairwise.wilcox.test(natlog_bshD$Current_Natural_Log_Ratio, natlog_bshD$condition,
                     p.adjust.method="fdr")

# NA     FA    
# FA 0.0006 -     
#   FT 0.0006 0.3095

natlog_bshFT<-natlog_bsh%>%filter(condition=="FT")
pairwise.wilcox.test(natlog_bshFT$Current_Natural_Log_Ratio, natlog_bshFT$phase,
                     p.adjust.method="fdr")

# light
# dark 0.39 

natlog_bshFA<-natlog_bsh%>%filter(condition=="FA")
pairwise.wilcox.test(natlog_bshFA$Current_Natural_Log_Ratio, natlog_bshFA$phase,
                     p.adjust.method="fdr")

# light
# dark 0.065

natlog_bshNA<-natlog_bsh%>%filter(condition=="NA")
pairwise.wilcox.test(natlog_bshNA$Current_Natural_Log_Ratio, natlog_bshNA$phase,
                     p.adjust.method="fdr")

# light
# dark 0.077
###########################################################

#plot the distribution of BSH by genus (Light)

mdT<-fread("metaT_metadata_ztcat_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))

selfeat<-fread("species_pfam/rpca_results_BSH_RPOB_rmzero/selected_features_speciespfam_BSH_RPOB.tsv")%>%
  dplyr::filter(grepl("PF02275.21",`Feature ID`))

bsh_dat<-fread("species_pfam/species_pfam_clean_noNT-TPM.tsv") %>%
  filter(FeatureID %in% selfeat$`Feature ID`)%>%
  mutate(FeatureID=gsub("\\|PF02275.21","",FeatureID))%>%
  gather(sample_name,TPM_counts,-FeatureID)%>%
  left_join(.,mdT,by="sample_name")%>%
  group_by(FeatureID,condition,phase)%>%summarise(mn_TPM=mean(TPM_counts), log_mn_TPM=log10(mean(TPM_counts)+1))%>%
  mutate(condition=factor(condition, levels = c("FT", "FA","NA")),
         phase=factor(phase, levels = c("light","dark")))%>%
  filter(phase=="light")%>%
  #filter(phase=="dark")%>%
  group_by(FeatureID)%>%mutate(sum_TPM=sum(mn_TPM), sum_log_TPM=sum(log_mn_TPM))%>%
  filter(mn_TPM>0,
         !grepl('[[:digit:]]+',FeatureID)),
         FeatureID %in% lst_keep)

lst_keep<-((bsh_dat%>%arrange(-sum_log_TPM))$FeatureID%>%unique())[1:15]

p<-ggplot(data=bsh_dat, aes(x=reorder(FeatureID, sum_log_TPM), y=log_mn_TPM, fill=condition)) +
  geom_bar(stat="identity") + coord_flip() + theme_classic()+
  scale_fill_manual(values=c("#009E73","#D55E00","#0072B2"))+
  #labs(y="log10(mean_TPM)",x="Species",title="Light")+
  labs(y="log10(mean_TPM)",x="Species",title="Dark")+
  theme(legend.position = "none")+
  scale_x_discrete(expand = c(0, 0)) +
  #scale_y_continuous(limits=c(0,70),expand = c(0, 0))
  scale_y_continuous(expand = c(0, 0))

ggsave("species_pfam/SFR23_0719_BSHlight_logTPM.pdf", height=3, width=3.5)
ggsave("species_pfam/SFR23_0719_BSHdark_logTPM.pdf", height=3, width=3.5)

###########################################################
#pull out BSH

md<-fread("metaT_metadata_ztcat_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))

bsh_dat<-fread("genus_pfam/genus_pfam_clean_noNT.tsv") %>%
  dplyr::filter(grepl("PF02275.21",FeatureID))%>%
  gather(sample_name,count,-FeatureID)%>%
  left_join(.,md,by="sample_name")%>%
  mutate(condition=factor(condition, levels = c("NA", "FA","FT")),
         phase=factor(phase, levels = c("light","dark")))%>%
  dplyr::filter(grepl("Kineothrix",FeatureID))

###########################################################
#load the differentials

#FAFTL_diff <- read_qza("genus_pfam/songbird_BSH_rmzero_light_FA/differentials.qza")$data%>%
FAFTL_diff <- read_qza("species_pfam/songbird_BSH_rmzero_light_FA/differentials.qza")$data%>%
  dplyr::rename(feature_ranking=`C(condition, Treatment('FA'))[T.FT]`)%>%
  mutate(condition=ifelse(feature_ranking<0,"FA","FT")) %>%
  mutate(condition=factor(condition,levels=c("FA","FT")),
         phase="light")%>%
  arrange(feature_ranking)

FAFTL_diff_min<-FAFTL_diff%>%slice_min(feature_ranking, n=10)
FAFTL_diff_max<-FAFTL_diff%>%slice_max(feature_ranking, n=10)

combFAFTL<-rbind(FAFTL_diff_min,FAFTL_diff_max)%>%arrange(feature_ranking)

#FAFTD_diff <- read_qza("genus_pfam/songbird_BSH_rmzero_dark_FA/differentials.qza")$data%>%
FAFTD_diff <- read_qza("species_pfam/songbird_BSH_rmzero_dark_FA/differentials.qza")$data%>%
  dplyr::rename(feature_ranking=`C(condition, Treatment('FA'))[T.FT]`)%>%
  mutate(condition=ifelse(feature_ranking<0,"FA","FT")) %>%
  mutate(condition=factor(condition,levels=c("FA","FT")),
         phase="dark")%>%
  filter(Feature.ID %in% combFAFTL$Feature.ID)

combFAFT<-rbind(combFAFTL,FAFTD_diff)%>%arrange(feature_ranking)%>%
  mutate(phase=factor(phase,levels=c("light","dark")))

#plot the differentials

combFAFT$Feature.ID <- factor(combFAFT$Feature.ID,levels = combFAFTL$Feature.ID )

ggplot(data=combFAFT, aes(x=Feature.ID, y=feature_ranking, fill=phase)) +
  geom_bar(stat="identity", position=position_dodge()) +coord_flip() +theme_pubr() +
  scale_fill_manual(values=c("gray70","gray10")) +
  scale_y_continuous(expand=c(0,0))

ggsave("species_pfam/songbird_BSH_rmzero_light_FA/SFR23_0714_BSHsbdiffFTFA_tb15.pdf",height=5, width=6)

FANAL_diff <- read_qza("species_pfam/songbird_BSH_rmzero_light/differentials.qza")$data%>%
  dplyr::rename(feature_ranking=`C(condition, Treatment('NA'))[T.FA]`)%>%
  mutate(condition=ifelse(feature_ranking<0,"NA","FA")) %>%
  mutate(condition=factor(condition,levels=c("NA","FA")),
         phase="light")%>%
  arrange(feature_ranking)

FANAL_diff_min<-FANAL_diff%>%slice_min(feature_ranking, n=10)
FANAL_diff_max<-FANAL_diff%>%slice_max(feature_ranking, n=10)

combFANAL<-rbind(FANAL_diff_min,FANAL_diff_max)%>%arrange(feature_ranking)

FANAD_diff <- read_qza("species_pfam/songbird_BSH_rmzero_dark/differentials.qza")$data%>%
  dplyr::rename(feature_ranking=`C(condition, Treatment('NA'))[T.FA]`)%>%
  mutate(condition=ifelse(feature_ranking<0,"NA","FA")) %>%
  mutate(condition=factor(condition,levels=c("NA","FA")),
         phase="dark")%>%
  filter(Feature.ID %in% combFANAL$Feature.ID)

combFANA<-rbind(combFANAL,FANAD_diff)%>%arrange(feature_ranking)%>%
  mutate(phase=factor(phase,levels=c("light","dark")))

#plot the differentials

combFANA$Feature.ID <- factor(combFANA$Feature.ID,levels = combFANAL$Feature.ID )

ggplot(data=combFANA, aes(x=Feature.ID, y=feature_ranking, fill=phase)) +
  geom_bar(stat="identity", position=position_dodge()) +coord_flip() +theme_pubr() +
  scale_fill_manual(values=c("gray70","gray10")) +
  scale_y_continuous(expand=c(0,0))

ggsave("species_pfam/songbird_BSH_rmzero_light/SFR23_0714_BSHsbdiffFANA_tb15.pdf",height=5, width=6)

p<-bsh_dat%>%
  ggplot(aes(x=condition, y=log10(count+1), fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
  facet_wrap(~phase)+
  theme_classic()+
  labs(x="condition",y="log10(counts)", title="")+
  theme(legend.position = "none")
ggsave(paste("species_func/BSH_lightdark/",FeatureID,".pdf",sep=""), height=3, width=3)

###########################################################
 
#load the birdman results (FA vs FT)
BSHlight<-fread("species_pfam/birdman/species_pfam_BSHonly_clean_rmzero_noNT_light.beta_var.tsv")%>%
  dplyr::rename(ratio=`C(condition, Treatment('FA'))[T.FT]_mean`,
                FeatureID=Feature)%>%
  mutate(`C(condition, Treatment('FA'))[T.FT]_hdi`=gsub("[(]|[)]","",`C(condition, Treatment('FA'))[T.FT]_hdi`))%>%
  separate(`C(condition, Treatment('FA'))[T.FT]_hdi`,c("min","max"), sep=",")%>%
  mutate(min=as.numeric(min),
         max=as.numeric(max),
         diff=max-min,
         phase="light")%>%
  filter(diff<10)%>%
  select(FeatureID, ratio, min, max, phase)%>%
  arrange(ratio)#%>%
  filter(FeatureID %in% BSHdark$FeatureID)


BSHdark<-fread("species_pfam/birdman/species_pfam_BSHonly_clean_rmzero_noNT_dark.beta_var.tsv")%>%
  dplyr::rename(ratio=`C(condition, Treatment('FA'))[T.FT]_mean`,
                FeatureID=Feature)%>%
  mutate(`C(condition, Treatment('FA'))[T.FT]_hdi`=gsub("[(]|[)]","",`C(condition, Treatment('FA'))[T.FT]_hdi`))%>%
  separate(`C(condition, Treatment('FA'))[T.FT]_hdi`,c("min","max"), sep=",")%>%
  mutate(min=as.numeric(min),
         max=as.numeric(max),
         diff=max-min,
         phase="dark")%>%
  select(FeatureID, ratio, min, max, phase)%>%
  filter(FeatureID %in% BSHlight$Feature)%>%
  arrange(ratio)

BSH<-rbind(BSHlight,BSHdark)%>%mutate(phase=factor(phase,levels=c("light","dark")))

BSH$Feature <- factor(BSH$Feature,levels = BSHlight$Feature )
BSH$Feature <- factor(BSH$Feature,levels = BSHdark$Feature )

ggplot(BSH, aes(x =Feature , y = ratio, ymin = min, ymax = max, color=phase, group=phase)) + 
  geom_linerange(position = position_dodge(width = 0.8)) + theme_pubr()+
  geom_pointrange(position = position_dodge(width = 0.8)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  coord_flip() +scale_color_manual(values=c("gray70","gray10"))
#ggsave("species_pfam/birdman/FAFTLD_full.pdf",height=10, width=6)
#ggsave("species_pfam/birdman/FAFTLD_smhdi.pdf",height=5, width=6)
ggsave("species_pfam/birdman/FAFTLD_darr_smhdi.pdf",height=5, width=6)
ggsave("species_pfam/birdman/FAFTLD_darr_smhdiL.pdf",height=5, width=6)

BSH<-fread("species_pfam/birdman/species_pfam_BSHonly_clean_rmzero_noNT.beta_var.tsv")%>%
  dplyr::rename(ratio=`C(condition, Treatment('FA'))[T.FT]_mean`,
                FeatureID=Feature)%>%
  mutate(`C(condition, Treatment('FA'))[T.FT]_hdi`=gsub("[(]|[)]","",`C(condition, Treatment('FA'))[T.FT]_hdi`))%>%
  separate(`C(condition, Treatment('FA'))[T.FT]_hdi`,c("min","max"), sep=",")%>%
  mutate(min=as.numeric(min),
         max=as.numeric(max),
         phase="light",
         cred=ifelse(min>0|max<0,"credible","not_credible"))%>%
  filter(cred=="credible")%>%
  select(FeatureID, ratio, min, max, phase)%>%
  arrange(ratio)

BSH$Feature <- factor(BSH$Feature,levels = BSH$Feature )

ggplot(BSH, aes(x =Feature , y = ratio, ymin = min, ymax = max)) + 
  geom_linerange(position = position_dodge(width = 0.8)) + theme_pubr()+
  geom_pointrange(position = position_dodge(width = 0.8)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  coord_flip()

#ggsave("species_pfam/birdman/FAFT_full.pdf",height=10, width=6)
ggsave("species_pfam/birdman/FAFT_justcred.pdf",height=2, width=6)

###########################################################

#load the birdman results (NA vs. FA)
BSHlight<-fread("species_pfam/birdman/species_pfam_BSHonly_clean_rmzero_noNT_lightNA.beta_var.tsv")%>%
  dplyr::rename(ratio=`C(condition, Treatment('NA'))[T.FA]_mean`,
                FeatureID=Feature)%>%
  mutate(`C(condition, Treatment('NA'))[T.FA]_hdi`=gsub("[(]|[)]","",`C(condition, Treatment('NA'))[T.FA]_hdi`))%>%
  separate(`C(condition, Treatment('NA'))[T.FA]_hdi`,c("min","max"), sep=",")%>%
  mutate(min=as.numeric(min),
         max=as.numeric(max),
         cred=ifelse(min>0|max<0,"credible","not_credible"),
         phase="light",
         diff=max-min)%>%
  #filter(diff<10)%>%
  select(FeatureID, ratio, min, max, phase)%>%
  #arrange(ratio)#%>%
  filter(FeatureID %in% BSHdark$Feature)


BSHdark<-fread("species_pfam/birdman/species_pfam_BSHonly_clean_rmzero_noNT_darkNA.beta_var.tsv")%>%
  dplyr::rename(ratio=`C(condition, Treatment('NA'))[T.FA]_mean`,
                FeatureID=Feature)%>%
  mutate(`C(condition, Treatment('NA'))[T.FA]_hdi`=gsub("[(]|[)]","",`C(condition, Treatment('NA'))[T.FA]_hdi`))%>%
  separate(`C(condition, Treatment('NA'))[T.FA]_hdi`,c("min","max"), sep=",")%>%
  mutate(min=as.numeric(min),
         max=as.numeric(max),
         cred=ifelse(min>0|max<0,"credible","not_credible"),
         phase="dark",
         diff=max-min)%>%
  filter(diff<10)%>%
  select(FeatureID, ratio, min, max, phase)%>%
  #filter(FeatureID %in% BSHlight$Feature)
  arrange(ratio)

BSH<-rbind(BSHlight,BSHdark)%>%mutate(phase=factor(phase,levels=c("light","dark")))

#BSH$Feature <- factor(BSH$Feature,levels = BSHlight$Feature )
BSH$Feature <- factor(BSH$Feature,levels = BSHdark$Feature )

ggplot(BSH, aes(x =Feature , y = ratio, ymin = min, ymax = max, color=phase, group=phase)) + 
  geom_linerange(position = position_dodge(width = 0.8)) + theme_pubr()+
  geom_pointrange(position = position_dodge(width = 0.8)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  coord_flip() +scale_color_manual(values=c("gray70","gray10"))
#ggsave("species_pfam/birdman/FANALD_full.pdf",height=10, width=6)
#ggsave("species_pfam/birdman/FANALD_smhdi.pdf",height=5, width=6)
ggsave("species_pfam/birdman/FANALD_darr_smhdi.pdf",height=5, width=6)


BSH<-fread("species_pfam/birdman/species_pfam_BSHonly_clean_rmzero_noNT_NA.beta_var.tsv")%>%
  dplyr::rename(ratio=`C(condition, Treatment('NA'))[T.FA]_mean`,
                FeatureID=Feature)%>%
  mutate(`C(condition, Treatment('NA'))[T.FA]_hdi`=gsub("[(]|[)]","",`C(condition, Treatment('NA'))[T.FA]_hdi`))%>%
  separate(`C(condition, Treatment('NA'))[T.FA]_hdi`,c("min","max"), sep=",")%>%
  mutate(min=as.numeric(min),
         max=as.numeric(max),
         phase="light",
         cred=ifelse(min>0|max<0,"credible","not_credible"))%>%
  filter(cred=="credible")%>%
  select(FeatureID, ratio, min, max, phase)%>%
  arrange(ratio)

BSH$Feature <- factor(BSH$Feature,levels = BSH$Feature )

ggplot(BSH, aes(x =Feature , y = ratio, ymin = min, ymax = max)) + 
  geom_linerange(position = position_dodge(width = 0.8)) + theme_pubr()+
  geom_pointrange(position = position_dodge(width = 0.8)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  coord_flip()

#ggsave("species_pfam/birdman/FANA_full.pdf",height=10, width=6)
ggsave("species_pfam/birdman/FANA_justcred.pdf",height=4, width=6)
