setwd("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metagenomic/woltka2_results")

library(tidyverse)
library(data.table)
library("Biostrings")
library(ggrepel)
library(ggpubr)
library(ggbreak)
library(ALDEx2)
library(RColorBrewer)
library(qiime2R)

###########################################################

md<-fread("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metagenomic/metaG_metadata.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  filter(condition!="NT")%>%
  mutate(cond_phase=paste(condition,lightdark,sep="_"),
         cond_zt=paste(condition,zt_time,sep="_"),
         fasted=ifelse(zt_time %in% c(5,9,13), "fasted","fed"),
         cond_fasted=paste(condition,fasted,sep="_"))
#write.table(md, "/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metagenomic/metaG_metadata_noNT.txt",sep = "\t",row.names = FALSE, quote=FALSE)
write.table(md, "/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metagenomic/metaG_metadata_noNT_wfasting.txt",sep = "\t",row.names = FALSE, quote=FALSE)


dat<-fread("filtered_metaG/go/function.tsv") %>%dplyr::select(c(1:43))
names(dat) <- c("FeatureID","FA1b","FA1c","FA2b","FA2c","FA3b","FA3c","FA4b","FA4c","FA5b","FA5c",
                         "FA6b","FA6c","FT1b","FT1c","FT2b","FT2c","FT3b","FT3c","FT4b","FT4c","FT5b",
                         "FT5c","FT6b","FT6c","NA1a","NA1b","NA1c","NA2a","NA2b","NA2c","NA3a","NA3b","NA3c",
                         "NA4a","NA4b","NA4c","NA5a","NA5b","NA5c","NA6a","NA6b","NA6c")
write.table(dat,"filtered_metaG/go/function_clean_noNT.txt",sep = "\t",row.names = FALSE, quote=FALSE)

###########################################################

# dat<-fread("filtered_metaG/genus_pfam.tsv") %>%dplyr::select(c(1:43))
# names(dat) <- c("FeatureID","FA1b","FA1c","FA2b","FA2c","FA3b","FA3c","FA4b","FA4c","FA5b","FA5c",
#                 "FA6b","FA6c","FT1b","FT1c","FT2b","FT2c","FT3b","FT3c","FT4b","FT4c","FT5b",
#                 "FT5c","FT6b","FT6c","NA1a","NA1b","NA1c","NA2a","NA2b","NA2c","NA3a","NA3b","NA3c",
#                 "NA4a","NA4b","NA4c","NA5a","NA5b","NA5c","NA6a","NA6b","NA6c")
# write.table(dat,"filtered_metaG/genus_pfam_clean_noNT.txt",sep = "\t",row.names = FALSE, quote=FALSE)
###########################################################

#RPCA
##for all the transcripts removing singletons/doubletons
ord <- read_qza("filtered_metaG/pfam_notnorm/metaG_rpca_results/ordination.qza")

samp_ord<-ord$data$Vectors
write.table(samp_ord,"filtered_metaG/pfam_notnorm/metaG_rpca_results/sample_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)
feat_ord<-ord$data$Species
write.table(feat_ord,"filtered_metaG/pfam_notnorm/metaG_rpca_results/feature_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)

md<-fread("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metagenomic/metaG_metadata_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))

rpca<-ord$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2,PC3)%>%
  dplyr::rename(sample_name=SampleID)%>%
  left_join(md,by="sample_name")%>%
  mutate(condition=factor(condition,levels=c("FT","FA","NA")),
         phase=factor(lightdark,levels=c("light","dark")))


pdf("filtered_metaG/pfam_notnorm/metaG_rpca_results/SFR23_0517_metaG_RPCA.pdf",width = 5, height = 5)

# Create a 3D PCA plot
colors <- c("#0072B2","#D55E00","#009E73")
colors <- colors[as.numeric(rpca$condition)]

scatterplot3d(rpca[,2:4],angle = 230,pch = 16,color=colors,grid=FALSE,
              xlab = "PC1", ylab = "PC2", zlab = "PC3",
              main = "MGX")
legend("bottom", legend = levels(rpca$condition),
       col =  c("#0072B2","#D55E00","#009E73"), pch = 16, 
       inset = -0.25, xpd = TRUE, horiz = TRUE)

dev.off()

p<-rpca %>%
  ggplot(aes(x=PC1, y=PC2, color=condition, shape=phase)) +
  geom_point(alpha=1.0) + #alpha controls transparency and helps when points are overlapping
  theme_pubr() +stat_ellipse(type = "t", linetype = 2,aes(group = condition))+
  scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+
  scale_shape_manual(values=c(3,16)) +
  labs(color="condition",
       x =paste("PC1 (",round(ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""))+ggtitle("metaG RPCA")+ theme(plot.title = element_text(face = "bold"))
ggsave("filtered_metaG/pfam_notnorm/metaG_rpca_results/SFR23_0327_metaG_RPCA.pdf", plot=p,height=3.5, width=3.5)

#plot the ZT time per condition
NA_rpca<-rpca%>%
  filter(condition=="NA")%>%
  mutate(zt_time=factor(zt_time,levels=c("1","5","9","13","17","21")))

p<-NA_rpca %>%
  ggplot(aes(x=PC1, y=PC2, fill=zt_time)) +
  geom_point(alpha=1.0,size=3,shape=21) + 
  theme_pubr() +
  scale_fill_manual(values=c("#67001f","#d6604d","#ffba92","#8ccff3","#4393c3","#053061"))+
  labs(color="ZT time",
       x =paste("PC1 (",round(ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""))+ggtitle("MGX NA")+ theme(plot.title = element_text(face = "bold"))

ggsave("filtered_metaG/pfam_notnorm/metaG_rpca_results/SFR24_0412_mgx_NA_RPCA.pdf", plot=p,height=3.5, width=3.5)

FA_rpca<-rpca%>%
  filter(condition=="FA")%>%
  mutate(zt_time=factor(zt_time,levels=c("1","5","9","13","17","21")))

p<-FA_rpca %>%
  ggplot(aes(x=PC1, y=PC2, fill=zt_time)) +
  geom_point(alpha=1.0,size=3,shape=21) + 
  theme_pubr() +
  scale_fill_manual(values=c("#67001f","#d6604d","#ffba92","#8ccff3","#4393c3","#053061"))+
  labs(color="ZT time",
       x =paste("PC1 (",round(ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""))+ggtitle("MGX FA")+ theme(plot.title = element_text(face = "bold"))

ggsave("filtered_metaG/pfam_notnorm/metaG_rpca_results/SFR24_0412_mgx_FA_RPCA.pdf", plot=p,height=3.5, width=3.5)

FT_rpca<-rpca%>%
  filter(condition=="FT")%>%
  mutate(zt_time=factor(zt_time,levels=c("1","5","9","13","17","21")))

p<-FT_rpca %>%
  ggplot(aes(x=PC1, y=PC2, fill=zt_time)) +
  geom_point(alpha=1.0,size=3,shape=21) + 
  theme_pubr() +
  scale_fill_manual(values=c("#67001f","#d6604d","#ffba92","#8ccff3","#4393c3","#053061"))+
  labs(color="ZT time",
       x =paste("PC1 (",round(ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""))+ggtitle("MGX FT")+ theme(plot.title = element_text(face = "bold"))

ggsave("filtered_metaG/pfam_notnorm/metaG_rpca_results/SFR24_0412_mgx_FT_RPCA.pdf", plot=p,height=3.5, width=3.5)

###########################################################

#alpha diversity -with rarefied table

md<-fread("~/scratch/metatranscript/metagenomic/metaG_metadata_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))

alpha <- read_qza("filtered_metaG/pfam/diversity-core-metrics12k/shannon_vector.qza")$data %>%
  rownames_to_column("sample_name") %>%
  left_join(.,md,by="sample_name") %>%
  mutate(condition=factor(condition,levels=c("NA","FA","FT")),
         lightdark=factor(lightdark,levels=c("light","dark")))%>%
  filter(condition!="NT")

p<-ggplot(alpha, aes(x=condition, y=shannon_entropy, fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  scale_fill_manual(values=c("#009E73","#D55E00","#0072B2"))+
  theme_classic()+
  labs(x="condition",y="shannon distance", title="metaG alpha diversity")+
  theme(legend.position = "none")

ggsave("filtered_metaG/pfam/diversity-core-metrics12k/shannon_metaG.pdf", plot=p,height=2.5, width=3.5)

pairwise.wilcox.test(alpha$shannon, alpha$condition,p.adjust.method="fdr")

# NA   FA  
# FA 0.26 -   
#   FT 0.26 0.77

p<-ggplot(alpha, aes(x=condition, y=shannon_entropy, fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  scale_fill_manual(values=c("#009E73","#D55E00","#0072B2"))+
  facet_wrap(~lightdark)+
  theme_classic()+
  labs(x="condition",y="shannon distance", title="metaG Alpha Diversity")+
  theme(legend.position = "none", plot.title = element_text(face = "bold"))
ggsave("filtered_metaG/pfam/diversity-core-metrics12k/shannon_metaG_lightdark.pdf", plot=p,height=2.5, width=3.5)

alphaL<-alpha%>%filter(lightdark=="light")
pairwise.wilcox.test(alphaL$shannon, alphaL$condition,
                     p.adjust.method="fdr")
# NA   FA  
# FA 0.88 -   
#   FT 0.88 0.88

alphaD<-alpha%>%filter(lightdark=="dark")
pairwise.wilcox.test(alphaD$shannon, alphaD$condition,
                     p.adjust.method="fdr")

# NA  FA 
# FA 0.1 -  
#   FT 0.1 1.0

alphaFT<-alpha%>%filter(condition=="FT")
pairwise.wilcox.test(alphaFT$shannon, alphaFT$lightdark,
                     p.adjust.method="fdr")
# light
# dark 0.69

alphaFA<-alpha%>%filter(condition=="FA")
pairwise.wilcox.test(alphaFA$shannon, alphaFA$lightdark,
                     p.adjust.method="fdr")
# light
# dark 0.18 

alphaNA<-alpha%>%filter(condition=="NA")
pairwise.wilcox.test(alphaNA$shannon, alphaNA$lightdark,
                     p.adjust.method="fdr")
# light
# dark 0.016

#rarefied visual

#plotting the rarefaction
# rarefac<-fread("filtered_metaG/pfam/alpha-rarefaction20k/shannon.csv")%>%
#   gather(key,values,-diet,-schedule,-replicate,-condition,-lightdark,-organ,-sample.id) %>%
#   mutate(condition=ifelse(is.na(condition),"NA",condition)) %>%
#   separate(key,c("depth","iteration"),sep='_') %>%
#   mutate_at("depth", str_replace, "depth.", "") %>%
#   mutate_at("iteration", str_replace, "iter.", "") %>%
#   filter(condition!="NT")%>%
#   mutate(depth=as.numeric(depth),
#          iteration=as.numeric(iteration),
#          condition=factor(condition,levels=c("FT","FA","NA")))
# 
# rfac_sub<-rarefac%>% group_by(condition,depth)%>%mutate(meanAbun=mean(values,na.rm=TRUE),str=sd(values/sqrt(length(values)),na.rm=TRUE),
#                                                         t.score = qt(p=0.05/2, df=length(values)-1,lower.tail=F),
#                                                         margin.error=t.score*str)
# 
# p<-rfac_sub%>%
#   ggplot(aes(x=depth, y=meanAbun, color=condition)) +
#   geom_point(alpha=1.0) + geom_line() +
#   theme_pubr() +
#   scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+
#   scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
#   geom_ribbon(aes(ymin = meanAbun-margin.error, ymax = meanAbun+margin.error, fill=condition),alpha=0.5,colour = NA)+
#   labs(color="condition",
#        y ="shannon distance",
#        x ="sequencing depth")
# ggsave("alpha-rarefaction20k/SFR22_0824_metaG_rarefaction.pdf", plot=p,height=4, width=5)
# 
# rarS<-fread("alpha-rarefaction20k/rarefaction_sampdrop_metaG.txt")%>%
#   mutate(condition=ifelse(is.na(condition),"NA",condition)) %>%
#   mutate(condition=factor(condition,levels=c("FT","FA","NA")))
# 
# p<-rarS%>%
#   ggplot(aes(x=depth, y=num_samps, color=condition)) +
#   geom_point(alpha=1.0) + geom_line() +
#   theme_pubr() +
#   scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+
#   scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
#   labs(color="condition",
#        y ="number of samples",
#        x ="sequencing depth")
# ggsave("alpha-rarefaction20k/SFR22_0824_metaG_rarefaction_samps.pdf", plot=p,height=4, width=5)
