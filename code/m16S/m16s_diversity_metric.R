setwd("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/micro16s")

library(tidyverse)
library(data.table)
library("qiime2R")
library("Biostrings")
library(ggrepel)
library(ggpubr)
library(ggbreak)
library(scatterplot3d)
library(car)

md<-fread("metadata.TRF_combined.tab")%>%
  mutate(phase=ifelse(zt<13,"light","dark"),
         cond_phase=paste(condition,phase,sep="_"),
         fasted=ifelse(zt %in% c(5,9,13), "fasted","fed"),
         cond_fasted=paste(condition,fasted,sep="_"))
write.table(md, "metadata.TRF_combined_wLD.tab",sep = "\t",row.names = FALSE, quote=FALSE)
write.table(md, "metadata.TRF_combined_wLD_wfasting.tab",sep = "\t",row.names = FALSE, quote=FALSE)


###shannon-rarefied
md<-fread("metadata.TRF_combined_wLD.tab")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))

alpha <- read_qza("diversity-core-metrics6k/shannon_vector.qza")$data %>%
  rownames_to_column("sample_name") %>%
  left_join(.,md,by="sample_name") %>%
  filter(condition!="NT")%>%
  mutate(condition=factor(condition,levels=c("NA","FA","FT")),
         lightdark=factor(phase,levels=c("light","dark")))

p<-ggplot(alpha, aes(x=condition, y=shannon_entropy, fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
  theme_classic()+
  labs(x="condition",y="shannon distance", title="16S alpha diversity")+
  theme(legend.position = "none")

ggsave("diversity-core-metrics6k/SFR23_0602_shannon_m16s.pdf", plot=p,height=2.5, width=3.5)

pairwise.wilcox.test(alpha$shannon, alpha$condition,p.adjust.method="fdr")

# NA   FA  
# FA 0.89 -   
#   FT 0.89 0.56

p<-ggplot(alpha, aes(x=condition, y=shannon_entropy, fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
  facet_wrap(~phase)+
  theme_classic()+
  labs(x="condition",y="shannon distance", title="16S Alpha Diversity")+
  theme(legend.position = "none", plot.title = element_text(face = "bold"))
ggsave("diversity-core-metrics6k/SFR23_0602_shannon_m16s_lightdark.pdf", plot=p,height=2.5, width=3.5)

alphaL<-alpha%>%filter(phase=="light")
pairwise.wilcox.test(alphaL$shannon, alphaL$condition,
                     p.adjust.method="fdr")
# NA   FA  
# FA 0.33 -   
#   FT 0.29 0.29

alphaD<-alpha%>%filter(phase=="dark")
pairwise.wilcox.test(alphaD$shannon, alphaD$condition,
                     p.adjust.method="fdr")

# NA   FA  
# FA 0.41 -   
#   FT 0.41 0.69

alphaFT<-alpha%>%filter(condition=="FT")
pairwise.wilcox.test(alphaFT$shannon, alphaFT$phase,
                     p.adjust.method="fdr")
# light
# dark 0.85

alphaFA<-alpha%>%filter(condition=="FA")
pairwise.wilcox.test(alphaFA$shannon, alphaFA$phase,
                     p.adjust.method="fdr")
# light
# dark 0.81

alphaNA<-alpha%>%filter(condition=="NA")
pairwise.wilcox.test(alphaNA$shannon, alphaNA$phase,
                     p.adjust.method="fdr")
# light
# dark 0.046

alphaZT<-alpha%>% group_by(condition,zt)%>%mutate(meanAbun=mean(shannon_entropy),str=sd(shannon_entropy)/sqrt(length(shannon_entropy)),
                                                       t.score = qt(p=0.05/2, df=length(shannon_entropy)-1,lower.tail=F),
                                                       margin.error=t.score*str)

# p<-alphaZT%>%
#   ggplot(aes(x=zt, y=meanAbun, color=condition)) +
#   geom_point(alpha=1.0) + geom_line() +
#   theme_pubr() +
#   scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+
#   scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
#   geom_ribbon(aes(ymin = meanAbun-margin.error, ymax = meanAbun+margin.error, fill=condition),alpha=0.5,colour = NA)+
#   labs(color="condition",
#        y =paste("shannon distance",sep=""),
#        x ="ZT time")+ theme(plot.title = element_text(face = "bold"))
# ggsave("diversity-core-metrics6k/SFR22_1028_shannon_m16s_ztime.pdf", plot=p,height=4, width=5)

#RPCA results
ord <- read_qza("ordination_deicode_asv16S.qza")

samp_ord<-ord$data$Vectors
write.table(samp_ord,"ordination_deicode_asv16s/sample_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)
feat_ord<-ord$data$Species
write.table(feat_ord,"ordination_deicode_asv16s/feature_ordination.txt",sep = "\t",row.names = FALSE,quote=FALSE)

md<-fread("metadata.TRF_combined_wLD.tab")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))

rpca<-ord$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2, PC3)%>%
  dplyr::rename(sample_name=SampleID)%>%
  left_join(md,by="sample_name")%>%
  filter(condition!="NT")%>%
  mutate(condition=factor(condition,levels=c("FT","FA","NA")),
         phase=factor(phase,levels=c("light","dark")))

pdf("SFR23_0602_m16s_RPCA.pdf",width = 5, height = 5)
pdf("SFR24_0411_m16s_RPCA.pdf",width = 5, height = 5)

# Create a 3D PCA plot
colors <- c("#009E73","#D55E00","#0072B2")
colors <- colors[as.numeric(rpca$condition)]

scatterplot3d(rpca[,2:4],angle = 50,pch = 16,color=colors,grid=FALSE, size=5,
              xlab = "PC1", ylab = "PC2", zlab = "PC3",
              main = "16S")
legend("bottom", legend = levels(rpca$condition),
       col =  c("#009E73","#D55E00","#0072B2"), pch = 16, 
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
       y =paste("PC2 (",round(ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""))+ggtitle("16S RPCA")+ theme(plot.title = element_text(face = "bold"))
ggsave("SFR22_1028_m16s_RPCA.pdf", plot=p,height=3.5, width=3.5)


#plot the ZT time per condition
NA_rpca<-rpca%>%
  filter(condition=="NA")%>%
  mutate(zt=factor(zt,levels=c("1","5","9","13","17","21")))

p<-NA_rpca %>%
  ggplot(aes(x=PC1, y=PC2, fill=zt)) +
  geom_point(alpha=1.0,size=3,shape=21) + 
  theme_pubr() +
  scale_fill_manual(values=c("#67001f","#d6604d","#ffba92","#8ccff3","#4393c3","#053061"))+
  labs(color="ZT time",
       x =paste("PC1 (",round(ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""))+ggtitle("16S NA")+ theme(plot.title = element_text(face = "bold"))

ggsave("SFR24_0412_m16s_NA_RPCA.pdf", plot=p,height=3.5, width=3.5)

FA_rpca<-rpca%>%
  filter(condition=="FA")%>%
  mutate(zt=factor(zt,levels=c("1","5","9","13","17","21")))

p<-FA_rpca %>%
  ggplot(aes(x=PC1, y=PC2, fill=zt)) +
  geom_point(alpha=1.0,size=3, shape=21) + 
  theme_pubr() +
  scale_fill_manual(values=c("#67001f","#d6604d","#ffba92","#8ccff3","#4393c3","#053061"))+
  labs(color="ZT time",
       x =paste("PC1 (",round(ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""))+ggtitle("16S FA")+ theme(plot.title = element_text(face = "bold"))

ggsave("SFR24_0412_m16s_FA_RPCA.pdf", plot=p,height=3.5, width=3.5)

FT_rpca<-rpca%>%
  filter(condition=="FT")%>%
  mutate(zt=factor(zt,levels=c("1","5","9","13","17","21")))

p<-FT_rpca %>%
  ggplot(aes(x=PC1, y=PC2, fill=zt)) +
  geom_point(alpha=1.0,size=3,shape=21) + 
  theme_pubr() +
  scale_fill_manual(values=c("#67001f","#d6604d","#ffba92","#8ccff3","#4393c3","#053061"))+
  labs(color="ZT time",
       x =paste("PC1 (",round(ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""))+ggtitle("16S FT")+ theme(plot.title = element_text(face = "bold"))

ggsave("SFR24_0412_m16s_FT_RPCA.pdf", plot=p,height=3.5, width=3.5)


