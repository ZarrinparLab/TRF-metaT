setwd("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/multiomics/")

library(tidyverse)
library(data.table)
library("qiime2R")
library("Biostrings")
library(ggrepel)
library(ggpubr)
library(ggbreak)
library(ggh4x)
library(rstatix)
library(multcomp)
#load 16S shannon

md16s<-fread("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/micro16s/metadata.TRF_combined_wLD.tab")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))

alpha16s <- read_qza("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/micro16s/diversity-core-metrics6k/shannon_vector.qza")$data %>%
  rownames_to_column("sample_name") %>%
  left_join(.,md16s,by="sample_name") %>%
  filter(condition!="NT")%>%
  dplyr::rename(zt_time=zt)%>%
  mutate(condition=factor(condition,levels=c("NA","FA","FT")),
         phase=factor(phase,levels=c("light","dark")),
         method="16S")%>%
  dplyr::select(sample_name,condition,phase,zt_time,method,shannon_entropy)

pairwise.wilcox.test(alpha16s$shannon, alpha16s$condition,p.adjust.method="fdr")

# NA   FA  
# FA 0.89 -   
#   FT 0.89 0.56

#load mgx shannon

mdG<-fread("~/scratch/TRF_multiomics/metagenomic/metaG_metadata_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  dplyr::rename(phase=lightdark)

alphaG <- read_qza("~/scratch/TRF_multiomics/metagenomic/woltka2_results/filtered_metaG/pfam_notnorm/diversity-core-metrics12k/shannon_vector.qza")$data %>%
  rownames_to_column("sample_name") %>%
  left_join(.,mdG,by="sample_name") %>%
  mutate(condition=factor(condition,levels=c("NA","FA","FT")),
         phase=factor(phase,levels=c("light","dark")),
         method="MGX")%>%
  filter(condition!="NT")%>%
  dplyr::select(sample_name,condition,phase,zt_time,method,shannon_entropy)

pairwise.wilcox.test(alphaG$shannon, alphaG$condition,p.adjust.method="fdr")

# NA   FA  
# FA 0.26 -   
#   FT 0.26 0.77

#load mtx shannon
mdT<-fread("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_m_results/metaT_metadata_ztcat_noNT.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))

alphaT <- read_qza("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_m_results/pfam/diversity-core-metrics2.3M/shannon_vector.qza")$data %>%
  rownames_to_column("sample_name") %>%
  left_join(.,mdT,by="sample_name") %>%
  mutate(condition=factor(condition,levels=c("NA","FA","FT")),
         phase=factor(phase,levels=c("light","dark")),
         method="MTX")%>%
  dplyr::select(sample_name,condition,phase,zt_time,method,shannon_entropy)

pairwise.wilcox.test(alphaT$shannon, alphaT$condition,p.adjust.method="fdr")

# NA      FA  
# FA 3.4e-06 -   
#   FT 6.9e-08 0.22

alpha<-rbind(alpha16s,alphaG,alphaT)

p<-ggplot(alpha, aes(x=condition, y=shannon_entropy, fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  facet_wrap(~method,scales = "free_y")+
  scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
  theme_pubr()+
  labs(x="condition",y="shannon distance", title="16S alpha diversity")+
  theme(legend.position = "none")

#ggsave("SFR23_0619_shannon_m16sMGXMTG_pfam.pdf", plot=p,height=4, width=6)
ggsave("SFR24_0412_shannon_m16sMGXMTG_pfam.pdf", plot=p,height=3, width=8)

stderror <- function(x) sd(x)/sqrt(length(x))

alpha_summ<-alpha%>%
  group_by(method,condition,zt_time)%>%
  summarise(mn_shannon_entropy=mean(shannon_entropy),sem=stderror(shannon_entropy))
#over time
p<-ggplot(alpha_summ, aes(x=zt_time, y=mn_shannon_entropy, color=condition)) +
  geom_point(alpha=1.0) + geom_line() +
  theme_pubr() +
  scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+
  scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
  scale_x_continuous(breaks=c(1,5,9,13,17,21))+
  geom_ribbon(aes(ymin = mn_shannon_entropy-sem, ymax = mn_shannon_entropy+sem, fill=condition),alpha=0.3,colour = NA)+
  facet_wrap(~method,scales = "free_y")+
  labs(x="condition",y="shannon distance")+
  theme(legend.position = "none")
ggsave("SFR24_0412_shannon_m16sMGXMTG_overZT_pfam.pdf", plot=p,height=2.5, width=8)


#stats
alpha16s<-alpha16s%>%mutate(zt_time=as.factor(zt_time))
res_aov <- aov(shannon_entropy ~ condition*zt_time,
               data = alpha16s)
summary(res_aov)

# Df Sum Sq Mean Sq F value Pr(>F)  
# condition          2   0.29  0.1449   0.230 0.7961  
# zt_time            5   8.27  1.6541   2.623 0.0441 *
#   condition:zt_time 10  13.25  1.3251   2.101 0.0568 .
# Residuals         30  18.92  0.6307    

pwc <- alpha16s %>%
  group_by(zt_time)%>%
  pairwise_t_test(
    shannon_entropy ~ condition, paired = FALSE,
    p.adjust.method = "fdr"
  )

write.table(pwc,"SFR24_0412_shannon_m16soverZT_pval.txt",sep = "\t",row.names = FALSE, quote=FALSE)

res_aov <- aov(shannon_entropy ~ condition*zt_time,
               data = alphaG)
summary(res_aov)

# Df Sum Sq Mean Sq F value  Pr(>F)   
# condition          2  0.272  0.1362   1.384 0.26429   
# zt_time            1  0.780  0.7804   7.930 0.00804 **
#   condition:zt_time  2  0.494  0.2472   2.512 0.09605 . 
# Residuals         34  3.346  0.0984 

pwc <- alphaG %>%
  group_by(zt_time)%>%
  pairwise_t_test(
    shannon_entropy ~ condition, paired = FALSE,
    p.adjust.method = "fdr"
  )
write.table(pwc,"SFR24_0412_shannon_mgxoverZT_pval.txt",sep = "\t",row.names = FALSE, quote=FALSE)


res_aov <- aov(shannon_entropy ~ condition*zt_time,
               data = alphaT)
summary(res_aov)

# Df Sum Sq Mean Sq F value  Pr(>F)    
# condition          2 15.575   7.787  21.825 6.2e-07 ***
#   zt_time            1  0.029   0.029   0.081   0.777    
# condition:zt_time  2  1.374   0.687   1.925   0.161    
# Residuals         36 12.845   0.357  

post_test <- glht(res_aov,
                  linfct = mcp(condition = "Tukey")
)

summary(post_test)

pwc <- alphaT %>%
  group_by(zt_time)%>%
  pairwise_t_test(
    shannon_entropy ~ condition, paired = FALSE,
    p.adjust.method = "fdr"
  )
write.table(pwc,"SFR24_0412_shannon_mtxoverZT_pval.txt",sep = "\t",row.names = FALSE, quote=FALSE)

