setwd("~/Notebooks/sfloresr/TRF-metaT/")

library(tidyverse)
library(data.table)
library("qiime2R")
library(ggpubr)
library(rstatix)

##########################################################
#paths
m16s_metadata<-"data/metadata.TRF_combined_wLD.tab"
mgx_metadata<-"data/metaG_metadata_noNT.txt"
mtx_metadata<-"data/metaT_metadata_ztcat_noNT.txt"
a_m16s<-"data/diversity_analysis/diversity-core-metrics6k_16s/shannon_vector.qza"
a_mgx<-"data/diversity_analysis/diversity-core-metrics12k_metaG/shannon_vector.qza"
a_mtx<-"data/diversity_analysis/diversity-core-metrics2.3M_metaT/shannon_vector.qza"
res_path<-"data/diversity_analysis/"
fig_path<-"figures/diversity_analysis/"
##########################################################
#functions

get_alpha_dat<-function(dt,md){
  alpha <- read_qza(dt)$data %>%
    rownames_to_column("sample_name") %>%
    left_join(.,md,by="sample_name") %>%
    mutate(condition=factor(condition,levels=c("NA","FA","FT")),
           phase=factor(phase,levels=c("light","dark")),
           method="MGX")%>%
    filter(condition!="NT")%>%
    dplyr::select(sample_name,condition,phase,zt_time,method,shannon_entropy)
  return(alpha)
}

stderror <- function(x) sd(x)/sqrt(length(x))


get_posthoc_t<-function(dt){
  pwc <- dt %>%
    group_by(zt_time)%>%
    pairwise_t_test(
      shannon_entropy ~ condition, paired = FALSE,
      p.adjust.method = "fdr")
  return(pwc)
}
##########################################################
#load shannon

md16s<-fread(m16s_metadata)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  dplyr::rename(zt_time=zt)

mdG<-fread(mgx_metadata)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  dplyr::rename(phase=lightdark)

mdT<-fread(mtx_metadata)%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))

alpha16s <- get_alpha_dat(a_m16s,md16s)
alphaG <- get_alpha_dat(a_mgx,mdG)
alphaT <- get_alpha_dat(a_mtx,mdT)

##########################################################
#plot by condition--Figure S1C

alpha<-rbind(alpha16s,alphaG,alphaT)

p<-ggplot(alpha, aes(x=condition, y=shannon_entropy, fill=condition)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  facet_wrap(~method,scales = "free_y")+
  scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
  theme_pubr()+
  labs(x="condition",y="shannon distance", title="16S alpha diversity")+
  theme(legend.position = "none")

ggsave(paste0(fig_path,"SFR24_0412_shannon_m16sMGXMTG_pfam.pdf"), plot=p,height=3, width=8)

pairwise.wilcox.test(alpha16s$shannon, alpha16s$condition,p.adjust.method="fdr")

# NA   FA  
# FA 0.89 -   
#   FT 0.89 0.56

pairwise.wilcox.test(alphaG$shannon, alphaG$condition,p.adjust.method="fdr")

# NA   FA  
# FA 0.55 -   
#   FT 0.55 0.55

pairwise.wilcox.test(alphaT$shannon, alphaT$condition,p.adjust.method="fdr")

# NA      FA  
# FA 3.4e-06 -   
#   FT 6.9e-08 0.22

##########################################################
#plot over ZT--Figure S1D

alpha_summ<-alpha%>%
  group_by(method,condition,zt_time)%>%
  summarise(mn_shannon_entropy=mean(shannon_entropy),sem=stderror(shannon_entropy))

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

ggsave(paste0(fig_path,"SFR24_0412_shannon_m16sMGXMTG_overZT_pfam.pdf"), plot=p,height=2.5, width=8)

#stats
alpha16s<-alpha16s%>%mutate(zt_time=as.factor(zt_time))

res_aov <- aov(shannon_entropy ~ condition*zt_time, data = alpha16s)
summary(res_aov)

# Df Sum Sq Mean Sq F value Pr(>F)  
# condition          2   0.29  0.1449   0.230 0.7961  
# zt_time            5   8.27  1.6541   2.623 0.0441 *
#   condition:zt_time 10  13.25  1.3251   2.101 0.0568 .
# Residuals         30  18.92  0.6307    

pwc<-get_posthoc_t(alpha16s)
write.table(pwc,paste0(res_path,"diversity-core-metrics6k_16s/SFR24_0412_shannon_m16soverZT_pval.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

res_aov <- aov(shannon_entropy ~ condition*zt_time, data = alphaG)
summary(res_aov)

# Df Sum Sq Mean Sq F value  Pr(>F)   
# condition          2  0.272  0.1362   1.384 0.26429   
# zt_time            1  0.780  0.7804   7.930 0.00804 **
#   condition:zt_time  2  0.494  0.2472   2.512 0.09605 . 
# Residuals         34  3.346  0.0984 

pwc<-get_posthoc_t(alphaG)
write.table(pwc,paste0(res_path,"diversity-core-metrics12k_metaG/SFR24_0412_shannon_mgxoverZT_pval.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

res_aov <- aov(shannon_entropy ~ condition*zt_time, data = alphaT)
summary(res_aov)

# Df Sum Sq Mean Sq F value  Pr(>F)    
# condition          2 15.575   7.787  21.825 6.2e-07 ***
#   zt_time            1  0.029   0.029   0.081   0.777    
# condition:zt_time  2  1.374   0.687   1.925   0.161    
# Residuals         36 12.845   0.357  

pwc<-get_posthoc_t(alphaT)
write.table(pwc,paste0(res_path,"diversity-core-metrics2.3M_metaT/SFR24_0412_shannon_mtxoverZT_pval.txt"),sep = "\t",row.names = FALSE, quote=FALSE)

