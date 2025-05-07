setwd("~/Notebooks/sfloresr/TRF-metaT/")

library(tidyverse)
library(data.table)
library(ggrepel)
library(ggpubr)
library(RColorBrewer)
library(ComplexHeatmap)
###########################################################
#paths
dat_library<-"data/ENB_invitro_metab/GNPS/GNPS_outputs/networking/library-results-merged_results_with_gnps_justBA_collapse.tsv"
dat_feattab<-"data/ENB_invitro_metab/GNPS/GNPS_outputs/clustering/featuretable_reformated_justBA_summcollapse_cln.csv"
dat_metadata<-"data/ENB_invitro_metab/metadata/ZT_BSH_metadata.tsv"
dat_path<-"data/ENB_invitro_metab/"
fig_path<-"figures/ENB_invitro_metab/"
###########################################################
#functions

TRF_unpr_ttest<-function(mtb){
  x<-pairwise.t.test(mtb$log_10_peakabun,mtb$timepoint,p.adjust.method = "fdr")
  df<-data.frame(t0hvt24h=x$p.value[1],
                 t0hvt48h=x$p.value[2],
                 t24hvt48h=x$p.value[4])
  return(df)
}

run_ttests_BA<-function(mtb,BA){
  stat.test<-mtb%>%
    group_by(BSH_strain,FeatureID)%>%
    nest()%>%
    mutate(pvals = map(data, ~ TRF_unpr_ttest(.x)))%>%
    dplyr::select(-data)%>%
    unnest()
  write.table(stat.test,paste(dat_path,"GNPS/ttest_unpr_results/",BA,"suppl_ttest_unprd_tmpt_pvals.txt",sep=""),
              sep = "\t",row.names = FALSE, quote=FALSE)  
}

get_signif_hits_allBA<-function(mtb,BA){
  mtb_summ<-mtb%>%
    gather(comparison,pval,-FeatureID,-BSH_strain)%>%
    filter(comparison=="t0hvt48h"&pval<0.05)%>%
    mutate(BA_suppl=BA)
  return(mtb_summ)
}
###########################################################
#load files

annot<-fread(dat_library)%>%
  mutate(Compound_Name = gsub("\\|.*", "", Compound_Name))%>%
  mutate(Compound_Name = gsub("\\(delta mass.*", "", Compound_Name))

mtb<-fread(dat_feattab)

md<-fread(dat_metadata)%>%
  filter(sample_name %in% colnames(mtb))%>%
  mutate(name=gsub("\\.mzML$","",filename))%>%
  mutate(name = sapply(strsplit(name, "_"), function(x) paste(x[3:length(x)], collapse = "_")))%>%
  dplyr::select(1:2,8,3:7)
write.table(md,paste0(dat_path,"metadata/ZT_BSH_metadata_justBA.tsv"),sep = "\t",row.names = FALSE, quote=FALSE)

###########################################################
#clean files

md_sub<-md%>%
  filter(timepoint!="")

mtb_sub<-mtb%>%
  dplyr::rename(FeatureID=`row.ID`)%>%
  dplyr::select(FeatureID, all_of(md_sub$sample_name))%>%
  gather(sample_name,peak_abun,-FeatureID)%>%
  left_join(.,annot,by="FeatureID")%>%
  left_join(.,md,by="sample_name")%>%
  filter(!(BA_suppl=="BHI-50DMSO"|BA_suppl=="BHI"))%>%
  mutate(log_10_peakabun=log10(peak_abun+1))
###########################################################
#run t-test

BA_list<-c("GCA","GCDCA","GDCA","GLCA","GUDCA","TCA","TCDCA","TDCA","TLCA","TUDCA")
for(i in BA_list){
  mtb_new<-mtb_sub%>%filter(BA_suppl==i)
  run_ttests_BA(mtb_new,i)
}

BA_t0ht48h_df<-data.frame(FeatureID=NA,BSH_strain=NA,comparison=NA,pval=NA,BA_suppl=NA)
for(i in BA_list){
  df_BA<-fread(paste(dat_path,"GNPS/ttest_unpr_results/",i,"suppl_ttest_unprd_tmpt_pvals.txt",sep=""))
  x<-get_signif_hits_allBA(df_BA,i)
  BA_t0ht48h_df<-rbind(BA_t0ht48h_df,x)
}

#summary of results for all BAs
BA_t0ht48h_df<-BA_t0ht48h_df%>%filter(!is.na(FeatureID))
write.table(BA_t0ht48h_df,paste0(dat_path,"GNPS/ttest_unpr_results/summ_hits_t0ht48h_allBAstrainsuppl_padj0.05.txt"),
            sep = "\t",row.names = FALSE, quote=FALSE)   

#filter out hits that were diff under control
list_BA_ctrl<-BA_t0ht48h_df%>%
  filter(BSH_strain=="ctrl"|BSH_strain=="EcAZ-1-cat")

BA_t0ht48h_rmctrl<-BA_t0ht48h_df%>%
  filter(BSH_strain!="ctrl")%>%
  filter(!(FeatureID %in% list_BA_ctrl$FeatureID))%>%
  mutate(unique_id=paste(FeatureID,BA_suppl,sep="_"))
#35 unique mtb diff that not diff in controls

###########################################################
#make heatmap of these 35 unique mtbs--Figure S5A

mtb_plt<-mtb%>%
  dplyr::rename(FeatureID=`row.ID`)%>%
  dplyr::select(FeatureID, all_of(md_sub$sample_name))%>%
  gather(sample_name,peak_abun,-FeatureID)%>%
  left_join(.,md,by="sample_name")%>%
  mutate(unique_id=paste(FeatureID,BA_suppl,sep="_"))%>%
  filter(unique_id %in% unique(BA_t0ht48h_rmctrl$unique_id))%>%
  group_by(FeatureID,name,BA_suppl,BSH_strain,timepoint)%>%
  summarise(mn_peak_abun=mean(peak_abun))%>%
  mutate(log_10_peakabun=log10(mn_peak_abun+1))%>%
  group_by(FeatureID,BA_suppl)%>%mutate(Zscore=(log_10_peakabun - mean(log_10_peakabun))/sd(log_10_peakabun))%>%
  left_join(.,annot,by="FeatureID")%>%
  mutate(label_name=paste(FeatureID,Compound_Name,sep=" "))%>%
  mutate(timepoint=factor(timepoint, levels = c("0h","24h","48h")),
         BSH_strain=factor(BSH_strain,levels=c("EcAZ-1-cat","AZ-52","Dny-BSH1","Dny-BSH2","Lgasseri-BSH",
                                               "LCAG-95-BSH1208","Ep-BSH101","25MeOH","50MeOH","ctrl")),
         BA_suppl=factor(BA_suppl,levels=c("BHI-50DMSO","BHI","GCA","GCDCA","GDCA","GLCA","GUDCA","TCA",
                                           "TCDCA","TDCA","TLCA","TUDCA")))%>%
  filter(!(BA_suppl=="BHI-50DMSO"|BA_suppl=="BHI"))%>%
  filter(timepoint!="24h")

plt<-ggplot(mtb_plt,aes(x=label_name, y=fct_rev(BSH_strain))) +theme_classic()+
  geom_tile(aes(fill=Zscore))+
  scale_x_discrete(expand = c(0, 0))+
  facet_grid(timepoint~BA_suppl,scales="free",space="free")+
  theme(axis.ticks.y=element_blank(),panel.spacing.x=unit(0.3, "lines"),axis.text.y = element_text(size = 8),
        panel.spacing.y=unit(0.3, "lines"),axis.text.x = element_text(angle=90,hjust = 1))+
  scale_fill_distiller(palette = "Spectral", direction = -1)+
  labs(y="Bile acids",x="BSH strains")

ggsave(paste0(fig_path,"ttest_unpr_results/SFR24_0501_sighitst0vt48_heatmap_zscore.pdf"), plt,width = 10, height = 6)
