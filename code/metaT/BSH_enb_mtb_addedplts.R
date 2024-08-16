setwd("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_m_results/")

library(tidyverse)
library(data.table)
library("qiime2R")
library("Biostrings")
library(ggrepel)
library(ggpubr)
library(gplots)
library(ggvenn)
library(viridis)
library(pheatmap)
library(ggpubfigs)
library(rstatix)
library(UpSetR)
library(ComplexHeatmap)
library(RColorBrewer)
###########################################################
annot<-fread("BSH/BSH_ENB/GNPS/nf_output/networking/library-results-merged_results_with_gnps_justBA_collapse.tsv")%>%
  mutate(Compound_Name = gsub("\\|.*", "", Compound_Name))%>%
  mutate(Compound_Name = gsub("\\(delta mass.*", "", Compound_Name))

mtb<-fread("BSH/BSH_ENB/GNPS/nf_output/clustering/featuretable_reformated_justBA_summcollapse_cln.csv")

md<-fread("BSH/BSH_ENB/GNPS/metadata_filename/ZT_BSH_metadata.tsv")%>%
  filter(sample_name %in% colnames(mtb))%>%
  mutate(name=gsub("\\.mzML$","",filename))%>%
  mutate(name = sapply(strsplit(name, "_"), function(x) paste(x[3:length(x)], collapse = "_")))%>%
  dplyr::select(1:2,8,3:7)
write.table(md,"BSH/BSH_ENB/GNPS/metadata_filename/ZT_BSH_metadata_justBA.tsv",sep = "\t",row.names = FALSE, quote=FALSE)  

# list_venn <- list(BSH_tab = colnames(mtb),
#                   BSH_md = md$sample_name)
# ItemsList <- venn(list_venn, show.plot = FALSE)
# all<-attributes(ItemsList)$intersections


summ_plate<-md%>%
  group_by(plate)%>%
  summarise(n=n())

plate_nest<-md%>%
  separate(plate_location, into = c("column", "row"), sep = "(?<=[A-Za-z])(?=[0-9])", remove = FALSE)%>%
  mutate(row=as.numeric(row))%>%
  group_by(plate)%>%
  nest()
###########################################################

#heatmap of all hits w/replicates 
md_sub<-md%>%
  filter(timepoint!="")

mtb_sub<-mtb%>%
  dplyr::rename(FeatureID=`row.ID`)%>%
  dplyr::select(FeatureID, all_of(md_sub$sample_name))%>%
  gather(sample_name,peak_abun,-FeatureID)%>%
  left_join(.,annot,by="FeatureID")%>%
  left_join(.,md,by="sample_name")%>%
  group_by(FeatureID,name,sample_name,BA_suppl,BSH_strain,timepoint)%>%
  summarise(mn_peak_abun=mean(peak_abun))%>%
  mutate(log_10_peakabun=log10(mn_peak_abun+1),
         timepoint=factor(timepoint, levels = c("0h","24h","48h")),
         BSH_strain=factor(BSH_strain,levels=c("EcAZ-1-cat","AZ-52","Dny-BSH1","Dny-BSH2","Lgasseri-BSH",
                                               "LCAG-95-BSH1208","Ep-BSH101","25MeOH","50MeOH","ctrl")),
         BA_suppl=factor(BA_suppl,levels=c("BHI-50DMSO","BHI","GCA","GCDCA","GDCA","GLCA","GUDCA","TCA",
                                           "TCDCA","TDCA","TLCA","TUDCA")))%>%
  filter(!(BA_suppl=="BHI-50DMSO"|BA_suppl=="BHI"))
#group_by(label_name,condition)%>%mutate(Zscore=(mn_TPM - mean(mn_TPM))/sd(mn_TPM))%>%

plt<-ggplot(mtb_sub,aes(x=FeatureID, y=BSH_strain)) +theme_classic()+
  geom_tile(aes(fill=log_10_peakabun))+
  scale_x_discrete(expand = c(0, 0))+
  facet_grid(timepoint~BA_suppl,scales="free",space="free")+
  theme(axis.ticks.y=element_blank(),panel.spacing.x=unit(0.3, "lines"),axis.text.y = element_text(size = 8),
        panel.spacing.y=unit(0.3, "lines"),strip.text.y = element_text(angle = 0),axis.text.x = element_blank())+
  scale_fill_distiller(palette = "Spectral", direction = -1)

ggsave("BSH/BSH_ENB/GNPS/SFR24_429_heatmap_all.pdf", plt,width = 20, height = 8)
ggsave("BSH/BSH_ENB/GNPS/SFR24_429_heatmap_all_rmBHI.pdf", plt,width = 18, height = 6)

#function to get a heatmap for each BA supplementation

df_forheatmap<-function(mtb,md,annot,BA_added){
  md_sub<-md%>%
    filter(BA_suppl==BA_added)
  mtb_sub<-mtb%>%
    dplyr::select(`row.ID`, all_of(md_sub$sample_name))%>%
    gather(sample_name,peak_abun,-row.ID)%>%
    left_join(.,md_sub,by="sample_name")%>%
    dplyr::rename(FeatureID=row.ID)%>%
    left_join(.,annot,by="FeatureID")%>%
    mutate(label_name=paste(FeatureID,Compound_Name,sep=" "))%>%
    group_by(name,label_name)%>%
    summarise(mn_peak_abun=mean(peak_abun))%>%
    spread(label_name,mn_peak_abun)%>%
    column_to_rownames("name")%>%as.matrix()%>%
    log10()
  mtb_sub[mtb_sub== -Inf] <- 0
  return(mtb_sub)
}

BHI<-df_forheatmap(mtb,md,annot,"BHI")
p<-pheatmap(BHI,cluster_rows=FALSE)
ggsave("BSH/BSH_ENB/GNPS/SFR24_429_heatmap_BHI.pdf", p,width = 15, height = 6)

BHIDMSO<-df_forheatmap(mtb,md,annot,"BHI-50DMSO")
p<-pheatmap(BHIDMSO,cluster_rows=FALSE)
ggsave("BSH/BSH_ENB/GNPS/SFR24_429_heatmap_BHI-50DMSO.pdf", p,width = 15, height = 6)

GCA<-df_forheatmap(mtb,md,annot,"GCA")
p<-pheatmap(GCA,cluster_rows=FALSE)
ggsave("BSH/BSH_ENB/GNPS/SFR24_429_heatmap_GCA.pdf", p,width = 15, height = 8)

GCDCA<-df_forheatmap(mtb,md,annot,"GCDCA")
p<-pheatmap(GCDCA,cluster_rows=FALSE)
ggsave("BSH/BSH_ENB/GNPS/SFR24_429_heatmap_GCDCA.pdf", p,width = 15, height = 8)

GDCA<-df_forheatmap(mtb,md,annot,"GDCA")
p<-pheatmap(GDCA,cluster_rows=FALSE)
ggsave("BSH/BSH_ENB/GNPS/SFR24_429_heatmap_GDCA.pdf", p,width = 15, height = 8)

GLCA<-df_forheatmap(mtb,md,annot,"GLCA")
p<-pheatmap(GLCA,cluster_rows=FALSE)
ggsave("BSH/BSH_ENB/GNPS/SFR24_429_heatmap_GLCA.pdf", p,width = 15, height = 8)

GUDCA<-df_forheatmap(mtb,md,annot,"GUDCA")
p<-pheatmap(GUDCA,cluster_rows=FALSE)
ggsave("BSH/BSH_ENB/GNPS/SFR24_429_heatmap_GUDCA.pdf", p,width = 15, height = 8)

TCA<-df_forheatmap(mtb,md,annot,"TCA")
p<-pheatmap(TCA,cluster_rows=FALSE)
ggsave("BSH/BSH_ENB/GNPS/SFR24_429_heatmap_TCA.pdf", p,width = 15, height = 8)

TCDCA<-df_forheatmap(mtb,md,annot,"TCDCA")
p<-pheatmap(TCDCA,cluster_rows=FALSE)
ggsave("BSH/BSH_ENB/GNPS/SFR24_429_heatmap_TCDCA.pdf", p,width = 15, height = 8)

TDCA<-df_forheatmap(mtb,md,annot,"TDCA")
p<-pheatmap(TDCA,cluster_rows=FALSE)
ggsave("BSH/BSH_ENB/GNPS/SFR24_429_heatmap_TDCA.pdf", p,width = 15, height = 8)

TLCA<-df_forheatmap(mtb,md,annot,"TLCA")
p<-pheatmap(TLCA,cluster_rows=FALSE)
ggsave("BSH/BSH_ENB/GNPS/SFR24_429_heatmap_TLCA.pdf", p,width = 15, height = 8)

TUDCA<-df_forheatmap(mtb,md,annot,"TUDCA")
p<-pheatmap(TUDCA,cluster_rows=FALSE)
ggsave("BSH/BSH_ENB/GNPS/SFR24_429_heatmap_TUDCA.pdf", p,width = 15, height = 8)


mtb_sub<-mtb%>%
  dplyr::rename(FeatureID=`row.ID`)%>%
  dplyr::select(FeatureID, all_of(md_sub$sample_name))%>%
  gather(sample_name,peak_abun,-FeatureID)%>%
  left_join(.,annot,by="FeatureID")%>%
  left_join(.,md,by="sample_name")%>%
  filter(!(BA_suppl=="BHI-50DMSO"|BA_suppl=="BHI"))%>%
  mutate(log_10_peakabun=log10(peak_abun+1))

#this shows we are missing some of our replicates, prob from prelim run
# summ_trip<-md_sub%>%
#   group_by(BA_suppl,BSH_strain,timepoint)%>%
#   summarise(n=n())
# 
# summ_plates<-md%>%
#   group_by(plate)%>%
#   summarise(n=n())

# TRF_unpr_wilcox<-function(mtb){
#   x<-pairwise.wilcox.test(mtb$peak_abun,mtb$timepoint,p.adjust.method = "none")
#   df<-data.frame(t0hvt24h=x$p.value[1],
#                  t0hvt48h=x$p.value[2],
#                  t24hvt48h=x$p.value[4])
#   return(df)
# }

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
  write.table(stat.test,paste("BSH/BSH_ENB/GNPS/ttest_unpr_results/",BA,"suppl_ttest_unprd_tmpt_pvals.txt",sep=""),sep = "\t",row.names = FALSE, quote=FALSE)  
}

BA_list<-c("GCA","GCDCA","GDCA","GLCA","GUDCA","TCA","TCDCA","TDCA","TLCA","TUDCA")
for(i in BA_list){
  mtb_new<-mtb_sub%>%filter(BA_suppl==i)
  run_ttests_BA(mtb_new,i)
}

# mtb_ttests <- mtb_sub%>% 
#   group_by(BA_suppl) %>% 
#   nest()%>%
#   mutate(pval = map2(data, BA_suppl,  ~ run_ttests_BA(.x,.y)))


#load the sig BAs from each BA suppl 0v48h

get_upset_data<-function(mtb){
  list_df<-mtb%>%
    gather(comparison,pval,-FeatureID,-BSH_strain)%>%
    filter(comparison=="t0hvt48h"&pval<0.05)%>%
    group_by(BSH_strain)%>%nest()
  return(list_df)
}

BA_list<-c("GCA","GCDCA","GDCA","GLCA","GUDCA","TCA","TCDCA","TDCA","TLCA","TUDCA")

mtb_upset<-fread(paste("BSH/BSH_ENB/GNPS/ttest_unpr_results/",BA_list[10],"suppl_ttest_unprd_tmpt_pvals.txt"))
list_df<-get_upset_data(mtb_upset)
list_venn <- list(AZ52 = list_df[[2]][[1]]$FeatureID,
                  #ctrl = list_df[[2]][[2]]$FeatureID,
                  DnyBSH1 = list_df[[2]][[2]]$FeatureID,
                  DnyBSH2 = list_df[[2]][[3]]$FeatureID,
                  EcAZ1cat = list_df[[2]][[4]]$FeatureID,
                  EpBSH101 = list_df[[2]][[5]]$FeatureID,
                  LCAG95BSH1208 = list_df[[2]][[6]]$FeatureID,
                  LgasseriBSH = list_df[[2]][[7]]$FeatureID)
m = make_comb_mat(list_venn)
pdf(file=paste("BSH/BSH_ENB/GNPS/ttest_unpr_results/",BA_list[10],"suppl_ttest_unprd_tmpt_upset.pdf",sep=""),width=8,height=4)
UpSet(m,set_order=c("AZ52","DnyBSH1","DnyBSH2","LgasseriBSH","EpBSH101","LCAG95BSH1208","EcAZ1cat"))#,,"ctrl"
dev.off()

#make a summary fig of num of mtb with time diff t0hvt48h

make_summplt_BA<-function(mtb,BA){
  mtb_summ<-mtb%>%
    gather(comparison,pval,-FeatureID,-BSH_strain)%>%
    filter(comparison=="t0hvt48h"&pval<0.05)%>%
    group_by(BSH_strain)%>%
    summarise(n=n())%>%
    mutate(BA_suppl=BA)
  return(mtb_summ)
  }

BA_summ_df<-data.frame(BSH_strain=NA,n=NA,BA_suppl=NA)
for(i in BA_list){
  df_BA<-fread(paste("BSH/BSH_ENB/GNPS/ttest_unpr_results/",i,"suppl_ttest_unprd_tmpt_pvals.txt",sep=""))
  x<-make_summplt_BA(df_BA,i)
  BA_summ_df<-rbind(BA_summ_df,x)
}

BA_summ_plt<-  BA_summ_df%>%filter(!is.na(n))%>%
  mutate(BSH_strain=factor(BSH_strain,levels=c("EcAZ-1-cat","AZ-52","Dny-BSH1","Dny-BSH2","Lgasseri-BSH",
                                    "LCAG-95-BSH1208","Ep-BSH101","ctrl")),
         BA_suppl=factor(BA_suppl,levels=c("GCA","GCDCA","GDCA","GLCA","GUDCA","TCA",
                                  "TCDCA","TDCA","TLCA","TUDCA")))

p<-ggplot(data=BA_summ_plt, aes(x=fct_reorder(BA_suppl,n), y=n, fill=BSH_strain)) + 
  geom_bar(stat="identity") + theme_classic()+
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_brewer(palette = "Paired")+
  coord_flip()+
  labs(y="num. of BA with changes over time (padj<0.05)",
       x="BA Supplemented to culture")+
  theme(legend.key.size = unit(0.5, 'cm'))

ggsave("BSH/BSH_ENB/GNPS/ttest_unpr_results/SFR24_0501_summofBAHwdifft0t48h_plt_BAord.pdf", p,width = 6, height = 4)

nb.cols <- 10
mycolors <- colorRampPalette(brewer.pal(8, "Paired"))(nb.cols)

p<-ggplot(data=BA_summ_plt, aes(x=fct_reorder(BSH_strain,n), y=n, fill=BA_suppl)) + 
  geom_bar(stat="identity") + theme_classic()+
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_manual(values = mycolors) +
  coord_flip()+
  labs(y="num. of BA with changes over time (padj<0.05)",
       x="BSH strains")+
  theme(legend.key.size = unit(0.5, 'cm'))

ggsave("BSH/BSH_ENB/GNPS/ttest_unpr_results/SFR24_0501_summofBAHwdifft0t48h_plt_strainord.pdf", p,width = 6, height = 4)

get_signif_hits_allBA<-function(mtb,BA){
  mtb_summ<-mtb%>%
    gather(comparison,pval,-FeatureID,-BSH_strain)%>%
    filter(comparison=="t0hvt48h"&pval<0.05)%>%
    mutate(BA_suppl=BA)
  return(mtb_summ)
}

BA_t0ht48h_df<-data.frame(FeatureID=NA,BSH_strain=NA,comparison=NA,pval=NA,BA_suppl=NA)
for(i in BA_list){
  df_BA<-fread(paste("BSH/BSH_ENB/GNPS/ttest_unpr_results/",i,"suppl_ttest_unprd_tmpt_pvals.txt",sep=""))
  x<-get_signif_hits_allBA(df_BA,i)
  BA_t0ht48h_df<-rbind(BA_t0ht48h_df,x)
}

BA_t0ht48h_df<-  BA_t0ht48h_df%>%filter(!is.na(FeatureID))
write.table(BA_t0ht48h_df,"BSH/BSH_ENB/GNPS/ttest_unpr_results/summ_hits_t0ht48h_allBAstrainsuppl_padj0.05.txt",sep = "\t",row.names = FALSE, quote=FALSE)   

#filter out hits that were diff under control
BA_t0ht48h_df<-fread("BSH/BSH_ENB/GNPS/ttest_unpr_results/summ_hits_t0ht48h_allBAstrainsuppl_padj0.05.txt")
list_BA_ctrl<-BA_t0ht48h_df%>%
  filter(BSH_strain=="ctrl"|BSH_strain=="EcAZ-1-cat")

BA_t0ht48h_rmctrl<-BA_t0ht48h_df%>% #35 unique mtb diff that not diff in controls 
  filter(BSH_strain!="ctrl")%>%
  filter(!(FeatureID %in% list_BA_ctrl$FeatureID))%>%
  mutate(unique_id=paste(FeatureID,BA_suppl,sep="_"))
#make heatmap of these mtbs

#heatmap of these hits w/replicates 
md_sub<-md%>%
  filter(timepoint!="")

mtb_sub<-mtb%>%
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


plt<-ggplot(mtb_sub,aes(x=label_name, y=fct_rev(BSH_strain))) +theme_classic()+
  #geom_tile(aes(fill=log_10_peakabun))+
  geom_tile(aes(fill=Zscore))+
  scale_x_discrete(expand = c(0, 0))+
  facet_grid(timepoint~BA_suppl,scales="free",space="free")+
  theme(axis.ticks.y=element_blank(),panel.spacing.x=unit(0.3, "lines"),axis.text.y = element_text(size = 8),
        panel.spacing.y=unit(0.3, "lines"),axis.text.x = element_text(angle=90,hjust = 1))+
  scale_fill_distiller(palette = "Spectral", direction = -1)+
  labs(y="Bile acids",x="BSH strains")

ggsave("BSH/BSH_ENB/GNPS/ttest_unpr_results/SFR24_0501_sighitst0vt48_heatmap.pdf", plt,width = 10, height = 6)
ggsave("BSH/BSH_ENB/GNPS/ttest_unpr_results/SFR24_0501_sighitst0vt48_heatmap_zscore.pdf", plt,width = 10, height = 6)
ggsave("BSH/BSH_ENB/GNPS/ttest_unpr_results/SFR24_0501_sighitst0vt48_heatmap_zscore_rmnotsig.pdf", plt,width = 10, height = 6)

#these hits individually
show_deconj<-function(mtb,mtb_id){
  BA_candi_x<-mtb%>%filter(BA_suppl=="GCA")
  p <- ggplot(BA_candi_x, aes(x=fct_rev(BSH_strain),y=log10(baseline_chng+0.000001),colour=condition)) +
    geom_boxplot()+
    #geom_point()+
    facet_wrap(~BA_suppl)+ coord_flip()+
    geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
    scale_colour_manual(values=c("#0072B2","#D55E00","#009E73"))+
    theme_pubr()+labs(title=unique(BA_candi_x$cmpd_name_new))+
    theme(legend.position = "none")
  
  ggsave(paste("BSH/BSH_ENB/GNPS/change_baseline/indv_mtb_plt/SFR24_0603_mtb",mtb_id,"GCA.pdf",sep="_"), p,width = 4, height = 3)
}

BA_suppl_subset<-BA_t0ht48h_rmctrl%>%filter(BA_suppl=="GCA")
mtb_sub<-mtb%>%
  mutate(condition=case_when(BSH_strain=="Dny-BSH1"|BSH_strain=="Dny-BSH2"~"FT",
                             BSH_strain=="LCAG-95-BSH1208"|BSH_strain=="Ep-BSH101"~"NA",
                             BSH_strain=="Lgasseri-BSH"~"FA",
                              .default = "none"))%>%
  mutate(baseline_chng=ifelse(`0h`+`48h`==0, NA,baseline_chng))%>%
  mutate(cmpd_name_new=paste(FeatureID,Compound_Name,sep="_"))%>%
  filter(FeatureID %in% BA_suppl_subset$FeatureID)%>%
  mutate(BSH_strain=factor(BSH_strain,levels=c("EcAZ-1-cat","AZ-52","Dny-BSH1","Dny-BSH2","Lgasseri-BSH",
                                        "LCAG-95-BSH1208","Ep-BSH101","25MeOH","50MeOH","ctrl")),
         condition=factor(condition,levels=c("NA","FA","FT")))
  

indv_plot <- mtb_sub %>%
  group_by(FeatureID) %>%
  nest() %>%
  mutate(pval = map2(data, FeatureID,  ~ show_deconj(.x,.y)))


#plot three examples of conjugated bile acids

mtb_sub<-mtb%>%
  dplyr::rename(FeatureID=`row.ID`)%>%
  dplyr::select(FeatureID, all_of(md_sub$sample_name))%>%
  gather(sample_name,peak_abun,-FeatureID)%>%
  left_join(.,md,by="sample_name")%>%
  mutate(unique_id=paste(FeatureID,BA_suppl,sep="_"))%>%
  filter(unique_id %in% unique(BA_t0ht48h_rmctrl$unique_id))%>%
  left_join(.,annot,by="FeatureID")%>%
  mutate(log_10_peakabun=log10(peak_abun+1),
         timepoint=factor(timepoint, levels = c("0h","24h","48h")),
         BSH_strain=factor(BSH_strain,levels=c("EcAZ-1-cat","AZ-52","Dny-BSH1","Dny-BSH2","Lgasseri-BSH",
                                               "LCAG-95-BSH1208","Ep-BSH101","25MeOH","50MeOH","ctrl")),
         BA_suppl=factor(BA_suppl,levels=c("BHI-50DMSO","BHI","GCA","GCDCA","GDCA","GLCA","GUDCA","TCA",
                                           "TCDCA","TDCA","TLCA","TUDCA")))%>%
  filter(!(BA_suppl=="BHI-50DMSO"|BA_suppl=="BHI"))%>%
  filter(timepoint!="24h")


plt_sig_mtbs<-function(mtb,mtb_id){
  cmp_name<-unique(mtb$Compound_Name)
  p<-ggplot(data=mtb, aes(x=timepoint, y=log_10_peakabun, colour=BA_suppl, fill=BA_suppl)) +
    geom_boxplot(color = "black",alpha=0.3)+
    geom_point(shape=21, color = "black", position=position_jitterdodge())+
    facet_wrap(~BSH_strain, nrow=1)+
    theme_classic()+
    labs(title=paste(cmp_name," (",mtb_id,")",sep=""),
         x="timepoint", y="log10(peak_area+1)")+
    theme(legend.position = "right", plot.title = element_text(size = 12),axis.title.x = element_text(size = 10))
  ggsave(paste("BSH/BSH_ENB/GNPS/ttest_unpr_results/indv_mtb_plt/SFR24_0501_mtb",mtb_id,"_notpaired.pdf",sep=""), p,width = 8, height = 2)
  
}

indv_plot<-mtb_sub%>%
  group_by(FeatureID)%>%
  nest()%>%
  mutate(pval = map2(data, FeatureID,  ~ plt_sig_mtbs(.x,.y)))

###########################################################
#get change from baseline
annot<-fread("BSH/BSH_ENB/GNPS/nf_output/networking/library-results-merged_results_with_gnps_justBA_collapse.tsv")%>%
  mutate(Compound_Name = gsub("\\|.*", "", Compound_Name))%>%
  mutate(Compound_Name = gsub("\\(delta mass.*", "", Compound_Name))%>%
  dplyr::select(FeatureID,Compound_Name)

md<-fread("BSH/BSH_ENB/GNPS/metadata_filename/ZT_BSH_metadata_justBA.tsv")%>%
  mutate(cult_plt=case_when(plate=="P4"|plate=="P5"|plate=="P6"~"P2",
                            plate=="P7"|plate=="P8"|plate=="P9"~"P3",
                            .default = "P1"))%>%
  filter(timepoint!="24h")
  
mtb<-fread("BSH/BSH_ENB/GNPS/nf_output/clustering/featuretable_reformated_justBA_summcollapse_cln.csv")%>%
  dplyr::rename(FeatureID=`row.ID`)%>%
  dplyr::select(FeatureID,any_of(md$sample_name))%>%
  gather(sample_name,peak_area,-FeatureID)%>%
  left_join(.,md,by="sample_name")%>%
  mutate(name_new=paste(BA_suppl,BSH_strain,sep="-"))%>%
  dplyr::select(FeatureID,name_new,BA_suppl,BSH_strain,cult_plt,plate_location,timepoint,peak_area)%>%
  filter(timepoint!=""& !(BA_suppl%in% c("BHI","BHI-50DMSO")))%>%
  pivot_wider(names_from = timepoint, values_from = peak_area)%>%
  #mutate(baseline_chng=`48h`/(`0h`+1))%>%
  mutate(baseline_chng=`48h`-`0h`)%>%
  right_join(.,annot,by="FeatureID")%>%
  filter(`0h`+`48h`>0)#%>%
  # filter(baseline_chng>0)

p <- ggplot(mtb, aes(x=log10(baseline_chng+0.000001), colour=BSH_strain)) + 
  geom_density()+ facet_wrap(~BA_suppl, nrow=2)+
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")+
  theme_pubr()+labs(title="all BAs")

ggsave("BSH/BSH_ENB/GNPS/change_baseline/allBA_density.pdf", p,width = 10, height = 4)

BA_amine<-annot%>% filter(str_detect(Compound_Name, "[^- ]-[^- ]"))

BA_candi<-mtb%>%filter(!(FeatureID %in% BA_amine$FeatureID))%>%
  mutate(condition=case_when(BSH_strain=="Dny-BSH1"|BSH_strain=="Dny-BSH2"~"FT",
                                      BSH_strain=="LCAG-95-BSH1208"|BSH_strain=="Ep-BSH101"~"NA",
                                      BSH_strain=="Lgasseri-BSH"~"FA",
                                      .default = "none"),
         BSH_strain=factor(BSH_strain,levels=c("EcAZ-1-cat","AZ-52","Dny-BSH1","Dny-BSH2","Lgasseri-BSH",
                                               "LCAG-95-BSH1208","Ep-BSH101","25MeOH","50MeOH","ctrl")))%>%
    mutate(condition=factor(condition,levels=c("NA","FA","FT")))

p <- ggplot(BA_candi, aes(x=log10(baseline_chng+0.000001), colour=BA_suppl)) + 
  geom_density()+ facet_wrap(~BSH_strain, nrow=2)+
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")+
  theme_pubr()+labs(title="candidate BA")
ggsave("BSH/BSH_ENB/GNPS/change_baseline/candBA_density.pdf", p,width = 8, height = 4)

show_deconj<-function(mtb, BA,mtb_id){
  BA_candi_x<-mtb%>%filter(BA_suppl==BA)%>%
    mutate(cmpd_name_new=paste(FeatureID,Compound_Name,sep="_"))%>%
    filter(FeatureID==mtb_id)
    #filter(`0h`>0)
  p <- ggplot(BA_candi_x, aes(x=fct_rev(BSH_strain),y=log10(baseline_chng+0.000001),colour=condition)) +
    geom_boxplot()+
    #geom_point()+
    facet_wrap(~BA_suppl)+ coord_flip()+
    scale_colour_manual(values=c("#0072B2","#D55E00","#009E73"))+
    geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
    theme_pubr()+labs(title=unique(BA_candi_x$cmpd_name_new))+theme(legend.position = "none")
  
  return(p)
}
BA_candi_GCA<-show_deconj(BA_candi,"GCA","238_i")
ggsave("BSH/BSH_ENB/GNPS/change_baseline/GCA_238_i_deconj.pdf", BA_candi_GCA,width =4, height = 3)
BA_candi_GCA<-show_deconj(BA_candi,"GCA","185_i")
ggsave("BSH/BSH_ENB/GNPS/change_baseline/GCA_185_i_deconj.pdf", BA_candi_GCA,width =4, height = 3)
BA_candi_GCA<-show_deconj(BA_candi,"GCA","43992")
ggsave("BSH/BSH_ENB/GNPS/change_baseline/GCA_43992_deconj.pdf", BA_candi_GCA,width =4, height = 3)

BA_candi_TLCA<-show_deconj(BA_candi,"TLCA","240_i")
ggsave("BSH/BSH_ENB/GNPS/change_baseline/TLCA_240_i_deconj.pdf", BA_candi_TLCA,width =4, height = 3)
BA_candi_TDCA<-show_deconj(BA_candi,"TDCA","88_i")
ggsave("BSH/BSH_ENB/GNPS/change_baseline/TDCA_88_i_deconj.pdf", BA_candi_TDCA,width =4, height = 3)
BA_candi_TCA<-show_deconj(BA_candi,"TCA","201_i")
ggsave("BSH/BSH_ENB/GNPS/change_baseline/TCA_201_i_deconj.pdf", BA_candi_TCA,width =4, height = 3)
# BA_candi_TCA<-show_deconj(BA_candi,"TCA","50529")




BA_amine<-mtb%>%filter(FeatureID %in% BA_amine$FeatureID)%>%
  mutate(BSH_strain=factor(BSH_strain,levels=c("EcAZ-1-cat","AZ-52","Dny-BSH1","Dny-BSH2","Lgasseri-BSH",
                                        "LCAG-95-BSH1208","Ep-BSH101","25MeOH","50MeOH","ctrl")),
         BA_suppl=factor(BA_suppl,levels=c("GCA","TCA","GCDCA","TCDCA","GDCA","TDCA","GLCA","TLCA","GUDCA","TUDCA")))%>%
  group_by(FeatureID,name_new,BSH_strain,BA_suppl)%>%
  summarise(mn_chg_baseline=mean(baseline_chng))%>%
  left_join(annot,by="FeatureID")%>%
  mutate(cmpd_name_new=paste(FeatureID,Compound_Name))

BA_amine <- BA_amine %>%
  mutate(cmpd_name_new = factor(cmpd_name_new, levels = unique(BA_amine$cmpd_name_new[order(BA_amine$Compound_Name)])))

p <- ggplot(BA_amine, aes(x=fct_rev(BSH_strain), y=log10(mn_chg_baseline), fill=cmpd_name_new)) + 
  geom_bar(stat="identity")+
  facet_wrap(~BA_suppl, nrow=5)+
  scale_fill_manual(values = c("#B4A7D6","#582C83","#A6CEE3","#3789AD","#87CA6A","#3DA533","#238B45","#979C62",
                               "#F98F8E","#E94330","#FDA440","#FF7F00","#783F04")) +
  scale_y_continuous(expand=c(0,0), breaks = seq(-10, 20,by = 5))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
  theme_bw()+labs(title="amine BAs")+coord_flip()+
  theme(legend.position = "right")

ggsave(paste("BSH/BSH_ENB/GNPS/change_baseline/amine_conj_summ_tall.pdf",sep=""), p,width = 6, height = 7)

###########################################################

#load the counts table from Ipsita to create heatmap

annot<-fread("BSH/BSH_ENB/BSH_Zarrinpar_libhit.csv")%>%
  dplyr::rename(FeatureID=Scan)%>%
  mutate(FeatureID=as.character(FeatureID))%>%
  dplyr::select(FeatureID,Compound_Name)

order<-c("TCA_EcAZ-1-cat_0h","TCA_EcAZ-1-cat_48h",
         "TCA_AZ-52_0h","TCA_AZ-52_48h",
         "TCA_Dny_BSH1_0h","TCA_Dny_BSH1_48h",
         "TCA_Dny_BSH2_0h","TCA_Dny_BSH2_48h",
         "TCA_Lgasseri_0h","TCA_Lgasseri_48h", 
         "TCA_Ep_BSH101_0h","TCA_Ep_BSH101_48h",
         "TCA_LCAG-95_0h","TCA_LCAG-95_48h",
         "TDCA_EcAZ-1-cat_0h","TDCA_EcAZ-1-cat_48h",
         "TDCA_AZ-52_0h", "TDCA_AZ-52_48h",
         "TDCA_Dny_BSH1_0h","TDCA_Dny_BSH1_48h",
         "TDCA_Dny_BSH2_0h","TDCA_Dny_BSH2_48h",
         "TDCA_Lgasseri_0h","TDCA_Lgasseri_48h",
         "TDCA_Ep_BSH101_0h","TDCA_Ep_BSH101_48h",
         "TDCA_LCAG-95_0h","TDCA_LCAG-95_48h",
         "Dny_BSH1_blank_0h","Dny_BSH1_blank_48h",
         "Dny_BSH2_blank_0h","Dny_BSH2_blank_48h",
         "Lgasseri_blank_0h", "Lgasseri_blank_48h",
         "Ep_BSH101_blank_0h","Ep_BSH101_blank_48h",
         "LCAG-95_blank_0h","LCAG-95_blank_48h",
         "Media_TCA_blank_nobacteria_0h","Media_TCA_blank_nobacteria_48h",
         "Media_TDCA_blank_nobacteria_0h","Media_TDCA_blank_nobacteria_48h",
         "Media_blank_nobacteria_0h","Media_blank_nobacteria_48h",
         "EcAZ-1-cat_blank_0h","EcAZ-1-cat_blank_48h",
         "AZ-52_blank_0h","AZ-52_blank_48h",
         "Blank")

mtb<-read_csv("BSH/BSH_ENB/BSH_feature_with_metadat_sub.csv")%>%
  group_by(ATTRIBUTE_name) %>%
  summarize_all(mean)

mtb <- mtb[order(match(mtb$ATTRIBUTE_name, order)), ]

mtb<-mtb%>%column_to_rownames("ATTRIBUTE_name")%>% as.matrix()%>%
  log10()

mtb[ mtb== -Inf] <- 0

# mtb<-read_csv("BSH/BSH_ENB/BSH_feature_with_metadat_sub.csv")%>%
#   gather(FeatureID,peak_area,-ATTRIBUTE_name)%>%
#   left_join(.,annot, by="FeatureID")

p<-pheatmap(mtb,cluster_rows=FALSE,row_order = order)
ggsave("BSH/BSH_ENB/SFR24_0307_heatmapTCATDCA.pdf", p,width = 10, height = 8)


mtb_annot<-mtb%>%as.data.frame()%>%
  rownames_to_column("ATTRIBUTE_name")%>%
  gather(FeatureID,log_peak,-ATTRIBUTE_name)%>%
  left_join(.,annot,by="FeatureID")%>%
  mutate(name=paste(FeatureID,Compound_Name, sep=" "))%>%
  dplyr::select(name,ATTRIBUTE_name,log_peak)%>%
  spread(key=name, value=log_peak)
  
mtb_annot <- mtb_annot[order(match(mtb_annot$ATTRIBUTE_name, order)), ]
mtb_annot<-mtb_annot%>%
  rownames_to_column()%>%
  dplyr::select(-rowname)%>%
  column_to_rownames("ATTRIBUTE_name")%>% as.matrix()

p<-pheatmap(mtb_annot,cluster_rows=FALSE,row_order = order)
ggsave("BSH/BSH_ENB/SFR24_0307_heatmapTCATDCA_annot.pdf", p,width = 10, height = 30)

#subset heatmap 

order_new<-c("TCA_EcAZ-1-cat_0h",
         "TCA_AZ-52_0h",
         "TCA_Dny_BSH1_0h",
         "TCA_Dny_BSH2_0h",
         "TCA_Lgasseri_0h",
         "TCA_Ep_BSH101_0h",
         "TCA_LCAG-95_0h",
         "TDCA_EcAZ-1-cat_0h",
         "TDCA_AZ-52_0h", 
         "TDCA_Dny_BSH1_0h",
         "TDCA_Dny_BSH2_0h",
         "TDCA_Lgasseri_0h",
         "TDCA_Ep_BSH101_0h",
         "TDCA_LCAG-95_0h",
         "TCA_EcAZ-1-cat_48h",
         "TCA_AZ-52_48h",
         "TCA_Dny_BSH1_48h",
         "TCA_Dny_BSH2_48h",
         "TCA_Lgasseri_48h",
         "TCA_Ep_BSH101_48h",
         "TCA_LCAG-95_48h",
         "TDCA_EcAZ-1-cat_48h",
         "TDCA_AZ-52_48h",
         "TDCA_Dny_BSH1_48h",
         "TDCA_Dny_BSH2_48h",
         "TDCA_Lgasseri_48h",
         "TDCA_Ep_BSH101_48h",
         "TDCA_LCAG-95_48h")

mtb_noblanks<-mtb[1:28,]

mtb_noblanks<-mtb_noblanks%>%as.data.frame()%>%
  rownames_to_column("ATTRIBUTE_name")

mtb_noblanks<- mtb_noblanks[order(match(mtb_noblanks$ATTRIBUTE_name, order_new)), ]

mtb_noblanks<-mtb_noblanks%>%
  rownames_to_column()%>%
  dplyr::select(-rowname)%>%
  column_to_rownames("ATTRIBUTE_name")%>% as.matrix()

p<-pheatmap(mtb_noblanks,cluster_rows=FALSE)
ggsave("BSH/BSH_ENB/SFR24_0307_heatmapTCATDCA_noblanks_048h.pdf", p,width = 10, height = 4)

mtb_noblanks_annot<-mtb_noblanks%>%as.data.frame()%>%
  rownames_to_column("ATTRIBUTE_name")%>%
  gather(FeatureID,log_peak,-ATTRIBUTE_name)%>%
  left_join(.,annot,by="FeatureID")%>%
  mutate(name=paste(FeatureID,Compound_Name, sep=" "))%>%
  dplyr::select(name,ATTRIBUTE_name,log_peak)%>%
  spread(key=name, value=log_peak)

mtb_noblanks_annot <- mtb_noblanks_annot[order(match(mtb_noblanks_annot$ATTRIBUTE_name, order_new)), ]
mtb_noblanks_annot<-mtb_noblanks_annot%>%
  rownames_to_column()%>%
  dplyr::select(-rowname)%>%
  column_to_rownames("ATTRIBUTE_name")%>% as.matrix()

p<-pheatmap(mtb_noblanks_annot,cluster_rows=FALSE,row_order = order)
ggsave("BSH/BSH_ENB/SFR24_0307_heatmapTCATDCA_noblanks_048h_annot.pdf", p,width = 10, height = 25)

#TCA and TDCA

TCA<-c("12997","12865","9195","12880","12814","12897","12891","12901","12878","12861","12922","12902",
       "12905","12875","11900","9161","9273")
TDCA<-c("16285","16162","16167","16231","16330","16226","16235")

mtb_forplt<-read_csv("BSH/BSH_ENB/BSH_feature_with_metadat_sub.csv")%>%
  gather(FeatureID,peak_area,-ATTRIBUTE_name)%>%
  filter(FeatureID %in% TCA)%>%
  #filter(FeatureID %in% TDCA)%>%
  filter(!grepl("blank", ATTRIBUTE_name,ignore.case = TRUE))%>%
  separate(ATTRIBUTE_name,c("BA","BSH_origin","BSH_name","timepoint"),sep="_",remove=FALSE)%>%
  mutate(timepoint=ifelse(is.na(timepoint),BSH_name,timepoint),
         BSH_name=ifelse(BSH_name=="0h"|BSH_name=="48h",NA,BSH_name),
         name=sub("_[^_]+$", "", ATTRIBUTE_name),
         peak_area=log10(peak_area+1))%>%
  filter(BA=="TCA")%>%
  #filter(BA=="TDCA")%>%
  mutate(name=factor(name,levels=c("TCA_EcAZ-1-cat","TCA_AZ-52","TCA_Dny_BSH1","TCA_Dny_BSH2",
                                   "TCA_Lgasseri","TCA_Ep_BSH101","TCA_LCAG-95")))%>%
  # mutate(name=factor(name,levels=c("TDCA_EcAZ-1-cat","TDCA_AZ-52","TDCA_Dny_BSH1","TDCA_Dny_BSH2",
  #                                  "TDCA_Lgasseri","TDCA_Ep_BSH101","TDCA_LCAG-95")))%>%
  dplyr::select(-ATTRIBUTE_name)%>%
  arrange(FeatureID)

mtb_forplt<-mtb_forplt%>%spread(timepoint,peak_area)

stat.test<-mtb_forplt %>%
  group_by(name) %>%
  wilcox_test(data =., peak_area ~ timepoint) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")

stat.test

p<-ggpaired(mtb_forplt, cond1="0h",cond2="48h",fill="name", id="FeatureID")+
  facet_wrap(~name, nrow=1)+
  labs(title="TDCA", x="timepoint", y="log10(peak_area+1)")+
  theme(legend.position = "none")
ggsave("BSH/BSH_ENB/SFR24_0403_TCA048h_quant_paired.pdf", p,width = 12, height = 3)
ggsave("BSH/BSH_ENB/SFR24_0403_TDCA048h_quant_paired.pdf", p,width = 12, height = 3)

p<-ggplot(data=mtb_forplt, aes(x=timepoint, y=peak_area, colour=name, fill=name)) +
  geom_boxplot(color = "black",alpha=0.3)+
  geom_point(shape=21, color = "black", position=position_jitterdodge())+
  facet_wrap(~name, nrow=1)+
  theme_classic()+
  labs(title="TCA", x="timepoint", y="log10(peak_area+1)")+
  theme(legend.position = "none", plot.title = element_text(size = 12),axis.title.x = element_text(size = 10))
ggsave("BSH/BSH_ENB/SFR24_0403_TCA048h_quant_notpaired.pdf", p,width = 8, height = 2)
ggsave("BSH/BSH_ENB/SFR24_0403_TDCA048h_quant_notpaired.pdf", p,width = 8, height = 2)
