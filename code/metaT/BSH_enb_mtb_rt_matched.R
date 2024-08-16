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
library(ComplexHeatmap)
library(RColorBrewer)
library(forcats)
library(ggbeeswarm)
library(ggsci)
###########################################################
md<-fread("BSH/BSH_ENB/GNPS/metadata_filename/ZT_BSH_metadata_justBA.tsv")%>%
  filter(!(plate=="P1"|plate=="P2"))
md_missing<-fread("BSH/BSH_ENB/GNPS_missingplates/nf_output/metadata/ZT_BSH_metadata_addplates.txt")
md<-rbind(md,md_missing)%>%
  mutate(cult_plt=case_when(plate=="P4"|plate=="P5"|plate=="P6"~"P2",
                            plate=="P7"|plate=="P8"|plate=="P9"~"P3",
                            .default = "P1"))
write.table(md,"BSH/BSH_ENB/GNPS_rt_matching/ZT_metadata_comb.txt",sep = "\t",row.names = FALSE, quote=FALSE)   


mtb_0h<-fread("BSH/BSH_ENB/GNPS_rt_matching/T_0hr_peak_areas_skyline_P4_P7.csv")
mtb_0h_mp<-fread("BSH/BSH_ENB/GNPS_rt_matching/T_0hr_peak_areas_skyline_P2_repeat.csv")
mtb_48h<-fread("BSH/BSH_ENB/GNPS_rt_matching/T_48hr_peak_areas_skyline_P6_P9.csv")
mtb_48h_mp<-fread("BSH/BSH_ENB/GNPS_rt_matching/T_48hr_peak_areas_skyline_P1_repeat.csv")

#rt matched quantification 
mtb<-rbind(mtb_0h,mtb_0h_mp,mtb_48h,mtb_48h_mp)%>%
  dplyr::rename(filename=`File Name`)%>%
  dplyr::select(-c(8:9))%>%
  left_join(.,md,by="filename")%>%
  filter(!(BA_suppl=="BHI-50DMSO"|BA_suppl=="BHI"))%>%
  mutate(log_10_peakabun=log10(Area+1))

###########################################################

TRF_unpr_ttest<-function(mtb){
  x<-pairwise.t.test(mtb$log_10_peakabun,mtb$timepoint,p.adjust.method = "none")
  df<-data.frame(t0hvt48h=x$p.value[1])
  return(df)
}

run_ttests_BA<-function(mtb,BA){
  stat.test<-mtb%>%
    group_by(BSH_strain,Molecule)%>%
    nest()%>%
    mutate(pvals = map(data, ~ TRF_unpr_ttest(.x)))%>%
    dplyr::select(-data)%>%
    unnest()
  write.table(stat.test,paste("BSH/BSH_ENB/GNPS_rt_matching/ttest_unpr_results/",BA,"suppl_ttest_unprd_tmpt_pvals_nofdr.txt",sep=""),sep = "\t",row.names = FALSE, quote=FALSE)  
}

BA_list<-c("GCA","GCDCA","GDCA","GLCA","GUDCA","TCA","TCDCA","TDCA","TLCA","TUDCA")
for(i in BA_list){
  mtb_new<-mtb%>%filter(BA_suppl==i)
  run_ttests_BA(mtb_new,i)
}
# mtb_sub<-mtb%>%filter(BA_suppl=="GUDCA"& Molecule=="Ala-UDCA" & BSH_strain=="Dny-BSH1")
# pairwise.t.test(mtb_sub$log_10_peakabun,mtb_sub$timepoint,p.adjust.method = "fdr")

mtb_bsln<-mtb%>%
  mutate(name_new=paste(BA_suppl,BSH_strain,sep="-"))%>%
  dplyr::select(Molecule,name_new,BA_suppl,BSH_strain,cult_plt,plate_location,timepoint,Area)%>%
  filter(timepoint!=""& !(BA_suppl%in% c("BHI","BHI-50DMSO")))%>%
  pivot_wider(names_from = timepoint, values_from = Area)%>%
  mutate(baseline_chng=`48h`/(`0h`+1))%>%
  #mutate(baseline_chng=`48h`-`0h`)%>%
  #filter(`0h`+`48h`>0)%>%
  mutate(condition=case_when(BSH_strain=="Dny-BSH1"|BSH_strain=="Dny-BSH2"~"FT",
                             BSH_strain=="LCAG-95-BSH1208"|BSH_strain=="Ep-BSH101"~"NA",
                             BSH_strain=="Lgasseri-BSH"~"FA",
                             .default = "none"))%>%
  mutate(baseline_chng=ifelse(`0h`+`48h`==0, NA,baseline_chng))%>%
  mutate(BSH_strain=factor(BSH_strain,levels=c("EcAZ-1-cat","AZ-52","Dny-BSH1","Dny-BSH2","Lgasseri-BSH",
                                               "LCAG-95-BSH1208","Ep-BSH101","25MeOH","50MeOH","ctrl")),
         condition=factor(condition,levels=c("NA","FA","FT")))%>%
  mutate(Molecule=ifelse(Molecule=="Ile/Leu-UDCA","Ile_Leu_UDCA",Molecule))

#these hits individually
show_deconj<-function(mtb,mtb_id){
  BA_candi_x<-mtb%>%filter(BA_suppl=="GCA")
  #BA_candi_x<-mtb%>%filter(BA_suppl=="GCDCA")
  #BA_candi_x<-mtb%>%filter(BA_suppl=="GDCA")
  #BA_candi_x<-mtb%>%filter(BA_suppl=="GLCA")
  #BA_candi_x<-mtb%>%filter(BA_suppl=="GUDCA")
  #BA_candi_x<-mtb%>%filter(BA_suppl=="TCA")
  #BA_candi_x<-mtb%>%filter(BA_suppl=="TCDCA")
  #BA_candi_x<-mtb%>%filter(BA_suppl=="TDCA")
  #BA_candi_x<-mtb%>%filter(BA_suppl=="TLCA")
  #BA_candi_x<-mtb%>%filter(BA_suppl=="TUDCA")
  p <- ggplot(BA_candi_x, aes(x=fct_rev(BSH_strain),y=log10(baseline_chng+0.000001),fill=condition,color=condition)) +
    #geom_boxplot(color = "black",alpha=0.3)+
    #geom_jitter()+
    geom_hline(yintercept = 0, linetype = "dashed", color = "red")+
    geom_beeswarm(shape=21,color="black",size=3) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.7, color = "black", linetype = "solid", size = 0.25)+
    facet_wrap(~BA_suppl)+ coord_flip()+
    scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
    scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+
    theme_pubr()+labs(title=mtb_id)+
    theme(legend.position = "none",panel.grid.major.y = element_line(color = "gray",size = 0.5,linetype = 2))
  
  ggsave(paste("BSH/BSH_ENB/GNPS_rt_matching/change_baseline/indv_mtb_plt/SFR24_0814_mtb",mtb_id,"GCA_v2.pdf",sep="_"), p,width = 4, height = 3)
  #ggsave(paste("BSH/BSH_ENB/GNPS_rt_matching/change_baseline/indv_mtb_plt/SFR24_0703_mtb",mtb_id,"GCDCA.pdf",sep="_"), p,width = 4, height = 3)
  #ggsave(paste("BSH/BSH_ENB/GNPS_rt_matching/change_baseline/indv_mtb_plt/SFR24_0703_mtb",mtb_id,"GDCA.pdf",sep="_"), p,width = 4, height = 3)
  #ggsave(paste("BSH/BSH_ENB/GNPS_rt_matching/change_baseline/indv_mtb_plt/SFR24_0703_mtb",mtb_id,"GLCA.pdf",sep="_"), p,width = 4, height = 3)
  #ggsave(paste("BSH/BSH_ENB/GNPS_rt_matching/change_baseline/indv_mtb_plt/SFR24_0703_mtb",mtb_id,"GUDCA.pdf",sep="_"), p,width = 4, height = 3)
  #ggsave(paste("BSH/BSH_ENB/GNPS_rt_matching/change_baseline/indv_mtb_plt/SFR24_0703_mtb",mtb_id,"TCA.pdf",sep="_"), p,width = 4, height = 3)
  #ggsave(paste("BSH/BSH_ENB/GNPS_rt_matching/change_baseline/indv_mtb_plt/SFR24_0703_mtb",mtb_id,"TCDCA.pdf",sep="_"), p,width = 4, height = 3)
  #ggsave(paste("BSH/BSH_ENB/GNPS_rt_matching/change_baseline/indv_mtb_plt/SFR24_0703_mtb",mtb_id,"TDCA.pdf",sep="_"), p,width = 4, height = 3)
  #ggsave(paste("BSH/BSH_ENB/GNPS_rt_matching/change_baseline/indv_mtb_plt/SFR24_0703_mtb",mtb_id,"TLCA.pdf",sep="_"), p,width = 4, height = 3)
  #ggsave(paste("BSH/BSH_ENB/GNPS_rt_matching/change_baseline/indv_mtb_plt/SFR24_0703_mtb",mtb_id,"TUDCA.pdf",sep="_"), p,width = 4, height = 3)
}


indv_plot <- mtb_bsln %>%
  group_by(Molecule) %>%
  nest() %>%
  mutate(pval = map2(data, Molecule,  ~ show_deconj(.x,.y)))

#get stats for change baseline
TRF_unpr_ttest_conj<-function(mtb){
  x<-pairwise.t.test(mtb$baseline_chng,mtb$BSH_strain,p.adjust.method = "none")
  df<-data.frame(AZ52 = x$p.value[1],
                 DnyBSH1 = x$p.value[2],
                 DnyBSH2 = x$p.value[3],
                 LgasseriBSH = x$p.value[4],
                 LCAG95BSH1208 = x$p.value[5],
                 EpBSH101 = x$p.value[6],
                 ctrl = x$p.value[7])
  return(df)
}

stat.test.uncong<-mtb_bsln%>%
  mutate(Molecule=ifelse(Molecule=="TCAA","TCA",Molecule))%>%
  filter(Molecule==BA_suppl)%>%
  group_by(Molecule)%>%
    nest()%>%
    mutate(pvals = map(data, ~ TRF_unpr_ttest_conj(.x)))%>%
    dplyr::select(-data)%>%
    unnest()

write.table(stat.test.uncong,"BSH/BSH_ENB/GNPS_rt_matching/ttest_unpr_results/BAsuppl_ttest_unprd_baselinechngbyaz-cat_fdr.txt",sep = "\t",row.names = FALSE, quote=FALSE)  
write.table(stat.test.uncong,"BSH/BSH_ENB/GNPS_rt_matching/ttest_unpr_results/BAsuppl_ttest_unprd_baselinechngbyaz-cat_nofdr.txt",sep = "\t",row.names = FALSE, quote=FALSE)  

mtb_bsln_tdca<-mtb_bsln%>%filter(Molecule=="TDCA"& BA_suppl=="TDCA")
pairwise.t.test(mtb_bsln_tdca$baseline_chng,mtb_bsln_tdca$BSH_strain,p.adjust.method = "none")

mtb_bsln_tlca<-mtb_bsln%>%filter(Molecule=="TLCA"& BA_suppl=="TLCA")
pairwise.t.test(mtb_bsln_tlca$baseline_chng,mtb_bsln_tlca$BSH_strain,p.adjust.method = "none")
###########################################################

p <- ggplot(mtb_bsln, aes(x=log10(baseline_chng+0.000001), colour=BA_suppl)) + 
  geom_density()+ facet_wrap(~BSH_strain, nrow=2)+
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")+
  theme_pubr()+labs(title="candidate BA")
ggsave("BSH/BSH_ENB/GNPS/change_baseline/candBA_density.pdf", p,width = 8, height = 4)

BA_amine<-mtb_bsln%>%
  mutate(Molecule=ifelse(Molecule=="Ile_Leu_UDCA","IleLeu-UDCA",Molecule))%>%
  filter(str_detect(Molecule, "[^- ]-[^- ]"))%>%
  mutate(BSH_strain=factor(BSH_strain,levels=c("EcAZ-1-cat","AZ-52","Dny-BSH1","Dny-BSH2","Lgasseri-BSH",
                                               "LCAG-95-BSH1208","Ep-BSH101","25MeOH","50MeOH","ctrl")),
         BA_suppl=factor(BA_suppl,levels=c("GCA","TCA","GCDCA","TCDCA","GDCA","TDCA","GLCA","TLCA","GUDCA","TUDCA")))%>%
  group_by(Molecule,BSH_strain,BA_suppl)%>%
  summarise(mn_chg_baseline=mean(baseline_chng))%>%
  filter(!(Molecule=="Phe-UDCA"|Molecule=="Lys-UDCA"))

p <- ggplot(BA_amine, aes(x=fct_rev(BSH_strain), y=log10(mn_chg_baseline), fill=Molecule)) + 
  geom_bar(stat="identity")+
  facet_wrap(~BA_suppl, nrow=5)+
  # scale_fill_manual(values=c("#2A6EBB","#F0AB00","#C50084","#F98F8E","#7D5CC6",
  #                            "#E37222","#69BE28","#238B45"))+
  scale_fill_manual(values=c("#2A6EBB","#69BE28","#7D5CC6",
                             "#C50084","#FDA440","#E37222"))+
  scale_y_continuous(expand=c(0,0), limits = c(0, 25), breaks = seq(-10, 30,by = 5))+
  #geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
  theme_bw()+labs(title="amine BAs")+coord_flip()+
  theme(legend.position = "right")

ggsave("BSH/BSH_ENB/GNPS_rt_matching/change_baseline/amine_conj_summ_tall.pdf", p,width = 6, height = 7)
ggsave("BSH/BSH_ENB/GNPS_rt_matching/change_baseline/amine_conj_summ_tall_rmpheudca.pdf", p,width = 5.5, height = 7)
ggsave("BSH/BSH_ENB/GNPS_rt_matching/change_baseline/amine_conj_summ_tall_rmpheudcalysudca.pdf", p,width = 5.5, height = 7)
ggsave("BSH/BSH_ENB/GNPS_rt_matching/change_baseline/amine_conj_summ_tall_rmpheudcalysudca_48h-0h.pdf", p,width = 5.5, height = 7)

###########################################################

#plot three examples of conjugated bile acids

plt_sig_mtbs<-function(mtb,mtb_id){
  p<-ggplot(data=mtb, aes(x=timepoint, y=log_10_peakabun, colour=BA_suppl, fill=BA_suppl)) +
    geom_boxplot(color = "black",alpha=0.3)+
    geom_point(shape=21, color = "black", position=position_jitterdodge())+
    facet_wrap(~BSH_strain, nrow=1)+
    theme_classic()+
    labs(title=paste(mtb_id," (",mtb_id,")",sep=""),
         x="timepoint", y="log10(peak_area+1)")+
    theme(legend.position = "right", plot.title = element_text(size = 12),axis.title.x = element_text(size = 10))
  #ggsave(paste("BSH/BSH_ENB/GNPS_rt_matching/ttest_unpr_results/indv_mtb_plt/SFR24_0501_mtb",mtb_id,"_notpaired.pdf",sep=""), p,width = 8, height = 2)
  #ggsave(paste("BSH/BSH_ENB/GNPS_rt_matching/ttest_unpr_results/indv_mtb_plt_just_gudcatudca/SFR24_0501_mtb",mtb_id,"_notpaired.pdf",sep=""), p,width = 8, height = 2)
  #ggsave(paste("BSH/BSH_ENB/GNPS_rt_matching/ttest_unpr_results/indv_mtb_plt_just_gcdcatcdca/SFR24_0501_mtb",mtb_id,"_notpaired.pdf",sep=""), p,width = 8, height = 2)
  ggsave(paste("BSH/BSH_ENB/GNPS_rt_matching/ttest_unpr_results/indv_mtb_plt_just_gcdatca/SFR24_0501_mtb",mtb_id,"_notpaired.pdf",sep=""), p,width = 8, height = 2)
  
}

mtb_sub<-mtb%>%
  mutate(Molecule=ifelse(Molecule=="Ile/Leu-UDCA","Ile_Leu_UDCA",Molecule))%>%
  group_by(BA_suppl,Molecule)%>%
  summarise(BA_summ=sum(log_10_peakabun))%>%
  filter(BA_summ>0)%>%
  mutate(suppl_mol=paste(BA_suppl,Molecule,sep="_"))

BA_mtb<-mtb%>%
  mutate(Molecule=ifelse(Molecule=="Ile/Leu-UDCA","Ile_Leu_UDCA",Molecule))%>%
  mutate(suppl_mol=paste(BA_suppl,Molecule,sep="_"))%>%
  filter(suppl_mol %in% mtb_sub$suppl_mol)%>%
  #filter(BA_suppl=="GUDCA"|BA_suppl=="TUDCA")%>%
  #filter(BA_suppl=="GCDCA"|BA_suppl=="TCDCA")%>%
  filter(BA_suppl=="GCA"|BA_suppl=="TCA")%>%
  mutate(BSH_strain=factor(BSH_strain,levels=c("EcAZ-1-cat","AZ-52","Dny-BSH1","Dny-BSH2","Lgasseri-BSH",
                                               "LCAG-95-BSH1208","Ep-BSH101","25MeOH","50MeOH","ctrl")),
         BA_suppl=factor(BA_suppl,levels=c("GCA","TCA","GCDCA","TCDCA","GDCA","TDCA","GLCA","TLCA","GUDCA","TUDCA")))

indv_plot<-BA_mtb%>%
  group_by(Molecule)%>%
  nest()%>%
  mutate(pval = map2(data, Molecule,  ~ plt_sig_mtbs(.x,.y)))
