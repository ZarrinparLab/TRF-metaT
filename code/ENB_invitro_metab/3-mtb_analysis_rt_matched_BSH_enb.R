setwd("~/Notebooks/sfloresr/TRF-metaT/")

library(tidyverse)
library(data.table)
library(ggrepel)
library(ggpubr)
library(ggbeeswarm)
library(rstatix)
library(forcats)
library(RColorBrewer)
###########################################################
#paths
dat_metadata<-"data/ENB_invitro_metab/metadata/ZT_BSH_metadata_justBA.tsv"
dat_metadata_mp<-"data/ENB_invitro_metab/metadata/ZT_BSH_metadata_addplates.txt"
dat_feattab_0h<-"data/ENB_invitro_metab/GNPS_rt_matching/T_0hr_peak_areas_skyline_P4_P7.csv"
dat_feattab_0h_mp<-"data/ENB_invitro_metab/GNPS_rt_matching/T_0hr_peak_areas_skyline_P2_repeat.csv"
dat_feattab_48h<-"data/ENB_invitro_metab/GNPS_rt_matching/T_48hr_peak_areas_skyline_P6_P9.csv"
dat_feattab_48h_mp<-"data/ENB_invitro_metab/GNPS_rt_matching/T_48hr_peak_areas_skyline_P1_repeat.csv"
dat_path<-"data/ENB_invitro_metab/"
fig_path<-"figures/ENB_invitro_metab/"
###########################################################
#functions

show_deconj<-function(mtb,mtb_id){
  p <- ggplot(mtb, aes(x=fct_rev(BSH_strain),y=log10(baseline_chng+0.000001),fill=condition,color=condition)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red")+
    geom_beeswarm(shape=21,color="black",size=3) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.7, color = "black", linetype = "solid", size = 0.25)+
    facet_wrap(~BA_suppl)+ coord_flip()+
    scale_fill_manual(values=c("#0072B2","#D55E00","#009E73"))+
    scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+
    theme_pubr()+labs(title=mtb_id)+
    theme(legend.position = "none")
  ggsave(paste(fig_path,"SFR24_0814_mtb_",mtb_id,".pdf",sep=""), p,width = 4, height = 3)
}

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

plt_sig_mtbs<-function(mtb,BA,mtb_id){
  mtb_sub<-mtb%>%filter((BA_suppl %in% BA) & Molecule==mtb_id)
  p<-ggplot(data=mtb_sub, aes(x=timepoint, y=log_10_peakabun, colour=BA_suppl, fill=BA_suppl)) +
    geom_boxplot(color = "black",alpha=0.3)+
    geom_point(shape=21, color = "black", position=position_jitterdodge())+
    facet_wrap(~BSH_strain, nrow=1)+
    theme_classic()+
    labs(title=paste(mtb_id," (",mtb_id,")",sep=""),
         x="timepoint", y="log10(peak_area+1)")+
    theme(legend.position = "right", plot.title = element_text(size = 12),axis.title.x = element_text(size = 10))
  return(p)
}

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
  write.table(stat.test,paste(dat_path,"GNPS_rt_matching/",BA,"suppl_ttest_unprd_tmpt_pvals.txt",sep=""),sep = "\t",row.names = FALSE, quote=FALSE)  
}

###########################################################
#load files
mtb_0h<-fread(dat_feattab_0h)
mtb_0h_mp<-fread(dat_feattab_0h_mp)
mtb_48h<-fread(dat_feattab_48h)
mtb_48h_mp<-fread(dat_feattab_48h_mp)

md<-fread(dat_metadata)%>%filter(!(plate=="P1"|plate=="P2"))
md_missing<-fread(dat_metadata_mp)

###########################################################
#clean files

md<-rbind(md,md_missing)%>%
  mutate(cult_plt=case_when(plate=="P4"|plate=="P5"|plate=="P6"~"P2",
                            plate=="P7"|plate=="P8"|plate=="P9"~"P3",
                            .default = "P1"))
write.table(md,paste0(dat_path,"metadata/ZT_metadata_comb.txt"),sep = "\t",row.names = FALSE, quote=FALSE)   

#rt matched quantification 
mtb<-rbind(mtb_0h,mtb_0h_mp,mtb_48h,mtb_48h_mp)%>%
  dplyr::rename(filename=`File Name`)%>%
  dplyr::select(-c(8:9))%>%
  left_join(.,md,by="filename")%>%
  filter(!(BA_suppl=="BHI-50DMSO"|BA_suppl=="BHI"))%>%
  mutate(log_10_peakabun=log10(Area+1))

###########################################################
#quantifying deconjugation of Gly and Tau BAs--Figure 5C

#get change from baseline 48h/0h
mtb_bsln<-mtb%>%
  mutate(name_new=paste(BA_suppl,BSH_strain,sep="-"))%>%
  dplyr::select(Molecule,name_new,BA_suppl,BSH_strain,cult_plt,plate_location,timepoint,Area)%>%
  filter(timepoint!=""& !(BA_suppl%in% c("BHI","BHI-50DMSO")))%>%
  pivot_wider(names_from = timepoint, values_from = Area)%>%
  mutate(baseline_chng=`48h`/(`0h`+1))%>%
  mutate(condition=case_when(BSH_strain=="Dny-BSH1"|BSH_strain=="Dny-BSH2"~"FT",
                             BSH_strain=="LCAG-95-BSH1208"|BSH_strain=="Ep-BSH101"~"NA",
                             BSH_strain=="Lgasseri-BSH"~"FA",
                             .default = "none"))%>%
  mutate(baseline_chng=ifelse(`0h`+`48h`==0, NA,baseline_chng))%>%
  mutate(BSH_strain=factor(BSH_strain,levels=c("EcAZ-1-cat","AZ-52","Dny-BSH1","Dny-BSH2","Lgasseri-BSH",
                                               "LCAG-95-BSH1208","Ep-BSH101","25MeOH","50MeOH","ctrl")),
         condition=factor(condition,levels=c("NA","FA","FT")))%>%
  mutate(Molecule=ifelse(Molecule=="Ile/Leu-UDCA","Ile_Leu_UDCA",Molecule))%>%
  mutate(Molecule=ifelse(Molecule=="TCAA","TCA",Molecule))

BA_list<-c("GCA","GCDCA","GDCA","GLCA","GUDCA","TCA","TCDCA","TDCA","TLCA","TUDCA")
for(i in BA_list){
  mtb_new<-mtb_bsln%>%filter(BA_suppl==i & Molecule==i)
  show_deconj(mtb_new,i)
}

#get stats for change baseline
stat.test.uncong<-mtb_bsln%>%
  filter(Molecule==BA_suppl)%>%
  group_by(Molecule)%>%
  nest()%>%
  mutate(pvals = map(data, ~ TRF_unpr_ttest_conj(.x)))%>%
  dplyr::select(-data)%>%
  unnest()

write.table(stat.test.uncong,paste0(dat_path,"GNPS_rt_matching/ttest_unpr_results/BAsuppl_ttest_unprd_baselinechngbyaz-cat.txt"),
            sep = "\t",row.names = FALSE, quote=FALSE)  

###########################################################
#quantifying BBAAs--Figure 5E

#change from baseline 48h-0h
mtb_cult<-mtb%>%
  mutate(name_new=paste(BA_suppl,BSH_strain,sep="-"))%>%
  dplyr::select(Molecule,name_new,BA_suppl,BSH_strain,cult_plt,plate_location,timepoint,Area)%>%
  filter(timepoint!=""& !(BA_suppl%in% c("BHI","BHI-50DMSO")))%>%
  filter(BSH_strain=="ctrl")%>%
  group_by(Molecule,BA_suppl,timepoint)%>%
  summarise(ctrl=mean(Area))%>%
  mutate(ctrl=as.numeric(ctrl))%>%
  mutate(ctrl=ifelse(ctrl==0,1,ctrl))

mtb_bsln<-mtb%>%
  mutate(name_new=paste(BA_suppl,BSH_strain,sep="-"))%>%
  dplyr::select(Molecule,name_new,BA_suppl,BSH_strain,cult_plt,plate_location,timepoint,Area)%>%
  filter(timepoint!=""& !(BA_suppl%in% c("BHI","BHI-50DMSO")))%>%
  filter(BSH_strain!="ctrl")%>%
  left_join(.,mtb_cult,by=c("Molecule","BA_suppl","timepoint"))%>%
  mutate(norm_abun=Area/ctrl)%>%
  dplyr::select(-ctrl,-Area)%>%
  pivot_wider(names_from = timepoint, values_from = norm_abun)%>%
  mutate(baseline_chng=`48h`-`0h`)%>%
  mutate(condition=case_when(BSH_strain=="Dny-BSH1"|BSH_strain=="Dny-BSH2"~"FT",
                             BSH_strain=="LCAG-95-BSH1208"|BSH_strain=="Ep-BSH101"~"NA",
                             BSH_strain=="Lgasseri-BSH"~"FA",
                             .default = "none"))%>%
  mutate(BSH_strain=factor(BSH_strain,levels=c("EcAZ-1-cat","AZ-52","Dny-BSH1","Dny-BSH2","Lgasseri-BSH",
                                               "LCAG-95-BSH1208","Ep-BSH101","25MeOH","50MeOH","ctrl")),
         condition=factor(condition,levels=c("NA","FA","FT")))%>%
  mutate(Molecule=ifelse(Molecule=="Ile/Leu-UDCA","Ile_Leu_UDCA",Molecule))

BA_amine<-mtb_bsln%>%
  mutate(Molecule=ifelse(Molecule=="Ile_Leu_UDCA","IleLeu-UDCA",Molecule))%>%
  filter(str_detect(Molecule, "[^- ]-[^- ]"))%>%
  mutate(BSH_strain=factor(BSH_strain,levels=c("EcAZ-1-cat","AZ-52","Dny-BSH1","Dny-BSH2","Lgasseri-BSH",
                                               "LCAG-95-BSH1208","Ep-BSH101","25MeOH","50MeOH","ctrl")),
         BA_suppl=factor(BA_suppl,levels=c("GCA","TCA","GCDCA","TCDCA","GDCA","TDCA","GLCA","TLCA","GUDCA","TUDCA")))%>%
  group_by(Molecule,BSH_strain,BA_suppl)%>%
  summarise(mn_chg_baseline=mean(baseline_chng))%>%
  filter(!(Molecule=="Phe-UDCA"|Molecule=="Lys-UDCA"))

p <- ggplot(BA_amine, aes(x=fct_rev(BSH_strain), y=mn_chg_baseline, fill=Molecule)) + 
  geom_bar(stat="identity")+
  facet_wrap(~BA_suppl, nrow=5)+
  scale_fill_manual(values=c("#2A6EBB","#69BE28","#7D5CC6",
                             "#C50084","#FDA440","#E37222"))+
  theme_bw()+labs(title="amine BAs")+coord_flip()+
  theme(legend.position = "right")

ggsave(paste0(fig_path,"amine_conj_summ_tall_rmpheudcalysudca_48h-0h_nolognormctrl.pdf"), p,width = 5.5, height = 7)

###########################################################
#plot BBAAs produced and determine differences--Figure 5F-I

mtb_sub<-mtb%>%
  mutate(Molecule=ifelse(Molecule=="Ile/Leu-UDCA","Ile_Leu_UDCA",Molecule))%>%
  group_by(BA_suppl,Molecule)%>%
  summarise(BA_summ=sum(log_10_peakabun))%>%
  filter(BA_summ>0)%>%
  mutate(suppl_mol=paste(BA_suppl,Molecule,sep="_"))

BA_mtb<-mtb%>%
  mutate(name_new=paste(BA_suppl,BSH_strain,sep="-"))%>%
  dplyr::select(Molecule,name_new,BA_suppl,BSH_strain,cult_plt,plate_location,timepoint,Area)%>%
  filter(timepoint!=""& !(BA_suppl%in% c("BHI","BHI-50DMSO")))%>%
  filter(BSH_strain!="ctrl")%>%
  left_join(.,mtb_cult,by=c("Molecule","BA_suppl","timepoint"))%>%
  mutate(norm_abun=Area/ctrl)%>%
  dplyr::select(-ctrl,-Area)%>%
  mutate(log_10_peakabun=log10(norm_abun+1))%>%
  mutate(Molecule=ifelse(Molecule=="Ile/Leu-UDCA","Ile_Leu_UDCA",Molecule))%>%
  mutate(suppl_mol=paste(BA_suppl,Molecule,sep="_"))%>%
  filter(suppl_mol %in% mtb_sub$suppl_mol)%>%
  mutate(BSH_strain=factor(BSH_strain,levels=c("EcAZ-1-cat","AZ-52","Dny-BSH1","Dny-BSH2","Lgasseri-BSH",
                                               "LCAG-95-BSH1208","Ep-BSH101","25MeOH","50MeOH","ctrl")),
         BA_suppl=factor(BA_suppl,levels=c("GCA","TCA","GCDCA","TCDCA","GDCA","TDCA","GLCA","TLCA","GUDCA","TUDCA")))

plt_sig_mtbs(BA_mtb,c("GUDCA","TUDCA"),"Ile_Leu_UDCA") #Fig5F
ggsave(paste0(fig_path,"SFR24_0501_gudcatudca_mtb_IleLeuUDCA_notpaired_norm.pdf"), width = 8, height = 2)
plt_sig_mtbs(BA_mtb,c("GUDCA","TUDCA"),"Ala-UDCA") #Fig5G
ggsave(paste0(fig_path,"SFR24_0501_gudcatudca_mtb_AlaUDCA_notpaired_norm.pdf"), width = 8, height = 2)
plt_sig_mtbs(BA_mtb,c("GCA","TCA"),"Lys-CA") #Fig5H
ggsave(paste0(fig_path,"SFR24_0501_gcatca_mtb_LysCA_notpaired_norm.pdf"), width = 8, height = 2)
plt_sig_mtbs(BA_mtb,c("GCDCA","TCDCA"),"Lys-CDCA") #Fig5I
ggsave(paste0(fig_path,"SFR24_0501_gcdcatcdca_mtb_LysCDCA_notpaired_norm.pdf"), width = 8, height = 2)

#run t-test, Table S5
BA_list<-c("GCA","GCDCA","GDCA","GLCA","GUDCA","TCA","TCDCA","TDCA","TLCA","TUDCA")
for(i in BA_list){
  mtb_new<-mtb%>%filter(BA_suppl==i)
  run_ttests_BA(mtb_new,i)
}

rt_matched_BA_results <- lapply(BA_list, function(BA) {
  file_path <- paste0(dat_path,"GNPS_rt_matching/ttest_unpr_results/", BA, "suppl_ttest_unprd_tmpt_pvals.txt")
  read.table(file_path, header = TRUE, sep = "\t") %>%
    dplyr::rename(rt_matched_detected_BA = Molecule) %>%
    mutate(culture_supplement_BA = BA)
}) %>% bind_rows()%>%arrange(t0hvt48h)

write.table(rt_matched_BA_results, file = paste0(dat_path,"GNPS_rt_matching/ttest_unpr_results/BA_invitro_results_combined.txt"), 
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
