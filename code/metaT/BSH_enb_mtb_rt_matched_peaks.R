setwd("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_m_results/BSH/BSH_ENB/GNPS_rt_matching/rt_matched_peaks/")

library(tidyverse)
library(data.table)
library(ggpubfigs)
library(ggpubr)
###########################################################
#suppl BA

skyline<-fread("MSV000094578_Skyline_list.csv")
md<-fread("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_m_results/BSH/BSH_ENB/GNPS_missingplates/nf_output/metadata/ZT_BSH_metadata_addplates.txt")%>%
  mutate(filename = sub("\\..*", "", filename))

make_df_rt_plt<-function(sel_BA,mz,rt_min,rt_max){
  BA_files<-md%>%
    filter(BA_suppl==sel_BA & plate =="P1")
  
  BA_std<-fread(paste(sel_BA,"_std/",sel_BA,"_std_",mz,".csv",sep=""))%>%mutate(filename=paste(sel_BA,"_std_",mz,sep=""))
  
  df<-data.frame()
  for (i in 1:nrow(BA_files)) {
    ba_tmp<-fread(paste(BA_files$filename[i],"/",BA_files$filename[i],"_",mz,".csv",sep=""))
    ba_tmp<-ba_tmp%>%mutate(filename=BA_files$filename[i])
    df <- rbind(df, ba_tmp)
  }
  df<-df%>%
    rbind(BA_std)%>%
    filter(`Time(min)`> rt_min & `Time(min)` < rt_max)%>%
    left_join(.,md,by="filename")%>%
    mutate(BSH_strain=ifelse(is.na(BSH_strain),"std",BSH_strain))%>%
    mutate(condition=case_when(BSH_strain=="Dny-BSH1"|BSH_strain=="Dny-BSH2"~"FT",
                               BSH_strain=="LCAG-95-BSH1208"|BSH_strain=="Ep-BSH101"~"NA",
                               BSH_strain=="Lgasseri-BSH"~"FA",
                               .default = "none"))%>%
    mutate(BSH_strain=factor(BSH_strain,levels=c("std","EcAZ-1-cat","AZ-52","Dny-BSH1","Dny-BSH2","Lgasseri-BSH",
                                                 "LCAG-95-BSH1208","Ep-BSH101","ctrl")),
           condition=factor(condition,levels=c("NA","FA","FT","none")))
  return(df)
}

rt_plot<-function(df,rt_max,rt_seq,rt_title){
  
  p <- ggplot(df, aes(x=`Time(min)`, y=`Relative Abundance`, colour=condition)) + 
    geom_line()+ facet_grid(BSH_strain ~ .)+
    scale_y_continuous(breaks = seq(0,rt_max, by = rt_seq)) +
    scale_color_manual(values=c("#0072B2","#D55E00","#009E73","gray30"))+
    theme_pubr()+labs(title=rt_title)+theme(legend.position = "none")
  return(p)
}

GCA_df<-make_df_rt_plt("GCA","466_3163",5,6.5)
p<-rt_plot(GCA_df,2.0E+08,1E+08,"GCA (48h)")
ggsave("rt_matched_figs/GCA_rt_matched.pdf", p,width = 2.5, height = 6)

GCDCA_df<-make_df_rt_plt("GCDCA","450_3214",6.5,9.5)
p<-rt_plot(GCDCA_df,2.0E+08,0.5E+08,"GCDCA (48h)")
ggsave("rt_matched_figs/GCDCA_rt_matched.pdf", p,width = 2.5, height = 6)

GDCA_df<-make_df_rt_plt("GDCA","450_3214",6.5,9.5)
p<-rt_plot(GDCA_df,5.0E+08,2E+08,"GDCA (48h)")+
  scale_x_continuous(breaks = seq(7,9, by = 1))
ggsave("rt_matched_figs/GDCA_rt_matched.pdf", p,width = 2.5, height = 6)

GLCA_df<-make_df_rt_plt("GLCA","434_3265",8,9.5)
p<-rt_plot(GLCA_df,2.5E+08,1E+08,"GLCA (48h)")
ggsave("rt_matched_figs/GLCA_rt_matched.pdf", p,width = 2.5, height = 6)

GUDCA_df<-make_df_rt_plt("GUDCA","450_3214",5,6.5)
p<-rt_plot(GUDCA_df,2E+08,1E+08,"GUDCA (48h)")
ggsave("rt_matched_figs/GUDCA_rt_matched.pdf", p,width = 2.5, height = 6)

TCA_df<-make_df_rt_plt("TCA","516_299",4,5.5)
p<-rt_plot(TCA_df,1E+07,0.5E+07,"TCA (48h)")
ggsave("rt_matched_figs/TCA_rt_matched.pdf", p,width = 2.5, height = 6)

TCDCA_df<-make_df_rt_plt("TCDCA","500_304",5.5,7.5)
p<-rt_plot(TCDCA_df,2.0E+07,0.5E+07,"TCDCA (48h)")
ggsave("rt_matched_figs/TCDCA_rt_matched.pdf", p,width = 2.5, height = 6)

TDCA_df<-make_df_rt_plt("TDCA","500_304",6,7.5)
p<-rt_plot(TDCA_df,7.0E+07,3E+07,"TDCA (48h)")
ggsave("rt_matched_figs/TDCA_rt_matched.pdf", p,width = 2.5, height = 6)

TLCA_df<-make_df_rt_plt("TLCA","484_3091",8,9.5)
p<-rt_plot(TLCA_df,4E+07,2E+07,"TLCA (48h)")+
  scale_x_continuous(breaks = seq(8,10, by = 0.5))
ggsave("rt_matched_figs/TLCA_rt_matched.pdf", p,width = 2.5, height = 6)

TUDCA_df<-make_df_rt_plt("TUDCA","500_304",4,5.5)
p<-rt_plot(TUDCA_df,3E+07,1.5E+07,"TUDCA (48h)")
ggsave("rt_matched_figs/TUDCA_rt_matched.pdf", p,width = 2.5, height = 6)

###########################################################
#amine BA

skyline<-fread("MSV000094578_Skyline_list.csv")
md<-fread("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_m_results/BSH/BSH_ENB/GNPS_missingplates/nf_output/metadata/ZT_BSH_metadata_addplates.txt")%>%
  mutate(filename = sub("\\..*", "", filename))%>%
  mutate(name_v2 = str_replace(name, "^[^_]*_", "")) %>%
  mutate(name_v2 = str_replace(name_v2,"_[^_]*$", ""))%>%
  mutate(name_v2=paste(name_v2,BA_suppl,sep="_"))

make_df_rt_aa_plt<-function(sel_BA,aminoacid,mz,rt_min,rt_max){
  BA_files<-md%>%
    filter(BA_suppl==sel_BA & plate =="P1")
  
  #BA_std<-fread(paste(aminoacid,"_aminoacids/",aminoacid,"_aminoacids_",mz,".csv",sep=""))%>%mutate(filename=paste(aminoacid,"_aminoacids_",mz,sep=""))
  BA_std<-fread(paste(aminoacid,"_UDCA/",aminoacid,"_UDCA_std.csv",sep=""),skip = 3, header = TRUE)%>%mutate(name_v2=paste(aminoacid,"_UDCA",sep=""))
  
  df<-data.frame()
  for (i in 1:nrow(BA_files)) {
    ba_tmp<-fread(paste(aminoacid,"_UDCA/",aminoacid,"_UDCA_",BA_files$name_v2[i],".csv",sep=""),skip = 3, header = TRUE)
    #ba_tmp<-fread(paste(BA_files$filename[i],"/",BA_files$filename[i],"_",mz,".csv",sep=""))
    #ba_tmp<-ba_tmp%>%mutate(filename=BA_files$filename[i])
    ba_tmp<-ba_tmp%>%mutate(name_v2=BA_files$name_v2[i])
    df <- rbind(df, ba_tmp)
  }
  df<-df%>%
    rbind(BA_std)%>%
    #filter(`Time(min)`> rt_min & `Time(min)` < rt_max)%>%
    filter(`Time` > rt_min & `Time` < rt_max)%>%
    #left_join(.,md,by="filename")%>%
    left_join(.,BA_files,by="name_v2")%>%
    mutate(BSH_strain=ifelse(is.na(BSH_strain),"std",BSH_strain))%>%
    mutate(condition=case_when(BSH_strain=="Dny-BSH1"|BSH_strain=="Dny-BSH2"~"FT",
                               BSH_strain=="LCAG-95-BSH1208"|BSH_strain=="Ep-BSH101"~"NA",
                               BSH_strain=="Lgasseri-BSH"~"FA",
                               .default = "none"))%>%
    mutate(BSH_strain=factor(BSH_strain,levels=c("std","EcAZ-1-cat","AZ-52","Dny-BSH1","Dny-BSH2","Lgasseri-BSH",
                                                 "LCAG-95-BSH1208","Ep-BSH101","ctrl")),
           condition=factor(condition,levels=c("NA","FA","FT","none")))
  return(df)
}

rt_plot<-function(df,rt_max,rt_seq,rt_title){
  
  p <- ggplot(df, aes(x=Time, y=Intensity, colour=BA_suppl)) + 
    #ggplot(df, aes(x=`Time(min)`, y=`Relative Abundance`, colour=BA_suppl)) + 
    geom_line()+ facet_grid(BSH_strain ~ .)+
    scale_y_continuous(breaks = seq(0,rt_max,by = rt_seq),limits = c(0, rt_max)) +
    theme_pubr()+labs(title=rt_title)+theme(legend.position = "top")
  return(p)
}

ileleuUDCA_gudca_df<-make_df_rt_aa_plt("GUDCA","Ile","506_3840",8.5,9)
ileleuUDCA_tudca_df<-make_df_rt_aa_plt("TUDCA","Ile","506_3840",8.5,9)
ileleuUDCA_df<-rbind(ileleuUDCA_gudca_df,ileleuUDCA_tudca_df)

p<-rt_plot(ileleuUDCA_df,1E+08,0.5E+08,"Ile/Leu-UDCA (48h)")+ #1.5E6
     scale_x_continuous(breaks = seq(8,10, by = 0.2))
ggsave("rt_matched_figs/ileleu_UDCA_rt_matched_sclstd.pdf", p,width = 2.5, height = 6)

p<-rt_plot(ileleuUDCA_df,2E+06,1E+06,"Ile/Leu-UDCA (48h)")+
  scale_x_continuous(breaks = seq(8,10, by = 0.2))
ggsave("rt_matched_figs/ileleu_UDCA_rt_matched_sclothers.pdf", p,width = 2.5, height = 6)

# ileleuUDCA_gudca_df<-make_df_rt_aa_plt("GUDCA","UDCA","506_3840",8,9.5)
# ileleuUDCA_tudca_df<-make_df_rt_aa_plt("TUDCA","UDCA","506_3840",8,9.5)
# ileleuUDCA_df<-rbind(ileleuUDCA_gudca_df,ileleuUDCA_tudca_df)
# 
# p<-rt_plot(ileleuUDCA_df,2E+06,0.5E+06,"Ile/Leu-UDCA (48h)")+
#   scale_x_continuous(breaks = seq(8,10, by = 0.5))
# ggsave("rt_matched_figs/ileleu_UDCA_rt_matched.pdf", p,width = 2.5, height = 6)

alaUDCA_gudca_df<-make_df_rt_aa_plt("GUDCA","Ala","464_3371",6.75,7.4)
alaUDCA_tudca_df<-make_df_rt_aa_plt("TUDCA","Ala","464_3371",6.75,7.4)
alaUDCA_df<-rbind(alaUDCA_gudca_df,alaUDCA_tudca_df)

p<-rt_plot(alaUDCA_df,1.3E+06,0.5E+06,"Ala-UDCA (48h)")
ggsave("rt_matched_figs/ala_UDCA_rt_matched_sclDny.pdf", p,width = 2.5, height = 6)

p<-rt_plot(alaUDCA_df,2E+07,1E+07,"Ala-UDCA (48h)")
ggsave("rt_matched_figs/ala_UDCA_rt_matched_sclstd.pdf", p,width = 2.5, height = 6)

p<-rt_plot(alaUDCA_df,6E+04,3E+04,"Ala-UDCA (48h)")
ggsave("rt_matched_figs/ala_UDCA_rt_matched_sclother.pdf", p,width = 2.5, height = 6)

# alaUDCA_gudca_df<-make_df_rt_aa_plt("GUDCA","UDCA","464_3371",5.5,7)
# alaUDCA_tudca_df<-make_df_rt_aa_plt("TUDCA","UDCA","464_3371",5.5,7)
# alaUDCA_df<-rbind(alaUDCA_gudca_df,alaUDCA_tudca_df)
# 
# p<-rt_plot(alaUDCA_df,6E+05,3E+05,"Ala-UDCA (48h)")
# ggsave("rt_matched_figs/ala_UDCA_rt_matched.pdf", p,width = 2.5, height = 6)

# lysUDCA_gudca_df<-make_df_rt_aa_plt("GUDCA","UDCA","521_3949",3,4.5)
# lysUDCA_tudca_df<-make_df_rt_aa_plt("TUDCA","UDCA","521_3949",3,4.5)
# lysUDCA_df<-rbind(lysUDCA_gudca_df,lysUDCA_tudca_df)
# 
# p<-rt_plot(lysUDCA_df,8E+05,4E+05,"Lys-UDCA (48h)")
# ggsave("rt_matched_figs/lys_UDCA_rt_matched.pdf", p,width = 2.5, height = 6)

lysCA_gca_df<-make_df_rt_aa_plt("GCA","CA","537_3898",3.3,4.1)
lysCA_tca_df<-make_df_rt_aa_plt("TCA","CA","537_3898",3.3,4.1)
lysCA_df<-rbind(lysCA_gca_df,lysCA_tca_df)

p<-rt_plot(lysCA_df,4E+05,1.5E+05,"Lys-CA (48h)")
ggsave("rt_matched_figs/lys_CA_rt_matched_smrt.pdf", p,width = 2.5, height = 6)

lysCDCA_gcdca_df<-make_df_rt_aa_plt("GCDCA","CDCA","521_3949",5,5.7)
lysCDCA_tcdca_df<-make_df_rt_aa_plt("TCDCA","CDCA","521_3949",5,5.7)
lysCDCA_df<-rbind(lysCDCA_gcdca_df,lysCDCA_tcdca_df)

p<-rt_plot(lysCDCA_df,5E+05,2.5E+05,"Lys-CDCA (48h)")
ggsave("rt_matched_figs/lys_CDCA_rt_matched_smrt.pdf", p,width = 2.5, height = 6)




