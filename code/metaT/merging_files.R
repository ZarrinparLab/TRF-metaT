setwd("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics")
library(tidyverse)
library(data.table)

gonames<-fread("metatranscript/woltka2_m_results/pfam/go_name.txt")
pfam2GO<-read.table("metatranscript/woltka2_m_results/pfam/pfam-to-go-process.map",header = FALSE, sep = "\t",
                    col.names = paste0("V",seq_len(4)), fill = TRUE)%>%
  gather(column,GO_Term,-V1)%>%
  dplyr::select(1,3)%>%
  dplyr::rename(FeatureID=V1)%>%
  left_join(.,gonames,by="GO_Term")%>%
  filter(!is.na(name))%>%
  group_by(FeatureID)%>%
  slice(1)%>%
  select(-GO_Term)%>%
  rename(GO_Term=name)

feat_annot<-fread("metatranscript/woltka2_m_results/pfam/pfam_annotationkey.csv")%>%
  left_join(.,pfam2GO,by="FeatureID")

#aldex no LD
dir1<-"metatranscript/woltka2_m_results/pfam/aldex/SFR23_0620_"
FAFT.effect<-fread(paste0(dir1,"FAFT_ald_effectwpval_wannot.txt"))%>%
  dplyr::rename(rab.win.A=rab.win.FA,rab.win.B=rab.win.FT)%>%
  mutate(comparison="FAvFT",datatype="metaT")
FANA.effect<-fread(paste0(dir2,"FANA_ald_effectwpval_wannot.txt"))%>%
dplyr::rename(rab.win.A=rab.win.FA,rab.win.B=rab.win.NA)%>%
  mutate(comparison="FAvNA",datatype="metaT")
FTNA.effect<-fread(paste0(dir1,"FTNA_ald_effectwpval_wannot.txt"))%>%
  dplyr::rename(rab.win.A=rab.win.FT,rab.win.B=rab.win.NA)%>%
  mutate(comparison="FTvNA",datatype="metaT")

dir2<-"metagenomic/woltka2_results/filtered_metaG/pfam_notnorm/aldex/SFR23_0606_"
FAFT.effect.g<-fread(paste0(dir2,"FAFT_ald_effectwpval_wannot.txt"))%>%
  dplyr::rename(rab.win.A=rab.win.FA,rab.win.B=rab.win.FT)%>%
  mutate(comparison="FAvFT",datatype="metaG")
FANA.effect.g<-fread(paste0(dir2,"FANA_ald_effectwpval_wannot.txt"))%>%
  dplyr::rename(rab.win.A=rab.win.FA,rab.win.B=rab.win.NA)%>%
  mutate(comparison="FAvNA",datatype="metaG")
FTNA.effect.g<-fread(paste0(dir1,"FTNA_ald_effectwpval_wannot.txt"))%>%
  dplyr::rename(rab.win.A=rab.win.FT,rab.win.B=rab.win.NA)%>%
  mutate(comparison="FTvNA",datatype="metaG")

comALDEx<-rbind(FANA.effect,FTNA.effect,FAFT.effect,FANA.effect.g,FTNA.effect.g,FAFT.effect.g)%>%
  left_join(.,feat_annot,by="FeatureID")%>%
  dplyr::rename(Name=Name.x)%>%
  dplyr::select(FeatureID,Name,GO_Term,comparison,datatype,diffexpr,everything(),-Name.y)%>%
  arrange(wi.eBH)

write.table(comALDEx, file = "multiomics/aldex2_results_combined.txt", #Table S1
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#aldex LD
dir1<-"metatranscript/woltka2_m_results/pfam/aldex/SFR23_0620_"
FAFTL.effect.annot<-fread(paste0(dir1,"FAFTL_ald_effectwpval_wannot.txt"))%>%
  dplyr::rename(rab.win.A=rab.win.FA,rab.win.B=rab.win.FT)%>%
  mutate(comparison="FAvFT",phase="Light",datatype="metaT")
FAFTD.effect.annot<-fread(paste0(dir1,"FAFTD_ald_effectwpval_wannot.txt"))%>%
  dplyr::rename(rab.win.A=rab.win.FA,rab.win.B=rab.win.FT)%>%
  mutate(comparison="FAvFT",phase="Dark",datatype="metaT")
FTNAL.effect.annot<-fread(paste0(dir1,"FTNAL_ald_effectwpval_wannot.txt"))%>%
  dplyr::rename(rab.win.A=rab.win.FT,rab.win.B=rab.win.NA)%>%
  mutate(comparison="FTvNA",phase="Light",datatype="metaT")
FTNAD.effect.annot<-fread(paste0(dir1,"FTNAD_ald_effectwpval_wannot.txt"))%>%
  dplyr::rename(rab.win.A=rab.win.FT,rab.win.B=rab.win.NA)%>%
  mutate(comparison="FTvNA",phase="Dark",datatype="metaT")
FANAL.effect.annot<-fread(paste0(dir1,"FANAL_ald_effectwpval_wannot.txt"))%>%
  dplyr::rename(rab.win.A=rab.win.FA,rab.win.B=rab.win.NA)%>%
  mutate(comparison="FAvNA",phase="Light",datatype="metaT")
FANAD.effect.annot<-fread(paste0(dir1,"FANAL_ald_effectwpval_wannot.txt"))%>%
  dplyr::rename(rab.win.A=rab.win.FA,rab.win.B=rab.win.NA)%>%
  mutate(comparison="FAvNA",phase="Dark",datatype="metaT")

dir2<-"metagenomic/woltka2_results/filtered_metaG/pfam_notnorm/aldex/SFR23_0606_"
FAFTL.effect.annot.g<-fread(paste0(dir2,"FAFTL_ald_effectwpval_wannot.txt"))%>%
  dplyr::rename(rab.win.A=rab.win.FA,rab.win.B=rab.win.FT)%>%
  mutate(comparison="FAvFT",phase="Light",datatype="metaG")
FAFTD.effect.annot.g<-fread(paste0(dir2,"FAFTD_ald_effectwpval_wannot.txt"))%>%
  dplyr::rename(rab.win.A=rab.win.FA,rab.win.B=rab.win.FT)%>%
  mutate(comparison="FAvFT",phase="Dark",datatype="metaG")
FTNAL.effect.annot.g<-fread(paste0(dir2,"FTNAL_ald_effectwpval_wannot.txt"))%>%
  dplyr::rename(rab.win.A=rab.win.FT,rab.win.B=rab.win.NA)%>%
  mutate(comparison="FTvNA",phase="Light",datatype="metaG")
FTNAD.effect.annot.g<-fread(paste0(dir2,"FTNAD_ald_effectwpval_wannot.txt"))%>%
  dplyr::rename(rab.win.A=rab.win.FT,rab.win.B=rab.win.NA)%>%
  mutate(comparison="FTvNA",phase="Dark",datatype="metaG")
FANAL.effect.annot.g<-fread(paste0(dir2,"FANAL_ald_effectwpval_wannot.txt"))%>%
  dplyr::rename(rab.win.A=rab.win.FA,rab.win.B=rab.win.NA)%>%
  mutate(comparison="FAvNA",phase="Light",datatype="metaG")
FANAD.effect.annot.g<-fread(paste0(dir2,"FANAL_ald_effectwpval_wannot.txt"))%>%
  dplyr::rename(rab.win.A=rab.win.FA,rab.win.B=rab.win.NA)%>%
  mutate(comparison="FAvNA",phase="Dark",datatype="metaG")

comALDExLD<-rbind(FANAL.effect.annot,FTNAL.effect.annot,FAFTL.effect.annot,FANAD.effect.annot,FTNAD.effect.annot,FAFTD.effect.annot,
                  FANAL.effect.annot.g,FTNAL.effect.annot.g,FAFTL.effect.annot.g,FANAD.effect.annot.g,FTNAD.effect.annot.g,FAFTD.effect.annot.g)%>%
  left_join(.,feat_annot,by="FeatureID")%>%
  dplyr::rename(Name=Name.x)%>%
  dplyr::select(FeatureID,Name,GO_Term,comparison,phase,datatype,diffexpr,everything(),-Name.y)%>%
  arrange(wi.eBH)

write.table(comALDExLD, file = "multiomics/aldex2LD_results_combined.txt", #Table S2
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)









