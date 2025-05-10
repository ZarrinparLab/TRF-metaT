setwd("~/Notebooks/sfloresr/TRF-metaT/")

library(tidyverse)
library(data.table)

##############################################################################
#paths
ani_data<-"data/ENB_invitro_metab/prot_seqn_analysis_foley/BSH_ani_subset.txt"
fig_path<-"figures/ENB_invitro_metab/"
##############################################################################
#plot ANI for BSH of interest--Figure S5F

ani<-fread(ani_data)%>%
  dplyr::select(-V1)

names(ani)<-c("Group1","G000364225_101","prot_1","prot_4","prot_2","KRM51566.1","KFOAOCCA_01485","LGAS_RS00260",
              "G000014425_51","tr|A0A1U7NKD7|A0A1U7NKD7_9FIRM","prot_5","prot_6","tr|A0A1U7NP31|A0A1U7NP31_9FIRM",
              "JIJPODMP_01796","AZ52_04447","SFE73109.1","G009917455_1208","WP_056959219.1","EFK28582.1","prot_3" )

ord_list<-c("G000364225_101","prot_1","prot_4","prot_2","KRM51566.1","KFOAOCCA_01485","LGAS_RS00260",
            "G000014425_51","tr|A0A1U7NKD7|A0A1U7NKD7_9FIRM","prot_5","prot_6","tr|A0A1U7NP31|A0A1U7NP31_9FIRM",
            "JIJPODMP_01796","AZ52_04447","SFE73109.1","G009917455_1208","WP_056959219.1","EFK28582.1","prot_3")

ani<-ani%>%
  gather(Group2,ani_values,-Group1)%>%
  mutate(Group1=factor(Group1,levels=ord_list),
         Group2=factor(Group2,levels=ord_list))


p<-ggplot(ani, aes(x = Group2, y = Group1, fill = ani_values)) +
  geom_tile() +
  geom_text(aes(label = signif(ani_values,digits=3))) +
  scale_fill_distiller(palette = "Reds",direction=1) +
  theme_minimal() +
  theme(panel.grid = element_blank(),legend.position = "top")+labs(title="Percent Indentity Matrix",x="",y="")
ggsave(paste0(fig_path,"SFR25_0207_anisub.pdf"),plot=p, width = 10, height = 10)