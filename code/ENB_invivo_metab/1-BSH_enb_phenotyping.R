setwd("~/Notebooks/sfloresr/TRF-metaT/")

library(tidyverse)
library(data.table)
library(ggpubr)

###########################################################################
#paths
load_dat<-"data/ENB_invivo_metab/BSH_enb_micro_load_wk6.txt"
food_dat<-"data/ENB_invivo_metab/BSH_enb_food_intake_cumul_g.txt"
weight_dat<-"data/ENB_invivo_metab/BSH_bodyweight_g.txt"
fglucose_dat<-"data/ENB_invivo_metab/BSH_enb_fasting_blood_glucose.txt"
pglucose_dat<-"data/ENB_invivo_metab/BSH_enb_postprandial_blood_glucose.txt"
insulin_dat<-"data/ENB_invivo_metab/BSH_enb_postprandial_serum_insulin.txt"
fat_dat<-"data/ENB_invivo_metab/BSH_enb_fatmass_percbodywgt.txt"
lean_dat<-"data/ENB_invivo_metab/BSH_enb_leanmass_percbodywgt.txt"
fig_path<-"figures/ENB_invivo_metab/"

###########################################################################
#microbial load--Fig 6B

df_load<-fread(load_dat)%>%
  mutate(BSH_enb=factor(BSH_enb,levels=c("EcAZ-2","EcAZ-2BSH+","EcAZ_DneBSH1","EcAZ_LgaBSH")))

load <- ggplot(df_load, aes(x=BSH_enb, y=log10_CFU_g_feces, group=BSH_enb, label=BSH_enb,fill=BSH_enb))+
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(0.75)) +
  theme_pubr()+
  scale_fill_manual(values = c("gray70","#CC79A7","#009E73","#D55E00")) +
  labs(title="Microbial load, Wk 6") +
  scale_y_continuous(limits = c(0,10))+
  ylab("log10(CFU/g feces) ")+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 

ggsave(paste0(fig_path,"BSH_enb_micro_load_wk6.pdf"), plot=load,height=3.5, width=3)

pairwise.t.test(df_load$log10_CFU_g_feces, df_load$BSH_enb,
                     p.adjust.method = "none")

# EcAZ-2 EcAZ-2BSH+ EcAZ_DneBSH1
# EcAZ-2BSH+   0.124  -          -           
#   EcAZ_DneBSH1 0.365  0.515      -           
#   EcAZ_LgaBSH  0.097  0.897      0.436       

#none are significantly different

##########################################################################
#cumulative food intake--Figure 6D

df_foodcons_cum<-fread(food_dat)%>%
  separate_rows(mouseID,sep=", ")%>%
  gather(key="week",value="foodcons_g",-mouseID,-cageID,-BSH_enb)%>%
  mutate(week = as.numeric(str_extract(week, "\\d+")),
         foodcons_kcal=foodcons_g*3.1)

stderror <- function(x) sd(x)/sqrt(length(x))

df_foodcons_cum_summary <- df_foodcons_cum %>%
  group_by(week, BSH_enb) %>%
  summarise(mn_foodcons_g = mean(foodcons_g, na.rm = TRUE),
            SEM_foodcons_g = stderror(foodcons_g),
            mn_foodcons_kcal = mean(foodcons_kcal, na.rm = TRUE),
            SEM_foodcons_kcal = stderror(foodcons_kcal),)%>%
  mutate(BSH_enb = factor(BSH_enb, levels = c("EcAZ-2", "EcAZ-2BSH+", "EcAZ_DneBSH1", "EcAZ_LgaBSH")))

foodcons_cum_kcal_plot <- ggplot(df_foodcons_cum_summary, aes(x= week, y= mn_foodcons_kcal))+
  geom_line(aes(color=BSH_enb), size=0.4)+
  geom_errorbar(aes(ymin=mn_foodcons_kcal-SEM_foodcons_kcal, ymax=mn_foodcons_kcal+SEM_foodcons_kcal),
                width=.2, padding=0.2)+
  geom_point(aes(fill=BSH_enb), colour="black",pch=21, size=2) +
  labs(title="Cumulative Food\nIntake vs. Time") +
  ylab("Food Intake (kcal/mouse)")+
  xlab("Week (after single gavage)")+
  theme_classic()+
  scale_x_continuous(breaks = seq(0, 6, by = 1))+
  scale_fill_manual(values = c("gray70","#CC79A7","#009E73","#D55E00")) +
  scale_color_manual(values = c("gray70","#CC79A7","#009E73","#D55E00")) +
  theme(plot.title = element_text(hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "top") 

ggsave(paste0(fig_path,"BSH_enb_food_intake_cumul_kcal.pdf"), plot=foodcons_cum_kcal_plot,height=3.2, width=3)

for (i in 1:6){
  bshenb_foodcons_t_test <- pairwise.t.test(x= subset(df_foodcons_cum, week==i)$foodcons_kcal, 
                                            g=subset(df_foodcons_cum, week==i)$BSH_enb,
                                            p.adjust.method = "fdr")
  print(i)
  print(bshenb_foodcons_t_test)
}#none are significantly different

############################################################################
#body weight g--Figure 6E

df_weight<-fread(weight_dat)%>%
  gather(key="week",value="weight_g",-mouseID,-cageID,-BSH_enb)%>%
  mutate(week = as.numeric(str_extract(week, "\\d+")))

stderror <- function(x) sd(x)/sqrt(length(x))

df_weight_summary <- df_weight %>%
  group_by(week, BSH_enb) %>%
  summarise(mn_weight_g = mean(weight_g, na.rm = TRUE),
            SEM = stderror(weight_g))%>%
  mutate(BSH_enb = factor(BSH_enb, levels = c("EcAZ-2", "EcAZ-2BSH+", "EcAZ_DneBSH1", "EcAZ_LgaBSH")))

weights_plot <- ggplot(df_weight_summary, aes(x= week, y= mn_weight_g))+
  geom_line(aes(color=BSH_enb), size=0.4)+
  geom_errorbar(aes(ymin=mn_weight_g-SEM, ymax=mn_weight_g+SEM),
                width=.2, padding=0.2)+
  geom_point(aes(fill=BSH_enb), colour="black",pch=21, size=2) +
  labs(title="Average Mouse\nWeight vs. Time") +
  ylab("Body weight (g)")+
  xlab("Week (after single gavage)")+
  theme_classic()+
  scale_x_continuous(breaks = seq(0, 6, by = 1))+
  scale_fill_manual(values = c("gray70","#CC79A7","#009E73","#D55E00")) +
  scale_color_manual(values = c("gray70","#CC79A7","#009E73","#D55E00")) +
  theme(plot.title = element_text(hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "top") 

ggsave(paste0(fig_path,"BSH_bodyweight_g.pdf"), plot=weights_plot,height=3.2, width=3)


for (i in 0:6){
  bshenb_weight_t_test <- pairwise.t.test(x= subset(df_weight, week==i)$weight_g, 
                                             g=subset(df_weight, week==i)$BSH_enb,
                                             p.adjust.method = "fdr")
  print(i)
  print(bshenb_weight_t_test)
} #none are significantly different 

###########################################################################
#fasting glucose--Figure S6A

df_fastglucose<-fread(fglucose_dat)%>%
  gather(key="week",value="fastglucose_mgdl",-mouseID,-cageID,-BSH_enb)%>%
  mutate(week = str_extract(week, "^[^_]+"))%>%
  mutate(BSH_enb=factor(BSH_enb,levels=c("EcAZ-2","EcAZ-2BSH+","EcAZ_DneBSH1","EcAZ_LgaBSH")))

fastglucose <- ggplot(df_fastglucose, aes(x=BSH_enb, y=fastglucose_mgdl, group=BSH_enb, label=BSH_enb,fill=BSH_enb))+
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(0.75)) +
  facet_wrap(~week)+
  theme_pubr()+
  scale_fill_manual(values = c("gray70","#CC79A7","#009E73","#D55E00")) +
  labs(title="Fasting Blood Glucose") +
  ylab("Fasting blood glucose (mg/dl)")+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 

ggsave(paste0(fig_path,"BSH_enb_fasting_blood_glucose.pdf"), plot=fastglucose,height=4, width=4)

#just wk3
df_fastglucose_wk3<-df_fastglucose%>%
  filter(week=="wk3")

#just wk6
df_fastglucose_wk6<-df_fastglucose%>%
  filter(week=="wk6")

pairwise.t.test(df_fastglucose_wk3$fastglucose_mgdl, df_fastglucose_wk3$BSH_enb,
                p.adjust.method = "none")

# EcAZ-2 EcAZ-2BSH+ EcAZ_DneBSH1
# EcAZ-2BSH+   0.854  -          -           
#   EcAZ_DneBSH1 0.425  0.538      -           
#   EcAZ_LgaBSH  0.257  0.190      0.058 

pairwise.t.test(df_fastglucose_wk6$fastglucose_mgdl, df_fastglucose_wk6$BSH_enb,
                     p.adjust.method = "none")

# EcAZ-2 EcAZ-2BSH+ EcAZ_DneBSH1
# EcAZ-2BSH+   0.3156 -          -           
#   EcAZ_DneBSH1 0.0282 0.2129     -           
#   EcAZ_LgaBSH  0.0044 0.0509     0.4572 

###########################################################################
#postprandial glucose--Figure 6F

df_prandglucose<-fread(pglucose_dat)%>%
  mutate(BSH_enb=factor(BSH_enb,levels=c("EcAZ-2","EcAZ-2BSH+","EcAZ_DneBSH1","EcAZ_LgaBSH")))

prandglucose <- ggplot(df_prandglucose, aes(x=BSH_enb, y=wk6_mgdl, group=BSH_enb, label=BSH_enb,fill=BSH_enb))+
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(0.75)) +
  theme_pubr()+
  scale_fill_manual(values = c("gray70","#CC79A7","#009E73","#D55E00")) +
  labs(title="Postprandial Blood Glucose, Wk6") +
  ylab("Postprandial blood glucose (mg/dl)")+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
ggsave(paste0(fig_path,"BSH_enb_postprandial_blood_glucose.pdf"), plot=prandglucose,height=4, width=3)

pairwise.t.test(df_prandglucose$wk6_mgdl, df_prandglucose$BSH_enb,
                p.adjust.method = "none")

# EcAZ-2 EcAZ-2BSH+ EcAZ_DneBSH1
# EcAZ-2BSH+   0.011  -          -           
#   EcAZ_DneBSH1 0.041  0.576      -           
#   EcAZ_LgaBSH  0.240  0.146      0.363  

###########################################################################
#postprandial insulin--Figure 6G

df_prandinsulin<-fread(insulin_dat)%>%
  mutate(BSH_enb=factor(BSH_enb,levels=c("EcAZ-2","EcAZ-2BSH+","EcAZ_DneBSH1","EcAZ_LgaBSH")))

prandinsulin <- ggplot(df_prandinsulin, aes(x=BSH_enb, y=wk6_ngdl, group=BSH_enb, label=BSH_enb,fill=BSH_enb))+
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(0.75)) +
  theme_pubr()+
  scale_fill_manual(values = c("gray70","#CC79A7","#009E73","#D55E00")) +
  labs(title="Postprandial Serum Insulin, Wk6") +
  ylab("Postprandial serum insulin (ng/dl)")+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
ggsave(paste0(fig_path,"BSH_enb_postprandial_serum_insulin.pdf"), plot=prandinsulin,height=4, width=3)

pairwise.t.test(df_prandinsulin$wk6_ngdl, df_prandinsulin$BSH_enb,
                p.adjust.method = "none")

# EcAZ-2 EcAZ-2BSH+ EcAZ_DneBSH1
# EcAZ-2BSH+   0.219  -          -           
#   EcAZ_DneBSH1 0.025  0.280      -           
#   EcAZ_LgaBSH  0.054  0.461      0.727 

###########################################################################
#fatmass diff--Figure 6H

df_fatmass_diff<-fread(fat_dat)%>%
  mutate(fatmass_diff=wk6-wk0)%>%
  mutate(BSH_enb=factor(BSH_enb,levels=c("EcAZ-2","EcAZ-2BSH+","EcAZ_DneBSH1","EcAZ_LgaBSH")))

fatmass_diff <- ggplot(df_fatmass_diff, aes(x=BSH_enb, y=fatmass_diff, group=BSH_enb, label=BSH_enb,fill=BSH_enb))+
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(0.75)) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed")+
  theme_pubr()+
  scale_fill_manual(values = c("gray70","#CC79A7","#009E73","#D55E00")) +
  labs(title="Fat mass") +
  scale_y_continuous(limits=c(-7,7),breaks = seq(-6, 6, by = 3))+
  ylab("Body weight change (%)")+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 

ggsave(paste0(fig_path,"BSH_enb_fatmass_percbodywgt_diff.pdf"), plot=fatmass_diff,height=4, width=3)

pairwise.t.test(df_fatmass_diff$fatmass_diff,df_fatmass_diff$BSH_enb,
                p.adjust.method = "none")

# EcAZ-2 EcAZ-2BSH+ EcAZ_DneBSH1
# EcAZ-2BSH+   0.0555 -          -           
#   EcAZ_DneBSH1 0.0075 0.3974     -           
#   EcAZ_LgaBSH  0.0677 0.9241     0.3473      

###########################################################################
#leanmass diff--Figure 6I

df_leanmass_diff<-fread(lean_dat)%>%
  mutate(leanmass_diff=wk6-wk0)%>%
  mutate(BSH_enb=factor(BSH_enb,levels=c("EcAZ-2","EcAZ-2BSH+","EcAZ_DneBSH1","EcAZ_LgaBSH")))

leanmass_diff <- ggplot(df_leanmass_diff, aes(x=BSH_enb, y=leanmass_diff, group=BSH_enb, label=BSH_enb,fill=BSH_enb))+
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(0.75)) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed")+
  theme_pubr()+
  scale_fill_manual(values = c("gray70","#CC79A7","#009E73","#D55E00")) +
  labs(title="Lean mass") +
  scale_y_continuous(limits=c(-5,10),breaks = seq(-10, 10, by = 2))+
  ylab("Body weight change (%)")+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(paste0(fig_path,"BSH_enb_leanmass_percbodywgt_diff.pdf"), plot=leanmass_diff,height=4, width=3)

pairwise.t.test(df_leanmass_diff$leanmass_diff,df_leanmass_diff$BSH_enb,
                p.adjust.method = "none")

# EcAZ-2 EcAZ-2BSH+ EcAZ_DneBSH1
# EcAZ-2BSH+   0.107  -          -           
#   EcAZ_DneBSH1 0.033  0.579      -           
#   EcAZ_LgaBSH  0.084  0.902      0.665       