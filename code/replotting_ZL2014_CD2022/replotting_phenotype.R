setwd("~/scratch/TRF_multiomics/phenotypic_fromTI")

library(tidyverse)
library(data.table)
library(ggpubfigs)

#################################################################
#plot weight
wgt_df<-fread("mouse_n_food_weights_TRF2014.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  group_by(week, condition)

wgt_df_mean<- summarise(wgt_df,weight= mean(weight,na.rm = T)) #average weight per group per week
wgt_df_SEM<- summarise(wgt_df, SEM= sd(weight,na.rm = T)/sqrt(length(weight))) #SEM per group per week
wgt_df_summarized <- merge(wgt_df_mean, wgt_df_SEM) %>% # merge mean weight & SEM tables
  filter(condition!="NT")%>%
  mutate(condition=factor(condition, levels = c("NA","FA","FT")))

#TRF cohort weights 
weights_plot <- ggplot(wgt_df_summarized, aes(x= week, y= weight))+
  geom_line(aes(color=condition), size=0.4)+
  geom_errorbar(aes(ymin=weight-SEM, ymax=weight+SEM),
                width=.2, padding=0.2)+
  geom_point(aes(fill=condition), colour="black",pch=21, size=2) +
  labs(title="Average Mouse\nWeight vs Time") +
  ylab("Weight (g)")+
  xlab("Weeks post intervention")+
  theme_classic()+
  scale_x_continuous(breaks = seq(0, 8, by = 1))+
  scale_fill_manual(values = friendly_pal("ito_seven")) +
  scale_color_manual(values = friendly_pal("ito_seven")) +
  theme(plot.title = element_text(hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "top") 

ggsave("avg_mouse_weight.pdf", plot=weights_plot,height=3.2, width=3)

for (i in 0:8){
  weight_t_test <- pairwise.wilcox.test(x= subset(wgt_df, week==i)$weight, 
                                             g=subset(wgt_df, week==i)$condition,
                                             p.adjust.method = "fdr")
  print(i)
  print(weight_t_test)
}

#################################################################
#plot glucose

glucose_df<-fread("microbiome_gtt_TRF2014.txt")%>%
  mutate(condition=ifelse(is.na(condition),"NA",condition))%>%
  group_by(time, condition)

glucose_df_mean<- summarise(glucose_df,glucose= mean(glucose,na.rm = T)) #average weight per group per week
glucose_df_SEM<- summarise(glucose_df, SEM= sd(glucose,na.rm = T)/sqrt(length(glucose))) #SEM per group per week
glucose_df_summarized <- merge(glucose_df_mean, glucose_df_SEM) %>% # merge mean weight & SEM tables
  filter(condition!="NT")%>%
  mutate(condition=factor(condition, levels = c("NA","FA","FT")))

glucose_plot <- ggplot(glucose_df_summarized, aes(x= time, y= glucose))+
  geom_line(aes(color=condition), size=0.4)+
  geom_errorbar(aes(ymin=glucose-SEM, ymax=glucose+SEM),
                width=.2, padding=0.2)+
  geom_point(aes(fill=condition), colour="black",pch=21, size=2) +
  labs(title="Blood Glucose vs Time") +
  ylab("Blood glucose (mg/dL)")+
  xlab("time")+
  theme_classic()+
  scale_x_continuous(breaks = seq(0, 120, by = 30))+
  scale_fill_manual(values = friendly_pal("ito_seven")) +
  scale_color_manual(values = friendly_pal("ito_seven")) +
  theme(plot.title = element_text(hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "top") 

ggsave("mouse_glucose.pdf", plot=glucose_plot,height=3.2, width=3)

time_list<-unique(glucose_df$time)
for (i in 1:5){
  
  glucose_t_test <- pairwise.wilcox.test(x= subset(glucose_df, time==time_list[i])$glucose, 
                                        g=subset(glucose_df, time==time_list[i])$condition,
                                        p.adjust.method = "fdr")
  print(time_list[i])
  print(glucose_t_test)
}
