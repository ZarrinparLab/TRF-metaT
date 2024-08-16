setwd("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_m_results/")

library(tidyverse)
library(data.table)
library(ggplot2)
###########################################################################
#run a density plot on the MTG data

cov <- read.table("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metagenomic/woltka2_results/metaG_coverage.tsv", 
                  sep="\t", 
                  dec=".", 
                  header = TRUE, 
                  stringsAsFactors = FALSE, 
                  na.strings = NA, 
                  check.names = FALSE, 
                  comment.char = "")
cov <- cov[order(cov$coverage_ratio, decreasing = TRUE), ]
cov



ggplot(cov, aes(x=coverage_ratio)) + 
  geom_density() + xlim(c(0,0.005))#filtered at 0.001 based on this results and 
###########################################################################
#since we use zebra filter on the MTG data, I will check that the zeroes from MTG are also zeroes in MTX

#pfam_mtx<-fread("pfam/pfam_clean_noNT.tsv")%>%
pfam_mtx<-fread("genome_clean_noNT.tsv")%>%
  gather(sample_id, counts, -FeatureID)

pfam_mtg<-fread("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metagenomic/woltka2_results/filtered_metaG/genome.tsv")
#pfam_mtg<-fread("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metagenomic/woltka2_results/filtered_metaG/pfam_notnorm/pfam_clean_noNT.txt")
names(pfam_mtg) <- c("FeatureID","cFA13b","cFA13c","cFA17b","cFA17c","cFA21b","cFA21c","cFA01b","cFA01c","cFA05b","cFA05c","cFA09b","cFA09c","cFT13b",
                     "cFT13c","cFT17b","cFT17c","cFT21b","cFT21c","cFT01b","cFT01c","cFT05b","cFT05c","cFT09b","cFT09c","cNA13a","cNA13b","cNA13c","cNA17a",
                     "cNA17b","cNA17c","cNA21a","cNA21b","cNA21c","cNA01a","cNA01b","cNA01c","cNA05a","cNA05b","cNA05c","cNA09a","cNA09b","cNA09c")

pfam_mtg<-pfam_mtg%>%
  gather(sample_id, counts, -FeatureID)%>%
  filter(counts==0)%>%
  left_join(.,pfam_mtx,by=c("FeatureID","sample_id"))%>%
  mutate(unmatched=ifelse(counts.x==counts.y, "matched", "not_matched"))

pfam_mtg %>%
  group_by(unmatched) %>%
  tally()

48234/90527
16802/90527
25491/90527

#a lot of the nonzero hits in MTX are zeroes in the MTG, will run zebra filter on MTX to see if this affects this outcome

###########################################################################

cov <- read.table("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/metatranscript/woltka2_m_results/zebra_filtered/metaT_coverage.tsv", 
                  sep="\t", 
                  dec=".", 
                  header = TRUE, 
                  stringsAsFactors = FALSE, 
                  na.strings = NA, 
                  check.names = FALSE, 
                  comment.char = "")
cov <- cov[order(cov$coverage_ratio, decreasing = TRUE), ]
cov



ggplot(cov, aes(x=coverage_ratio)) + 
  geom_density() + xlim(c(0,0.0025)) #choosing 0.001, will run filtered sam thru pipeline to see how drastic of a difference using zebra filtering makes on the results
