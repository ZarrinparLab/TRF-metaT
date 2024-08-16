setwd("/mnt/zarrinpar/scratch/sfloresr/TRF_multiomics/bileacid_fromTI/")

library(tidyverse)
library(MetaCycle)
library(data.table)
library("qiime2R")
library("Biostrings")
library(ggrepel)
library(ggpubr)
library(gplots)
library(ggvenn)
library(ggbreak)
library(viridis)
###########################################################
infile <- "TI_bileacids_abs_quant_mgtissue_20211122.csv"
salksamplesheet <- "TI_bileacids_sample_sheet_20211122.csv"
mdatf <- "TI_metadata_final.csv"
fig_pre <- "bileacids_"

# Set colors/shapes
mycols <- c("#0072B2","#D55E00","#009E73")
cols_by_group <- c(`NA`="#0072B2", FA="#D55E00", FT="#009E73")

# Functions

merge_log_dat_w_mdat <- function (df, x, metadata=mdat) {
  
  df <- merge(df, metadata, by.x="sample_name", by.y="sample", all.x=T, all.y=F)
  df <- pivot_longer(df, cols=x, names_to="bile_acid", values_to="log2")
  df$Acronym <-  factor(df$Acronym, levels=c("NA", "FA", "FT"))
  df$phase <- factor(df$phase, levels = c("light", "dark"))  
  return(df)
}

# Read and process input data ------------------------

# Read input files
dat <- read.csv(infile, row.names = 1, check.names = F)
samples <- read.csv(salksamplesheet)
mdat <- read.csv(mdatf, na.strings = "")
mdat$Acronym <- as.factor(mdat$Acronym)

# Pre-process main bile acid names and data
bas <- colnames(dat)
dat$sample_name <- gsub("ileum_.*\\s", "", rownames(dat))
fdat <- merge(dat, mdat, by.x="sample_name", by.y="sample", all.x=T, all.y=F)
fdat <- pivot_longer(fdat, cols=bas, names_to="bile_acid", values_to="pmol/mg")
fdat$Acronym <-  factor(fdat$Acronym, levels=c("NA", "FA", "FT"))
fdat$phase <- factor(fdat$phase, levels = c("light", "dark"))

# Get ratios from original data table 
dat$`CA/TCA` <- log(dat$CA/dat$TCA, base = 2)
dat$`b-MCA/T-b-MCA` <- log(dat$`b-MCA`/dat$`T-b-MCA`, base = 2)
dat$`a-MCA/T-a-MCA` <- log(dat$`a-MCA`/dat$`T-a-MCA`, base = 2)
dat$`o-MCA/T-o-MCA` <- log(dat$`o-MCA`/dat$`T-o-MCA`, base = 2)
dat$`DCA/TDCA` <- log(dat$`DCA`/dat$`TDCA`, base = 2)
dat$`CDCA/TCDCA` <- log(dat$`CDCA`/dat$`TCDCA`, base = 2)


# Get data with ratios only
ba_ratios <- c("CA/TCA", "a-MCA/T-a-MCA", "b-MCA/T-b-MCA", "o-MCA/T-o-MCA", "DCA/TDCA", "CDCA/TCDCA")
bdat <- merge_log_dat_w_mdat(df = dat[,c(ba_ratios, "sample_name")], x=ba_ratios)



# Plot BAs ---------------------------


# Plot ratios (Figure 7C)
bdat %>%
  ggplot(aes(x=Acronym, y=log2, color = Acronym, fill=Acronym, shape=phase)) +
  theme_bw()+
  theme(legend.position = "bottom",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  geom_boxplot(alpha=0.2,position = "dodge") +
  geom_point(position = position_jitterdodge()) +
  stat_compare_means(method="wilcox", size=2, label = "p.signif" ) +
  facet_wrap(~bile_acid, scales = "free_y", nrow=1) +
  scale_shape_manual(values=c(3,16)) +
  scale_fill_manual(values = mycols) +
  scale_color_manual(values = mycols)
ggsave(filename ="boxplot_phase_ratios.pdf", width = 9, height = 2.5)


# Plot unconjugated - S7
subset(fdat, bile_acid %in% c("a-MCA", "b-MCA", "CA", "CDCA", "DCA", "o-MCA")) %>%
  ggplot(aes(y=`pmol/mg`, x=Acronym,  color=Acronym, fill=phase, shape=phase)) +
  geom_boxplot(alpha=0) +
  theme_bw()+
  theme(legend.position = "bottom") +
  geom_point(position = position_jitterdodge()) +
  facet_wrap(~bile_acid, nrow=1, scales = "free_y") +
  stat_compare_means(method="wilcox", size=2, label = "p.signif" ) +
  scale_shape_manual(values=c(5,8)) +
  scale_fill_manual(values = mycols) +
  scale_color_manual(values = mycols)
ggsave(filename = paste0(fig_pre, "boxplot_unconj_group_phase_comparison.pdf" ), width = 9, height = 2.5)

# Plot conjugated - S7
subset(fdat, bile_acid %in% c("T-a-MCA", "T-b-MCA", "TCA", "TCDCA", "TDCA", "T-o-MCA")) %>%
  ggplot(aes(y=`pmol/mg`, x=Acronym,  color=Acronym, fill=phase, shape=phase)) +
  geom_boxplot(alpha=0) +
  theme_bw()+
  theme(legend.position = "bottom") +
  geom_point(position = position_jitterdodge()) +
  facet_wrap(~bile_acid, nrow=1, scales = "free_y") +
  stat_compare_means(method="wilcox", size=2, label = "p.signif" ) +
  scale_shape_manual(values=c(5,8)) +
  scale_fill_manual(values = mycols) +
  scale_color_manual(values = mycols)
ggsave(filename = paste0(fig_pre, "boxplot_conj_group_phase_comparison.pdf" ), width = 9, height = 2.5)



# Check significant hits
dl_bygroup <- compare_means(data = fdat, formula = `pmol/mg` ~ phase,
                            p.adjust.method = "fdr",
                            group.by = c("bile_acid", "Acronym"), method="wilcox") 

# Subset dl significant results
dl_sig <- subset(dl_bygroup, p < 0.05, select=c("bile_acid", "Acronym"))
dl_sig