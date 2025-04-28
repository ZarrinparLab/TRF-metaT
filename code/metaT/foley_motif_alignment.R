setwd("~/scratch/TRF_multiomics/metatranscript/woltka2_m_results")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("msa")