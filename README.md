# TRF-metaT

This is the code for the analysis of the microbial composition (16S), functional potential (functional metagenomics), and gene expression (metatranscriptomics) of mice exposed to high-fat diet (HFD) and treated with the intervetion time-restricted feeding (TRF). The goal was to characterize the microbial activity of microbes under HFD and TRF and use what we learn of microbial expression to guide the engineering of new bacterial therapeutics.

## Instructions on Data Access

### Multi-omics

We used samples collected and phenotyped from [Zarrinpar, et al. 2014](https://www.cell.com/cell-metabolism/fulltext/S1550-4131(14)00505-1). The 16S data was sequenced as part of this earlier study. Shotgun metagenomics and metatranscriptomics are available in ENA under the project ID PRJEB89098. Untargeted metabolomics data collected on the new engineered microbes is available in the MassIVE database [https://massive.ucsd.edu/] under the MassIVE ID MSV000094578 (culture) and MSV000097414 (fecal). Feature based molecular networking data is available in GNPS2 under https://gnps2.org/status?task=f67a6015e87942fab913beedd63a36aa and https://gnps2.org/status?task=092d720b607e4957ab6b2cd67d15d98f.

## Organization of files

The main folder are `code`, `data`, and `figures`. Within each folder, it is organized as follows:

- `diversity_analysis`: contains the code, data, or figures for the diversity analysis comparing 16S, metagenomics genes, and metatranscriptomic transcripts. Related to Figure 1.
- `DE_analysis`: contains the code, data, or figures for the differential abundance analysis for the metagenomics genes and metatranscriptomic transcripts. Related to Figure 2.
- `cycling_analysis`: contains the code, data, or figures for the cycling analysis for the metagenomics genes and metatranscriptomic transcripts. Related to Figure 3.
- `bsh_analysis`: contains the code, data, or figures for the bsh differential expression analysis using the targeted and untargeted search for bsh in the metatranscriptomics. Related to Figure 4.
- `ENB_invitro_metab`: contains the code, data, or figures the bile acid quantification from the cultures of ENBs. Related to Figure 5.
- `ENB_invivo_metab`: contains the code, data, or figures for the phenotypic measurements of mice gavaged with ENBs. Also contains the fecal bile acid quantification of these mice. Related to Figure 6.
- `replotting_AZ2014_CD2022`: contains the code, data, or figures for weight and glucose measurements from Zarrinpar 2014 and ileum bile acid quantification from Dantas Machado, 2022 we replotted in Figure S1 and S4, respectively.

Under `code` we also include:
- `processing_metaG`: contains the code used to process the metagenomic data from fastq to counts table.
- `processing_metaT`: contains the code used to process the metatranscriptomics data from fastq to counts table.

Under `data` we also include:
- `m16s`: contains the processed 16s counts data from Zarrinpar 2014
- `genome_metaG`: contains woltka processed metagenomic counts data at the genome level
- `genome_metaT`: contains woltka processed metatranscriptomics counts data at the genome level
- `pfam_metaG`: contains woltka processed metagenomic counts data collapsed to pfam functions
- `pfam_metaT`: contains woltka processed metatranscriptomics counts data collapsed to pfam functions

## Citation

Flores Ramos, S., Siguenza, N., Zhong, W., Mohanty, I., Lingaraju, A., Richter, R.A., Karthikeyan, S., Lukowski, A.L., Zhu, Q., Nunes, W.D.G., Zemlin, J., Xu, Z.Z., Hasty, J., Dorrestein, P.C., Panda, S., Knight, R., Zarrinpar, A. “Metatranscriptomics Uncover Diurnal Functional Shifts in Bacterial Transgenes with Profound Metabolic Effects”. Cell Host & Microbe (2025)
