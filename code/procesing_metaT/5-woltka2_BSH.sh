#!/bin/bash
#SBATCH -J TRFmetat_woltka2_BSH
#SBATCH -n 8
#SBATCH -t 0-48:00:00
#SBATCH --mem=96G
#SBATCH --mail-user=sfloresr@ucsd.edu
#SBATCH --mail-type=all
#SBATCH --export=all
#SBATCH -o /home/sfloresr/Output_Files/%x-%N-%j.out 
#SBATCH -e /home/sfloresr/Error_Files/%x-%N-%j.err

dbdir=/projects/wol/qiyun/wol2
ext=tsv

# genome-based mappings
woltka classify \
	-i /mnt/zarrinpar/Pynchon/Notebooks/sfloresr/TRF-metaT/data/processing_files/trf_metaT_woltka_results/diamond_BSH/ \
	-e .cb.merged.out.xz \
        --format b6o \
        --no-demux \
	-o /mnt/zarrinpar/Pynchon/Notebooks/sfloresr/TRF-metaT/data/bsh_analysis/BSH_proteindb_metaT/genome.$ext
