#!/bin/bash
#SBATCH -J TRFmetaG_woltka2_BSH
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
	-i /projects/zl_trf_metat/trf_metaG_woltka_results/diamond_BSH/ \
	-e .merged.out.xz \
        --format b6o \
        --no-demux \
	-o /mnt/zarrinpar/Pynchon/Notebooks/sfloresr/TRF-metaT/data/bsh_analysis/BSH_proteindb_metaG/genome.$ext
