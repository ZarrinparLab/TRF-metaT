#!/bin/bash
#SBATCH -J TRFmetaG_woltka2_notnorm
#SBATCH -n 8
#SBATCH -t 0-48:00:00
#SBATCH --mem=64G
#SBATCH --mail-user=sfloresr@ucsd.edu
#SBATCH --mail-type=all
#SBATCH --export=all
#SBATCH -o /home/sfloresr/Output_Files/%x-%N-%j.out 
#SBATCH -e /home/sfloresr/Error_Files/%x-%N-%j.err

dbdir=/projects/wol/qiyun/wol2
ext=tsv
path=/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/TRF-metaT/data/processing_files/trf_metaG_woltka_results/
path2=/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/TRF-metaT/data/genome_metaG/

# genome-based mappings
woltka classify \
	-i $path/zebra_filtered/sam_filtered/ \
	-e .bowtie2_filtered.sam.xz \
	-o $path2/genome.$ext \
	-c $dbdir/proteins/coords.txt.xz
