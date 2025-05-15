#!/bin/bash
#SBATCH -J TRFmetat_woltka2_notnorm
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
path=/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/TRF-metaT/data/processing_files/trf_metaT_woltka_results/
path2=/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/TRF-metaT/data/geome_metaT/

# genome-based mappings
woltka classify \
	-i $path/bowtie2/ \
	-e .merged.bowtie2.sam.xz \
	-o $path2/genome.$ext \
	-c $dbdir/proteins/coords.txt.xz
