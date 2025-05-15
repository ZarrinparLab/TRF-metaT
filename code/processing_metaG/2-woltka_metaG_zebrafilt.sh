#!/bin/bash
#BATCH -J TRFmetaG_zebrafilt
#SBATCH -n 8
#SBATCH -t 0-48:00:00
#SBATCH --mem=240G
#SBATCH --mail-user=sfloresr@ucsd.edu
#SBATCH --mail-type=all
#SBATCH --export=all
#SBATCH -o /home/sfloresr/Output_Files/%x-%N-%j.out
#SBATCH -e /home/sfloresr/Error_Files/%x-%N-%j.err

db=/home/sfloresr/zebra_filter/databases/WoL/metadata.tsv
path=/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/TRF-metaT/data/processing_files/trf_metaG_woltka_results/

#python /home/sfloresr/zebra_filter/calculate_coverages.py \
# -i $path/bowtie2/ \
# -o $path/zebra_filtered/metaG_coverage.tsv \
# -d $db

python /home/sfloresr/zebra_filter/filter_sam.py \
  -i $path/zebra_filtered/metaG_coverage.tsv \
  -s $path/bowtie2/ \
  -c .001 \
  -o $path/zebra_filtered/sam_filtered/
