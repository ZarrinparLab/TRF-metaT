#!/bin/bash
#SBATCH -J BSH_diamond_metaG
#SBATCH -n 8
#SBATCH -t 0-72:00:00
#SBATCH --mem=320G
#SBATCH --mail-user=sfloresr@ucsd.edu
#SBATCH --mail-type=all
#SBATCH --export=all
#SBATCH -o /home/sfloresr/Output_Files/%x-%N-%j.out
#SBATCH -e /home/sfloresr/Error_Files/%x-%N-%j.err

#diamond makedb --in /projects/zl_trf_metat/trf_metaT_woltka_results/wol2_wrep_notnorm/Uniprot_6K_BSH_genes.fasta.gz -d /projects/zl_trf_metat/trf_metaG_woltka_results/wol2/filtered_notnorm/BSH_db

bt2db=/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/TRF-metaT/data/processing_files/trf_metaG_woltka_results/wol2/filtered_notnorm/BSH_db/BSH_db.dmnd
path=/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/TRF-metaT/data/processing_files/trf_metaG_woltka_results

sed -n "${SLURM_ARRAY_TASK_ID}p" /mnt/zarrinpar/Pynchon/Notebooks/sfloresr/TRF-metaT/data/processing_files/trf_metaG_woltka_results/zebra_filtered/sam_filtered/summ_fastq_zf.txt | while read id fq1 fq2 id1 id2; do

diamond blastx -d $bt2db -q $path/zebra_filtered/sam_filtered/$fq1 -o $path/diamond_BSH_zebra/${id1}.out
diamond blastx -d $bt2db -q $path/zebra_filtered/sam_filtered/$fq2 -o $path/diamond_BSH_zebra/${id2}.out
done
