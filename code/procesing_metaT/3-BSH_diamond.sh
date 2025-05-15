#!/bin/bash
#SBATCH -J BSH_diamond
#SBATCH -n 8
#SBATCH -t 0-72:00:00
#SBATCH --mem=320G
#SBATCH --mail-user=sfloresr@ucsd.edu
#SBATCH --mail-type=all
#SBATCH --export=all
#SBATCH -o /home/sfloresr/Output_Files/%x-%N-%j.out
#SBATCH -e /home/sfloresr/Error_Files/%x-%N-%j.err

path=/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/TRF-metaT/data/processing_files/trf_metaT_woltka_results/
path2=/projects/zl_trf_metat/20201015_trf_metatranscriptome

#diamond makedb --in $path/wol2_wrep_notnorm/Uniprot_6K_BSH_genes.fasta.gz -d $path/wol2_wrep_notnorm/BSH_db

bt2db=$path/wol2_wrep_notnorm/BSH_db/BSH_db.dmnd

sed -n "${SLURM_ARRAY_TASK_ID}p" $path2/cleaned/cleaned_files/fastq_files_clean.txt | while read id fq1 fq2 id1 id2; do

diamond blastx -d $bt2db -q $path2/cleaned/cleaned_files/$fq1 -o $path/diamond_BSH/${id1}.out
diamond blastx -d $bt2db -q $path2/cleaned/cleaned_files/$fq2 -o $path/diamond_BSH/${id2}.out
done
