#!/bin/bash
#SBATCH -J TRFmetat_bowtie2_single
#SBATCH -n 8
#SBATCH -t 0-48:00:00
#SBATCH --mem=240G
#SBATCH --mail-user=sfloresr@ucsd.edu
#SBATCH --mail-type=all
#SBATCH --export=all
#SBATCH -o /home/sfloresr/Output_Files/%x-%N-%j.out
#SBATCH -e /home/sfloresr/Error_Files/%x-%N-%j.err

#bowtie2=/home/drz/Programs/Bowtie2/2.4.5/bowtie2
db=/projects/wol/qiyun/wol2/databases
bt2db=$db/bowtie2/WoLr2
path=/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/TRF-metaT/data/processing_files/20201015_trf_metatranscriptome/
path2=/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/TRF-metaT/data/processing_files/trf_metaT_woltka_results/

sed -n "${SLURM_ARRAY_TASK_ID}p" $path/cleaned/cleaned_files/fastq_files_v2.txt | while read id fq1 fq2; do

bowtie2 -p 8 -x $bt2db --very-sensitive -1 $path/cleaned/cleaned_files/$fq1 -2 $path/cleaned/cleaned_files/$fq2 --seed 42 --no-head --no-unal | cut -f1-9 | sed 's/$/\t*\t*/' | xz -9 > $path2/bowtie2/${id}.bowtie2.sam.xz | 2> $path2/bowtie2/${id}.bowtie2.log

bowtie2 -p 8 -x $bt2db --very-sensitive -1 $path/cleaned/cleaned_files/$fq1 -2 $path/cleaned/cleaned_files/$fq2 --seed 42 --no-head --no-unal -k 16 --np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" --score-min "L,0,-0.05" | cut -f1-9 | sed 's/$/\t*\t*/' | xz -9 > $path2/bowtie2/${id}.bt2sho.sam.xz
done
