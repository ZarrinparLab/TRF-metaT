#!/bin/bash
#SBATCH -J TRFmetaG_bowtie2
#SBATCH -n 8
#SBATCH -t 0-48:00:00
#SBATCH --mem=96G
#SBATCH --mail-user=sfloresr@ucsd.edu
#SBATCH --mail-type=all
#SBATCH --export=all
#SBATCH -o /home/sfloresr/Output_Files/%x-%N-%j.out
#SBATCH -e /home/sfloresr/Error_Files/%x-%N-%j.err

#bowtie2=/home/drz/Programs/Bowtie2/2.4.5/bowtie2
db=/projects/wol/qiyun/wol2/databases
bt2db=$db/bowtie2/WoLr2
path=/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/TRF-metaT/data/processing_files/20210829_trf_metaG_assembly/
path2=/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/TRF-metaT/data/processing_files/trf_metaG_woltka_results/

#script to run shogun on the metaG data 
read_list='/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/TRF-metaT/data/processing_files/trf_metaG_woltka_results/filtered_fastq_list.txt'
while read line; do
echo "'$line'"

bowtie2 -p 8 -x $bt2db --very-sensitive -1 $path/filtered.fastq/$line'_L001_R1_001.filtered.fastq.gz' -2 $path/filtered.fastq/$line'_L001_R2_001.filtered.fastq.gz' --seed 42 --no-head --no-unal | cut -f1-9 | sed 's/$/\t*\t*/' | xz -9 > $path2/bowtie2/$line'.bowtie2.sam.xz'

bowtie2 -p 8 -x $bt2db --very-sensitive -1 $path/filtered.fastq/$line'_L001_R1_001.filtered.fastq.gz' -2 $path/filtered.fastq/$line'_L001_R2_001.filtered.fastq.gz' --seed 42 --no-head --no-unal -k 16 --np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" --score-min "L,0,-0.05" | cut -f1-9 | sed 's/$/\t*\t*/' | xz -9 > $path2/bowtie2/$line'.bt2sho.sam.xz'
done < $read_list
