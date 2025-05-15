#!/bin/bash
#SBATCH -J BSH_convertsamtofastq_metaG
#SBATCH -n 8
#SBATCH -t 0-72:00:00
#SBATCH --mem=320G
#SBATCH --mail-user=sfloresr@ucsd.edu
#SBATCH --mail-type=all
#SBATCH --export=all
#SBATCH -o /home/sfloresr/Output_Files/%x-%N-%j.out
#SBATCH -e /home/sfloresr/Error_Files/%x-%N-%j.err

convertsamtofastq_files() {
    inpath="/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/TRF-metaT/data/processing_files/trf_metaG_woltka_results/zebra_filtered/sam_filtered/"
    line=$1
    infile=$inpath$line".bowtie2_filtered.sam.xz"
    outfile1=$inpath$line".bowtie2_filtered.R1.fastq.gz"
    outfile2=$inpath$line".bowtie2_filtered.R2.fastq.gz"
    echo "Processing: $infile -> $outfile"
    xzcat "$infile" | samtools view -h -t "/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/TRF-metaT/data/processing_files/trf_metaG_woltka_results/zebra_filtered/sam_filtered/all.fna.faidx" -F 0x0C | samtools fastq - -1 "$outfile1" -2 "$outfile2"
}

export -f convertsamtofastq_files

parallel convertsamtofastq_files :::: /mnt/zarrinpar/Pynchon/Notebooks/sfloresr/TRF-metaT/data/processing_files/trf_metaG_woltka_results/file_names.txt
