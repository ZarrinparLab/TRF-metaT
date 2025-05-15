#!/bin/bash
#SBATCH -J TRFmetat_mergeBSH
#SBATCH -n 8
#SBATCH -t 0-48:00:00
#SBATCH --mem=120G
#SBATCH --mail-user=sfloresr@ucsd.edu
#SBATCH --mail-type=all
#SBATCH --export=all
#SBATCH -o /home/sfloresr/Output_Files/%x-%N-%j.out 
#SBATCH -e /home/sfloresr/Error_Files/%x-%N-%j.err

 
# Read input files .out  matching id for technical replicates
# Output combined technical replicates into one .out file

merge_files() {
    inpath="/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/TRF-metaT/data/processing_files/trf_metaT_woltka_results/diamond_BSH/"
    line=$1
    #infile1=$inpath$line"_L001_R2.out"
    #infile2=$inpath$line"_L002_R2.out"
    #outfile=$inpath$line".R2.merged.out.xz"
    infile1=$inpath$line".R1.merged.out.xz"
    infile2=$inpath$line".R2.merged.out.xz"
    outfile=$inpath$line".cb.merged.out.xz"
    echo $infile1
    echo $infile2
    echo $outfile
    xzcat $infile1 $infile2 | xz -c - > $outfile
    #cat $infile1 $infile2 | xz -c - > $outfile
}

export -f merge_files

while read -r line
do
    parallel merge_files ::: "$line" &
done < /mnt/zarrinpar/Pynchon/Notebooks/sfloresr/TRF-metaT/data/processing_files/trf_metaT_woltka_results/file_names.txt

wait 
