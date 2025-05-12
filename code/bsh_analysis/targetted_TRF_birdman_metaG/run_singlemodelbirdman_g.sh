#!/bin/bash
#SBATCH --chdir=/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/TRF-metaT/data/bsh_analysis/BSH_proteindb_metaG/birdman_outputs/
#SBATCH --output=/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/TRF-metaT/data/bsh_analysis/BSH_proteindb_metaG/birdman_outputs/slurm/%x.%a.out
#SBATCH --partition=short
#SBATCH --mail-user="sfloresr@ucsd.edu"
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mem=64G
#SBATCH --nodes=1
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --time=6:00:00
#SBATCH --array=1-10

pwd; hostname; date

set -e

source ~/anaconda3/bin/activate birdman

echo Chunk $SLURM_ARRAY_TASK_ID / $SLURM_ARRAY_TASK_MAX

TABLEID="genome_noNT_rmdbton_NA"
TABLE="/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/TRF-metaT/data/bsh_analysis/BSH_proteindb_metaG/"$TABLEID".biom"
SLURMS="/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/TRF-metaT/data/bsh_analysis/BSH_proteindb_metaG/birdman_outputs/slurm_out/"$TABLEID
OUTDIR="/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/TRF-metaT/data/bsh_analysis/BSH_proteindb_metaG/birdman_outputs/inferences/"$TABLEID
LOGDIR="/mnt/zarrinpar/Pynchon/Notebooks/sfloresr/TRF-metaT/data/bsh_analysis/BSH_proteindb_metaG/birdman_outputs/logs/"$TABLEID
mkdir -p $SLURMS
mkdir -p $OUTDIR
mkdir -p $LOGDIR

echo Starting Python script...
time python /mnt/zarrinpar/Pynchon/Notebooks/sfloresr/TRF-metaT/code/bsh_analysis/targetted_TRF_birdman_metaG/TRF_BSH_birdman_chunked.py \
    --table-path $TABLE \
    --inference-dir $OUTDIR \
    --num-chunks $SLURM_ARRAY_TASK_MAX \
    --chunk-num $SLURM_ARRAY_TASK_ID \
    --logfile "${LOGDIR}/chunk_${SLURM_ARRAY_TASK_ID}.log" && echo Finished Python script!

