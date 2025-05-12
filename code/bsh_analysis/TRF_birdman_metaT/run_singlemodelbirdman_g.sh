#!/bin/bash
#SBATCH --chdir=/projects/zl_trf_metat/
#SBATCH --output=/projects/zl_trf_metat/scripts/woltka/TRF_birdman/slurm/%x.%a.out
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

TABLEID="genomeL_noNT_rmdbton"
TABLE="/projects/zl_trf_metat/trf_metaT_woltka_results/wol2_wrep_notnorm/BSH_results/targeted_BSH/"$TABLEID".biom"
SLURMS="/projects/zl_trf_metat/trf_metaT_woltka_results/wol2_wrep_notnorm/BSH_results/targeted_BSH/slurm_out/"$TABLEID
OUTDIR="/projects/zl_trf_metat/trf_metaT_woltka_results/wol2_wrep_notnorm/BSH_results/targeted_BSH/inferences/"$TABLEID
LOGDIR="/projects/zl_trf_metat/trf_metaT_woltka_results/wol2_wrep_notnorm/BSH_results/targeted_BSH/logs/"$TABLEID
mkdir -p $SLURMS
mkdir -p $OUTDIR
mkdir -p $LOGDIR

echo Starting Python script...
time python /projects/zl_trf_metat/scripts/woltka/TRF_birdman/TRF_BSH_birdman_chunked.py \
    --table-path $TABLE \
    --inference-dir $OUTDIR \
    --num-chunks $SLURM_ARRAY_TASK_MAX \
    --chunk-num $SLURM_ARRAY_TASK_ID \
    --logfile "${LOGDIR}/chunk_${SLURM_ARRAY_TASK_ID}.log" && echo Finished Python script!

