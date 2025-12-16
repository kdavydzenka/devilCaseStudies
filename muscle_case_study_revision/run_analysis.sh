#!/bin/bash
#SBATCH --partition=GENOA
#SBATCH --nodes=1
#SBATCH --job-name=downstream_array
#SBATCH --output=logs/downstream_%A_%a.out
#SBATCH --error=logs/downstream_%A_%a.err
#SBATCH --array=0-7
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=24GB

# -------------------------
# Define parameter grids
# -------------------------
SF_LIST=("subsampled" "full")
COND_LIST=("age_only" "age_type1" "age_type2" "interaction")

# -------------------------
# Compute combinations
# -------------------------
# sf index = array_id % 2
# condition index = array_id / 2
SF_INDEX=$(( SLURM_ARRAY_TASK_ID % 2 ))
COND_INDEX=$(( SLURM_ARRAY_TASK_ID / 2 ))

SF=${SF_LIST[$SF_INDEX]}
COND=${COND_LIST[$COND_INDEX]}

echo "----------------------------------------"
echo " SLURM ARRAY JOB: $SLURM_ARRAY_TASK_ID"
echo " Running with:"
echo "   sf       = $SF"
echo "   condition = $COND"
echo "----------------------------------------"

# -------------------------
# Run the R script
# Pass parameters as command-line args
# -------------------------

source ~/.bashrc
conda activate process

Rscript 02_per_type_analysis.R --sf "$SF" --condition "$COND"
