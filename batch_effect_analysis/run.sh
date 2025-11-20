#!/usr/bin/env bash
#SBATCH --job-name=batch_sim
#SBATCH --output=logs/%A_%a.out
#SBATCH --error=logs/%A_%a.err
#SBATCH --partition=THIN
#SBATCH --mem=128GB
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --array=1-320%40

source ~/.bashrc
conda activate process

echo "SLURM_ARRAY_TASK_ID=\$SLURM_ARRAY_TASK_ID"

Rscript fit_data.R
