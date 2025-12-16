#!/bin/bash
#SBATCH --partition=GENOA
#SBATCH --job-name=muscle_de
#SBATCH --array=1-12
#SBATCH --output=logs/muscle_de_%A_%a.out
#SBATCH --error=logs/muscle_de_%A_%a.err
#SBATCH --mem=200GB
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4

source ~/.bashrc
conda activate process

# Define the combinations
methods=("devil" "glmGamPoi" "nebula")
designs=("age_only" "interaction")
subsamples=("FALSE" "TRUE")

# Total combinations
n_methods=${#methods[@]}      # 3
n_designs=${#designs[@]}      # 2
n_subsamples=${#subsamples[@]} # 2

# Convert SLURM_ARRAY_TASK_ID (1..12) to zero-based index
idx=$((SLURM_ARRAY_TASK_ID - 1))

# Decode index into (subsample, design_test, method)
combo_per_sub=$((n_methods * n_designs))   # 6
sub_idx=$((idx / combo_per_sub))
rem=$((idx % combo_per_sub))
des_idx=$((rem / n_methods))
met_idx=$((rem % n_methods))

METHOD=${methods[$met_idx]}
DESIGN_TEST=${designs[$des_idx]}
SUBSAMPLE=${subsamples[$sub_idx]}

echo "Running combination:"
echo "  METHOD      = ${METHOD}"
echo "  DESIGN_TEST = ${DESIGN_TEST}"
echo "  SUBSAMPLE   = ${SUBSAMPLE}"

# Run the R script with arguments
Rscript 01_fitData.R "${METHOD}" "${DESIGN_TEST}" "${SUBSAMPLE}"
