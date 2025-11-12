#!/usr/bin/env bash
#SBATCH --job-name=PBMC_run
#SBATCH --output=logs/pbmc.out
#SBATCH --error=logs/pbmc.err
#SBATCH --partition=GENOA
#SBATCH --mem=128GB
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4

source ~/.bashrc
conda activate process
Rscript run_purifiedPBMC.R