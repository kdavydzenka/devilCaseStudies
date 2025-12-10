#!/usr/bin/env bash
#SBATCH --job-name=de-sim-plotting
#SBATCH --output=logs/plotting.out
#SBATCH --error=logs/plotting.err
#SBATCH --partition=GENOA
#SBATCH --mem=50GB
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1

source ~/.bashrc
conda activate process

Rscript 02_summarise_results.R
Rscript 03_plot_results.R
