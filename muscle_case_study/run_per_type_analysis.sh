#!/bin/bash
#SBATCH --partition=GENOA
#SBATCH --mem=50GB
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --output=analysis.log
#SBATCH --job-name=muscle_analysis

source ~/.bashrc
conda activate process

Rscript per_type_clustering_analysis_revision.R
