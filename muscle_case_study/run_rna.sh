#!/bin/bash
#SBATCH --partition=THIN
#SBATCH --mem=200GB
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --output=out.log
#SBATCH --job-name=de_analysis

module load R/4.4.1

Rscript 01_fitData.R

