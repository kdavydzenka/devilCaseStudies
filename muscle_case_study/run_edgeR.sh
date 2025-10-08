#!/bin/bash
#SBATCH --partition=cpuq
#SBATCH --mem=200GB
#SBATCH --time=200:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --output=out.log
#SBATCH --job-name=de_analysis

module load R/4.2.3

LC_ALL=C.UTF-8 Rscript extra.R
