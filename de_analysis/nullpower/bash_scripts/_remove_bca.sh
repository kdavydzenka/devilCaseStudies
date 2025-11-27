#!/bin/bash
#SBATCH --partition=EPYC
#SBATCH --mem=128GB
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --output=out/bca.log
#SBATCH --job-name=bca

conda activate process
LC_ALL=C.UTF-8 Rscript run_models.R bca
