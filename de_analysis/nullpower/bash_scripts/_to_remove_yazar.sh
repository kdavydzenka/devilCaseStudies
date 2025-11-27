#!/bin/bash
#SBATCH --partition=EPYC
#SBATCH --mem=200GB
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --output=out/yazar.log
#SBATCH --job-name=yazar

conda activate process
LC_ALL=C.UTF-8 Rscript run_models.R yazar
