#!/bin/bash
#SBATCH --mem=20GB
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=GENOA
#SBATCH --output=out/small_cpu_all.log
#SBATCH --job-name=cpu_small_all

source ~/.bashrc
conda activate process

Rscript scripts/BaronPancreas/cpu.R
# Rscript scripts/BaronPancreas/cpu_nebula.R
