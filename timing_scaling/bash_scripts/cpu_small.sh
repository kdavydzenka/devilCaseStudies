#!/bin/bash
#SBATCH --mem=20GB
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=112
#SBATCH --partition=H100
#SBATCH --output=out/small_cpu.log
#SBATCH --job-name=cpu_small

module use /u/area/ntosato/scratch/timing/time_scale_final/software/modules/

module load R/4.3.3-h100
module load cutensor/2.2.0.0
module load openBLAS/0.3.29-h100
module load cuda/12.6
#source env.sh
cd ..
Rscript scripts/BaronPancreas/cpu.R
