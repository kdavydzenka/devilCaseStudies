#!/bin/bash
#SBATCH --partition=H100
#SBATCH --mem=500GB
#SBATCH --time=18:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=112
#SBATCH --output=gpu_macaque_h100.log
#SBATCH --job-name=gpu_macaque
#SBATCH --gpus=8


#export R_LIBS_USER=/u/area/ntosato/scratch/devil2025/r-pkg
module use /u/area/ntosato/scratch/timing/time_scale_final/software/modules/

module load R/4.3.3-h100
module load cutensor/2.2.0.0
module load openBLAS/0.3.29-h100
module load cuda/12.6
#source env.sh
#export R_LIBS_USER=~/scratch/r_package_dgx/

cd ..

Rscript scripts/install.R
Rscript scripts/MacaqueBrain/gpu.R
