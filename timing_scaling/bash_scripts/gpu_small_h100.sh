#!/bin/bash
#SBATCH --partition=H100
#SBATCH --mem=40GB
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=112
#SBATCH --output=out/small_gpu_h100_1gpus.log
#SBATCH --job-name=gpu_small
#SBATCH --gpus=1
#module load R/4.3.3

module use /u/area/ntosato/scratch/timing/time_scale_final/software/modules/

module load R/4.3.3-h100
module load cutensor/2.2.0.0
module load openBLAS/0.3.29-h100
module load cuda/12.6
#source env.sh
#export R_LIBS_USER=~/scratch/r_package_dgx/

cd ..

Rscript scripts/install.R
Rscript scripts/BaronPancreas/gpu.R
