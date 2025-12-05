#!/bin/bash
#SBATCH --partition=DGX
#SBATCH --mem=400GB
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=112
#SBATCH --output=out/gpu_macaque_a100.log
#SBATCH --job-name=gpu_macaque
#SBATCH --gres=gpu:A100:8
#module load R/4.3.3

source /etc/profile.d/lmod.sh
module use /u/area/ntosato/scratch/timing/time_scale_final/software/modules/

module load R/4.3.3-a100
module load cutensor/2.2.0.0
module load openBLAS/0.3.29-a100
module load cuda/12.6
#source env.sh
#export R_LIBS_USER=~/scratch/r_package_dgx/

cd ..

#Rscript scripts/install.R
Rscript scripts/MacaqueBrain/gpu.R
