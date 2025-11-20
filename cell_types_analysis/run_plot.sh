#!/bin/bash
#SBATCH --partition=THIN
#SBATCH --mem=100GB
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --output=out/plot.log
#SBATCH --job-name=cell_type_plotting

source ~/.bashrc
conda activate process

LC_ALL=C.UTF-8 Rscript 03_supp_figure.R BaronPancreasData pancreas
LC_ALL=C.UTF-8 Rscript 03_supp_figure.R pbmc blood
LC_ALL=C.UTF-8 Rscript 03_supp_figure.R liver liver
LC_ALL=C.UTF-8 Rscript 03_supp_figure.R BigBloodData blood
LC_ALL=C.UTF-8 Rscript 03_supp_figure.R BigLiverData liver
