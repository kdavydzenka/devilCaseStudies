#!/usr/bin/env bash
set -euo pipefail

AUTHOR="${1:-bca}"   # pass author on the command line: bca|yazar|hsc|kumar

mkdir -p logs

sbatch <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=de-sim-${AUTHOR}
#SBATCH --output=logs/plotting.out
#SBATCH --error=logs/plotting.err
#SBATCH --partition=EPYC
#SBATCH --mem=200GB
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1

source ~/.bashrc
conda activate process

Rscript plot_results.R 
EOF
