#!/usr/bin/env bash
set -euo pipefail

AUTHOR="${1:-bca}"   # pass author on the command line: bca|yazar|hsc|kumar

mkdir -p logs

sbatch <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=de-sim-${AUTHOR}
#SBATCH --output=logs/prepping_%x.out
#SBATCH --error=logs/prepping_%x.err
#SBATCH --partition=EPYC
#SBATCH --mem=200GB
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1

source ~/.bashrc
conda activate process

echo "Author: ${AUTHOR}"

Rscript 00_prep_data_v2.R "${AUTHOR}"
EOF
