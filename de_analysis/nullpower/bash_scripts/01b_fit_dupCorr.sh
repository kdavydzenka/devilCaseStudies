#!/usr/bin/env bash
set -euo pipefail

AUTHOR="${1:-bca}"   # pass author on the command line: bca|yazar|hsc|kumar

mkdir -p logs

sbatch <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=de-sim-${AUTHOR}
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --partition=GENOA
#SBATCH --mem=32GB
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --array=1-360%60 

source ~/.bashrc
conda activate process

echo "Author: ${AUTHOR}"
echo "SLURM_ARRAY_TASK_ID=\$SLURM_ARRAY_TASK_ID"

Rscript 01b_fit_data_dupCorr.R "${AUTHOR}"
EOF
