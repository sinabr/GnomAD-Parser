#!/bin/bash
#SBATCH --job-name=missense_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=170G
#SBATCH --time=10:00:00
#SBATCH --array=1-24%4
#SBATCH --output=logs/missense_gpu_%A_%a.out
#SBATCH --error=logs/missense_gpu_%A_%a.err

source ~/miniconda3/etc/profile.d/conda.sh
conda activate protein

set -euo pipefail

mkdir -p logs
mkdir -p gnomad_missense_validated

TASK=${SLURM_ARRAY_TASK_ID}
if   [ "$TASK" -ge 1 ] && [ "$TASK" -le 22 ]; then CHROM="$TASK"
elif [ "$TASK" -eq 23 ]; then CHROM="X"
elif [ "$TASK" -eq 24 ]; then CHROM="Y"
else
  echo "Bad SLURM_ARRAY_TASK_ID=$TASK"
  exit 2
fi

echo "JobID=$SLURM_JOB_ID ArrayJobID=$SLURM_ARRAY_JOB_ID TaskID=$SLURM_ARRAY_TASK_ID"
echo "Node=$SLURM_NODELIST CPUs=$SLURM_CPUS_PER_TASK MEM=$SLURM_MEM_PER_NODE"
echo "Chrom=$CHROM Start=$(date)"

python -u process_all_dataframes.py process-chrom \
  --chrom "$CHROM" \
  --input_dir gnomad_all_genes/chromosome_dataframes \
  --output_dir gnomad_missense_validated \
  --reference_dir /projects/lugoteam/protein_graphs/GnomAD-Parser/ref_data \
  --batch_size 500000

echo "Done Chrom=$CHROM End=$(date)"
