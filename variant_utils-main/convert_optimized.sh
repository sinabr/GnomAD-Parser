#!/bin/bash
#SBATCH --job-name=gnomad_all
#SBATCH --partition=pool1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=3-00:00:00
#SBATCH --array=1-24%10
#SBATCH --output=logs/gnomad_%A_%a.out
#SBATCH --error=logs/gnomad_%A_%a.err



source ~/miniconda3/etc/profile.d/conda.sh
conda activate protein

set -euo pipefail

mkdir -p logs
mkdir -p gnomad_all_genes/chromosome_dataframes

# Map array index -> chromosome
# 1..22 => chr1..chr22, 23 => X, 24 => Y
TASK=${SLURM_ARRAY_TASK_ID}
if [ "$TASK" -le 22 ]; then
  CHR="$TASK"
elif [ "$TASK" -eq 23 ]; then
  CHR="X"
else
  CHR="Y"
fi

echo "============================================================"
echo "Job: $SLURM_JOB_ID  Task: $SLURM_ARRAY_TASK_ID  Node: $SLURM_NODELIST"
echo "CPUs: $SLURM_CPUS_PER_TASK  Mem: ${SLURM_MEM_PER_NODE:-NA}"
echo "Chromosome: chr${CHR}"
echo "Start: $(date)"
echo "============================================================"

python -u convert_optimized.py \
  --chrom "${CHR}" \
  --chunk-rows 500000 \
  --threads "${SLURM_CPUS_PER_TASK}"

echo "============================================================"
echo "Done chr${CHR} at $(date)"
echo "Output:"
ls -lh "gnomad_all_genes/chromosome_dataframes/chr${CHR}_variants.parquet" || true
echo "============================================================"
