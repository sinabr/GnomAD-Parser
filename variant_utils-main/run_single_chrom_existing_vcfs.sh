#!/bin/bash
# SLURM job wrapper to run a single chromosome using pre-split VCFs.
# Adjust CHROM, GENE_WORKERS, MAX_GENES as needed.

#SBATCH --job-name=gnomad_one_chrom
#SBATCH --partition=model4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64GB
#SBATCH --time=06:00:00
#SBATCH --output=./logs/one_chrom_%j.out
#SBATCH --error=./logs/one_chrom_%j.err

set -euo pipefail

CHROM="${CHROM:-22}"
GENE_WORKERS="${GENE_WORKERS:-8}"
MAX_GENES="${MAX_GENES:-0}"   # 0 = all genes on chromosome

# Activate environment
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate protein

PROJ_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PROJ_DIR"

echo "Running single-chromosome existing-VCF job"
echo "Chromosome: $CHROM"
echo "Gene workers: $GENE_WORKERS"
echo "Max genes: $MAX_GENES"
echo "Start: $(date)"

python run_single_chrom_existing_vcfs.py \
    --chrom "$CHROM" \
    --gene-workers "$GENE_WORKERS" \
    --max-genes "$MAX_GENES"

echo "Done: $(date)"
