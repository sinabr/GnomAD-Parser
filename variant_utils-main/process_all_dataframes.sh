#!/bin/bash
#SBATCH --job-name=process_missense
#SBATCH --partition=model4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=logs/process_missense_%j.out
#SBATCH --error=logs/process_missense_%j.err

# ============================================================================
# Process Chromosome Parquet Files → Validated Missense Variants
# ============================================================================
# This script reads the chromosome-level Parquet files and produces
# validated missense variants with protein sequences.
#
# Memory: 64GB is sufficient for loading all chromosome data in batches
# Time: ~2-3 hours for all chromosomes
# ============================================================================

source ~/miniconda3/etc/profile.d/conda.sh
conda activate protein

set -euo pipefail

mkdir -p logs
mkdir -p gnomad_missense_validated

echo "============================================================================"
echo "GNOMAD MISSENSE VARIANT PROCESSING"
echo "============================================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: $SLURM_MEM_PER_NODE MB"
echo "Start Time: $(date)"
echo "============================================================================"
echo ""

# Run the processing script
python -u process_missense_from_parquet.py

EXIT_CODE=$?

echo ""
echo "============================================================================"
echo "PROCESSING COMPLETE"
echo "============================================================================"
echo "Exit code: $EXIT_CODE"
echo "End time: $(date)"
echo ""

if [ $EXIT_CODE -eq 0 ]; then
    echo "✅ Success! Results saved to gnomad_missense_validated/"
    echo ""
    echo "Output files:"
    ls -lh gnomad_missense_validated/*.parquet 2>/dev/null || echo "  (no parquet files found)"
    echo ""
    echo "Statistics:"
    cat gnomad_missense_validated/statistics.json 2>/dev/null || echo "  (statistics not found)"
else
    echo "❌ Processing failed. Check logs/process_missense_${SLURM_JOB_ID}.err for details."
fi

echo "============================================================================"

exit $EXIT_CODE