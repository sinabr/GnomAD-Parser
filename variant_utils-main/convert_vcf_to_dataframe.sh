#!/bin/bash
#SBATCH --job-name=gnomad_direct_convert
#SBATCH --partition=model4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100GB                 # Enough for chr1 (~10GB) + processing
#SBATCH --time=8:00:00              # ~4-6 hours for all chromosomes
#SBATCH --output=./logs/direct_convert_%j.out
#SBATCH --error=./logs/direct_convert_%j.err

# ============================================================================
# SLURM Job: Direct VCF to DataFrame Conversion
# ============================================================================
# 
# Streamlined approach:
# - Read source VCFs directly (no bcftools filtering)
# - Filter in Python (SNPs + PASS)
# - Parse VEP annotations
# - Save as Parquet DataFrames
# - Optionally save VCF cache too (for original pipeline)
#
# Everything happens in one Python pass - simpler and cleaner!
# ============================================================================

source ~/miniconda3/etc/profile.d/conda.sh
conda activate protein

# Create output directories
mkdir -p logs
mkdir -p gnomad_all_genes/chromosome_dataframes

echo "============================================================================"
echo "JOB INFORMATION"
echo "============================================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: $SLURM_MEM_PER_NODE"
echo "Start Time: $(date)"
echo "============================================================================"
echo ""

# Check Python and pysam
echo "Checking environment..."
python --version
python -c "import pysam; print(f'pysam version: {pysam.__version__}')"
python -c "import pandas; print(f'pandas version: {pandas.__version__}')"
echo ""

# ============================================================================
# RUN CONVERSION
# ============================================================================

echo "============================================================================"
echo "DIRECT VCF TO DATAFRAME CONVERSION"
echo "============================================================================"
echo ""
echo "Strategy:"
echo "  - Read source VCFs directly (skip bcftools)"
echo "  - Filter in Python (SNPs + PASS)"
echo "  - Parse VEP annotations"
echo "  - Save as Parquet DataFrames"
echo ""
echo "Start time: $(date)"
echo ""

# Run conversion
# Add --save-vcf-cache if you also want VCF cache for original pipeline
python convert_vcf_to_dataframe.py

EXIT_CODE=$?

echo ""
echo "Conversion completed with exit code: $EXIT_CODE"
echo "End time: $(date)"
echo ""

if [ $EXIT_CODE -ne 0 ]; then
    echo "ERROR: Conversion failed. Check logs for details."
    exit 1
fi

# ============================================================================
# SUMMARY
# ============================================================================

echo "============================================================================"
echo "JOB SUMMARY"
echo "============================================================================"
echo ""

echo "DataFrames created:"
du -sh gnomad_all_genes/chromosome_dataframes
ls -lh gnomad_all_genes/chromosome_dataframes/*.parquet 2>/dev/null | wc -l | xargs echo "  Files:"
echo ""

echo "Sample files:"
ls -lh gnomad_all_genes/chromosome_dataframes/*.parquet 2>/dev/null | head -5
echo ""

if [ -f "gnomad_all_genes/chromosome_dataframes/conversion_summary.csv" ]; then
    echo "Conversion summary:"
    cat gnomad_all_genes/chromosome_dataframes/conversion_summary.csv
    echo ""
fi

echo "============================================================================"
echo "JOB COMPLETED SUCCESSFULLY"
echo "============================================================================"
echo "End Time: $(date)"
echo ""
echo "Next steps:"
echo "  1. Load DataFrames: pd.read_parquet('gnomad_all_genes/chromosome_dataframes/chr1_variants.parquet')"
echo "  2. Query genes from memory (fast!)"
echo "============================================================================"

exit 0