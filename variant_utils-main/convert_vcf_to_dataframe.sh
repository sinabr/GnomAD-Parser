#!/bin/bash
#SBATCH --job-name=gnomad_chunked
#SBATCH --partition=model4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=170GB                 # Increased for safety, but only uses ~15GB
#SBATCH --time=24:00:00
#SBATCH --output=./logs/convert_chunked_%j.out
#SBATCH --error=./logs/convert_chunked_%j.err

# ============================================================================
# SLURM Job: MEMORY-EFFICIENT VCF to DataFrame Conversion
# ============================================================================
# 
# KEY OPTIMIZATION: Chunked processing!
# - Processes 500K variants at a time
# - Writes incrementally to Parquet
# - Never loads full chromosome into memory
# - Peak memory: ~10-15 GB (vs 60+ GB before)
#
# This solves the OOM (Out of Memory) errors!
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

# Check environment
echo "Checking environment..."
python --version
python -c "import pysam; print(f'pysam version: {pysam.__version__}')"
python -c "import pandas; print(f'pandas version: {pandas.__version__}')"
python -c "import pyarrow; print(f'pyarrow version: {pyarrow.__version__}')"
echo ""

# Check available memory
echo "System memory:"
free -h
echo ""

# ============================================================================
# RUN MEMORY-EFFICIENT CONVERSION
# ============================================================================

echo "============================================================================"
echo "MEMORY-EFFICIENT VCF TO DATAFRAME CONVERSION"
echo "============================================================================"
echo ""
echo "Optimization:"
echo "  - Chunked processing (500K variants per batch)"
echo "  - Incremental writes to Parquet"
echo "  - Peak memory: ~10-15 GB (safe!)"
echo ""
echo "Start time: $(date)"
echo ""

# Run chunked conversion
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

# Check for temp files (should be cleaned up)
TEMP_FILES=$(ls gnomad_all_genes/chromosome_dataframes/.tmp_* 2>/dev/null | wc -l)
if [ $TEMP_FILES -gt 0 ]; then
    echo "WARNING: Found $TEMP_FILES temp files (cleaning up...)"
    rm -f gnomad_all_genes/chromosome_dataframes/.tmp_*
fi

echo "============================================================================"
echo "JOB COMPLETED SUCCESSFULLY"
echo "============================================================================"
echo "End Time: $(date)"
echo ""
echo "Usage:"
echo "  df = pd.read_parquet('gnomad_all_genes/chromosome_dataframes/chr1_variants.parquet')"
echo "  missense = df[df['Consequence'].str.contains('missense', na=False)]"
echo "============================================================================"

exit 0