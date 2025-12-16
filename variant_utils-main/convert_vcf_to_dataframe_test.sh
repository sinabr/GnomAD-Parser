#!/bin/bash
#SBATCH --job-name=gnomad_test_chr22
#SBATCH --partition=model4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32GB                  # Chr22 is small, 32GB should be plenty
#SBATCH --time=1:00:00              # 1 hour should be enough for chr22
#SBATCH --output=./logs/test_convert_%j.out
#SBATCH --error=./logs/test_convert_%j.err

# ============================================================================
# TEST SLURM Job: Direct VCF to DataFrame Conversion (Chr22 Only)
# ============================================================================
#
# Quick test on chromosome 22 (smallest autosomal chromosome) to validate:
# - VCF reading works
# - Filtering logic (SNPs + PASS + MISSENSE) is correct
# - VEP parsing works
# - Output format is correct
#
# Usage:
#   sbatch convert_vcf_to_dataframe_test.sh
#
# Optional - for very quick test (first 1000 variants per VCF):
#   Edit the python command below to add: --max-variants 1000
# ============================================================================

source ~/miniconda3/etc/profile.d/conda.sh
conda activate protein

# Create output directories
mkdir -p logs
mkdir -p gnomad_test/chromosome_dataframes

echo "============================================================================"
echo "TEST JOB INFORMATION"
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
# RUN TEST CONVERSION
# ============================================================================

echo "============================================================================"
echo "TEST: DIRECT VCF TO DATAFRAME CONVERSION (CHR22 - MISSENSE ONLY)"
echo "============================================================================"
echo ""
echo "Strategy:"
echo "  - Read source VCFs directly (skip bcftools)"
echo "  - Filter: SNPs + PASS + MISSENSE variants only"
echo "  - Parse VEP annotations"
echo "  - Save as Parquet DataFrames"
echo ""
echo "Test chromosome: 22 (smallest autosomal chromosome)"
echo ""
echo "Start time: $(date)"
echo ""

# Run test conversion
# Add --max-variants 1000 for a very quick test (processes only first 1000 variants per VCF)
python convert_vcf_to_dataframe_test.py

EXIT_CODE=$?

echo ""
echo "Test conversion completed with exit code: $EXIT_CODE"
echo "End time: $(date)"
echo ""

if [ $EXIT_CODE -ne 0 ]; then
    echo "ERROR: Test conversion failed. Check logs for details."
    exit 1
fi

# ============================================================================
# VERIFICATION
# ============================================================================

echo "============================================================================"
echo "TEST VERIFICATION"
echo "============================================================================"
echo ""

echo "Output directory contents:"
ls -lh gnomad_test/chromosome_dataframes/
echo ""

if [ -f "gnomad_test/chromosome_dataframes/chr22_missense_variants.parquet" ]; then
    echo "✅ Test output file created successfully!"
    echo ""
    echo "File details:"
    ls -lh gnomad_test/chromosome_dataframes/chr22_missense_variants.parquet
    echo ""

    echo "Quick data inspection with Python:"
    python << 'EOF'
import pandas as pd
from pathlib import Path

parquet_file = Path("gnomad_test/chromosome_dataframes/chr22_missense_variants.parquet")
if parquet_file.exists():
    df = pd.read_parquet(parquet_file)
    print(f"Shape: {df.shape}")
    print(f"\nColumns ({len(df.columns)}):")
    for i, col in enumerate(df.columns, 1):
        print(f"  {i}. {col}")

    print(f"\nSample data (first 3 rows):")
    print(df.head(3))

    print(f"\nConsequence field check:")
    if 'Consequence' in df.columns:
        print(df['Consequence'].value_counts().head(10))

        # Verify all are missense
        has_missense = df['Consequence'].str.contains('missense_variant', na=False)
        print(f"\nAll rows contain 'missense_variant': {has_missense.all()}")
        print(f"Rows with missense: {has_missense.sum()} / {len(df)}")

    print(f"\nSource distribution:")
    if 'source' in df.columns:
        print(df['source'].value_counts())
EOF

else
    echo "❌ Test output file NOT created!"
    exit 1
fi

echo ""
echo "============================================================================"
echo "TEST SUMMARY"
echo "============================================================================"
echo ""

if [ -f "gnomad_test/chromosome_dataframes/test_conversion_summary.csv" ]; then
    echo "Test summary:"
    cat gnomad_test/chromosome_dataframes/test_conversion_summary.csv
    echo ""
fi

echo "============================================================================"
echo "TEST COMPLETED SUCCESSFULLY"
echo "============================================================================"
echo "End Time: $(date)"
echo ""
echo "Next steps:"
echo "  1. Review the output above to verify missense filtering worked"
echo "  2. If everything looks good, run full version:"
echo "     sbatch convert_vcf_to_dataframe.sh"
echo "============================================================================"

exit 0
