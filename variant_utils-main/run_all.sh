#!/bin/bash
#SBATCH --job-name=gnomad_chrom_batch
#SBATCH --partition=RM
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=200GB
#SBATCH --time=12:00:00
#SBATCH --output=./logs/fetch_missense_%j.out
#SBATCH --error=./logs/fetch_missense_%j.err

# ============================================================================
# SLURM Job: Fetch All gnomAD Missense Variants (Chromosome-Batched)
# ============================================================================
# 
# OPTIMIZATIONS:
# - Process ONE chromosome at a time (load VCF once per chromosome)
# - Within each chromosome, process multiple genes in parallel
# - Aggressive caching to avoid redundant I/O
# - Memory-efficient: 32GB for chromosome VCF + gene processing
#
# Resource justification:
# - 8 CPUs: 4 gene workers + GATK overhead
# - 32GB RAM: Chromosome VCF (~5-10GB) + gene processing (4x 2GB = 8GB) + buffer
# - 48 hours: Conservative for ~19,000 genes with checkpointing
# - RM-shared: Good balance of cost and memory
#
# Memory breakdown per chromosome:
#   - Chromosome VCF extraction: ~5-10GB (chr1 largest)
#   - Gene processing (4 parallel): 4 x 2GB = 8GB
#   - Python overhead: ~4GB
#   - Total: ~20GB typical, 32GB safe
#
# ============================================================================

source ~/miniconda3/etc/profile.d/conda.sh
conda activate protein

echo "============================================================================"
echo "JOB INFORMATION"
echo "============================================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Job Name: $SLURM_JOB_NAME"
echo "Node: $SLURM_NODELIST"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: 32GB"
echo "Start Time: $(date)"
echo "Strategy: Chromosome-batched processing"
echo "============================================================================"

# Load modules
export JAVA_HOME=/jet/home/barazand/NEWOCEAN/java/jdk-21.0.9
export PATH=$JAVA_HOME/bin:$PATH

echo ""
echo "Checking Java installation..."
java -version
echo ""

# Set thread limits for GATK
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export JAVA_OPTS="-Xmx28g -XX:+UseParallelGC -XX:ParallelGCThreads=2"

echo "Environment configured:"
echo "  OMP_NUM_THREADS: $OMP_NUM_THREADS"
echo "  JAVA_OPTS: $JAVA_OPTS"
echo ""

# Navigate to working directory
cd $SLURM_SUBMIT_DIR
echo "Working directory: $(pwd)"
echo ""

# Create logs directory if it doesn't exist
mkdir -p logs

# ============================================================================
# STEP 1: FETCH GNOMAD DATA (CHROMOSOME-BATCHED)
# ============================================================================

echo "============================================================================"
echo "STEP 1: FETCHING GNOMAD DATA (CHROMOSOME-BATCHED)"
echo "============================================================================"
echo ""
echo "Strategy:"
echo "  - Process ONE chromosome at a time (reduces memory)"
echo "  - Extract chromosome VCF once, use for ALL genes on that chromosome"
echo "  - Process 4 genes in parallel per chromosome"
echo "  - Full checkpointing and resume capability"
echo ""
echo "Start time: $(date)"
echo ""

# Run with optimized settings
python /jet/home/barazand/NEWOCEAN/proteins/gnomAD/variant_utils-main/download_all_gnomad_parallel.py \
    --chrom-workers 2 \
    --gene-workers 32

STEP1_EXIT=$?

echo ""
echo "Step 1 completed with exit code: $STEP1_EXIT"
echo "End time: $(date)"
echo ""

# Check if Step 1 succeeded
if [ $STEP1_EXIT -ne 0 ]; then
    echo "ERROR: Step 1 failed. Check logs for details."
    echo "Progress is saved - you can resume by re-running this job"
    exit 1
fi

# ============================================================================
# STEP 2: VALIDATE AND ENRICH ALL VARIANTS
# ============================================================================

echo "============================================================================"
echo "STEP 2: VALIDATING AND ENRICHING VARIANTS"
echo "============================================================================"
echo "Start time: $(date)"
echo ""

python /jet/home/barazand/NEWOCEAN/proteins/gnomAD/variant_utils-main/validate_all_missenses.py

STEP2_EXIT=$?

echo ""
echo "Step 2 completed with exit code: $STEP2_EXIT"
echo "End time: $(date)"
echo ""

if [ $STEP2_EXIT -ne 0 ]; then
    echo "ERROR: Step 2 failed. Check logs for details."
    exit 1
fi

# ============================================================================
# SUMMARY AND CLEANUP
# ============================================================================

echo "============================================================================"
echo "PIPELINE SUMMARY"
echo "============================================================================"

# Display statistics
if [ -f "gnomad_all_genes_validated/statistics.json" ]; then
    echo ""
    echo "Validation Statistics:"
    cat gnomad_all_genes_validated/statistics.json
    echo ""
fi

# Show output files
echo "Output files:"
echo ""
echo "Raw gnomAD results:"
ls -lh gnomad_all_genes/all_gnomad_missense_variants.* 2>/dev/null || echo "  (not found)"
echo ""

echo "Chromosome cache:"
du -sh gnomad_all_genes/chromosome_cache 2>/dev/null || echo "  (not found)"
echo ""

echo "Validated datasets:"
ls -lh gnomad_all_genes_validated/*.parquet 2>/dev/null || echo "  (not found)"
echo ""

# Display disk usage
echo "Total disk usage:"
du -sh gnomad_all_genes 2>/dev/null
du -sh gnomad_all_genes_validated 2>/dev/null
echo ""

# Check for failed genes
if [ -f "gnomad_all_genes/failed_genes.txt" ]; then
    FAILED_COUNT=$(wc -l < gnomad_all_genes/failed_genes.txt)
    echo "WARNING: $FAILED_COUNT genes failed. See gnomad_all_genes/failed_genes.txt"
    echo ""
fi

# Chromosome cache info
if [ -d "gnomad_all_genes/chromosome_cache" ]; then
    echo "Chromosome cache contains:"
    ls -lh gnomad_all_genes/chromosome_cache/*.vcf.gz 2>/dev/null | wc -l | xargs echo "  VCF files:"
    echo ""
fi

echo "============================================================================"
echo "JOB COMPLETED SUCCESSFULLY"
echo "============================================================================"
echo "End Time: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo ""
echo "NOTE: Chromosome VCF cache is preserved for future runs."
echo "      To clear cache: rm -rf gnomad_all_genes/chromosome_cache"
echo "============================================================================"

exit 0