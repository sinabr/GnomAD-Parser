#!/bin/bash
#SBATCH --job-name=gnomad_chrom_batch
#SBATCH --partition=model4          # 64 CPUs, ~192GB RAM, fast scratch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32          # leave headroom below 64 for GC/IO
#SBATCH --mem=170GB                 # below the ~192GB node cap for safety
#SBATCH --time=24:00:00
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

# Stage to fast scratch
SOURCE_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
SCRATCH_ROOT="/scratch/barazand"
WORKDIR="${SCRATCH_ROOT}/gnomad_run_${SLURM_JOB_ID:-manual}"
mkdir -p "$WORKDIR"
echo "Staging to scratch: $WORKDIR"
rsync -a --delete "$SOURCE_DIR"/ "$WORKDIR"/
cd "$WORKDIR"
# Determine project directory (works whether submitted from repo root or within subdir)
if [ -d "$WORKDIR/variant_utils-main" ]; then
    PROJ_DIR="$WORKDIR/variant_utils-main"
else
    PROJ_DIR="$WORKDIR"
fi

# Optionally rewrite external_tools.json to point to staged copies (default: skip).
# Set REWRITE_EXTERNAL=1 only if gatk/picard/gnomAD/spliceAI are present under the staged tree.
REWRITE_EXTERNAL="${REWRITE_EXTERNAL:-0}"
if [ "$REWRITE_EXTERNAL" -eq 1 ]; then
python - <<'PY'
import json, os, pathlib
proj_dir = pathlib.Path(os.environ["PROJ_DIR"])
cfg_path = proj_dir / "external_tools.json"
with cfg_path.open() as f:
    cfg = json.load(f)
java_path = pathlib.Path(os.path.expanduser(cfg.get("java", "")))
base = proj_dir
v4_root = base / "data"
v2_root = base / "data"
splice_root = base / "data" / "dbs" / "spliceAI"
if not v4_root.exists() or not v2_root.exists():
    print("Scratch VCF roots missing; keeping external_tools.json unchanged.")
    raise SystemExit(0)
new_cfg = {
    "java": str(java_path),
    "gatk": str(base / "gatk-4.6.1.0" / "gatk"),
    "picard_filepath": str(base / "picard.jar"),
    "gnomad_v4_vcf_root": str(v4_root),
    "gnomad_v2_vcf_root": str(v2_root),
    "spliceAIRoot": str(splice_root)
}
cfg_path.write_text(json.dumps(new_cfg, indent=4))
print("external_tools.json rewritten for scratch staging:")
print(json.dumps(new_cfg, indent=2))
PY
else
    echo "Using external_tools.json paths as-is (cluster absolute paths expected)."
fi

# Allow overriding parallelism and temp location without editing the script
# OPTIMIZED: Increased gene workers to fully utilize CPU (ProcessPoolExecutor has no GIL)
CHROM_WORKERS="${CHROM_WORKERS:-1}"   # chromosome-level parallelism (memory-heavy)
GENE_WORKERS="${GENE_WORKERS:-32}"    # gene-level parallelism (CPU/I/O-bound); increased for ProcessPool
export TMPDIR="${TMPDIR:-$PROJ_DIR/gnomad_all_genes/tmp}"

echo "============================================================================"
echo "JOB INFORMATION"
echo "============================================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Job Name: $SLURM_JOB_NAME"
echo "Node: $SLURM_NODELIST"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: $SLURM_MEM_PER_NODE"
echo "Start Time: $(date)"
echo "Strategy: Chromosome-batched processing"
echo "Chromosome workers: $CHROM_WORKERS"
echo "Gene workers: $GENE_WORKERS"
echo "TMPDIR: $TMPDIR"
echo "============================================================================"

# Load modules
export JAVA_HOME=$HOME/software/jdk-21
export PATH=$JAVA_HOME/bin:$PATH

echo ""
echo "Checking Java installation..."
java -version
echo ""

# OPTIMIZATION: bcftools replaces GATK for chromosome extraction (10-50x faster)
# GATK/Java only needed if you fall back to old method
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export JAVA_OPTS="-Xmx96g -XX:+UseParallelGC -XX:ParallelGCThreads=4"

# Note: With bcftools + pysam optimization:
# - bcftools: Fast chromosome extraction (minutes instead of hours)
# - pysam: Fast gene-level queries (milliseconds instead of seconds)

echo "Environment configured:"
echo "  OMP_NUM_THREADS: $OMP_NUM_THREADS"
echo "  JAVA_OPTS: $JAVA_OPTS"
echo ""

# Create logs directory if it doesn't exist
mkdir -p "$PROJ_DIR/logs"
mkdir -p "$TMPDIR"

# ============================================================================
# STEP 1: FETCH GNOMAD DATA (CHROMOSOME-BATCHED)
# ============================================================================

echo "============================================================================"
echo "STEP 1: FETCHING GNOMAD DATA (CHROMOSOME-BATCHED)"
echo "============================================================================"
echo ""
echo "Strategy (OPTIMIZED):"
echo "  - Process ONE chromosome at a time (reduces memory)"
echo "  - Extract chromosome VCF once with bcftools (10-50x faster than GATK, cached)"
echo "  - Use pysam for fast gene-level extraction (100-500x faster than GATK)"
echo "  - ProcessPoolExecutor for true parallel CPU execution (no GIL)"
echo "  - Process multiple genes in parallel per chromosome (${GENE_WORKERS} workers)"
echo "  - Full checkpointing and resume capability"
echo ""
echo "Start time: $(date)"
echo ""

# Run with optimized settings (scratch copy)
(cd "$PROJ_DIR" && python download_all_gnomad_parallel.py \
    --chrom-workers "$CHROM_WORKERS" \
    --gene-workers "$GENE_WORKERS")

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

(cd "$PROJ_DIR" && python validate_all_missenses.py)

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

# Sync results back to submit directory
if [ -n "$SLURM_SUBMIT_DIR" ]; then
    echo "Syncing results back to $SLURM_SUBMIT_DIR"
    rsync -a "$PROJ_DIR/gnomad_all_genes"/ "$SLURM_SUBMIT_DIR/gnomad_all_genes"/ 2>/dev/null
    rsync -a "$PROJ_DIR/gnomad_all_genes_validated"/ "$SLURM_SUBMIT_DIR/gnomad_all_genes_validated"/ 2>/dev/null
    rsync -a "$PROJ_DIR/logs"/ "$SLURM_SUBMIT_DIR/logs"/ 2>/dev/null
fi

exit 0
