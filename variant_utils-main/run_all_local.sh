#!/bin/bash
# Lightweight local runner for quick testing on an interactive node (no SLURM).
# Assumes only 2 CPU cores are available.

source ~/miniconda3/etc/profile.d/conda.sh
conda activate protein

set -euo pipefail

# Activate environment (adjust if your conda is elsewhere). Temporarily disable
# nounset to avoid issues in conda activate/deactivate scripts that expect
# unset backup vars.
if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
    set +u
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
    conda activate protein
    set -u
fi

# Resolve project directory (this script lives in variant_utils-main/)
PROJ_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PROJ_DIR"

# Parallelism tuned for 2 cores
CHROM_WORKERS="${CHROM_WORKERS:-1}"
GENE_WORKERS="${GENE_WORKERS:-2}"
# Optional single-chromosome smoke test
SIMPLE_CHROM_TEST="${SIMPLE_CHROM_TEST:-0}"
TEST_CHROM="${TEST_CHROM:-4}"
MAX_GENES="${MAX_GENES:-5}"

# Temp/log locations
export TMPDIR="${TMPDIR:-$PROJ_DIR/gnomad_all_genes/tmp_local}"
mkdir -p "$TMPDIR" "$PROJ_DIR/logs"

echo "============================================================"
echo "Local gnomAD run (interactive, no SLURM)"
echo "Project dir: $PROJ_DIR"
echo "Chromosome workers: $CHROM_WORKERS"
echo "Gene workers: $GENE_WORKERS"
echo "Simple chrom test: $SIMPLE_CHROM_TEST (chrom=$TEST_CHROM, max_genes=$MAX_GENES)"
echo "TMPDIR: $TMPDIR"
echo "Start: $(date)"
echo "============================================================"

# Step 1: Fetch gnomAD data
if [ "$SIMPLE_CHROM_TEST" -eq 1 ]; then
    # Lightweight single-chromosome test to verify bcftools + parsing
    python download_one_chrom_test.py \
        --chrom "$TEST_CHROM" \
        --gene-workers "$GENE_WORKERS" \
        --max-genes "$MAX_GENES"
else
    # Full multi-chromosome run
    python download_all_gnomad_parallel.py \
        --chrom-workers "$CHROM_WORKERS" \
        --gene-workers "$GENE_WORKERS"
fi

# Step 2: Validate/enrich all variants
python validate_all_missenses.py

echo "============================================================"
echo "Local run complete at $(date)"
echo "============================================================"
