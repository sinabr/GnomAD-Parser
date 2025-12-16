#!/bin/bash
# Install bcftools via conda for optimized VCF extraction

echo "=================================="
echo "Installing bcftools via conda"
echo "=================================="

# Activate your conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate protein

# Install bcftools (includes htslib)
echo "Installing bcftools..."
conda install -y -c bioconda bcftools

# Verify installation
echo ""
echo "Verifying installation..."
bcftools --version

echo ""
echo "✅ bcftools installed successfully!"
echo ""
echo "Path: $(which bcftools)"
echo ""
