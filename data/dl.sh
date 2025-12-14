#!/bin/bash
#SBATCH --job-name=download_gnomad
#SBATCH --partition=RM-shared
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --output=logs/download_gnomad_dl-%j.out
#SBATCH --error=logs/download_gnomad-%j.err

# Script to download gnomAD files and genomic reference data
# Usage: sbatch download_gnomad.sh

set -e  # Exit on error

# Create logs directory if it doesn't exist
mkdir -p logs

# Print job information
echo "Job ID: $SLURM_JOB_ID"
echo "Job Name: $SLURM_JOB_NAME"
echo "Node: $SLURM_NODELIST"
echo "Start Time: $(date)"
echo "Working Directory: $(pwd)"
echo ""

# Create directory structure
echo "Creating directory structure..."
mkdir -p exomes
mkdir -p genomes

# Base URLs
GNOMAD_BASE="https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf"
REFSEQ_BASE="https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers"
GRCH38_BASE="https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest"
UNIPROT_BASE="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete"

# Chromosomes to download
CHROMOSOMES=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y)

echo "================================"
echo "Downloading gnomAD Exomes v4.1"
echo "================================"

for chr in "${CHROMOSOMES[@]}"; do
    echo "Downloading chromosome ${chr} exomes..."
    
    # Download VCF file
    wget -c -P exomes \
        "${GNOMAD_BASE}/exomes/gnomad.exomes.v4.1.sites.chr${chr}.vcf.bgz" || {
        echo "Warning: Failed to download chr${chr} exomes VCF"
    }
    
    # Download index file
    wget -c -P exomes \
        "${GNOMAD_BASE}/exomes/gnomad.exomes.v4.1.sites.chr${chr}.vcf.bgz.tbi" || {
        echo "Warning: Failed to download chr${chr} exomes index"
    }
done

echo ""
echo "================================"
echo "Downloading gnomAD Genomes v4.1"
echo "================================"

for chr in "${CHROMOSOMES[@]}"; do
    echo "Downloading chromosome ${chr} genomes..."
    
    # Download VCF file
    wget -c -P genomes \
        "${GNOMAD_BASE}/genomes/gnomad.genomes.v4.1.sites.chr${chr}.vcf.bgz" || {
        echo "Warning: Failed to download chr${chr} genomes VCF"
    }
    
    # Download index file
    wget -c -P genomes \
        "${GNOMAD_BASE}/genomes/gnomad.genomes.v4.1.sites.chr${chr}.vcf.bgz.tbi" || {
        echo "Warning: Failed to download chr${chr} genomes index"
    }
done

echo ""
echo "================================"
echo "Downloading Reference Files"
echo "================================"

# Download RefSeq gene files
echo "Downloading RefSeq gene files..."
wget -c "${REFSEQ_BASE}/refseqgene.1.genomic.fna.gz" || {
    echo "Warning: Failed to download refseqgene.1.genomic.fna.gz"
}

wget -c "${REFSEQ_BASE}/refseqgene.1.genomic.gbff.gz" || {
    echo "Warning: Failed to download refseqgene.1.genomic.gbff.gz"
}

# Download GRCh38 annotation files
echo "Downloading GRCh38 annotation files..."
wget -c "${GRCH38_BASE}/GRCh38_latest_genomic.gff.gz" && \
    gunzip -f GRCh38_latest_genomic.gff.gz || {
    echo "Warning: Failed to download GRCh38_latest_genomic.gff"
}

wget -c "${GRCH38_BASE}/GRCh38_latest_rna.fna.gz" || {
    echo "Warning: Failed to download GRCh38_latest_rna.fna.gz"
}

# Download UniProt Swiss-Prot
echo "Downloading UniProt Swiss-Prot database..."
wget -c "${UNIPROT_BASE}/uniprot_sprot.fasta.gz" || {
    echo "Warning: Failed to download uniprot_sprot.fasta.gz"
}

echo ""
echo "================================"
echo "Download Summary"
echo "================================"
echo "Exomes directory: $(ls -1 exomes | wc -l) files"
echo "Genomes directory: $(ls -1 genomes | wc -l) files"
echo "Reference files in current directory: $(ls -1 *.gz *.gff 2>/dev/null | wc -l) files"
echo ""
echo "Download complete!"
echo "Note: Some files may have failed - check warnings above"
echo ""
echo "End Time: $(date)"
echo "Job completed successfully!"