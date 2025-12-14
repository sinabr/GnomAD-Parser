#!/usr/bin/env python3

"""
Download gnomAD files and genomic reference data
Usage: python download_gnomad.py [--dry-run] [--chromosomes 1,2,X]
"""

import os
import sys
import argparse
import subprocess
from pathlib import Path
from typing import List, Tuple

class FileDownloader:
    def __init__(self, dry_run=False):
        self.dry_run = dry_run
        self.failed_downloads = []
        
    def download_file(self, url: str, output_path: str) -> bool:
        """Download a file using wget with progress bar"""
        if self.dry_run:
            print(f"[DRY RUN] Would download: {url}")
            return True
            
        # Check if file already exists
        if os.path.exists(output_path):
            file_size = os.path.getsize(output_path)
            if file_size > 0:
                print(f"File exists ({file_size:,} bytes): {output_path}")
                return True
        
        print(f"Downloading: {url}")
        try:
            cmd = [
                "wget",
                "-c",  # Continue partial downloads
                "--progress=bar:force",
                "-O", output_path,
                url
            ]
            result = subprocess.run(cmd, check=True)
            return True
        except subprocess.CalledProcessError as e:
            print(f"ERROR: Failed to download {url}")
            self.failed_downloads.append(url)
            return False
        except FileNotFoundError:
            print("ERROR: wget not found. Please install wget.")
            sys.exit(1)

def setup_directories():
    """Create directory structure"""
    Path("exomes").mkdir(exist_ok=True)
    Path("genomes").mkdir(exist_ok=True)
    print("✓ Directory structure created")

def download_gnomad_files(downloader: FileDownloader, chromosomes: List[str]):
    """Download gnomAD exomes and genomes VCF files"""
    
    base_url = "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf"
    
    # Download exomes
    print("\n" + "="*60)
    print("DOWNLOADING GNOMAD EXOMES v4.1")
    print("="*60)
    
    for chr_num in chromosomes:
        print(f"\n--- Chromosome {chr_num} (Exomes) ---")
        
        # VCF file
        vcf_file = f"gnomad.exomes.v4.1.sites.chr{chr_num}.vcf.bgz"
        vcf_url = f"{base_url}/exomes/{vcf_file}"
        downloader.download_file(vcf_url, f"exomes/{vcf_file}")
        
        # Index file
        tbi_file = f"{vcf_file}.tbi"
        tbi_url = f"{vcf_url}.tbi"
        downloader.download_file(tbi_url, f"exomes/{tbi_file}")
    
    # Download genomes
    print("\n" + "="*60)
    print("DOWNLOADING GNOMAD GENOMES v4.1")
    print("="*60)
    
    for chr_num in chromosomes:
        print(f"\n--- Chromosome {chr_num} (Genomes) ---")
        
        # VCF file
        vcf_file = f"gnomad.genomes.v4.1.sites.chr{chr_num}.vcf.bgz"
        vcf_url = f"{base_url}/genomes/{vcf_file}"
        downloader.download_file(vcf_url, f"genomes/{vcf_file}")
        
        # Index file
        tbi_file = f"{vcf_file}.tbi"
        tbi_url = f"{vcf_url}.tbi"
        downloader.download_file(tbi_url, f"genomes/{tbi_file}")

def download_reference_files(downloader: FileDownloader):
    """Download reference genome and annotation files"""
    
    print("\n" + "="*60)
    print("DOWNLOADING REFERENCE FILES")
    print("="*60)
    
    files_to_download = [
        # RefSeq gene files
        (
            "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/refseqgene.1.genomic.fna.gz",
            "refseqgene.1.genomic.fna.gz"
        ),
        (
            "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/refseqgene.1.genomic.gbff.gz",
            "refseqgene.1.genomic.gbff.gz"
        ),
        # GRCh38 annotation files
        (
            "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/GRCh38_latest_genomic.gff.gz",
            "GRCh38_latest_genomic.gff.gz"
        ),
        (
            "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/GRCh38_latest_rna.fna.gz",
            "GRCh38_latest_rna.fna.gz"
        ),
        # UniProt Swiss-Prot
        (
            "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz",
            "uniprot_sprot.fasta.gz"
        ),
    ]
    
    for url, filename in files_to_download:
        print(f"\n--- {filename} ---")
        downloader.download_file(url, filename)
    
    # Optionally decompress GFF file
    if os.path.exists("GRCh38_latest_genomic.gff.gz") and not downloader.dry_run:
        print("\nDecompressing GRCh38_latest_genomic.gff.gz...")
        try:
            subprocess.run(["gunzip", "-f", "GRCh38_latest_genomic.gff.gz"], check=True)
            print("✓ Decompressed successfully")
        except subprocess.CalledProcessError:
            print("Warning: Failed to decompress GFF file")

def print_summary(downloader: FileDownloader):
    """Print download summary"""
    
    print("\n" + "="*60)
    print("DOWNLOAD SUMMARY")
    print("="*60)
    
    if not downloader.dry_run:
        exomes_count = len(list(Path("exomes").glob("*")))
        genomes_count = len(list(Path("genomes").glob("*")))
        ref_count = len([f for f in os.listdir(".") if f.endswith(('.gz', '.gff'))])
        
        print(f"Exomes directory:  {exomes_count:3d} files")
        print(f"Genomes directory: {genomes_count:3d} files")
        print(f"Reference files:   {ref_count:3d} files")
    
    if downloader.failed_downloads:
        print(f"\n⚠ {len(downloader.failed_downloads)} downloads failed:")
        for url in downloader.failed_downloads:
            print(f"  - {url}")
    else:
        print("\n✓ All downloads completed successfully!")

def main():
    parser = argparse.ArgumentParser(
        description="Download gnomAD and reference genomic data"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be downloaded without actually downloading"
    )
    parser.add_argument(
        "--chromosomes",
        type=str,
        default="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y",
        help="Comma-separated list of chromosomes to download (default: all)"
    )
    parser.add_argument(
        "--skip-gnomad",
        action="store_true",
        help="Skip gnomAD files, only download reference files"
    )
    parser.add_argument(
        "--skip-reference",
        action="store_true",
        help="Skip reference files, only download gnomAD"
    )
    
    args = parser.parse_args()
    
    # Parse chromosomes
    chromosomes = [c.strip() for c in args.chromosomes.split(",")]
    
    print("gnomAD and Reference Data Downloader")
    print(f"Chromosomes: {', '.join(chromosomes)}")
    print(f"Dry run: {args.dry_run}")
    print("")
    
    # Setup
    setup_directories()
    downloader = FileDownloader(dry_run=args.dry_run)
    
    # Download files
    if not args.skip_gnomad:
        download_gnomad_files(downloader, chromosomes)
    
    if not args.skip_reference:
        download_reference_files(downloader)
    
    # Summary
    print_summary(downloader)

if __name__ == "__main__":
    main()