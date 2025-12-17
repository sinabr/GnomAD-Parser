#!/usr/bin/env python3
"""
Validate Setup - Check all paths and data before running full pipeline

This script checks:
1. All reference data files exist
2. Chromosome Parquet files exist
3. Can load a sample of data
4. VEP fields are present

Run this before submitting the SLURM job to catch any issues early.
"""

import sys
from pathlib import Path
import pandas as pd

# Colors for output
GREEN = '\033[92m'
RED = '\033[91m'
YELLOW = '\033[93m'
RESET = '\033[0m'

def check_file(path: Path, description: str) -> bool:
    """Check if a file exists."""
    if path.exists():
        size = path.stat().st_size / (1024**3)  # GB
        print(f"{GREEN}✓{RESET} {description}: {path} ({size:.2f} GB)")
        return True
    else:
        print(f"{RED}✗{RESET} {description}: {path} NOT FOUND")
        return False

def check_directory(path: Path, description: str, pattern: str = "*") -> tuple:
    """Check if a directory exists and count files."""
    if not path.exists():
        print(f"{RED}✗{RESET} {description}: {path} NOT FOUND")
        return False, 0
    
    files = list(path.glob(pattern))
    print(f"{GREEN}✓{RESET} {description}: {path} ({len(files)} files)")
    return True, len(files)

def main():
    print("="*80)
    print("GNOMAD MISSENSE PROCESSOR - SETUP VALIDATION")
    print("="*80)
    print()
    
    all_ok = True
    
    # ========================================================================
    # 1. CHECK REFERENCE DATA
    # ========================================================================
    print("1. REFERENCE DATA FILES")
    print("-" * 80)
    
    ref_dir = Path("/projects/lugoteam/protein_graphs/GnomAD-Parser/ref_data")
    
    reference_files = {
        "Ensembl GTF": ref_dir / "ensembl" / "Homo_sapiens.GRCh38.112.gtf.gz",
        "Ensembl Peptides": ref_dir / "ensembl" / "Homo_sapiens.GRCh38.pep.all.fa.gz",
        "MANE Summary": ref_dir / "mane" / "MANE.GRCh38.v1.3.summary.txt.gz",
        "UniProt FASTA": ref_dir / "uniprot" / "uniprot_sprot.fasta.gz",
        "UniProt ID Mapping": ref_dir / "uniprot" / "HUMAN_9606_idmapping.dat.gz",
    }
    
    for desc, path in reference_files.items():
        if not check_file(path, desc):
            all_ok = False
    
    print()
    
    # ========================================================================
    # 2. CHECK CHROMOSOME PARQUET FILES
    # ========================================================================
    print("2. CHROMOSOME PARQUET FILES")
    print("-" * 80)
    
    input_dir = Path("gnomad_all_genes/chromosome_dataframes")
    exists, count = check_directory(input_dir, "Chromosome dataframes", "chr*.parquet")
    
    if not exists:
        all_ok = False
    elif count == 0:
        print(f"{RED}✗{RESET} No chromosome Parquet files found!")
        all_ok = False
    else:
        # List all chromosome files
        expected_chroms = [str(i) for i in range(1, 23)] + ['X', 'Y']
        found_chroms = []
        missing_chroms = []
        
        for chrom in expected_chroms:
            chrom_file = input_dir / f"chr{chrom}_variants.parquet"
            if chrom_file.exists():
                found_chroms.append(chrom)
            else:
                missing_chroms.append(chrom)
        
        print(f"   Found: {len(found_chroms)}/24 chromosomes")
        if missing_chroms:
            print(f"{YELLOW}!{RESET} Missing chromosomes: {', '.join(missing_chroms)}")
    
    print()
    
    # ========================================================================
    # 3. TEST LOADING A SAMPLE
    # ========================================================================
    print("3. TEST DATA LOADING")
    print("-" * 80)
    
    # Try to load chromosome 22 (smallest)
    test_file = input_dir / "chr22_variants.parquet"
    
    if test_file.exists():
        try:
            print(f"Loading {test_file.name} as test...")
            df = pd.read_parquet(test_file)
            
            print(f"{GREEN}✓{RESET} Successfully loaded {len(df):,} rows")
            print(f"   Columns: {len(df.columns)}")
            print()
            
            # Check for required fields
            print("   Checking required fields:")
            required_fields = ['CHROM', 'POS', 'REF', 'ALT', 'source', 'Consequence', 'Feature', 'HGVSp', 'SYMBOL']
            
            missing_fields = []
            for field in required_fields:
                if field in df.columns:
                    non_null = df[field].notna().sum()
                    print(f"   {GREEN}✓{RESET} {field}: {non_null:,}/{len(df):,} non-null ({100*non_null/len(df):.1f}%)")
                else:
                    print(f"   {RED}✗{RESET} {field}: MISSING")
                    missing_fields.append(field)
                    all_ok = False
            
            if not missing_fields:
                # Check for missense variants
                print()
                print("   Checking for missense variants:")
                missense_count = df['Consequence'].str.contains('missense_variant', case=False, na=False).sum()
                print(f"   Found {missense_count:,} missense variants ({100*missense_count/len(df):.1f}% of total)")
                
                if missense_count == 0:
                    print(f"{YELLOW}!{RESET} WARNING: No missense variants found in chr22")
                    print("   This might be normal if chr22 has few missense variants")
            
        except Exception as e:
            print(f"{RED}✗{RESET} Failed to load test file: {e}")
            all_ok = False
    else:
        print(f"{YELLOW}!{RESET} chr22 not found, skipping sample test")
    
    print()
    
    # ========================================================================
    # 4. CHECK OUTPUT DIRECTORY
    # ========================================================================
    print("4. OUTPUT DIRECTORY")
    print("-" * 80)
    
    output_dir = Path("gnomad_missense_validated")
    if output_dir.exists():
        print(f"{GREEN}✓{RESET} Output directory exists: {output_dir}")
        existing_files = list(output_dir.glob("*.parquet"))
        if existing_files:
            print(f"{YELLOW}!{RESET} WARNING: Output directory contains {len(existing_files)} existing files")
            print("   These will be overwritten when the pipeline runs")
    else:
        print(f"{GREEN}✓{RESET} Output directory will be created: {output_dir}")
    
    print()
    
    # ========================================================================
    # SUMMARY
    # ========================================================================
    print("="*80)
    if all_ok:
        print(f"{GREEN}✓ ALL CHECKS PASSED{RESET}")
        print()
        print("You're ready to run the pipeline:")
        print("  sbatch run_process_missense.sh")
    else:
        print(f"{RED}✗ SOME CHECKS FAILED{RESET}")
        print()
        print("Please fix the issues above before running the pipeline.")
        return 1
    print("="*80)
    
    return 0

if __name__ == "__main__":
    sys.exit(main())