#!/usr/bin/env python3
"""
OPTIONAL: Preprocess gnomAD VCFs to Parquet for ultra-fast gene queries.

This script is NOT needed for the main pipeline. Use it only if:
- You plan to run gene queries many times
- You want instant gene lookups (seconds instead of minutes)
- You have 50-200GB storage available

One-time preprocessing: 3-6 hours
Subsequent gene queries: 5-15 minutes for all genes, seconds for individual genes

Usage:
    python prep_gnomad_parquet.py --chromosomes 1,2,3,X,Y
    python prep_gnomad_parquet.py --all-chromosomes
"""

import pandas as pd
import argparse
from pathlib import Path
from pysam import VariantFile
import logging
from typing import List
import time

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Configuration
EXTERNAL_TOOLS_CONFIG = "external_tools.json"
OUTPUT_DIR = Path("gnomad_parquet_cache")
OUTPUT_DIR.mkdir(exist_ok=True)
ASSEMBLY = "GRCh38"


def expand_vep_annotations(vep_str: str, vep_columns: List[str]) -> dict:
    """Expand VEP pipe-delimited string into dictionary."""
    if not vep_str or vep_str == '':
        return {}

    # VEP can have multiple annotations separated by comma
    # Take the first one for simplicity (you can modify this)
    first_annotation = vep_str.split(',')[0] if ',' in vep_str else vep_str
    fields = first_annotation.split('|')

    return {col: fields[i] if i < len(fields) else None
            for i, col in enumerate(vep_columns)}


def preprocess_chromosome(
    chrom: str,
    assembly: str,
    external_config: dict
) -> None:
    """
    Convert chromosome VCF to Parquet with expanded VEP annotations.

    This makes subsequent queries 10-50x faster than reading VCF.
    """
    logger.info(f"\n{'='*80}")
    logger.info(f"PREPROCESSING CHROMOSOME {chrom}")
    logger.info(f"{'='*80}")

    start_time = time.time()

    # Setup paths
    release_version = "v4.1" if assembly == "GRCh38" else "r2.1.1"
    chr_prefix = "chr" if release_version == "v4.1" else ""

    gnomad_vcf_root = Path(
        external_config['gnomad_v4_vcf_root'] if release_version == "v4.1"
        else external_config['gnomad_v2_vcf_root']
    )

    exomes_vcf = gnomad_vcf_root / f"exomes/gnomad.exomes.{release_version}.sites.{chr_prefix}{chrom}.vcf.bgz"
    genomes_vcf = gnomad_vcf_root / f"genomes/gnomad.genomes.{release_version}.sites.{chr_prefix}{chrom}.vcf.bgz"

    if not exomes_vcf.exists() or not genomes_vcf.exists():
        logger.warning(f"Skipping chr{chrom}: VCF files not found")
        return

    # Get VEP columns
    vcf = VariantFile(str(exomes_vcf))
    vep_columns = vcf.header.info['vep'].description.split("Format: ")[1].split("|")
    vcf.close()

    logger.info(f"VEP columns: {len(vep_columns)}")

    # Process exomes
    logger.info("Processing exomes...")
    exomes_variants = []
    vcf_exomes = VariantFile(str(exomes_vcf))

    for i, variant in enumerate(vcf_exomes):
        if i > 0 and i % 100000 == 0:
            logger.info(f"  Processed {i:,} exome variants...")

        # Only PASS SNPs
        if variant.filter.keys() != ['PASS']:
            continue
        if len(variant.ref) != 1:
            continue

        for alt_idx, alt in enumerate(variant.alts or []):
            if len(str(alt)) != 1:
                continue

            # Extract basic fields
            ac = variant.info.get('AC', [None])
            af = variant.info.get('AF', [None])

            ac_val = ac[alt_idx] if isinstance(ac, (list, tuple)) and alt_idx < len(ac) else ac
            af_val = af[alt_idx] if isinstance(af, (list, tuple)) and alt_idx < len(af) else af

            # Extract and expand VEP
            vep_data = variant.info.get('vep', [])
            vep_str = ','.join(vep_data) if isinstance(vep_data, (list, tuple)) else str(vep_data)
            vep_dict = expand_vep_annotations(vep_str, vep_columns)

            # Combine all fields
            row = {
                'CHROM': variant.chrom,
                'POS': variant.pos,
                'ID': variant.id or '.',
                'REF': variant.ref,
                'ALT': str(alt),
                'AC': ac_val,
                'AF': af_val,
                'source': 'exomes',
                **vep_dict  # Expand VEP columns
            }
            exomes_variants.append(row)

    vcf_exomes.close()
    logger.info(f"  Collected {len(exomes_variants):,} exome variants")

    # Process genomes
    logger.info("Processing genomes...")
    genomes_variants = []
    vcf_genomes = VariantFile(str(genomes_vcf))

    for i, variant in enumerate(vcf_genomes):
        if i > 0 and i % 100000 == 0:
            logger.info(f"  Processed {i:,} genome variants...")

        if variant.filter.keys() != ['PASS']:
            continue
        if len(variant.ref) != 1:
            continue

        for alt_idx, alt in enumerate(variant.alts or []):
            if len(str(alt)) != 1:
                continue

            ac = variant.info.get('AC', [None])
            af = variant.info.get('AF', [None])

            ac_val = ac[alt_idx] if isinstance(ac, (list, tuple)) and alt_idx < len(ac) else ac
            af_val = af[alt_idx] if isinstance(af, (list, tuple)) and alt_idx < len(af) else af

            vep_data = variant.info.get('vep', [])
            vep_str = ','.join(vep_data) if isinstance(vep_data, (list, tuple)) else str(vep_data)
            vep_dict = expand_vep_annotations(vep_str, vep_columns)

            row = {
                'CHROM': variant.chrom,
                'POS': variant.pos,
                'ID': variant.id or '.',
                'REF': variant.ref,
                'ALT': str(alt),
                'AC': ac_val,
                'AF': af_val,
                'source': 'genomes',
                **vep_dict
            }
            genomes_variants.append(row)

    vcf_genomes.close()
    logger.info(f"  Collected {len(genomes_variants):,} genome variants")

    # Combine and save
    logger.info("Combining and saving to Parquet...")
    all_variants = exomes_variants + genomes_variants
    df = pd.DataFrame(all_variants)

    # Save as Parquet with compression
    output_file = OUTPUT_DIR / f"chr{chrom}.parquet"
    df.to_parquet(
        output_file,
        engine='pyarrow',
        compression='snappy',
        index=False
    )

    elapsed = time.time() - start_time
    file_size_mb = output_file.stat().st_size / 1e6

    logger.info(f"✅ Chromosome {chrom} complete:")
    logger.info(f"   Variants: {len(df):,}")
    logger.info(f"   File size: {file_size_mb:.1f} MB")
    logger.info(f"   Time: {elapsed/60:.1f} minutes")


def main():
    """Main execution."""
    from variant_utils.utils import read_external_config

    parser = argparse.ArgumentParser(
        description='Preprocess gnomAD VCFs to Parquet for fast queries'
    )
    parser.add_argument('--chromosomes', type=str,
                       help='Comma-separated list of chromosomes (e.g., 1,2,3,X,Y)')
    parser.add_argument('--all-chromosomes', action='store_true',
                       help='Process all chromosomes (1-22, X, Y)')
    args = parser.parse_args()

    if not args.chromosomes and not args.all_chromosomes:
        parser.error("Must specify --chromosomes or --all-chromosomes")

    # Determine chromosomes to process
    if args.all_chromosomes:
        chromosomes = [str(i) for i in range(1, 23)] + ['X', 'Y']
    else:
        chromosomes = args.chromosomes.split(',')

    logger.info("\n" + "="*80)
    logger.info("GNOMAD PARQUET PREPROCESSING")
    logger.info("="*80)
    logger.info(f"Assembly: {ASSEMBLY}")
    logger.info(f"Chromosomes: {', '.join(chromosomes)}")
    logger.info(f"Output directory: {OUTPUT_DIR}")
    logger.info("")
    logger.info("NOTE: This is optional preprocessing for repeated queries.")
    logger.info("      Not needed if you only run the pipeline once.")
    logger.info("="*80)

    # Load external config
    external_config = read_external_config(EXTERNAL_TOOLS_CONFIG)

    # Process each chromosome
    total_start = time.time()

    for chrom in chromosomes:
        try:
            preprocess_chromosome(chrom, ASSEMBLY, external_config)
        except Exception as e:
            logger.error(f"Failed to process chromosome {chrom}: {e}")
            continue

    total_elapsed = time.time() - total_start

    logger.info("\n" + "="*80)
    logger.info("PREPROCESSING COMPLETE")
    logger.info("="*80)
    logger.info(f"Total time: {total_elapsed/3600:.2f} hours")
    logger.info(f"Output directory: {OUTPUT_DIR}")
    logger.info("")
    logger.info("To use preprocessed data:")
    logger.info("  - Modify query code to read from Parquet instead of VCF")
    logger.info("  - Gene queries will be 10-50x faster")
    logger.info("="*80)


if __name__ == "__main__":
    main()
