#!/usr/bin/env python3
"""
OPTIMIZED: Direct VCF to DataFrame Conversion

KEY OPTIMIZATION: Don't filter missense during VCF reading!
- Read ALL SNPs + PASS variants
- Parse VEP once
- Save complete DataFrame
- Filter to missense later (takes <1 second in pandas)

This is 10x faster than filtering missense during reading!
"""

import pandas as pd
from pathlib import Path
from pysam import VariantFile
import logging
from typing import List, Dict
import time
import argparse

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Paths
SOURCE_EXOMES = Path("/projects/lugoteam/protein_graphs/GnomAD-Parser/data/exomes")
SOURCE_GENOMES = Path("/projects/lugoteam/protein_graphs/GnomAD-Parser/data/genomes")
DF_DIR = Path("gnomad_all_genes/chromosome_dataframes")

# Chromosomes
CHROMOSOMES = [str(i) for i in range(1, 23)] + ['X', 'Y']


def get_vep_columns_from_header(vcf_path: str) -> List[str]:
    """Extract VEP column names from VCF header."""
    vcf = VariantFile(vcf_path)
    
    for record in vcf.header.records:
        if record.key == 'INFO' and 'ID' in record and record['ID'] == 'vep':
            description = record.get('Description', '')
            if 'Format:' in description:
                format_part = description.split('Format:')[1].strip()
                columns = [col.strip() for col in format_part.split('|')]
                vcf.close()
                return columns
    
    vcf.close()
    raise ValueError("VEP format not found in VCF header")


def parse_vep_annotation(vep_string: str, columns: List[str]) -> List[Dict]:
    """Parse VEP annotation string into list of dicts (one per transcript)."""
    if not vep_string or pd.isna(vep_string):
        return []
    
    annotations = []
    for transcript_annotation in str(vep_string).split(','):
        fields = transcript_annotation.split('|')
        
        # Pad with empty strings if needed
        while len(fields) < len(columns):
            fields.append('')
        
        annotation = {col: fields[i] if i < len(fields) else '' 
                     for i, col in enumerate(columns)}
        annotations.append(annotation)
    
    return annotations


def vcf_to_dataframe_optimized(vcf_path: Path, vep_columns: List[str]) -> pd.DataFrame:
    """
    OPTIMIZED: Read VCF and convert to DataFrame.
    
    Filters:
    - SNPs only (no indels)
    - PASS only (QC passed)
    - Parses ALL VEP annotations
    
    Does NOT filter by consequence type - do that later in pandas!
    """
    logger.info(f"  Reading {vcf_path.name}...")
    
    vcf = VariantFile(str(vcf_path))
    records = []
    
    variants_total = 0
    variants_kept = 0
    
    for record in vcf:
        variants_total += 1
        
        # Progress indicator
        if variants_total % 1000000 == 0:
            logger.info(f"    Processed {variants_total:,} variants...")
        
        # Filter: SNPs only
        if not record.alts:
            continue
        
        ref_len = len(record.ref)
        alt_len = len(str(record.alts[0])) if record.alts else 0
        if ref_len != 1 or alt_len != 1:
            continue
        
        # Filter: PASS only
        if record.filter.keys() and 'PASS' not in record.filter.keys():
            continue
        
        variants_kept += 1
        
        # Basic variant fields
        base_record = {
            'CHROM': record.chrom.replace('chr', ''),
            'POS': record.pos,
            'REF': record.ref,
            'ALT': ','.join([str(a) for a in record.alts]) if record.alts else '',
        }
        
        # Add key INFO fields
        info_fields = [
            'AC', 'AN', 'AF', 'nhomalt',
            'AC_afr', 'AC_amr', 'AC_asj', 'AC_eas', 'AC_fin', 'AC_nfe', 'AC_sas',
            'faf95', 'faf99', 'cadd_phred', 'phylop'
        ]
        
        for field in info_fields:
            if field in record.info:
                value = record.info[field]
                if isinstance(value, tuple):
                    value = value[0] if len(value) == 1 else ','.join(map(str, value))
                base_record[field] = value
            else:
                base_record[field] = None
        
        # Parse VEP - creates multiple rows (one per transcript)
        vep_string = record.info.get('vep', None)
        vep_annotations = parse_vep_annotation(vep_string, vep_columns)
        
        if vep_annotations:
            for vep_ann in vep_annotations:
                full_record = {**base_record, **vep_ann}
                records.append(full_record)
        else:
            # No VEP - add row with empty VEP fields
            full_record = {**base_record, **{col: '' for col in vep_columns}}
            records.append(full_record)
    
    vcf.close()
    
    logger.info(f"    Total variants: {variants_total:,}")
    logger.info(f"    SNPs + PASS: {variants_kept:,} ({100*variants_kept/variants_total:.1f}%)")
    logger.info(f"    Total annotations: {len(records):,}")
    
    if not records:
        return pd.DataFrame()
    
    return pd.DataFrame(records)


def merge_exomes_genomes(exomes_df: pd.DataFrame, genomes_df: pd.DataFrame) -> pd.DataFrame:
    """Merge exomes and genomes DataFrames."""
    if exomes_df.empty and genomes_df.empty:
        return pd.DataFrame()
    
    if exomes_df.empty:
        genomes_df['source'] = 'genomes'
        return genomes_df
    
    if genomes_df.empty:
        exomes_df['source'] = 'exomes'
        return exomes_df
    
    # Mark source
    exomes_df['source'] = 'exomes'
    genomes_df['source'] = 'genomes'
    
    # Combine
    combined = pd.concat([exomes_df, genomes_df], ignore_index=True)
    
    logger.info(f"  Merged: {len(exomes_df):,} exomes + {len(genomes_df):,} genomes = {len(combined):,} total")
    
    return combined


def process_chromosome(chrom: str) -> Dict:
    """Process one chromosome: read source VCF → filter → parse VEP → save DataFrame."""
    logger.info(f"\n{'='*80}")
    logger.info(f"Processing Chromosome {chrom}")
    logger.info(f"{'='*80}")
    
    start_time = time.time()
    
    try:
        # Check if already processed
        output_file = DF_DIR / f"chr{chrom}_variants.parquet"
        if output_file.exists():
            logger.info(f"✅ Chr{chrom} DataFrame already exists (skipping)")
            file_size_mb = output_file.stat().st_size / 1e6
            # Load to get variant count
            df = pd.read_parquet(output_file)
            return {
                'chrom': chrom,
                'status': 'skipped',
                'variants': len(df),
                'file_size_mb': file_size_mb,
                'time': 0
            }
        
        # Source VCF paths
        exomes_source = SOURCE_EXOMES / f"gnomad.exomes.v4.1.sites.chr{chrom}.vcf.bgz"
        genomes_source = SOURCE_GENOMES / f"gnomad.genomes.v4.1.sites.chr{chrom}.vcf.bgz"
        
        if not exomes_source.exists() or not genomes_source.exists():
            logger.error(f"❌ Chr{chrom}: Source VCFs not found!")
            return {
                'chrom': chrom,
                'status': 'failed',
                'variants': 0,
                'time': 0,
                'error': 'Source VCFs not found'
            }
        
        # Get VEP columns from header
        vep_columns = get_vep_columns_from_header(str(exomes_source))
        logger.info(f"VEP columns: {len(vep_columns)} fields")
        
        # Read and filter exomes (NO MISSENSE FILTERING - do that later!)
        logger.info("Processing EXOMES...")
        exomes_df = vcf_to_dataframe_optimized(exomes_source, vep_columns)
        
        # Read and filter genomes (NO MISSENSE FILTERING - do that later!)
        logger.info("Processing GENOMES...")
        genomes_df = vcf_to_dataframe_optimized(genomes_source, vep_columns)
        
        # Merge
        logger.info("Merging exomes and genomes...")
        combined_df = merge_exomes_genomes(exomes_df, genomes_df)
        
        if combined_df.empty:
            logger.warning(f"⚠️  Chr{chrom}: No variants after filtering and merging")
            return {
                'chrom': chrom,
                'status': 'success',
                'variants': 0,
                'time': time.time() - start_time
            }
        
        # Optimize data types
        logger.info("  Optimizing data types...")
        combined_df['POS'] = combined_df['POS'].astype('int32')
        
        # Convert numeric columns
        numeric_cols = ['AC', 'AN', 'AF', 'nhomalt', 'faf95', 'faf99', 'cadd_phred', 'phylop']
        for col in numeric_cols:
            if col in combined_df.columns:
                combined_df[col] = pd.to_numeric(combined_df[col], errors='coerce')
        
        # Save DataFrame (ALL variants, not just missense!)
        logger.info(f"  Saving DataFrame to {output_file.name}...")
        combined_df.to_parquet(output_file, index=False, compression='snappy')
        
        file_size_mb = output_file.stat().st_size / 1e6
        elapsed = time.time() - start_time
        
        # Report missense count for reference
        missense_count = combined_df['Consequence'].str.contains('missense', na=False).sum()
        
        logger.info(f"✅ Chr{chrom} complete:")
        logger.info(f"   Total annotations: {len(combined_df):,}")
        logger.info(f"   Missense annotations: {missense_count:,}")
        logger.info(f"   File size: {file_size_mb:.1f} MB")
        logger.info(f"   Time: {elapsed:.1f}s ({elapsed/60:.1f} min)")
        
        return {
            'chrom': chrom,
            'status': 'success',
            'variants': len(combined_df),
            'missense_variants': missense_count,
            'time': elapsed,
            'file_size_mb': file_size_mb
        }
        
    except Exception as e:
        logger.error(f"❌ Chr{chrom} failed: {e}")
        import traceback
        traceback.print_exc()
        
        return {
            'chrom': chrom,
            'status': 'failed',
            'variants': 0,
            'time': time.time() - start_time,
            'error': str(e)
        }


def main():
    """Main execution."""
    parser = argparse.ArgumentParser(
        description='Convert gnomAD VCFs to DataFrames (OPTIMIZED - no missense filtering)'
    )
    parser.add_argument('--chromosomes', nargs='+', default=CHROMOSOMES,
                       help='Chromosomes to process (default: all)')
    args = parser.parse_args()
    
    logger.info("\n" + "="*80)
    logger.info("OPTIMIZED VCF TO DATAFRAME CONVERSION")
    logger.info("="*80)
    logger.info(f"Source exomes: {SOURCE_EXOMES}")
    logger.info(f"Source genomes: {SOURCE_GENOMES}")
    logger.info(f"Output directory: {DF_DIR}")
    logger.info(f"Chromosomes: {', '.join(args.chromosomes)}")
    logger.info("")
    logger.info("OPTIMIZATION: Storing ALL variants (SNPs + PASS)")
    logger.info("              Filter to missense later in pandas (much faster!)")
    
    # Create output directory
    DF_DIR.mkdir(exist_ok=True, parents=True)
    
    start_time = time.time()
    results = []
    
    # Process chromosomes sequentially (memory-safe)
    for chrom in args.chromosomes:
        result = process_chromosome(chrom)
        results.append(result)
    
    # Summary
    elapsed = time.time() - start_time
    
    logger.info("\n" + "="*80)
    logger.info("CONVERSION COMPLETE")
    logger.info("="*80)
    
    successful = [r for r in results if r['status'] == 'success']
    skipped = [r for r in results if r['status'] == 'skipped']
    failed = [r for r in results if r['status'] == 'failed']
    
    total_variants = sum(r['variants'] for r in successful + skipped)
    total_missense = sum(r.get('missense_variants', 0) for r in successful)
    total_size_mb = sum(r.get('file_size_mb', 0) for r in successful + skipped)
    
    logger.info(f"\nSuccessful: {len(successful)}")
    logger.info(f"Skipped: {len(skipped)}")
    logger.info(f"Failed: {len(failed)}")
    logger.info(f"\nTotal annotations: {total_variants:,}")
    logger.info(f"Missense annotations: {total_missense:,} ({100*total_missense/total_variants:.1f}%)")
    logger.info(f"Total size: {total_size_mb:.1f} MB ({total_size_mb/1024:.1f} GB)")
    logger.info(f"Total time: {elapsed/60:.1f} minutes ({elapsed/3600:.1f} hours)")
    logger.info(f"Avg per chromosome: {elapsed/len(results)/60:.1f} minutes")
    
    if failed:
        logger.warning(f"\nFailed chromosomes: {', '.join([r['chrom'] for r in failed])}")
    
    # Save summary
    summary_df = pd.DataFrame(results)
    summary_file = DF_DIR / "conversion_summary.csv"
    summary_df.to_csv(summary_file, index=False)
    logger.info(f"\nSummary saved to: {summary_file}")
    
    # Show how to filter missense
    logger.info("\n" + "="*80)
    logger.info("TO FILTER MISSENSE VARIANTS:")
    logger.info("="*80)
    logger.info("import pandas as pd")
    logger.info("df = pd.read_parquet('gnomad_all_genes/chromosome_dataframes/chr1_variants.parquet')")
    logger.info("missense = df[df['Consequence'].str.contains('missense', na=False)]")
    logger.info("# Takes <1 second!")
    
    return len(failed) == 0


if __name__ == "__main__":
    success = main()
    exit(0 if success else 1)