#!/usr/bin/env python3
"""
Direct VCF to DataFrame Conversion (No Intermediate Filtering)

Streamlined approach:
1. Read source VCF directly
2. Filter in Python (SNPs + PASS only)
3. Parse VEP annotations
4. Save as DataFrame (Parquet)
5. Optionally save filtered VCF cache too

This skips the bcftools filtering step - everything happens in Python.
"""

import pandas as pd
from pathlib import Path
from pysam import VariantFile
import logging
from typing import List, Dict, Optional
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
CACHE_DIR = Path("gnomad_all_genes/chromosome_cache")

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


def vcf_to_dataframe_direct(
    vcf_path: Path,
    vep_columns: List[str],
    filter_snps: bool = True,
    filter_pass: bool = True,
    filter_missense: bool = True
) -> pd.DataFrame:
    """
    Read VCF directly and convert to DataFrame with filtering.

    Parameters:
    -----------
    vcf_path : Path
        Source VCF file path
    vep_columns : List[str]
        VEP column names
    filter_snps : bool
        Only keep SNPs (filter out indels)
    filter_pass : bool
        Only keep variants that passed QC filters
    filter_missense : bool
        Only keep missense variants
    """
    logger.info(f"  Reading {vcf_path.name}...")
    logger.info(f"    Filters: SNPs={filter_snps}, PASS={filter_pass}, Missense={filter_missense}")
    
    vcf = VariantFile(str(vcf_path))
    records = []
    
    variants_total = 0
    variants_kept = 0
    
    for record in vcf:
        variants_total += 1
        
        # Filter: SNPs only
        if filter_snps and not record.alts:
            continue
        if filter_snps:
            # Check if it's a SNP (single nucleotide change)
            ref_len = len(record.ref)
            alt_len = len(str(record.alts[0])) if record.alts else 0
            if ref_len != 1 or alt_len != 1:
                continue
        
        # Filter: PASS only
        if filter_pass:
            # Keep if no filters (empty = PASS) OR if PASS is explicitly listed
            if record.filter.keys() and 'PASS' not in record.filter.keys():
                continue
        
        variants_kept += 1
        
        # Basic variant fields
        base_record = {
            'CHROM': record.chrom.replace('chr', ''),
            'POS': record.pos,
            'REF': record.ref,
            'ALT': ','.join([str(a) for a in record.alts]) if record.alts else '',
            'FILTER': ','.join(record.filter.keys()) if record.filter.keys() else 'PASS',
        }
        
        # Add key INFO fields (gnomAD-specific)
        info_fields = [
            'AC', 'AN', 'AF', 'nhomalt',
            'AC_afr', 'AN_afr', 'AF_afr',
            'AC_amr', 'AN_amr', 'AF_amr',
            'AC_asj', 'AN_asj', 'AF_asj',
            'AC_eas', 'AN_eas', 'AF_eas',
            'AC_fin', 'AN_fin', 'AF_fin',
            'AC_nfe', 'AN_nfe', 'AF_nfe',
            'AC_sas', 'AN_sas', 'AF_sas',
            'faf95', 'faf99',
            'cadd_phred', 'cadd_raw_score',
            'phylop', 'pangolin_largest_ds'
        ]
        
        for field in info_fields:
            if field in record.info:
                value = record.info[field]
                # Handle tuples (some fields return tuples)
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
                # Filter: Missense variants only
                if filter_missense:
                    consequence = vep_ann.get('Consequence', '')
                    if consequence:
                        consequences = consequence.split('&')
                        if 'missense_variant' not in consequences:
                            continue
                    else:
                        # No consequence annotation, skip
                        continue

                full_record = {**base_record, **vep_ann}
                records.append(full_record)
        else:
            # No VEP annotation - skip if filtering for missense
            if not filter_missense:
                full_record = {**base_record, **{col: '' for col in vep_columns}}
                records.append(full_record)
    
    vcf.close()
    
    logger.info(f"    Total variants: {variants_total:,}")
    logger.info(f"    Kept after filtering: {variants_kept:,} ({100*variants_kept/variants_total:.1f}%)")
    logger.info(f"    Annotations (with VEP): {len(records):,}")
    
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


def save_filtered_vcf_cache(
    source_vcf: Path,
    cache_vcf: Path,
    filter_snps: bool = True,
    filter_pass: bool = True
):
    """
    Optionally save filtered VCF to cache for compatibility with original pipeline.
    This is a bonus - not required for DataFrame approach.
    """
    if cache_vcf.exists():
        logger.info(f"  Cache VCF already exists: {cache_vcf.name}")
        return
    
    logger.info(f"  Creating cache VCF: {cache_vcf.name}")
    
    vcf_in = VariantFile(str(source_vcf))
    vcf_out = VariantFile(str(cache_vcf), 'w', header=vcf_in.header)
    
    for record in vcf_in:
        # Same filters as above
        if filter_snps and not record.alts:
            continue
        if filter_snps:
            ref_len = len(record.ref)
            alt_len = len(str(record.alts[0])) if record.alts else 0
            if ref_len != 1 or alt_len != 1:
                continue
        
        if filter_pass:
            if record.filter.keys() and 'PASS' not in record.filter.keys():
                continue
        
        vcf_out.write(record)
    
    vcf_in.close()
    vcf_out.close()
    
    # Index
    import pysam
    pysam.tabix_index(str(cache_vcf), preset='vcf', force=True)
    logger.info(f"  Cached and indexed: {cache_vcf.name}")


def process_chromosome(chrom: str, save_vcf_cache: bool = False) -> Dict:
    """
    Process one chromosome: read source VCF → filter → parse VEP → save DataFrame.
    
    Parameters:
    -----------
    chrom : str
        Chromosome (e.g., "1", "X")
    save_vcf_cache : bool
        If True, also save filtered VCF cache (for compatibility with original pipeline)
    """
    logger.info(f"\n{'='*80}")
    logger.info(f"Processing Chromosome {chrom}")
    logger.info(f"{'='*80}")
    
    start_time = time.time()
    
    try:
        # Check if already processed
        output_file = DF_DIR / f"chr{chrom}_variants.parquet"
        if output_file.exists():
            logger.info(f"✅ Chr{chrom} DataFrame already exists (skipping)")
            return {
                'chrom': chrom,
                'status': 'skipped',
                'variants': 0,
                'time': 0
            }
        
        # Source VCF paths
        exomes_source = SOURCE_EXOMES / f"gnomad.exomes.v4.1.sites.chr{chrom}.vcf.bgz"
        genomes_source = SOURCE_GENOMES / f"gnomad.genomes.v4.1.sites.chr{chrom}.vcf.bgz"

        # Check if VCF files exist
        if not exomes_source.exists() or not genomes_source.exists():
            logger.error(f"❌ Chr{chrom}: Source VCFs not found!")
            return {
                'chrom': chrom,
                'status': 'failed',
                'variants': 0,
                'time': 0,
                'error': 'Source VCFs not found'
            }

        # Check if index files exist (.tbi or .csi)
        for vcf_path, label in [(exomes_source, 'exomes'), (genomes_source, 'genomes')]:
            tbi_index = Path(str(vcf_path) + '.tbi')
            csi_index = Path(str(vcf_path) + '.csi')
            if not tbi_index.exists() and not csi_index.exists():
                logger.error(f"❌ Chr{chrom}: Index missing for {label} VCF: {vcf_path}")
                logger.error(f"   Please run: tabix -p vcf {vcf_path}")
                return {
                    'chrom': chrom,
                    'status': 'failed',
                    'variants': 0,
                    'time': 0,
                    'error': f'Index missing for {label} VCF'
                }
        
        # Get VEP columns from header
        vep_columns = get_vep_columns_from_header(str(exomes_source))
        logger.info(f"VEP columns: {len(vep_columns)} fields")
        
        # Read and filter exomes
        logger.info("Processing EXOMES...")
        exomes_df = vcf_to_dataframe_direct(
            exomes_source,
            vep_columns,
            filter_snps=True,
            filter_pass=True,
            filter_missense=True
        )

        # Read and filter genomes
        logger.info("Processing GENOMES...")
        genomes_df = vcf_to_dataframe_direct(
            genomes_source,
            vep_columns,
            filter_snps=True,
            filter_pass=True,
            filter_missense=True
        )
        
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
        numeric_cols = ['AC', 'AN', 'AF', 'nhomalt', 'faf95', 'faf99', 
                       'cadd_phred', 'phylop']
        for col in numeric_cols:
            if col in combined_df.columns:
                combined_df[col] = pd.to_numeric(combined_df[col], errors='coerce')
        
        # Save DataFrame
        logger.info(f"  Saving DataFrame to {output_file.name}...")
        combined_df.to_parquet(output_file, index=False, compression='snappy')
        
        file_size_mb = output_file.stat().st_size / 1e6
        
        # Optionally save VCF cache too (for original pipeline compatibility)
        if save_vcf_cache:
            logger.info("  Creating VCF cache (optional, for original pipeline)...")
            CACHE_DIR.mkdir(exist_ok=True, parents=True)
            
            cache_exomes = CACHE_DIR / f"gnomad.exomes.GRCh38.chr{chrom}.vcf.gz"
            cache_genomes = CACHE_DIR / f"gnomad.genomes.GRCh38.chr{chrom}.vcf.gz"
            
            save_filtered_vcf_cache(exomes_source, cache_exomes)
            save_filtered_vcf_cache(genomes_source, cache_genomes)
        
        elapsed = time.time() - start_time
        
        logger.info(f"✅ Chr{chrom} complete: {len(combined_df):,} annotations, {file_size_mb:.1f} MB, {elapsed:.1f}s")
        
        return {
            'chrom': chrom,
            'status': 'success',
            'variants': len(combined_df),
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
        description='Convert gnomAD VCFs to DataFrames (direct, no intermediate filtering)'
    )
    parser.add_argument('--save-vcf-cache', action='store_true',
                       help='Also save filtered VCF cache (for original pipeline compatibility)')
    parser.add_argument('--chromosomes', nargs='+', default=CHROMOSOMES,
                       help='Chromosomes to process (default: all)')
    args = parser.parse_args()
    
    logger.info("\n" + "="*80)
    logger.info("DIRECT VCF TO DATAFRAME CONVERSION")
    logger.info("="*80)
    logger.info(f"Source exomes: {SOURCE_EXOMES}")
    logger.info(f"Source genomes: {SOURCE_GENOMES}")
    logger.info(f"Output directory: {DF_DIR}")
    logger.info(f"Chromosomes: {', '.join(args.chromosomes)}")
    logger.info(f"Save VCF cache: {args.save_vcf_cache}")
    
    # Create output directory
    DF_DIR.mkdir(exist_ok=True, parents=True)
    
    start_time = time.time()
    results = []
    
    # Process chromosomes sequentially (memory-safe)
    for chrom in args.chromosomes:
        result = process_chromosome(chrom, save_vcf_cache=args.save_vcf_cache)
        results.append(result)
    
    # Summary
    elapsed = time.time() - start_time
    
    logger.info("\n" + "="*80)
    logger.info("CONVERSION COMPLETE")
    logger.info("="*80)
    
    successful = [r for r in results if r['status'] == 'success']
    skipped = [r for r in results if r['status'] == 'skipped']
    failed = [r for r in results if r['status'] == 'failed']
    
    total_variants = sum(r['variants'] for r in successful)
    total_size_mb = sum(r.get('file_size_mb', 0) for r in successful)
    
    logger.info(f"\nSuccessful: {len(successful)}")
    logger.info(f"Skipped: {len(skipped)}")
    logger.info(f"Failed: {len(failed)}")
    logger.info(f"\nTotal annotations: {total_variants:,}")
    logger.info(f"Total size: {total_size_mb:.1f} MB ({total_size_mb/1024:.1f} GB)")
    logger.info(f"Total time: {elapsed/60:.1f} minutes ({elapsed/3600:.1f} hours)")
    
    if failed:
        logger.warning(f"\nFailed chromosomes: {', '.join([r['chrom'] for r in failed])}")
    
    # Save summary
    summary_df = pd.DataFrame(results)
    summary_file = DF_DIR / "conversion_summary.csv"
    summary_df.to_csv(summary_file, index=False)
    logger.info(f"\nSummary saved to: {summary_file}")
    
    return len(failed) == 0


if __name__ == "__main__":
    success = main()
    exit(0 if success else 1)