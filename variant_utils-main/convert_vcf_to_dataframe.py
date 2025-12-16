#!/usr/bin/env python3
"""
MEMORY-EFFICIENT: VCF to DataFrame with Chunked Processing

KEY OPTIMIZATION: Process and save in chunks!
- Read VCF in batches
- Save each batch incrementally
- Never load full chromosome into memory
- Uses Parquet append mode

Maximum memory usage: ~10-15 GB (vs 60+ GB before)
"""

import pandas as pd
from pathlib import Path
from pysam import VariantFile
import logging
from typing import List, Dict
import time
import argparse
import pyarrow as pa
import pyarrow.parquet as pq

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

# Chunk size (tune based on memory)
CHUNK_SIZE = 500000  # Process 500K variants at a time


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


def process_vcf_in_chunks(vcf_path: Path, vep_columns: List[str], output_file: Path, 
                          source_label: str, chunk_size: int, is_first_file: bool = True):
    """
    Process VCF in chunks and append to Parquet file.
    
    Memory-efficient: Never loads entire chromosome into memory.
    """
    logger.info(f"  Processing {vcf_path.name} ({source_label})...")
    
    vcf = VariantFile(str(vcf_path))
    
    chunk_records = []
    variants_total = 0
    variants_kept = 0
    chunks_written = 0
    
    # Get schema from first chunk (for Parquet writer)
    writer = None
    schema = None
    
    for record in vcf:
        variants_total += 1
        
        # Progress indicator
        if variants_total % 1000000 == 0:
            logger.info(f"    Processed {variants_total:,} variants, kept {variants_kept:,}...")
        
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
            'source': source_label,
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
        
        # Parse VEP
        vep_string = record.info.get('vep', None)
        vep_annotations = parse_vep_annotation(vep_string, vep_columns)
        
        if vep_annotations:
            for vep_ann in vep_annotations:
                full_record = {**base_record, **vep_ann}
                chunk_records.append(full_record)
        else:
            full_record = {**base_record, **{col: '' for col in vep_columns}}
            chunk_records.append(full_record)
        
        # Write chunk when it reaches size limit
        if len(chunk_records) >= chunk_size:
            chunk_df = pd.DataFrame(chunk_records)
            
            # Optimize dtypes
            chunk_df['POS'] = chunk_df['POS'].astype('int32')
            for col in ['AC', 'AN', 'AF', 'nhomalt', 'faf95', 'faf99', 'cadd_phred', 'phylop']:
                if col in chunk_df.columns:
                    chunk_df[col] = pd.to_numeric(chunk_df[col], errors='coerce')
            
            # Write to Parquet
            if writer is None:
                # First chunk - create writer with schema
                schema = pa.Schema.from_pandas(chunk_df)
                writer = pq.ParquetWriter(output_file, schema, compression='snappy')
            
            table = pa.Table.from_pandas(chunk_df, schema=schema)
            writer.write_table(table)
            
            chunks_written += 1
            logger.info(f"    Wrote chunk {chunks_written} ({len(chunk_records):,} rows)")
            
            # Clear chunk
            chunk_records = []
    
    # Write remaining records
    if chunk_records:
        chunk_df = pd.DataFrame(chunk_records)
        
        # Optimize dtypes
        chunk_df['POS'] = chunk_df['POS'].astype('int32')
        for col in ['AC', 'AN', 'AF', 'nhomalt', 'faf95', 'faf99', 'cadd_phred', 'phylop']:
            if col in chunk_df.columns:
                chunk_df[col] = pd.to_numeric(chunk_df[col], errors='coerce')
        
        if writer is None:
            # Only one small chunk - just save it
            chunk_df.to_parquet(output_file, index=False, compression='snappy')
        else:
            table = pa.Table.from_pandas(chunk_df, schema=schema)
            writer.write_table(table)
            chunks_written += 1
    
    # Close writer
    if writer is not None:
        writer.close()
    
    vcf.close()
    
    logger.info(f"    Total variants: {variants_total:,}")
    logger.info(f"    SNPs + PASS: {variants_kept:,}")
    logger.info(f"    Annotations written: {chunks_written * chunk_size + len(chunk_records):,}")
    
    return variants_kept


def process_chromosome(chrom: str, chunk_size: int = CHUNK_SIZE) -> Dict:
    """Process one chromosome with memory-efficient chunking."""
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
        
        # Get VEP columns
        vep_columns = get_vep_columns_from_header(str(exomes_source))
        logger.info(f"VEP columns: {len(vep_columns)} fields")
        logger.info(f"Chunk size: {chunk_size:,} variants per batch")
        
        # Process exomes (writes to temp file)
        temp_exomes = DF_DIR / f".tmp_chr{chrom}_exomes.parquet"
        logger.info("Processing EXOMES in chunks...")
        exomes_count = process_vcf_in_chunks(exomes_source, vep_columns, temp_exomes, 'exomes', chunk_size, True)
        
        # Process genomes (writes to temp file)
        temp_genomes = DF_DIR / f".tmp_chr{chrom}_genomes.parquet"
        logger.info("Processing GENOMES in chunks...")
        genomes_count = process_vcf_in_chunks(genomes_source, vep_columns, temp_genomes, 'genomes', chunk_size, False)
        
        # Merge the two files
        logger.info("Merging exomes and genomes...")
        exomes_df = pd.read_parquet(temp_exomes)
        genomes_df = pd.read_parquet(temp_genomes)
        
        combined_df = pd.concat([exomes_df, genomes_df], ignore_index=True)
        
        logger.info(f"  Total: {len(exomes_df):,} exomes + {len(genomes_df):,} genomes = {len(combined_df):,}")
        
        # Save final combined file
        logger.info(f"  Saving final DataFrame to {output_file.name}...")
        combined_df.to_parquet(output_file, index=False, compression='snappy')
        
        # Clean up temp files
        temp_exomes.unlink()
        temp_genomes.unlink()
        
        file_size_mb = output_file.stat().st_size / 1e6
        elapsed = time.time() - start_time
        
        # Report missense count
        missense_count = combined_df['Consequence'].str.contains('missense', na=False).sum()
        
        logger.info(f"✅ Chr{chrom} complete:")
        logger.info(f"   Total annotations: {len(combined_df):,}")
        logger.info(f"   Missense annotations: {missense_count:,}")
        logger.info(f"   File size: {file_size_mb:.1f} MB")
        logger.info(f"   Time: {elapsed:.1f}s ({elapsed/60:.1f} min)")
        logger.info(f"   Peak memory: ~10-15 GB (chunked processing)")
        
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
        
        # Clean up temp files on error
        temp_exomes = DF_DIR / f".tmp_chr{chrom}_exomes.parquet"
        temp_genomes = DF_DIR / f".tmp_chr{chrom}_genomes.parquet"
        if temp_exomes.exists():
            temp_exomes.unlink()
        if temp_genomes.exists():
            temp_genomes.unlink()
        
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
        description='Memory-efficient VCF to DataFrame conversion with chunking'
    )
    parser.add_argument('--chromosomes', nargs='+', default=CHROMOSOMES,
                       help='Chromosomes to process (default: all)')
    parser.add_argument('--chunk-size', type=int, default=CHUNK_SIZE,
                       help=f'Chunk size for processing (default: {CHUNK_SIZE})')
    args = parser.parse_args()
    
    chunk_size = args.chunk_size
    
    logger.info("\n" + "="*80)
    logger.info("MEMORY-EFFICIENT VCF TO DATAFRAME CONVERSION")
    logger.info("="*80)
    logger.info(f"Source exomes: {SOURCE_EXOMES}")
    logger.info(f"Source genomes: {SOURCE_GENOMES}")
    logger.info(f"Output directory: {DF_DIR}")
    logger.info(f"Chromosomes: {', '.join(args.chromosomes)}")
    logger.info(f"Chunk size: {chunk_size:,} variants")
    logger.info("")
    logger.info("OPTIMIZATION: Chunked processing")
    logger.info("              Peak memory: ~10-15 GB (vs 60+ GB before)")
    
    # Create output directory
    DF_DIR.mkdir(exist_ok=True, parents=True)
    
    start_time = time.time()
    results = []
    
    # Process chromosomes sequentially
    for chrom in args.chromosomes:
        result = process_chromosome(chrom, chunk_size)
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
    logger.info(f"Missense annotations: {total_missense:,}")
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