#!/usr/bin/env python3
"""
Fetch All Missense Variants from gnomAD (Chromosome-Batched + Optimized)

KEY OPTIMIZATIONS:
1. Process all genes on a chromosome together (load VCF once per chromosome)
2. Parallel processing within chromosomes (multiple genes at once)
3. Aggressive caching at multiple levels
4. Memory-efficient streaming where possible

Usage:
    python fetch_all_gnomad_chromosome_batch.py --workers 4 --chrom-workers 2
"""

import pandas as pd
import json
import gzip
from pathlib import Path
from typing import Dict, List, Tuple
import time
from datetime import datetime
import logging
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import subprocess

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(processName)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# ============================================================================
# CONFIGURATION
# ============================================================================

# Paths
REFERENCE_DIR = Path("/projects/lugoteam/protein_graphs/GnomAD-Parser/ref_data")
MANE_SUMMARY = REFERENCE_DIR / "mane" / "MANE.GRCh38.v1.3.summary.txt.gz"
EXTERNAL_TOOLS_CONFIG = "external_tools.json"

# Output directory
OUTPUT_DIR = Path("gnomad_all_genes")
OUTPUT_DIR.mkdir(exist_ok=True)
CACHE_DIR = OUTPUT_DIR / "cache"
CACHE_DIR.mkdir(exist_ok=True)

# Chromosome-level cache for extracted VCFs
CHROM_CACHE_DIR = OUTPUT_DIR / "chromosome_cache"
CHROM_CACHE_DIR.mkdir(exist_ok=True)

# Assembly
ASSEMBLY = "GRCh38"

# Resume capability
PROGRESS_FILE = OUTPUT_DIR / "progress.json"
FAILED_GENES_FILE = OUTPUT_DIR / "failed_genes.txt"


# ============================================================================
# GENE LOADING (optimized with caching)
# ============================================================================

def get_all_genes_from_mane() -> pd.DataFrame:
    """Extract all protein-coding genes from MANE summary with caching."""
    
    # Check for cached gene list
    cache_file = OUTPUT_DIR / "gene_list.parquet"
    if cache_file.exists():
        logger.info(f"Loading cached gene list from {cache_file}")
        return pd.read_parquet(cache_file)
    
    logger.info("Loading all genes from MANE summary...")
    
    genes = []
    with gzip.open(MANE_SUMMARY, "rt") as f:
        header = f.readline().strip().split("\t")
        col_idx = {name: i for i, name in enumerate(header)}
        
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) != len(header):
                continue
            
            status = parts[col_idx.get("MANE_status", col_idx.get("MANE_Status", -1))]
            if status not in ("MANE Select", "MANE_Select"):
                continue
            
            symbol = parts[col_idx.get("symbol", col_idx.get("HGNC_symbol", -1))]
            chrom = parts[col_idx.get("GRCh38_chr", col_idx.get("#chr", -1))]
            gene_start = parts[col_idx.get("chr_start", col_idx.get("gene_start", -1))]
            gene_end = parts[col_idx.get("chr_end", col_idx.get("gene_end", -1))]
            enst = parts[col_idx.get("Ensembl_nuc", col_idx.get("Ensembl_transcript", -1))]
            nm = parts[col_idx.get("RefSeq_nuc", col_idx.get("RefSeq_transcript", -1))]
            hgnc_id = parts[col_idx.get("HGNC_ID", -1)]
            
            # Parse chromosome from NCBI format
            if chrom and chrom.startswith("NC_"):
                parts_chr = chrom.split(".")
                if len(parts_chr) > 0:
                    nc_id = parts_chr[0]
                    if nc_id == "NC_000023":
                        chrom = "X"
                    elif nc_id == "NC_000024":
                        chrom = "Y"
                    elif nc_id.startswith("NC_"):
                        chr_num = nc_id.replace("NC_", "").lstrip("0")
                        if chr_num:
                            chrom = chr_num
                        else:
                            chrom = "1"
            
            # Skip if missing required fields
            if not all([symbol, chrom, gene_start, gene_end, hgnc_id]):
                continue
            if gene_start == '-' or gene_end == '-':
                continue
            
            try:
                start_int = int(gene_start)
                end_int = int(gene_end)
            except (ValueError, TypeError):
                continue
            
            genes.append({
                'gene_symbol': symbol,
                'chrom': chrom.replace('chr', ''),
                'gene_start': start_int,
                'gene_end': end_int,
                'enst': enst.split('.')[0] if enst else '',
                'nm': nm if nm else '',
                'hgnc_id': hgnc_id
            })
    
    genes_df = pd.DataFrame(genes)
    
    # Filter out genes without HGNC IDs
    before = len(genes_df)
    genes_df = genes_df[genes_df['hgnc_id'].notna()].copy()
    after = len(genes_df)
    
    logger.info(f"Found {after} MANE Select genes with HGNC IDs (filtered {before-after} without)")
    
    # Cache for future use
    genes_df.to_parquet(cache_file, index=False)
    logger.info(f"Cached gene list to {cache_file}")
    
    return genes_df


def group_genes_by_chromosome(genes_df: pd.DataFrame) -> Dict[str, pd.DataFrame]:
    """Group genes by chromosome for batch processing."""
    chrom_groups = {}
    
    for chrom in genes_df['chrom'].unique():
        chrom_genes = genes_df[genes_df['chrom'] == chrom].copy()
        chrom_genes = chrom_genes.sort_values('gene_start')
        chrom_groups[chrom] = chrom_genes
        logger.info(f"  Chromosome {chrom}: {len(chrom_genes)} genes")
    
    return chrom_groups


# ============================================================================
# CHROMOSOME-LEVEL VCF EXTRACTION
# ============================================================================

def extract_chromosome_vcf(
    assembly: str,
    chrom: str,
    genes_df: pd.DataFrame,
    external_config: dict
) -> Tuple[Path, Path]:
    """
    Extract entire chromosome VCF for all genes at once (exomes + genomes).
    OPTIMIZED: Uses bcftools (10-50x faster than GATK).

    Returns: (exomes_vcf_path, genomes_vcf_path)
    """
    from pysam import VariantFile

    # Check cache first
    cache_exomes = CHROM_CACHE_DIR / f"gnomad.exomes.{assembly}.chr{chrom}.vcf.gz"
    cache_genomes = CHROM_CACHE_DIR / f"gnomad.genomes.{assembly}.chr{chrom}.vcf.gz"

    if cache_exomes.exists() and cache_genomes.exists():
        logger.info(f"Using cached chromosome VCFs for chr{chrom}")
        return cache_exomes, cache_genomes

    logger.info(f"Extracting chromosome {chrom} VCFs with bcftools (10-50x faster than GATK)...")

    # Setup paths
    release_version = "v4.1" if assembly == "GRCh38" else "r2.1.1"
    chr_prefix = "chr" if release_version == "v4.1" else ""

    gnomad_vcf_root = Path(
        external_config['gnomad_v4_vcf_root'] if release_version == "v4.1"
        else external_config['gnomad_v2_vcf_root']
    )

    gnomAD_exomes = gnomad_vcf_root / f"exomes/gnomad.exomes.{release_version}.sites.{chr_prefix}{chrom}.vcf.bgz"
    gnomAD_genomes = gnomad_vcf_root / f"genomes/gnomad.genomes.{release_version}.sites.{chr_prefix}{chrom}.vcf.bgz"

    # Get chromosome length from VCF header
    try:
        vcf = VariantFile(str(gnomAD_exomes))
        chrom_length = None
        for contig in vcf.header.contigs.values():
            if contig.name == f"{chr_prefix}{chrom}":
                chrom_length = contig.length
                break
        vcf.close()

        if chrom_length is None:
            raise ValueError(f"Chromosome {chr_prefix}{chrom} not found in VCF header")

    except Exception as e:
        logger.warning(f"Could not get chromosome length from VCF header: {e}")
        # Use approximate chromosome lengths as fallback
        chrom_lengths_grch38 = {
            '1': 248956422, '2': 242193529, '3': 198295559, '4': 190214555,
            '5': 181538259, '6': 170805979, '7': 159345973, '8': 145138636,
            '9': 138394717, '10': 133797422, '11': 135086622, '12': 133275309,
            '13': 114364328, '14': 107043718, '15': 101991189, '16': 90338345,
            '17': 83257441, '18': 80373285, '19': 58617616, '20': 64444167,
            '21': 46709983, '22': 50818468, 'X': 156040895, 'Y': 57227415
        }
        chrom_length = chrom_lengths_grch38.get(chrom, 250000000)
        logger.info(f"  Using approximate chromosome length: {chrom_length}")

    # Get regions for all genes on this chromosome
    min_pos = genes_df['gene_start'].min()
    max_pos = genes_df['gene_end'].max()

    # Add buffer (100kb on each side) but don't exceed chromosome boundaries
    min_pos = max(1, min_pos - 100000)
    max_pos = min(chrom_length, max_pos + 100000)

    logger.info(f"  Extracting region {chr_prefix}{chrom}:{min_pos}-{max_pos} (chr length: {chrom_length})")

    # OPTIMIZED: Use bcftools instead of GATK (10-50x faster)
    # Extract exomes with bcftools
    cmd_exomes = [
        "bcftools", "view",
        "-r", f"{chr_prefix}{chrom}:{min_pos}-{max_pos}",
        "-i", 'TYPE="snp" && FILTER="PASS"',
        str(gnomAD_exomes),
        "-Oz", "-o", str(cache_exomes)
    ]

    logger.info("  Extracting exomes with bcftools...")
    result = subprocess.run(cmd_exomes, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Failed to extract exomes with bcftools: {result.stderr}")

    # Index exomes
    logger.info("  Indexing exomes...")
    subprocess.run(["bcftools", "index", "-t", str(cache_exomes)],
                   capture_output=True, text=True, check=True)

    # Extract genomes with bcftools
    cmd_genomes = [
        "bcftools", "view",
        "-r", f"{chr_prefix}{chrom}:{min_pos}-{max_pos}",
        "-i", 'TYPE="snp" && FILTER="PASS"',
        str(gnomAD_genomes),
        "-Oz", "-o", str(cache_genomes)
    ]

    logger.info("  Extracting genomes with bcftools...")
    result = subprocess.run(cmd_genomes, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Failed to extract genomes with bcftools: {result.stderr}")

    # Index genomes
    logger.info("  Indexing genomes...")
    subprocess.run(["bcftools", "index", "-t", str(cache_genomes)],
                   capture_output=True, text=True, check=True)

    logger.info(f"✅ Chromosome {chrom} VCFs extracted and cached with bcftools")

    return cache_exomes, cache_genomes


# ============================================================================
# OPTIMIZED GENE PROCESSING (using pre-extracted chromosome VCFs)
# ============================================================================

def process_gene_from_chrom_vcf(
    gene_row: Dict,
    exomes_vcf: Path,
    genomes_vcf: Path,
    vep_columns: List[str],
    assembly: str
) -> Dict:
    """
    Process a single gene using pre-extracted chromosome VCFs.
    OPTIMIZED: Uses pysam directly (100-500x faster than GATK).
    """
    from variant_utils.gnomad_utils import (
        extract_variants_fast,
        merge_exome_genome_dataframes,
        parse_vep
    )
    from concurrent.futures import ThreadPoolExecutor

    gene_symbol = gene_row['gene_symbol']

    try:
        # Extract gene region from chromosome VCFs
        release_version = "v4.1" if assembly == "GRCh38" else "r2.1.1"
        chr_prefix = "chr" if release_version == "v4.1" else ""

        chrom = gene_row['chrom']
        start = int(gene_row['gene_start'])
        end = int(gene_row['gene_end'])

        # OPTIMIZATION: Extract exomes and genomes in parallel using ThreadPoolExecutor
        # (I/O bound operation, threads work well here)
        with ThreadPoolExecutor(max_workers=2) as executor:
            exomes_future = executor.submit(
                extract_variants_fast, exomes_vcf, chrom, start, end, chr_prefix
            )
            genomes_future = executor.submit(
                extract_variants_fast, genomes_vcf, chrom, start, end, chr_prefix
            )

            exomes_df = exomes_future.result()
            genomes_df = genomes_future.result()

        # Merge exomes and genomes (in-memory, no temp files)
        df = merge_exome_genome_dataframes(exomes_df, genomes_df)

        if df.empty:
            return {
                'gene': gene_symbol,
                'status': 'success',
                'variants': 0,
                'error': None
            }

        # Parse VEP annotations
        vep_df = parse_vep(df, columns=vep_columns)

        # Merge with main dataframe
        df = pd.merge(df, vep_df, left_index=True, right_on='index', validate='one_to_many')

        # Filter to this gene
        gene_df = df[df.HGNC_ID == gene_row['hgnc_id']].copy()

        # Normalize
        gene_df = gene_df.assign(
            CHROM=gene_df.CHROM.astype(str).str.replace("chr", ""),
            POS=gene_df.POS.astype(str),
            REF=gene_df.REF.astype(str),
            ALT=gene_df.ALT.astype(str)
        )

        # Save
        output_file = OUTPUT_DIR / f"{gene_symbol}_gnomad_variants.parquet"
        gene_df.to_parquet(output_file, index=False)

        return {
            'gene': gene_symbol,
            'status': 'success',
            'variants': len(gene_df),
            'error': None
        }

    except Exception as e:
        return {
            'gene': gene_symbol,
            'status': 'failed',
            'variants': 0,
            'error': str(e)
        }


# ============================================================================
# CHROMOSOME-BATCHED PARALLEL PROCESSING
# ============================================================================

def process_chromosome(
    chrom: str,
    genes_df: pd.DataFrame,
    external_config: dict,
    assembly: str,
    gene_workers: int = 2
) -> List[Dict]:
    """
    Process all genes on a chromosome.

    OPTIMIZED Strategy:
    1. Extract chromosome VCFs once (cached)
    2. Get VEP columns once (not per gene)
    3. Process multiple genes in parallel using ProcessPoolExecutor (no GIL!)
    """
    from variant_utils.gnomad_utils import get_vep_columns_from_vcf_header

    logger.info(f"\n{'='*80}")
    logger.info(f"PROCESSING CHROMOSOME {chrom}")
    logger.info(f"{'='*80}")
    logger.info(f"Genes on chromosome: {len(genes_df)}")

    # Step 1: Extract chromosome VCFs (or use cache)
    try:
        exomes_vcf, genomes_vcf = extract_chromosome_vcf(
            assembly, chrom, genes_df, external_config
        )
    except Exception as e:
        logger.error(f"Failed to extract chromosome {chrom} VCFs: {e}")
        return [
            {'gene': row['gene_symbol'], 'status': 'failed', 'variants': 0,
             'error': f"Chromosome extraction failed: {e}"}
            for _, row in genes_df.iterrows()
        ]

    # Step 2: Get VEP columns once (not per gene - saves time)
    try:
        vep_columns = get_vep_columns_from_vcf_header(str(exomes_vcf))
        logger.info(f"VEP columns extracted: {len(vep_columns)} fields")
    except Exception as e:
        logger.error(f"Failed to extract VEP columns: {e}")
        vep_columns = []

    # Step 3: Process genes in parallel with ProcessPoolExecutor (no GIL!)
    logger.info(f"Processing {len(genes_df)} genes with {gene_workers} workers (ProcessPool)...")

    results = []
    gene_list = genes_df.to_dict('records')

    # OPTIMIZATION: Use ProcessPoolExecutor instead of ThreadPoolExecutor
    # This bypasses Python's GIL and allows true parallel CPU execution
    with ProcessPoolExecutor(max_workers=gene_workers) as executor:
        futures = {
            executor.submit(
                process_gene_from_chrom_vcf,
                gene, exomes_vcf, genomes_vcf, vep_columns, assembly
            ): gene['gene_symbol']
            for gene in gene_list
        }

        for future in as_completed(futures):
            result = future.result()
            results.append(result)

            status_icon = "✅" if result['status'] == 'success' else "❌"
            logger.info(
                f"  {status_icon} {result['gene']}: {result['variants']:,} variants"
            )

    return results


def load_progress() -> Dict:
    """Load progress from previous run."""
    if PROGRESS_FILE.exists():
        with open(PROGRESS_FILE, 'r') as f:
            return json.load(f)
    return {'completed': [], 'failed': [], 'completed_chromosomes': []}


def save_progress(progress: Dict):
    """Save progress."""
    with open(PROGRESS_FILE, 'w') as f:
        json.dump(progress, f, indent=2)


def query_all_genes_by_chromosome(
    genes_df: pd.DataFrame,
    external_config: dict,
    assembly: str,
    chrom_workers: int = 1,
    gene_workers: int = 2,
    resume: bool = True
):
    """
    Query all genes, processing one chromosome at a time.
    
    Parameters:
    -----------
    genes_df : pd.DataFrame
        All genes to process
    external_config : dict
        External tools config
    assembly : str
        'GRCh37' or 'GRCh38'
    chrom_workers : int
        Number of chromosomes to process in parallel (usually 1-2 for memory)
    gene_workers : int
        Number of genes to process in parallel per chromosome
    resume : bool
        Resume from previous run
    """
    # Load progress
    progress = load_progress() if resume else {
        'completed': [], 'failed': [], 'completed_chromosomes': []
    }
    
    # Group genes by chromosome
    logger.info("\nGrouping genes by chromosome...")
    chrom_groups = group_genes_by_chromosome(genes_df)
    
    # Filter out completed genes
    completed_set = set(progress['completed'])
    completed_chroms = set(progress.get('completed_chromosomes', []))
    
    # Determine which chromosomes to process
    # Only process standard chromosomes (1-22, X, Y) - skip alternate contigs
    valid_chroms = set([str(i) for i in range(1, 23)] + ['X', 'Y'])
    
    chromosomes_to_process = []
    skipped_contigs = []
    
    for chrom in sorted(chrom_groups.keys(), key=lambda x: (x not in ['X', 'Y'], x)):
        # Skip alternate contigs (NT_, NW_, etc.) - no gnomAD data
        if chrom not in valid_chroms:
            skipped_contigs.append(chrom)
            continue
            
        if chrom in completed_chroms:
            logger.info(f"Chromosome {chrom}: already completed (skipping)")
            continue
        
        chrom_genes = chrom_groups[chrom]
        remaining = chrom_genes[~chrom_genes['gene_symbol'].isin(completed_set)]
        
        if len(remaining) == 0:
            logger.info(f"Chromosome {chrom}: all genes completed")
            progress['completed_chromosomes'].append(chrom)
            continue
        
        chromosomes_to_process.append((chrom, remaining))
        logger.info(f"Chromosome {chrom}: {len(remaining)} genes to process")
    
    if skipped_contigs:
        logger.info(f"\nSkipped {len(skipped_contigs)} alternate contigs (no gnomAD data available):")
        logger.info(f"  {', '.join(skipped_contigs[:10])}" + 
                   (f"... and {len(skipped_contigs)-10} more" if len(skipped_contigs) > 10 else ""))
    
    if not chromosomes_to_process:
        logger.info("\n✅ All chromosomes already processed!")
        return progress
    
    logger.info(f"\nWill process {len(chromosomes_to_process)} chromosomes")
    logger.info(f"  Chromosome workers: {chrom_workers}")
    logger.info(f"  Gene workers per chromosome: {gene_workers}")
    
    # Process chromosomes
    start_time = time.time()
    
    if chrom_workers == 1:
        # Sequential chromosome processing (memory-safe)
        for chrom, genes in chromosomes_to_process:
            results = process_chromosome(
                chrom, genes, external_config, assembly, gene_workers
            )
            
            # Update progress
            for result in results:
                if result['status'] == 'success':
                    progress['completed'].append(result['gene'])
                else:
                    progress['failed'].append({
                        'gene': result['gene'],
                        'error': result['error'],
                        'timestamp': datetime.now().isoformat()
                    })
            
            # Mark chromosome as complete
            progress['completed_chromosomes'].append(chrom)
            save_progress(progress)
            
            logger.info(f"✅ Chromosome {chrom} complete")
    else:
        # Parallel chromosome processing (use with caution - memory intensive)
        with ProcessPoolExecutor(max_workers=chrom_workers) as executor:
            futures = {
                executor.submit(
                    process_chromosome,
                    chrom, genes, external_config, assembly, gene_workers
                ): chrom
                for chrom, genes in chromosomes_to_process
            }
            
            for future in as_completed(futures):
                chrom = futures[future]
                results = future.result()
                
                # Update progress
                for result in results:
                    if result['status'] == 'success':
                        progress['completed'].append(result['gene'])
                    else:
                        progress['failed'].append({
                            'gene': result['gene'],
                            'error': result['error'],
                            'timestamp': datetime.now().isoformat()
                        })
                
                progress['completed_chromosomes'].append(chrom)
                save_progress(progress)
    
    # Final summary
    elapsed = time.time() - start_time
    logger.info(f"\n{'='*80}")
    logger.info("QUERY COMPLETE")
    logger.info(f"{'='*80}")
    logger.info(f"  Successful: {len(progress['completed'])}")
    logger.info(f"  Failed: {len(progress['failed'])}")
    logger.info(f"  Total time: {elapsed/3600:.2f} hours")
    
    return progress


# ============================================================================
# COMBINING RESULTS
# ============================================================================

def process_all_missense_variants():
    """Combine all gene results into single dataset."""
    logger.info("\n" + "="*80)
    logger.info("PROCESSING MISSENSE VARIANTS")
    logger.info("="*80)
    
    gene_files = list(OUTPUT_DIR.glob("*_gnomad_variants.parquet"))
    logger.info(f"Found {len(gene_files)} gene result files")
    
    if not gene_files:
        logger.error("No gene result files found!")
        return
    
    all_missense = []
    
    for gene_file in gene_files:
        gene_symbol = gene_file.stem.replace("_gnomad_variants", "")
        
        try:
            df = pd.read_parquet(gene_file)
            missense = df[df['Consequence'].str.contains('missense', case=False, na=False)].copy()
            
            if len(missense) > 0:
                missense['gene_symbol'] = gene_symbol
                all_missense.append(missense)
                logger.info(f"  {gene_symbol}: {len(missense)} missense variants")
            
        except Exception as e:
            logger.error(f"  Failed to process {gene_symbol}: {e}")
            continue
    
    if not all_missense:
        logger.error("No missense variants found!")
        return
    
    logger.info("\nCombining all missense variants...")
    combined = pd.concat(all_missense, ignore_index=True)
    
    logger.info(f"Total missense variants: {len(combined):,}")
    logger.info(f"Unique genes: {combined['gene_symbol'].nunique():,}")
    
    # Save results
    output_parquet = OUTPUT_DIR / "all_gnomad_missense_variants.parquet"
    output_csv = OUTPUT_DIR / "all_gnomad_missense_variants.csv.gz"
    
    combined.to_parquet(output_parquet, index=False)
    combined.to_csv(output_csv, index=False, compression='gzip')
    
    logger.info(f"\n✅ Saved combined results:")
    logger.info(f"  {output_parquet} ({output_parquet.stat().st_size / 1e6:.1f} MB)")
    logger.info(f"  {output_csv} ({output_csv.stat().st_size / 1e6:.1f} MB)")
    
    return combined


# ============================================================================
# MAIN
# ============================================================================

def main():
    """Main execution."""
    from variant_utils.utils import read_external_config
    
    parser = argparse.ArgumentParser(
        description='Fetch all gnomAD missense variants (chromosome-batched)'
    )
    parser.add_argument('--chrom-workers', type=int, default=1,
                       help='Number of chromosomes to process in parallel (default: 1, recommend 1-2)')
    parser.add_argument('--gene-workers', type=int, default=32,
                       help='Number of genes to process in parallel per chromosome (default: 32, optimized for ProcessPool)')
    parser.add_argument('--no-resume', action='store_true',
                       help='Start fresh (ignore previous progress)')
    args = parser.parse_args()
    
    logger.info("\n" + "="*80)
    logger.info("FETCH ALL MISSENSE VARIANTS (CHROMOSOME-BATCHED)")
    logger.info("="*80)
    logger.info(f"\nAssembly: {ASSEMBLY}")
    logger.info(f"Output directory: {OUTPUT_DIR}")
    logger.info(f"Chromosome workers: {args.chrom_workers}")
    logger.info(f"Gene workers: {args.gene_workers}")
    logger.info(f"Resume: {not args.no_resume}")
    
    start_time = time.time()
    
    try:
        # Load external config
        external_config = read_external_config(EXTERNAL_TOOLS_CONFIG)
        
        # Step 1: Get all genes
        logger.info("\n" + "="*80)
        logger.info("STEP 1: LOADING GENES")
        logger.info("="*80)
        
        genes_df = get_all_genes_from_mane()
        
        # Step 2: Query by chromosome
        logger.info("\n" + "="*80)
        logger.info("STEP 2: QUERYING GNOMAD (CHROMOSOME-BATCHED)")
        logger.info("="*80)
        
        progress = query_all_genes_by_chromosome(
            genes_df=genes_df,
            external_config=external_config,
            assembly=ASSEMBLY,
            chrom_workers=args.chrom_workers,
            gene_workers=args.gene_workers,
            resume=not args.no_resume
        )
        
        # Step 3: Combine results
        logger.info("\n" + "="*80)
        logger.info("STEP 3: COMBINING RESULTS")
        logger.info("="*80)
        
        combined_df = process_all_missense_variants()
        
        # Final summary
        elapsed = time.time() - start_time
        logger.info("\n" + "="*80)
        logger.info("PIPELINE COMPLETE")
        logger.info("="*80)
        logger.info(f"\nTotal time: {elapsed/3600:.2f} hours")
        logger.info(f"Genes processed: {len(progress['completed'])}")
        logger.info(f"Genes failed: {len(progress['failed'])}")
        
        if combined_df is not None:
            logger.info(f"Total missense variants: {len(combined_df):,}")
        
    except Exception as e:
        logger.error(f"\n❌ PIPELINE FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    return True


if __name__ == "__main__":
    success = main()
    exit(0 if success else 1)
