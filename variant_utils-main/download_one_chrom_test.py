#!/usr/bin/env python3
"""
Lightweight single-chromosome test runner.
Use this before launching the full pipeline to verify bcftools + parsing work.

Example:
    python download_one_chrom_test.py --chrom 4 --gene-workers 2 --max-genes 10
"""

import argparse
import logging
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import List

import pandas as pd

from variant_utils.utils import read_external_config
from variant_utils.gnomad_utils import get_vep_columns_from_vcf_header
from download_all_gnomad_parallel import (
    get_all_genes_from_mane,
    extract_chromosome_vcf,
    process_gene_from_chrom_vcf,
    resolve_bcftools,
)

ASSEMBLY = "GRCh38"
EXTERNAL_TOOLS_CONFIG = "external_tools.json"

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(processName)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


def run_single_chrom(
    chrom: str,
    gene_workers: int,
    max_genes: int,
    bcftools_threads: int,
) -> List[dict]:
    external_config = read_external_config(EXTERNAL_TOOLS_CONFIG)
    bcftools_cmd = resolve_bcftools(external_config)

    logger.info(f"Assembly: {ASSEMBLY}")
    logger.info(f"Chromosome: {chrom}")
    logger.info(f"Gene workers: {gene_workers}")
    logger.info(f"Max genes: {max_genes if max_genes > 0 else 'all'}")
    logger.info(f"bcftools: {bcftools_cmd} (threads={bcftools_threads})")

    genes_df = get_all_genes_from_mane()
    genes_df = genes_df[genes_df["CHROM"].astype(str) == str(chrom)]

    if max_genes > 0:
        genes_df = genes_df.head(max_genes)

    if genes_df.empty:
        raise ValueError(f"No genes found for chromosome {chrom}")

    logger.info(f"Genes to process on chr{chrom}: {len(genes_df)}")

    exomes_vcf, genomes_vcf = extract_chromosome_vcf(
        ASSEMBLY,
        str(chrom),
        genes_df,
        external_config,
        bcftools_cmd,
        bcftools_threads,
    )

    vep_columns = get_vep_columns_from_vcf_header(str(exomes_vcf))
    logger.info(f"VEP columns extracted: {len(vep_columns)}")

    logger.info("Processing genes...")
    results = []
    with ProcessPoolExecutor(max_workers=gene_workers) as executor:
        futures = {
            executor.submit(
                process_gene_from_chrom_vcf,
                gene,
                exomes_vcf,
                genomes_vcf,
                vep_columns,
                ASSEMBLY,
            ): gene["gene_symbol"]
            for gene in genes_df.to_dict("records")
        }

        for future in as_completed(futures):
            result = future.result()
            results.append(result)
            status = "✅" if result["status"] == "success" else "❌"
            logger.info(f"  {status} {result['gene']}: {result['variants']:,} variants")

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Single-chromosome test run (bcftools + gene processing)"
    )
    parser.add_argument("--chrom", required=True, help="Chromosome to process (e.g., 4, X)")
    parser.add_argument(
        "--gene-workers",
        type=int,
        default=2,
        help="ProcessPool workers for genes (default: 2)",
    )
    parser.add_argument(
        "--max-genes",
        type=int,
        default=5,
        help="Limit number of genes for a quick test (<=0 means all).",
    )
    parser.add_argument(
        "--bcftools-threads",
        type=int,
        default=int(os.environ.get("BCFTOOLS_THREADS", "2")),
        help="Threads for bcftools view/index (default: env BCFTOOLS_THREADS or 2)",
    )
    args = parser.parse_args()

    results = run_single_chrom(
        chrom=args.chrom,
        gene_workers=args.gene_workers,
        max_genes=args.max_genes,
        bcftools_threads=max(1, args.bcftools_threads),
    )

    successes = sum(1 for r in results if r["status"] == "success")
    failures = [r for r in results if r["status"] != "success"]
    total_variants = sum(r.get("variants", 0) for r in results)

    logger.info("============================================================")
    logger.info("Single-chromosome test complete")
    logger.info(f"Genes processed: {len(results)} (success: {successes}, failed: {len(failures)})")
    logger.info(f"Total variants across processed genes: {total_variants:,}")
    if failures:
        logger.info("Failures:")
        for f in failures:
            logger.info(f"  {f['gene']}: {f['error']}")
    logger.info("============================================================")


if __name__ == "__main__":
    main()
