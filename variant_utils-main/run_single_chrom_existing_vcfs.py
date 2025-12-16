#!/usr/bin/env python3
"""
Single-chromosome runner that uses existing per-chromosome exome/genome VCFs
without re-extracting or re-indexing. Useful for quick validation runs when
VCFs are already split by chromosome and indexed.

Example:
    python run_single_chrom_existing_vcfs.py --chrom 22 --gene-workers 4 --max-genes 10
"""

import argparse
import logging
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List

import pandas as pd

from variant_utils.utils import read_external_config
from variant_utils.gnomad_utils import (
    get_vep_columns_from_vcf_header,
    extract_variants_fast,
    merge_exome_genome_dataframes,
    parse_vep,
)
from download_all_gnomad_parallel import get_all_genes_from_mane

ASSEMBLY = "GRCh38"
OUTPUT_DIR = Path("gnomad_single_chrom_test")
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(processName)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


def build_vcf_paths(external_config: Dict, chrom: str) -> Dict[str, Path]:
    release_version = "v4.1" if ASSEMBLY == "GRCh38" else "r2.1.1"
    chr_prefix = "chr" if release_version == "v4.1" else ""
    root = Path(
        external_config["gnomad_v4_vcf_root"]
        if release_version == "v4.1"
        else external_config["gnomad_v2_vcf_root"]
    )
    exomes = root / f"exomes/gnomad.exomes.{release_version}.sites.{chr_prefix}{chrom}.vcf.bgz"
    genomes = root / f"genomes/gnomad.genomes.{release_version}.sites.{chr_prefix}{chrom}.vcf.bgz"
    return {"exomes": exomes, "genomes": genomes, "chr_prefix": chr_prefix}


def process_gene(
    gene_row: Dict,
    exomes_vcf: Path,
    genomes_vcf: Path,
    vep_columns: List[str],
    chr_prefix: str,
) -> Dict:
    gene_symbol = gene_row["gene_symbol"]
    chrom = gene_row["chrom"]
    start = int(gene_row["gene_start"])
    end = int(gene_row["gene_end"])
    try:
        # Extract exome/genome slices in parallel (I/O bound)
        with ThreadPoolExecutor(max_workers=2) as pool:
            exomes_future = pool.submit(
                extract_variants_fast, exomes_vcf, chrom, start, end, chr_prefix
            )
            genomes_future = pool.submit(
                extract_variants_fast, genomes_vcf, chrom, start, end, chr_prefix
            )
            exomes_df = exomes_future.result()
            genomes_df = genomes_future.result()

        df = merge_exome_genome_dataframes(exomes_df, genomes_df)
        if df.empty:
            return {"gene": gene_symbol, "status": "success", "variants": 0, "error": None}

        vep_df = parse_vep(df, columns=vep_columns)
        df = pd.merge(df, vep_df, left_index=True, right_on="index", validate="one_to_many")
        gene_df = df[df.HGNC_ID == gene_row["hgnc_id"]].copy()
        gene_df = gene_df.assign(
            CHROM=gene_df.CHROM.astype(str).str.replace("chr", ""),
            POS=gene_df.POS.astype(str),
            REF=gene_df.REF.astype(str),
            ALT=gene_df.ALT.astype(str),
        )
        out_file = OUTPUT_DIR / f"{gene_symbol}_gnomad_variants.parquet"
        gene_df.to_parquet(out_file, index=False)

        return {"gene": gene_symbol, "status": "success", "variants": len(gene_df), "error": None}
    except Exception as exc:  # noqa: BLE001
        return {"gene": gene_symbol, "status": "failed", "variants": 0, "error": str(exc)}


def main():
    parser = argparse.ArgumentParser(
        description="Process one chromosome using existing per-chromosome VCFs (no extraction)."
    )
    parser.add_argument("--chrom", required=True, help="Chromosome (e.g., 22, X)")
    parser.add_argument("--gene-workers", type=int, default=4, help="ProcessPool workers (default: 4)")
    parser.add_argument(
        "--max-genes",
        type=int,
        default=0,
        help="Limit number of genes for a quick test (<=0 means all genes on the chromosome).",
    )
    args = parser.parse_args()

    external_config = read_external_config("external_tools.json")
    vcf_info = build_vcf_paths(external_config, args.chrom)

    for label, path in [("exomes", vcf_info["exomes"]), ("genomes", vcf_info["genomes"])]:
        if not path.exists():
            raise FileNotFoundError(f"{label} VCF missing: {path}")
        if not Path(str(path) + ".tbi").exists() and not Path(str(path) + ".csi").exists():
            raise FileNotFoundError(
                f"Index missing for {label} VCF: {path}. Please run tabix/bcftools index once."
            )

    logger.info(f"Assembly: {ASSEMBLY}")
    logger.info(f"Chromosome: {args.chrom}")
    logger.info(f"Gene workers: {args.gene_workers}")
    logger.info(f"Max genes: {args.max_genes if args.max_genes > 0 else 'all'}")
    logger.info(f"Exomes VCF: {vcf_info['exomes']}")
    logger.info(f"Genomes VCF: {vcf_info['genomes']}")

    genes_df = get_all_genes_from_mane()
    genes_df = genes_df[genes_df["chrom"].astype(str) == str(args.chrom)]
    if args.max_genes > 0:
        genes_df = genes_df.head(args.max_genes)
    if genes_df.empty:
        raise ValueError(f"No genes found for chromosome {args.chrom}")

    vep_columns = get_vep_columns_from_vcf_header(str(vcf_info["exomes"]))
    logger.info(f"VEP columns: {len(vep_columns)}")
    logger.info(f"Genes to process: {len(genes_df)}")

    results: List[Dict] = []
    with ProcessPoolExecutor(max_workers=args.gene_workers) as executor:
        futures = {
            executor.submit(
                process_gene,
                gene,
                vcf_info["exomes"],
                vcf_info["genomes"],
                vep_columns,
                vcf_info["chr_prefix"],
            ): gene["gene_symbol"]
            for gene in genes_df.to_dict("records")
        }
        for future in as_completed(futures):
            result = future.result()
            results.append(result)
            icon = "✅" if result["status"] == "success" else "❌"
            logger.info(f"  {icon} {result['gene']}: {result['variants']:,} variants")

    successes = sum(1 for r in results if r["status"] == "success")
    failures = [r for r in results if r["status"] != "success"]
    total_variants = sum(r.get("variants", 0) for r in results)

    combined = OUTPUT_DIR / f"chr{args.chrom}_combined_results.parquet"
    # Concatenate successful gene outputs
    frames = []
    for r in results:
        if r["status"] != "success":
            continue
        fpath = OUTPUT_DIR / f"{r['gene']}_gnomad_variants.parquet"
        if fpath.exists():
            frames.append(pd.read_parquet(fpath))
    if frames:
        pd.concat(frames, ignore_index=True).to_parquet(combined, index=False)

    logger.info("============================================================")
    logger.info("Single-chromosome existing-VCF run complete")
    logger.info(f"Genes processed: {len(results)} (success: {successes}, failed: {len(failures)})")
    logger.info(f"Total variants across processed genes: {total_variants:,}")
    logger.info(f"Combined Parquet (if any): {combined}")
    if failures:
        logger.info("Failures:")
        for f in failures:
            logger.info(f"  {f['gene']}: {f['error']}")
    logger.info("============================================================")


if __name__ == "__main__":
    main()
