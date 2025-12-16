#!/usr/bin/env python3
"""
FAST: Stream gnomAD VCF -> Parquet per chromosome (no temp files, no pandas concat)

Key properties:
- Processes VCF sequentially
- Expands VEP annotations: 1 row per transcript annotation
- Writes in large chunks directly to a single Parquet file
- Memory is bounded by --chunk-rows
- Skips chromosome if output already exists (unless --force)

Example:
  python convert_vcf_to_dataframe_fast.py --chrom 1 --chunk-rows 500000 --threads 8
"""

from __future__ import annotations

import argparse
import logging
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pyarrow as pa
import pyarrow.parquet as pq
from pysam import VariantFile


# ----------------------------
# CONFIG: paths
# ----------------------------
SOURCE_EXOMES = Path("/projects/lugoteam/protein_graphs/GnomAD-Parser/data/exomes")
SOURCE_GENOMES = Path("/projects/lugoteam/protein_graphs/GnomAD-Parser/data/genomes")
OUT_DIR = Path("gnomad_all_genes/chromosome_dataframes")

DEFAULT_CHUNK_ROWS = 500_000  # rows (after VEP expansion)
DEFAULT_THREADS = 8

# INFO fields to pull (some may be missing depending on file)
INFO_FIELDS_INT = ["AC", "AN", "nhomalt"]
INFO_FIELDS_FLOAT = ["AF", "faf95", "faf99", "cadd_phred", "phylop"]
INFO_FIELDS_POP_AC = ["AC_afr", "AC_amr", "AC_asj", "AC_eas", "AC_fin", "AC_nfe", "AC_sas"]


# ----------------------------
# logging
# ----------------------------
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
log = logging.getLogger("gnomad-fast")


def get_vep_columns_from_header(vcf_path: Path, threads: int) -> List[str]:
    vcf = VariantFile(str(vcf_path), threads=threads)
    try:
        for record in vcf.header.records:
            if record.key == "INFO" and "ID" in record and record["ID"] == "vep":
                desc = record.get("Description", "")
                if "Format:" in desc:
                    fmt = desc.split("Format:", 1)[1].strip()
                    cols = [c.strip() for c in fmt.split("|")]
                    if cols and cols[0] != "":
                        return cols
        raise RuntimeError("Could not find INFO/vep Format: ... in VCF header.")
    finally:
        vcf.close()


def _scalar_info(v, default=None):
    """Convert pysam INFO value into scalar where possible."""
    if v is None:
        return default
    if isinstance(v, tuple):
        if len(v) == 0:
            return default
        # If allele-specific with 1 alt allele, often len=1. Keep first.
        if len(v) == 1:
            return v[0]
        # Otherwise keep as string to avoid type chaos
        return ",".join(map(str, v))
    return v


def build_schema(vep_cols: List[str]) -> pa.Schema:
    fields = [
        pa.field("CHROM", pa.string()),
        pa.field("POS", pa.int32()),
        pa.field("REF", pa.string()),
        pa.field("ALT", pa.string()),
        pa.field("source", pa.string()),
    ]

    for c in INFO_FIELDS_INT + INFO_FIELDS_POP_AC:
        fields.append(pa.field(c, pa.int64()))
    for c in INFO_FIELDS_FLOAT:
        fields.append(pa.field(c, pa.float64()))

    # VEP columns as strings (safe; you can cast later downstream if desired)
    for c in vep_cols:
        fields.append(pa.field(c, pa.string()))

    return pa.schema(fields)


def make_empty_buffers(schema: pa.Schema) -> Dict[str, list]:
    return {f.name: [] for f in schema}


def append_row(
    buffers: Dict[str, list],
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    source: str,
    info_vals: Dict[str, Optional[object]],
    vep_map: Optional[Dict[str, str]],
    vep_cols: List[str],
):
    buffers["CHROM"].append(chrom)
    buffers["POS"].append(pos)
    buffers["REF"].append(ref)
    buffers["ALT"].append(alt)
    buffers["source"].append(source)

    for c in INFO_FIELDS_INT + INFO_FIELDS_POP_AC:
        buffers[c].append(info_vals.get(c))
    for c in INFO_FIELDS_FLOAT:
        buffers[c].append(info_vals.get(c))

    if vep_map is None:
        for c in vep_cols:
            buffers[c].append("")
    else:
        for c in vep_cols:
            buffers[c].append(vep_map.get(c, ""))


def flush(writer: pq.ParquetWriter, buffers: Dict[str, list], schema: pa.Schema) -> int:
    if not buffers["POS"]:
        return 0
    table = pa.Table.from_pydict(buffers, schema=schema)
    writer.write_table(table)
    n = table.num_rows
    for k in buffers:
        buffers[k].clear()
    return n


def stream_one_vcf(
    vcf_path: Path,
    writer: pq.ParquetWriter,
    schema: pa.Schema,
    vep_cols: List[str],
    source_label: str,
    chunk_rows: int,
    threads: int,
) -> Tuple[int, int]:
    """
    Returns: (variants_total_seen, rows_written_after_vep_expansion)
    """
    vcf = VariantFile(str(vcf_path), threads=threads)
    buffers = make_empty_buffers(schema)

    variants_seen = 0
    rows_written = 0
    rows_in_buffer = 0
    last_log = time.time()

    try:
        for rec in vcf:
            variants_seen += 1

            # periodic progress
            if variants_seen % 2_000_000 == 0 or (time.time() - last_log) > 120:
                log.info(f"    {vcf_path.name}: seen {variants_seen:,} variants, wrote {rows_written:,} rows")
                last_log = time.time()

            # Require biallelic SNP
            if not rec.alts or len(rec.alts) != 1:
                continue
            if len(rec.ref) != 1 or len(rec.alts[0]) != 1:
                continue

            # PASS only
            if rec.filter.keys() and "PASS" not in rec.filter.keys():
                continue

            chrom = rec.chrom.replace("chr", "")
            pos = int(rec.pos)
            ref = rec.ref
            alt = str(rec.alts[0])

            info_vals: Dict[str, Optional[object]] = {}
            for c in INFO_FIELDS_INT + INFO_FIELDS_POP_AC:
                val = _scalar_info(rec.info.get(c))
                if isinstance(val, str):
                    # if it became CSV string, keep None to avoid type issues in int column
                    info_vals[c] = None
                else:
                    info_vals[c] = int(val) if val is not None else None

            for c in INFO_FIELDS_FLOAT:
                val = _scalar_info(rec.info.get(c))
                if isinstance(val, str):
                    info_vals[c] = None
                else:
                    info_vals[c] = float(val) if val is not None else None

            vep_str = rec.info.get("vep")
            if vep_str:
                # Expand 1 row per transcript annotation
                for ann in str(vep_str).split(","):
                    parts = ann.split("|")
                    vmap = {vep_cols[i]: (parts[i] if i < len(parts) else "") for i in range(len(vep_cols))}
                    append_row(buffers, chrom, pos, ref, alt, source_label, info_vals, vmap, vep_cols)
                    rows_in_buffer += 1
            else:
                append_row(buffers, chrom, pos, ref, alt, source_label, info_vals, None, vep_cols)
                rows_in_buffer += 1

            if rows_in_buffer >= chunk_rows:
                rows_written += flush(writer, buffers, schema)
                rows_in_buffer = 0

        # final flush
        rows_written += flush(writer, buffers, schema)

    finally:
        vcf.close()

    return variants_seen, rows_written


def process_chromosome(
    chrom: str,
    chunk_rows: int,
    threads: int,
    force: bool,
) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_file = OUT_DIR / f"chr{chrom}_variants.parquet"

    if out_file.exists() and not force:
        # Fast skip without reading whole file
        pf = pq.ParquetFile(out_file)
        log.info(f"✅ chr{chrom} exists, skipping (rows={pf.metadata.num_rows:,})")
        return

    exomes = SOURCE_EXOMES / f"gnomad.exomes.v4.1.sites.chr{chrom}.vcf.bgz"
    genomes = SOURCE_GENOMES / f"gnomad.genomes.v4.1.sites.chr{chrom}.vcf.bgz"

    if not exomes.exists():
        raise FileNotFoundError(f"Missing exomes VCF: {exomes}")
    if not genomes.exists():
        raise FileNotFoundError(f"Missing genomes VCF: {genomes}")

    log.info(f"{'='*80}")
    log.info(f"Processing chr{chrom}")
    log.info(f"  exomes:  {exomes.name}")
    log.info(f"  genomes: {genomes.name}")
    log.info(f"  chunk_rows={chunk_rows:,} threads={threads}")
    log.info(f"{'='*80}")

    # Use exomes header to get VEP columns
    vep_cols = get_vep_columns_from_header(exomes, threads=threads)
    log.info(f"VEP fields: {len(vep_cols)}")

    schema = build_schema(vep_cols)

    # Write to a temp name then move into place (atomic-ish)
    tmp_file = OUT_DIR / f".tmp_chr{chrom}_{int(time.time())}.parquet"
    start = time.time()

    writer = pq.ParquetWriter(tmp_file, schema=schema, compression="snappy")
    try:
        seen_e, rows_e = stream_one_vcf(exomes, writer, schema, vep_cols, "exomes", chunk_rows, threads)
        log.info(f"  exomes done: seen={seen_e:,} variants, wrote={rows_e:,} rows")

        seen_g, rows_g = stream_one_vcf(genomes, writer, schema, vep_cols, "genomes", chunk_rows, threads)
        log.info(f"  genomes done: seen={seen_g:,} variants, wrote={rows_g:,} rows")
    finally:
        writer.close()

    tmp_file.replace(out_file)

    elapsed = time.time() - start
    pf = pq.ParquetFile(out_file)
    size_gb = out_file.stat().st_size / (1024**3)

    log.info(f"✅ chr{chrom} complete: rows={pf.metadata.num_rows:,} size={size_gb:.2f} GB time={elapsed/60:.1f} min")


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", required=True, help="Chromosome: 1-22, X, Y")
    ap.add_argument("--chunk-rows", type=int, default=DEFAULT_CHUNK_ROWS, help="Rows per write chunk")
    ap.add_argument("--threads", type=int, default=DEFAULT_THREADS, help="Threads for bgzip decode")
    ap.add_argument("--force", action="store_true", help="Recompute even if output exists")
    return ap.parse_args()


def main():
    args = parse_args()
    chrom = str(args.chrom).replace("chr", "")
    if chrom not in [str(i) for i in range(1, 23)] + ["X", "Y"]:
        raise ValueError(f"Invalid chrom: {chrom}")

    process_chromosome(chrom, args.chunk_rows, args.threads, args.force)


if __name__ == "__main__":
    main()
