#!/usr/bin/env python3
"""
missense_pipeline_final.py

Memory-safe missense processing for gnomAD chromosome Parquets.

Modes:
  1) process-chrom  : stream one chr{chrom}_variants.parquet -> chr{chrom}_missense_enriched.parquet
  2) combine        : stream-combine per-chrom enriched files -> all/verified/mane/gold + statistics.json

Key properties:
- DOES NOT load an entire chromosome into RAM
- DOES NOT accumulate chromosomes in memory
- DOES NOT store full protein sequences per variant row
- Uses Arrow Dataset scanning + batch processing
"""

from __future__ import annotations

import argparse
import gzip
import json
import logging
import re
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.dataset as ds
import pyarrow.parquet as pq


# -----------------------------
# Logging
# -----------------------------
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger("missense_pipeline_final")


# -----------------------------
# Reference loaders
# -----------------------------
def load_enst_to_ensp(gtf_path: Path) -> Dict[str, str]:
    logger.info(f"   Parsing GTF: {gtf_path.name}")
    enst_to_ensp: Dict[str, str] = {}
    with gzip.open(gtf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue

            attr = {}
            for field in parts[8].split(";"):
                field = field.strip()
                if not field or " " not in field:
                    continue
                key, value = field.split(" ", 1)
                attr[key] = value.strip('"')

            tid = attr.get("transcript_id")
            pid = attr.get("protein_id")
            if tid and pid:
                enst_to_ensp[tid] = pid
                enst_to_ensp[tid.split(".")[0]] = pid

    logger.info(f"      Found {len(enst_to_ensp):,} ENST→ENSP mappings")
    return enst_to_ensp


def load_ensp_to_seq(pep_path: Path) -> Dict[str, str]:
    logger.info(f"   Parsing peptide FASTA: {pep_path.name}")
    ensp_to_seq: Dict[str, str] = {}
    current_id: Optional[str] = None
    seq_chunks: List[str] = []

    with gzip.open(pep_path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if current_id:
                    seq = "".join(seq_chunks)
                    ensp_to_seq[current_id] = seq
                    ensp_to_seq[current_id.split(".")[0]] = seq
                current_id = line[1:].strip().split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())

        if current_id:
            seq = "".join(seq_chunks)
            ensp_to_seq[current_id] = seq
            ensp_to_seq[current_id.split(".")[0]] = seq

    logger.info(f"      Found {len(ensp_to_seq):,} ENSP→sequence mappings")
    return ensp_to_seq


def build_enst_to_seq(gtf_path: Path, pep_path: Path) -> Dict[str, str]:
    enst_to_ensp = load_enst_to_ensp(gtf_path)
    ensp_to_seq = load_ensp_to_seq(pep_path)

    enst_to_seq: Dict[str, str] = {}
    for enst, ensp in enst_to_ensp.items():
        seq = ensp_to_seq.get(ensp) or ensp_to_seq.get(ensp.split(".")[0])
        if seq:
            enst_to_seq[enst] = seq
            enst_to_seq[enst.split(".")[0]] = seq

    logger.info(f"   ✅ Built {len(enst_to_seq):,} ENST→sequence mappings")
    return enst_to_seq


def load_mane_transcripts(mane_path: Path) -> Dict[str, Dict]:
    logger.info(f"   Parsing MANE summary: {mane_path.name}")
    mane_data: Dict[str, Dict] = {}

    with gzip.open(mane_path, "rt") as f:
        header = f.readline().strip().split("\t")
        col_idx = {name: i for i, name in enumerate(header)}

        def idx(*names: str) -> int:
            for n in names:
                if n in col_idx:
                    return col_idx[n]
            return -1

        i_status = idx("MANE_status", "MANE_Status")
        i_symbol = idx("symbol", "HGNC_symbol")
        i_enst = idx("Ensembl_nuc", "Ensembl_transcript")
        i_nm = idx("RefSeq_nuc", "RefSeq_transcript")

        for line in f:
            parts = line.strip().split("\t")
            if len(parts) != len(header):
                continue

            status = parts[i_status] if i_status >= 0 else ""
            if status not in ("MANE Select", "MANE_Select"):
                continue

            symbol = parts[i_symbol] if i_symbol >= 0 else ""
            enst = parts[i_enst] if i_enst >= 0 else ""
            nm = parts[i_nm] if i_nm >= 0 else ""

            if symbol and enst:
                mane_data[symbol] = {"enst": enst.split(".")[0], "refseq": nm or ""}

    logger.info(f"      Found {len(mane_data):,} MANE Select transcripts")
    return mane_data


def load_uniprot_id_mapping(idmapping_path: Path) -> Dict[str, str]:
    logger.info(f"   Parsing UniProt ID mapping: {idmapping_path.name}")
    gene_to_uniprot: Dict[str, str] = {}

    with gzip.open(idmapping_path, "rt") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) == 3 and parts[1] == "Gene_Name":
                if parts[2] not in gene_to_uniprot:
                    gene_to_uniprot[parts[2]] = parts[0]

    logger.info(f"      Found {len(gene_to_uniprot):,} gene→UniProt mappings")
    return gene_to_uniprot


# -----------------------------
# Fast HGVSp parsing (vectorized)
# -----------------------------
AA_3TO1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Ter": "*", "Stop": "*", "Sec": "U", "Pyl": "O",
}

def parse_hgvsp_cols(hgvsp: pd.Series) -> Tuple[pd.Series, pd.Series, pd.Series]:
    s = hgvsp.astype("string")
    # normalize "ENSP...:p." -> "p."
    s = s.str.replace(r"^.*?:p\.", "p.", regex=True)

    m3 = s.str.extract(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|\=)")
    ref3 = m3[0].map(AA_3TO1)
    pos3 = pd.to_numeric(m3[1], errors="coerce")
    alt3_raw = m3[2]
    alt3 = alt3_raw.map(AA_3TO1)
    alt3 = alt3.where(alt3_raw != "=", ref3)

    m1 = s.str.extract(r"p\.([A-Z\*])(\d+)([A-Z\*\=])")
    ref1 = m1[0]
    pos1 = pd.to_numeric(m1[1], errors="coerce")
    alt1 = m1[2].where(m1[2] != "=", ref1)

    ref = ref3.combine_first(ref1)
    pos = pos3.combine_first(pos1)
    alt = alt3.combine_first(alt1)
    return ref, pos, alt


# -----------------------------
# Fast seq verification (tight loop on valid subset)
# -----------------------------
def verify_seq_fast(feature: pd.Series, pos: pd.Series, ref: pd.Series, enst_to_seq: Dict[str, str]) -> pd.Series:
    feat_u = feature.astype("string").str.split(".", n=1).str[0]
    ref_s = ref.astype("string")

    valid = feat_u.notna() & (feat_u != "") & pos.notna() & ref_s.notna() & (ref_s != "")
    out = np.full(len(feature), pd.NA, dtype=object)

    idx = np.flatnonzero(valid.to_numpy())
    feats = feat_u.iloc[idx].to_numpy(dtype=object)
    poss = pos.iloc[idx].to_numpy()
    refs = ref_s.iloc[idx].to_numpy(dtype=object)

    for j in range(len(idx)):
        seq = enst_to_seq.get(feats[j])
        if not seq:
            continue
        p = int(poss[j]) - 1
        if p < 0 or p >= len(seq):
            out[idx[j]] = False
        else:
            out[idx[j]] = (seq[p] == refs[j])

    return pd.Series(out, index=feature.index, dtype="object")


# -----------------------------
# Streaming Parquet batch iterator (Arrow)
# -----------------------------
INPUT_COLS = [
    "CHROM", "POS", "REF", "ALT", "source",
    "Consequence", "SYMBOL", "Gene", "Feature", "HGVSc", "HGVSp",
    "CANONICAL", "MANE_SELECT", "MANE_PLUS_CLINICAL",
    "BIOTYPE", "IMPACT", "PolyPhen", "SIFT",
]

OUTPUT_EXTRA_COLS = [
    "ref_aa", "protein_pos", "alt_aa", "mutation",
    "seq_verified", "uniprot_id", "alphafold_id", "is_mane_select",
]

def iter_missense_batches(parquet_path: Path, batch_size: int):
    dataset = ds.dataset(str(parquet_path), format="parquet")

    # Filter at scan-time when possible (saves I/O)
    # Some files may not be string type; cast defensively.
    try:
        filt = pc.match_substring(pc.cast(ds.field("Consequence"), pa.string()), "missense_variant")
    except Exception:
        filt = None

    scanner = dataset.scanner(columns=INPUT_COLS, filter=filt, batch_size=batch_size)
    for rb in scanner.to_batches():
        table = pa.Table.from_batches([rb])
        df = table.to_pandas(types_mapper=pd.ArrowDtype)

        # Safety filter (Consequence can be comma-separated)
        mask = df["Consequence"].astype("string").str.contains("missense_variant", case=False, na=False)
        df = df.loc[mask]
        if len(df):
            yield df


def append_parquet(out_path: Path, df: pd.DataFrame, writer: Optional[pq.ParquetWriter]) -> pq.ParquetWriter:
    table = pa.Table.from_pandas(df, preserve_index=False)
    if writer is None:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        writer = pq.ParquetWriter(str(out_path), table.schema, compression="snappy")
    writer.write_table(table)
    return writer


# -----------------------------
# Process one chromosome
# -----------------------------
def process_chrom(args) -> int:
    chrom = args.chrom
    input_path = args.input_dir / f"chr{chrom}_variants.parquet"
    if not input_path.exists():
        logger.error(f"Missing input parquet: {input_path}")
        return 2

    out_dir = args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    out_parquet = out_dir / f"chr{chrom}_missense_enriched.parquet"
    out_stats = out_dir / f"chr{chrom}_stats.json"

    t0 = time.time()

    # Load references once per job
    ref_dir = args.reference_dir
    logger.info("Loading reference caches...")
    enst_to_seq = build_enst_to_seq(
        ref_dir / "ensembl" / "Homo_sapiens.GRCh38.112.gtf.gz",
        ref_dir / "ensembl" / "Homo_sapiens.GRCh38.pep.all.fa.gz",
    )
    mane = load_mane_transcripts(ref_dir / "mane" / "MANE.GRCh38.v1.3.summary.txt.gz")
    gene_to_uniprot = load_uniprot_id_mapping(ref_dir / "uniprot" / "HUMAN_9606_idmapping.dat.gz")

    writer: Optional[pq.ParquetWriter] = None

    total = 0
    with_seq_info = 0
    verified_true = 0
    with_uniprot = 0
    mane_true = 0

    logger.info(f"chr{chrom}: streaming batches from {input_path.name}")
    for df in iter_missense_batches(input_path, batch_size=args.batch_size):
        total += len(df)

        # HGVSp parse
        ref_aa, prot_pos, alt_aa = parse_hgvsp_cols(df["HGVSp"])
        df["ref_aa"] = ref_aa
        df["protein_pos"] = prot_pos
        df["alt_aa"] = alt_aa
        df["mutation"] = (
            df["ref_aa"].astype("string")
            + df["protein_pos"].astype("Int64").astype("string")
            + df["alt_aa"].astype("string")
        )

        # Sequence verification (no sequences stored)
        df["seq_verified"] = verify_seq_fast(df["Feature"], df["protein_pos"], df["ref_aa"], enst_to_seq)
        with_seq_info += int(df["seq_verified"].notna().sum())
        verified_true += int((df["seq_verified"] == True).sum())

        # UniProt + AF ids
        df["uniprot_id"] = df["SYMBOL"].astype("string").map(gene_to_uniprot)
        df["alphafold_id"] = df["uniprot_id"].apply(lambda x: f"AF-{x}-F1" if pd.notna(x) and x else pd.NA)
        with_uniprot += int(df["uniprot_id"].notna().sum())

        # MANE select flag
        feat_u = df["Feature"].astype("string").str.split(".", n=1).str[0]
        mane_enst = df["SYMBOL"].astype("string").map(lambda s: mane.get(s, {}).get("enst") if pd.notna(s) else None)
        df["is_mane_select"] = mane_enst.notna() & (feat_u == mane_enst.astype("string"))
        mane_true += int(df["is_mane_select"].sum())

        # Write only necessary columns
        keep_cols = [c for c in INPUT_COLS + OUTPUT_EXTRA_COLS if c in df.columns]
        df_out = df[keep_cols].copy()

        writer = append_parquet(out_parquet, df_out, writer)

        # Encourage early release
        del df, df_out

    if writer is not None:
        writer.close()

    stats = {
        "chrom": chrom,
        "total_missense": total,
        "with_seq_info": with_seq_info,
        "verified_true": verified_true,
        "with_uniprot": with_uniprot,
        "mane_select_true": mane_true,
        "minutes": round((time.time() - t0) / 60, 2),
        "out_file": str(out_parquet),
    }
    with open(out_stats, "w") as f:
        json.dump(stats, f, indent=2)

    logger.info(f"chr{chrom}: ✅ done missense={total:,} time={stats['minutes']} min")
    return 0


# -----------------------------
# Combine per-chrom outputs (memory-safe streaming)
# -----------------------------
def stream_write_dataset(input_files: List[Path], out_path: Path, kind: str, batch_size: int) -> int:
    dataset = ds.dataset([str(p) for p in input_files], format="parquet")

    f_verified = (ds.field("seq_verified") == True)
    f_mane = (ds.field("is_mane_select") == True)

    if kind == "all":
        filt = None
    elif kind == "verified":
        filt = f_verified
    elif kind == "mane":
        filt = f_mane
    elif kind == "gold":
        filt = f_verified & f_mane
    else:
        raise ValueError(kind)

    scanner = dataset.scanner(filter=filt, batch_size=batch_size)

    writer: Optional[pq.ParquetWriter] = None
    n = 0
    try:
        for rb in scanner.to_batches():
            table = pa.Table.from_batches([rb])
            if table.num_rows == 0:
                continue
            if writer is None:
                out_path.parent.mkdir(parents=True, exist_ok=True)
                writer = pq.ParquetWriter(str(out_path), table.schema, compression="snappy")
            writer.write_table(table)
            n += table.num_rows
    finally:
        if writer is not None:
            writer.close()
    return n


def combine(args) -> int:
    out_dir = args.output_dir
    inputs = sorted(out_dir.glob("chr*_missense_enriched.parquet"))
    if not inputs:
        logger.error(f"No per-chrom files found in {out_dir}")
        return 2

    t0 = time.time()
    out_all = out_dir / "all_missense_variants.parquet"
    out_verified = out_dir / "verified_missense_variants.parquet"
    out_mane = out_dir / "mane_select_missense_variants.parquet"
    out_gold = out_dir / "gold_standard_missense_variants.parquet"

    n_all = stream_write_dataset(inputs, out_all, "all", args.batch_size)
    n_ver = stream_write_dataset(inputs, out_verified, "verified", args.batch_size)
    n_man = stream_write_dataset(inputs, out_mane, "mane", args.batch_size)
    n_gol = stream_write_dataset(inputs, out_gold, "gold", args.batch_size)

    stats = {
        "per_chrom_files": len(inputs),
        "combined_rows": {"all": n_all, "verified": n_ver, "mane": n_man, "gold": n_gol},
        "minutes": round((time.time() - t0) / 60, 2),
    }
    with open(out_dir / "statistics.json", "w") as f:
        json.dump(stats, f, indent=2)

    logger.info(f"✅ combine done in {stats['minutes']} min (all={n_all:,})")
    return 0


# -----------------------------
# CLI
# -----------------------------
def main() -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("process-chrom")
    p.add_argument("--chrom", required=True, help="1..22, X, Y")
    p.add_argument("--input_dir", type=Path, required=True)
    p.add_argument("--output_dir", type=Path, required=True)
    p.add_argument("--reference_dir", type=Path, required=True)
    p.add_argument("--batch_size", type=int, default=200_000)

    c = sub.add_parser("combine")
    c.add_argument("--output_dir", type=Path, required=True)
    c.add_argument("--batch_size", type=int, default=500_000)

    args = ap.parse_args()

    if args.cmd == "process-chrom":
        return process_chrom(args)
    else:
        return combine(args)


if __name__ == "__main__":
    raise SystemExit(main())
