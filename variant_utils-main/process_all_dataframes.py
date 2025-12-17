#!/usr/bin/env python3
"""
missense_pipeline_robust_v2.py

Fixes schema mismatch across batches by writing seq_verified as int8:
  1=True, 0=False, -1=NA

Everything else remains streaming + schema-robust.
"""

from __future__ import annotations

import argparse
import gzip
import json
import logging
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.dataset as ds
import pyarrow.parquet as pq


logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger("missense_pipeline_robust_v2")


# =========================
# Reference loaders
# =========================
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


# =========================
# HGVSp parsing + seq verify
# =========================
AA_3TO1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Ter": "*", "Stop": "*", "Sec": "U", "Pyl": "O",
}

def parse_hgvsp_cols(hgvsp: pd.Series) -> Tuple[pd.Series, pd.Series, pd.Series]:
    s = hgvsp.astype("string")
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


def seq_verified_to_int8(s: pd.Series) -> pd.Series:
    # 1=True, 0=False, -1=NA
    out = pd.Series(np.full(len(s), -1, dtype=np.int8), index=s.index)
    out.loc[s == True] = 1
    out.loc[s == False] = 0
    return out


# =========================
# Streaming (schema-robust)
# =========================
CORE_COLS = [
    "CHROM", "POS", "REF", "ALT",
    "Consequence", "SYMBOL", "Gene", "Feature", "HGVSc", "HGVSp",
    "CANONICAL", "MANE_SELECT", "MANE_PLUS_CLINICAL", "BIOTYPE", "IMPACT",
]
OPTIONAL_COLS = ["source", "FILTER", "SIFT", "PolyPhen"]

OUTPUT_EXTRA = [
    "ref_aa", "protein_pos", "alt_aa", "mutation",
    "seq_verified_i8", "uniprot_id", "alphafold_id", "is_mane_select",
]

def available_columns(parquet_path: Path) -> List[str]:
    return pq.ParquetFile(str(parquet_path)).schema.names


def iter_missense_batches(parquet_path: Path, batch_size: int):
    cols_in_file = set(available_columns(parquet_path))
    needed = [c for c in (CORE_COLS + OPTIONAL_COLS) if c in cols_in_file]

    min_required = {"Consequence", "SYMBOL", "Feature", "HGVSp"}
    if not min_required.issubset(cols_in_file):
        missing = sorted(list(min_required - cols_in_file))
        raise RuntimeError(f"{parquet_path.name} missing required columns: {missing}")

    dataset = ds.dataset(str(parquet_path), format="parquet")

    try:
        filt = pc.match_substring(pc.cast(ds.field("Consequence"), pa.string()), "missense_variant")
    except Exception:
        filt = None

    scanner = dataset.scanner(columns=needed, filter=filt, batch_size=batch_size)

    for rb in scanner.to_batches():
        table = pa.Table.from_batches([rb])
        df = table.to_pandas(types_mapper=pd.ArrowDtype)

        mask = df["Consequence"].astype("string").str.contains("missense_variant", case=False, na=False)
        df = df.loc[mask]
        if len(df) == 0:
            continue

        for c in (CORE_COLS + OPTIONAL_COLS):
            if c not in df.columns:
                df[c] = pd.NA

        yield df


# =========================
# Process one chromosome
# =========================
def process_chrom(chrom: str, input_dir: Path, output_dir: Path, reference_dir: Path, batch_size: int) -> int:
    in_path = input_dir / f"chr{chrom}_variants.parquet"
    if not in_path.exists():
        logger.error(f"Missing: {in_path}")
        return 2

    output_dir.mkdir(parents=True, exist_ok=True)
    out_parquet = output_dir / f"chr{chrom}_missense_enriched.parquet"
    out_stats = output_dir / f"chr{chrom}_stats.json"

    t0 = time.time()

    logger.info("Loading reference caches...")
    enst_to_seq = build_enst_to_seq(
        reference_dir / "ensembl" / "Homo_sapiens.GRCh38.112.gtf.gz",
        reference_dir / "ensembl" / "Homo_sapiens.GRCh38.pep.all.fa.gz",
    )
    mane = load_mane_transcripts(reference_dir / "mane" / "MANE.GRCh38.v1.3.summary.txt.gz")
    gene_to_uniprot = load_uniprot_id_mapping(reference_dir / "uniprot" / "HUMAN_9606_idmapping.dat.gz")

    writer: Optional[pq.ParquetWriter] = None

    total = 0
    verified_true = 0
    verified_false = 0
    verified_na = 0
    mane_true = 0
    with_uniprot = 0

    logger.info(f"chr{chrom}: streaming batches from {in_path.name}")
    for df in iter_missense_batches(in_path, batch_size=batch_size):
        total += len(df)

        ref_aa, prot_pos, alt_aa = parse_hgvsp_cols(df["HGVSp"])
        df["ref_aa"] = ref_aa
        df["protein_pos"] = prot_pos.astype("Int64")
        df["alt_aa"] = alt_aa
        df["mutation"] = (
            df["ref_aa"].astype("string")
            + df["protein_pos"].astype("string")
            + df["alt_aa"].astype("string")
        )

        sv = verify_seq_fast(df["Feature"], df["protein_pos"], df["ref_aa"], enst_to_seq)
        df["seq_verified_i8"] = seq_verified_to_int8(sv)

        verified_true += int((df["seq_verified_i8"] == 1).sum())
        verified_false += int((df["seq_verified_i8"] == 0).sum())
        verified_na += int((df["seq_verified_i8"] == -1).sum())

        df["uniprot_id"] = df["SYMBOL"].astype("string").map(gene_to_uniprot)
        df["alphafold_id"] = df["uniprot_id"].apply(lambda x: f"AF-{x}-F1" if pd.notna(x) and x else pd.NA)
        with_uniprot += int(df["uniprot_id"].notna().sum())

        feat_u = df["Feature"].astype("string").str.split(".", n=1).str[0]
        mane_enst = df["SYMBOL"].astype("string").map(lambda s: mane.get(s, {}).get("enst") if pd.notna(s) else None)
        df["is_mane_select"] = mane_enst.notna() & (feat_u == mane_enst.astype("string"))
        mane_true += int(df["is_mane_select"].sum())

        keep = [c for c in (CORE_COLS + OPTIONAL_COLS + OUTPUT_EXTRA) if c in df.columns]
        df_out = df[keep].copy()

        table = pa.Table.from_pandas(df_out, preserve_index=False)
        if writer is None:
            writer = pq.ParquetWriter(str(out_parquet), table.schema, compression="snappy")
        writer.write_table(table)

        del df, df_out, table

    if writer is not None:
        writer.close()

    stats = {
        "chrom": chrom,
        "total_missense": total,
        "seq_verified_true": verified_true,
        "seq_verified_false": verified_false,
        "seq_verified_na": verified_na,
        "mane_select_true": mane_true,
        "with_uniprot": with_uniprot,
        "minutes": round((time.time() - t0) / 60, 2),
        "out_file": str(out_parquet),
    }
    with open(out_stats, "w") as f:
        json.dump(stats, f, indent=2)

    logger.info(f"chr{chrom}: ✅ done missense={total:,} time={stats['minutes']} min -> {out_parquet.name}")
    return 0


# =========================
# Combine step (streaming)
# =========================
def stream_write_dataset(input_files: List[Path], out_path: Path, kind: str, batch_size: int) -> int:
    dataset = ds.dataset([str(p) for p in input_files], format="parquet")
    f_verified = (ds.field("seq_verified_i8") == 1)
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


def combine(output_dir: Path, batch_size: int) -> int:
    files = sorted(output_dir.glob("chr*_missense_enriched.parquet"))
    if not files:
        logger.error(f"No per-chrom outputs found in {output_dir}")
        return 2

    t0 = time.time()
    out_all = output_dir / "all_missense_variants.parquet"
    out_ver = output_dir / "verified_missense_variants.parquet"
    out_mane = output_dir / "mane_select_missense_variants.parquet"
    out_gold = output_dir / "gold_standard_missense_variants.parquet"

    n_all = stream_write_dataset(files, out_all, "all", batch_size)
    n_ver = stream_write_dataset(files, out_ver, "verified", batch_size)
    n_man = stream_write_dataset(files, out_mane, "mane", batch_size)
    n_gol = stream_write_dataset(files, out_gold, "gold", batch_size)

    stats = {
        "per_chrom_files": len(files),
        "combined_rows": {"all": n_all, "verified": n_ver, "mane": n_man, "gold": n_gol},
        "minutes": round((time.time() - t0) / 60, 2),
    }
    with open(output_dir / "statistics.json", "w") as f:
        json.dump(stats, f, indent=2)

    logger.info(f"✅ combine done in {stats['minutes']} min (all={n_all:,})")
    return 0


# =========================
# CLI
# =========================
def main() -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("process-chrom")
    p.add_argument("--chrom", required=True)
    p.add_argument("--input_dir", type=Path, required=True)
    p.add_argument("--output_dir", type=Path, required=True)
    p.add_argument("--reference_dir", type=Path, required=True)
    p.add_argument("--batch_size", type=int, default=200_000)

    c = sub.add_parser("combine")
    c.add_argument("--output_dir", type=Path, required=True)
    c.add_argument("--batch_size", type=int, default=500_000)

    args = ap.parse_args()

    if args.cmd == "process-chrom":
        return process_chrom(args.chrom, args.input_dir, args.output_dir, args.reference_dir, args.batch_size)
    else:
        return combine(args.output_dir, args.batch_size)


if __name__ == "__main__":
    raise SystemExit(main())
