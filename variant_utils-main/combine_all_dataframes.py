#!/usr/bin/env python3
"""
Robust combiner for per-chrom missense enriched parquet outputs.

Fixes schema mismatch across chr files (e.g., FILTER null vs string) by
dropping problematic columns during combine (FILTER, SIFT, PolyPhen).

Creates:
  - all_missense_variants.parquet
  - verified_missense_variants.parquet
  - mane_select_missense_variants.parquet
  - gold_standard_missense_variants.parquet
and writes statistics.json

Usage:
  python -u combine_missense_outputs.py \
      --output_dir gnomad_missense_validated \
      --batch_size 500000
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.compute as pc
import pyarrow.parquet as pq


DROP_COLS = ("FILTER", "SIFT", "PolyPhen")


def stream_write(files: list[Path], out_path: Path, kind: str, batch_size: int) -> int:
    """
    Stream through all input files and write a single combined parquet.

    kind in {"all","verified","mane","gold"}:
      - all: no filter
      - verified: seq_verified_i8 == 1
      - mane: is_mane_select == True
      - gold: both
    """
    writer = None
    n_written = 0
    t0 = time.time()

    for fp in files:
        # scan file independently (avoids schema unification casting issues)
        dataset = ds.dataset(str(fp), format="parquet")
        scanner = dataset.scanner(batch_size=batch_size)

        for rb in scanner.to_batches():
            table = pa.Table.from_batches([rb])

            # Drop problematic columns to ensure identical schema across batches/files
            drop = [c for c in DROP_COLS if c in table.schema.names]
            if drop:
                table = table.drop_columns(drop)

            # Ensure filter columns exist with stable types
            # seq_verified_i8: int8, is_mane_select: bool
            if "seq_verified_i8" not in table.schema.names:
                table = table.append_column("seq_verified_i8", pa.array([-1] * table.num_rows, type=pa.int8()))
            else:
                col = table["seq_verified_i8"]
                if pa.types.is_null(col.type):
                    table = table.set_column(
                        table.schema.get_field_index("seq_verified_i8"),
                        "seq_verified_i8",
                        pa.array([-1] * table.num_rows, type=pa.int8()),
                    )
                elif col.type != pa.int8():
                    table = table.set_column(
                        table.schema.get_field_index("seq_verified_i8"),
                        "seq_verified_i8",
                        pc.cast(col, pa.int8(), safe=False),
                    )

            if "is_mane_select" not in table.schema.names:
                table = table.append_column("is_mane_select", pa.array([False] * table.num_rows, type=pa.bool_()))
            else:
                col = table["is_mane_select"]
                if pa.types.is_null(col.type):
                    table = table.set_column(
                        table.schema.get_field_index("is_mane_select"),
                        "is_mane_select",
                        pa.array([False] * table.num_rows, type=pa.bool_()),
                    )
                elif col.type != pa.bool_():
                    table = table.set_column(
                        table.schema.get_field_index("is_mane_select"),
                        "is_mane_select",
                        pc.cast(col, pa.bool_(), safe=False),
                    )

            # Apply selection filter
            if kind == "verified":
                mask = pc.equal(table["seq_verified_i8"], pa.scalar(1, pa.int8()))
                table = table.filter(mask)
            elif kind == "mane":
                mask = pc.equal(table["is_mane_select"], pa.scalar(True, pa.bool_()))
                table = table.filter(mask)
            elif kind == "gold":
                m1 = pc.equal(table["seq_verified_i8"], pa.scalar(1, pa.int8()))
                m2 = pc.equal(table["is_mane_select"], pa.scalar(True, pa.bool_()))
                table = table.filter(pc.and_(m1, m2))
            # kind == "all": no filter

            if table.num_rows == 0:
                continue

            # Initialize writer with the first non-empty batch schema
            if writer is None:
                out_path.parent.mkdir(parents=True, exist_ok=True)
                writer = pq.ParquetWriter(str(out_path), table.schema, compression="snappy")

            # IMPORTANT: enforce same column order as writer schema
            # (parquet writer is strict)
            if table.schema.names != writer.schema.names:
                # re-order columns to match writer
                table = table.select(writer.schema.names)

            writer.write_table(table)
            n_written += table.num_rows

    if writer is not None:
        writer.close()

    minutes = (time.time() - t0) / 60
    print(f"[{kind}] wrote {n_written:,} rows -> {out_path} in {minutes:.2f} min")
    return n_written


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output_dir", type=str, required=True, help="Directory containing chr*_missense_enriched.parquet")
    ap.add_argument("--batch_size", type=int, default=500_000)
    ap.add_argument("--keep_existing", action="store_true", help="Do not delete existing combined outputs")
    args = ap.parse_args()

    outdir = Path(args.output_dir)
    files = sorted(outdir.glob("chr*_missense_enriched.parquet"))
    if not files:
        raise SystemExit(f"No chr*_missense_enriched.parquet found in {outdir.resolve()}")

    # Output paths
    out_all = outdir / "all_missense_variants.parquet"
    out_ver = outdir / "verified_missense_variants.parquet"
    out_mane = outdir / "mane_select_missense_variants.parquet"
    out_gold = outdir / "gold_standard_missense_variants.parquet"

    if not args.keep_existing:
        for p in (out_all, out_ver, out_mane, out_gold):
            if p.exists():
                p.unlink()

    t0 = time.time()
    n_all = stream_write(files, out_all, "all", args.batch_size)
    n_ver = stream_write(files, out_ver, "verified", args.batch_size)
    n_man = stream_write(files, out_mane, "mane", args.batch_size)
    n_gol = stream_write(files, out_gold, "gold", args.batch_size)

    stats = {
        "per_chrom_files": len(files),
        "combined_rows": {"all": n_all, "verified": n_ver, "mane": n_man, "gold": n_gol},
        "batch_size": args.batch_size,
        "dropped_columns": list(DROP_COLS),
        "minutes": round((time.time() - t0) / 60, 2),
    }
    with open(outdir / "statistics.json", "w") as f:
        json.dump(stats, f, indent=2)

    print("\n✅ combine complete")
    print(json.dumps(stats, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
