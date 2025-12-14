#!/usr/bin/env python3
"""
Validate and Enrich All gnomAD Missense Variants

This script takes the raw gnomAD results and:
1. Parses protein-level information (HGVSp)
2. Attaches protein sequences from local Ensembl
3. Validates sequences
4. Adds UniProt and AlphaFold IDs
5. Creates analysis-ready dataset

Usage:
    python validate_all_missense.py
"""

import pandas as pd
import re
import gzip
from pathlib import Path
from typing import Dict, Optional
import logging
from collections import defaultdict

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# ============================================================================
# CONFIGURATION
# ============================================================================

# Paths
REFERENCE_DIR = Path("/jet/home/barazand/NEWOCEAN/ref_data")
ENSEMBL_GTF = REFERENCE_DIR / "ensembl" / "Homo_sapiens.GRCh38.112.gtf.gz"
ENSEMBL_PEP = REFERENCE_DIR / "ensembl" / "Homo_sapiens.GRCh38.pep.all.fa.gz"
MANE_SUMMARY = REFERENCE_DIR / "mane" / "MANE.GRCh38.v1.3.summary.txt.gz"
UNIPROT_FASTA = REFERENCE_DIR / "uniprot" / "uniprot_sprot.fasta.gz"
UNIPROT_IDMAPPING = REFERENCE_DIR / "uniprot" / "HUMAN_9606_idmapping.dat.gz"

# Input/Output
INPUT_DIR = Path("gnomad_all_genes")
INPUT_FILE = INPUT_DIR / "all_gnomad_missense_variants.parquet"
OUTPUT_DIR = Path("gnomad_all_genes_validated")
OUTPUT_DIR.mkdir(exist_ok=True)

# Global caches
ENST_TO_SEQ = None
GENE_TO_UNIPROT = None
UNIPROT_TO_SEQ = None


# ============================================================================
# REFERENCE DATA LOADERS (same as before)
# ============================================================================

def load_enst_to_ensp(gtf_path: Path) -> Dict[str, str]:
    """Parse Ensembl GTF to build ENST→ENSP mapping."""
    logger.info(f"   Parsing GTF: {gtf_path.name}")
    enst_to_ensp = {}
    
    with gzip.open(gtf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            
            attr_dict = {}
            for field in parts[8].split(";"):
                field = field.strip()
                if not field or " " not in field:
                    continue
                key, value = field.split(" ", 1)
                attr_dict[key] = value.strip('"')
            
            tid = attr_dict.get("transcript_id")
            pid = attr_dict.get("protein_id")
            if tid and pid:
                enst_to_ensp.setdefault(tid, pid)
                enst_to_ensp.setdefault(tid.split(".")[0], pid)
    
    logger.info(f"      Found {len(enst_to_ensp)} ENST→ENSP mappings")
    return enst_to_ensp


def load_ensp_to_seq(pep_path: Path) -> Dict[str, str]:
    """Parse Ensembl peptide FASTA to build ENSP→sequence mapping."""
    logger.info(f"   Parsing peptide FASTA: {pep_path.name}")
    ensp_to_seq = {}
    current_id = None
    seq_chunks = []
    
    with gzip.open(pep_path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if current_id:
                    seq = "".join(seq_chunks)
                    ensp_to_seq.setdefault(current_id, seq)
                    ensp_to_seq.setdefault(current_id.split(".")[0], seq)
                
                current_id = line[1:].strip().split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        
        if current_id:
            seq = "".join(seq_chunks)
            ensp_to_seq.setdefault(current_id, seq)
            ensp_to_seq.setdefault(current_id.split(".")[0], seq)
    
    logger.info(f"      Found {len(ensp_to_seq)} ENSP→sequence mappings")
    return ensp_to_seq


def build_enst_to_seq(gtf_path: Path, pep_path: Path) -> Dict[str, str]:
    """Build ENST→sequence mapping."""
    enst_to_ensp = load_enst_to_ensp(gtf_path)
    ensp_to_seq = load_ensp_to_seq(pep_path)
    
    enst_to_seq = {}
    for enst, ensp in enst_to_ensp.items():
        seq = ensp_to_seq.get(ensp) or ensp_to_seq.get(ensp.split(".")[0])
        if seq:
            enst_to_seq[enst] = seq
            enst_to_seq.setdefault(enst.split(".")[0], seq)
    
    logger.info(f"   ✅ Built {len(enst_to_seq)} ENST→sequence mappings")
    return enst_to_seq


def load_uniprot_id_mapping(idmapping_path: Path) -> Dict[str, str]:
    """Load gene→UniProt mapping."""
    logger.info(f"   Parsing UniProt ID mapping: {idmapping_path.name}")
    gene_to_uniprot = {}
    
    with gzip.open(idmapping_path, "rt") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) == 3 and parts[1] == "Gene_Name":
                if parts[2] not in gene_to_uniprot:
                    gene_to_uniprot[parts[2]] = parts[0]
    
    logger.info(f"      Found {len(gene_to_uniprot)} gene→UniProt mappings")
    return gene_to_uniprot


def load_uniprot_sequences(fasta_path: Path) -> Dict[str, str]:
    """Load UniProt→sequence mapping."""
    logger.info(f"   Parsing UniProt FASTA: {fasta_path.name}")
    uniprot_to_seq = {}
    current_id = None
    seq_chunks = []
    
    with gzip.open(fasta_path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if current_id:
                    uniprot_to_seq[current_id] = "".join(seq_chunks)
                
                parts = line[1:].strip().split("|")
                current_id = parts[1] if len(parts) >= 2 else line[1:].strip().split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        
        if current_id:
            uniprot_to_seq[current_id] = "".join(seq_chunks)
    
    logger.info(f"      Found {len(uniprot_to_seq)} UniProt→sequence mappings")
    return uniprot_to_seq


def initialize_reference_data():
    """Load all reference data."""
    global ENST_TO_SEQ, GENE_TO_UNIPROT, UNIPROT_TO_SEQ
    
    logger.info("\n" + "="*80)
    logger.info("LOADING REFERENCE DATA")
    logger.info("="*80)
    
    logger.info("\n📚 Loading Ensembl data...")
    ENST_TO_SEQ = build_enst_to_seq(ENSEMBL_GTF, ENSEMBL_PEP)
    
    logger.info("\n📚 Loading UniProt data...")
    GENE_TO_UNIPROT = load_uniprot_id_mapping(UNIPROT_IDMAPPING)
    UNIPROT_TO_SEQ = load_uniprot_sequences(UNIPROT_FASTA)
    
    logger.info("\n✅ All reference data loaded")


# ============================================================================
# PROCESSING FUNCTIONS
# ============================================================================

def parse_hgvsp_notation(hgvsp: str) -> Optional[Dict]:
    """Parse HGVSp notation."""
    if pd.isna(hgvsp) or hgvsp == '':
        return None
    
    aa_3to1 = {
        'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
        'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
        'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
        'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V',
        'Ter': '*', 'Stop': '*', 'Sec': 'U', 'Pyl': 'O'
    }
    
    hgvsp_str = str(hgvsp)
    if ':p.' in hgvsp_str:
        hgvsp_str = 'p.' + hgvsp_str.split(':p.')[1]
    
    pattern = r'p\.([A-Z][a-z]{2}|[A-Z\*])(\d+)([A-Z][a-z]{2}|[A-Z\*\=])'
    match = re.match(pattern, hgvsp_str)
    
    if match:
        ref_aa = aa_3to1.get(match.group(1), match.group(1))
        pos = int(match.group(2))
        alt_aa = match.group(3)
        alt_aa = aa_3to1.get(alt_aa, alt_aa) if alt_aa != '=' else ref_aa
        
        return {'ref_aa': ref_aa, 'pos': pos, 'alt_aa': alt_aa}
    
    return None


def get_protein_sequence(transcript_id: str) -> Optional[str]:
    """Get protein sequence from cache."""
    if not transcript_id or ENST_TO_SEQ is None:
        return None
    return ENST_TO_SEQ.get(transcript_id) or ENST_TO_SEQ.get(transcript_id.split(".")[0])


def get_uniprot_id(gene_symbol: str) -> Optional[str]:
    """Get UniProt ID for gene."""
    if GENE_TO_UNIPROT is None:
        return None
    return GENE_TO_UNIPROT.get(gene_symbol)


def process_batch(df: pd.DataFrame, batch_num: int, total_batches: int) -> pd.DataFrame:
    """Process a batch of variants."""
    logger.info(f"\nProcessing batch {batch_num}/{total_batches} ({len(df):,} variants)")
    
    # Parse HGVSp
    df['protein_change'] = df['HGVSp'].apply(parse_hgvsp_notation)
    df['ref_aa'] = df['protein_change'].apply(lambda x: x['ref_aa'] if x else None)
    df['protein_pos'] = df['protein_change'].apply(lambda x: x['pos'] if x else None)
    df['alt_aa'] = df['protein_change'].apply(lambda x: x['alt_aa'] if x else None)
    df['mutation'] = df.apply(
        lambda row: f"{row['ref_aa']}{row['protein_pos']}{row['alt_aa']}" if row['ref_aa'] else None,
        axis=1
    )
    
    # Attach protein sequences
    df['protein_seq'] = df['Feature'].apply(get_protein_sequence)
    
    # Verify sequences
    def verify_match(row):
        if pd.isna(row['protein_seq']) or pd.isna(row['protein_pos']):
            return None
        try:
            return row['protein_seq'][row['protein_pos'] - 1] == row['ref_aa']
        except:
            return None
    
    df['seq_verified'] = df.apply(verify_match, axis=1)
    
    # Add UniProt IDs
    df['uniprot_id'] = df['gene_symbol'].apply(get_uniprot_id)
    df['alphafold_id'] = df['uniprot_id'].apply(
        lambda x: f"AF-{x}-F1" if pd.notna(x) else None
    )
    
    verified = df['seq_verified'].sum()
    total_with_seq = df['protein_seq'].notna().sum()
    logger.info(f"  Sequences: {total_with_seq:,}/{len(df):,}")
    logger.info(f"  Verified: {verified:,}/{total_with_seq:,}")
    
    return df


def process_all_variants():
    """Process all variants with batching for memory efficiency."""
    logger.info("\n" + "="*80)
    logger.info("PROCESSING ALL VARIANTS")
    logger.info("="*80)
    
    if not INPUT_FILE.exists():
        raise FileNotFoundError(f"Input file not found: {INPUT_FILE}")
    
    logger.info(f"\nReading: {INPUT_FILE}")
    df = pd.read_parquet(INPUT_FILE)
    logger.info(f"Total variants: {len(df):,}")
    
    # Process in batches to manage memory
    batch_size = 100000
    total_batches = (len(df) + batch_size - 1) // batch_size
    
    processed_batches = []
    for i in range(0, len(df), batch_size):
        batch = df.iloc[i:i+batch_size].copy()
        batch_num = (i // batch_size) + 1
        
        processed = process_batch(batch, batch_num, total_batches)
        processed_batches.append(processed)
    
    # Combine all batches
    logger.info("\nCombining batches...")
    final_df = pd.concat(processed_batches, ignore_index=True)
    
    return final_df


def generate_statistics(df: pd.DataFrame):
    """Generate comprehensive statistics."""
    logger.info("\n" + "="*80)
    logger.info("STATISTICS")
    logger.info("="*80)
    
    stats = {
        'total_variants': len(df),
        'unique_genes': df['gene_symbol'].nunique(),
        'unique_mutations': df['mutation'].nunique(),
        'with_sequence': df['protein_seq'].notna().sum(),
        'sequence_verified': df['seq_verified'].sum(),
        'with_uniprot': df['uniprot_id'].notna().sum(),
        'canonical_only': (df['CANONICAL'] == 'YES').sum(),
    }
    
    logger.info(f"\nTotal variants: {stats['total_variants']:,}")
    logger.info(f"Unique genes: {stats['unique_genes']:,}")
    logger.info(f"Unique mutations: {stats['unique_mutations']:,}")
    logger.info(f"With sequence: {stats['with_sequence']:,} ({100*stats['with_sequence']/stats['total_variants']:.1f}%)")
    logger.info(f"Sequence verified: {stats['sequence_verified']:,} ({100*stats['sequence_verified']/stats['with_sequence']:.1f}%)")
    logger.info(f"With UniProt ID: {stats['with_uniprot']:,} ({100*stats['with_uniprot']/stats['total_variants']:.1f}%)")
    logger.info(f"Canonical transcripts: {stats['canonical_only']:,}")
    
    # Save statistics
    stats_file = OUTPUT_DIR / "statistics.json"
    import json
    with open(stats_file, 'w') as f:
        json.dump(stats, f, indent=2)
    
    logger.info(f"\n✅ Statistics saved to {stats_file}")
    
    return stats


def save_results(df: pd.DataFrame):
    """Save processed results."""
    logger.info("\n" + "="*80)
    logger.info("SAVING RESULTS")
    logger.info("="*80)
    
    # Full dataset
    full_parquet = OUTPUT_DIR / "all_missense_validated.parquet"
    full_csv = OUTPUT_DIR / "all_missense_validated.csv.gz"
    
    df.to_parquet(full_parquet, index=False)
    df.to_csv(full_csv, index=False, compression='gzip')
    
    logger.info(f"\n✅ Full dataset:")
    logger.info(f"  {full_parquet} ({full_parquet.stat().st_size / 1e6:.1f} MB)")
    logger.info(f"  {full_csv} ({full_csv.stat().st_size / 1e6:.1f} MB)")
    
    # Canonical only (high-confidence)
    canonical = df[df['CANONICAL'] == 'YES'].copy()
    canonical_parquet = OUTPUT_DIR / "canonical_missense_validated.parquet"
    canonical.to_parquet(canonical_parquet, index=False)
    
    logger.info(f"\n✅ Canonical only:")
    logger.info(f"  {canonical_parquet} ({canonical_parquet.stat().st_size / 1e6:.1f} MB)")
    logger.info(f"  {len(canonical):,} variants")
    
    # Verified sequences only (highest confidence)
    verified = df[df['seq_verified'] == True].copy()
    verified_parquet = OUTPUT_DIR / "verified_missense.parquet"
    verified.to_parquet(verified_parquet, index=False)
    
    logger.info(f"\n✅ Verified sequences only:")
    logger.info(f"  {verified_parquet} ({verified_parquet.stat().st_size / 1e6:.1f} MB)")
    logger.info(f"  {len(verified):,} variants")


def main():
    """Main execution."""
    logger.info("\n" + "="*80)
    logger.info("VALIDATE ALL GNOMAD MISSENSE VARIANTS")
    logger.info("="*80)
    
    try:
        # Load reference data
        initialize_reference_data()
        
        # Process variants
        df = process_all_variants()
        
        # Generate statistics
        generate_statistics(df)
        
        # Save results
        save_results(df)
        
        logger.info("\n" + "="*80)
        logger.info("✅ PIPELINE COMPLETE")
        logger.info("="*80)
        
        return True
        
    except Exception as e:
        logger.error(f"\n❌ PIPELINE FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == "__main__":
    success = main()
    exit(0 if success else 1)