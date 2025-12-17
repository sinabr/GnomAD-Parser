#!/usr/bin/env python3
"""
Process Chromosome Parquet Files → Validated Missense Variants

Takes the chromosome-level Parquet files created by convert_vcf_to_dataframe_fast.py
and produces a validated missense variant dataset with protein sequences.

Usage:
    python process_missense_from_parquet.py
"""

import pandas as pd
import re
import gzip
from pathlib import Path
from typing import Dict, Optional, List
import logging
from collections import defaultdict
import time

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# ============================================================================
# CONFIGURATION
# ============================================================================

# Paths - UPDATE THESE TO MATCH YOUR SYSTEM
REFERENCE_DIR = Path("/projects/lugoteam/protein_graphs/GnomAD-Parser/ref_data")
ENSEMBL_GTF = REFERENCE_DIR / "ensembl" / "Homo_sapiens.GRCh38.112.gtf.gz"
ENSEMBL_PEP = REFERENCE_DIR / "ensembl" / "Homo_sapiens.GRCh38.pep.all.fa.gz"
MANE_SUMMARY = REFERENCE_DIR / "mane" / "MANE.GRCh38.v1.3.summary.txt.gz"
UNIPROT_FASTA = REFERENCE_DIR / "uniprot" / "uniprot_sprot.fasta.gz"
UNIPROT_IDMAPPING = REFERENCE_DIR / "uniprot" / "HUMAN_9606_idmapping.dat.gz"

# Input/Output
INPUT_DIR = Path("gnomad_all_genes/chromosome_dataframes")
OUTPUT_DIR = Path("gnomad_missense_validated")
OUTPUT_DIR.mkdir(exist_ok=True)

# Global caches
ENST_TO_SEQ = None
ENSP_TO_SEQ = None
GENE_TO_UNIPROT = None
UNIPROT_TO_SEQ = None
MANE_TRANSCRIPTS = None

# VEP fields we care about for missense variant identification
REQUIRED_VEP_FIELDS = [
    'Consequence',      # Must contain 'missense_variant'
    'SYMBOL',           # Gene symbol
    'Gene',             # Ensembl gene ID
    'Feature',          # Transcript ID (ENST)
    'HGVSc',            # DNA-level change
    'HGVSp',            # Protein-level change
    'CANONICAL',        # Is this the canonical transcript?
    'MANE_SELECT',      # MANE Select transcript
    'MANE_PLUS_CLINICAL', # MANE Plus Clinical
    'BIOTYPE',          # protein_coding, etc.
    'IMPACT',           # HIGH, MODERATE, LOW, MODIFIER
    'PolyPhen',         # PolyPhen prediction
    'SIFT',             # SIFT prediction
]


# ============================================================================
# REFERENCE DATA LOADERS
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
                # Store both versioned and unversioned
                enst_to_ensp[tid] = pid
                enst_to_ensp[tid.split(".")[0]] = pid
    
    logger.info(f"      Found {len(enst_to_ensp):,} ENST→ENSP mappings")
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
                    # Store both versioned and unversioned
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
    """Build ENST→sequence mapping."""
    enst_to_ensp = load_enst_to_ensp(gtf_path)
    ensp_to_seq = load_ensp_to_seq(pep_path)
    
    enst_to_seq = {}
    for enst, ensp in enst_to_ensp.items():
        # Try versioned first, then unversioned
        seq = ensp_to_seq.get(ensp) or ensp_to_seq.get(ensp.split(".")[0])
        if seq:
            enst_to_seq[enst] = seq
            # Also store unversioned
            enst_to_seq[enst.split(".")[0]] = seq
    
    logger.info(f"   ✅ Built {len(enst_to_seq):,} ENST→sequence mappings")
    return enst_to_seq


def load_mane_transcripts(mane_path: Path) -> Dict[str, Dict]:
    """Load MANE Select transcripts for gene→transcript mapping."""
    logger.info(f"   Parsing MANE summary: {mane_path.name}")
    mane_data = {}
    
    with gzip.open(mane_path, "rt") as f:
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
            enst = parts[col_idx.get("Ensembl_nuc", col_idx.get("Ensembl_transcript", -1))]
            nm = parts[col_idx.get("RefSeq_nuc", col_idx.get("RefSeq_transcript", -1))]
            
            if symbol and enst:
                mane_data[symbol] = {
                    'enst': enst.split('.')[0],
                    'refseq': nm if nm else ''
                }
    
    logger.info(f"      Found {len(mane_data):,} MANE Select transcripts")
    return mane_data


def load_uniprot_id_mapping(idmapping_path: Path) -> Dict[str, str]:
    """Load gene→UniProt mapping."""
    logger.info(f"   Parsing UniProt ID mapping: {idmapping_path.name}")
    gene_to_uniprot = {}
    
    with gzip.open(idmapping_path, "rt") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) == 3 and parts[1] == "Gene_Name":
                # Keep first mapping (canonical)
                if parts[2] not in gene_to_uniprot:
                    gene_to_uniprot[parts[2]] = parts[0]
    
    logger.info(f"      Found {len(gene_to_uniprot):,} gene→UniProt mappings")
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
    
    logger.info(f"      Found {len(uniprot_to_seq):,} UniProt→sequence mappings")
    return uniprot_to_seq


def initialize_reference_data():
    """Load all reference data into global caches."""
    global ENST_TO_SEQ, ENSP_TO_SEQ, GENE_TO_UNIPROT, UNIPROT_TO_SEQ, MANE_TRANSCRIPTS
    
    logger.info("\n" + "="*80)
    logger.info("LOADING REFERENCE DATA")
    logger.info("="*80)
    
    logger.info("\n📚 Loading Ensembl data...")
    ENST_TO_SEQ = build_enst_to_seq(ENSEMBL_GTF, ENSEMBL_PEP)
    
    logger.info("\n📚 Loading MANE data...")
    MANE_TRANSCRIPTS = load_mane_transcripts(MANE_SUMMARY)
    
    logger.info("\n📚 Loading UniProt data...")
    GENE_TO_UNIPROT = load_uniprot_id_mapping(UNIPROT_IDMAPPING)
    UNIPROT_TO_SEQ = load_uniprot_sequences(UNIPROT_FASTA)
    
    logger.info("\n✅ All reference data loaded")


# ============================================================================
# VEP FIELD PARSING
# ============================================================================

def parse_hgvsp_notation(hgvsp: str) -> Optional[Dict]:
    """
    Parse HGVSp notation to extract ref_aa, position, alt_aa.
    
    Examples:
        p.Arg123Cys → {'ref_aa': 'R', 'pos': 123, 'alt_aa': 'C'}
        ENSP00000123:p.R123C → {'ref_aa': 'R', 'pos': 123, 'alt_aa': 'C'}
    """
    if pd.isna(hgvsp) or hgvsp == '':
        return None
    
    # 3-letter to 1-letter AA code mapping
    aa_3to1 = {
        'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
        'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
        'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
        'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V',
        'Ter': '*', 'Stop': '*', 'Sec': 'U', 'Pyl': 'O'
    }
    
    hgvsp_str = str(hgvsp)
    
    # Strip protein ID prefix if present (e.g., "ENSP00000123:p.Arg123Cys" → "p.Arg123Cys")
    if ':p.' in hgvsp_str:
        hgvsp_str = 'p.' + hgvsp_str.split(':p.')[1]
    
    # Match 3-letter AA code pattern: p.Arg123Cys
    pattern_3letter = r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|\=)'
    match = re.match(pattern_3letter, hgvsp_str)
    
    if match:
        ref_aa_3 = match.group(1)
        pos = int(match.group(2))
        alt_aa_3 = match.group(3)
        
        ref_aa = aa_3to1.get(ref_aa_3, ref_aa_3)
        alt_aa = aa_3to1.get(alt_aa_3, alt_aa_3) if alt_aa_3 != '=' else ref_aa
        
        return {'ref_aa': ref_aa, 'pos': pos, 'alt_aa': alt_aa}
    
    # Match 1-letter AA code pattern: p.R123C
    pattern_1letter = r'p\.([A-Z\*])(\d+)([A-Z\*\=])'
    match = re.match(pattern_1letter, hgvsp_str)
    
    if match:
        ref_aa = match.group(1)
        pos = int(match.group(2))
        alt_aa = match.group(3) if match.group(3) != '=' else ref_aa
        
        return {'ref_aa': ref_aa, 'pos': pos, 'alt_aa': alt_aa}
    
    return None


# ============================================================================
# CHROMOSOME PARQUET PROCESSING
# ============================================================================

def load_chromosome_parquet(chrom: str) -> pd.DataFrame:
    """Load a single chromosome Parquet file."""
    parquet_file = INPUT_DIR / f"chr{chrom}_variants.parquet"
    
    if not parquet_file.exists():
        logger.warning(f"  ⚠️  chr{chrom}: file not found")
        return pd.DataFrame()
    
    logger.info(f"  Loading chr{chrom}...")
    df = pd.read_parquet(parquet_file)
    
    logger.info(f"    Total variants: {len(df):,}")
    return df


def filter_missense_variants(df: pd.DataFrame) -> pd.DataFrame:
    """Filter to missense variants only."""
    if df.empty:
        return df
    
    # Check if Consequence field exists
    if 'Consequence' not in df.columns:
        logger.error("    ❌ No 'Consequence' field found in VEP annotations!")
        return pd.DataFrame()
    
    # Filter for missense variants
    # Consequence can be a comma-separated list, so check if 'missense_variant' is in it
    missense_mask = df['Consequence'].str.contains('missense_variant', case=False, na=False)
    missense_df = df[missense_mask].copy()
    
    logger.info(f"    Missense variants: {len(missense_df):,}")
    
    return missense_df


def enrich_with_protein_data(df: pd.DataFrame, chrom: str) -> pd.DataFrame:
    """Enrich variants with protein sequences and validation."""
    if df.empty:
        return df
    
    logger.info(f"  Enriching chr{chrom} with protein data...")
    
    # Parse HGVSp notation
    df['protein_change'] = df['HGVSp'].apply(parse_hgvsp_notation)
    df['ref_aa'] = df['protein_change'].apply(lambda x: x['ref_aa'] if x else None)
    df['protein_pos'] = df['protein_change'].apply(lambda x: x['pos'] if x else None)
    df['alt_aa'] = df['protein_change'].apply(lambda x: x['alt_aa'] if x else None)
    
    # Create mutation string (e.g., "R123C")
    df['mutation'] = df.apply(
        lambda row: f"{row['ref_aa']}{row['protein_pos']}{row['alt_aa']}" 
        if row['ref_aa'] else None,
        axis=1
    )
    
    # Get protein sequence from Ensembl via transcript ID
    def get_ensembl_seq(feature):
        if pd.isna(feature) or not feature:
            return None
        # Try both versioned and unversioned
        return ENST_TO_SEQ.get(feature) or ENST_TO_SEQ.get(feature.split(".")[0])
    
    df['protein_seq_ensembl'] = df['Feature'].apply(get_ensembl_seq)
    
    # Verify sequence match
    def verify_sequence(row):
        if pd.isna(row['protein_seq_ensembl']) or pd.isna(row['protein_pos']) or pd.isna(row['ref_aa']):
            return None
        try:
            pos_idx = int(row['protein_pos']) - 1  # 0-indexed
            if pos_idx < 0 or pos_idx >= len(row['protein_seq_ensembl']):
                return False
            return row['protein_seq_ensembl'][pos_idx] == row['ref_aa']
        except:
            return None
    
    df['seq_verified'] = df.apply(verify_sequence, axis=1)
    
    # Add UniProt IDs based on gene symbol
    def get_uniprot_id(symbol):
        if pd.isna(symbol):
            return None
        return GENE_TO_UNIPROT.get(symbol)
    
    df['uniprot_id'] = df['SYMBOL'].apply(get_uniprot_id)
    
    # Add AlphaFold IDs
    df['alphafold_id'] = df['uniprot_id'].apply(
        lambda x: f"AF-{x}-F1" if pd.notna(x) else None
    )
    
    # Get UniProt sequence for comparison
    df['protein_seq_uniprot'] = df['uniprot_id'].apply(
        lambda x: UNIPROT_TO_SEQ.get(x) if pd.notna(x) else None
    )
    
    # Flag MANE Select transcripts
    df['is_mane_select'] = df.apply(
        lambda row: (
            pd.notna(row['SYMBOL']) and 
            row['SYMBOL'] in MANE_TRANSCRIPTS and
            row['Feature'].split('.')[0] == MANE_TRANSCRIPTS[row['SYMBOL']]['enst']
        ) if pd.notna(row['Feature']) else False,
        axis=1
    )
    
    # Stats
    with_seq = df['protein_seq_ensembl'].notna().sum()
    verified = df['seq_verified'].sum()
    with_uniprot = df['uniprot_id'].notna().sum()
    mane_select = df['is_mane_select'].sum()
    
    logger.info(f"    With Ensembl sequence: {with_seq:,}/{len(df):,} ({100*with_seq/len(df):.1f}%)")
    if with_seq > 0:
        logger.info(f"    Sequence verified: {verified:,}/{with_seq:,} ({100*verified/with_seq:.1f}%)")
    logger.info(f"    With UniProt ID: {with_uniprot:,}/{len(df):,} ({100*with_uniprot/len(df):.1f}%)")
    logger.info(f"    MANE Select: {mane_select:,}/{len(df):,} ({100*mane_select/len(df):.1f}%)")
    
    return df


def process_all_chromosomes() -> pd.DataFrame:
    """Process all chromosome Parquet files."""
    logger.info("\n" + "="*80)
    logger.info("PROCESSING CHROMOSOME PARQUET FILES")
    logger.info("="*80)
    
    chromosomes = [str(i) for i in range(1, 23)] + ['X', 'Y']
    
    all_missense = []
    
    for chrom in chromosomes:
        logger.info(f"\n📊 Processing chromosome {chrom}")
        
        # Load chromosome data
        df = load_chromosome_parquet(chrom)
        if df.empty:
            continue
        
        # Filter to missense
        missense_df = filter_missense_variants(df)
        if missense_df.empty:
            continue
        
        # Enrich with protein data
        enriched_df = enrich_with_protein_data(missense_df, chrom)
        
        all_missense.append(enriched_df)
    
    if not all_missense:
        logger.error("\n❌ No missense variants found!")
        return pd.DataFrame()
    
    logger.info("\n" + "="*80)
    logger.info("COMBINING ALL CHROMOSOMES")
    logger.info("="*80)
    
    combined = pd.concat(all_missense, ignore_index=True)
    logger.info(f"Total missense variants: {len(combined):,}")
    
    return combined


# ============================================================================
# RESULTS AND STATISTICS
# ============================================================================

def generate_statistics(df: pd.DataFrame) -> Dict:
    """Generate comprehensive statistics."""
    logger.info("\n" + "="*80)
    logger.info("STATISTICS")
    logger.info("="*80)
    
    stats = {
        'total_variants': len(df),
        'unique_genes': df['SYMBOL'].nunique() if 'SYMBOL' in df.columns else 0,
        'unique_mutations': df['mutation'].nunique(),
        'by_source': {
            'exomes': (df['source'] == 'exomes').sum(),
            'genomes': (df['source'] == 'genomes').sum(),
        },
        'with_ensembl_seq': df['protein_seq_ensembl'].notna().sum(),
        'sequence_verified': df['seq_verified'].sum(),
        'with_uniprot': df['uniprot_id'].notna().sum(),
        'canonical': (df['CANONICAL'] == 'YES').sum() if 'CANONICAL' in df.columns else 0,
        'mane_select': df['is_mane_select'].sum(),
        'by_impact': {},
    }
    
    # Impact distribution
    if 'IMPACT' in df.columns:
        for impact in df['IMPACT'].unique():
            if pd.notna(impact):
                stats['by_impact'][impact] = (df['IMPACT'] == impact).sum()
    
    # Log statistics
    logger.info(f"\nTotal missense variants: {stats['total_variants']:,}")
    logger.info(f"Unique genes: {stats['unique_genes']:,}")
    logger.info(f"Unique mutations: {stats['unique_mutations']:,}")
    
    logger.info(f"\nBy source:")
    logger.info(f"  Exomes: {stats['by_source']['exomes']:,}")
    logger.info(f"  Genomes: {stats['by_source']['genomes']:,}")
    
    logger.info(f"\nProtein data:")
    logger.info(f"  With Ensembl sequence: {stats['with_ensembl_seq']:,} ({100*stats['with_ensembl_seq']/stats['total_variants']:.1f}%)")
    logger.info(f"  Sequence verified: {stats['sequence_verified']:,} ({100*stats['sequence_verified']/stats['with_ensembl_seq']:.1f}% of those with seq)")
    logger.info(f"  With UniProt ID: {stats['with_uniprot']:,} ({100*stats['with_uniprot']/stats['total_variants']:.1f}%)")
    
    logger.info(f"\nTranscript quality:")
    logger.info(f"  Canonical: {stats['canonical']:,} ({100*stats['canonical']/stats['total_variants']:.1f}%)")
    logger.info(f"  MANE Select: {stats['mane_select']:,} ({100*stats['mane_select']/stats['total_variants']:.1f}%)")
    
    if stats['by_impact']:
        logger.info(f"\nBy VEP impact:")
        for impact, count in sorted(stats['by_impact'].items(), key=lambda x: x[1], reverse=True):
            logger.info(f"  {impact}: {count:,}")
    
    return stats


def save_results(df: pd.DataFrame, stats: Dict):
    """Save processed results in multiple formats."""
    logger.info("\n" + "="*80)
    logger.info("SAVING RESULTS")
    logger.info("="*80)
    
    # Save statistics
    import json
    stats_file = OUTPUT_DIR / "statistics.json"
    with open(stats_file, 'w') as f:
        json.dump(stats, f, indent=2)
    logger.info(f"\n✅ Statistics: {stats_file}")
    
    # 1. Full dataset (all missense variants)
    full_parquet = OUTPUT_DIR / "all_missense_variants.parquet"
    full_csv = OUTPUT_DIR / "all_missense_variants.csv.gz"
    
    df.to_parquet(full_parquet, index=False)
    df.to_csv(full_csv, index=False, compression='gzip')
    
    logger.info(f"\n✅ Full dataset:")
    logger.info(f"  {full_parquet} ({full_parquet.stat().st_size / 1e6:.1f} MB)")
    logger.info(f"  {full_csv} ({full_csv.stat().st_size / 1e6:.1f} MB)")
    logger.info(f"  Variants: {len(df):,}")
    
    # 2. Sequence-verified only (highest confidence)
    verified = df[df['seq_verified'] == True].copy()
    if len(verified) > 0:
        verified_parquet = OUTPUT_DIR / "verified_missense_variants.parquet"
        verified_csv = OUTPUT_DIR / "verified_missense_variants.csv.gz"
        
        verified.to_parquet(verified_parquet, index=False)
        verified.to_csv(verified_csv, index=False, compression='gzip')
        
        logger.info(f"\n✅ Sequence-verified only (highest confidence):")
        logger.info(f"  {verified_parquet} ({verified_parquet.stat().st_size / 1e6:.1f} MB)")
        logger.info(f"  {verified_csv} ({verified_csv.stat().st_size / 1e6:.1f} MB)")
        logger.info(f"  Variants: {len(verified):,}")
    
    # 3. MANE Select only (standard clinical transcripts)
    mane = df[df['is_mane_select'] == True].copy()
    if len(mane) > 0:
        mane_parquet = OUTPUT_DIR / "mane_select_missense_variants.parquet"
        mane_csv = OUTPUT_DIR / "mane_select_missense_variants.csv.gz"
        
        mane.to_parquet(mane_parquet, index=False)
        mane.to_csv(mane_csv, index=False, compression='gzip')
        
        logger.info(f"\n✅ MANE Select only (clinical standard):")
        logger.info(f"  {mane_parquet} ({mane_parquet.stat().st_size / 1e6:.1f} MB)")
        logger.info(f"  {mane_csv} ({mane_csv.stat().st_size / 1e6:.1f} MB)")
        logger.info(f"  Variants: {len(mane):,}")
    
    # 4. MANE Select + Verified (gold standard)
    gold = df[(df['is_mane_select'] == True) & (df['seq_verified'] == True)].copy()
    if len(gold) > 0:
        gold_parquet = OUTPUT_DIR / "gold_standard_missense_variants.parquet"
        gold_csv = OUTPUT_DIR / "gold_standard_missense_variants.csv.gz"
        
        gold.to_parquet(gold_parquet, index=False)
        gold.to_csv(gold_csv, index=False, compression='gzip')
        
        logger.info(f"\n✅ Gold standard (MANE Select + Sequence Verified):")
        logger.info(f"  {gold_parquet} ({gold_parquet.stat().st_size / 1e6:.1f} MB)")
        logger.info(f"  {gold_csv} ({gold_csv.stat().st_size / 1e6:.1f} MB)")
        logger.info(f"  Variants: {len(gold):,}")
    
    # 5. Canonical transcripts
    if 'CANONICAL' in df.columns:
        canonical = df[df['CANONICAL'] == 'YES'].copy()
        if len(canonical) > 0:
            canonical_parquet = OUTPUT_DIR / "canonical_missense_variants.parquet"
            canonical.to_parquet(canonical_parquet, index=False)
            
            logger.info(f"\n✅ Canonical transcripts:")
            logger.info(f"  {canonical_parquet} ({canonical_parquet.stat().st_size / 1e6:.1f} MB)")
            logger.info(f"  Variants: {len(canonical):,}")


# ============================================================================
# MAIN
# ============================================================================

def main():
    """Main execution."""
    logger.info("\n" + "="*80)
    logger.info("GNOMAD MISSENSE VARIANT PROCESSOR")
    logger.info("="*80)
    logger.info(f"\nInput:  {INPUT_DIR}")
    logger.info(f"Output: {OUTPUT_DIR}")
    
    start_time = time.time()
    
    try:
        # Step 1: Load reference data
        initialize_reference_data()
        
        # Step 2: Process all chromosomes
        df = process_all_chromosomes()
        
        if df.empty:
            logger.error("\n❌ No data to process!")
            return False
        
        # Step 3: Generate statistics
        stats = generate_statistics(df)
        
        # Step 4: Save results
        save_results(df, stats)
        
        elapsed = time.time() - start_time
        logger.info("\n" + "="*80)
        logger.info("✅ PIPELINE COMPLETE")
        logger.info("="*80)
        logger.info(f"\nTotal time: {elapsed/60:.1f} minutes")
        logger.info(f"Total variants: {len(df):,}")
        logger.info(f"Verified variants: {stats['sequence_verified']:,}")
        
        return True
        
    except Exception as e:
        logger.error(f"\n❌ PIPELINE FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == "__main__":
    success = main()
    exit(0 if success else 1)