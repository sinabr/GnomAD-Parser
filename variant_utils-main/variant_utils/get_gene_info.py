import pandas as pd
from pathlib import Path
from typing import Optional

# Global cache to avoid re-reading MANE file for every gene
_MANE_CACHE = None

def get_gene_info(
    geneSymbol: str,
    mane_summary_file: Optional[str] = None
) -> dict:
    """
    Get information about a gene from the LOCAL MANE summary file (offline).
    
    This version uses a local MANE file instead of downloading.
    The MANE data is cached in memory after first load for fast lookups.

    Parameters
    ----------
    geneSymbol : str
        Gene symbol (e.g., "HBB", "TP53")
    mane_summary_file : str, optional
        Path to local MANE summary file. If None, uses default location.

    Returns
    -------
    dict
        NCBI_GeneID
        Ensembl_Gene
        HGNC_ID
        symbol
        name
        RefSeq_nuc
        RefSeq_prot
        Ensembl_nuc
        Ensembl_prot
        MANE_status
        GRCh38_chr
        chr_start
        chr_end
        chr_strand
        CHROM (parsed chromosome: 1-22, X, Y)
    """
    global _MANE_CACHE
    
    # Default to local MANE file
    if mane_summary_file is None:
        mane_summary_file = "/jet/home/barazand/NEWOCEAN/ref_data/mane/MANE.GRCh38.v1.3.summary.txt.gz"
    
    mane_path = Path(mane_summary_file)
    
    if not mane_path.exists():
        raise FileNotFoundError(
            f"MANE summary file not found: {mane_path}\n"
            f"Please ensure the file exists or provide correct path."
        )
    
    # Load MANE data into cache (only once)
    if _MANE_CACHE is None:
        _MANE_CACHE = pd.read_csv(mane_path, sep="\t", compression='gzip')
    
    # Look up gene
    record = _MANE_CACHE[_MANE_CACHE.symbol == geneSymbol]
    
    if record.empty:
        raise ValueError(f"Gene {geneSymbol} not found in MANE summary file")
    
    record = record.iloc[0]
    
    # Parse chromosome from NCBI format
    CHROM = record['GRCh38_chr']
    
    try:
        # Handle NC_000001.11 -> 1, NC_000023.11 -> X, etc.
        chrom_num = int(CHROM.split(".")[0].replace("NC_", "").replace("000", ""))
        
        if chrom_num == 23:
            chrom = "X"
        elif chrom_num == 24:
            chrom = "Y"
        else:
            chrom = str(chrom_num)
    except (ValueError, AttributeError):
        # Handle alternate contigs (NW_, NT_, etc.) or other formats
        # These will fail downstream in gnomAD query, which is expected
        raise ValueError(
            f"Gene {geneSymbol} is on alternate contig {CHROM} - "
            f"not supported for gnomAD queries"
        )
    
    # Add parsed chromosome to record
    record = record.to_dict()
    record['CHROM'] = chrom
    
    return record