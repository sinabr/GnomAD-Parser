from variant_utils.utils import read_external_config
from pathlib import Path
import subprocess
from datetime import datetime
import pandas as pd
from pysam import VariantFile
from typing import List, Optional
import tempfile
import shutil
import logging

logger = logging.getLogger(__name__)


def queryGnomAD(
    assembly: str,
    CHROM: str,
    START: int,
    STOP: int,
    HGNC_ID: str,
    external_config_filepath: str,
    **kwargs
) -> pd.DataFrame:
    """
    Query gnomAD for missense variants in a gene; if assembly is 'GRCh37' gnomAD v2.1.1 is used, 
    otherwise gnomAD v4.1 is used
    
    IMPROVEMENTS:
    - Parallel GATK execution for exomes/genomes
    - Caching support to avoid re-querying same regions
    - Better error handling and validation
    - Optional cleanup of intermediate files
    - Better memory management
    
    Dependencies:
    - GATK
    - picard
    - java
    - gnomAD exomes and genomes VCF files

    Parameters:
    -----------
    assembly : str
        The genome assembly to use for the query, either 'GRCh37' or 'GRCh38'
    CHROM : str
        The chromosome for which to query gnomAD
    START : int
        The minimum position in the chromosome for which to query gnomAD
    STOP : int
        The maximum position in the chromosome for which to query gnomAD
    HGNC_ID : str
        The HGNC ID for filtering variants to the gene
    external_config_filepath : str
        Path to external tools configuration file

    Optional kwargs:
    - write_dir: str : Path to the directory where the output files will be written : default "/tmp"
    - use_cache: bool : Whether to use cached results if available : default False
    - cleanup: bool : Whether to delete intermediate files : default True
    - parallel: bool : Whether to run exomes/genomes queries in parallel : default False
    - gene_symbol: str : Gene symbol for cache naming (optional, uses HGNC_ID if not provided)

    Returns:
    - gene_df: pd.DataFrame : A DataFrame containing parsed VEP annotations for matched 
                              missense variants in gnomAD exomes and genomes
    """
    # Read configuration
    external_tools = read_external_config(external_config_filepath)
    
    # Parse kwargs
    write_dir = Path(kwargs.get("write_dir", "/tmp"))
    use_cache = kwargs.get("use_cache", False)
    cleanup = kwargs.get("cleanup", True)
    parallel = kwargs.get("parallel", False)
    gene_symbol = kwargs.get("gene_symbol", HGNC_ID)
    
    write_dir.mkdir(exist_ok=True)
    
    # Check for cached results
    if use_cache:
        cache_file = write_dir / f"gnomAD_{gene_symbol}_{assembly}_{CHROM}_{START}_{STOP}.parquet"
        if cache_file.exists():
            logger.info(f"Using cached gnomAD data from {cache_file}")
            return pd.read_parquet(cache_file)
    
    # Validate tools
    java = Path(external_tools.get("java"))
    picard_filepath = Path(external_tools.get("picard_filepath"))
    gatk = external_tools.get("gatk")
    
    if not picard_filepath.exists():
        raise FileNotFoundError(f"picard_filepath does not exist: {picard_filepath}")
    
    # Setup version-specific parameters
    release_version = "v4.1" if assembly == "GRCh38" else "r2.1.1"
    chr_prefix = "chr" if release_version == "v4.1" else ""
    
    gnomad_vcf_root = Path(
        external_tools['gnomad_v4_vcf_root'] if release_version == "v4.1" 
        else external_tools['gnomad_v2_vcf_root']
    )
    
    if not gnomad_vcf_root.exists():
        raise FileNotFoundError(f"gnomad_vcf_root does not exist: {gnomad_vcf_root}")
    
    # Setup VCF paths
    gnomAD_exomes_filepath = (
        gnomad_vcf_root / f"exomes/gnomad.exomes.{release_version}.sites.{chr_prefix}{CHROM}.vcf.bgz"
    )
    gnomAD_genomes_filepath = (
        gnomad_vcf_root / f"genomes/gnomad.genomes.{release_version}.sites.{chr_prefix}{CHROM}.vcf.bgz"
    )
    
    # Validate input VCFs exist
    if not gnomAD_exomes_filepath.exists():
        raise FileNotFoundError(f"Exomes VCF not found: {gnomAD_exomes_filepath}")
    if not gnomAD_genomes_filepath.exists():
        raise FileNotFoundError(f"Genomes VCF not found: {gnomAD_genomes_filepath}")
    
    # Generate output filenames
    timestamp = str(datetime.now()).replace(' ', '_')
    exomes_output_file = write_dir / f"selectvariants_{timestamp}.exomes.vcf"
    genomes_output_file = write_dir / f"selectvariants_{timestamp}.genomes.vcf"
    combined_output_file = write_dir / f"combinevariants_{timestamp}.vcf"
    tsv_output_file = write_dir / f"combinevariants_{timestamp}.tsv"
    
    try:
        # Step 1: Extract variants from exomes and genomes
        if parallel:
            # Run in parallel using subprocess.Popen
            exomes_proc = run_select_variants_async(
                gatk, gnomAD_exomes_filepath, chr_prefix, CHROM, START, STOP, exomes_output_file
            )
            genomes_proc = run_select_variants_async(
                gatk, gnomAD_genomes_filepath, chr_prefix, CHROM, START, STOP, genomes_output_file
            )
            
            # Wait for both to complete
            exomes_proc.wait()
            genomes_proc.wait()
            
            if exomes_proc.returncode != 0:
                raise RuntimeError(f"Exomes SelectVariants failed with code {exomes_proc.returncode}")
            if genomes_proc.returncode != 0:
                raise RuntimeError(f"Genomes SelectVariants failed with code {genomes_proc.returncode}")
        else:
            # Run sequentially
            run_select_variants(
                gatk, gnomAD_exomes_filepath, chr_prefix, CHROM, START, STOP, exomes_output_file
            )
            run_select_variants(
                gatk, gnomAD_genomes_filepath, chr_prefix, CHROM, START, STOP, genomes_output_file
            )
        
        # Step 2: Merge VCFs
        run_merge_vcfs(java, picard_filepath, exomes_output_file, genomes_output_file, combined_output_file)
        
        # Step 3: Convert to table
        run_variants_to_table(gatk, combined_output_file, tsv_output_file)
        
        # Step 4: Parse results
        gnomAD_df = pd.read_csv(tsv_output_file, delimiter='\t')
        
        # Step 5: Parse VEP annotations
        vep_columns = get_vep_columns_from_vcf_header(str(combined_output_file))
        vep_df = parse_vep(gnomAD_df, columns=vep_columns)
        
        # Step 6: Merge with main dataframe
        gnomAD_df = pd.merge(
            gnomAD_df, vep_df, 
            left_index=True, right_on='index', 
            validate='one_to_many'
        )
        
        # Step 7: Filter for the specific gene
        gene_df = gnomAD_df[gnomAD_df.HGNC_ID == HGNC_ID].copy()
        
        # Step 8: Normalize columns
        gene_df = gene_df.assign(
            CHROM=gene_df.CHROM.astype(str).str.replace("chr", ""),
            POS=gene_df.POS.astype(str),
            REF=gene_df.REF.astype(str),
            ALT=gene_df.ALT.astype(str)
        )
        
        # Step 9: Cache results if requested
        if use_cache:
            cache_file = write_dir / f"gnomAD_{gene_symbol}_{assembly}_{CHROM}_{START}_{STOP}.parquet"
            gene_df.to_parquet(cache_file, index=False)
            logger.info(f"Cached gnomAD data to {cache_file}")
        
        return gene_df
        
    finally:
        # Cleanup intermediate files if requested
        if cleanup:
            for f in [exomes_output_file, genomes_output_file, 
                     combined_output_file, tsv_output_file,
                     Path(str(combined_output_file) + ".idx")]:
                if f.exists():
                    f.unlink()


def run_select_variants(
    gatk: str, 
    input_vcf: Path, 
    chr_prefix: str, 
    chrom: str, 
    start: int, 
    stop: int, 
    output_file: Path
) -> None:
    """Run GATK SelectVariants synchronously."""
    cmd = [
        gatk, "SelectVariants",
        "-V", str(input_vcf),
        "-L", f"{chr_prefix}{chrom}:{start}-{stop}",
        "--select-type-to-include", "SNP",
        "--exclude-filtered",
        "--output", str(output_file)
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"SelectVariants failed:\n{result.stderr}")


def run_select_variants_async(
    gatk: str,
    input_vcf: Path,
    chr_prefix: str,
    chrom: str,
    start: int,
    stop: int,
    output_file: Path
) -> subprocess.Popen:
    """Run GATK SelectVariants asynchronously."""
    cmd = [
        gatk, "SelectVariants",
        "-V", str(input_vcf),
        "-L", f"{chr_prefix}{chrom}:{start}-{stop}",
        "--select-type-to-include", "SNP",
        "--exclude-filtered",
        "--output", str(output_file)
    ]
    
    return subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )


def run_merge_vcfs(
    java: Path,
    picard_filepath: Path,
    exomes_vcf: Path,
    genomes_vcf: Path,
    output_file: Path
) -> None:
    """Run Picard MergeVcfs."""
    cmd = [
        str(java), "-jar", str(picard_filepath),
        "MergeVcfs",
        f"I={exomes_vcf}",
        f"I={genomes_vcf}",
        f"O={output_file}"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"MergeVcfs failed:\n{result.stderr}")


def run_variants_to_table(
    gatk: str,
    input_vcf: Path,
    output_tsv: Path
) -> None:
    """Run GATK VariantsToTable."""
    cmd = [
        gatk, "VariantsToTable",
        "-V", str(input_vcf),
        "-F", "CHROM",
        "-F", "POS",
        "-F", "ID",
        "-F", "REF",
        "-F", "ALT",
        "-F", "QUAL",
        "-F", "FILTER",
        "-ASF", "AC",
        "-ASF", "AF",
        "-ASF", "vep",
        "-O", str(output_tsv)
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"VariantsToTable failed:\n{result.stderr}")


def get_vep_columns_from_vcf_header(vcf_file: str) -> List[str]:
    """
    Read a gnomAD vcf file and extract the VEP columns from the header

    Parameters:
    -----------
    vcf_file : str
        The path to the gnomAD VCF file

    Returns:
    --------
    list : A list of VEP columns
    """
    vcf = VariantFile(vcf_file)
    return vcf.header.info['vep'].description.split("Format: ")[1].split("|")


def parse_vep(df: pd.DataFrame, columns: List[str]) -> pd.DataFrame:
    """
    Parse the 'vep' column of the gnomAD dataframe into its own dataframe

    Parameters:
    -----------
    df : pd.DataFrame
        The gnomAD dataframe
    columns : list
        The VEP columns
    
    Returns:
    --------
    pd.DataFrame : A DataFrame containing the parsed VEP annotations
    """
    vep_series = df.vep.apply(
        lambda r: list(map(lambda s: dict(zip(columns, s.split('|'))), r.split(",")))
    )
    vep_df = pd.DataFrame(vep_series, index=df.index).explode('vep')
    vep_df = pd.DataFrame.from_records(vep_df.vep.values, index=vep_df.index).reset_index()
    return vep_df


def query_gnomad_batch(
    genes: List[dict],
    assembly: str,
    external_config_filepath: str,
    **kwargs
) -> pd.DataFrame:
    """
    Query gnomAD for multiple genes in batch.
    
    This is more efficient than calling queryGnomAD repeatedly because it:
    - Reuses the same write directory
    - Can cache intermediate results
    - Provides progress tracking
    
    Parameters:
    -----------
    genes : List[dict]
        List of gene dictionaries, each containing:
        - 'symbol': gene symbol
        - 'chrom': chromosome
        - 'start': start position
        - 'stop': stop position
        - 'hgnc_id': HGNC ID
    assembly : str
        'GRCh37' or 'GRCh38'
    external_config_filepath : str
        Path to external tools config
    **kwargs : additional arguments passed to queryGnomAD
    
    Returns:
    --------
    pd.DataFrame : Combined results for all genes
    """
    all_results = []
    
    for i, gene in enumerate(genes, 1):
        logger.info(f"Processing gene {i}/{len(genes)}: {gene['symbol']}")
        
        try:
            result = queryGnomAD(
                assembly=assembly,
                CHROM=gene['chrom'],
                START=gene['start'],
                STOP=gene['stop'],
                HGNC_ID=gene['hgnc_id'],
                external_config_filepath=external_config_filepath,
                gene_symbol=gene['symbol'],
                **kwargs
            )
            
            # Add gene identifier
            result['query_gene'] = gene['symbol']
            all_results.append(result)
            
        except Exception as e:
            logger.error(f"Failed to process {gene['symbol']}: {e}")
            continue
    
    if not all_results:
        raise RuntimeError("No genes processed successfully")
    
    return pd.concat(all_results, ignore_index=True)