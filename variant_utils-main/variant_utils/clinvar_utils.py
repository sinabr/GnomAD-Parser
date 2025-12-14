from variant_utils.utils import read_external_config
from pathlib import Path
import subprocess
from typing import Optional
import pandas as pd

def annotate_variants(df):
    low_qual_rev_stat = {'no_assertion_criteria_provided','no_classification_provided'}
    df = df.assign(CHROM=df.CHROM.astype(str),
                   POS=df.POS.astype(str),
                   REF=df.REF.astype(str),
                   ALT=df.ALT.astype(str),
                   is_pathogenic=(df.CLNSIG.isin({"Pathogenic","Likely_pathogenic","Pathogenic/Likely_pathogenic"})) & (~df.CLNREVSTAT.isin(low_qual_rev_stat)),
                   is_benign=(df.CLNSIG.isin({"Benign","Likely_benign","Benign/Likely_benign"})) & (~df.CLNREVSTAT.isin(low_qual_rev_stat)),
                   is_conflicting=(df.CLNSIG == "Conflicting_classifications_of_pathogenicity") & (~df.CLNREVSTAT.isin(low_qual_rev_stat)),
                   is_VUS=(df.CLNSIG == "Uncertain_significance") & (~df.CLNREVSTAT.isin(low_qual_rev_stat)))
    return df

def queryClinVarVCF(clinvar_filepath:str|Path, CHROM : str,START: int,STOP : int,external_tools_filepath:str|Path,gene_name:Optional[str]=None, gene_id:Optional[str]=None,**kwargs)->pd.DataFrame:
    """
    Query ClinVar for variants in a region of the genome

    This function uses GATK SelectVariants to extract variants in a region of the genome from the ClinVar VCF file
    ensure it is installed the path is specified in the external tools configuration file


    Parameters:
    -----------
    clinvar_filepath : str
        The path to the ClinVar VCF file
        downloadable from https://ftp.ncbi.nlm.nih.gov/pub/clinvar/
    CHROM : str
        The chromosome for which to query ClinVar
    START : int
        The minimum position in the chromosome for which to query ClinVar
    STOP : int
        The maximum position in the chromosome for which to query ClinVar
    external_tools_filepath : str|Path
        The path to the external tools configuration file

    Optional Parameters:
    --------------------
    gene_name : str|None
        The gene name for which to filter on
    gene_id : str|None
        NCBI Gene ID for which to filter on

    Returns:
    --------
    clinVar_df : pd.DataFrame
        A DataFrame containing the ClinVar variants in the specified region
    """
    filepath = Path(clinvar_filepath)
    assert filepath.exists(), "filepath {} does not exist".format(filepath)
    write_dir = Path(kwargs.get("write_dir","/tmp"))
    write_dir.mkdir(exist_ok=True)
    external_tools = read_external_config(external_tools_filepath)
    output_file = write_dir / f"ClinVar_selectvariants_chr{CHROM}:{START}-{STOP}.vcf"
    cmd = f"{external_tools['gatk']} SelectVariants -V {filepath} -L {CHROM}:{START}-{STOP} --exclude-filtered --output {output_file}"
    subprocess.run(cmd.split(" "))

    tsvout = str(output_file).replace('.vcf','.tsv')
    variants2table = f"{external_tools['gatk']} VariantsToTable -V {output_file} -O {tsvout}"
    subprocess.run(variants2table.split(" "))

    clinVar_df = pd.read_csv(tsvout,delimiter='\t')
    if gene_name is not None and gene_id is not None:
        clinVar_df = clinVar_df[clinVar_df.GENEINFO == f"{gene_name}:{gene_id}"]

    clinVar_df = annotate_variants(clinVar_df)
    clinVar_df.rename(columns={"#CHROM":'CHROM'},inplace=True)
    return clinVar_df