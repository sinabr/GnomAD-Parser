# Variant Utils
> Utilities for querying external variant databases

## Required External Tools
- [gatk](https://github.com/broadinstitute/gatk/releases) must be installed
- [Picard](https://broadinstitute.github.io/picard/) must be installed
- java (see Picard install instructions)

## Required External Data
- gnomAD v4 (if using GRCh38 mappings)
    directory structure:
    ```
    /data/dbs/gnomad/release/v4.1.0/
    ├── exomes
    │   ├── gnomad.exomes.v4.1.sites.chr10.vcf.bgz
    │   ├── gnomad.exomes.v4.1.sites.chr10.vcf.bgz.tbi
    │   ├── gnomad.exomes.v4.1.sites.chr11.vcf.bgz
    │   ├── gnomad.exomes.v4.1.sites.chr11.vcf.bgz.tbi
    │   ├── gnomad.exomes.v4.1.sites.chr12.vcf.bgz
    │   ├── gnomad.exomes.v4.1.sites.chr12.vcf.bgz.tbi
    │   ├── gnomad.exomes.v4.1.sites.chr13.vcf.bgz
    │   ├── gnomad.exomes.v4.1.sites.chr13.vcf.bgz.tbi
    │   ├── gnomad.exomes.v4.1.sites.chr14.vcf.bgz
    │   ├── gnomad.exomes.v4.1.sites.chr14.vcf.bgz.tbi
    │   ├── gnomad.exomes.v4.1.sites.chr15.vcf.bgz
    │   ├── gnomad.exomes.v4.1.sites.chr15.vcf.bgz.tbi
    │   .
        .
        .
    ├── genomes
    │   ├── gnomad.genomes.v4.1.sites.chr10.vcf.bgz
    │   ├── gnomad.genomes.v4.1.sites.chr10.vcf.bgz.tbi
    │   ├── gnomad.genomes.v4.1.sites.chr11.vcf.bgz
    │   ├── gnomad.genomes.v4.1.sites.chr11.vcf.bgz.tbi
    │   ├── gnomad.genomes.v4.1.sites.chr12.vcf.bgz
    │   ├── gnomad.genomes.v4.1.sites.chr12.vcf.bgz.tbi
    │   ├── gnomad.genomes.v4.1.sites.chr13.vcf.bgz
    │   ├── gnomad.genomes.v4.1.sites.chr13.vcf.bgz.tbi
    │   ├── gnomad.genomes.v4.1.sites.chr14.vcf.bgz
    │   ├── gnomad.genomes.v4.1.sites.chr14.vcf.bgz.tbi
    │   ├── gnomad.genomes.v4.1.sites.chr15.vcf.bgz
    │   ├── gnomad.genomes.v4.1.sites.chr15.vcf.bgz.tbi
    │   .
        .
        .
    ```
- gnomAD v2 (if using GRCh37 mappings)
    ```
    /data/dbs/gnomad/release/v2.1.1/
    ├── exomes
    │   ├── gnomad.exomes.r2.1.1.sites.10.vcf.bgz
    │   ├── gnomad.exomes.r2.1.1.sites.10.vcf.bgz.tbi
    │   ├── gnomad.exomes.r2.1.1.sites.11.vcf.bgz
        .
        .
        .
    ├── genomes
    │   ├── gnomad.genomes.r2.1.1.sites.10.vcf.bgz
    │   ├── gnomad.genomes.r2.1.1.sites.10.vcf.bgz.tbi
    │   ├── gnomad.genomes.r2.1.1.sites.11.vcf.bgz
    │   ├── gnomad.genomes.r2.1.1.sites.11.vcf.bgz.tbi
    │   ├── gnomad.genomes.r2.1.1.sites.12.vcf.bgz
    │   .
        .
        .
- spliceAI scores (if interested in filtering splice variants)
    ```
    /data/dbs/spliceAI/
    ├── spliceai_scores.raw.snv.hg19.vcf.gz
    ├── spliceai_scores.raw.snv.hg19.vcf.gz.tbi
    ├── spliceai_scores.raw.snv.hg38.vcf.gz
    └── spliceai_scores.raw.snv.hg38.vcf.gz.tbi
    ```

## External Config File
Create `external_tools.json` with the following example structure:

```json
{
    "java" : "/usr/bin/java",
    "gatk" : "/path/to/gatk-4.4.0.0/gatk",
    "picard_filepath": "/path/to/picard/build/libs/picard.jar",
    "gnomad_v4_vcf_root": "/data/dbs/gnomad/release/v4.1.0/",
    "gnomad_v2_vcf_root": "/data/dbs/gnomad/release/v2.1.1/",
    "spliceAIRoot" : "/data/dbs/spliceAI/"
}
```

## Install
```bash
cd variant_utils
pip install -e .
```

## Usage
### Fetch information on a gene
```python
from variant_utils.get_gene_info import get_gene_info

brca1_info = get_gene_info("BRCA1")
```

### Fetch gnomAD variants for a gene
```python
from variant_utils.get_gene_info import get_gene_info
from variant_utils.gnomad_utils import queryGnomAD

brca1_info = get_gene_info("BRCA1")
brca1_gnomad_variants = queryGnomAD("GRCh38",brca1_info.CHROM, brca1_info.chr_start, brca1_info.chr_end, brca1_info.HGNC_ID,"external_tools.json")
```

### Fetch ClinVar variants for a gene
```python
from pathlib import Path
from variant_utils.get_gene_info import get_gene_info
from variant_utils.clinvar_utils import queryClinVarVCF
import urllib.request

brca1_info = get_gene_info("BRCA1")
# set destination to save/reload ClinVar
cache_dir = Path(".cache")
cache_dir.mkdir(exist_ok=True)
# download ClinVar release (e.g., 2018-12-17) if file is not present
clinvar_filepath = cache_dir / "clinvar_20181217.vcf.gz"
idx_filepath = cache_dir / "clinvar_20181217.vcf.gz.tbi"
if not clinvar_filepath.exists():
    urllib.request.urlretrieve(f"https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/weekly/clinvar_20181217.vcf.gz",str(clinvar_filepath))
if not idx_filepath.exists():
    urllib.request.urlretrieve(f"https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/weekly/clinvar_20181217.vcf.gz.tbi",str(idx_filepath))

brca1_clinvar_variants = queryClinVarVCF(str(clinvar_filepath), brca1_info.CHROM, brca1_info.chr_start, brca1_info.chr_end, "external_tools.json",write_dir=".cache")
```