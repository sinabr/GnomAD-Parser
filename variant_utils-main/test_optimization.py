#!/usr/bin/env python3
"""
Quick test to verify optimization changes work correctly.
Tests the new pysam-based extraction against a known gene.
"""

import sys
from pathlib import Path
import pandas as pd
from variant_utils.gnomad_utils import extract_variants_fast, get_vep_columns_from_vcf_header, parse_vep

def test_pysam_extraction():
    """Test that pysam extraction works."""

    print("Testing pysam-based variant extraction...")
    print("="*80)

    # Test with a small chromosome region
    # You'll need to adjust these paths based on your setup
    test_vcf = Path("/projects/lugoteam/protein_graphs/GnomAD-Parser/data/exomes/gnomad.exomes.v4.1.sites.chr22.vcf.bgz")

    if not test_vcf.exists():
        print(f"❌ Test VCF not found: {test_vcf}")
        print("Please run this script on the cluster where gnomAD data is available.")
        return False

    print(f"✅ Test VCF found: {test_vcf}")

    # Test extraction on a small region (chr22, first 1MB)
    print("\nExtracting variants from chr22:1-1000000...")

    try:
        df = extract_variants_fast(
            vcf_path=test_vcf,
            chrom="22",
            start=1,
            end=1000000,
            chr_prefix="chr"
        )

        print(f"✅ Extracted {len(df)} variants")

        if len(df) > 0:
            print("\nSample variant:")
            print(df.iloc[0].to_dict())

            # Test VEP parsing
            print("\nTesting VEP column extraction...")
            vep_columns = get_vep_columns_from_vcf_header(str(test_vcf))
            print(f"✅ Found {len(vep_columns)} VEP columns")

            print("\nTesting VEP parsing...")
            vep_df = parse_vep(df, columns=vep_columns)
            print(f"✅ Parsed {len(vep_df)} VEP annotations")

            if len(vep_df) > 0:
                print("\nSample VEP annotation:")
                print(vep_df.iloc[0].to_dict())

        print("\n" + "="*80)
        print("✅ ALL TESTS PASSED")
        print("="*80)
        return True

    except Exception as e:
        print(f"\n❌ ERROR: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == "__main__":
    success = test_pysam_extraction()
    sys.exit(0 if success else 1)
