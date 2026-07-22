#!/usr/bin/env python3
"""
Annotate ClinVar variants with gnomAD v4.1 JOINT allele frequencies

This script queries gnomAD v4.1 JOINT VCF files to extract allele counts
and frequencies for all ancestries for ClinVar P/LP and VUS variants.

Extracted from: Reproducible_Code_Appendix.Rmd
Updated: March 20, 2026

Requirements:
- pysam (pip install pysam)
- gnomAD v4.1 JOINT VCF files downloaded locally
"""

import pandas as pd
import pysam
import argparse
from pathlib import Path
import sys

def parse_args():
    """Parse command-line arguments"""
    parser = argparse.ArgumentParser(
        description='Annotate variants with gnomAD v4.1 frequencies'
    )
    parser.add_argument(
        '--variants',
        required=True,
        help='Input variant file (CSV with Chromosome, Position, Ref, Alt columns)'
    )
    parser.add_argument(
        '--gnomad-dir',
        required=True,
        help='Directory containing gnomAD v4.1 JOINT VCF files'
    )
    parser.add_argument(
        '--output',
        default='variants_annotated_gnomAD.csv',
        help='Output file with gnomAD annotations'
    )
    parser.add_argument(
        '--chromosomes',
        default='1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X',
        help='Comma-separated list of chromosomes to process'
    )
    
    return parser.parse_args()

def create_variant_lookup(variants_df):
    """Create lookup set of variant keys"""
    variants_df['lookup_key'] = (
        variants_df['Chromosome'].astype(str) + ':' +
        variants_df['Position'].astype(str) + ':' +
        variants_df['Ref'] + ':' +
        variants_df['Alt']
    )
    
    return set(variants_df['lookup_key']), variants_df

def extract_gnomad_frequencies(gnomad_dir, chromosomes, variant_lookup):
    """
    Extract allele frequencies from gnomAD VCF files
    
    Parameters:
    -----------
    gnomad_dir : str
        Directory containing gnomAD VCF files
    chromosomes : list
        List of chromosomes to process
    variant_lookup : set
        Set of variant keys to extract
    
    Returns:
    --------
    DataFrame with gnomAD annotations
    """
    
    gnomad_data = []
    total_variants_found = 0
    
    for chrom in chromosomes:
        print(f"\nProcessing chromosome {chrom}...")
        
        # gnomAD v4.1 JOINT filename pattern
        vcf_file = Path(gnomad_dir) / f'gnomad.genomes_exomes.v4.1.sites.chr{chrom}.vcf.bgz'
        
        if not vcf_file.exists():
            print(f"  ⚠️  File not found: {vcf_file}")
            print(f"      Skipping chromosome {chrom}")
            continue
        
        try:
            vcf = pysam.VariantFile(str(vcf_file))
        except Exception as e:
            print(f"  ❌ Error opening VCF: {e}")
            continue
        
        chr_variants = 0
        
        for record in vcf:
            # Create lookup key
            # gnomAD uses 'chr' prefix
            chrom_clean = record.contig.replace('chr', '')
            key = f"{chrom_clean}:{record.pos}:{record.ref}:{record.alts[0]}"
            
            if key not in variant_lookup:
                continue
            
            chr_variants += 1
            
            # Extract INFO fields
            info = record.info
            
            variant_data = {
                'Chromosome': chrom_clean,
                'Position': record.pos,
                'Ref': record.ref,
                'Alt': record.alts[0],
                'lookup_key': key,
                
                # Total (all ancestries combined)
                'gnomAD_AC_total': info.get('AC', [0])[0] if 'AC' in info else 0,
                'gnomAD_AN_total': info.get('AN', 0),
                'gnomAD_AF_total': info.get('AF', [0])[0] if 'AF' in info else 0.0,
                'gnomAD_nhomalt_total': info.get('nhomalt', [0])[0] if 'nhomalt' in info else 0,
                
                # African/African American
                'gnomAD_AC_afr': info.get('AC_afr', [0])[0] if 'AC_afr' in info else 0,
                'gnomAD_AN_afr': info.get('AN_afr', 0),
                'gnomAD_AF_afr': info.get('AF_afr', [0])[0] if 'AF_afr' in info else 0.0,
                
                # Latino/Admixed American
                'gnomAD_AC_amr': info.get('AC_amr', [0])[0] if 'AC_amr' in info else 0,
                'gnomAD_AN_amr': info.get('AN_amr', 0),
                'gnomAD_AF_amr': info.get('AF_amr', [0])[0] if 'AF_amr' in info else 0.0,
                
                # Ashkenazi Jewish
                'gnomAD_AC_asj': info.get('AC_asj', [0])[0] if 'AC_asj' in info else 0,
                'gnomAD_AN_asj': info.get('AN_asj', 0),
                'gnomAD_AF_asj': info.get('AF_asj', [0])[0] if 'AF_asj' in info else 0.0,
                
                # East Asian
                'gnomAD_AC_eas': info.get('AC_eas', [0])[0] if 'AC_eas' in info else 0,
                'gnomAD_AN_eas': info.get('AN_eas', 0),
                'gnomAD_AF_eas': info.get('AF_eas', [0])[0] if 'AF_eas' in info else 0.0,
                
                # Finnish
                'gnomAD_AC_fin': info.get('AC_fin', [0])[0] if 'AC_fin' in info else 0,
                'gnomAD_AN_fin': info.get('AN_fin', 0),
                'gnomAD_AF_fin': info.get('AF_fin', [0])[0] if 'AF_fin' in info else 0.0,
                
                # Non-Finnish European (NFE)
                'gnomAD_AC_nfe': info.get('AC_nfe', [0])[0] if 'AC_nfe' in info else 0,
                'gnomAD_AN_nfe': info.get('AN_nfe', 0),
                'gnomAD_AF_nfe': info.get('AF_nfe', [0])[0] if 'AF_nfe' in info else 0.0,
                
                # Remaining (other)
                'gnomAD_AC_remaining': info.get('AC_remaining', [0])[0] if 'AC_remaining' in info else 0,
                'gnomAD_AN_remaining': info.get('AN_remaining', 0),
                'gnomAD_AF_remaining': info.get('AF_remaining', [0])[0] if 'AF_remaining' in info else 0.0,
            }
            
            gnomad_data.append(variant_data)
        
        vcf.close()
        
        total_variants_found += chr_variants
        print(f"  ✅ Found {chr_variants} variants on chr{chrom}")
    
    print(f"\n✅ Total variants found in gnomAD: {total_variants_found:,}")
    
    return pd.DataFrame(gnomad_data)

def main():
    """Main execution"""
    args = parse_args()
    
    print("="*80)
    print("GNOMAD v4.1 JOINT ANNOTATION")
    print("="*80)
    
    # Load variant file
    print(f"\nLoading variants from: {args.variants}")
    
    if args.variants.endswith('.xlsx'):
        variants_df = pd.read_excel(args.variants)
    else:
        variants_df = pd.read_csv(args.variants)
    
    print(f"✅ Loaded {len(variants_df):,} variants")
    
    # Verify required columns
    required_cols = ['Chromosome', 'Position', 'Ref', 'Alt']
    missing_cols = [col for col in required_cols if col not in variants_df.columns]
    
    # Try alternative column names
    if 'Reference' in variants_df.columns:
        variants_df = variants_df.rename(columns={'Reference': 'Ref'})
    if 'Alternate' in variants_df.columns:
        variants_df = variants_df.rename(columns={'Alternate': 'Alt'})
    
    if missing_cols:
        print(f"❌ Error: Missing required columns: {missing_cols}")
        print(f"   Available columns: {list(variants_df.columns)}")
        sys.exit(1)
    
    # Create lookup
    variant_lookup, variants_df = create_variant_lookup(variants_df)
    print(f"✅ Created lookup for {len(variant_lookup):,} unique variants")
    
    # Parse chromosomes
    chromosomes = args.chromosomes.split(',')
    print(f"✅ Will process {len(chromosomes)} chromosomes: {', '.join(chromosomes)}")
    
    # Extract gnomAD data
    print(f"\nExtracting data from gnomAD directory: {args.gnomad_dir}")
    gnomad_df = extract_gnomad_frequencies(
        args.gnomad_dir,
        chromosomes,
        variant_lookup
    )
    
    if len(gnomad_df) == 0:
        print("\n⚠️  Warning: No variants found in gnomAD!")
        print("   Check that:")
        print("   1. gnomAD VCF files are in the correct directory")
        print("   2. Chromosome names match (with/without 'chr' prefix)")
        print("   3. Variant coordinates match exactly")
    
    # Merge with input variants
    print(f"\nMerging with input variants...")
    
    final_data = variants_df.merge(
        gnomad_df,
        on='lookup_key',
        how='left',
        suffixes=('', '_gnomad')
    )
    
    # Fill missing gnomAD data with 0/NaN
    gnomad_cols = [col for col in final_data.columns if col.startswith('gnomAD_')]
    for col in gnomad_cols:
        if 'AF' in col:
            final_data[col] = final_data[col].fillna(0.0)
        elif 'AC' in col or 'AN' in col or 'nhomalt' in col:
            final_data[col] = final_data[col].fillna(0)
    
    # Drop duplicate coordinate columns
    cols_to_drop = [col for col in final_data.columns if col.endswith('_gnomad')]
    final_data = final_data.drop(columns=cols_to_drop)
    
    # Save
    print(f"\n💾 Saving annotated variants to: {args.output}")
    final_data.to_csv(args.output, index=False)
    
    # Summary
    print("\n" + "="*80)
    print("ANNOTATION SUMMARY")
    print("="*80)
    print(f"Total input variants: {len(variants_df):,}")
    print(f"Variants found in gnomAD: {len(gnomad_df):,}")
    print(f"Variants not in gnomAD: {len(variants_df) - len(gnomad_df):,}")
    print(f"Coverage: {len(gnomad_df)/len(variants_df)*100:.1f}%")
    
    # Ancestry-specific coverage
    if len(gnomad_df) > 0:
        print(f"\ngnomAD v4.1 JOINT sample sizes:")
        print(f"  Total: {gnomad_df['gnomAD_AN_total'].max():,} alleles")
        print(f"  AFR: {gnomad_df['gnomAD_AN_afr'].max():,} alleles")
        print(f"  AMR: {gnomad_df['gnomAD_AN_amr'].max():,} alleles")
        print(f"  ASJ: {gnomad_df['gnomAD_AN_asj'].max():,} alleles")
        print(f"  EAS: {gnomad_df['gnomAD_AN_eas'].max():,} alleles")
        print(f"  FIN: {gnomad_df['gnomAD_AN_fin'].max():,} alleles")
        print(f"  NFE: {gnomad_df['gnomAD_AN_nfe'].max():,} alleles")
    
    print("\n✅ ANNOTATION COMPLETE!")
    print("="*80)

if __name__ == '__main__':
    main()

# Usage:
# python annotate_with_gnomad_STANDALONE.py \
#     --variants ClinVar_PLP_2star.csv \
#     --gnomad-dir /path/to/gnomad/vcf/files \
#     --output ClinVar_PLP_with_gnomAD.csv
#
# Note: gnomAD v4.1 JOINT files can be downloaded from:
# https://gnomad.broadinstitute.org/downloads#v4-variants
