"""
================================================================================
STRUCTURAL VARIANT COVERAGE ANALYSIS - gnomAD v4.1 JOINT
================================================================================

Analyzes how well gnomAD v4.1 JOINT covers structural variants (large deletions)
in ClinVar P/LP and VUS variants.

Input Files:
  - any_PLP_MASTER_JOINT.txt  (P/LP variants with gnomAD JOINT annotation)
  - VUS_2plus_stars_JOINT.txt (VUS variants with gnomAD JOINT annotation)

Analysis:
  - Identifies variants missing from gnomAD (gnomAD_AF_total is NA)
  - Classifies large deletions (>50bp) vs small variants
  - Calculates structural variant coverage gap

================================================================================
"""

import pandas as pd
import numpy as np

print("=" * 80)
print("STRUCTURAL VARIANT COVERAGE ANALYSIS - gnomAD v4.1 JOINT")
print("=" * 80)

# ============================================================================
# CONFIGURATION
# ============================================================================

LARGE_DELETION_THRESHOLD = 50  # bp

# ============================================================================
# FUNCTION: Analyze SV coverage
# ============================================================================

def analyze_sv_coverage(df, variant_type):
    """
    Analyze structural variant coverage in gnomAD
    
    Args:
        df: DataFrame with gnomAD annotation
        variant_type: 'P/LP' or 'VUS'
    """
    print(f"\n{'=' * 80}")
    print(f"{variant_type} ANALYSIS")
    print(f"{'=' * 80}")
    
    # Total variants
    total = len(df)
    print(f"\n📊 Total variants: {total:,}")
    
    # Check for gnomAD annotation
    af_col = 'gnomAD_AF_total'
    if af_col not in df.columns:
        print(f"❌ ERROR: {af_col} column not found")
        print(f"Available columns: {', '.join(df.columns[:10])}...")
        return
    
    # Convert AF to numeric
    df[af_col] = pd.to_numeric(df[af_col], errors='coerce')
    
    # Variants missing from gnomAD (AF is NA)
    missing = df[df[af_col].isna()].copy()
    n_missing = len(missing)
    pct_missing = (n_missing / total * 100) if total > 0 else 0
    
    print(f"\n🔍 Missing from gnomAD: {n_missing:,} ({pct_missing:.1f}%)")
    
    if n_missing == 0:
        print("✅ All variants have gnomAD annotation!")
        return
    
    # Analyze missing variants by size
    print(f"\n📏 Analyzing variant sizes for missing variants...")
    
    # Calculate reference allele length
    if 'Ref' in missing.columns:
        ref_col = 'Ref'
    elif 'ref_allele' in missing.columns:
        ref_col = 'ref_allele'
    else:
        print(f"❌ ERROR: Cannot find Ref allele column")
        return
    
    missing['ref_len'] = missing[ref_col].astype(str).str.len()
    
    # Classify by size
    small_variants = missing[missing['ref_len'] <= LARGE_DELETION_THRESHOLD]
    large_deletions = missing[missing['ref_len'] > LARGE_DELETION_THRESHOLD]
    
    n_small = len(small_variants)
    n_large = len(large_deletions)
    
    print(f"\n   Small variants (≤{LARGE_DELETION_THRESHOLD}bp): {n_small:,}")
    print(f"   Large deletions (>{LARGE_DELETION_THRESHOLD}bp): {n_large:,}")
    
    # Structural variant gap
    sv_gap = (n_large / total * 100) if total > 0 else 0
    print(f"\n📊 Structural variant gap: {n_large:,}/{total:,} = {sv_gap:.2f}%")
    
    # Size distribution of missing large deletions
    if n_large > 0:
        print(f"\n📏 Large deletion size distribution:")
        print(f"   Min: {large_deletions['ref_len'].min()} bp")
        print(f"   Median: {large_deletions['ref_len'].median():.0f} bp")
        print(f"   Mean: {large_deletions['ref_len'].mean():.0f} bp")
        print(f"   Max: {large_deletions['ref_len'].max()} bp")
        
        # Show size bins
        bins = [50, 100, 500, 1000, 5000, 10000, float('inf')]
        labels = ['51-100bp', '101-500bp', '501-1000bp', '1001-5000bp', '5001-10000bp', '>10000bp']
        large_deletions['size_bin'] = pd.cut(large_deletions['ref_len'], bins=bins, labels=labels)
        
        print(f"\n   Size distribution:")
        for label in labels:
            count = (large_deletions['size_bin'] == label).sum()
            if count > 0:
                print(f"      {label}: {count:,}")
    
    # Gene distribution of missing variants
    if 'Gene' in missing.columns:
        gene_col = 'Gene'
    elif 'gene_symbol' in missing.columns:
        gene_col = 'gene_symbol'
    else:
        gene_col = None
    
    if gene_col:
        print(f"\n🧬 Top 10 genes with missing variants:")
        gene_counts = missing[gene_col].value_counts().head(10)
        for gene, count in gene_counts.items():
            print(f"   {gene}: {count:,}")
    
    # Star rating distribution (for P/LP)
    if 'StarRating_Numeric' in missing.columns or 'StarRating' in missing.columns:
        star_col = 'StarRating_Numeric' if 'StarRating_Numeric' in missing.columns else 'StarRating'
        print(f"\n⭐ Star rating distribution of missing variants:")
        star_counts = missing[star_col].value_counts().sort_index()
        for stars, count in star_counts.items():
            print(f"   {stars} stars: {count:,}")
    
    return {
        'variant_type': variant_type,
        'total': total,
        'missing': n_missing,
        'missing_pct': pct_missing,
        'small_missing': n_small,
        'large_deletions': n_large,
        'sv_gap_pct': sv_gap
    }

# ============================================================================
# LOAD AND ANALYZE P/LP VARIANTS
# ============================================================================

print(f"\n📂 Loading P/LP file (gnomAD JOINT)...")
try:
    plp_df = pd.read_csv('any_PLP_MASTER_JOINT.txt', sep='\t')
    print(f"✅ Loaded: {len(plp_df):,} variants")
    
    plp_results = analyze_sv_coverage(plp_df, 'P/LP')
except FileNotFoundError:
    print(f"❌ ERROR: any_PLP_MASTER_JOINT.txt not found")
    plp_results = None

# ============================================================================
# LOAD AND ANALYZE VUS VARIANTS
# ============================================================================

print(f"\n📂 Loading VUS file (gnomAD JOINT)...")
try:
    vus_df = pd.read_csv('VUS_2plus_stars_JOINT.txt', sep='\t')
    print(f"✅ Loaded: {len(vus_df):,} variants")
    
    vus_results = analyze_sv_coverage(vus_df, 'VUS')
except FileNotFoundError:
    print(f"❌ ERROR: VUS_2plus_stars_JOINT.txt not found")
    vus_results = None

# ============================================================================
# COMBINED SUMMARY
# ============================================================================

print(f"\n{'=' * 80}")
print("COMBINED SUMMARY")
print(f"{'=' * 80}")

if plp_results and vus_results:
    total_vars = plp_results['total'] + vus_results['total']
    total_missing = plp_results['missing'] + vus_results['missing']
    total_large_dels = plp_results['large_deletions'] + vus_results['large_deletions']
    
    print(f"\n📊 Overall Statistics (P/LP + VUS):")
    print(f"   Total variants: {total_vars:,}")
    print(f"   Missing from gnomAD: {total_missing:,} ({total_missing/total_vars*100:.1f}%)")
    print(f"   Large deletions (>{LARGE_DELETION_THRESHOLD}bp) missing: {total_large_dels:,}")
    print(f"   Structural variant gap: {total_large_dels/total_vars*100:.2f}%")
    
    print(f"\n📋 Breakdown by variant type:")
    print(f"\n   P/LP:")
    print(f"      Total: {plp_results['total']:,}")
    print(f"      Missing: {plp_results['missing']:,} ({plp_results['missing_pct']:.1f}%)")
    print(f"      Large deletions missing: {plp_results['large_deletions']:,} ({plp_results['sv_gap_pct']:.2f}%)")
    
    print(f"\n   VUS:")
    print(f"      Total: {vus_results['total']:,}")
    print(f"      Missing: {vus_results['missing']:,} ({vus_results['missing_pct']:.1f}%)")
    print(f"      Large deletions missing: {vus_results['large_deletions']:,} ({vus_results['sv_gap_pct']:.2f}%)")

print(f"\n{'=' * 80}")
print("✅ ANALYSIS COMPLETE")
print(f"{'=' * 80}")

print(f"""
Key Findings:
  • gnomAD v4.1 JOINT dataset (807,162 individuals)
  • Large deletions = variants with Ref allele >{LARGE_DELETION_THRESHOLD}bp
  • Missing variants are typically ultra-rare (not in gnomAD cohort)
  • SV gap represents limitation of short-read sequencing for structural variants
""")
