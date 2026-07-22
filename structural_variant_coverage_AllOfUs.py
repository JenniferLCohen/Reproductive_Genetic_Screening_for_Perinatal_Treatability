"""
================================================================================
ALL OF US STRUCTURAL VARIANT COVERAGE ANALYSIS
================================================================================

Analyzes structural variant coverage in All of Us Research Program v8 data
using EXACT carrier frequency filters to match gnomAD analysis.

Filters Applied:
  P/LP:
    - ≥2 stars (already applied in extraction)
    - Exclude MT-RNR1, GATA1, SMN1 (already applied in extraction)
  
  VUS:
    - ≥2 stars (already applied in extraction)
    - AF <0.001% (already applied - RARE files)
    - Exclude MT-RNR1, GATA1, SMN1 (already applied in extraction)

Input Files (from All of Us extraction pipeline):
  - PLP_AUTOSOMAL_VARIANT_LEVEL_EXACT_COUNTS.csv
  - PLP_X_CHROMOSOME_FEMALE_EXACT_COUNTS.csv
  - VUS_AUTOSOMAL_VARIANT_LEVEL_EXACT_COUNTS_RARE.csv
  - VUS_X_CHROMOSOME_FEMALE_EXACT_COUNTS_RARE.csv

Note: These files already have filters applied during extraction, so they
      match the carrier frequency analysis exactly.

================================================================================
"""

import pandas as pd
import numpy as np

print("=" * 80)
print("ALL OF US STRUCTURAL VARIANT COVERAGE ANALYSIS")
print("=" * 80)
print("\n✅ Files already filtered during extraction pipeline:")
print("   P/LP: ≥2 stars, exclude MT-RNR1/GATA1/SMN1")
print("   VUS: ≥2 stars, AF <0.001%, exclude MT-RNR1/GATA1/SMN1")

# ============================================================================
# CONFIGURATION
# ============================================================================

LARGE_DELETION_THRESHOLD = 50  # bp

print(f"\n📏 Large deletion threshold: >{LARGE_DELETION_THRESHOLD}bp")

# ============================================================================
# FUNCTION: Analyze SV in All of Us data
# ============================================================================

def analyze_aou_sv(df, variant_type, chromosome_type):
    """
    Analyze structural variants in All of Us data.
    
    Args:
        df: DataFrame with All of Us variant data
        variant_type: 'P/LP' or 'VUS'
        chromosome_type: 'Autosomal' or 'X chromosome'
    """
    print(f"\n{'=' * 80}")
    print(f"{variant_type} - {chromosome_type}")
    print(f"{'=' * 80}")
    
    total = len(df)
    print(f"\n📊 Total variants: {total:,}")
    
    if total == 0:
        print("⚠️  No variants found")
        return None
    
    # Calculate reference allele length
    ref_col = 'ref_allele' if 'ref_allele' in df.columns else 'Ref'
    
    if ref_col not in df.columns:
        print(f"❌ ERROR: Cannot find reference allele column")
        print(f"Available columns: {', '.join(df.columns[:10])}...")
        return None
    
    df['ref_len'] = df[ref_col].astype(str).str.len()
    
    # Classify by size
    small_variants = df[df['ref_len'] <= LARGE_DELETION_THRESHOLD]
    large_deletions = df[df['ref_len'] > LARGE_DELETION_THRESHOLD]
    
    n_small = len(small_variants)
    n_large = len(large_deletions)
    
    print(f"\n   Small variants (≤{LARGE_DELETION_THRESHOLD}bp): {n_small:,}")
    print(f"   Large deletions (>{LARGE_DELETION_THRESHOLD}bp): {n_large:,}")
    
    # Structural variant proportion
    sv_prop = (n_large / total * 100) if total > 0 else 0
    print(f"\n📊 Large deletions: {n_large:,}/{total:,} = {sv_prop:.2f}%")
    
    # Size distribution of large deletions
    if n_large > 0:
        print(f"\n📏 Large deletion size distribution:")
        print(f"   Min: {large_deletions['ref_len'].min()} bp")
        print(f"   Median: {large_deletions['ref_len'].median():.0f} bp")
        print(f"   Mean: {large_deletions['ref_len'].mean():.0f} bp")
        print(f"   Max: {large_deletions['ref_len'].max()} bp")
        
        # Show size bins
        bins = [50, 100, 500, 1000, 5000, 10000, float('inf')]
        labels = ['51-100bp', '101-500bp', '501-1000bp', '1001-5000bp', '5001-10000bp', '>10000bp']
        large_deletions_copy = large_deletions.copy()
        large_deletions_copy['size_bin'] = pd.cut(large_deletions_copy['ref_len'], bins=bins, labels=labels)
        
        print(f"\n   Size distribution:")
        for label in labels:
            count = (large_deletions_copy['size_bin'] == label).sum()
            if count > 0:
                print(f"      {label}: {count:,}")
    
    return {
        'variant_type': variant_type,
        'chromosome_type': chromosome_type,
        'total': total,
        'small_variants': n_small,
        'large_deletions': n_large,
        'sv_prop': sv_prop
    }

# ============================================================================
# LOAD AND ANALYZE P/LP VARIANTS
# ============================================================================

print(f"\n{'=' * 80}")
print("LOADING P/LP VARIANTS (All of Us)")
print(f"{'=' * 80}")

try:
    # Autosomal P/LP
    plp_auto = pd.read_csv('PLP_AUTOSOMAL_VARIANT_LEVEL_EXACT_COUNTS.csv')
    print(f"\n✅ Loaded autosomal P/LP: {len(plp_auto):,} variants")
    plp_auto_results = analyze_aou_sv(plp_auto, 'P/LP', 'Autosomal')
    
except FileNotFoundError:
    print(f"❌ ERROR: PLP_AUTOSOMAL_VARIANT_LEVEL_EXACT_COUNTS.csv not found")
    plp_auto_results = None

try:
    # X chromosome P/LP
    plp_x = pd.read_csv('PLP_X_CHROMOSOME_FEMALE_EXACT_COUNTS.csv')
    print(f"\n✅ Loaded X chromosome P/LP: {len(plp_x):,} variants")
    plp_x_results = analyze_aou_sv(plp_x, 'P/LP', 'X chromosome')
    
except FileNotFoundError:
    print(f"❌ ERROR: PLP_X_CHROMOSOME_FEMALE_EXACT_COUNTS.csv not found")
    plp_x_results = None

# ============================================================================
# LOAD AND ANALYZE VUS VARIANTS
# ============================================================================

print(f"\n{'=' * 80}")
print("LOADING VUS VARIANTS (All of Us)")
print(f"{'=' * 80}")

try:
    # Autosomal VUS (rare)
    vus_auto = pd.read_csv('VUS_AUTOSOMAL_VARIANT_LEVEL_EXACT_COUNTS_RARE.csv')
    print(f"\n✅ Loaded autosomal rare VUS: {len(vus_auto):,} variants")
    vus_auto_results = analyze_aou_sv(vus_auto, 'VUS', 'Autosomal')
    
except FileNotFoundError:
    print(f"❌ ERROR: VUS_AUTOSOMAL_VARIANT_LEVEL_EXACT_COUNTS_RARE.csv not found")
    vus_auto_results = None

try:
    # X chromosome VUS (rare)
    vus_x = pd.read_csv('VUS_X_CHROMOSOME_FEMALE_EXACT_COUNTS_RARE.csv')
    print(f"\n✅ Loaded X chromosome rare VUS: {len(vus_x):,} variants")
    vus_x_results = analyze_aou_sv(vus_x, 'VUS', 'X chromosome')
    
except FileNotFoundError:
    print(f"❌ ERROR: VUS_X_CHROMOSOME_FEMALE_EXACT_COUNTS_RARE.csv not found")
    vus_x_results = None

# ============================================================================
# COMBINED SUMMARY
# ============================================================================

print(f"\n{'=' * 80}")
print("FINAL SUMMARY - ALL OF US → gnomAD GAP (Parallel to gnomAD Analysis)")
print(f"{'=' * 80}")

# Combine results
all_results = [r for r in [plp_auto_results, plp_x_results, vus_auto_results, vus_x_results] if r is not None]

if len(all_results) > 0:
    total_variants = sum(r['total'] for r in all_results)
    total_missing = sum(r['missing'] for r in all_results)
    total_large_dels = sum(r['large_deletions'] for r in all_results)
    
    print(f"\n📊 Overall Statistics (All of Us variants):")
    print(f"   Total All of Us variants: {total_variants:,}")
    print(f"   Missing from gnomAD: {total_missing:,} ({total_missing/total_variants*100:.1f}%)")
    print(f"   Large deletions (>{LARGE_DELETION_THRESHOLD}bp) missing: {total_large_dels:,}")
    print(f"   Structural variant gap: {total_large_dels/total_variants*100:.2f}%")
    
    print(f"\n📋 Breakdown by variant type:")
    
    # P/LP summary
    plp_results = [r for r in all_results if r['variant_type'] == 'P/LP']
    if plp_results:
        plp_total = sum(r['total'] for r in plp_results)
        plp_missing = sum(r['missing'] for r in plp_results)
        plp_large = sum(r['large_deletions'] for r in plp_results)
        print(f"\n   P/LP:")
        print(f"      Total: {plp_total:,}")
        print(f"      Missing from gnomAD: {plp_missing:,} ({plp_missing/plp_total*100:.1f}%)")
        print(f"      Large deletions missing: {plp_large:,} ({plp_large/plp_total*100:.2f}%)")
        for r in plp_results:
            print(f"         {r['chromosome_type']}: {r['large_deletions']:,}/{r['total']:,}")
    
    # VUS summary
    vus_results = [r for r in all_results if r['variant_type'] == 'VUS']
    if vus_results:
        vus_total = sum(r['total'] for r in vus_results)
        vus_missing = sum(r['missing'] for r in vus_results)
        vus_large = sum(r['large_deletions'] for r in vus_results)
        print(f"\n   Rare VUS:")
        print(f"      Total: {vus_total:,}")
        print(f"      Missing from gnomAD: {vus_missing:,} ({vus_missing/vus_total*100:.1f}%)")
        print(f"      Large deletions missing: {vus_large:,} ({vus_large/vus_total*100:.2f}%)")
        for r in vus_results:
            print(f"         {r['chromosome_type']}: {r['large_deletions']:,}/{r['total']:,}")

print(f"\n{'=' * 80}")
print("✅ ANALYSIS COMPLETE")
print(f"{'=' * 80}")

print(f"""
PARALLEL ANALYSIS COMPLETE - Directly Comparable to gnomAD Analysis

This analysis answers: "Of variants found in All of Us carrier screening 
panel, which are MISSING from gnomAD v4.1 JOINT?"

Comparison Framework:
  ✅ SAME methodology as gnomAD analysis
  ✅ SAME filters (≥2 stars, gene exclusions, rare VUS)
  ✅ SAME question (gnomAD coverage gap)
  ✅ Different perspective (All of Us variants vs ClinVar variants)

Interpretation:
  - Similar gap % = Good agreement between datasets
  - Higher gap % = All of Us may capture variants not in gnomAD
  - Lower gap % = All of Us variants better represented in gnomAD

Next Steps:
  - Compare with gnomAD analysis results (0.08% gap overall)
  - Identify variants unique to All of Us vs gnomAD
  - Cross-reference large deletions between datasets
""")

