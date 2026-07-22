"""
Identify genes with carrier frequency >= 1/200 (0.5%)
Analyzes All of Us and gnomAD gene-level data
"""

import pandas as pd
import numpy as np

# Threshold: 1 in 200 = 0.005 = 0.5%
CARRIER_FREQ_THRESHOLD = 0.005

print("=" * 80)
print("HIGH CARRIER FREQUENCY GENES ANALYSIS")
print("=" * 80)
print(f"\nThreshold: ≥ 1/200 (carrier frequency ≥ {CARRIER_FREQ_THRESHOLD} = {CARRIER_FREQ_THRESHOLD*100}%)")
print("\nAnalyzing P/LP variants (≥2 stars)")

# ============================================================================
# Load gene-level data
# ============================================================================

print("\n" + "=" * 80)
print("LOADING DATA")
print("=" * 80)

# All of Us
aou_plp = pd.read_csv('AoU_PLP_GENE_LEVEL.csv')
print(f"\n✅ All of Us P/LP: {len(aou_plp):,} rows, {aou_plp['gene_symbol'].nunique()} unique genes")

# gnomAD
gnomad_plp = pd.read_csv('gnomAD_PLP_GENE_LEVEL.csv')
print(f"✅ gnomAD P/LP: {len(gnomad_plp):,} rows, {gnomad_plp['gene_symbol'].nunique()} unique genes")

# ============================================================================
# All of Us: Find high carrier frequency genes
# ============================================================================

print("\n" + "=" * 80)
print("ALL OF US - HIGH CARRIER FREQUENCY GENES")
print("=" * 80)

# Filter for combined ancestry (ALL)
aou_all = aou_plp[aou_plp['ancestry'] == 'ALL'].copy()
aou_high = aou_all[aou_all['carrier_frequency'] >= CARRIER_FREQ_THRESHOLD].copy()
aou_high = aou_high.sort_values('carrier_frequency', ascending=False)

print(f"\n📊 Genes with CF ≥ 1/200 in ALL ancestry: {len(aou_high)}")

if len(aou_high) > 0:
    print("\nGene Rankings (ALL ancestry):")
    print("-" * 80)
    for i, row in aou_high.iterrows():
        cf_pct = row['carrier_frequency'] * 100
        nns = 1 / row['carrier_frequency'] if row['carrier_frequency'] > 0 else np.inf
        print(f"  {row['gene_symbol']:15} CF = {cf_pct:6.2f}%  (1 in {nns:5.0f})")

# Check by ancestry
print("\n" + "-" * 80)
print("BY ANCESTRY:")
print("-" * 80)

ancestry_summary = []
for ancestry in ['ALL', 'AFR', 'AMR', 'EAS', 'EUR']:
    aou_anc = aou_plp[aou_plp['ancestry'] == ancestry].copy()
    aou_anc_high = aou_anc[aou_anc['carrier_frequency'] >= CARRIER_FREQ_THRESHOLD].copy()
    
    genes_high = aou_anc_high['gene_symbol'].unique()
    
    print(f"\n{ancestry:12} {len(aou_anc_high):3} genes with CF ≥ 1/200")
    
    if len(aou_anc_high) > 0:
        aou_anc_high = aou_anc_high.sort_values('carrier_frequency', ascending=False)
        top_5 = aou_anc_high.head(5)
        for _, row in top_5.iterrows():
            cf_pct = row['carrier_frequency'] * 100
            print(f"             {row['gene_symbol']:15} {cf_pct:6.2f}%")
    
    ancestry_summary.append({
        'ancestry': ancestry,
        'n_genes_high_cf': len(aou_anc_high),
        'genes': list(genes_high)
    })

# ============================================================================
# gnomAD: Find high carrier frequency genes
# ============================================================================

print("\n" + "=" * 80)
print("gnomAD - HIGH CARRIER FREQUENCY GENES")
print("=" * 80)

# Filter for total ancestry (TOTAL)
gnomad_total = gnomad_plp[gnomad_plp['ancestry'] == 'TOTAL'].copy()
gnomad_high = gnomad_total[gnomad_total['carrier_frequency'] >= CARRIER_FREQ_THRESHOLD].copy()
gnomad_high = gnomad_high.sort_values('carrier_frequency', ascending=False)

print(f"\n📊 Genes with CF ≥ 1/200 in TOTAL ancestry: {len(gnomad_high)}")

if len(gnomad_high) > 0:
    print("\nGene Rankings (TOTAL ancestry):")
    print("-" * 80)
    for i, row in gnomad_high.iterrows():
        cf_pct = row['carrier_frequency'] * 100
        nns = 1 / row['carrier_frequency'] if row['carrier_frequency'] > 0 else np.inf
        print(f"  {row['gene_symbol']:15} CF = {cf_pct:6.2f}%  (1 in {nns:5.0f})")

# Check by ancestry
print("\n" + "-" * 80)
print("BY ANCESTRY:")
print("-" * 80)

gnomad_ancestry_summary = []
for ancestry in ['TOTAL', 'AFR', 'AMR', 'EAS', 'NFE', 'ASJ']:
    gnomad_anc = gnomad_plp[gnomad_plp['ancestry'] == ancestry].copy()
    gnomad_anc_high = gnomad_anc[gnomad_anc['carrier_frequency'] >= CARRIER_FREQ_THRESHOLD].copy()
    
    genes_high = gnomad_anc_high['gene_symbol'].unique()
    
    print(f"\n{ancestry:12} {len(gnomad_anc_high):3} genes with CF ≥ 1/200")
    
    if len(gnomad_anc_high) > 0:
        gnomad_anc_high = gnomad_anc_high.sort_values('carrier_frequency', ascending=False)
        top_5 = gnomad_anc_high.head(5)
        for _, row in top_5.iterrows():
            cf_pct = row['carrier_frequency'] * 100
            print(f"             {row['gene_symbol']:15} {cf_pct:6.2f}%")
    
    gnomad_ancestry_summary.append({
        'ancestry': ancestry,
        'n_genes_high_cf': len(gnomad_anc_high),
        'genes': list(genes_high)
    })

# ============================================================================
# Compare: AoU vs gnomAD
# ============================================================================

print("\n" + "=" * 80)
print("COMPARISON: ALL OF US vs gnomAD")
print("=" * 80)

aou_genes_high = set(aou_high['gene_symbol'].unique())
gnomad_genes_high = set(gnomad_high['gene_symbol'].unique())

genes_both = aou_genes_high & gnomad_genes_high
genes_aou_only = aou_genes_high - gnomad_genes_high
genes_gnomad_only = gnomad_genes_high - aou_genes_high

print(f"\n📊 SUMMARY (Combined/Total ancestry):")
print(f"   Genes ≥ 1/200 in BOTH datasets:     {len(genes_both):3}")
print(f"   Genes ≥ 1/200 in All of Us ONLY:    {len(genes_aou_only):3}")
print(f"   Genes ≥ 1/200 in gnomAD ONLY:       {len(genes_gnomad_only):3}")
print(f"   Total unique genes ≥ 1/200:         {len(aou_genes_high | gnomad_genes_high):3}")

if len(genes_both) > 0:
    print(f"\nGenes in BOTH datasets (n={len(genes_both)}):")
    # Get carrier frequencies for comparison
    both_comparison = []
    for gene in sorted(genes_both):
        aou_cf = aou_high[aou_high['gene_symbol'] == gene]['carrier_frequency'].iloc[0]
        gnomad_cf = gnomad_high[gnomad_high['gene_symbol'] == gene]['carrier_frequency'].iloc[0]
        both_comparison.append({
            'gene': gene,
            'aou_cf': aou_cf,
            'gnomad_cf': gnomad_cf,
            'avg_cf': (aou_cf + gnomad_cf) / 2
        })
    
    both_df = pd.DataFrame(both_comparison).sort_values('avg_cf', ascending=False)
    print("  " + "-" * 76)
    print(f"  {'Gene':<15} {'AoU CF':>10} {'gnomAD CF':>12} {'Avg NNS':>12}")
    print("  " + "-" * 76)
    for _, row in both_df.iterrows():
        aou_pct = row['aou_cf'] * 100
        gnomad_pct = row['gnomad_cf'] * 100
        avg_nns = 1 / row['avg_cf']
        print(f"  {row['gene']:<15} {aou_pct:9.2f}% {gnomad_pct:11.2f}%  {avg_nns:11.0f}")

if len(genes_aou_only) > 0:
    print(f"\nGenes in All of Us ONLY (n={len(genes_aou_only)}):")
    aou_only_sorted = aou_high[aou_high['gene_symbol'].isin(genes_aou_only)].sort_values('carrier_frequency', ascending=False)
    for _, row in aou_only_sorted.iterrows():
        cf_pct = row['carrier_frequency'] * 100
        print(f"  {row['gene_symbol']:<15} {cf_pct:6.2f}%")

if len(genes_gnomad_only) > 0:
    print(f"\nGenes in gnomAD ONLY (n={len(genes_gnomad_only)}):")
    gnomad_only_sorted = gnomad_high[gnomad_high['gene_symbol'].isin(genes_gnomad_only)].sort_values('carrier_frequency', ascending=False)
    for _, row in gnomad_only_sorted.iterrows():
        cf_pct = row['carrier_frequency'] * 100
        print(f"  {row['gene_symbol']:<15} {cf_pct:6.2f}%")

# ============================================================================
# Save detailed results
# ============================================================================

print("\n" + "=" * 80)
print("SAVING RESULTS")
print("=" * 80)

# All of Us high CF genes (all ancestries)
aou_all_high = aou_plp[aou_plp['carrier_frequency'] >= CARRIER_FREQ_THRESHOLD].copy()
aou_all_high['1_in_N'] = 1 / aou_all_high['carrier_frequency']
aou_all_high = aou_all_high.sort_values(['ancestry', 'carrier_frequency'], ascending=[True, False])
aou_all_high.to_csv('AoU_High_CF_Genes_1in200.csv', index=False)
print(f"\n✅ Saved: AoU_High_CF_Genes_1in200.csv ({len(aou_all_high)} rows)")

# gnomAD high CF genes (all ancestries)
gnomad_all_high = gnomad_plp[gnomad_plp['carrier_frequency'] >= CARRIER_FREQ_THRESHOLD].copy()
gnomad_all_high['1_in_N'] = 1 / gnomad_all_high['carrier_frequency']
gnomad_all_high = gnomad_all_high.sort_values(['ancestry', 'carrier_frequency'], ascending=[True, False])
gnomad_all_high.to_csv('gnomAD_High_CF_Genes_1in200.csv', index=False)
print(f"✅ Saved: gnomAD_High_CF_Genes_1in200.csv ({len(gnomad_all_high)} rows)")

# Comparison summary
comparison_summary = []
all_genes = aou_genes_high | gnomad_genes_high

for gene in sorted(all_genes):
    aou_row = aou_high[aou_high['gene_symbol'] == gene]
    gnomad_row = gnomad_high[gnomad_high['gene_symbol'] == gene]
    
    aou_cf = aou_row['carrier_frequency'].iloc[0] if len(aou_row) > 0 else 0
    gnomad_cf = gnomad_row['carrier_frequency'].iloc[0] if len(gnomad_row) > 0 else 0
    
    in_aou = len(aou_row) > 0
    in_gnomad = len(gnomad_row) > 0
    
    comparison_summary.append({
        'gene_symbol': gene,
        'in_AoU': in_aou,
        'in_gnomAD': in_gnomad,
        'in_both': in_aou and in_gnomad,
        'AoU_carrier_freq': aou_cf,
        'gnomAD_carrier_freq': gnomad_cf,
        'AoU_1_in_N': 1/aou_cf if aou_cf > 0 else np.inf,
        'gnomAD_1_in_N': 1/gnomad_cf if gnomad_cf > 0 else np.inf,
        'avg_carrier_freq': (aou_cf + gnomad_cf) / 2 if (in_aou and in_gnomad) else (aou_cf if in_aou else gnomad_cf),
    })

comparison_df = pd.DataFrame(comparison_summary)
comparison_df = comparison_df.sort_values('avg_carrier_freq', ascending=False)
comparison_df.to_csv('High_CF_Genes_Comparison_1in200.csv', index=False)
print(f"✅ Saved: High_CF_Genes_Comparison_1in200.csv ({len(comparison_df)} unique genes)")

# ============================================================================
# Summary
# ============================================================================

print("\n" + "=" * 80)
print("✅ ANALYSIS COMPLETE")
print("=" * 80)

print(f"""
OUTPUT FILES:
  ✅ AoU_High_CF_Genes_1in200.csv
     - All genes with CF ≥ 1/200 in All of Us
     - All ancestries included
     - {len(aou_all_high)} rows
  
  ✅ gnomAD_High_CF_Genes_1in200.csv
     - All genes with CF ≥ 1/200 in gnomAD
     - All ancestries included
     - {len(gnomad_all_high)} rows
  
  ✅ High_CF_Genes_Comparison_1in200.csv
     - Comparison of AoU vs gnomAD
     - Combined/Total ancestry only
     - {len(comparison_df)} unique genes

KEY FINDINGS (Combined/Total ancestry):
  • {len(genes_both)} genes ≥ 1/200 in BOTH datasets
  • {len(genes_aou_only)} genes ≥ 1/200 in All of Us ONLY
  • {len(genes_gnomad_only)} genes ≥ 1/200 in gnomAD ONLY
  • {len(aou_genes_high | gnomad_genes_high)} total unique genes ≥ 1/200

Threshold: 1 in 200 = 0.5% carrier frequency
Data: P/LP variants with ≥2 star rating
""")
