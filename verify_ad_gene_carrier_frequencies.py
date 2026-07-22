"""
Verify: Do all AD genes have carrier frequency < 0.5%?
Checks All of Us and gnomAD data for AD gene carrier frequencies
"""

import pandas as pd
import numpy as np

print("=" * 80)
print("VERIFICATION: AD GENE CARRIER FREQUENCIES")
print("=" * 80)

THRESHOLD = 0.005  # 0.5% = 1/200

# Load gene classifications
gene_class = pd.read_csv('GENE_CLASSIFICATIONS_293.csv')
ad_genes = gene_class[gene_class['inheritance'] == 'AD']['gene_symbol'].tolist()

print(f"\n📊 AD genes in panel: {len(ad_genes)}")
print(f"   Threshold: ≥{THRESHOLD*100}% (1 in {1/THRESHOLD:.0f})")

# Load gene-level data
aou_plp = pd.read_csv('AoU_PLP_GENE_LEVEL.csv')
gnomad_plp = pd.read_csv('gnomAD_PLP_GENE_LEVEL.csv')

# Filter to AD genes, combined/total ancestry
aou_ad = aou_plp[
    (aou_plp['gene_symbol'].isin(ad_genes)) & 
    (aou_plp['ancestry'] == 'ALL')
].copy()

gnomad_ad = gnomad_plp[
    (gnomad_plp['gene_symbol'].isin(ad_genes)) & 
    (gnomad_plp['ancestry'] == 'TOTAL')
].copy()

print(f"\n📂 Data loaded:")
print(f"   All of Us: {len(aou_ad)} AD genes with data")
print(f"   gnomAD: {len(gnomad_ad)} AD genes with data")

# Check for genes ≥ threshold
print(f"\n{'=' * 80}")
print(f"CHECKING FOR AD GENES WITH CF ≥ {THRESHOLD*100}%")
print(f"{'=' * 80}")

# All of Us
aou_high = aou_ad[aou_ad['carrier_frequency'] >= THRESHOLD].copy()
aou_high = aou_high.sort_values('carrier_frequency', ascending=False)

print(f"\n📊 All of Us - AD genes with CF ≥ {THRESHOLD*100}%:")
if len(aou_high) > 0:
    print(f"   ❌ FOUND {len(aou_high)} AD genes ≥ {THRESHOLD*100}%")
    print(f"\n   {'Gene':<15} {'CF':<12} {'1 in N':<12}")
    print(f"   {'-'*40}")
    for _, row in aou_high.iterrows():
        cf_pct = row['carrier_frequency'] * 100
        nns = 1 / row['carrier_frequency'] if row['carrier_frequency'] > 0 else np.inf
        print(f"   {row['gene_symbol']:<15} {cf_pct:>6.3f}%      1 in {nns:>5.0f}")
else:
    print(f"   ✅ NONE - All AD genes have CF < {THRESHOLD*100}%")

# gnomAD
gnomad_high = gnomad_ad[gnomad_ad['carrier_frequency'] >= THRESHOLD].copy()
gnomad_high = gnomad_high.sort_values('carrier_frequency', ascending=False)

print(f"\n📊 gnomAD - AD genes with CF ≥ {THRESHOLD*100}%:")
if len(gnomad_high) > 0:
    print(f"   ❌ FOUND {len(gnomad_high)} AD genes ≥ {THRESHOLD*100}%")
    print(f"\n   {'Gene':<15} {'CF':<12} {'1 in N':<12}")
    print(f"   {'-'*40}")
    for _, row in gnomad_high.iterrows():
        cf_pct = row['carrier_frequency'] * 100
        nns = 1 / row['carrier_frequency'] if row['carrier_frequency'] > 0 else np.inf
        print(f"   {row['gene_symbol']:<15} {cf_pct:>6.3f}%      1 in {nns:>5.0f}")
else:
    print(f"   ✅ NONE - All AD genes have CF < {THRESHOLD*100}%")

# Summary statistics
print(f"\n{'=' * 80}")
print(f"SUMMARY STATISTICS - AD GENES")
print(f"{'=' * 80}")

print(f"\nAll of Us (n={len(aou_ad)}):")
if len(aou_ad) > 0:
    print(f"   Maximum CF: {aou_ad['carrier_frequency'].max()*100:.3f}%")
    print(f"   Mean CF: {aou_ad['carrier_frequency'].mean()*100:.3f}%")
    print(f"   Median CF: {aou_ad['carrier_frequency'].median()*100:.3f}%")
    
    # Top 5
    top5_aou = aou_ad.nlargest(5, 'carrier_frequency')
    print(f"\n   Top 5 AD genes by carrier frequency:")
    for _, row in top5_aou.iterrows():
        cf_pct = row['carrier_frequency'] * 100
        print(f"      {row['gene_symbol']:<15} {cf_pct:>6.3f}%")

print(f"\ngnomAD (n={len(gnomad_ad)}):")
if len(gnomad_ad) > 0:
    print(f"   Maximum CF: {gnomad_ad['carrier_frequency'].max()*100:.3f}%")
    print(f"   Mean CF: {gnomad_ad['carrier_frequency'].mean()*100:.3f}%")
    print(f"   Median CF: {gnomad_ad['carrier_frequency'].median()*100:.3f}%")
    
    # Top 5
    top5_gnomad = gnomad_ad.nlargest(5, 'carrier_frequency')
    print(f"\n   Top 5 AD genes by carrier frequency:")
    for _, row in top5_gnomad.iterrows():
        cf_pct = row['carrier_frequency'] * 100
        print(f"      {row['gene_symbol']:<15} {cf_pct:>6.3f}%")

# Final verdict
print(f"\n{'=' * 80}")
print(f"FINAL VERDICT")
print(f"{'=' * 80}")

any_high_aou = len(aou_high) > 0
any_high_gnomad = len(gnomad_high) > 0

if not any_high_aou and not any_high_gnomad:
    print(f"\n✅ STATEMENT IS TRUE:")
    print(f"   All autosomal dominant genes on this list have carrier")
    print(f"   frequencies < {THRESHOLD*100}% in both All of Us and gnomAD databases.")
    print(f"\n   Maximum observed:")
    if len(aou_ad) > 0 and len(gnomad_ad) > 0:
        max_aou = aou_ad['carrier_frequency'].max() * 100
        max_gnomad = gnomad_ad['carrier_frequency'].max() * 100
        print(f"     All of Us: {max_aou:.3f}%")
        print(f"     gnomAD: {max_gnomad:.3f}%")
elif any_high_aou or any_high_gnomad:
    print(f"\n❌ STATEMENT IS FALSE:")
    print(f"   Some AD genes have carrier frequencies ≥ {THRESHOLD*100}%")
    if any_high_aou:
        print(f"\n   All of Us: {len(aou_high)} genes ≥ {THRESHOLD*100}%")
    if any_high_gnomad:
        print(f"   gnomAD: {len(gnomad_high)} genes ≥ {THRESHOLD*100}%")
    
    print(f"\n   SUGGESTED REVISED STATEMENT:")
    if len(aou_ad) > 0 and len(gnomad_ad) > 0:
        max_aou = aou_ad['carrier_frequency'].max() * 100
        max_gnomad = gnomad_ad['carrier_frequency'].max() * 100
        max_overall = max(max_aou, max_gnomad)
        
        # Suggest appropriate threshold
        if max_overall < 0.1:
            suggested = 0.1
        elif max_overall < 0.2:
            suggested = 0.2
        elif max_overall < 0.3:
            suggested = 0.3
        else:
            suggested = round(max_overall + 0.05, 1)
        
        print(f"   'Autosomal dominant genes on this list have carrier")
        print(f"    frequencies < {suggested}% in All of Us and gnomAD databases'")

# Save detailed results
results = []

for _, row in aou_ad.iterrows():
    gnomad_row = gnomad_ad[gnomad_ad['gene_symbol'] == row['gene_symbol']]
    gnomad_cf = gnomad_row['carrier_frequency'].iloc[0] if len(gnomad_row) > 0 else np.nan
    
    results.append({
        'gene_symbol': row['gene_symbol'],
        'AoU_carrier_freq': row['carrier_frequency'],
        'AoU_CF_pct': row['carrier_frequency'] * 100,
        'gnomAD_carrier_freq': gnomad_cf,
        'gnomAD_CF_pct': gnomad_cf * 100 if not pd.isna(gnomad_cf) else np.nan,
        'max_CF_pct': max(row['carrier_frequency'] * 100, 
                         gnomad_cf * 100 if not pd.isna(gnomad_cf) else 0),
        'exceeds_0.5pct': (row['carrier_frequency'] >= THRESHOLD) or 
                          (gnomad_cf >= THRESHOLD if not pd.isna(gnomad_cf) else False)
    })

results_df = pd.DataFrame(results)
results_df = results_df.sort_values('max_CF_pct', ascending=False)
results_df.to_csv('AD_Gene_Carrier_Frequencies.csv', index=False)

print(f"\n💾 Detailed results saved: AD_Gene_Carrier_Frequencies.csv")
print(f"   {len(results_df)} AD genes with carrier frequency data")

print(f"\n{'=' * 80}")
