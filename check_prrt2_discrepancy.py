"""
Investigate PRRT2 carrier frequency discrepancy between AoU and gnomAD.

AoU CF: ~0.05% (ranks 5th in AD panel)
gnomAD CF: ~0.46% (ranks 1st in AD panel)
~10-fold difference

Checks:
1. Which PRRT2 variants are in each dataset
2. Whether the same variants exist in both
3. Whether AF differs for shared variants
4. Whether gnomAD has additional high-AF variants not in AoU
"""

import pandas as pd
import os

# ── Paths ──────────────────────────────────────────────────────────────────
BASE = os.path.expanduser(
    '~/Desktop/Figures_for_Aim1_publication 2/send_to_coauthors/raw_files'
)
AOU_FILE    = os.path.join(BASE, 'raw_AoU',
              'PLP_AUTOSOMAL_VARIANT_LEVEL_EXACT_COUNTS_renamedPLEC1.csv')
GNOMAD_FILE = os.path.join(BASE, 'raw_gnomAD',
              'any_PLP_MASTER_JOINT_renamedPLEC1.csv')

print("=" * 70)
print("PRRT2 DISCREPANCY INVESTIGATION")
print("=" * 70)

# ── Load and filter to PRRT2 ───────────────────────────────────────────────
aou   = pd.read_csv(AOU_FILE,    low_memory=False)
gno   = pd.read_csv(GNOMAD_FILE, low_memory=False)

aou_p = aou[aou['gene_symbol'].astype(str).str.upper() == 'PRRT2'].copy()
gno_p = gno[gno['Gene'].astype(str).str.upper()        == 'PRRT2'].copy()

# gnomAD: apply >=2 star filter (same as pipeline)
if 'StarRating' in gno_p.columns:
    gno_p = gno_p[pd.to_numeric(gno_p['StarRating'], errors='coerce') >= 2]
elif 'StarRating_Numeric' in gno_p.columns:
    gno_p = gno_p[pd.to_numeric(gno_p['StarRating_Numeric'], errors='coerce') >= 2]

print(f"\nAoU  PRRT2 variants (after star filter): {len(aou_p)}")
print(f"gnomAD PRRT2 variants (after star filter): {len(gno_p)}")

# ── AoU variant-level summary ─────────────────────────────────────────────
print("\n── AoU PRRT2 variants ──────────────────────────────────────────────")
aou_cols = ['variant_id' if 'variant_id' in aou_p.columns else
            'NC_Variant_ID' if 'NC_Variant_ID' in aou_p.columns else
            aou_p.columns[0]]

# Find AC and AN columns
ac_col_aou = next((c for c in aou_p.columns
                   if c.lower() in ['ac','ac_total','allele_count']), None)
an_col_aou = next((c for c in aou_p.columns
                   if c.lower() in ['an','an_total','allele_number']), None)
af_col_aou = next((c for c in aou_p.columns
                   if c.lower() in ['af','af_total','allele_frequency']), None)

id_col_aou = next((c for c in ['NC_Variant_ID','variant_id','Position']
                   if c in aou_p.columns), aou_p.columns[0])

show_cols = [id_col_aou]
if ac_col_aou: show_cols.append(ac_col_aou)
if an_col_aou: show_cols.append(an_col_aou)
if af_col_aou: show_cols.append(af_col_aou)

print(aou_p[show_cols].to_string(index=False))

if ac_col_aou and an_col_aou:
    total_ac = pd.to_numeric(aou_p[ac_col_aou], errors='coerce').sum()
    mean_an  = pd.to_numeric(aou_p[an_col_aou], errors='coerce').mean()
    af_gene  = total_ac / mean_an if mean_an > 0 else 0
    cf_gene  = 2 * af_gene * (1 - af_gene)
    print(f"\nAoU gene-level: AC={total_ac:.0f}, mean AN={mean_an:.0f}, "
          f"AF={af_gene:.6f}, CF={cf_gene*100:.4f}%")

# ── gnomAD variant-level summary ──────────────────────────────────────────
print("\n── gnomAD PRRT2 variants ───────────────────────────────────────────")

id_col_gno = next((c for c in ['NC_Variant_ID','Position']
                   if c in gno_p.columns), gno_p.columns[0])

gno_show = [id_col_gno, 'ClinicalSignificance',
            'gnomAD_AC_total', 'gnomAD_AN_total', 'gnomAD_AF_total']
gno_show = [c for c in gno_show if c in gno_p.columns]

gno_p_sorted = gno_p.sort_values('gnomAD_AF_total', ascending=False)
print(gno_p_sorted[gno_show].to_string(index=False))

# Gene-level gnomAD CF
ac_gno = pd.to_numeric(gno_p['gnomAD_AC_total'], errors='coerce').fillna(0)
an_gno = pd.to_numeric(gno_p['gnomAD_AN_total'], errors='coerce').dropna()
if len(an_gno) > 0:
    total_ac_gno = ac_gno.sum()
    mean_an_gno  = an_gno.mean()
    af_gno       = total_ac_gno / mean_an_gno if mean_an_gno > 0 else 0
    cf_gno       = 2 * af_gno * (1 - af_gno)
    print(f"\ngnomAD gene-level: AC={total_ac_gno:.0f}, mean AN={mean_an_gno:.0f}, "
          f"AF={af_gno:.6f}, CF={cf_gno*100:.4f}%")

# ── Overlap analysis ───────────────────────────────────────────────────────
print("\n── Variant overlap ─────────────────────────────────────────────────")

# Standardise position-based key for matching
def make_key(df, chrom_col, pos_col, ref_col, alt_col):
    return (df[chrom_col].astype(str).str.replace('chr','') + ':' +
            df[pos_col].astype(str) + ':' +
            df[ref_col].astype(str) + ':' +
            df[alt_col].astype(str))

aou_key_cols  = ['chromosome','position','ref','alt'] if 'chromosome' in aou_p.columns \
                else ['Chromosome','Position','Ref','Alt']
gno_key_cols  = ['Chromosome','Position','Ref','Alt']

try:
    aou_p['_key'] = make_key(aou_p, *aou_key_cols)
    gno_p['_key'] = make_key(gno_p, *gno_key_cols)

    aou_keys = set(aou_p['_key'])
    gno_keys = set(gno_p['_key'])

    shared     = aou_keys & gno_keys
    aou_only   = aou_keys - gno_keys
    gno_only   = gno_keys - aou_keys

    print(f"Shared variants (same position/ref/alt): {len(shared)}")
    print(f"AoU only:    {len(aou_only)}")
    print(f"gnomAD only: {len(gno_only)}")

    if gno_only:
        print(f"\ngnomAD-only PRRT2 variants (not in AoU):")
        gno_unique = gno_p[gno_p['_key'].isin(gno_only)].sort_values(
                     'gnomAD_AF_total', ascending=False)
        print(gno_unique[gno_show].to_string(index=False))
        print(f"\ngnomAD-only variants AC total: "
              f"{pd.to_numeric(gno_unique['gnomAD_AC_total'], errors='coerce').sum():.0f}")
        print(f"Contribution to gnomAD CF: "
              f"these variants explain the majority of the gnomAD-AoU discrepancy "
              f"if their AC >> AoU AC")

    if shared:
        print(f"\nShared variants — AF comparison:")
        aou_shared  = aou_p[aou_p['_key'].isin(shared)][
                      [id_col_aou, '_key'] + ([af_col_aou] if af_col_aou else [])
                      ].rename(columns={af_col_aou: 'AoU_AF'} if af_col_aou else {})
        gno_shared  = gno_p[gno_p['_key'].isin(shared)][
                      ['_key','gnomAD_AF_total']
                      ].rename(columns={'gnomAD_AF_total': 'gnomAD_AF'})
        merged = aou_shared.merge(gno_shared, on='_key')
        if 'AoU_AF' in merged.columns:
            merged['AF_ratio_gnomAD_over_AoU'] = (
                pd.to_numeric(merged['gnomAD_AF'], errors='coerce') /
                pd.to_numeric(merged['AoU_AF'],   errors='coerce')
            )
        print(merged.to_string(index=False))

except KeyError as e:
    print(f"Column mismatch for overlap analysis: {e}")
    print(f"AoU columns:   {list(aou_p.columns[:15])}")
    print(f"gnomAD columns:{list(gno_p.columns[:15])}")

print("\n" + "=" * 70)
print("INTERPRETATION GUIDE")
print("=" * 70)
print("""
If gnomAD has many variants NOT in AoU:
  → AoU's star-rating filter or variant extraction captured fewer PRRT2 variants
  → Check ClinVar star ratings of gnomAD-only variants

If shared variants have similar AFs:
  → Discrepancy is driven by gnomAD having more variants, not different AFs

If shared variants have much higher gnomAD AF:
  → True population frequency difference (ancestry composition)

If gnomAD AC >> AoU AC for same variants:
  → PRRT2 variants are genuinely more common in gnomAD reference population
""")
