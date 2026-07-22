"""
================================================================================
CARRIER SCREENING ANALYSIS - COMPLETE PIPELINE
================================================================================

Datasets:
  - All of Us (AoU): ~250,071 females with srWGS
  - gnomAD v4.1 JOINT: 807,162 individuals

Scenarios:
  1. AoU P/LP only (star rating ≥2)
  2. AoU P/LP + LoF
  3. AoU P/LP + 5% rare VUS (≥2 stars, AF <0.001%)
  4. gnomAD P/LP only (star rating ≥2)
  5. gnomAD P/LP + 5% rare VUS (≥2 stars, AF <0.001%)

Gene-level carrier frequency formula:
  AF = sum(AC across variants) / mean(AN across variants)
  Autosomal: carrier_freq = 2 * AF * (1 - AF)   [Hardy-Weinberg]
  X chromosome: carrier_freq = 2 * AF * (1 - AF)  [Hardy-Weinberg, females]
  NOTE: AN for X chromosome = n_females * 2 (allele count, not person count),
  so AF is a per-allele frequency and HWE applies identically to autosomal.

5 Gene Panels:
  1. Full 293-gene list
  2. AR + XL genes only
  3. AD-only genes
  4. Full list minus FAM111A, SNAP25, SCN8A, KCNQ3
  5. AD genes minus FAM111A, SNAP25, SCN8A, KCNQ3

Input Files Required (place in same directory as script):
  AoU:
    - PLP_AUTOSOMAL_VARIANT_LEVEL_EXACT_COUNTS.csv
    - PLP_X_CHROMOSOME_FEMALE_EXACT_COUNTS.csv
    - LOF_AUTOSOMAL_VARIANT_LEVEL_EXACT_COUNTS.csv
    - LOF_X_CHROMOSOME_FEMALE_EXACT_COUNTS.csv
    - VUS_AUTOSOMAL_VARIANT_LEVEL_EXACT_COUNTS_RARE.csv
    - VUS_X_CHROMOSOME_FEMALE_EXACT_COUNTS_RARE.csv
  gnomAD:
    - any_PLP_MASTER_JOINT.txt
    - VUS_2plus_stars_JOINT.txt
  Gene list:
    - mmc2_unmerged.xlsx

Output Files:
  - GENE_CLASSIFICATIONS_293.csv
  - AoU_PLP_GENE_LEVEL.csv
  - AoU_LOF_GENE_LEVEL.csv
  - AoU_VUS_GENE_LEVEL_RARE.csv
  - AoU_PLP_LOF_COMBINED_GENE_LEVEL.csv
  - AoU_PLP_VUS_COMBINED_GENE_LEVEL.csv
  - gnomAD_PLP_GENE_LEVEL.csv
  - gnomAD_VUS_GENE_LEVEL_RARE.csv
  - gnomAD_PLP_VUS_COMBINED_GENE_LEVEL.csv
  - NNS_ALL_SCENARIOS.csv
  - NNS_SUMMARY.csv
================================================================================
"""

import pandas as pd
import numpy as np
from scipy import stats
import warnings
warnings.filterwarnings('ignore')
from pathlib import Path

print("=" * 80)
print("CARRIER SCREENING ANALYSIS - COMPLETE PIPELINE")
print("=" * 80)

# ============================================================================
# CONFIGURATION
# ============================================================================

STAR_RATING_THRESHOLD = 2
RARE_VUS_AF_THRESHOLD = 0.00001    # <0.001%
VUS_RECLASSIFICATION_RATE = 0.05   # 5% of rare VUS assumed pathogenic
N_BOOTSTRAP = 10000
GENES_TO_EXCLUDE_FROM_PANELS = ['FAM111A', 'SNAP25', 'SCN8A', 'KCNQ3']

# AoU ancestry columns (from GVS in VAT)
AOU_ANCESTRIES = {
    'ALL': ('gvs_all_ac',  'gvs_all_an'),
    'AFR': ('gvs_afr_ac',  'gvs_afr_an'),
    'AMR': ('gvs_amr_ac',  'gvs_amr_an'),
    'EAS': ('gvs_eas_ac',  'gvs_eas_an'),
    'EUR': ('gvs_eur_ac',  'gvs_eur_an'),
    'MID': ('gvs_mid_ac',  'gvs_mid_an'),
    'OTH': ('gvs_oth_ac',  'gvs_oth_an'),
    'SAS': ('gvs_sas_ac',  'gvs_sas_an'),
}

# gnomAD ancestry columns
GNOMAD_ANCESTRIES = {
    'TOTAL':     ('gnomAD_AC_total',     'gnomAD_AN_total'),
    'AFR':       ('gnomAD_AC_afr',       'gnomAD_AN_afr'),
    'AMR':       ('gnomAD_AC_amr',       'gnomAD_AN_amr'),
    'ASJ':       ('gnomAD_AC_asj',       'gnomAD_AN_asj'),
    'EAS':       ('gnomAD_AC_eas',       'gnomAD_AN_eas'),
    'FIN':       ('gnomAD_AC_fin',       'gnomAD_AN_fin'),
    'NFE':       ('gnomAD_AC_nfe',       'gnomAD_AN_nfe'),
    'REMAINING': ('gnomAD_AC_remaining', 'gnomAD_AN_remaining'),
}

print("\n✅ Configuration loaded")
print(f"   Star rating threshold: ≥{STAR_RATING_THRESHOLD}")
print(f"   Rare VUS AF threshold: <{RARE_VUS_AF_THRESHOLD} ({RARE_VUS_AF_THRESHOLD*100}%)")
print(f"   VUS reclassification rate: {VUS_RECLASSIFICATION_RATE*100}%")

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def wilson_ci(ac, an, confidence=0.95):
    """
    Wilson score confidence interval for allele frequency.
    Returns: (af_estimate, ci_lower, ci_upper)
    """
    if an == 0 or ac < 0:
        return 0.0, 0.0, 0.0
    p = ac / an
    z = stats.norm.ppf((1 + confidence) / 2)
    denominator = 1 + z**2 / an
    center = (p + z**2 / (2 * an)) / denominator
    margin = z * np.sqrt((p * (1 - p) / an) + (z**2 / (4 * an**2))) / denominator
    ci_lower = max(0, center - margin)
    ci_upper = min(1, center + margin)
    return p, ci_lower, ci_upper


def autosomal_carrier_freq(af):
    """Hardy-Weinberg carrier frequency for autosomal gene."""
    return 2 * af * (1 - af)


def calc_gene_level_autosomal(variants_df, gene_col, ancestry_cols, label):
    """
    Calculate gene-level carrier frequencies for autosomal variants.

    Formula:
      AF = sum(AC) / mean(AN)
      carrier_freq = 2 * AF * (1 - AF)
      CI from Wilson score on AC/mean_AN, then converted to carrier freq CI

    Parameters:
      variants_df: DataFrame of variants
      gene_col: column name for gene symbol
      ancestry_cols: dict of {ancestry_label: (ac_col, an_col)}
      label: dataset label for printing

    Returns: DataFrame with gene-level frequencies
    """
    results = []
    genes = variants_df[gene_col].dropna().unique()

    for gene in genes:
        gene_vars = variants_df[variants_df[gene_col] == gene]

        for ancestry, (ac_col, an_col) in ancestry_cols.items():
            if ac_col not in gene_vars.columns or an_col not in gene_vars.columns:
                continue

            ac_vals = pd.to_numeric(gene_vars[ac_col], errors='coerce').fillna(0)
            an_vals = pd.to_numeric(gene_vars[an_col], errors='coerce').dropna()

            if len(an_vals) == 0 or an_vals.mean() == 0:
                continue

            total_ac = ac_vals.sum()
            mean_an = an_vals.mean()

            af, af_ci_lower, af_ci_upper = wilson_ci(total_ac, mean_an)

            carrier_freq = autosomal_carrier_freq(af)
            ci_lower = autosomal_carrier_freq(af_ci_lower)
            ci_upper = autosomal_carrier_freq(af_ci_upper)

            results.append({
                'gene_symbol':      gene,
                'ancestry':         ancestry,
                'n_variants':       len(gene_vars),
                'total_ac':         int(total_ac),
                'mean_an':          round(mean_an, 1),
                'allele_frequency': af,
                'carrier_frequency': carrier_freq,
                'ci_95_lower':      ci_lower,
                'ci_95_upper':      ci_upper,
                'chromosome_type':  'autosomal'
            })

    df = pd.DataFrame(results)
    if len(df) > 0:
        print(f"   {label}: {df['gene_symbol'].nunique()} autosomal genes")
    return df


def calc_gene_level_x(variants_df, gene_col, ancestry_cols, ac_col_override=None,
                       an_col_override=None, label=''):
    """
    Calculate gene-level carrier frequencies for X chromosome variants.

    For X-linked genes, carrier_freq = 2 * AF * (1 - AF) [Hardy-Weinberg].
    AN = n_females * 2 (allele-based denominator) in both AoU and gnomAD,
    so AF is a per-allele frequency and HWE applies identically to autosomal.
    Uses ac_female/an_female columns if available (AoU),
    otherwise falls back to ancestry_cols pattern.

    Returns: DataFrame with gene-level frequencies
    """
    results = []

    if len(variants_df) == 0:
        return pd.DataFrame()

    genes = variants_df[gene_col].dropna().unique()

    # Check if female-specific columns exist (AoU)
    has_female_cols = (ac_col_override is not None and
                       ac_col_override in variants_df.columns and
                       an_col_override in variants_df.columns)

    for gene in genes:
        gene_vars = variants_df[variants_df[gene_col] == gene]

        if has_female_cols:
            # AoU: use female-specific ac/an for ALL ancestry
            ac_vals = pd.to_numeric(gene_vars[ac_col_override], errors='coerce').fillna(0)
            an_vals = pd.to_numeric(gene_vars[an_col_override], errors='coerce').dropna()

            if len(an_vals) > 0 and an_vals.mean() > 0:
                total_ac = ac_vals.sum()
                mean_an = an_vals.mean()
                af, af_ci_lower, af_ci_upper = wilson_ci(total_ac, mean_an)
                # AN = n_females * 2 (allele-based), so HWE applies like autosomal
                carrier_freq = autosomal_carrier_freq(af)
                ci_lower     = autosomal_carrier_freq(af_ci_lower)
                ci_upper     = autosomal_carrier_freq(af_ci_upper)

                results.append({
                    'gene_symbol':      gene,
                    'ancestry':         'ALL',
                    'n_variants':       len(gene_vars),
                    'total_ac':         int(total_ac),
                    'mean_an':          round(mean_an, 1),
                    'allele_frequency': af,
                    'carrier_frequency': carrier_freq,
                    'ci_95_lower':      ci_lower,
                    'ci_95_upper':      ci_upper,
                    'chromosome_type':  'X'
                })

        # Also try ancestry-specific columns (gnomAD, or AoU using GVS)
        for ancestry, (ac_col, an_col) in ancestry_cols.items():
            if has_female_cols and ancestry == 'ALL':
                continue  # Already handled above for AoU
            if ac_col not in gene_vars.columns or an_col not in gene_vars.columns:
                continue

            ac_vals = pd.to_numeric(gene_vars[ac_col], errors='coerce').fillna(0)
            an_vals = pd.to_numeric(gene_vars[an_col], errors='coerce').dropna()

            if len(an_vals) == 0 or an_vals.mean() == 0:
                continue

            total_ac = ac_vals.sum()
            mean_an = an_vals.mean()

            if total_ac == 0:
                continue

            af, af_ci_lower, af_ci_upper = wilson_ci(total_ac, mean_an)
            # AN = n_females * 2 (allele-based), so HWE applies like autosomal
            carrier_freq = autosomal_carrier_freq(af)
            ci_lower     = autosomal_carrier_freq(af_ci_lower)
            ci_upper     = autosomal_carrier_freq(af_ci_upper)

            results.append({
                'gene_symbol':      gene,
                'ancestry':         ancestry,
                'n_variants':       len(gene_vars),
                'total_ac':         int(total_ac),
                'mean_an':          round(mean_an, 1),
                'allele_frequency': af,
                'carrier_frequency': carrier_freq,
                'ci_95_lower':      ci_lower,
                'ci_95_upper':      ci_upper,
                'chromosome_type':  'X'
            })

    df = pd.DataFrame(results)
    if len(df) > 0:
        print(f"   {label}: {df['gene_symbol'].nunique()} X-linked genes")
    return df


def combine_gene_frequencies(df_primary, df_secondary, label='secondary'):
    """
    Combine two gene-level frequency DataFrames.
    For each gene+ancestry: combined_cf = 1 - (1-cf1)*(1-cf2)
    CI: propagated approximately.
    """
    merged = df_primary.merge(
        df_secondary[['gene_symbol', 'ancestry', 'carrier_frequency',
                      'ci_95_lower', 'ci_95_upper', 'n_variants']],
        on=['gene_symbol', 'ancestry'],
        how='outer',
        suffixes=('', f'_{label}')
    )

    cf_col   = f'carrier_frequency_{label}'
    ci_l_col = f'ci_95_lower_{label}'
    ci_u_col = f'ci_95_upper_{label}'
    n_col    = f'n_variants_{label}'

    for col in [cf_col, ci_l_col, ci_u_col, n_col]:
        if col not in merged.columns:
            merged[col] = 0.0

    merged['carrier_frequency'].fillna(0, inplace=True)
    merged[cf_col].fillna(0, inplace=True)
    merged['ci_95_lower'].fillna(0, inplace=True)
    merged['ci_95_upper'].fillna(0, inplace=True)
    merged[ci_l_col].fillna(0, inplace=True)
    merged[ci_u_col].fillna(0, inplace=True)
    merged['n_variants'].fillna(0, inplace=True)
    merged[n_col].fillna(0, inplace=True)

    # Combined using product rule
    merged['carrier_frequency_combined'] = (
        1 - (1 - merged['carrier_frequency']) * (1 - merged[cf_col])
    )
    merged['ci_95_lower_combined'] = (
        1 - (1 - merged['ci_95_lower']) * (1 - merged[ci_l_col])
    )
    merged['ci_95_upper_combined'] = (
        1 - (1 - merged['ci_95_upper']) * (1 - merged[ci_u_col])
    )
    merged['n_variants_combined'] = merged['n_variants'] + merged[n_col]

    # Fill missing gene/ancestry info
    for col in ['gene_symbol', 'ancestry', 'chromosome_type']:
        if col in df_secondary.columns:
            merged[col] = merged[col].combine_first(
                merged.get(f'{col}_{label}', merged[col])
            )

    result = merged[[
        'gene_symbol', 'ancestry', 'n_variants_combined',
        'carrier_frequency_combined', 'ci_95_lower_combined', 'ci_95_upper_combined',
        'chromosome_type'
    ]].rename(columns={
        'n_variants_combined':     'n_variants',
        'carrier_frequency_combined': 'carrier_frequency',
        'ci_95_lower_combined':    'ci_95_lower',
        'ci_95_upper_combined':    'ci_95_upper'
    })

    return result


def calculate_combined_carrier_freq(carrier_freqs):
    """Combined carrier frequency across genes using product rule."""
    valid = [f for f in carrier_freqs if pd.notna(f) and f > 0]
    if not valid:
        return 0.0
    return 1 - np.prod([1 - f for f in valid])


def calculate_nns_with_ci(carrier_freqs, ci_lowers, ci_uppers, n_bootstrap=N_BOOTSTRAP):
    """
    Calculate NNS with 95% bootstrap confidence intervals.
    Returns: (NNS, NNS_95CI_lower, NNS_95CI_upper, combined_carrier_freq)
    """
    combined_freq = calculate_combined_carrier_freq(carrier_freqs)
    if combined_freq == 0:
        return None, None, None, 0.0
    nns = 1 / combined_freq

    bootstrap_combined = []
    for _ in range(n_bootstrap):
        sampled = []
        for cf, cl, cu in zip(carrier_freqs, ci_lowers, ci_uppers):
            if pd.notna(cf) and cf > 0 and pd.notna(cl) and pd.notna(cu):
                if cu > cl:
                    se = (cu - cl) / (2 * 1.96)
                    if se > 0:
                        s = np.random.normal(cf, se)
                        s = max(0.0, min(1.0, s))
                        sampled.append(s)
                    else:
                        sampled.append(cf)
                else:
                    sampled.append(cf)
        if sampled:
            c = calculate_combined_carrier_freq(sampled)
            if c > 0:
                bootstrap_combined.append(c)

    if len(bootstrap_combined) > 100:
        ci_lower_cf = np.percentile(bootstrap_combined, 2.5)
        ci_upper_cf = np.percentile(bootstrap_combined, 97.5)
        nns_lower = 1 / ci_upper_cf if ci_upper_cf > 0 else None
        nns_upper = 1 / ci_lower_cf if ci_lower_cf > 0 else None
    else:
        nns_lower = None
        nns_upper = None

    return nns, nns_lower, nns_upper, combined_freq


# ============================================================================
# PART 1: GENE CLASSIFICATIONS FROM mmc2_unmerged.xlsx
# ============================================================================

print("\n" + "=" * 80)
print("PART 1: GENE CLASSIFICATIONS")
print("=" * 80)

gene_df = pd.read_excel('../raw_files/mmc2_unmerged.xlsx', header=None)
print(f"✅ Loaded mmc2_unmerged.xlsx: {gene_df.shape}")

# Column 0 = gene name, column 4 = inheritance
gene_raw = gene_df[[2, 3]].copy()
gene_raw.columns = ['gene_symbol', 'inheritance']
gene_raw['gene_symbol'] = gene_raw['gene_symbol'].astype(str).str.strip()
gene_raw['inheritance'] = gene_raw['inheritance'].astype(str).str.strip()

# Remove non-gene rows (NaN, empty, numeric)
gene_raw = gene_raw[
    gene_raw['gene_symbol'].notna() &
    (gene_raw['gene_symbol'] != 'nan') &
    (gene_raw['gene_symbol'] != '')
]

# Clean up inheritance: strip whitespace and newlines
gene_raw['inheritance'] = gene_raw['inheritance'].str.replace(r'\n', '', regex=True).str.strip()

# Keep only valid inheritance patterns
# Normalize compound and malformed entries
gene_raw['inheritance'] = gene_raw['inheritance'].str.replace('ARAR', 'AR', regex=False)
gene_raw.loc[gene_raw['inheritance'].isin(['AD;AR', 'AR;AD']), 'inheritance'] = 'AD;AR'

print(f"\nRaw gene entries: {len(gene_raw)}")
print(f"Unique genes: {gene_raw['gene_symbol'].nunique()}")
print(f"\nInheritance distribution (raw):")
print(gene_raw['inheritance'].value_counts().head(10))

# Resolve inheritance per gene
# Rule: if gene has multiple inheritance patterns, AR > XL > AD
def resolve_inheritance(inheritances):
    inh_set = set()
    for i in inheritances:
        for part in str(i).split(';'):
            p = part.strip()
            if p in ['AD', 'AR', 'XL']:
                inh_set.add(p)
    if 'AR' in inh_set:
        return 'AR'
    elif 'XL' in inh_set:
        return 'XL'
    elif 'AD' in inh_set:
        return 'AD'
    return None

gene_class_list = []
for gene, grp in gene_raw.groupby('gene_symbol'):
    inheritances = grp['inheritance'].tolist()
    resolved = resolve_inheritance(inheritances)
    if resolved:
        gene_class_list.append({
            'gene_symbol': gene,
            'inheritance': resolved,
            'original_inheritance': ';'.join(sorted(set(inheritances)))
        })

gene_class = pd.DataFrame(gene_class_list).sort_values('gene_symbol').reset_index(drop=True)

print(f"\n✅ Final gene classifications: {len(gene_class)} genes")
print(f"   AD: {(gene_class['inheritance']=='AD').sum()}")
print(f"   AR: {(gene_class['inheritance']=='AR').sum()}")
print(f"   XL: {(gene_class['inheritance']=='XL').sum()}")

gene_class.to_csv('GENE_CLASSIFICATIONS_293.csv', index=False)
print(f"✅ Saved: GENE_CLASSIFICATIONS_293.csv")

# ============================================================================
# PART 2: DEFINE GENE PANELS
# ============================================================================

print("\n" + "=" * 80)
print("PART 2: GENE PANELS")
print("=" * 80)

all_genes = gene_class['gene_symbol'].tolist()
ar_xl_genes = gene_class[gene_class['inheritance'].isin(['AR', 'XL'])]['gene_symbol'].tolist()
ad_only_genes = gene_class[gene_class['inheritance'] == 'AD']['gene_symbol'].tolist()
all_minus4 = [g for g in all_genes if g not in GENES_TO_EXCLUDE_FROM_PANELS]
ad_minus4 = [g for g in ad_only_genes if g not in GENES_TO_EXCLUDE_FROM_PANELS]

panels = {
    '1_All_293_genes':  all_genes,
    '2_AR_and_XL':      ar_xl_genes,
    '3_AD_only':        ad_only_genes,
    '4_All_minus_4':    all_minus4,
    '5_AD_minus_4':     ad_minus4,
}

print(f"\n✅ Gene panels defined:")
for name, genes in panels.items():
    print(f"   {name}: {len(genes)} genes")

# ============================================================================
# PART 3: ALL OF US - P/LP GENE-LEVEL CARRIER FREQUENCIES
# ============================================================================



def convert_stars_to_numeric(df):
    """Convert star symbols (⭐) to numeric ratings."""
    if 'StarRating' in df.columns and 'StarRating_Numeric' in df.columns:
        # Count stars in StarRating column
        df['StarRating_Numeric'] = df['StarRating'].astype(str).str.count('⭐')
        # Fill NaN with 0
        df['StarRating_Numeric'].fillna(0, inplace=True)
    return df

def filter_star_rating(df, col_name, threshold, label):
    """Filter DataFrame to star rating >= threshold. Handles NaN and string values."""
    if col_name not in df.columns:
        print(f"   ⚠️  No {col_name} column found in {label} - keeping all variants")
        return df
    
    df = df.copy()
    numeric_vals = pd.to_numeric(df[col_name], errors='coerce')
    
    n_nan = numeric_vals.isna().sum()
    n_total = len(df)
    print(f"   {label} star rating: {numeric_vals.value_counts(dropna=False).sort_index().to_dict()}")
    
    if numeric_vals.notna().sum() == 0:
        print(f"   ⚠️  All star rating values are NaN in {label} - keeping all variants")
        return df
    
    filtered = df[numeric_vals >= threshold].copy()
    print(f"   After ≥{threshold} stars: {len(filtered):,} kept, {n_total - len(filtered):,} removed")
    return filtered

print("\n" + "=" * 80)
print("PART 3: ALL OF US - P/LP GENE-LEVEL CARRIER FREQUENCIES")
print("=" * 80)

print("\n📂 Loading AoU P/LP files...")
aou_plp_auto = pd.read_csv('../raw_files/raw_AoU/PLP_AUTOSOMAL_VARIANT_LEVEL_EXACT_COUNTS_renamedPLEC1.csv')
aou_plp_x    = pd.read_csv('../raw_files/raw_AoU/PLP_X_CHROMOSOME_FEMALE_EXACT_COUNTS.csv')

# Convert star symbols to numeric ratings
aou_plp_auto = convert_stars_to_numeric(aou_plp_auto)
aou_plp_x    = convert_stars_to_numeric(aou_plp_x)

print(f"   Autosomal P/LP: {len(aou_plp_auto):,} variants")
print(f"   X chromosome P/LP: {len(aou_plp_x):,} variants")

# Filter to star rating ≥2
star_col = next((c for c in ['StarRating_Numeric', 'star_rating_numeric', 'StarRating', 'star_rating', 'clinvar_star_rating'] if c in aou_plp_auto.columns), None)
if star_col:
    aou_plp_auto = filter_star_rating(aou_plp_auto, star_col, STAR_RATING_THRESHOLD, 'AoU P/LP autosomal')
    aou_plp_x    = filter_star_rating(aou_plp_x,    star_col if star_col in aou_plp_x.columns else next((c for c in ['StarRating_Numeric','star_rating_numeric','StarRating'] if c in aou_plp_x.columns), None) or star_col, STAR_RATING_THRESHOLD, 'AoU P/LP X chromosome')
else:
    print(f"   ⚠️  No star rating column found - print columns to diagnose:")
    print(f"   {list(aou_plp_auto.columns[:15])}")

# Calculate gene-level frequencies
print("\n🧬 Calculating AoU P/LP gene-level carrier frequencies...")
aou_plp_auto_gene = calc_gene_level_autosomal(
    aou_plp_auto, 'gene_symbol', AOU_ANCESTRIES, 'AoU P/LP autosomal'
)
aou_plp_x_gene = calc_gene_level_x(
    aou_plp_x, 'gene_symbol', AOU_ANCESTRIES,
    ac_col_override='ac_female', an_col_override='an_female',
    label='AoU P/LP X chromosome'
)

aou_plp_gene = pd.concat([aou_plp_auto_gene, aou_plp_x_gene], ignore_index=True)
aou_plp_gene.to_csv('AoU_PLP_GENE_LEVEL.csv', index=False)
print(f"✅ Saved: AoU_PLP_GENE_LEVEL.csv ({aou_plp_gene['gene_symbol'].nunique()} genes)")

# ============================================================================
# PART 4: ALL OF US - LOF GENE-LEVEL CARRIER FREQUENCIES
# ============================================================================

print("\n" + "=" * 80)
print("PART 4: ALL OF US - LOF GENE-LEVEL CARRIER FREQUENCIES")
print("=" * 80)

print("\n📂 Loading AoU LoF files...")
aou_lof_auto = pd.read_csv('../raw_files/raw_AoU/LOF_AUTOSOMAL_VARIANT_LEVEL_EXACT_COUNTS.csv')
aou_lof_x    = pd.read_csv('../raw_files/raw_AoU/LOF_X_CHROMOSOME_FEMALE_EXACT_COUNTS.csv')

# Convert star symbols to numeric ratings
aou_lof_auto = convert_stars_to_numeric(aou_lof_auto)
aou_lof_x    = convert_stars_to_numeric(aou_lof_x)

print(f"   Autosomal LoF: {len(aou_lof_auto):,} variants")
print(f"   X chromosome LoF: {len(aou_lof_x):,} variants")

print("\n🧬 Calculating AoU LoF gene-level carrier frequencies...")
aou_lof_auto_gene = calc_gene_level_autosomal(
    aou_lof_auto, 'gene_symbol', AOU_ANCESTRIES, 'AoU LoF autosomal'
)
aou_lof_x_gene = calc_gene_level_x(
    aou_lof_x, 'gene_symbol', AOU_ANCESTRIES,
    ac_col_override='ac_female', an_col_override='an_female',
    label='AoU LoF X chromosome'
)

aou_lof_gene = pd.concat([aou_lof_auto_gene, aou_lof_x_gene], ignore_index=True)
aou_lof_gene.to_csv('AoU_LOF_GENE_LEVEL.csv', index=False)
print(f"✅ Saved: AoU_LOF_GENE_LEVEL.csv ({aou_lof_gene['gene_symbol'].nunique()} genes)")

# ============================================================================
# PART 5: ALL OF US - VUS GENE-LEVEL CARRIER FREQUENCIES
# ============================================================================

print("\n" + "=" * 80)
print("PART 5: ALL OF US - VUS GENE-LEVEL CARRIER FREQUENCIES")
print("=" * 80)

print("\n📂 Loading AoU VUS files (already filtered to ≥2 stars and AF <0.001%)...")
aou_vus_auto = pd.read_csv('../raw_files/raw_AoU/VUS_AUTOSOMAL_VARIANT_LEVEL_EXACT_COUNTS_RARE_renamedPLEC1.csv')
aou_vus_x    = pd.read_csv('../raw_files/raw_AoU/VUS_X_CHROMOSOME_FEMALE_EXACT_COUNTS_RARE.csv')

# Convert star symbols to numeric ratings
aou_vus_auto = convert_stars_to_numeric(aou_vus_auto)
aou_vus_x    = convert_stars_to_numeric(aou_vus_x)

print(f"   Autosomal VUS (rare): {len(aou_vus_auto):,} variants")
print(f"   X chromosome VUS (rare): {len(aou_vus_x):,} variants")

# Double-check star rating filter (should already be applied from AoU extraction)
vus_star_col = next((c for c in ['StarRating_Numeric', 'star_rating_numeric', 'StarRating'] if c in aou_vus_auto.columns), None)
if vus_star_col:
    n_before = len(aou_vus_auto)
    aou_vus_auto = filter_star_rating(aou_vus_auto, vus_star_col, STAR_RATING_THRESHOLD, 'AoU VUS autosomal')
    if n_before != len(aou_vus_auto):
        print(f"   Additional star filter applied")

print("\n🧬 Calculating AoU VUS gene-level carrier frequencies...")
aou_vus_auto_gene = calc_gene_level_autosomal(
    aou_vus_auto, 'gene_symbol', AOU_ANCESTRIES, 'AoU VUS autosomal'
)
aou_vus_x_gene = calc_gene_level_x(
    aou_vus_x, 'gene_symbol', AOU_ANCESTRIES,
    ac_col_override='ac_female', an_col_override='an_female',
    label='AoU VUS X chromosome'
)

aou_vus_gene = pd.concat([aou_vus_auto_gene, aou_vus_x_gene], ignore_index=True)

# Apply 5% reclassification rate
aou_vus_gene['carrier_frequency_vus_5pct'] = (
    aou_vus_gene['carrier_frequency'] * VUS_RECLASSIFICATION_RATE
)
aou_vus_gene.to_csv('AoU_VUS_GENE_LEVEL_RARE.csv', index=False)
print(f"✅ Saved: AoU_VUS_GENE_LEVEL_RARE.csv ({aou_vus_gene['gene_symbol'].nunique()} genes)")

# ============================================================================
# PART 6: ALL OF US - COMBINED GENE LEVELS
# ============================================================================

print("\n" + "=" * 80)
print("PART 6: ALL OF US - COMBINED GENE-LEVEL FREQUENCIES")
print("=" * 80)

# P/LP + LoF
print("\n🔗 Combining P/LP + LoF...")
aou_plp_lof_gene = combine_gene_frequencies(aou_plp_gene, aou_lof_gene, label='lof')
aou_plp_lof_gene.to_csv('AoU_PLP_LOF_COMBINED_GENE_LEVEL.csv', index=False)
print(f"✅ Saved: AoU_PLP_LOF_COMBINED_GENE_LEVEL.csv ({aou_plp_lof_gene['gene_symbol'].nunique()} genes)")

# P/LP + 5% VUS
print("\n🔗 Combining P/LP + 5% VUS...")
# Create VUS gene-level with 5% applied
aou_vus_5pct = aou_vus_gene.copy()
aou_vus_5pct['carrier_frequency'] = aou_vus_5pct['carrier_frequency_vus_5pct']
aou_vus_5pct['ci_95_lower'] = aou_vus_5pct['ci_95_lower'] * VUS_RECLASSIFICATION_RATE
aou_vus_5pct['ci_95_upper'] = aou_vus_5pct['ci_95_upper'] * VUS_RECLASSIFICATION_RATE

aou_plp_vus_gene = combine_gene_frequencies(aou_plp_gene, aou_vus_5pct, label='vus')
aou_plp_vus_gene.to_csv('AoU_PLP_VUS_COMBINED_GENE_LEVEL.csv', index=False)
print(f"✅ Saved: AoU_PLP_VUS_COMBINED_GENE_LEVEL.csv ({aou_plp_vus_gene['gene_symbol'].nunique()} genes)")

# ============================================================================
# PART 7: gnomAD - P/LP GENE-LEVEL CARRIER FREQUENCIES
# ============================================================================

print("\n" + "=" * 80)
print("PART 7: gnomAD - P/LP GENE-LEVEL CARRIER FREQUENCIES")
print("=" * 80)

print("\n📂 Loading gnomAD P/LP file...")
gnomad_plp_all = pd.read_csv('../raw_files/raw_gnomAD/any_PLP_MASTER_JOINT_renamedPLEC1.csv')
print(f"   Loaded: {len(gnomad_plp_all):,} variants")

# Filter to star rating ≥2
if 'StarRating_Numeric' in gnomad_plp_all.columns:
    gnomad_plp = gnomad_plp_all[
        pd.to_numeric(gnomad_plp_all['StarRating_Numeric'], errors='coerce') >= STAR_RATING_THRESHOLD
    ].copy()
elif 'StarRating' in gnomad_plp_all.columns:
    gnomad_plp_all['StarRating_Numeric'] = pd.to_numeric(
        gnomad_plp_all['StarRating'].str.extract(r'(\d+)', expand=False), errors='coerce'
    )
    gnomad_plp = gnomad_plp_all[gnomad_plp_all['StarRating_Numeric'] >= STAR_RATING_THRESHOLD].copy()
else:
    print("   ⚠️  No StarRating_Numeric column found - using all variants")
    gnomad_plp = gnomad_plp_all.copy()

print(f"   After star rating ≥{STAR_RATING_THRESHOLD}: {len(gnomad_plp):,} variants "
      f"({len(gnomad_plp_all)-len(gnomad_plp):,} removed)")


# DEBUG: Print columns to help identify correct names
print(f"   Available columns: {list(gnomad_plp.columns[:20])}...")

# Identify gene column
# Identify gene column (try multiple variations)
gene_col_gnomad = None
for col_name in ['Gene', 'gene', 'gene_symbol', 'GENE', 'Gene_Symbol', 'gene_name']:
    if col_name in gnomad_plp.columns:
        gene_col_gnomad = col_name
        break

if gene_col_gnomad is None:
    print("   ❌ ERROR: No gene column found!")
    print(f"   Available columns: {list(gnomad_plp.columns)}")
    raise ValueError("Cannot find gene column in gnomAD file")

print(f"   Using gene column: {gene_col_gnomad}")

# Separate autosomal and X chromosome
# Identify chromosome column (try multiple variations)
chrom_col = None
for col_name in ['Chromosome', 'chromosome', 'CHROMOSOME', 'Chr', 'chr', 'CHROM', 'chrom']:
    if col_name in gnomad_plp.columns:
        chrom_col = col_name
        break

if chrom_col is None:
    print("   ⚠️  No chromosome column found - assuming all variants are autosomal")
    gnomad_plp['is_x'] = False
else:
    print(f"   Using chromosome column: {chrom_col}")
    gnomad_plp['is_x'] = gnomad_plp[chrom_col].astype(str).isin(['X', 'chrX', 'x', '23'])
gnomad_plp_auto = gnomad_plp[~gnomad_plp['is_x']].copy()
gnomad_plp_x    = gnomad_plp[gnomad_plp['is_x']].copy()
print(f"   Autosomal: {len(gnomad_plp_auto):,}, X chromosome: {len(gnomad_plp_x):,}")

# Check for XX-specific columns for X chromosome
has_xx_cols = any('XX' in col for col in gnomad_plp_x.columns)
if has_xx_cols:
    print("   ✅ XX-specific columns found for X chromosome")
    gnomad_x_ancestry = {
        anc: (f'gnomAD_AC_XX_{anc.lower()}' if anc != 'TOTAL' else 'gnomAD_AC_XX',
              f'gnomAD_AN_XX_{anc.lower()}' if anc != 'TOTAL' else 'gnomAD_AN_XX')
        for anc in GNOMAD_ANCESTRIES.keys()
    }
    # Only keep columns that exist
    gnomad_x_ancestry = {
        k: v for k, v in gnomad_x_ancestry.items()
        if v[0] in gnomad_plp_x.columns
    }
    if not gnomad_x_ancestry:
        gnomad_x_ancestry = GNOMAD_ANCESTRIES
else:
    print("   ℹ️  No XX columns found - using standard gnomAD counts for X chromosome")
    gnomad_x_ancestry = GNOMAD_ANCESTRIES

print("\n🧬 Calculating gnomAD P/LP gene-level carrier frequencies...")
gnomad_plp_auto_gene = calc_gene_level_autosomal(
    gnomad_plp_auto, gene_col_gnomad, GNOMAD_ANCESTRIES, 'gnomAD P/LP autosomal'
)
gnomad_plp_auto_gene.rename(columns={'gene_symbol': 'gene_symbol'}, inplace=True)

gnomad_plp_x_gene = calc_gene_level_x(
    gnomad_plp_x, gene_col_gnomad, gnomad_x_ancestry,
    label='gnomAD P/LP X chromosome'
)

gnomad_plp_gene = pd.concat([gnomad_plp_auto_gene, gnomad_plp_x_gene], ignore_index=True)
# Standardize gene column name
if gene_col_gnomad != 'gene_symbol':
    gnomad_plp_gene.rename(columns={gene_col_gnomad: 'gene_symbol'}, inplace=True)

gnomad_plp_gene.to_csv('gnomAD_PLP_GENE_LEVEL.csv', index=False)
print(f"✅ Saved: gnomAD_PLP_GENE_LEVEL.csv ({gnomad_plp_gene['gene_symbol'].nunique()} genes)")

# ============================================================================
# PART 8: gnomAD - VUS GENE-LEVEL CARRIER FREQUENCIES
# ============================================================================

print("\n" + "=" * 80)
print("PART 8: gnomAD - VUS GENE-LEVEL CARRIER FREQUENCIES")
print("=" * 80)

print("\n📂 Loading gnomAD VUS file (already ≥2 stars)...")
gnomad_vus_all = pd.read_csv('../raw_files/raw_gnomAD/VUS_2plus_stars_JOINT_renamedPLEC1.csv')
print(f"   Loaded: {len(gnomad_vus_all):,} variants")

# Filter to rare VUS (AF < 0.001%)
af_col_gnomad_vus = 'gnomAD_AF_total' if 'gnomAD_AF_total' in gnomad_vus_all.columns else 'gnomAD_AF_TOTAL'
if af_col_gnomad_vus in gnomad_vus_all.columns:
    gnomad_vus_all[af_col_gnomad_vus] = pd.to_numeric(
        gnomad_vus_all[af_col_gnomad_vus], errors='coerce'
    ).fillna(0)
    gnomad_vus = gnomad_vus_all[
        gnomad_vus_all[af_col_gnomad_vus] < RARE_VUS_AF_THRESHOLD
    ].copy()
    print(f"   After AF <{RARE_VUS_AF_THRESHOLD}: {len(gnomad_vus):,} rare VUS "
          f"({len(gnomad_vus_all)-len(gnomad_vus):,} removed)")
else:
    print(f"   ⚠️  AF column not found ({af_col_gnomad_vus}) - using all VUS")
    gnomad_vus = gnomad_vus_all.copy()

# Identify gene column (try multiple variations)
gene_col_gnomad_vus = None
for col_name in ['Gene', 'gene', 'gene_symbol', 'GENE', 'Gene_Symbol', 'gene_name']:
    if col_name in gnomad_vus.columns:
        gene_col_gnomad_vus = col_name
        break
if gene_col_gnomad_vus is None:
    raise ValueError("Cannot find gene column in gnomAD VUS file")
print(f"   Using gene column: {gene_col_gnomad_vus}")

# Identify chromosome column
chrom_col_vus = None
for col_name in ['Chromosome', 'chromosome', 'CHROMOSOME', 'Chr', 'chr', 'CHROM', 'chrom']:
    if col_name in gnomad_vus.columns:
        chrom_col_vus = col_name
        break
if chrom_col_vus is None:
    print("   ⚠️  No chromosome column - assuming all autosomal")
    gnomad_vus['is_x'] = False
else:
    print(f"   Using chromosome column: {chrom_col_vus}")
    gnomad_vus['is_x'] = gnomad_vus[chrom_col_vus].astype(str).isin(['X', 'chrX', 'x', '23'])
gnomad_vus_auto = gnomad_vus[~gnomad_vus['is_x']].copy()
gnomad_vus_x    = gnomad_vus[gnomad_vus['is_x']].copy()
print(f"   Autosomal rare VUS: {len(gnomad_vus_auto):,}, X chromosome: {len(gnomad_vus_x):,}")

print("\n🧬 Calculating gnomAD VUS gene-level carrier frequencies...")
gnomad_vus_auto_gene = calc_gene_level_autosomal(
    gnomad_vus_auto, gene_col_gnomad_vus, GNOMAD_ANCESTRIES, 'gnomAD VUS autosomal'
)
gnomad_vus_x_gene = calc_gene_level_x(
    gnomad_vus_x, gene_col_gnomad_vus, gnomad_x_ancestry,
    label='gnomAD VUS X chromosome'
)

gnomad_vus_gene = pd.concat([gnomad_vus_auto_gene, gnomad_vus_x_gene], ignore_index=True)
if gene_col_gnomad_vus != 'gene_symbol' and gene_col_gnomad_vus in gnomad_vus_gene.columns:
    gnomad_vus_gene.rename(columns={gene_col_gnomad_vus: 'gene_symbol'}, inplace=True)

gnomad_vus_gene['carrier_frequency_vus_5pct'] = (
    gnomad_vus_gene['carrier_frequency'] * VUS_RECLASSIFICATION_RATE
)
gnomad_vus_gene.to_csv('gnomAD_VUS_GENE_LEVEL_RARE.csv', index=False)
print(f"✅ Saved: gnomAD_VUS_GENE_LEVEL_RARE.csv ({gnomad_vus_gene['gene_symbol'].nunique()} genes)")

# P/LP + 5% VUS combined
print("\n🔗 Combining gnomAD P/LP + 5% VUS...")
gnomad_vus_5pct = gnomad_vus_gene.copy()
gnomad_vus_5pct['carrier_frequency'] = gnomad_vus_5pct['carrier_frequency_vus_5pct']
gnomad_vus_5pct['ci_95_lower'] = gnomad_vus_5pct['ci_95_lower'] * VUS_RECLASSIFICATION_RATE
gnomad_vus_5pct['ci_95_upper'] = gnomad_vus_5pct['ci_95_upper'] * VUS_RECLASSIFICATION_RATE

gnomad_plp_vus_gene = combine_gene_frequencies(gnomad_plp_gene, gnomad_vus_5pct, label='vus')
gnomad_plp_vus_gene.to_csv('gnomAD_PLP_VUS_COMBINED_GENE_LEVEL.csv', index=False)
print(f"✅ Saved: gnomAD_PLP_VUS_COMBINED_GENE_LEVEL.csv ({gnomad_plp_vus_gene['gene_symbol'].nunique()} genes)")

# ============================================================================
# PART 9: NNS CALCULATIONS FOR ALL SCENARIOS AND PANELS
# ============================================================================

print("\n" + "=" * 80)
print("PART 9: NNS CALCULATIONS")
print("=" * 80)

# Define all scenarios
scenarios = [
    # (label, dataset_label, gene_level_df, ancestry_list)
    ('AoU_PLP_only',       'All of Us',  aou_plp_gene,     list(AOU_ANCESTRIES.keys())),
    ('AoU_PLP_plus_LoF',   'All of Us',  aou_plp_lof_gene, list(AOU_ANCESTRIES.keys())),
    ('AoU_PLP_plus_5pctVUS','All of Us', aou_plp_vus_gene, list(AOU_ANCESTRIES.keys())),
    ('gnomAD_PLP_only',    'gnomAD',     gnomad_plp_gene,  list(GNOMAD_ANCESTRIES.keys())),
    ('gnomAD_PLP_plus_5pctVUS','gnomAD', gnomad_plp_vus_gene, list(GNOMAD_ANCESTRIES.keys())),
]

all_nns_results = []

for scenario_label, dataset, gene_level_df, ancestries in scenarios:
    print(f"\n📊 Scenario: {scenario_label}")

    # Merge with inheritance
    gene_with_inh = gene_level_df.merge(
        gene_class[['gene_symbol', 'inheritance']],
        on='gene_symbol',
        how='left'
    )

    for panel_name, panel_genes in panels.items():
        for ancestry in ancestries:
            panel_data = gene_with_inh[
                (gene_with_inh['gene_symbol'].isin(panel_genes)) &
                (gene_with_inh['ancestry'] == ancestry) &
                (gene_with_inh['carrier_frequency'].notna()) &
                (gene_with_inh['carrier_frequency'] > 0)
            ]

            carrier_freqs = panel_data['carrier_frequency'].tolist()
            ci_lowers = panel_data['ci_95_lower'].tolist()
            ci_uppers = panel_data['ci_95_upper'].tolist()

            nns, nns_lower, nns_upper, combined_freq = calculate_nns_with_ci(
                carrier_freqs, ci_lowers, ci_uppers
            )

            all_nns_results.append({
                'scenario':           scenario_label,
                'dataset':            dataset,
                'panel':              panel_name,
                'ancestry':           ancestry,
                'n_genes_in_panel':   len(panel_genes),
                'n_genes_with_data':  len(panel_data),
                'combined_carrier_freq': combined_freq,
                'NNS':                nns,
                'NNS_95CI_lower':     nns_lower,
                'NNS_95CI_upper':     nns_upper,
                '1_in_N':             f"1 in {nns:,.0f}" if nns else "N/A"
            })

    # Summary for this scenario
    scenario_data = [r for r in all_nns_results if r['scenario'] == scenario_label]
    print(f"   Calculated {len(scenario_data)} NNS values ({len(panels)} panels × {len(ancestries)} ancestries)")

# Save results
nns_df = pd.DataFrame(all_nns_results)
nns_df.to_csv('NNS_ALL_SCENARIOS.csv', index=False)
print(f"\n✅ Saved: NNS_ALL_SCENARIOS.csv ({len(nns_df):,} rows)")

# ============================================================================
# PART 10: SUMMARY TABLE
# ============================================================================

print("\n" + "=" * 80)
print("PART 10: SUMMARY TABLE")
print("=" * 80)

# Create readable summary focusing on ALL/TOTAL ancestry
summary_ancestry = {
    'AoU_PLP_only':          'ALL',
    'AoU_PLP_plus_LoF':      'ALL',
    'AoU_PLP_plus_5pctVUS':  'ALL',
    'gnomAD_PLP_only':       'TOTAL',
    'gnomAD_PLP_plus_5pctVUS':'TOTAL',
}

summary_rows = []
for scenario_label, ancestry in summary_ancestry.items():
    for panel_name in panels.keys():
        row = nns_df[
            (nns_df['scenario'] == scenario_label) &
            (nns_df['panel'] == panel_name) &
            (nns_df['ancestry'] == ancestry)
        ]
        if len(row) > 0:
            r = row.iloc[0]
            summary_rows.append({
                'Scenario':             scenario_label,
                'Panel':                panel_name,
                'Ancestry':             ancestry,
                'N_Genes_With_Data':    r['n_genes_with_data'],
                'Combined_CF':          round(r['combined_carrier_freq'], 6) if r['combined_carrier_freq'] else None,
                'Combined_CF_Pct':      f"{r['combined_carrier_freq']*100:.2f}%" if r['combined_carrier_freq'] else None,
                'NNS':                  round(r['NNS'], 1) if r['NNS'] else None,
                'NNS_95CI_lower':       round(r['NNS_95CI_lower'], 1) if r['NNS_95CI_lower'] else None,
                'NNS_95CI_upper':       round(r['NNS_95CI_upper'], 1) if r['NNS_95CI_upper'] else None,
                '1_in_N':               r['1_in_N']
            })

summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv('NNS_SUMMARY.csv', index=False)

# Print summary
print("\n📊 NNS SUMMARY (All scenarios, combined ancestry, Panel 2: AR+XL):")
print("-" * 80)
ar_xl_summary = summary_df[summary_df['Panel'] == '2_AR_and_XL']
for _, row in ar_xl_summary.iterrows():
    ci_str = ""
    if row['NNS_95CI_lower'] and row['NNS_95CI_upper']:
        ci_str = f" (95% CI: {row['NNS_95CI_lower']:,.1f} - {row['NNS_95CI_upper']:,.1f})"
    nns_str = f"{row['NNS']:,.1f}" if row['NNS'] else "N/A"
    cf_str = row['Combined_CF_Pct'] if row['Combined_CF_Pct'] else "N/A"
    print(f"  {row['Scenario']:<30} NNS={nns_str}{ci_str} | CF={cf_str}")

print(f"\n✅ Saved: NNS_SUMMARY.csv ({len(summary_df):,} rows)")

# ============================================================================
# FINAL SUMMARY
# ============================================================================

print("\n" + "=" * 80)
print("✅ ANALYSIS COMPLETE!")
print("=" * 80)

print(f"""
📊 OUTPUT FILES:

Gene Classifications:
  ✅ GENE_CLASSIFICATIONS_293.csv

All of Us Gene-Level Data:
  ✅ AoU_PLP_GENE_LEVEL.csv
  ✅ AoU_LOF_GENE_LEVEL.csv
  ✅ AoU_VUS_GENE_LEVEL_RARE.csv
  ✅ AoU_PLP_LOF_COMBINED_GENE_LEVEL.csv
  ✅ AoU_PLP_VUS_COMBINED_GENE_LEVEL.csv

gnomAD Gene-Level Data:
  ✅ gnomAD_PLP_GENE_LEVEL.csv
  ✅ gnomAD_VUS_GENE_LEVEL_RARE.csv
  ✅ gnomAD_PLP_VUS_COMBINED_GENE_LEVEL.csv

NNS Results:
  ✅ NNS_ALL_SCENARIOS.csv  ({len(nns_df):,} rows: 5 scenarios × 5 panels × ancestries)
  ✅ NNS_SUMMARY.csv        (Summary for combined ancestry)

Gene Panels:
  1. All 293 genes:     {len(panels['1_All_293_genes'])} genes
  2. AR + XL:           {len(panels['2_AR_and_XL'])} genes
  3. AD only:           {len(panels['3_AD_only'])} genes
  4. All minus 4 genes: {len(panels['4_All_minus_4'])} genes
  5. AD minus 4 genes:  {len(panels['5_AD_minus_4'])} genes

Gene-Level Formula: sum(AC) / mean(AN) → AF → HWE carrier freq (autosomal & X-linked)
                    X-linked: AN = n_females × 2, so HWE applies identically to autosomal
Wilson score 95% CI, Bootstrap NNS CI (10,000 iterations)
""")
