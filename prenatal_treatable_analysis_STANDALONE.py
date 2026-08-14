#!/usr/bin/env python3
"""
Prenatal Treatable Genes - Carrier Frequency and NNS Analysis

Analyzes carrier frequencies for 52 genes with prenatal treatable conditions (GATA1 excluded — somatic variant)
and calculates combined panel-level NNS with bootstrap CI.

Prenatal treatable = conditions with in utero interventions that improve outcomes.

Gene list matches Table S13 (52 genes; SMN1 and GATA1 excluded — GATA1 is a somatic variant not appropriate for germline carrier screening).
For per-gene carrier frequency comparison table, see build_prenatal_gene_table_v3.py.

Updated: July 2026
"""

import pandas as pd
import numpy as np
from scipy import stats
import argparse

# 53 Prenatal Treatable Genes (matches Table S13 gene list; SMN1 excluded)
PRENATAL_TREATABLE_GENES = [
    # Connective tissue
    "COL1A1", "COL1A2",

    # Lysosomal storage / metabolic
    "GAA", "LIPA", "IDUA", "IDS", "GALNS", "ARSB", "GUSB", "GBA1",

    # Ectodermal / X-linked
    "EDA",

    # Cardiac channelopathies
    "KCNH2", "SCN5A", "CACNA1C",

    # Thyroid
    "TPO",

    # X-linked endocrine / transporter
    "SLC16A2",

    # Diamond-Blackfan anaemia / ribosomal
    "RPS19", "RPL5", "RPS10", "RPS24",

    # Coagulation — fibrinogen
    "FGB", "FGA", "FGG",

    # Coagulation — factors
    "F2", "F5", "F7", "F8", "F9", "F10", "F11", "F12", "F13A1", "F13B",

    # Haemoglobin
    "HBA1", "HBA2",

    # X-linked immunodeficiency
    "IL2RG",

    # Cobalamin / organic acid metabolism
    "MMACHC", "HLCS", "MMAA", "MMAB", "MMADHC",

    # Rare metabolic
    "NANS", "PHGDH", "DHCR7", "TRMU", "ALDH7A1",

    # Epilepsy / channelopathy
    "SCN8A",

    # TSC
    "TSC1", "TSC2",

    # X-linked / other
    "MAGED2", "ACE",

    # Cystic fibrosis
    "CFTR",
]

def parse_args():
    """Parse command-line arguments"""
    parser = argparse.ArgumentParser(
        description='Analyze prenatal treatable genes'
    )
    parser.add_argument(
        '--aou-plp',
        default='AoU_PLP_GENE_LEVEL.csv',
        help='All of Us P/LP gene-level file'
    )
    parser.add_argument(
        '--gnomad-plp',
        default='gnomAD_PLP_GENE_LEVEL.csv',
        help='gnomAD P/LP gene-level file'
    )
    parser.add_argument(
        '--output',
        default='Prenatal_Treatable_Analysis.csv',
        help='Output file'
    )
    parser.add_argument(
        '--n-bootstrap',
        type=int,
        default=10000,
        help='Number of bootstrap iterations (default: 10000)'
    )
    
    return parser.parse_args()

def calculate_combined_cf_with_ci(carrier_freqs, ci_lowers, ci_uppers, n_bootstrap=10000):
    """
    Calculate combined carrier frequency with bootstrap CI
    
    Parameters:
    -----------
    carrier_freqs : array
        Individual gene carrier frequencies
    ci_lowers : array
        Lower 95% CI bounds
    ci_uppers : array
        Upper 95% CI bounds
    n_bootstrap : int
        Number of bootstrap iterations
    
    Returns:
    --------
    dict with combined_cf, ci_lower, ci_upper, nns, nns_ci_lower, nns_ci_upper
    """
    
    # Combined carrier frequency = 1 - product(1 - individual CFs)
    combined_cf = 1 - np.prod(1 - np.array(carrier_freqs))
    
    # Bootstrap CI
    boot_cfs = []
    for _ in range(n_bootstrap):
        # Sample from CI for each gene
        boot_cf_genes = [
            np.random.uniform(low, high)
            for low, high in zip(ci_lowers, ci_uppers)
        ]
        boot_combined = 1 - np.prod(1 - np.array(boot_cf_genes))
        boot_cfs.append(boot_combined)
    
    ci_lower = np.percentile(boot_cfs, 2.5)
    ci_upper = np.percentile(boot_cfs, 97.5)
    
    # NNS = 1/CF
    nns = 1 / combined_cf if combined_cf > 0 else np.inf
    nns_ci_lower = 1 / ci_upper if ci_upper > 0 else np.inf
    nns_ci_upper = 1 / ci_lower if ci_lower > 0 else np.inf
    
    return {
        'combined_carrier_freq': combined_cf,
        'cf_95ci_lower': ci_lower,
        'cf_95ci_upper': ci_upper,
        'nns': nns,
        'nns_95ci_lower': nns_ci_lower,
        'nns_95ci_upper': nns_ci_upper
    }

def main():
    """Main execution"""
    args = parse_args()
    
    print("="*80)
    print("PRENATAL TREATABLE GENES ANALYSIS")
    print("="*80)
    
    print(f"\nPrenatal treatable genes: {len(PRENATAL_TREATABLE_GENES)}")
    print(f"  (Conditions with in utero interventions)")
    
    # Load gene-level data
    print(f"\nLoading gene-level data...")
    print(f"  All of Us: {args.aou_plp}")
    print(f"  gnomAD: {args.gnomad_plp}")
    
    aou_plp = pd.read_csv(args.aou_plp)
    gnomad_plp = pd.read_csv(args.gnomad_plp)
    
    print(f"✅ Loaded gene-level carrier frequencies")
    
    # Filter to prenatal genes, combined/total ancestry
    aou_prenatal = aou_plp[
        (aou_plp['gene_symbol'].isin(PRENATAL_TREATABLE_GENES)) &
        (aou_plp['ancestry'] == 'ALL')
    ].copy()
    
    gnomad_prenatal = gnomad_plp[
        (gnomad_plp['gene_symbol'].isin(PRENATAL_TREATABLE_GENES)) &
        (gnomad_plp['ancestry'] == 'TOTAL')
    ].copy()
    
    print(f"\nGenes with carrier frequency data:")
    print(f"  All of Us: {len(aou_prenatal)} / {len(PRENATAL_TREATABLE_GENES)} genes")
    print(f"  gnomAD: {len(gnomad_prenatal)} / {len(PRENATAL_TREATABLE_GENES)} genes")
    
    # Identify missing genes
    aou_genes = set(aou_prenatal['gene_symbol'])
    gnomad_genes = set(gnomad_prenatal['gene_symbol'])
    missing_aou = set(PRENATAL_TREATABLE_GENES) - aou_genes
    missing_gnomad = set(PRENATAL_TREATABLE_GENES) - gnomad_genes
    
    if missing_aou:
        print(f"\n⚠️  All of Us missing genes: {', '.join(sorted(missing_aou))}")
    if missing_gnomad:
        print(f"⚠️  gnomAD missing genes: {', '.join(sorted(missing_gnomad))}")
    
    # Calculate combined NNS for All of Us
    print(f"\n{'='*80}")
    print("COMBINED CARRIER FREQUENCY AND NNS")
    print("="*80)
    
    if len(aou_prenatal) > 0:
        aou_result = calculate_combined_cf_with_ci(
            aou_prenatal['carrier_frequency'].values,
            aou_prenatal['ci_95_lower'].values,
            aou_prenatal['ci_95_upper'].values,
            args.n_bootstrap
        )
        
        print(f"\n📊 All of Us ({len(aou_prenatal)} genes):")
        print(f"   Combined Carrier Frequency: {aou_result['combined_carrier_freq']*100:.2f}%")
        print(f"   95% CI: {aou_result['cf_95ci_lower']*100:.2f}% - {aou_result['cf_95ci_upper']*100:.2f}%")
        print(f"   NNS: {aou_result['nns']:.1f}")
        print(f"   95% CI: {aou_result['nns_95ci_lower']:.1f} - {aou_result['nns_95ci_upper']:.1f}")
    else:
        aou_result = None
        print(f"\n⚠️  No All of Us data available")
    
    # Calculate combined NNS for gnomAD
    if len(gnomad_prenatal) > 0:
        gnomad_result = calculate_combined_cf_with_ci(
            gnomad_prenatal['carrier_frequency'].values,
            gnomad_prenatal['ci_95_lower'].values,
            gnomad_prenatal['ci_95_upper'].values,
            args.n_bootstrap
        )
        
        print(f"\n📊 gnomAD ({len(gnomad_prenatal)} genes):")
        print(f"   Combined Carrier Frequency: {gnomad_result['combined_carrier_freq']*100:.2f}%")
        print(f"   95% CI: {gnomad_result['cf_95ci_lower']*100:.2f}% - {gnomad_result['cf_95ci_upper']*100:.2f}%")
        print(f"   NNS: {gnomad_result['nns']:.1f}")
        print(f"   95% CI: {gnomad_result['nns_95ci_lower']:.1f} - {gnomad_result['nns_95ci_upper']:.1f}")
    else:
        gnomad_result = None
        print(f"\n⚠️  No gnomAD data available")
    
    # Save results
    print(f"\n💾 Saving results to: {args.output}")
    
    results = []
    
    if aou_result:
        results.append({
            'dataset': 'All of Us',
            'n_prenatal_genes_total': len(PRENATAL_TREATABLE_GENES),
            'n_genes_with_data': len(aou_prenatal),
            'combined_carrier_freq': aou_result['combined_carrier_freq'],
            'cf_95ci_lower': aou_result['cf_95ci_lower'],
            'cf_95ci_upper': aou_result['cf_95ci_upper'],
            'nns': aou_result['nns'],
            'nns_95ci_lower': aou_result['nns_95ci_lower'],
            'nns_95ci_upper': aou_result['nns_95ci_upper']
        })
    
    if gnomad_result:
        results.append({
            'dataset': 'gnomAD',
            'n_prenatal_genes_total': len(PRENATAL_TREATABLE_GENES),
            'n_genes_with_data': len(gnomad_prenatal),
            'combined_carrier_freq': gnomad_result['combined_carrier_freq'],
            'cf_95ci_lower': gnomad_result['cf_95ci_lower'],
            'cf_95ci_upper': gnomad_result['cf_95ci_upper'],
            'nns': gnomad_result['nns'],
            'nns_95ci_lower': gnomad_result['nns_95ci_lower'],
            'nns_95ci_upper': gnomad_result['nns_95ci_upper']
        })
    
    results_df = pd.DataFrame(results)
    results_df.to_csv(args.output, index=False)
    
    # Also save gene-level details
    gene_level_output = args.output.replace('.csv', '_gene_level.csv')
    
    gene_level_data = []
    
    for _, row in aou_prenatal.iterrows():
        gene_level_data.append({
            'gene_symbol': row['gene_symbol'],
            'dataset': 'All of Us',
            'carrier_frequency': row['carrier_frequency'],
            'ci_95_lower': row['ci_95_lower'],
            'ci_95_upper': row['ci_95_upper'],
            'n_variants': row.get('n_variants', np.nan)
        })
    
    for _, row in gnomad_prenatal.iterrows():
        gene_level_data.append({
            'gene_symbol': row['gene_symbol'],
            'dataset': 'gnomAD',
            'carrier_frequency': row['carrier_frequency'],
            'ci_95_lower': row['ci_95_lower'],
            'ci_95_upper': row['ci_95_upper'],
            'n_variants': row.get('n_variants', np.nan)
        })
    
    gene_level_df = pd.DataFrame(gene_level_data)
    gene_level_df = gene_level_df.sort_values(['gene_symbol', 'dataset'])
    gene_level_df.to_csv(gene_level_output, index=False)
    
    print(f"✅ Saved summary: {args.output}")
    print(f"✅ Saved gene-level: {gene_level_output}")
    
    # Top genes by carrier frequency
    if len(gene_level_df) > 0:
        print(f"\n📊 Top 10 prenatal treatable genes by carrier frequency:")
        print(f"{'Gene':<12} {'Dataset':<12} {'CF':<12} {'1 in N':<10}")
        print("-"*50)
        
        top_genes = gene_level_df.nlargest(10, 'carrier_frequency')
        for _, row in top_genes.iterrows():
            cf_pct = row['carrier_frequency'] * 100
            nns = 1 / row['carrier_frequency'] if row['carrier_frequency'] > 0 else np.inf
            print(f"{row['gene_symbol']:<12} {row['dataset']:<12} {cf_pct:>6.3f}%      1 in {nns:>5.0f}")
    
    print("\n" + "="*80)
    print("✅ ANALYSIS COMPLETE!")
    print("="*80)

if __name__ == '__main__':
    main()

# Usage:
# python prenatal_treatable_analysis_STANDALONE.py \
#     --aou-plp AoU_PLP_GENE_LEVEL.csv \
#     --gnomad-plp gnomAD_PLP_GENE_LEVEL.csv \
#     --output Prenatal_Treatable_NNS.csv \
#     --n-bootstrap 10000
