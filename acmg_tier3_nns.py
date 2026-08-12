"""
ACMG Tier 3 — NNS Analysis (P/LP, all ancestries combined)

Reads pre-computed gene-level carrier frequency CSVs from the main pipeline
and calculates combined carrier frequency and NNS for the subset of ACMG tier 3
genes that overlap with the 293-gene carrier screening panel (30 of 112 genes).

Note: 82 of 112 ACMG tier 3 genes are not in the 293-gene panel and therefore
have no carrier frequency data. Results reflect the 30 overlapping genes only
and should be interpreted as a partial estimate of ACMG tier 3 panel yield.

Usage:
  python acmg_tier3_nns.py \
      --aou-plp AoU_PLP_GENE_LEVEL.csv \
      --gnomad-plp gnomAD_PLP_GENE_LEVEL.csv

Output:
  ACMG_Tier3_NNS.csv              — panel-level summary (30 overlapping genes)
  ACMG_Tier3_Gene_Level.csv       — per-gene carrier frequencies
"""

import argparse
import pandas as pd
import numpy as np

# ── ACMG Tier 3 gene list — 113 genes minus SMN1 = 112 ────────────────────
ACMG_TIER3 = [
    "ABCA3", "ABCC8", "ABCD1", "ACADM", "ACADVL", "ACAT1", "AFF2",
    "AGA", "AGXT", "AHI1", "AIRE", "ALDOB", "ALPL", "ANO10", "ARSA",
    "ARX", "ASL", "ASPA", "ATP7B", "BBS1", "BBS2", "BCKDHB", "BLM",
    "BTD", "CBS", "CC2D2A", "CCDC88C", "CEP290", "CFTR", "CHRNE",
    "CLCN1", "CLRN1", "CNGB3", "COL7A1", "CPT2", "CYP11A1", "CYP21A2",
    "CYP27A1", "CYP27B1", "DHCR7", "DHDDS", "DLD", "DMD", "DYNC2H1",
    "ELP1", "ERCC2", "EVC2", "F8", "F9", "FAH", "FANCC", "FKRP",
    "FKTN", "FMO3", "FMR1", "FXN", "G6PC", "GAA", "GALT",
    "GBA1",    # listed as GBA in ACMG; standardised to GBA1 in this pipeline
    "GBE1", "GJB2", "GLA", "GNPTAB", "GRIP1", "HBA1", "HBA2", "HBB",
    "HEXA", "HPS1", "HPS3", "IDUA", "L1CAM", "LRP2", "MCCC2", "MCOLN1",
    "MCPH1", "MID1", "MLC1", "MMACHC",
    "MMUT",    # MUT in some databases; standardised to MMUT in this pipeline
    "MVK", "NAGA", "NEB", "NPHS1", "NR0B1", "OCA2", "OTC", "PAH",
    "PCDH15", "PKHD1", "PLP1", "PMM2", "POLG", "PRF1", "RARS2",
    "RNASEH2B", "RPGR", "RS1", "SCO2", "SLC19A3", "SLC26A2", "SLC26A4",
    "SLC37A4", "SLC6A8", "SMPD1", "TF", "TMEM216", "TNXB", "TYR",
    "USH2A", "XPC",
    # SMN1 excluded (copy-number based screening; incompatible with SNV pipeline)
]

assert len(ACMG_TIER3) == 112, f"Expected 112 genes, got {len(ACMG_TIER3)}"
assert len(set(ACMG_TIER3)) == 112, "Duplicate genes detected"

def parse_args():
    p = argparse.ArgumentParser(description='ACMG Tier 3 NNS analysis')
    p.add_argument('--aou-plp',    default='AoU_PLP_GENE_LEVEL.csv')
    p.add_argument('--gnomad-plp', default='gnomAD_PLP_GENE_LEVEL.csv')
    p.add_argument('--n-bootstrap', type=int, default=10000)
    p.add_argument('--output',     default='ACMG_Tier3_NNS.csv')
    return p.parse_args()

def product_rule(cfs):
    """Combined carrier frequency via product rule."""
    return 1.0 - np.prod(1.0 - np.array(cfs))

def bootstrap_nns_ci(cfs, ci_lowers, ci_uppers, n_boot=10000, seed=42):
    """Bootstrap 95% CI for NNS from individual gene CFs and their CIs."""
    rng = np.random.default_rng(seed)
    ses = (np.array(ci_uppers) - np.array(ci_lowers)) / (2 * 1.96)
    boot_nns = []
    for _ in range(n_boot):
        sampled = rng.normal(cfs, ses)
        sampled = np.clip(sampled, 0, 1)
        cf_combined = product_rule(sampled)
        if cf_combined > 0:
            boot_nns.append(1.0 / cf_combined)
    boot_nns = np.array(boot_nns)
    return np.percentile(boot_nns, 2.5), np.percentile(boot_nns, 97.5)

def analyse(df_raw, ancestry_col, ancestry_val, dataset_label, n_boot):
    df = df_raw[df_raw[ancestry_col] == ancestry_val].copy()
    df = df[df['gene_symbol'].isin(ACMG_TIER3)]

    genes_with_data = df['gene_symbol'].unique().tolist()
    genes_no_data   = [g for g in ACMG_TIER3 if g not in genes_with_data]

    cfs        = df['carrier_frequency'].values
    ci_lo      = df['ci_95_lower'].values
    ci_hi      = df['ci_95_upper'].values

    cf_combined = product_rule(cfs)
    nns         = 1.0 / cf_combined if cf_combined > 0 else float('inf')
    nns_lo, nns_hi = bootstrap_nns_ci(cfs, ci_lo, ci_hi, n_boot)

    print(f"\n{'='*60}")
    print(f"{dataset_label}  (P/LP only, all ancestries combined)")
    print(f"NNS based on {len(genes_with_data)} of 30 overlapping ACMG genes with data")
    print(f"{'='*60}")
    print(f"ACMG genes with carrier frequency data: {len(genes_with_data)} / {len(ACMG_TIER3)}")
    print(f"Combined carrier frequency: {cf_combined*100:.2f}%")
    print(f"NNS: {nns:.1f}  (95% CI {nns_lo:.1f} – {nns_hi:.1f})")
    print(f"1 in N: 1 in {round(nns)}")

    if genes_no_data:
        print(f"\nGenes with no qualifying P/LP variants (≥2 stars):")
        print("  " + ", ".join(sorted(genes_no_data)))

    print(f"\nTop 15 genes by carrier frequency:")
    top = df.sort_values('carrier_frequency', ascending=False).head(15)
    print(top[['gene_symbol','carrier_frequency','ci_95_lower','ci_95_upper']].to_string(index=False))

    return {
        'dataset':            dataset_label,
        'n_acmg_genes':       len(ACMG_TIER3),
        'n_genes_with_data':  len(genes_with_data),
        'n_genes_no_data':    len(genes_no_data),
        'genes_no_data':      "; ".join(sorted(genes_no_data)),
        'combined_cf':        cf_combined,
        'combined_cf_pct':    round(cf_combined * 100, 4),
        'nns':                round(nns, 2),
        'nns_95ci_lower':     round(nns_lo, 2),
        'nns_95ci_upper':     round(nns_hi, 2),
        '1_in_n':             f'1 in {round(nns)}',
    }, df.sort_values('carrier_frequency', ascending=False)

def main():
    args = parse_args()

    print("="*60)
    print("ACMG TIER 3 — NNS ANALYSIS (P/LP, all ancestries combined)")
    print("="*60)
    print(f"ACMG tier 3 gene list: 112 genes (113 minus SMN1)")
    print(f"Genes overlapping with 293-panel: 30 of 112")
    print(f"Note: NNS reflects the 30 overlapping genes only.")
    print(f"      GBA → GBA1 and MUT → MMUT standardisation applied.")

    aou    = pd.read_csv(args.aou_plp)
    gnomad = pd.read_csv(args.gnomad_plp)

    # Determine ancestry column values
    aou_anc    = 'ALL'   if 'ALL'   in aou['ancestry'].values else aou['ancestry'].unique()[0]
    gno_anc    = 'TOTAL' if 'TOTAL' in gnomad['ancestry'].values else gnomad['ancestry'].unique()[0]

    aou_result,    aou_genes    = analyse(aou,    'ancestry', aou_anc,
                                          'All of Us', args.n_bootstrap)
    gnomad_result, gnomad_genes = analyse(gnomad, 'ancestry', gno_anc,
                                          'gnomAD v4.1', args.n_bootstrap)

    # ── Save summary ──────────────────────────────────────────────────────
    summary = pd.DataFrame([aou_result, gnomad_result])
    summary.to_csv(args.output, index=False)
    print(f"\n✅ Saved summary: {args.output}")

    # ── Save gene-level ───────────────────────────────────────────────────
    aou_genes['dataset']    = 'All of Us'
    gnomad_genes['dataset'] = 'gnomAD v4.1'
    gene_out = pd.concat([aou_genes, gnomad_genes], ignore_index=True)
    gene_csv = args.output.replace('.csv', '_gene_level.csv')
    gene_out[['dataset','gene_symbol','carrier_frequency',
              'ci_95_lower','ci_95_upper','chromosome_type']].to_csv(gene_csv, index=False)
    print(f"✅ Saved gene-level: {gene_csv}")

if __name__ == '__main__':
    main()
