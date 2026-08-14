"""
Build per-gene comparison table for 52 prenatal treatable genes (GATA1 excluded — somatic variant):
Gene, gnomAD CF (%), AoU CF (%), Difference (%), AoU/gnomAD CF Ratio

Uses corrected (HWE-fixed) gene-level CSVs:
  AoU_PLP_GENE_LEVEL.csv
  gnomAD_PLP_GENE_LEVEL.csv
"""

import pandas as pd
import argparse

PRENATAL_TREATABLE_GENES = [
    "COL1A1", "COL1A2",
    "GAA",
    "LIPA",
    "IDUA",
    "IDS",
    "GALNS",
    "ARSB",
    "GUSB",
    "GBA1",
    "EDA",
    "KCNH2", "SCN5A", "CACNA1C",
    "TPO",
    "SLC16A2",
    "RPS19", "RPL5", "RPS10", "RPS24",
    "FGB", "FGA", "FGG",
    "F2", "F5", "F7", "F8", "F9", "F10", "F11", "F12", "F13A1", "F13B",
    "HBA1", "HBA2",
    "IL2RG",
    "MMACHC",
    "HLCS",
    "MMAA", "MMAB", "MMADHC",
    "NANS",
    "PHGDH",
    "DHCR7",
    "TRMU",
    "ALDH7A1",
    "SCN8A",
    "TSC1", "TSC2",
    "MAGED2",
    "ACE",
    "CFTR"
]

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--aou-plp', default='AoU_PLP_GENE_LEVEL.csv')
    p.add_argument('--gnomad-plp', default='gnomAD_PLP_GENE_LEVEL.csv')
    p.add_argument('--output', default='Prenatal_Treatable_Gene_Level_Comparison.csv')
    return p.parse_args()

def main():
    args = parse_args()

    aou = pd.read_csv(args.aou_plp)
    gnomad = pd.read_csv(args.gnomad_plp)

    # AoU uses 'ALL', gnomAD uses 'TOTAL' for the overall ancestry row
    aou_all = aou[aou['ancestry'] == 'ALL'].set_index('gene_symbol')
    gnomad_total = gnomad[gnomad['ancestry'] == 'TOTAL'].set_index('gene_symbol')

    rows = []
    missing_aou = []
    missing_gnomad = []

    for gene in PRENATAL_TREATABLE_GENES:
        aou_cf = aou_all.loc[gene, 'carrier_frequency'] if gene in aou_all.index else None
        gnomad_cf = gnomad_total.loc[gene, 'carrier_frequency'] if gene in gnomad_total.index else None

        if aou_cf is None:
            missing_aou.append(gene)
        if gnomad_cf is None:
            missing_gnomad.append(gene)

        aou_pct = aou_cf * 100 if aou_cf is not None else None
        gnomad_pct = gnomad_cf * 100 if gnomad_cf is not None else None

        if aou_pct is not None and gnomad_pct is not None:
            difference = aou_pct - gnomad_pct
            ratio = aou_pct / gnomad_pct if gnomad_pct != 0 else None
        else:
            difference = None
            ratio = None

        rows.append({
            'Gene': gene,
            'gnomAD Carrier Frequency (%)': round(gnomad_pct, 4) if gnomad_pct is not None else 'no qualifying variants',
            'AoU Carrier Frequency (%)': round(aou_pct, 4) if aou_pct is not None else 'no qualifying variants',
            'Difference (%)': round(difference, 4) if difference is not None else 'N/A',
            'AoU/gnomAD CF Ratio': round(ratio, 3) if ratio is not None else 'N/A'
        })

    out = pd.DataFrame(rows)
    out.to_csv(args.output, index=False)

    print(f"Saved: {args.output}")
    print(f"\nTotal genes: {len(PRENATAL_TREATABLE_GENES)}")
    print(f"Genes missing from AoU: {len(missing_aou)} -> {missing_aou}")
    print(f"Genes missing from gnomAD: {len(missing_gnomad)} -> {missing_gnomad}")
    print(f"\n{out.to_string(index=False)}")

if __name__ == '__main__':
    main()
