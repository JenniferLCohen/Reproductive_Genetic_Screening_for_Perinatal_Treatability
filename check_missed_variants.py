"""
Check missed variants from GENEINFO secondary gene parsing
Filters to >=2 star rating and P/LP classification
to determine true impact on carrier frequency analysis
"""

import gzip
import pandas as pd

CLINVAR_VCF = '/Users/jc745/clinvar.vcf.gz'
GENE_CLASS_CSV = '../Figures/v2_XL_corrected/GENE_CLASSIFICATIONS_293.csv'
MIN_STARS = 2

def parse_star_rating(revstat):
    star_map = {
        'practice_guideline': 4,
        'reviewed_by_expert_panel': 3,
        'criteria_provided,_multiple_submitters,_no_conflicts': 2,
        'criteria_provided,_single_submitter': 1,
        'no_assertion_criteria_provided': 0,
        'no_assertion_provided': 0
    }
    max_stars = 0
    for status, rating in star_map.items():
        if status in revstat.lower():
            max_stars = max(max_stars, rating)
    return max_stars

def is_plp(clnsig):
    if 'Benign' in clnsig or 'Conflicting' in clnsig:
        return False
    return 'Pathogenic' in clnsig or 'Likely_pathogenic' in clnsig

def is_vus(clnsig):
    return 'Uncertain_significance' in clnsig

genes_293 = set(pd.read_csv(GENE_CLASS_CSV)['gene_symbol'])

missed_plp = []
missed_vus = []
missed_low_star = 0

with gzip.open(CLINVAR_VCF, 'rt') as f:
    for line in f:
        if line.startswith('#'): continue
        if 'GENEINFO=' not in line: continue

        fields = line.strip().split('\t')
        info = {k: v for item in fields[7].split(';')
                if '=' in item for k, v in [item.split('=', 1)]}

        if 'GENEINFO' not in info: continue

        all_genes = [g.split(':')[0] for g in info['GENEINFO'].split('|')]
        primary = all_genes[0]
        secondary = all_genes[1:]

        # Only care about records where panel gene is ONLY secondary
        if primary in genes_293: continue
        panel_secondary = [g for g in secondary if g in genes_293]
        if not panel_secondary: continue

        gene = panel_secondary[0]

        # Star rating filter
        revstat = info.get('CLNREVSTAT', '')
        stars = parse_star_rating(revstat)
        if stars < MIN_STARS:
            missed_low_star += 1
            continue

        clnsig = info.get('CLNSIG', '')
        chrom = fields[0]
        pos = fields[1]
        ref = fields[3]
        alt = fields[4]

        record = {
            'Gene': gene,
            'Chromosome': chrom,
            'Position': pos,
            'Ref': ref,
            'Alt': alt,
            'ClinicalSignificance': clnsig,
            'StarRating': stars,
            'Primary_gene_in_record': primary,
            'All_genes_in_record': '|'.join(all_genes)
        }

        if is_plp(clnsig):
            missed_plp.append(record)
        elif is_vus(clnsig):
            missed_vus.append(record)

print("=" * 60)
print("MISSED VARIANT ANALYSIS (secondary GENEINFO only)")
print("=" * 60)
print(f"\nTotal missed (all star ratings): 1,661")
print(f"Filtered out by <{MIN_STARS} stars:  {missed_low_star:,}")
print(f"Remaining ≥{MIN_STARS} stars:         {len(missed_plp) + len(missed_vus):,}")
print(f"  P/LP (would affect CF):   {len(missed_plp):,}")
print(f"  VUS (would affect VUS CF): {len(missed_vus):,}")

if missed_plp:
    plp_df = pd.DataFrame(missed_plp)
    print(f"\nP/LP missed variants by gene:")
    print(plp_df.groupby('Gene').size().sort_values(ascending=False).to_string())
    plp_df.to_csv('missed_plp_variants.csv', index=False)
    print(f"\nSaved: missed_plp_variants.csv")

if missed_vus:
    vus_df = pd.DataFrame(missed_vus)
    print(f"\nVUS missed variants by gene:")
    print(vus_df.groupby('Gene').size().sort_values(ascending=False).to_string())
    vus_df.to_csv('missed_vus_variants.csv', index=False)
    print(f"\nSaved: missed_vus_variants.csv")
