#!/usr/bin/env python3
"""
Extract ClinVar P/LP and VUS variants for 293-gene carrier screening panel

This script extracts variants by gene name from ClinVar VCF using the GENEINFO field.
Filters to ≥2 star rating and creates separate P/LP and VUS files.

Extracted from: Reproducible_Code_Appendix.Rmd
Updated: March 20, 2026
"""

import pandas as pd
import gzip
import argparse
from pathlib import Path

def parse_args():
    """Parse command-line arguments"""
    parser = argparse.ArgumentParser(
        description='Extract ClinVar variants by gene name'
    )
    parser.add_argument(
        '--clinvar-vcf',
        default='clinvar.vcf.gz',
        help='Path to ClinVar VCF file (default: clinvar.vcf.gz)'
    )
    parser.add_argument(
        '--gene-list',
        required=True,
        help='Gene list file (CSV/Excel with gene_symbol column)'
    )
    parser.add_argument(
        '--min-stars',
        type=int,
        default=2,
        help='Minimum star rating (default: 2)'
    )
    parser.add_argument(
        '--output-plp',
        default='ClinVar_PLP_2star_for_analysis.csv',
        help='Output file for P/LP variants'
    )
    parser.add_argument(
        '--output-vus',
        default='ClinVar_VUS_2star_for_analysis.csv',
        help='Output file for VUS variants'
    )
    
    return parser.parse_args()

def load_gene_list(gene_file):
    """Load gene list from CSV or Excel file"""
    if gene_file.endswith('.xlsx'):
        gene_df = pd.read_excel(gene_file)
    else:
        gene_df = pd.read_csv(gene_file)
    
    if 'gene_symbol' in gene_df.columns:
        genes = set(gene_df['gene_symbol'].tolist())
    elif 'Gene' in gene_df.columns:
        genes = set(gene_df['Gene'].tolist())
    else:
        raise ValueError("Gene list must have 'gene_symbol' or 'Gene' column")
    
    return genes

def parse_star_rating(revstat):
    """Convert review status to star rating"""
    star_map = {
        'practice_guideline': 4,
        'reviewed_by_expert_panel': 3,
        'criteria_provided,_multiple_submitters,_no_conflicts': 2,
        'criteria_provided,_single_submitter': 1,
        'no_assertion_criteria_provided': 0,
        'no_assertion_provided': 0
    }
    
    # Handle multiple review statuses (take maximum)
    max_stars = 0
    for status, rating in star_map.items():
        if status in revstat.lower():
            max_stars = max(max_stars, rating)
    
    return max_stars

def extract_clinvar_variants(clinvar_vcf, genes_of_interest, min_stars=2):
    """
    Extract P/LP and VUS variants from ClinVar VCF
    
    Parameters:
    -----------
    clinvar_vcf : str
        Path to ClinVar VCF file
    genes_of_interest : set
        Set of gene symbols to extract
    min_stars : int
        Minimum star rating threshold
    
    Returns:
    --------
    tuple: (plp_variants, vus_variants) as DataFrames
    """
    
    print(f"Processing {clinvar_vcf}...")
    print(f"Genes of interest: {len(genes_of_interest)}")
    print(f"Minimum star rating: {min_stars}")
    
    clinvar_plp = []
    clinvar_vus = []
    
    total_lines = 0
    variants_in_genes = 0
    
    with gzip.open(clinvar_vcf, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            total_lines += 1
            if total_lines % 100000 == 0:
                print(f"  Processed {total_lines:,} lines... "
                      f"(P/LP: {len(clinvar_plp)}, VUS: {len(clinvar_vus)})")
            
            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = fields[1]
            ref = fields[3]
            alt = fields[4]
            info = fields[7]
            
            # Parse INFO field
            info_dict = {}
            for item in info.split(';'):
                if '=' in item:
                    key, value = item.split('=', 1)
                    info_dict[key] = value
            
            # Extract gene name from GENEINFO
            if 'GENEINFO' not in info_dict:
                continue
            
            gene_info = info_dict['GENEINFO']
            # GENEINFO format: GENE:GENEID or GENE1:ID1|GENE2:ID2
            gene = gene_info.split(':')[0].split('|')[0]
            
            # Check if gene is in our list
            if gene not in genes_of_interest:
                continue
            
            variants_in_genes += 1
            
            # Extract clinical significance
            if 'CLNSIG' not in info_dict:
                continue
            
            clnsig = info_dict['CLNSIG']
            
            # Extract and filter star rating
            if 'CLNREVSTAT' in info_dict:
                revstat = info_dict['CLNREVSTAT']
                stars = parse_star_rating(revstat)
            else:
                stars = 0
            
            if stars < min_stars:
                continue
            
            # Create variant record
            variant_record = {
                'Chromosome': chrom,
                'Position': pos,
                'Ref': ref,
                'Alt': alt,
                'Gene': gene,
                'ClinicalSignificance': clnsig,
                'ReviewStatus': info_dict.get('CLNREVSTAT', ''),
                'StarRating': stars,
                'ClinVar_VariationID': info_dict.get('ALLELEID', ''),
                'ClinVar_Condition': info_dict.get('CLNDN', ''),
                'NC_Variant_ID': f"{chrom}_{pos}_{ref}_{alt}"
            }
            
            # Classify as P/LP or VUS
            if 'Pathogenic' in clnsig and 'Benign' not in clnsig and 'Conflicting' not in clnsig:
                clinvar_plp.append(variant_record)
            elif 'Likely_pathogenic' in clnsig and 'Benign' not in clnsig and 'Conflicting' not in clnsig:
                clinvar_plp.append(variant_record)
            elif 'Uncertain_significance' in clnsig:
                clinvar_vus.append(variant_record)
    
    print(f"\n✅ Processing complete!")
    print(f"   Total ClinVar variants: {total_lines:,}")
    print(f"   Variants in panel genes: {variants_in_genes:,}")
    print(f"   P/LP variants (≥{min_stars} stars): {len(clinvar_plp):,}")
    print(f"   VUS variants (≥{min_stars} stars): {len(clinvar_vus):,}")
    
    # Create DataFrames
    plp_df = pd.DataFrame(clinvar_plp)
    vus_df = pd.DataFrame(clinvar_vus)
    
    return plp_df, vus_df

def main():
    """Main execution"""
    args = parse_args()
    
    print("="*80)
    print("CLINVAR VARIANT EXTRACTION BY GENE NAME")
    print("="*80)
    
    # Load gene list
    print(f"\nLoading gene list from: {args.gene_list}")
    genes_of_interest = load_gene_list(args.gene_list)
    print(f"✅ Loaded {len(genes_of_interest)} genes")
    
    # Extract variants
    plp_df, vus_df = extract_clinvar_variants(
        args.clinvar_vcf,
        genes_of_interest,
        args.min_stars
    )
    
    # Save results
    print(f"\n💾 Saving results...")
    
    if len(plp_df) > 0:
        plp_df.to_csv(args.output_plp, index=False)
        print(f"   ✅ P/LP variants: {args.output_plp}")
        print(f"      Genes with P/LP: {plp_df['Gene'].nunique()}")
    else:
        print(f"   ⚠️  No P/LP variants found")
    
    if len(vus_df) > 0:
        vus_df.to_csv(args.output_vus, index=False)
        print(f"   ✅ VUS variants: {args.output_vus}")
        print(f"      Genes with VUS: {vus_df['Gene'].nunique()}")
    else:
        print(f"   ⚠️  No VUS variants found")
    
    print("\n" + "="*80)
    print("✅ EXTRACTION COMPLETE!")
    print("="*80)

if __name__ == '__main__':
    main()

# Usage:
# python extract_clinvar_variants_by_gene_STANDALONE.py \
#     --gene-list mmc2_unmerged.xlsx \
#     --clinvar-vcf clinvar.vcf.gz \
#     --min-stars 2 \
#     --output-plp ClinVar_PLP_2star.csv \
#     --output-vus ClinVar_VUS_2star.csv
