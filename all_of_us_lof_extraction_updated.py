"""
All of Us - LoF Variant Extraction using Hail
==============================================

This script identifies Loss of Function (LoF) variants from the VAT dataset
using Hail to filter by consequence terms, then saves them for downstream processing.

LoF consequences include:
- transcript_ablation
- splice_acceptor_variant
- splice_donor_variant
- stop_gained
- frameshift_variant
- stop_lost
- start_lost

Gene filtering:
- Excludes SMN1, GATA1, MT-RNR1 (excluded from all variant types)
- Excludes 15 autosomal dominant gain-of-function genes:
  GLUD1, FAM111A, CYP11B1, CLCN2, KCNJ5, KCNJ11, STAT3, TSHR,
  CACNA1D, GUCY2C, NLRC4, NLRP3, RAC2, CHRNB1, CACNA1C

After running this script, use the main extraction script to process LoF variants
through canonical transcript processing, X chromosome identification, etc.

Author: LoF extraction component
Date: 2026-02-13
"""

from datetime import datetime
import os
import pandas as pd
import hail as hl

################################################################################
# SETUP
################################################################################

bucket = os.getenv('WORKSPACE_BUCKET')
genomic_location = os.getenv("CDR_STORAGE_PATH")

start_time = datetime.now()

print("=" * 80)
print("LOF VARIANT EXTRACTION USING HAIL")
print("=" * 80)
print(f"✅ Start time: {start_time}")
print(f"✅ Workspace bucket: {bucket}")
print("=" * 80)

# Gene exclusion lists
GENES_TO_EXCLUDE_ALL = ['SMN1', 'GATA1', 'MT-RNR1']  # Exclude from ALL variant types

GENES_TO_EXCLUDE_LOF = [  # Exclude from LoF (autosomal dominant gain-of-function)
    'GLUD1', 'FAM111A', 'CYP11B1', 'CLCN2', 'KCNJ5', 'KCNJ11', 'STAT3', 'TSHR',
    'CACNA1D', 'GUCY2C', 'NLRC4', 'NLRP3', 'RAC2', 'CHRNB1', 'CACNA1C'
]

print(f"\n📋 Gene filtering configuration:")
print(f"   Genes excluded from all variant types: {', '.join(GENES_TO_EXCLUDE_ALL)}")
print(f"   Additional LoF exclusions (gain-of-function): {len(GENES_TO_EXCLUDE_LOF)} genes")

################################################################################
# INITIALIZE HAIL
################################################################################

print(f"\n🔧 Initializing Hail...")
hl.default_reference(new_default_reference = "GRCh38")
print(f"✅ Hail initialized with GRCh38 reference")

################################################################################
# LOAD VAT DATASET
################################################################################

vat_path = "gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/aux/vat/vat_complete.bgz.tsv.gz"

print("\n" + "=" * 80)
print(f"📊 Loading VAT dataset")
print("⏰ Expected time: 1-2 hours")
print("=" * 80)

vat_load_start = datetime.now()

vat_ht = hl.import_table(
    vat_path,
    delimiter="\t",
    force_bgz=True,
    source_file_field="source",
    types={
        "vid": hl.tstr,
        "gene_symbol": hl.tstr,
        "consequence": hl.tstr,
        "clinvar_classification": hl.tstr
    },
    min_partitions=128
)

print(f"✅ VAT loaded in {datetime.now() - vat_load_start}")
print(f"   Partitions: {vat_ht.n_partitions()}")

################################################################################
# LOAD GENE LIST FROM EXCEL
################################################################################

print("\n" + "=" * 80)
print("📋 Loading gene list from Excel")
print("=" * 80)

# Load gene list from Excel file
gene_file = "/home/jupyter/mmc2_unmerged.xlsx"

try:
    # Read Excel file (no header, genes in column 2)
    gene_df = pd.read_excel(gene_file, header=None)
    
    print(f"✅ Loaded Excel file: {gene_file}")
    print(f"   Shape: {gene_df.shape}")
    
    # Genes are in column 2 (index 2), clean whitespace/newlines
    gene_list = gene_df.iloc[:, 2].dropna().astype(str).str.strip().unique().tolist()
    
    print(f"\n✅ Extracted {len(gene_list)} unique genes from column 2")
    print(f"\n📊 First 10 genes:")
    for gene in gene_list[:10]:
        print(f"   - {gene}")
    
    # Filter out excluded genes
    print(f"\n🔧 Filtering excluded genes...")
    n_before = len(gene_list)
    
    # Remove genes excluded from all variant types
    all_exclusions = set(GENES_TO_EXCLUDE_ALL + GENES_TO_EXCLUDE_LOF)
    gene_list = [g for g in gene_list if g not in all_exclusions]
    
    n_excluded = n_before - len(gene_list)
    print(f"   Original gene count: {n_before}")
    print(f"   Excluded genes: {n_excluded}")
    print(f"   Final gene count: {len(gene_list)}")
    
    if n_excluded > 0:
        excluded_found = [g for g in all_exclusions if g in gene_df.iloc[:, 2].astype(str).str.strip().values]
        if excluded_found:
            print(f"   Genes removed: {', '.join(excluded_found)}")
    
except Exception as e:
    print(f"❌ Error loading Excel file: {e}")
    print(f"\n💡 Troubleshooting tips:")
    print(f"   1. Make sure file is uploaded to /home/jupyter/")
    print(f"   2. Check filename is exactly: mmc2_unmerged.xlsx")
    print(f"   3. Verify openpyxl is installed: pip install openpyxl --break-system-packages")
    raise

################################################################################
# DEFINE LOF CONSEQUENCES
################################################################################

print("\n" + "=" * 80)
print("🧬 Defining LoF consequence terms")
print("=" * 80)

lof_consequences = [
    'transcript_ablation',
    'splice_acceptor_variant',
    'splice_donor_variant',
    'stop_gained',
    'frameshift_variant',
    'stop_lost',
    'start_lost'
]

print("LoF consequences:")
for consequence in lof_consequences:
    print(f"   - {consequence}")

################################################################################
# FILTER VAT FOR LOF VARIANTS
################################################################################

print("\n" + "=" * 80)
print("🔍 Filtering VAT for LoF variants")
print("=" * 80)

filter_start = datetime.now()

# Filter by gene
gene_filter = hl.literal(set(gene_list)).contains(vat_ht.gene_symbol)

# Filter by LoF consequence
consequence_filter = hl.literal(set(lof_consequences)).contains(vat_ht.consequence)

# Filter by ClinVar classification (empty or not provided - we want non-P/LP variants)
clinvar_filter = (
    hl.is_missing(vat_ht.clinvar_classification) |
    (vat_ht.clinvar_classification == '') |
    (vat_ht.clinvar_classification == 'not_provided')
)

# Apply all filters
lof_filtered = vat_ht.filter(gene_filter & consequence_filter & clinvar_filter)

print(f"✅ Filters applied in {datetime.now() - filter_start}")

################################################################################
# EXPORT RESULTS
################################################################################

print("\n" + "=" * 80)
print("💾 Exporting LoF variants")
print("=" * 80)

export_start = datetime.now()

output_path = f"{bucket}/matched_variants_lof.tsv"
lof_filtered.export(output_path)

print(f"✅ Exported in {datetime.now() - export_start}")
print(f"   Output: {output_path}")

################################################################################
# DOWNLOAD AND PROCESS
################################################################################

print("\n📥 Downloading results...")
os.system(f"gsutil cp {output_path} .")

# Load and inspect
lof_df = pd.read_csv('matched_variants_lof.tsv', sep='\t')

print(f"\n📊 LoF Variants Summary:")
print(f"   Total records: {len(lof_df):,} (all transcripts)")
print(f"   Unique variants (VID): {lof_df['vid'].nunique():,}")
print(f"   Genes: {lof_df['gene_symbol'].nunique()}")

# Top genes by variant count
print(f"\n🧬 Top 20 genes by LoF variant count:")
gene_counts = lof_df['gene_symbol'].value_counts().head(20)
print(gene_counts)

# Consequence distribution
print(f"\n📊 LoF consequence distribution:")
consequence_counts = lof_df['consequence'].value_counts()
print(consequence_counts)

# Save as CSV
lof_df.to_csv('matched_variants_lof.csv', index=False)
print(f"\n✅ Saved: matched_variants_lof.csv")

# Create file for next step with Chromosome, Position, Ref, Alt, VID
# Extract chromosome, position, ref, alt from VID
lof_df['Chromosome'] = lof_df['vid'].str.split('-').str[0]
lof_df['Position'] = lof_df['vid'].str.split('-').str[1]
lof_df['Ref'] = lof_df['vid'].str.split('-').str[2]
lof_df['Alt'] = lof_df['vid'].str.split('-').str[3]

# Save file with standard columns for main extraction script
lof_for_extraction = lof_df[['vid', 'Chromosome', 'Position', 'Ref', 'Alt', 'gene_symbol', 'consequence']].drop_duplicates(subset='vid')
lof_for_extraction.rename(columns={'vid': 'VID', 'gene_symbol': 'Gene'}, inplace=True)
lof_for_extraction.to_csv('LoF_with_VID.txt', sep='\t', index=False)

print(f"✅ Saved: LoF_with_VID.txt ({len(lof_for_extraction):,} unique variants)")
print(f"\n📋 File format for main extraction script:")
print(lof_for_extraction.head(10))

################################################################################
# FINAL SUMMARY
################################################################################

total_time = datetime.now() - start_time

print("\n" + "=" * 80)
print("✅ LOF EXTRACTION COMPLETE!")
print("=" * 80)
print(f"\n⏱️  Total time: {total_time}")

print(f"\n📝 NEXT STEPS:")
print(f"   1. Review matched_variants_lof.csv")
print(f"   2. Use LoF_with_VID.txt as input to the main extraction script")
print(f"   3. In main script, uncomment the LoF section:")
print(f"      lof_files = process_variant_list(")
print(f"          variant_file='/home/jupyter/LoF_with_VID.txt',")
print(f"          variant_type='LOF',")
print(f"          file_separator='\\t'")
print(f"      )")
print(f"   4. This will create:")
print(f"      - LOF_CANONICAL_TRANSCRIPT_WITH_STARRATING.csv")
print(f"      - LOF_AUTOSOMAL_VARIANT_LEVEL_EXACT_COUNTS.csv")
print(f"      - LOF_X_CHROMOSOME_VARIANTS.csv")
print(f"      - LOF_X_CHROMOSOME_FEMALE_EXACT_COUNTS.csv")

print("\n" + "=" * 80)
