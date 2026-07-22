"""
All of Us - Clean Variant Extraction Script
============================================

Extracts 3 variant lists from All of Us VAT with exact allele counts (no suppression):
1. P/LP (Pathogenic/Likely Pathogenic)
2. VUS (Variants of Uncertain Significance - Star Rating >= 2)
3. LoF (Loss of Function)

For each list:
- Creates VID from Chromosome-Position-Ref-Alt
- Canonical transcript processing
- Gene filtering:
  * Excludes SMN1, GATA1, MT-RNR1 from ALL variant types
  * LoF only: Also excludes 15 gain-of-function genes (GLUD1, FAM111A, etc.)
- X chromosome identification
- Female cohort extraction (n=250,071) for X variants
- Exact allele counts and allele numbers (internal use - no small cell suppression)
- VUS only: Filters to AF <= 0.001% after allele frequencies extracted

Final outputs (4-5 files per variant type):
- {type}_CANONICAL_TRANSCRIPT_WITH_STARRATING.csv
- {type}_AUTOSOMAL_VARIANT_LEVEL_EXACT_COUNTS.csv (or _RARE.csv for VUS)
- {type}_X_CHROMOSOME_VARIANTS.csv
- {type}_X_CHROMOSOME_FEMALE_EXACT_COUNTS.csv (or _RARE.csv for VUS)

Gene-level carrier frequency calculations done separately using MEAN(AN):
  carrier_frequency = sum(AC_across_variants) / mean(AN_across_variants)

Author: Clean extraction for reproducible analysis
Date: 2026-02-13
"""

from datetime import datetime
import os
import pandas as pd
import hail as hl

################################################################################
# SETUP AND CONFIGURATION
################################################################################

# Get workspace configuration
bucket = os.getenv('WORKSPACE_BUCKET')
genomic_location = os.getenv("CDR_STORAGE_PATH")
dataset = os.getenv('WORKSPACE_CDR')

start_time = datetime.now()

print("=" * 80)
print("ALL OF US - CLEAN VARIANT EXTRACTION")
print("=" * 80)
print(f"✅ Start time: {start_time}")
print(f"✅ Workspace bucket: {bucket}")
print(f"✅ Genomic location: {genomic_location}")
print(f"✅ Dataset: {dataset}")
print("=" * 80)

# Female cohort configuration (for X chromosome)
FEMALE_COHORT_SIZE = 250071
FEMALE_COHORT_ID = 137526  # Females_XX_with_WGS cohort

# Gene exclusion lists
GENES_TO_EXCLUDE_ALL = ['SMN1', 'GATA1', 'MT-RNR1']  # Exclude from ALL variant types

GENES_TO_EXCLUDE_LOF = [  # Exclude from LoF only (autosomal dominant gain-of-function)
    'GLUD1', 'FAM111A', 'CYP11B1', 'CLCN2', 'KCNJ5', 'KCNJ11', 'STAT3', 'TSHR',
    'CACNA1D', 'GUCY2C', 'NLRC4', 'NLRP3', 'RAC2', 'CHRNB1', 'CACNA1C'
]

# VUS allele frequency threshold
VUS_AF_THRESHOLD = 0.00001  # 0.001% = 0.00001

################################################################################
# INITIALIZE HAIL
################################################################################

hail_start = datetime.now()
print(f"\n🔧 Initializing Hail at {hail_start}...")

hl.default_reference(new_default_reference = "GRCh38")

print(f"✅ Hail initialized")
print(f"⏱️  Time: {datetime.now() - hail_start}")

################################################################################
# LOAD VAT DATASET
################################################################################

vat_path = "gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/aux/vat/vat_complete.bgz.tsv.gz"

vat_load_start = datetime.now()
print("\n" + "=" * 80)
print(f"📊 Loading VAT table at {vat_load_start}")
print("⏰ Expected time: 1-2 hours")
print("💡 This is a good time to grab lunch or coffee!")
print("=" * 80)

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

print("=" * 80)
print("✅ VAT TABLE LOADED SUCCESSFULLY!")
print(f"   Partitions: {vat_ht.n_partitions()}")
print(f"⏱️  Time elapsed: {datetime.now() - vat_load_start}")
print("=" * 80)


################################################################################
# HELPER FUNCTION: PROCESS VARIANT LIST
################################################################################

def process_variant_list(variant_file, variant_type, file_separator='\t'):
    """
    Process a variant list through the complete workflow:
    1. Load file and create VID from Chromosome-Position-Ref-Alt
    2. Extract from VAT by VID
    3. Canonical transcript processing
    4. X chromosome identification
    5. Female cohort extraction for X variants
    
    Args:
        variant_file: Path to input file with Chromosome, Position, Ref, Alt columns
        variant_type: 'PLP', 'VUS', or 'LOF'
        file_separator: Separator for input file (default: tab)
    
    Returns:
        Dictionary with paths to output files
    """
    
    print("\n" + "=" * 80)
    print(f"PROCESSING {variant_type} VARIANTS")
    print("=" * 80)
    
    # Step 1: Load variant list and create VID
    print(f"\n📂 Loading {variant_file}...")
    df = pd.read_csv(variant_file, sep=file_separator)
    print(f"   Loaded {len(df):,} rows")
    
    # Check for required columns
    required_cols = ['Chromosome', 'Position', 'Ref', 'Alt']
    missing_cols = [col for col in required_cols if col not in df.columns]
    
    if missing_cols:
        print(f"❌ Missing required columns: {missing_cols}")
        print(f"   Available columns: {df.columns.tolist()}")
        return None
    
    # Create VID if not present
    if 'VID' not in df.columns:
        print(f"\n🔧 Creating VID from Chromosome-Position-Ref-Alt...")
        df['VID'] = (
            df['Chromosome'].astype(str) + '-' + 
            df['Position'].astype(str) + '-' + 
            df['Ref'].astype(str) + '-' + 
            df['Alt'].astype(str)
        )
        print(f"   ✅ Created VID for {len(df):,} variants")
    
    # Get VID list
    vid_list = df['VID'].dropna().unique().tolist()
    print(f"   Unique VIDs: {len(vid_list):,}")
    
    # Step 2: Filter VAT and extract
    filter_start = datetime.now()
    print(f"\n🔍 Filtering VAT for {variant_type} variants...")
    
    vid_filter_expression = hl.literal(set(vid_list)).contains(vat_ht.vid)
    filtered_ht = vat_ht.filter(vid_filter_expression)
    
    # Export to workspace bucket
    output_path = f"{bucket}/matched_variants_{variant_type.lower()}.tsv"
    filtered_ht.export(output_path)
    
    print(f"✅ Filtered VAT")
    print(f"⏱️  Time: {datetime.now() - filter_start}")
    
    # Step 3: Download results
    print(f"\n📥 Downloading {variant_type} results...")
    local_file = f"matched_variants_{variant_type.lower()}.tsv"
    os.system(f"gsutil cp {output_path} {local_file}")
    
    # Load results
    results = pd.read_csv(local_file, sep='\t')
    print(f"   Retrieved {len(results):,} records (all transcripts)")
    
    # Save as CSV
    csv_file = f"matched_variants_{variant_type.lower()}.csv"
    results.to_csv(csv_file, index=False)
    print(f"✅ Saved: {csv_file}")
    
    # Step 4: Canonical transcript processing
    print(f"\n🧬 Processing canonical transcripts...")
    
    # Deduplicate to canonical (MANE_SELECT or MANE_PLUS_CLINICAL)
    results['is_canonical'] = (
        (results['mane_select'].notna()) | 
        (results['mane_plus_clinical'].notna())
    )
    
    canonical = results[results['is_canonical']].copy()
    
    # If multiple canonical transcripts per variant, take first
    canonical_dedup = canonical.drop_duplicates(subset='vid', keep='first')
    
    print(f"   All transcripts: {len(results):,}")
    print(f"   Canonical transcripts: {len(canonical_dedup):,}")
    
    # Merge with original variant info (for star ratings, gene names, etc.)
    print(f"\n⭐ Merging with original variant data...")
    
    # Merge on VID
    canonical_with_orig = canonical_dedup.merge(
        df,
        left_on='vid',
        right_on='VID',
        how='left',
        suffixes=('', '_orig')
    )
    
    # Prefer original gene symbol if available (corrects VAT synonyms)
    if 'Gene' in df.columns:
        print(f"   Using original Gene column to correct VAT gene symbols")
        canonical_with_orig['gene_symbol_VAT'] = canonical_with_orig['gene_symbol']  # Keep VAT name for reference
        canonical_with_orig['gene_symbol'] = canonical_with_orig['Gene'].fillna(canonical_with_orig['gene_symbol'])
    
    # GENE SYNONYM STANDARDIZATION
    # Some genes have multiple names - standardize to preferred name
    print(f"\n🔧 Standardizing gene synonyms...")
    gene_name_map = {
        'GBA': 'GBA1',      # Glucocerebrosidase
        'MUT': 'MMUT',      # Methylmalonyl-CoA mutase
        'TAZ': 'TAFAZZIN'   # Tafazzin
    }
    
    n_before_map = canonical_with_orig['gene_symbol'].nunique()
    canonical_with_orig['gene_symbol'] = canonical_with_orig['gene_symbol'].replace(gene_name_map)
    n_after_map = canonical_with_orig['gene_symbol'].nunique()
    
    # Report changes
    for old_name, new_name in gene_name_map.items():
        n_changed = (canonical_with_orig['gene_symbol'] == new_name).sum()
        if n_changed > 0:
            print(f"   {old_name} → {new_name}: {n_changed} variants")
    
    print(f"   Genes before mapping: {n_before_map}")
    print(f"   Genes after mapping: {n_after_map}")
    
    # Add numeric star rating if StarRating present
    if 'StarRating' in canonical_with_orig.columns:
        star_map = {'0': 0, '1': 1, '2': 2, '3': 3, '4': 4}
        canonical_with_orig['StarRating_Numeric'] = canonical_with_orig['StarRating'].map(star_map)
        print(f"   Merged with star ratings")
    
    canonical_final = canonical_with_orig
    
    # FILTER 1: Remove excluded genes (applies to all variant types)
    print(f"\n🔧 Filtering excluded genes...")
    n_before_filter = len(canonical_final)
    
    # Remove SMN1, GATA1, MT-RNR1 from all variant types
    canonical_final = canonical_final[~canonical_final['gene_symbol'].isin(GENES_TO_EXCLUDE_ALL)].copy()
    
    n_excluded_all = n_before_filter - len(canonical_final)
    print(f"   Removed {n_excluded_all:,} variants in excluded genes: {', '.join(GENES_TO_EXCLUDE_ALL)}")
    
    # FILTER 2: For LoF only, also remove gain-of-function genes
    if variant_type == 'LOF':
        n_before_lof_filter = len(canonical_final)
        canonical_final = canonical_final[~canonical_final['gene_symbol'].isin(GENES_TO_EXCLUDE_LOF)].copy()
        n_excluded_lof = n_before_lof_filter - len(canonical_final)
        print(f"   LoF-specific: Removed {n_excluded_lof:,} variants in gain-of-function genes")
        print(f"   Excluded genes: {', '.join(GENES_TO_EXCLUDE_LOF[:5])}... ({len(GENES_TO_EXCLUDE_LOF)} total)")
    
    print(f"   Variants after filtering: {len(canonical_final):,}")
    
    # Save canonical file
    canonical_file = f"{variant_type}_CANONICAL_TRANSCRIPT_WITH_STARRATING.csv"
    canonical_final.to_csv(canonical_file, index=False)
    print(f"✅ Saved: {canonical_file}")
    
    # Step 5: Identify X chromosome variants
    print(f"\n🧬 Identifying X chromosome variants...")
    
    x_variants = canonical_final[
        (canonical_final['contig'] == 'X') | 
        (canonical_final['contig'] == 'chrX') |
        (canonical_final['contig'] == 'x') |
        (canonical_final['contig'] == 'chr23')
    ].copy()
    
    print(f"   X chromosome variants: {len(x_variants):,}")
    
    if len(x_variants) > 0:
        x_genes = x_variants['gene_symbol'].value_counts()
        print(f"   X-linked genes: {len(x_genes)}")
        print(f"\n   Top 10 genes:")
        print(x_genes.head(10))
        
        # Save X chromosome variants
        x_file = f"{variant_type}_X_CHROMOSOME_VARIANTS.csv"
        x_variants.to_csv(x_file, index=False)
        print(f"\n✅ Saved: {x_file}")
        
        # Extract X chromosome VIDs for BigQuery
        x_vids = x_variants['vid'].tolist()
        vid_list_str = "', '".join(x_vids)
        
        # Step 6: Query female cohort for X chromosome exact counts
        print(f"\n👥 Querying female cohort for X chromosome variants...")
        print(f"   Female cohort size: {FEMALE_COHORT_SIZE:,}")
        print(f"   X chromosome VIDs to query: {len(x_vids):,}")
        
        query = f"""
        WITH female_ids AS (
          -- Get all female person_ids from cohort
          SELECT person_id
          FROM `{dataset}.cb_search_person`
          WHERE person_id IN (
            SELECT person_id
            FROM `{dataset}.cb_search_cohort`
            WHERE cohort_definition_id = {FEMALE_COHORT_ID}
          )
        ),
        all_x_variants AS (
          -- Get all X chromosome VIDs
          SELECT DISTINCT vid
          FROM `{dataset}.cb_variant_to_person`
          WHERE vid IN ('{vid_list_str}')
        ),
        x_variants_unnested AS (
          -- Unnest person_ids array
          SELECT 
            vid,
            pid
          FROM `{dataset}.cb_variant_to_person`
          CROSS JOIN UNNEST(person_ids) AS pid
          WHERE vid IN ('{vid_list_str}')
        ),
        female_variant_counts AS (
          -- Count females with each variant
          SELECT 
            vid,
            COUNT(DISTINCT x.pid) AS n_females_with_variant
          FROM x_variants_unnested x
          INNER JOIN female_ids f ON x.pid = f.person_id
          GROUP BY vid
        )
        -- Include all variants (even 0 counts)
        SELECT 
          a.vid,
          COALESCE(c.n_females_with_variant, 0) AS n_females_with_variant
        FROM all_x_variants a
        LEFT JOIN female_variant_counts c ON a.vid = c.vid
        ORDER BY n_females_with_variant DESC
        """
        
        print("\n   Executing BigQuery (may take 5-10 minutes)...")
        query_start = datetime.now()
        
        try:
            female_counts = pd.read_gbq(query, dialect='standard')
            
            print(f"✅ Query complete in {datetime.now() - query_start}")
            print(f"   Retrieved counts for {len(female_counts):,} variants")
            
            # Calculate exact allele counts and numbers (NO SUPPRESSION)
            female_counts['ac_female'] = female_counts['n_females_with_variant']  # Assume heterozygous
            female_counts['an_female'] = FEMALE_COHORT_SIZE * 2  # 2 X chromosomes per female
            female_counts['af_female'] = female_counts['ac_female'] / female_counts['an_female']
            
            # Merge with X variant info
            x_final = x_variants.merge(female_counts, left_on='vid', right_on='vid', how='left')
            
            # Save X chromosome female exact counts
            x_female_file = f"{variant_type}_X_CHROMOSOME_FEMALE_EXACT_COUNTS.csv"
            x_final.to_csv(x_female_file, index=False)
            print(f"✅ Saved: {x_female_file}")
            
            # FILTER 3: For VUS only, filter by allele frequency threshold
            if variant_type == 'VUS':
                print(f"\n🔧 Filtering VUS X chromosome variants by allele frequency...")
                n_before_af_filter = len(x_final)
                
                # Filter to AF <= 0.001% (0.00001) using female AF
                if 'af_female' in x_final.columns:
                    x_filtered = x_final[
                        (x_final['af_female'] <= VUS_AF_THRESHOLD) |
                        (x_final['af_female'].isna()) |  # Keep variants with missing AF
                        (x_final['af_female'] == 0)  # Keep variants not observed
                    ].copy()
                    
                    n_excluded_af = n_before_af_filter - len(x_filtered)
                    print(f"   Removed {n_excluded_af:,} X chromosome variants with AF > {VUS_AF_THRESHOLD} (0.001%)")
                    print(f"   Variants retained: {len(x_filtered):,}")
                    
                    # Save filtered X chromosome file
                    x_filtered_file = f"{variant_type}_X_CHROMOSOME_FEMALE_EXACT_COUNTS_RARE.csv"
                    x_filtered.to_csv(x_filtered_file, index=False)
                    print(f"✅ Saved filtered file: {x_filtered_file}")
                    
                    # Update x_final for summary
                    x_final = x_filtered
                    x_female_file = x_filtered_file
                else:
                    print(f"   ⚠️  Warning: af_female column not found, skipping AF filtering")
            
            print(f"\n📊 X chromosome summary:")
            print(f"   Variants observed in females: {(female_counts['n_females_with_variant'] > 0).sum():,}")
            print(f"   Variants NOT observed: {(female_counts['n_females_with_variant'] == 0).sum():,}")
            
        except Exception as e:
            print(f"❌ Error querying BigQuery: {e}")
            x_female_file = None
    else:
        print("   No X chromosome variants found")
        x_file = None
        x_female_file = None
    
    # Step 7: Extract autosomal variants with exact GVS counts
    print(f"\n🧬 Extracting autosomal variants with GVS allele counts...")
    
    autosomal = canonical_final[canonical_final['contig'] != 'chrX'].copy()
    print(f"   Autosomal variants: {len(autosomal):,}")
    
    # Get GVS columns
    gvs_cols = [col for col in autosomal.columns if col.startswith('gvs_')]
    print(f"   GVS columns present: {len(gvs_cols)}")
    
    # Reorder columns for easier viewing
    key_cols = [
        'vid', 'contig', 'position', 'ref_allele', 'alt_allele',
        'gene_symbol', 'consequence', 'aa_change'
    ]
    
    # Add GVS columns in organized order
    gvs_order = []
    for anc in ['all', 'afr', 'amr', 'eas', 'eur', 'mid', 'oth', 'sas']:
        for metric in ['ac', 'an', 'af']:
            col = f'gvs_{anc}_{metric}'
            if col in autosomal.columns:
                gvs_order.append(col)
    
    # Add star rating columns if present
    if 'StarRating' in autosomal.columns:
        key_cols.extend(['StarRating', 'StarRating_Numeric', 'ClinicalSignificance'])
    
    key_cols_existing = [col for col in key_cols if col in autosomal.columns]
    other_cols = [col for col in autosomal.columns if col not in key_cols_existing and col not in gvs_order]
    
    autosomal_reordered = autosomal[key_cols_existing + gvs_order + other_cols]
    
    # Save autosomal file
    autosomal_file = f"{variant_type}_AUTOSOMAL_VARIANT_LEVEL_EXACT_COUNTS.csv"
    autosomal_reordered.to_csv(autosomal_file, index=False)
    print(f"✅ Saved: {autosomal_file}")
    
    # FILTER 3: For VUS only, filter by allele frequency threshold
    if variant_type == 'VUS':
        print(f"\n🔧 Filtering VUS variants by allele frequency...")
        n_before_af_filter = len(autosomal_reordered)
        
        # Filter to AF <= 0.001% (0.00001) using GVS all ancestry
        if 'gvs_all_af' in autosomal_reordered.columns:
            autosomal_reordered['gvs_all_af_numeric'] = pd.to_numeric(autosomal_reordered['gvs_all_af'], errors='coerce')
            autosomal_filtered = autosomal_reordered[
                (autosomal_reordered['gvs_all_af_numeric'] <= VUS_AF_THRESHOLD) |
                (autosomal_reordered['gvs_all_af_numeric'].isna())  # Keep variants with missing AF
            ].copy()
            
            n_excluded_af = n_before_af_filter - len(autosomal_filtered)
            print(f"   Removed {n_excluded_af:,} variants with AF > {VUS_AF_THRESHOLD} (0.001%)")
            print(f"   Variants retained: {len(autosomal_filtered):,}")
            
            # Save filtered autosomal file
            autosomal_filtered_file = f"{variant_type}_AUTOSOMAL_VARIANT_LEVEL_EXACT_COUNTS_RARE.csv"
            autosomal_filtered.to_csv(autosomal_filtered_file, index=False)
            print(f"✅ Saved filtered file: {autosomal_filtered_file}")
            
            # Update autosomal for final summary
            autosomal_reordered = autosomal_filtered
            autosomal_file = autosomal_filtered_file
        else:
            print(f"   ⚠️  Warning: gvs_all_af column not found, skipping AF filtering")
    
    # Summary
    print(f"\n📊 {variant_type} Summary:")
    print(f"   Total variants (canonical): {len(canonical_final):,}")
    print(f"   Autosomal variants: {len(autosomal):,}")
    print(f"   X chromosome variants: {len(x_variants):,}")
    print(f"   Genes: {canonical_final['gene_symbol'].nunique()}")
    
    return {
        'canonical': canonical_file,
        'autosomal': autosomal_file,
        'x_variants': x_file,
        'x_female_counts': x_female_file
    }


################################################################################
# PROCESS ALL THREE VARIANT LISTS
################################################################################

print("\n" + "=" * 80)
print("PROCESSING ALL VARIANT LISTS")
print("=" * 80)

# 1. P/LP variants
plp_files = process_variant_list(
    variant_file='/home/jupyter/any_PLP_final_GENENAME_01_13_2026.txt',
    variant_type='PLP',
    file_separator='\t'
)

# 2. VUS variants (2+ stars only)
vus_files = process_variant_list(
    variant_file='/home/jupyter/VUS_2plus_stars_annotated.txt',
    variant_type='VUS',
    file_separator='\t'
)

# 3. LoF variants (uncomment after running LoF extraction script)
# lof_files = process_variant_list(
#     variant_file='/home/jupyter/LoF_with_VID.txt',
#     variant_type='LOF',
#     file_separator='\t'
# )


################################################################################
# FINAL SUMMARY
################################################################################

total_time = datetime.now() - start_time

print("\n" + "=" * 80)
print("🎉 EXTRACTION COMPLETE!")
print("=" * 80)
print(f"\n⏱️  Total processing time: {total_time}")

print(f"\n✅ OUTPUT FILES CREATED:")
print(f"\nP/LP Variants:")
for key, path in plp_files.items():
    if path:
        print(f"   {key}: {path}")

print(f"\nVUS Variants:")
for key, path in vus_files.items():
    if path:
        print(f"   {key}: {path}")

print(f"\n📝 NEXT STEPS:")
print(f"   1. Review output files for quality")
print(f"   2. VUS files will have _RARE.csv suffix (AF <= 0.001%)")
print(f"   3. All files exclude: SMN1, GATA1, MT-RNR1")
print(f"   4. LoF files also exclude 15 gain-of-function genes")
print(f"   5. Use variant-level files for gene-level aggregation with MEAN(AN)")
print(f"   6. Gene-level carrier frequency = sum(AC) / mean(AN)")

print("\n" + "=" * 80)
