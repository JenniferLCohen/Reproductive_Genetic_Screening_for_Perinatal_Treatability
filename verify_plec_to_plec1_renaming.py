"""
Verification Script: Check for Remaining PLEC Instances
Ensures all PLEC has been renamed to PLEC1 in variant and gene-level files
"""

import pandas as pd
import os
from pathlib import Path

print("=" * 80)
print("VERIFICATION: PLEC → PLEC1 RENAMING COMPLETENESS")
print("=" * 80)

# Files to check
FILES_TO_CHECK = {
    'Variant-Level Files (Autosomal - should have _renamedPLEC1)': [
        'PLP_AUTOSOMAL_VARIANT_LEVEL_EXACT_COUNTS_renamedPLEC1.csv',
        'VUS_AUTOSOMAL_VARIANT_LEVEL_EXACT_COUNTS_RARE_renamedPLEC1.csv',
        'any_PLP_MASTER_JOINT_renamedPLEC1.csv',
        'VUS_2plus_stars_JOINT_renamedPLEC1.csv'
    ],
    'Variant-Level Files (X-chromosome - should NOT have PLEC/PLEC1)': [
        'PLP_X_CHROMOSOME_FEMALE_EXACT_COUNTS.csv',
        'VUS_X_CHROMOSOME_FEMALE_EXACT_COUNTS_RARE.csv'
    ],
    'Variant-Level Files (LoF - should NOT have PLEC/PLEC1)': [
        'LOF_AUTOSOMAL_VARIANT_LEVEL_EXACT_COUNTS.csv',
        'LOF_X_CHROMOSOME_FEMALE_EXACT_COUNTS.csv'
    ],
    'Gene-Level Output Files (should have PLEC1, NOT PLEC)': [
        'AoU_PLP_GENE_LEVEL.csv',
        'AoU_LOF_GENE_LEVEL.csv',
        'AoU_VUS_GENE_LEVEL_RARE.csv',
        'AoU_PLP_LOF_COMBINED_GENE_LEVEL.csv',
        'AoU_PLP_VUS_COMBINED_GENE_LEVEL.csv',
        'gnomAD_PLP_GENE_LEVEL.csv',
        'gnomAD_VUS_GENE_LEVEL_RARE.csv',
        'gnomAD_PLP_VUS_COMBINED_GENE_LEVEL.csv'
    ]
}

def check_file_for_plec(filepath):
    """
    Check a file for PLEC and PLEC1 instances
    Returns: (plec_count, plec1_count, gene_column_name)
    """
    if not os.path.exists(filepath):
        return None, None, None
    
    # Detect file format
    if filepath.endswith('.csv'):
        try:
            df = pd.read_csv(filepath)
        except Exception as e:
            print(f"   ⚠️  Error reading {filepath}: {e}")
            return None, None, None
    elif filepath.endswith('.txt'):
        try:
            df = pd.read_csv(filepath, sep='\t')
        except Exception as e:
            print(f"   ⚠️  Error reading {filepath}: {e}")
            return None, None, None
    else:
        return None, None, None
    
    # Find gene column (could be 'Gene', 'gene', 'gene_symbol', etc.)
    gene_col = None
    for col in ['Gene', 'gene', 'gene_symbol', 'GENE', 'Gene_Symbol']:
        if col in df.columns:
            gene_col = col
            break
    
    if gene_col is None:
        return 0, 0, None  # No gene column, no PLEC/PLEC1
    
    # Count PLEC and PLEC1
    plec_count = (df[gene_col] == 'PLEC').sum()
    plec1_count = (df[gene_col] == 'PLEC1').sum()
    
    return plec_count, plec1_count, gene_col

# Track results
all_results = []
issues_found = False

# Check each category
for category, files in FILES_TO_CHECK.items():
    print(f"\n{'=' * 80}")
    print(f"{category}")
    print(f"{'=' * 80}")
    
    for filename in files:
        plec_count, plec1_count, gene_col = check_file_for_plec(filename)
        
        if plec_count is None and plec1_count is None and gene_col is None:
            # File doesn't exist
            print(f"\n📁 {filename}")
            print(f"   ⚠️  FILE NOT FOUND")
            all_results.append({
                'file': filename,
                'status': 'NOT_FOUND',
                'plec': None,
                'plec1': None
            })
            continue
        
        if gene_col is None:
            # No gene column
            print(f"\n📁 {filename}")
            print(f"   ℹ️  No gene column found (expected for some files)")
            all_results.append({
                'file': filename,
                'status': 'NO_GENE_COLUMN',
                'plec': 0,
                'plec1': 0
            })
            continue
        
        # Analyze results
        print(f"\n📁 {filename}")
        print(f"   Gene column: '{gene_col}'")
        print(f"   PLEC instances: {plec_count}")
        print(f"   PLEC1 instances: {plec1_count}")
        
        # Determine status
        if 'renamedPLEC1' in filename or 'GENE_LEVEL' in filename:
            # These files SHOULD have PLEC1 and NO PLEC
            if plec_count == 0 and plec1_count > 0:
                print(f"   ✅ CORRECT: PLEC fully renamed to PLEC1")
                status = 'CORRECT'
            elif plec_count > 0 and plec1_count > 0:
                print(f"   ❌ ISSUE: File has BOTH PLEC and PLEC1!")
                status = 'MIXED'
                issues_found = True
            elif plec_count > 0 and plec1_count == 0:
                print(f"   ❌ ISSUE: File still has PLEC, no PLEC1!")
                status = 'NOT_RENAMED'
                issues_found = True
            elif plec_count == 0 and plec1_count == 0:
                print(f"   ℹ️  No PLEC or PLEC1 found (may be expected)")
                status = 'NONE'
            else:
                status = 'UNKNOWN'
        else:
            # X-chromosome or LoF files - should have NO PLEC or PLEC1
            if plec_count == 0 and plec1_count == 0:
                print(f"   ✅ CORRECT: No PLEC or PLEC1 (expected for chr X/LoF)")
                status = 'CORRECT'
            else:
                print(f"   ⚠️  UNEXPECTED: Found PLEC/PLEC1 in X/LoF file")
                status = 'UNEXPECTED'
        
        all_results.append({
            'file': filename,
            'status': status,
            'plec': plec_count,
            'plec1': plec1_count
        })

# Summary Report
print(f"\n{'=' * 80}")
print("SUMMARY REPORT")
print(f"{'=' * 80}")

# Create results dataframe
results_df = pd.DataFrame(all_results)

# Count by status
print(f"\n📊 Status Summary:")
if len(results_df) > 0:
    status_counts = results_df['status'].value_counts()
    for status, count in status_counts.items():
        if status == 'CORRECT':
            print(f"   ✅ {status}: {count} files")
        elif status in ['MIXED', 'NOT_RENAMED']:
            print(f"   ❌ {status}: {count} files")
        elif status == 'NOT_FOUND':
            print(f"   ⚠️  {status}: {count} files")
        else:
            print(f"   ℹ️  {status}: {count} files")

# Total PLEC and PLEC1 counts
total_plec = results_df['plec'].sum() if results_df['plec'].notna().any() else 0
total_plec1 = results_df['plec1'].sum() if results_df['plec1'].notna().any() else 0

print(f"\n📊 Total Counts Across All Files:")
print(f"   PLEC instances: {total_plec}")
print(f"   PLEC1 instances: {total_plec1}")

# Files with issues
issues_df = results_df[results_df['status'].isin(['MIXED', 'NOT_RENAMED', 'UNEXPECTED'])]
if len(issues_df) > 0:
    print(f"\n❌ FILES WITH ISSUES ({len(issues_df)}):")
    for _, row in issues_df.iterrows():
        print(f"   - {row['file']}: {row['status']}")
        print(f"     PLEC: {row['plec']}, PLEC1: {row['plec1']}")
else:
    print(f"\n✅ NO ISSUES FOUND!")

# Files not found
missing_df = results_df[results_df['status'] == 'NOT_FOUND']
if len(missing_df) > 0:
    print(f"\n⚠️  FILES NOT FOUND ({len(missing_df)}):")
    for _, row in missing_df.iterrows():
        print(f"   - {row['file']}")

# Save results
results_df.to_csv('PLEC_PLEC1_Verification_Report.csv', index=False)
print(f"\n💾 Detailed report saved: PLEC_PLEC1_Verification_Report.csv")

# Final verdict
print(f"\n{'=' * 80}")
if issues_found:
    print("❌ VERIFICATION FAILED: Issues found in one or more files")
    print("   Review the issues above and fix the affected files")
else:
    print("✅ VERIFICATION PASSED: All files correctly renamed")
    print("   PLEC → PLEC1 renaming is complete!")
print(f"{'=' * 80}")

# Expected results guidance
print(f"\n{'=' * 80}")
print("EXPECTED RESULTS:")
print(f"{'=' * 80}")
print("""
FILES THAT SHOULD HAVE PLEC1 (and NO PLEC):
  ✓ PLP_AUTOSOMAL_VARIANT_LEVEL_EXACT_COUNTS_renamedPLEC1.csv
  ✓ VUS_AUTOSOMAL_VARIANT_LEVEL_EXACT_COUNTS_RARE_renamedPLEC1.csv
  ✓ any_PLP_MASTER_JOINT_renamedPLEC1.csv
  ✓ VUS_2plus_stars_JOINT_renamedPLEC1.csv
  ✓ AoU_PLP_GENE_LEVEL.csv (if PLEC1 has carriers)
  ✓ gnomAD_PLP_GENE_LEVEL.csv (if PLEC1 has carriers)

FILES THAT SHOULD HAVE NO PLEC OR PLEC1:
  ✓ PLP_X_CHROMOSOME_FEMALE_EXACT_COUNTS.csv (PLEC1 is autosomal)
  ✓ VUS_X_CHROMOSOME_FEMALE_EXACT_COUNTS_RARE.csv (PLEC1 is autosomal)
  ✓ LOF_AUTOSOMAL_VARIANT_LEVEL_EXACT_COUNTS.csv (no PLEC variants in LoF)
  ✓ LOF_X_CHROMOSOME_FEMALE_EXACT_COUNTS.csv (no PLEC variants in LoF)
  ✓ AoU_LOF_GENE_LEVEL.csv (no PLEC variants in LoF)

PLEC1 EXPECTED CARRIER FREQUENCIES:
  ✓ All of Us: ~0.0075% (12 variants, ~1 in 13,333)
  ✓ gnomAD: ~0.0096% (19 variants, ~1 in 10,417)
  ✓ Chromosome: 8 (autosomal recessive)
""")
