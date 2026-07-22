"""
Quick Fix: Rename PLEC → PLEC1 in All of Us Variant Files
Fixes the 2 files identified by verification script
"""

import pandas as pd
import os

print("=" * 80)
print("QUICK FIX: RENAME PLEC → PLEC1 IN ALL OF US FILES")
print("=" * 80)

# Files that need fixing (identified by verification)
FILES_TO_FIX = [
    'PLP_AUTOSOMAL_VARIANT_LEVEL_EXACT_COUNTS_renamedPLEC1.csv',
    'VUS_AUTOSOMAL_VARIANT_LEVEL_EXACT_COUNTS_RARE_renamedPLEC1.csv'
]

for filename in FILES_TO_FIX:
    print(f"\n{'=' * 80}")
    print(f"Processing: {filename}")
    print(f"{'=' * 80}")
    
    if not os.path.exists(filename):
        print(f"❌ File not found: {filename}")
        continue
    
    # Read the file
    print(f"📂 Reading file...")
    df = pd.read_csv(filename, low_memory=False)
    
    print(f"   Total rows: {len(df):,}")
    
    # Check for Gene column
    if 'Gene' not in df.columns:
        print(f"❌ No 'Gene' column found in {filename}")
        continue
    
    # Count PLEC instances before
    plec_before = (df['Gene'] == 'PLEC').sum()
    plec1_before = (df['Gene'] == 'PLEC1').sum()
    
    print(f"\n📊 Before renaming:")
    print(f"   PLEC instances: {plec_before}")
    print(f"   PLEC1 instances: {plec1_before}")
    
    if plec_before == 0:
        print(f"✅ No PLEC instances to rename!")
        continue
    
    # Rename PLEC → PLEC1
    print(f"\n🔄 Renaming PLEC → PLEC1...")
    df.loc[df['Gene'] == 'PLEC', 'Gene'] = 'PLEC1'
    
    # Count after
    plec_after = (df['Gene'] == 'PLEC').sum()
    plec1_after = (df['Gene'] == 'PLEC1').sum()
    
    print(f"\n📊 After renaming:")
    print(f"   PLEC instances: {plec_after}")
    print(f"   PLEC1 instances: {plec1_after}")
    print(f"   ✅ Renamed {plec_before} instances")
    
    # Create backup
    backup_name = filename.replace('.csv', '_BACKUP_BEFORE_PLEC1_FIX.csv')
    if not os.path.exists(backup_name):
        print(f"\n💾 Creating backup: {backup_name}")
        os.rename(filename, backup_name)
    
    # Save the fixed file
    print(f"💾 Saving fixed file: {filename}")
    df.to_csv(filename, index=False)
    
    print(f"✅ Successfully fixed {filename}")

print(f"\n{'=' * 80}")
print("FIX COMPLETE")
print(f"{'=' * 80}")

print("""
Next Steps:
1. Run verification script again to confirm all files are correct:
   python verify_plec_to_plec1_renaming.py

2. If verification passes, rerun the main analysis:
   python carrier_screening_analysis_PLEC1.py

3. PLEC1 should now be included in all calculations

Backup files created:
  - PLP_AUTOSOMAL_VARIANT_LEVEL_EXACT_COUNTS_renamedPLEC1_BACKUP_BEFORE_PLEC1_FIX.csv
  - VUS_AUTOSOMAL_VARIANT_LEVEL_EXACT_COUNTS_RARE_renamedPLEC1_BACKUP_BEFORE_PLEC1_FIX.csv

If you need to revert:
  1. Delete the fixed files
  2. Rename the backup files (remove _BACKUP_BEFORE_PLEC1_FIX)
""")
