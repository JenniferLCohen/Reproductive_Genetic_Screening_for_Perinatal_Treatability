"""
All of Us - Setup and File Verification
========================================

Run this script FIRST to:
1. Install required Python packages
2. Verify all required files are uploaded
3. Check file formats

Author: Setup script
Date: 2026-02-13
"""

import os
import sys

print("=" * 80)
print("ALL OF US EXTRACTION - SETUP & VERIFICATION")
print("=" * 80)

################################################################################
# INSTALL REQUIRED PACKAGES
################################################################################

print("\n📦 Installing required packages...")

packages = [
    'pandas',
    'openpyxl',  # For reading Excel files
    'google-cloud-bigquery',
    'pandas-gbq'
]

for package in packages:
    print(f"\n   Installing {package}...")
    os.system(f"pip install {package} --break-system-packages --quiet")

print("\n✅ All packages installed")

################################################################################
# VERIFY FILE UPLOADS
################################################################################

print("\n" + "=" * 80)
print("📂 VERIFYING UPLOADED FILES")
print("=" * 80)

required_files = {
    'P/LP variants': '/home/jupyter/any_PLP_final_GENENAME_01_13_2026.txt',
    'VUS variants (2+ stars)': '/home/jupyter/VUS_2plus_stars_annotated.txt',
    'Gene list (Excel)': '/home/jupyter/mmc2_unmerged.xlsx'
}

all_present = True

for description, filepath in required_files.items():
    if os.path.exists(filepath):
        size_mb = os.path.getsize(filepath) / (1024 * 1024)
        print(f"✅ {description}: {filepath}")
        print(f"   Size: {size_mb:.2f} MB")
    else:
        print(f"❌ MISSING: {description}")
        print(f"   Expected location: {filepath}")
        all_present = False

################################################################################
# CHECK FILE FORMATS
################################################################################

if all_present:
    print("\n" + "=" * 80)
    print("🔍 CHECKING FILE FORMATS")
    print("=" * 80)
    
    import pandas as pd
    
    # Check P/LP file
    print("\n📋 P/LP File:")
    try:
        plp_df = pd.read_csv(required_files['P/LP variants'], sep='\t', nrows=5)
        print(f"   Columns: {plp_df.columns.tolist()}")
        print(f"   First row preview:")
        print(plp_df.head(1).to_dict('records'))
        
        # Check for required columns
        required_cols = ['Chromosome', 'Position', 'Ref', 'Alt']
        missing = [c for c in required_cols if c not in plp_df.columns]
        if missing:
            print(f"   ⚠️  Missing columns: {missing}")
        else:
            print(f"   ✅ All required columns present")
    except Exception as e:
        print(f"   ❌ Error reading file: {e}")
    
    # Check VUS file
    print("\n📋 VUS File (2+ stars):")
    try:
        vus_df = pd.read_csv(required_files['VUS variants (2+ stars)'], sep='\t', nrows=5)
        print(f"   Columns: {vus_df.columns.tolist()}")
        print(f"   First row preview:")
        print(vus_df.head(1).to_dict('records'))
        
        # Check for required columns
        required_cols = ['Chromosome', 'Position', 'Ref', 'Alt']
        missing = [c for c in required_cols if c not in vus_df.columns]
        if missing:
            print(f"   ⚠️  Missing columns: {missing}")
        else:
            print(f"   ✅ All required columns present")
    except Exception as e:
        print(f"   ❌ Error reading file: {e}")
    
    # Check Excel file
    print("\n📋 Gene List (Excel):")
    try:
        gene_df = pd.read_excel(required_files['Gene list (Excel)'])
        print(f"   Shape: {gene_df.shape}")
        print(f"   Columns: {gene_df.columns.tolist()}")
        
        # Show first column (likely contains genes)
        first_col = gene_df.columns[0]
        print(f"   First column '{first_col}' preview:")
        print(f"   {gene_df[first_col].head(10).tolist()}")
        print(f"   Total unique values: {gene_df[first_col].nunique()}")
        
    except Exception as e:
        print(f"   ❌ Error reading Excel file: {e}")
        print(f"   💡 You may need to install openpyxl:")
        print(f"      pip install openpyxl --break-system-packages")

################################################################################
# FINAL INSTRUCTIONS
################################################################################

print("\n" + "=" * 80)
print("📝 NEXT STEPS")
print("=" * 80)

if all_present:
    print("\n✅ All files are present!")
    print("\n🚀 Ready to run extraction scripts:")
    print("\n1. Extract LoF variants:")
    print("   %run all_of_us_lof_extraction_updated.py")
    print("\n2. Extract all variant types (P/LP, VUS, LoF):")
    print("   %run all_of_us_variant_extraction_updated.py")
    print("\n⏰ Total expected runtime: ~2-4 hours")
    print("   (Most time is VAT loading, which happens once)")
else:
    print("\n❌ Some files are missing!")
    print("\n📤 Please upload the following to /home/jupyter/:")
    for description, filepath in required_files.items():
        if not os.path.exists(filepath):
            filename = os.path.basename(filepath)
            print(f"   - {filename} ({description})")
    
    print("\n💡 How to upload:")
    print("   1. In Jupyter, click the Upload button (⬆️)")
    print("   2. Select your files")
    print("   3. Upload to /home/jupyter/")
    print("   4. Re-run this setup script to verify")

print("\n" + "=" * 80)
