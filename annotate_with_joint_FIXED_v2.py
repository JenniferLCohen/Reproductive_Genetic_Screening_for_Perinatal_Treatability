import pandas as pd
import pysam
import sys
import os
from datetime import datetime
import subprocess

print("="*70)
print("gnomAD v4.1 JOINT ANNOTATION - FIXED + X CHROMOSOME HANDLING")
print("="*70)
print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

PLP_FILE = 'any_PLP_MASTER.txt'
VUS_FILE = 'VUS_2plus_stars_annotated.txt'
VCF_DIR = 'gnomad_joint_vcfs'
CHECKPOINT_DIR = 'gnomad_checkpoints'
BASE_URL = 'https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/joint/'

print(f"\nLoading {PLP_FILE}...")
df_plp = pd.read_csv(PLP_FILE, sep='\t')
print(f"  Total P/LP variants: {len(df_plp):,}")

print(f"\nLoading {VUS_FILE}...")
df_vus = pd.read_csv(VUS_FILE, sep='\t')
print(f"  Total VUS variants: {len(df_vus):,}")

os.makedirs(VCF_DIR, exist_ok=True)
os.makedirs(CHECKPOINT_DIR, exist_ok=True)

chromosomes = [str(i) for i in range(1, 23)] + ['X', 'Y']

def validate_vcf(vcf_file):
    try:
        vcf = pysam.VariantFile(vcf_file)
        vcf.close()
        return True
    except:
        return False

def download_chromosome(chrom):
    vcf_file = f'{VCF_DIR}/gnomad.joint.v4.1.sites.chr{chrom}.vcf.bgz'
    tbi_file = f'{vcf_file}.tbi'
    
    if os.path.exists(vcf_file) and os.path.exists(tbi_file) and validate_vcf(vcf_file):
        print(f"  Valid VCF exists: chr{chrom}")
        return vcf_file
    
    if os.path.exists(vcf_file):
        print(f"  Removing corrupted chr{chrom} VCF...")
        os.remove(vcf_file)
    if os.path.exists(tbi_file):
        os.remove(tbi_file)
    
    print(f"  Downloading chr{chrom} VCF...")
    
    vcf_url = f'{BASE_URL}gnomad.joint.v4.1.sites.chr{chrom}.vcf.bgz'
    tbi_url = f'{BASE_URL}gnomad.joint.v4.1.sites.chr{chrom}.vcf.bgz.tbi'
    
    subprocess.run(['curl', '-L', '-o', vcf_file, vcf_url], check=True)
    subprocess.run(['curl', '-s', '-L', '-o', tbi_file, tbi_url], check=True)
    
    if not validate_vcf(vcf_file):
        print(f"  ERROR: chr{chrom} download failed validation!")
        os.remove(vcf_file)
        if os.path.exists(tbi_file):
            os.remove(tbi_file)
        raise Exception(f"Failed to download valid VCF for chr{chrom}")
    
    print(f"  ✓ chr{chrom} download validated")
    return vcf_file

def get_info_value(rec, field):
    if field not in rec.info:
        return None
    val = rec.info[field]
    if isinstance(val, tuple) and len(val) == 1:
        return val[0]
    return val

def annotate_chromosome(df, chrom, vcf_file):
    df_chr = df[(df['Chromosome'] == chrom) & (~df['Chromosome'].isin(['MT', 'M']))].copy()
    
    if len(df_chr) == 0:
        return pd.DataFrame()
    
    use_xx = (chrom == 'X')
    
    if use_xx:
        print(f"    Annotating {len(df_chr):,} variants on chr{chrom} (using XX-specific frequencies for carrier screening)...")
    else:
        print(f"    Annotating {len(df_chr):,} variants on chr{chrom}...")
    
    vcf = pysam.VariantFile(vcf_file)
    results = []
    found_count = 0
    
    for idx, row in df_chr.iterrows():
        pos = int(row['Position'])
        ref = row['Ref']
        alt = row['Alt']
        
        gnomad_data = {
            'gnomAD_AC_total': None, 'gnomAD_AN_total': None, 'gnomAD_AF_total': None,
            'gnomAD_AC_afr': None, 'gnomAD_AN_afr': None, 'gnomAD_AF_afr': None,
            'gnomAD_AC_amr': None, 'gnomAD_AN_amr': None, 'gnomAD_AF_amr': None,
            'gnomAD_AC_asj': None, 'gnomAD_AN_asj': None, 'gnomAD_AF_asj': None,
            'gnomAD_AC_eas': None, 'gnomAD_AN_eas': None, 'gnomAD_AF_eas': None,
            'gnomAD_AC_fin': None, 'gnomAD_AN_fin': None, 'gnomAD_AF_fin': None,
            'gnomAD_AC_nfe': None, 'gnomAD_AN_nfe': None, 'gnomAD_AF_nfe': None,
            'gnomAD_AC_remaining': None, 'gnomAD_AN_remaining': None, 'gnomAD_AF_remaining': None,
            'gnomAD_nhomalt_total': None,
        }
        
        try:
            for rec in vcf.fetch(f'chr{chrom}', pos-1, pos):
                if rec.pos == pos and rec.ref == ref and alt in rec.alts:
                    if use_xx:
                        gnomad_data['gnomAD_AC_total'] = get_info_value(rec, 'AC_joint_XX')
                        gnomad_data['gnomAD_AN_total'] = get_info_value(rec, 'AN_joint_XX')
                        gnomad_data['gnomAD_AF_total'] = get_info_value(rec, 'AF_joint_XX')
                        gnomad_data['gnomAD_nhomalt_total'] = get_info_value(rec, 'nhomalt_joint_XX')
                        gnomad_data['gnomAD_AC_afr'] = get_info_value(rec, 'AC_joint_afr_XX')
                        gnomad_data['gnomAD_AN_afr'] = get_info_value(rec, 'AN_joint_afr_XX')
                        gnomad_data['gnomAD_AF_afr'] = get_info_value(rec, 'AF_joint_afr_XX')
                        gnomad_data['gnomAD_AC_amr'] = get_info_value(rec, 'AC_joint_amr_XX')
                        gnomad_data['gnomAD_AN_amr'] = get_info_value(rec, 'AN_joint_amr_XX')
                        gnomad_data['gnomAD_AF_amr'] = get_info_value(rec, 'AF_joint_amr_XX')
                        gnomad_data['gnomAD_AC_asj'] = get_info_value(rec, 'AC_joint_asj_XX')
                        gnomad_data['gnomAD_AN_asj'] = get_info_value(rec, 'AN_joint_asj_XX')
                        gnomad_data['gnomAD_AF_asj'] = get_info_value(rec, 'AF_joint_asj_XX')
                        gnomad_data['gnomAD_AC_eas'] = get_info_value(rec, 'AC_joint_eas_XX')
                        gnomad_data['gnomAD_AN_eas'] = get_info_value(rec, 'AN_joint_eas_XX')
                        gnomad_data['gnomAD_AF_eas'] = get_info_value(rec, 'AF_joint_eas_XX')
                        gnomad_data['gnomAD_AC_fin'] = get_info_value(rec, 'AC_joint_fin_XX')
                        gnomad_data['gnomAD_AN_fin'] = get_info_value(rec, 'AN_joint_fin_XX')
                        gnomad_data['gnomAD_AF_fin'] = get_info_value(rec, 'AF_joint_fin_XX')
                        gnomad_data['gnomAD_AC_nfe'] = get_info_value(rec, 'AC_joint_nfe_XX')
                        gnomad_data['gnomAD_AN_nfe'] = get_info_value(rec, 'AN_joint_nfe_XX')
                        gnomad_data['gnomAD_AF_nfe'] = get_info_value(rec, 'AF_joint_nfe_XX')
                        gnomad_data['gnomAD_AC_remaining'] = get_info_value(rec, 'AC_joint_remaining_XX')
                        gnomad_data['gnomAD_AN_remaining'] = get_info_value(rec, 'AN_joint_remaining_XX')
                        gnomad_data['gnomAD_AF_remaining'] = get_info_value(rec, 'AF_joint_remaining_XX')
                    else:
                        gnomad_data['gnomAD_AC_total'] = get_info_value(rec, 'AC_joint')
                        gnomad_data['gnomAD_AN_total'] = get_info_value(rec, 'AN_joint')
                        gnomad_data['gnomAD_AF_total'] = get_info_value(rec, 'AF_joint')
                        gnomad_data['gnomAD_nhomalt_total'] = get_info_value(rec, 'nhomalt_joint')
                        gnomad_data['gnomAD_AC_afr'] = get_info_value(rec, 'AC_joint_afr')
                        gnomad_data['gnomAD_AN_afr'] = get_info_value(rec, 'AN_joint_afr')
                        gnomad_data['gnomAD_AF_afr'] = get_info_value(rec, 'AF_joint_afr')
                        gnomad_data['gnomAD_AC_amr'] = get_info_value(rec, 'AC_joint_amr')
                        gnomad_data['gnomAD_AN_amr'] = get_info_value(rec, 'AN_joint_amr')
                        gnomad_data['gnomAD_AF_amr'] = get_info_value(rec, 'AF_joint_amr')
                        gnomad_data['gnomAD_AC_asj'] = get_info_value(rec, 'AC_joint_asj')
                        gnomad_data['gnomAD_AN_asj'] = get_info_value(rec, 'AN_joint_asj')
                        gnomad_data['gnomAD_AF_asj'] = get_info_value(rec, 'AF_joint_asj')
                        gnomad_data['gnomAD_AC_eas'] = get_info_value(rec, 'AC_joint_eas')
                        gnomad_data['gnomAD_AN_eas'] = get_info_value(rec, 'AN_joint_eas')
                        gnomad_data['gnomAD_AF_eas'] = get_info_value(rec, 'AF_joint_eas')
                        gnomad_data['gnomAD_AC_fin'] = get_info_value(rec, 'AC_joint_fin')
                        gnomad_data['gnomAD_AN_fin'] = get_info_value(rec, 'AN_joint_fin')
                        gnomad_data['gnomAD_AF_fin'] = get_info_value(rec, 'AF_joint_fin')
                        gnomad_data['gnomAD_AC_nfe'] = get_info_value(rec, 'AC_joint_nfe')
                        gnomad_data['gnomAD_AN_nfe'] = get_info_value(rec, 'AN_joint_nfe')
                        gnomad_data['gnomAD_AF_nfe'] = get_info_value(rec, 'AF_joint_nfe')
                        gnomad_data['gnomAD_AC_remaining'] = get_info_value(rec, 'AC_joint_remaining')
                        gnomad_data['gnomAD_AN_remaining'] = get_info_value(rec, 'AN_joint_remaining')
                        gnomad_data['gnomAD_AF_remaining'] = get_info_value(rec, 'AF_joint_remaining')
                    
                    found_count += 1
                    break
        except:
            pass
        
        results.append(gnomad_data)
    
    for col in results[0].keys():
        df_chr[col] = [r[col] for r in results]
    
    vcf.close()
    print(f"    Found {found_count:,}/{len(df_chr):,} variants in gnomAD")
    return df_chr

print("\n" + "="*70)
print("ANNOTATING P/LP VARIANTS")
print("="*70)

for chrom in chromosomes:
    checkpoint_file = f'{CHECKPOINT_DIR}/plp_chr{chrom}.txt'
    
    if os.path.exists(checkpoint_file):
        print(f"\n✓ chr{chrom} already done")
        continue
    
    print(f"\nProcessing chr{chrom}...")
    
    try:
        vcf_file = download_chromosome(chrom)
        result = annotate_chromosome(df_plp, chrom, vcf_file)
        
        if len(result) > 0:
            result.to_csv(checkpoint_file, sep='\t', index=False)
            print(f"  ✓ Checkpoint saved")
        
        print(f"  Cleaning up...")
        if os.path.exists(vcf_file):
            os.remove(vcf_file)
        if os.path.exists(f'{vcf_file}.tbi'):
            os.remove(f'{vcf_file}.tbi')
    except Exception as e:
        print(f"  ERROR: {e}")
        continue

print("\nMerging P/LP checkpoints...")
plp_results = []
for chrom in chromosomes:
    checkpoint_file = f'{CHECKPOINT_DIR}/plp_chr{chrom}.txt'
    if os.path.exists(checkpoint_file):
        df_chr = pd.read_csv(checkpoint_file, sep='\t')
        plp_results.append(df_chr)
        print(f"  Loaded chr{chrom}: {len(df_chr):,} variants")

if plp_results:
    df_plp_annotated = pd.concat(plp_results, ignore_index=True)
    df_plp_no_annot = df_plp[~df_plp.index.isin(df_plp_annotated.index)]
    df_plp_final = pd.concat([df_plp_annotated, df_plp_no_annot], ignore_index=True)
else:
    df_plp_final = df_plp

output_plp = 'any_PLP_MASTER_JOINT.txt'
df_plp_final.to_csv(output_plp, sep='\t', index=False)

variants_in_gnomad = df_plp_final['gnomAD_AF_total'].notna().sum()
print(f"\n✓ P/LP complete: {variants_in_gnomad:,}/{len(df_plp_final):,} in gnomAD")
print(f"  Saved: {output_plp}")

print("\n" + "="*70)
print("ANNOTATING VUS VARIANTS")
print("="*70)

for chrom in chromosomes:
    checkpoint_file = f'{CHECKPOINT_DIR}/vus_chr{chrom}.txt'
    
    if os.path.exists(checkpoint_file):
        print(f"\n✓ chr{chrom} already done")
        continue
    
    print(f"\nProcessing chr{chrom}...")
    
    try:
        vcf_file = download_chromosome(chrom)
        result = annotate_chromosome(df_vus, chrom, vcf_file)
        
        if len(result) > 0:
            result.to_csv(checkpoint_file, sep='\t', index=False)
            print(f"  ✓ Checkpoint saved")
        
        print(f"  Cleaning up...")
        if os.path.exists(vcf_file):
            os.remove(vcf_file)
        if os.path.exists(f'{vcf_file}.tbi'):
            os.remove(f'{vcf_file}.tbi')
    except Exception as e:
        print(f"  ERROR: {e}")
        continue

print("\nMerging VUS checkpoints...")
vus_results = []
for chrom in chromosomes:
    checkpoint_file = f'{CHECKPOINT_DIR}/vus_chr{chrom}.txt'
    if os.path.exists(checkpoint_file):
        df_chr = pd.read_csv(checkpoint_file, sep='\t')
        vus_results.append(df_chr)
        print(f"  Loaded chr{chrom}: {len(df_chr):,} variants")

if vus_results:
    df_vus_annotated = pd.concat(vus_results, ignore_index=True)
    df_vus_no_annot = df_vus[~df_vus.index.isin(df_vus_annotated.index)]
    df_vus_final = pd.concat([df_vus_annotated, df_vus_no_annot], ignore_index=True)
else:
    df_vus_final = df_vus

output_vus = 'VUS_2plus_stars_JOINT.txt'
df_vus_final.to_csv(output_vus, sep='\t', index=False)

variants_in_gnomad = df_vus_final['gnomAD_AF_total'].notna().sum()
print(f"\n✓ VUS complete: {variants_in_gnomad:,}/{len(df_vus_final):,} in gnomAD")
print(f"  Saved: {output_vus}")

print("\n" + "="*70)
print("COMPLETE!")
print("="*70)
print(f"Finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
