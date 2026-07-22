"""
Generate ALL figures for PLEC1-corrected carrier screening analysis
Creates 7 publication-quality figures using NNS_ALL_SCENARIOS.csv
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# EXACT ORIGINAL COLORS
COLOR_AOU = '#4ECDC4'       # Teal for All of Us
COLOR_GNOMAD = '#FF6B6B'    # Coral for gnomAD

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica']
plt.rcParams['font.size'] = 12

print("=" * 80)
print("GENERATING ALL FIGURES - PLEC1 CORRECTED DATA")
print("=" * 80)

# Load data
nns_all = pd.read_csv('NNS_ALL_SCENARIOS.csv')

print(f"\n✅ Loaded NNS_ALL_SCENARIOS.csv: {len(nns_all)} rows")

# Normalize column names to lowercase for consistency
nns_all.columns = [col.lower() for col in nns_all.columns]

print(f"   Columns: {list(nns_all.columns)}")
print(f"   Panels: {sorted(nns_all['panel'].unique())}")
print(f"   Scenarios: {sorted(nns_all['scenario'].unique())}")
print(f"   Ancestries: {sorted(nns_all['ancestry'].unique())}")

# Convert carrier frequency from decimal to percentage
nns_all['combined_carrier_freq'] = nns_all['combined_carrier_freq'] * 100
print(f"   ✅ Converted carrier frequency to percentage (sample: {nns_all['combined_carrier_freq'].iloc[0]:.2f}%)")

# ============================================================================
# FIGURE 1 & 2: AoU vs gnomAD across ancestries (All Minus 4 and AR+XL)
# ============================================================================

def create_ancestry_comparison_figure(panel_code, panel_name, filename_base):
    """Create side-by-side ancestry comparison figure"""
    
    print(f"\n📊 Creating ancestry comparison: {panel_name}...")
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # Define ancestries to show
    ancestries = ['ALL', 'AFR', 'AMR', 'EAS', 'EUR']
    ancestry_labels = ['All Combined', 'African', 'Latino/Admixed', 'East Asian', 'European']
    
    # Ancestry mapping: All of Us → gnomAD
    ancestry_map = {
        'ALL': 'TOTAL',
        'AFR': 'AFR',
        'AMR': 'AMR',
        'EAS': 'EAS',
        'EUR': 'NFE'  # gnomAD uses NFE (Non-Finnish European) not EUR
    }
    
    # Left panel: P/LP only
    aou_plp = []
    aou_plp_ci_lower = []
    aou_plp_ci_upper = []
    gno_plp = []
    gno_plp_ci_lower = []
    gno_plp_ci_upper = []
    
    for anc in ancestries:
        # All of Us P/LP only
        aou_row = nns_all[
            (nns_all['panel'] == panel_code) &
            (nns_all['scenario'] == 'AoU_PLP_only') &
            (nns_all['ancestry'] == anc)
        ]
        if len(aou_row) > 0:
            nns = aou_row['nns'].iloc[0]
            aou_plp.append(nns)
            aou_plp_ci_lower.append(nns - aou_row['nns_95ci_lower'].iloc[0])
            aou_plp_ci_upper.append(aou_row['nns_95ci_upper'].iloc[0] - nns)
        else:
            aou_plp.append(0)
            aou_plp_ci_lower.append(0)
            aou_plp_ci_upper.append(0)
        
        # gnomAD P/LP only (use ancestry mapping)
        gno_anc = ancestry_map[anc]
        gno_row = nns_all[
            (nns_all['panel'] == panel_code) &
            (nns_all['scenario'] == 'gnomAD_PLP_only') &
            (nns_all['ancestry'] == gno_anc)
        ]
        if len(gno_row) > 0:
            nns = gno_row['nns'].iloc[0]
            gno_plp.append(nns)
            gno_plp_ci_lower.append(nns - gno_row['nns_95ci_lower'].iloc[0])
            gno_plp_ci_upper.append(gno_row['nns_95ci_upper'].iloc[0] - nns)
        else:
            gno_plp.append(0)
            gno_plp_ci_lower.append(0)
            gno_plp_ci_upper.append(0)
    
    # Plot left panel (P/LP only)
    x = np.arange(len(ancestry_labels))
    width = 0.35
    
    bars1 = ax1.bar(x - width/2, aou_plp, width, label='All of Us',
                    color=COLOR_AOU, alpha=0.8,
                    yerr=[aou_plp_ci_lower, aou_plp_ci_upper],
                    capsize=5, error_kw={'linewidth': 2, 'alpha': 0.3})
    bars2 = ax1.bar(x + width/2, gno_plp, width, label='gnomAD',
                    color=COLOR_GNOMAD, alpha=0.8,
                    yerr=[gno_plp_ci_lower, gno_plp_ci_upper],
                    capsize=5, error_kw={'linewidth': 2, 'alpha': 0.3})
    
    # Add value labels on bars
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            if height > 0:
                ax1.text(bar.get_x() + bar.get_width()/2., height,
                        f'{height:.1f}', ha='center', va='bottom',
                        fontsize=11, fontweight='bold')
    
    ax1.set_xlabel('ancestry', fontweight='bold', fontsize=13)
    ax1.set_ylabel('NNS', fontweight='bold', fontsize=13)
    ax1.set_title('P/LP only (≥2 stars)', fontweight='bold', fontsize=14)
    ax1.set_xticks(x)
    ax1.set_xticklabels(ancestry_labels)
    ax1.legend(loc='upper right', framealpha=0.95, fontsize=12)
    ax1.grid(axis='y', alpha=0.3, linestyle='--')
    ax1.set_ylim(0, max(max(aou_plp), max(gno_plp)) * 1.2)
    
    # Right panel: P/LP + 5% VUS
    aou_vus = []
    aou_vus_ci_lower = []
    aou_vus_ci_upper = []
    gno_vus = []
    gno_vus_ci_lower = []
    gno_vus_ci_upper = []
    
    for anc in ancestries:
        # All of Us P/LP + VUS
        aou_row = nns_all[
            (nns_all['panel'] == panel_code) &
            (nns_all['scenario'] == 'AoU_PLP_plus_5pctVUS') &
            (nns_all['ancestry'] == anc)
        ]
        if len(aou_row) > 0:
            nns = aou_row['nns'].iloc[0]
            aou_vus.append(nns)
            aou_vus_ci_lower.append(nns - aou_row['nns_95ci_lower'].iloc[0])
            aou_vus_ci_upper.append(aou_row['nns_95ci_upper'].iloc[0] - nns)
        else:
            aou_vus.append(0)
            aou_vus_ci_lower.append(0)
            aou_vus_ci_upper.append(0)
        
        # gnomAD P/LP + VUS (use ancestry mapping)
        gno_anc = ancestry_map[anc]
        gno_row = nns_all[
            (nns_all['panel'] == panel_code) &
            (nns_all['scenario'] == 'gnomAD_PLP_plus_5pctVUS') &
            (nns_all['ancestry'] == gno_anc)
        ]
        if len(gno_row) > 0:
            nns = gno_row['nns'].iloc[0]
            gno_vus.append(nns)
            gno_vus_ci_lower.append(nns - gno_row['nns_95ci_lower'].iloc[0])
            gno_vus_ci_upper.append(gno_row['nns_95ci_upper'].iloc[0] - nns)
        else:
            gno_vus.append(0)
            gno_vus_ci_lower.append(0)
            gno_vus_ci_upper.append(0)
    
    # Plot right panel
    bars1 = ax2.bar(x - width/2, aou_vus, width, label='All of Us',
                    color=COLOR_AOU, alpha=0.8,
                    yerr=[aou_vus_ci_lower, aou_vus_ci_upper],
                    capsize=5, error_kw={'linewidth': 2, 'alpha': 0.3})
    bars2 = ax2.bar(x + width/2, gno_vus, width, label='gnomAD',
                    color=COLOR_GNOMAD, alpha=0.8,
                    yerr=[gno_vus_ci_lower, gno_vus_ci_upper],
                    capsize=5, error_kw={'linewidth': 2, 'alpha': 0.3})
    
    # Add value labels
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            if height > 0:
                ax2.text(bar.get_x() + bar.get_width()/2., height,
                        f'{height:.1f}', ha='center', va='bottom',
                        fontsize=11, fontweight='bold')
    
    ax2.set_xlabel('ancestry', fontweight='bold', fontsize=13)
    ax2.set_ylabel('NNS', fontweight='bold', fontsize=13)
    ax2.set_title('P/LP + 5% rare VUS', fontweight='bold', fontsize=14)
    ax2.set_xticks(x)
    ax2.set_xticklabels(ancestry_labels)
    ax2.legend(loc='upper right', framealpha=0.95, fontsize=12)
    ax2.grid(axis='y', alpha=0.3, linestyle='--')
    ax2.set_ylim(0, max(max(aou_vus), max(gno_vus)) * 1.2)
    
    # Overall title
    title_text = f'{panel_name}: All of Us vs gnomAD'
    if 'All Minus 4' in panel_name or '289' in panel_name:
        title_text += '\nExcluded: FAM111A, SNAP25, SCN8A, KCNQ3'
    
    fig.suptitle(title_text, fontweight='bold', fontsize=15, y=0.98)
    
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(f'{filename_base}.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(f'{filename_base}.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"   ✅ Saved: {filename_base}.pdf/png")

# Create ancestry comparison figures
create_ancestry_comparison_figure(
    '4_All_minus_4',
    'All Minus 4 Genes Panel (289 genes)',
    'Fig_Compare_All_Minus4_AoU_vs_gnomAD_PLEC1'
)

create_ancestry_comparison_figure(
    '2_AR_and_XL',
    'AR + XL Panel (247 genes)',
    'Fig_Compare_ARplusXL_AoU_vs_gnomAD_PLEC1'
)

create_ancestry_comparison_figure(
    '1_All_293_genes',
    'All 293 Genes Panel',
    'Fig_Compare_All293_AoU_vs_gnomAD_PLEC1'
)

# ============================================================================
# FIGURE 3 & 4: Impact of LoF and VUS (horizontal bars)
# ============================================================================

def create_scenario_comparison_figure(panel_code, panel_name, filename_base):
    """Create horizontal bar chart comparing scenarios"""
    
    print(f"\n📊 Creating scenario comparison: {panel_name}...")
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Get data for ALL/TOTAL ancestry
    scenarios = [
        ('AoU_PLP_only', 'ALL', 'AoU P/LP only'),
        ('AoU_PLP_plus_LoF', 'ALL', 'AoU P/LP + LoF'),
        ('AoU_PLP_plus_5pctVUS', 'ALL', 'AoU P/LP + 5% VUS'),
        ('gnomAD_PLP_only', 'TOTAL', 'gnomAD P/LP only'),
        ('gnomAD_PLP_plus_5pctVUS', 'TOTAL', 'gnomAD P/LP + 5% VUS'),
    ]
    
    nns_values = []
    ci_lower = []
    ci_upper = []
    cf_values = []
    labels = []
    
    for scenario, ancestry, label in scenarios:
        row = nns_all[
            (nns_all['panel'] == panel_code) &
            (nns_all['scenario'] == scenario) &
            (nns_all['ancestry'] == ancestry)
        ]
        if len(row) > 0:
            nns = row['nns'].iloc[0]
            cf = row['combined_carrier_freq'].iloc[0]  # Already converted to percentage
            nns_values.append(nns)
            ci_lower.append(nns - row['nns_95ci_lower'].iloc[0])
            ci_upper.append(row['nns_95ci_upper'].iloc[0] - nns)
            cf_values.append(cf)
            labels.append(label)
    
    # Create horizontal bars (bottom to top)
    y_pos = np.arange(len(labels))
    
    bars = ax.barh(y_pos, nns_values, height=0.6,
                   color=[COLOR_GNOMAD if 'gnomAD' in l else COLOR_AOU for l in labels],
                   alpha=0.8,
                   xerr=[ci_lower, ci_upper],
                   capsize=5, error_kw={'linewidth': 2, 'alpha': 0.3})
    
    # Add NNS and CF labels on bars
    for i, (bar, nns, cf) in enumerate(zip(bars, nns_values, cf_values)):
        width = bar.get_width()
        ax.text(width + 0.1, bar.get_y() + bar.get_height()/2.,
               f'NNS = {nns:.1f}\nCF = {cf:.1f}%',
               ha='left', va='center', fontsize=11, fontweight='bold')
    
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=12)
    ax.set_xlabel('Number Needed to Screen (NNS)', fontweight='bold', fontsize=13)
    ax.grid(axis='x', alpha=0.3, linestyle='--')
    ax.set_xlim(0, max(nns_values) + 1)
    
    # Title
    title_text = f'{panel_name}: Impact of LoF and VUS\n(Combined Ancestry)'
    if 'All Minus 4' in panel_name or '289' in panel_name:
        title_text += '\nExcluded: FAM111A, SNAP25, SCN8A, KCNQ3'
    
    ax.set_title(title_text, fontweight='bold', fontsize=14, pad=20)
    
    plt.tight_layout()
    plt.savefig(f'{filename_base}.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(f'{filename_base}.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"   ✅ Saved: {filename_base}.pdf/png")

# Create scenario comparison figures
create_scenario_comparison_figure(
    '4_All_minus_4',
    'All Minus 4 Panel (289 genes)',
    'Fig_Compare_PLP_LoF_VUS_AllMinus4_289genes_PLEC1'
)

create_scenario_comparison_figure(
    '2_AR_and_XL',
    'AR + XL Panel (247 genes)',
    'Fig_Compare_PLP_LoF_VUS_ARplusXL_247genes_PLEC1'
)

create_scenario_comparison_figure(
    '1_All_293_genes',
    'All 293 Genes Panel',
    'Fig_Compare_PLP_LoF_VUS_All293_PLEC1'
)

# ============================================================================
# FIGURE 5: Three-panel comparison (requested earlier)
# ============================================================================

print(f"\n📊 Creating three-panel comparison...")

fig, ax = plt.subplots(figsize=(14, 6))

# Define panels (bottom to top)
panels_of_interest = {
    '5_AD_minus_4': 'AD genes',
    '2_AR_and_XL': 'AR + XL genes',
    '4_All_minus_4': 'All genes'
}
panel_order = ['5_AD_minus_4', '2_AR_and_XL', '4_All_minus_4']

# Extract data
aou_nns, aou_ci_lower, aou_ci_upper = [], [], []
gnomad_nns, gnomad_ci_lower, gnomad_ci_upper = [], [], []
panel_labels = []

for panel in panel_order:
    # All of Us
    aou_row = nns_all[
        (nns_all['panel']==panel) &
        (nns_all['scenario']=='AoU_PLP_only') &
        (nns_all['ancestry']=='ALL')
    ]
    if len(aou_row) > 0:
        nns = aou_row['nns'].iloc[0]
        aou_nns.append(nns)
        aou_ci_lower.append(nns - aou_row['nns_95ci_lower'].iloc[0])
        aou_ci_upper.append(aou_row['nns_95ci_upper'].iloc[0] - nns)
    else:
        aou_nns.append(0)
        aou_ci_lower.append(0)
        aou_ci_upper.append(0)
    
    # gnomAD
    gnomad_row = nns_all[
        (nns_all['panel']==panel) &
        (nns_all['scenario']=='gnomAD_PLP_only') &
        (nns_all['ancestry']=='TOTAL')
    ]
    if len(gnomad_row) > 0:
        nns = gnomad_row['nns'].iloc[0]
        gnomad_nns.append(nns)
        gnomad_ci_lower.append(nns - gnomad_row['nns_95ci_lower'].iloc[0])
        gnomad_ci_upper.append(gnomad_row['nns_95ci_upper'].iloc[0] - nns)
    else:
        gnomad_nns.append(0)
        gnomad_ci_lower.append(0)
        gnomad_ci_upper.append(0)
    
    panel_labels.append(panels_of_interest[panel])

# Create plot
y_pos = np.arange(len(panel_labels))
height = 0.35

bars1 = ax.barh(y_pos - height/2, aou_nns, height,
                label='All of Us', color=COLOR_AOU, alpha=0.8,
                xerr=[aou_ci_lower, aou_ci_upper],
                capsize=5, error_kw={'linewidth': 2, 'alpha': 0.3})

bars2 = ax.barh(y_pos + height/2, gnomad_nns, height,
                label='gnomAD', color=COLOR_GNOMAD, alpha=0.8,
                xerr=[gnomad_ci_lower, gnomad_ci_upper],
                capsize=5, error_kw={'linewidth': 2, 'alpha': 0.3})

# Add value labels OUTSIDE bars
for bars in [bars1, bars2]:
    for bar in bars:
        width = bar.get_width()
        if width > 0:
            ax.text(width + 0.5, bar.get_y() + bar.get_height()/2.,
                   f'{width:.1f}', ha='left', va='center',
                   fontsize=13, fontweight='bold')

ax.set_yticks(y_pos)
ax.set_yticklabels(panel_labels, fontweight='bold', fontsize=14)
ax.xaxis.tick_top()
ax.xaxis.set_label_position('top')
ax.set_xlabel('Number Needed to Screen (NNS)', fontweight='bold', fontsize=16, labelpad=10)
ax.tick_params(axis='x', labelsize=14)
fig.suptitle('Number Needed to Screen: All of Us vs gnomAD\n' +
             'P/LP ≥2 stars, Combined Ancestry',
             fontweight='bold', fontsize=14, y=0.995)
ax.legend(loc='upper right', framealpha=0.95, fontsize=13)
ax.grid(axis='x', alpha=0.3, linestyle='--')

max_nns = max(max(aou_nns), max(gnomad_nns))
max_ci = max(max([u for u in aou_ci_upper if u > 0]),
             max([u for u in gnomad_ci_upper if u > 0]))
ax.set_xlim(0, max_nns + max_ci + 5)

plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig('Fig_NNS_Three_Panel_Comparison_PLEC1.pdf', dpi=300, bbox_inches='tight')
plt.savefig('Fig_NNS_Three_Panel_Comparison_PLEC1.png', dpi=300, bbox_inches='tight')
plt.close()

print(f"   ✅ Saved: Fig_NNS_Three_Panel_Comparison_PLEC1.pdf/png")

# ============================================================================
# SUMMARY
# ============================================================================

print(f"\n{'=' * 80}")
print("✅ ALL FIGURES GENERATED")
print("=" * 80)
print(f"""
ANCESTRY COMPARISONS (AoU vs gnomAD across ancestries):
  ✅ Fig_Compare_All_Minus4_AoU_vs_gnomAD_PLEC1.pdf/png
  ✅ Fig_Compare_ARplusXL_AoU_vs_gnomAD_PLEC1.pdf/png
  ✅ Fig_Compare_All293_AoU_vs_gnomAD_PLEC1.pdf/png

SCENARIO COMPARISONS (Impact of LoF and VUS):
  ✅ Fig_Compare_PLP_LoF_VUS_AllMinus4_289genes_PLEC1.pdf/png
  ✅ Fig_Compare_PLP_LoF_VUS_ARplusXL_247genes_PLEC1.pdf/png
  ✅ Fig_Compare_PLP_LoF_VUS_All293_PLEC1.pdf/png

THREE-PANEL COMPARISON:
  ✅ Fig_NNS_Three_Panel_Comparison_PLEC1.pdf/png

DATA FILES (already created by analysis script):
  ✅ NNS_ALL_SCENARIOS.csv (200 rows)
  ✅ NNS_SUMMARY.csv (25 rows)

All figures use PLEC1-corrected data!
Colors: Teal (#4ECDC4) for All of Us, Coral (#FF6B6B) for gnomAD
""")
