"""
Cumulative carrier frequency contribution plot
Shows share of total panel qualifying heterozygote yield explained by
top 1, 3, 5, 10, and 20 genes — for both AoU and gnomAD, AR+XL panel.

Output: Fig_Cumulative_Gene_Contribution_PLEC1.pdf/.png
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import os

# ── Paths ──────────────────────────────────────────────────────────────────
FIGURES_DIR = os.path.join(
    os.path.expanduser('~'),
    'Desktop',
    'Figures_for_Aim1_publication 2',
    'send_to_coauthors',
    'Figures',
    'v2_XL_corrected'
)
AOU_FILE    = os.path.join(FIGURES_DIR, 'AoU_PLP_GENE_LEVEL.csv')
GNOMAD_FILE = os.path.join(FIGURES_DIR, 'gnomAD_PLP_GENE_LEVEL.csv')
GC_FILE     = os.path.join(FIGURES_DIR, 'GENE_CLASSIFICATIONS_293.csv')

# ── Colors (matching existing figures) ────────────────────────────────────
COLOR_AOU    = '#4ECDC4'
COLOR_GNOMAD = '#FF6B6B'
GRAY_LIGHT   = '#F5F5F5'
GRAY_MID     = '#CCCCCC'
TEXT_COLOR   = '#2D2D2D'

MILESTONES = [1, 3, 5, 10, 20]

plt.rcParams['font.family']     = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica']
plt.rcParams['font.size']       = 12

# ── Load data ──────────────────────────────────────────────────────────────
aou    = pd.read_csv(AOU_FILE)
gnomad = pd.read_csv(GNOMAD_FILE)
gc     = pd.read_csv(GC_FILE)

ar_xl_genes = set(gc[gc['inheritance'].isin(['AR', 'XL'])]['gene_symbol'])
all_genes = set(gc['gene_symbol'])

# Filter to ALL / TOTAL ancestry only — include ALL inheritance patterns
aou_arxl = (aou[(aou['ancestry'] == 'ALL') &
                (aou['gene_symbol'].isin(all_genes))]
            .sort_values('carrier_frequency', ascending=False)
            .reset_index(drop=True))

gno_arxl = (gnomad[(gnomad['ancestry'] == 'TOTAL') &
                   (gnomad['gene_symbol'].isin(all_genes))]
            .sort_values('carrier_frequency', ascending=False)
            .reset_index(drop=True))

def cumulative_contribution(df, cf_col='carrier_frequency', max_genes=20):
    """
    Returns arrays of (gene_rank, cumulative_pct_of_total)
    using the product rule for combined CF.
    """
    cfs = df[cf_col].values[:max_genes]
    total_cf = 1 - np.prod(1 - df[cf_col].values)

    cum_cf   = np.array([1 - np.prod(1 - cfs[:n+1]) for n in range(len(cfs))])
    cum_pct  = cum_cf / total_cf * 100
    ranks    = np.arange(1, len(cfs) + 1)
    return ranks, cum_pct, total_cf

aou_ranks, aou_pct, aou_total   = cumulative_contribution(aou_arxl)
gno_ranks, gno_pct, gno_total   = cumulative_contribution(gno_arxl)

# ── Figure ─────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 6.5))
fig.patch.set_facecolor('white')
ax.set_facecolor('white')

# Shade milestone columns
for x in MILESTONES:
    ax.axvline(x, color=GRAY_MID, linewidth=0.8, linestyle='--', zorder=1)

# Main lines
ax.plot(aou_ranks, aou_pct,
        color=COLOR_AOU, linewidth=2.5, zorder=3,
        label=f'All of Us  (total panel CF = {aou_total*100:.1f}%)')
ax.plot(gno_ranks, gno_pct,
        color=COLOR_GNOMAD, linewidth=2.5, zorder=3,
        label=f'gnomAD  (total panel CF = {gno_total*100:.1f}%)')

# Milestone dots + annotations
for n in MILESTONES:
    if n <= len(aou_pct):
        y_aou = aou_pct[n-1]
        ax.scatter(n, y_aou, color=COLOR_AOU,   s=60, zorder=4)
    if n <= len(gno_pct):
        y_gno = gno_pct[n-1]
        ax.scatter(n, y_gno, color=COLOR_GNOMAD, s=60, zorder=4)

# Annotate top gene names
for i, row in aou_arxl.head(5).iterrows():
    ax.annotate(row['gene_symbol'],
                xy=(i+1, aou_pct[i]),
                xytext=(4, 4), textcoords='offset points',
                fontsize=8, color=COLOR_AOU, fontweight='bold')

for i, row in gno_arxl.head(5).iterrows():
    ax.annotate(row['gene_symbol'],
                xy=(i+1, gno_pct[i]),
                xytext=(4, -12), textcoords='offset points',
                fontsize=8, color=COLOR_GNOMAD, fontweight='bold')

# Reference lines at 50% and 80%
for y_ref, label in [(50, '50%'), (80, '80%')]:
    ax.axhline(y_ref, color=GRAY_MID, linewidth=0.8,
               linestyle=':', zorder=1)
    ax.text(20.3, y_ref, label, va='center',
            fontsize=8, color='#888888')

# Milestone x-tick table annotation
milestone_data = []
for n in MILESTONES:
    aou_v = f"{aou_pct[n-1]:.0f}%" if n <= len(aou_pct) else '—'
    gno_v = f"{gno_pct[n-1]:.0f}%" if n <= len(gno_pct) else '—'
    milestone_data.append((n, aou_v, gno_v))

# Table below x-axis
table_y = -18
col_aou_color  = COLOR_AOU
col_gno_color  = COLOR_GNOMAD

ax.text(0.02, -0.14, 'Top N genes:',
        transform=ax.transAxes, fontsize=8.5,
        color=TEXT_COLOR, fontweight='bold')
ax.text(0.02, -0.19, 'AoU cumulative %:',
        transform=ax.transAxes, fontsize=8.5, color=COLOR_AOU)
ax.text(0.02, -0.24, 'gnomAD cumulative %:',
        transform=ax.transAxes, fontsize=8.5, color=COLOR_GNOMAD)

x_positions = [0.20, 0.32, 0.44, 0.60, 0.78]
for (n, aou_v, gno_v), xpos in zip(milestone_data, x_positions):
    ax.text(xpos, -0.14, f'Top {n}',
            transform=ax.transAxes, fontsize=8.5,
            color=TEXT_COLOR, ha='center', fontweight='bold')
    ax.text(xpos, -0.19, aou_v,
            transform=ax.transAxes, fontsize=8.5,
            color=COLOR_AOU, ha='center')
    ax.text(xpos, -0.24, gno_v,
            transform=ax.transAxes, fontsize=8.5,
            color=COLOR_GNOMAD, ha='center')

# Axes formatting
ax.set_xlim(0.5, 20.5)
ax.set_ylim(0, 105)
ax.set_xticks(MILESTONES)
ax.set_xticklabels([str(n) for n in MILESTONES])
ax.yaxis.set_major_formatter(mticker.FormatStrFormatter('%g%%'))
ax.set_xlabel('Number of genes (ranked by carrier frequency)', fontsize=12)
ax.set_ylabel('Cumulative % of total 293-gene panel carrier frequency', fontsize=12)
ax.set_title(
    'Cumulative Contribution of Top Genes to 293-Gene Panel Carrier Frequency\n'
    '(P/LP variants ≥2 stars, combined ancestry)',
    fontsize=13, fontweight='bold', pad=14, color=TEXT_COLOR
)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='both', labelsize=10)

legend = ax.legend(frameon=True, fontsize=10, loc='lower right',
                   framealpha=0.9, edgecolor=GRAY_MID)

plt.tight_layout()
plt.subplots_adjust(bottom=0.28)

# ── Save ───────────────────────────────────────────────────────────────────
for ext in ['pdf', 'png']:
    outpath = os.path.join(
        FIGURES_DIR,
        f'Fig_Cumulative_Gene_Contribution_PLEC1.{ext}'
    )
    plt.savefig(outpath, dpi=300, bbox_inches='tight',
                facecolor='white')
    print(f'Saved: {outpath}')

plt.close()

# ── Console summary ────────────────────────────────────────────────────────
print('\n── Cumulative contribution summary (293-gene panel) ──────────────────')
print(f'{"Top N":>6}  {"AoU %":>8}  {"gnomAD %":>10}')
for n in MILESTONES:
    aou_v = f"{aou_pct[n-1]:.1f}%" if n <= len(aou_pct) else '—'
    gno_v = f"{gno_pct[n-1]:.1f}%" if n <= len(gno_pct) else '—'
    print(f'{n:>6}  {aou_v:>8}  {gno_v:>10}')

print(f'\nAoU  top 5 genes: {", ".join(aou_arxl["gene_symbol"].head(5).tolist())}')
print(f'gnomAD top 5 genes: {", ".join(gno_arxl["gene_symbol"].head(5).tolist())}')
