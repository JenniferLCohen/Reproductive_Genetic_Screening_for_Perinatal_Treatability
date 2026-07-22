# Carrier Screening Analysis — Script Reference
## 293-Gene Panel: All of Us v8 and gnomAD v4.1 JOINT

**Publication:** Carrier screening yield estimates across ancestries using All of Us and gnomAD  
**Last updated:** July 2026

---

## Overview

This repository contains all Python scripts for the carrier screening panel analysis. Scripts are organized by pipeline stage. The main analysis flow is:

```
ClinVar VCF → extract_clinvar → annotate_with_joint → AoU extraction → carrier_screening_analysis → figures
```

---

## Pipeline Stages and Scripts

### Stage 0 — Environment Setup (All of Us Researcher Workbench only)

#### `all_of_us_setup.py`
**Run first in AoU Jupyter environment.**  
Installs required packages (`pandas`, `openpyxl`, `google-cloud-bigquery`, `pandas-gbq`) and verifies that required input files are present at `/home/jupyter/`.

**Required input files:**
- `any_PLP_final_GENENAME_01_13_2026.txt` — ClinVar P/LP variants
- `VUS_2plus_stars_annotated.txt` — ClinVar VUS (≥2 stars)
- `mmc2_unmerged.xlsx` — 293-gene panel list

---

### Stage 1 — ClinVar Variant Extraction (Local)

#### `extract_clinvar_variants_by_gene_STANDALONE.py`
Extracts P/LP and VUS variants from a ClinVar VCF file for the 293-gene panel.

**Key logic:**
- Gene assignment uses the primary (first-listed) gene in the GENEINFO field
- P/LP: `'Pathogenic' in CLNSIG AND 'Benign' not in CLNSIG AND 'Conflicting' not in CLNSIG`
- VUS: `'Uncertain_significance' in CLNSIG`
- Star filter: ≥2 stars (applied after classification)
- Of 249,031 panel-relevant ClinVar records, 1,661 (0.67%) had a panel gene only as a secondary annotation; after star filtering, 23 P/LP and 95 VUS variants were not captured — unlikely to materially affect results

**Usage:**
```bash
python extract_clinvar_variants_by_gene_STANDALONE.py \
    --gene-list mmc2_unmerged.xlsx \
    --clinvar-vcf clinvar.vcf.gz \
    --min-stars 2 \
    --output-plp ClinVar_PLP_2star.csv \
    --output-vus ClinVar_VUS_2star.csv
```

**Outputs:** `ClinVar_PLP_2star_for_analysis.csv`, `ClinVar_VUS_2star_for_analysis.csv`

---

### Stage 2 — gnomAD v4.1 JOINT Annotation (Local)

#### `annotate_with_joint_FIXED_v2.py`
**Primary gnomAD annotation script used in the analysis.**  
Annotates ClinVar P/LP and VUS variants with gnomAD v4.1 JOINT allele frequencies by downloading and querying per-chromosome VCF files. Implements checkpoint-based processing (one file per chromosome) for robustness.

**Key feature:** For X chromosome variants, uses XX-specific INFO fields (`AC_joint_XX`, `AN_joint_XX`, `AF_joint_XX`) rather than total counts. For autosomal variants, uses standard joint fields (`AC_joint`, `AN_joint`, etc.). This means `gnomAD_AN_total` for X-linked variants reflects female-only allele number.

**Requirements:** `pysam`, gnomAD v4.1 JOINT VCF files (downloaded per-chromosome from `https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/joint/`)

**Outputs:** `any_PLP_MASTER_JOINT.txt`, `VUS_2plus_stars_JOINT.txt`

#### `annotate_with_gnomad_STANDALONE.py`
Standalone/simplified version of the gnomAD annotation script for reference. Uses standard (non-XX-specific) gnomAD fields for all chromosomes. **Note:** This was not the script used for the primary analysis — use `annotate_with_joint_FIXED_v2.py` instead.

---

### Stage 3 — All of Us Variant Extraction (AoU Workbench)

#### `all_of_us_lof_extraction_updated.py`
**Run in AoU Jupyter environment (Step 1 of AoU extraction).**  
Uses Hail to filter the AoU Variant Annotation Table (VAT, v8) for loss-of-function consequences across all transcripts, then exports to Google Cloud Storage.

**LoF consequences included:**
- `transcript_ablation`, `splice_acceptor_variant`, `splice_donor_variant`
- `stop_gained`, `frameshift_variant`, `stop_lost`, `start_lost`

**Gene exclusions:**
- All variant types: SMN1, GATA1, MT-RNR1
- LoF-specific (gain-of-function genes): GLUD1, FAM111A, CYP11B1, CLCN2, KCNJ5, KCNJ11, STAT3, TSHR, CACNA1D, GUCY2C, NLRC4, NLRP3, RAC2, CHRNB1, CACNA1C

**Output:** `LoF_with_VID.txt` (passed to `all_of_us_variant_extraction_updated.py`)

#### `all_of_us_variant_extraction_updated.py`
**Run in AoU Jupyter environment (Step 2 of AoU extraction).**  
Main AoU extraction pipeline. Processes P/LP, VUS, and LoF variants through canonical transcript deduplication (MANE Select or MANE Plus Clinical), female-specific X chromosome allele counting, and ancestry-stratified allele frequency calculation.

**Key implementation details:**
- MANE Select/MANE Plus Clinical deduplication applied to all variant types
- X chromosome: AC = count of distinct female participants carrying each variant (assumes heterozygosity); AN = n_females × 2 = 500,142
- VUS AF filter: <0.001% (applied using AoU female AF for X chromosome, AoU total AF for autosomal)
- Cohort: Females_XX_with_WGS (cohort ID 137526, n=250,071) for X chromosome

**Outputs (6 files):**
- `PLP_AUTOSOMAL_VARIANT_LEVEL_EXACT_COUNTS_renamedPLEC1.csv`
- `PLP_X_CHROMOSOME_FEMALE_EXACT_COUNTS.csv`
- `LOF_AUTOSOMAL_VARIANT_LEVEL_EXACT_COUNTS.csv`
- `LOF_X_CHROMOSOME_FEMALE_EXACT_COUNTS.csv`
- `VUS_AUTOSOMAL_VARIANT_LEVEL_EXACT_COUNTS_RARE_renamedPLEC1.csv`
- `VUS_X_CHROMOSOME_FEMALE_EXACT_COUNTS_RARE.csv`

---

### Stage 4 — Main Carrier Frequency and NNS Analysis (Local)

#### `carrier_screening_analysis_PLEC1_v2.py`
**Primary analysis script.** Calculates gene-level carrier frequencies and panel-level NNS for all 5 scenarios across 8 ancestries.

**Key formulas:**
- AF_gene = sum(AC_variants) / mean(AN_variants)
- Autosomal CF = 2 × AF × (1 − AF) [Hardy-Weinberg]
- X-linked CF = 2 × AF × (1 − AF) [HWE; AN is allele-based in both AoU and gnomAD]
- Panel CF_combined = 1 − ∏(1 − CF_gene) [product rule]
- NNS = 1 / CF_combined

**5 scenarios:**
1. AoU P/LP only (≥2 stars)
2. AoU P/LP + LoF
3. AoU P/LP + 5% rare VUS (AF <0.001%)
4. gnomAD P/LP only (≥2 stars)
5. gnomAD P/LP + 5% rare VUS (AF <0.001%)

**Input files (9):**
- `PLP_AUTOSOMAL_VARIANT_LEVEL_EXACT_COUNTS_renamedPLEC1.csv`
- `PLP_X_CHROMOSOME_FEMALE_EXACT_COUNTS.csv`
- `LOF_AUTOSOMAL_VARIANT_LEVEL_EXACT_COUNTS.csv`
- `LOF_X_CHROMOSOME_FEMALE_EXACT_COUNTS.csv`
- `VUS_AUTOSOMAL_VARIANT_LEVEL_EXACT_COUNTS_RARE_renamedPLEC1.csv`
- `VUS_X_CHROMOSOME_FEMALE_EXACT_COUNTS_RARE.csv`
- `any_PLP_MASTER_JOINT_renamedPLEC1.csv` (gnomAD P/LP)
- `VUS_2plus_stars_JOINT_renamedPLEC1.csv` (gnomAD VUS)
- `mmc2_unmerged.xlsx` (gene panel list)

**Key outputs:**
- `GENE_CLASSIFICATIONS_293.csv`
- `AoU_PLP_GENE_LEVEL.csv`, `AoU_LOF_GENE_LEVEL.csv`, `AoU_VUS_GENE_LEVEL_RARE.csv`
- `AoU_PLP_LOF_COMBINED_GENE_LEVEL.csv`, `AoU_PLP_VUS_COMBINED_GENE_LEVEL.csv`
- `gnomAD_PLP_GENE_LEVEL.csv`, `gnomAD_VUS_GENE_LEVEL_RARE.csv`, `gnomAD_PLP_VUS_COMBINED_GENE_LEVEL.csv`
- `NNS_ALL_SCENARIOS.csv`, `NNS_SUMMARY.csv`

**Usage:**
```bash
python carrier_screening_analysis_PLEC1_v2.py 2>&1 | tee run_log_v2.txt
```
Run time: ~5–10 minutes (bootstrap iterations).

---

### Stage 5 — Figure Generation (Local)

#### `create_ALL_figures_PLEC1.py`
Generates all 7 main publication figures from `NNS_ALL_SCENARIOS.csv`.

**Figures produced (PDF + PNG each):**
1. `Fig_Compare_All293_AoU_vs_gnomAD_PLEC1` — Ancestry comparison, 293-gene panel
2. `Fig_Compare_ARplusXL_AoU_vs_gnomAD_PLEC1` — Ancestry comparison, AR+XL panel
3. `Fig_Compare_All_Minus4_AoU_vs_gnomAD_PLEC1` — Ancestry comparison, All-minus-4 panel
4. `Fig_Compare_PLP_LoF_VUS_All293_PLEC1` — Scenario comparison, 293-gene panel
5. `Fig_Compare_PLP_LoF_VUS_ARplusXL_247genes_PLEC1` — Scenario comparison, AR+XL
6. `Fig_Compare_PLP_LoF_VUS_AllMinus4_289genes_PLEC1` — Scenario comparison, All-minus-4
7. `Fig_NNS_Three_Panel_Comparison_PLEC1` — Three-panel NNS summary

**Colors:** All of Us = `#4ECDC4` (teal), gnomAD = `#FF6B6B` (coral)

**Usage:**
```bash
python create_ALL_figures_PLEC1.py
```

#### `create_cumulative_contribution_plot.py`
Generates cumulative carrier frequency contribution figure (`Fig_Cumulative_Gene_Contribution_PLEC1.pdf/.png`). Shows the share of total 293-gene panel carrier frequency explained by the top 1, 3, 5, 10, and 20 genes for both AoU and gnomAD.

**Input files (from `v2_XL_corrected/`):**
- `AoU_PLP_GENE_LEVEL.csv`
- `gnomAD_PLP_GENE_LEVEL.csv`
- `GENE_CLASSIFICATIONS_293.csv`

---

### Stage 6 — Supplementary Analyses (Local)

#### `analyze_high_carrier_freq_genes.py`
Identifies genes with carrier frequency ≥ 1/200 (0.5%) in AoU and gnomAD, across all ancestries. Confirmed: no X-linked genes meet this threshold.

**Outputs:** `AoU_High_CF_Genes_1in200.csv`, `gnomAD_High_CF_Genes_1in200.csv`, `High_CF_Genes_Comparison_1in200.csv`

#### `build_prenatal_gene_table_v3.py`
**Correct script for Table S13** — builds the per-gene comparison table for the 53 prenatal-treatable genes with columns: Gene, gnomAD CF (%), AoU CF (%), Difference (%), AoU/gnomAD CF Ratio. Uses the correct gene list (COL1A1, COL1A2, GAA, LIPA, IDUA, IDS, ... CFTR).

**Usage:**
```bash
python build_prenatal_gene_table_v3.py \
    --aou-plp AoU_PLP_GENE_LEVEL.csv \
    --gnomad-plp gnomAD_PLP_GENE_LEVEL.csv \
    --output Prenatal_Treatable_Gene_Level_Comparison_v3.csv
```

#### `prenatal_treatable_analysis_STANDALONE-corrected.py`
Calculates combined panel-level NNS for prenatal-treatable genes (bootstrap CI). 

#### `structural_variant_coverage_analysis_JOINT.py` and `structural_variant_coverage_AllOfUs.py`
These scripts are retained for reference but were not part of the primary published analysis. The published approach to structural variant limitation uses ClinGen dosage sensitivity flagging: genes with ClinGen haploinsufficiency or triplosensitivity score = 3 ("Sufficient Evidence") were designated as `SV_UNDERCOUNT_LIKELY`, reflecting the expectation that P/LP copy-number or dosage-altering variants may be systematically underrepresented in short-read population datasets. Results are reported in Supplemental Tables 9 and 10.

---

### Stage 7 — QC, Verification and Investigation Scripts (Local)

#### `verify_ad_gene_carrier_frequencies.py`
Verifies that all autosomal dominant (AD) genes have carrier frequency < 0.5% (1/200) in both AoU and gnomAD. Confirms the AD panel yield is driven by a small number of genes.

**Outputs:** `AD_Gene_Carrier_Frequencies.csv`

#### `fix_plec_to_plec1.py`
One-time utility: renamed `PLEC` → `PLEC1` in AoU autosomal P/LP and VUS variant files to correct gene symbol. Creates backups before modifying.

**Files fixed:**
- `PLP_AUTOSOMAL_VARIANT_LEVEL_EXACT_COUNTS_renamedPLEC1.csv`
- `VUS_AUTOSOMAL_VARIANT_LEVEL_EXACT_COUNTS_RARE_renamedPLEC1.csv`

#### `check_missed_variants.py`
Quantifies the impact of GENEINFO primary-gene-only parsing in ClinVar. Scans `clinvar.vcf.gz` to find variants where a panel gene appears only as a secondary GENEINFO annotation, then applies star and P/LP filters to determine the true number missed.

**Key finding:** 23 P/LP and 95 VUS variants across 25 panel genes missed after ≥2 star filter — unlikely to materially affect results.

**Input:** `clinvar.vcf.gz`, `GENE_CLASSIFICATIONS_293.csv`  
**Outputs:** `missed_plp_variants.csv`, `missed_vus_variants.csv`

#### `check_prrt2_discrepancy.py`
Investigates the ~9-fold PRRT2 carrier frequency discrepancy between AoU (CF = 0.05%) and gnomAD (CF = 0.46%). Performs variant-level overlap analysis and AF comparison.

**Finding:** The discrepancy is driven by two VCF representations of the PRRT2 c.649dup founder variant (dbSNP rs587778771) at chr16:29813694, which has 7.3× and 31.6× higher AF in gnomAD than AoU due to European ancestry enrichment in gnomAD.

**Input files:** `PLP_AUTOSOMAL_VARIANT_LEVEL_EXACT_COUNTS_renamedPLEC1.csv`, `any_PLP_MASTER_JOINT_renamedPLEC1.csv`

---

## Input Data Requirements Summary

| File | Source | Stage needed |
|------|--------|-------------|
| `clinvar.vcf.gz` | ClinVar FTP | Stage 1 |
| `mmc2_unmerged.xlsx` | Supplementary Table 2 | Stages 1, 3, 4 |
| `any_PLP_final_GENENAME_01_13_2026.txt` | ClinVar P/LP extraction | Stage 3 (AoU) |
| `VUS_2plus_stars_annotated.txt` | ClinVar VUS extraction | Stage 3 (AoU) |
| gnomAD v4.1 JOINT VCF files (per chromosome) | gnomAD downloads | Stage 2 |

---

## Key Methodological Notes

**Carrier frequency formula (both autosomal and X-linked):**
- CF = 2 × AF × (1 − AF) [Hardy-Weinberg]
- This applies to X-linked genes because AN is allele-based (n_females × 2) in both datasets

**gnomAD X-linked variants:**
- `gnomAD_AN_total` for X chromosome = XX-specific allele number (`AN_joint_XX`) 
- This is populated by `annotate_with_joint_FIXED_v2.py`, not `annotate_with_gnomad_STANDALONE.py`

**AoU X-linked female AC assumption:**
- AC = count of distinct female carriers (heterozygosity assumed)
- Homozygous females would have AC underestimated by 1 per variant

**VUS rarity filter:**
- AoU: applied using AoU female AF (`af_female`) for X chromosome during extraction
- gnomAD: applied using `gnomAD_AF_total` in `carrier_screening_analysis_PLEC1_v2.py`

**Panel-level NNS:**
- Combines AR, XL, and AD gene carrier frequencies via the product rule
- CF_combined = 1 − ∏(1 − CF_gene)

---

## Python Dependencies

```
pandas >= 2.1.0
numpy >= 1.24.0
scipy >= 1.11.0
matplotlib >= 3.7.0
openpyxl
pysam          # for annotate_with_joint_FIXED_v2.py
hail           # for AoU scripts (AoU Workbench only)
google-cloud-bigquery  # for AoU scripts
pandas-gbq     # for AoU scripts
```

Install locally:
```bash
pip install pandas numpy scipy matplotlib openpyxl pysam
```
