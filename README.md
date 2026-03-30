# Replication Package

**Paper:** "Maternal Income, Education, and Children's Epigenetic Age Acceleration."
**Manuscript ID:** DEM-2024-0306
**Journal:** Demography

---

## Data Access

This study uses restricted-use data from the **Future of Families and Child Wellbeing Study (FFCWS)**. The data are available to qualified researchers through the Office of Population Research (OPR) data archive at Princeton University.

**To request access:**
1. Visit https://fragilefamilies.princeton.edu/
2. Apply for restricted-use data access through the OPR Data Archive
3. Required data modules:
   - Baseline, Year 1, Year 3, Year 5 (age 9), Year 6 (age 15), and Year 7 (age 22) surveys
   - Epigenetic clock measures from Illumina 450K and EPIC arrays at ages 9 and 15
   - Cell composition estimates (EpiDISH: epithelial and immune cell fractions) at ages 9 and 15

**Data file:** After obtaining access, construct the analytic dataset using `00_variable_creation.do` (Stata) and save as `data_2026_01_05.dta` in this directory.

---

## Software Requirements

**Stata** (for variable creation only — script 00)

**R version 4.x** with the following packages:

| Package | Purpose |
|---------|---------|
| tidyverse | Data manipulation and visualization |
| haven | Reading Stata .dta files |
| estimatr | Robust standard errors (HC2) via `lm_robust()` |
| broom | Tidying model output |
| flextable | Formatted tables for Word export |
| officer | Word document creation |
| gtsummary | Descriptive statistics tables |
| patchwork | Multi-panel figure layouts |
| splines | Natural spline models for non-linearity tests |
| openxlsx | Excel file operations |

Install all packages:
```r
install.packages(c("tidyverse", "haven", "estimatr", "broom", "flextable",
                    "officer", "gtsummary", "patchwork", "splines", "openxlsx"))
```

---

## Execution Order

Set your R working directory to this folder before running scripts. Scripts should be run in numerical order. Scripts 01-12 produce tables; scripts 13-16 produce figures.

| # | Script | Description | Produces |
|---|--------|-------------|----------|
| 00 | `00_variable_creation.do` | Construct combined analytic dataset (Stata) | `data_2026_01_05.dta` |
| 01 | `01_descriptives.R` | Table 1: Descriptive statistics | Table 1 |
| 02 | `02_income_models.R` | Tables 2-3: Income and EAA models (ages 9, 15) | Tables 2-3 |
| 03 | `03_education_models.R` | Table 4: Education and EAA models | Table 4 |
| 04 | `04_height_analysis.R` | Table S5: EAA at age 9 and self-reported height at 15 | Table S5 |
| 05 | `05_height_growth_adult.R` | Height growth (+ height at 9 control) and adult height | Supporting analysis |
| 06 | `06_menarche_correlation.R` | Table S3: EAA-menarche correlations | Table S3 |
| 07 | `07_menarche_nonlinearity.R` | Table S7: Non-linearity tests for menarche | Table S7 |
| 08 | `08_cell_composition_sensitivity.R` | Table S10: With/without cell composition adjustment | Table S10 |
| 09 | `09_cell_composition_tables.R` | Cell composition sensitivity detail tables | Table S10 detail |
| 10 | `10_maternal_covariate_sensitivity.R` | Table S11: Maternal covariate sensitivity | Table S11 |
| 11 | `11_fdr_correction.R` | Table S8: FDR correction across 12 tests | Table S8 |
| 12 | `12_supplementary_tables.R` | Assembles all supplementary tables (S1-S11) | Supplement |
| 13 | `13_menarche_figures.R` | Figure 1: Menarche-EAA panels | Figure 1 |
| 14 | `14_menarche_figure_wave4.R` | Figure 2: Refined menarche figure (Wave 4) | Figure 2 |
| 15 | `15_extremes_figures.R` | Supplementary extremes figures | Supp. Figures |
| 16 | `16_menarche_trajectory.R` | Figure 3: Menarche-EAA trajectory | Figure 3 |

---

## Notes

- **Height at age 15** is based on youth self-report (`ck6chtcm` in FFCWS), standardized to mean = 0, SD = 1.
- **Epigenetic clocks** (GrimAge, PhenoAge, DunedinPACE) are residualized on cell composition (EpiDISH epithelial + immune cell fractions) and array type (450K vs EPIC) before analysis.
- **Robust standard errors** (HC2) are used throughout via the `estimatr` package.
- Script 11 (`11_fdr_correction.R`) applies FDR correction across all 12 developmental outcome tests (3 clocks × 4 outcomes: menarche linear, menarche quadratic, height growth, adult height). It contains hardcoded p-values from scripts 05 (height growth and adult height models), 06 (menarche correlation), and 07 (menarche non-linearity). Height growth values are from the fully adjusted growth model in script 05 (Model 4: controlling for height at age 9, maternal height, sex, and race). If rerunning from scratch, verify these values match your output before running the FDR correction.
- All output is saved to an `output/` subdirectory, which is created automatically.
