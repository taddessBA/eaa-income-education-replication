# Replication Package

**Paper:** "Maternal Income, Education, and Children's Epigenetic Age Acceleration."
**Manuscript ID:** DEM-2024-0306
**Journal:** Demography

---

## Data Access

This study uses restricted-use data from the **Future of Families and Child Wellbeing Study (FFCWS)**. The data are available to qualified researchers through the Office of Population Research (OPR) data archive at Princeton University.

**To request access:**
1. Visit [https://fragilefamilies.princeton.edu/](https://ffcws.princeton.edu/documentation)
2. Apply for restricted-use data access through the OPR Data Archive
3. Required data modules:
   - Baseline, Year 1, Year 3, Year 5 (age 9), Year 6 (age 15), and Year 7 (age 22) surveys
   - Epigenetic clock measures from Illumina 450K and EPIC arrays at ages 9 and 15
   - Cell composition estimates (EpiDISH: epithelial and immune cell fractions) at ages 9 and 15

---

## Software Requirements

- Stata: StataCorp. 2023. *Stata Statistical Software: Release 18*. College Station, TX: StataCorp LLC.
- R: R Core Team. 2024. *R: A Language and Environment for Statistical Computing* (version 4.4.1). Vienna, Austria: R Foundation for Statistical Computing. https://www.R-project.org/

---

## Execution Order

Set your R working directory to this folder before running scripts. Run scripts in numerical order. Scripts 01-04 produce tables; script 16 produces a figure.

| # | Script | Description | Produces |
|---|--------|-------------|----------|
| 00 | `00_variable_creation.do` | Construct combined analytic dataset (Stata) | Analytic dataset |
| 01 | `01_descriptives.R` | Table 1: Descriptive statistics | Table 1 |
| 02 | `02_income_models.R` | Tables 2-3: Income and EAA models (ages 9, 15) | Tables 2-3 |
| 03 | `03_education_models.R` | Table 4: Education and EAA models | Table 4 |
| 04 | `04_height_analysis.R` | Table S5: EAA at age 9 and self-reported height at 15 | Table S5 |
| 16 | `16_menarche_trajectory.R` | Figure 3: Menarche-EAA trajectory | Figure 3 |

---

## Notes

- **Epigenetic clocks** (GrimAge, PhenoAge, DunedinPACE) are residualized on cell composition (EpiDISH epithelial + immune cell fractions) and array type (450K vs EPIC) before analysis.
- **Robust standard errors** (HC2) are used throughout.
