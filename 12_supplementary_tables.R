# Replication code for: "Maternal Income, Education, and Children's Epigenetic Age Acceleration"
# Manuscript ID: DEM-2024-0306 | Journal: Demography
# Input: data_2026_01_05.dta (FFCWS restricted-use data)
# =============================================================================
# Supplementary Tables: All tables in one Word document
# =============================================================================

library(tidyverse)
library(haven)
library(estimatr)
library(broom)
library(flextable)
library(officer)

# Load data
data <- read_dta("data_2026_01_05.dta")

cat("=== DATA LOADED ===\n")

# Create output directory
if (!dir.exists("output")) dir.create("output")

# Initialize Word document
doc <- read_docx()

# =============================================================================
# Helper functions
# =============================================================================

format_coef <- function(estimate, se, p_value) {
  stars <- case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    p_value < 0.10 ~ "\u2020",
    TRUE ~ ""
  )
  coef_str <- sprintf("%.3f%s", estimate, stars)
  se_str <- sprintf("(%.3f)", se)
  return(list(coef = coef_str, se = se_str))
}

run_reg <- function(formula, data) {
  lm_robust(as.formula(formula),
            data = data,
            clusters = data$idnum,
            se_type = "stata")
}

# =============================================================================
# TABLE S1: Detailed Description of Epigenetic Clocks
# =============================================================================

cat("\n=== TABLE S1 ===\n")

s1_data <- data.frame(
  Clock = c("DNAmGrimAge", "PhenoAge", "DunedinPACE"),
  `Original Citation` = c(
    "Lu et al. (2019)",
    "Levine et al. (2018)",
    "Belsky et al. (2022)"
  ),
  `Training Outcome` = c(
    "Plasma protein surrogates associated with mortality; composite of DNAm-based surrogates for adrenomedullin, beta-2 microglobulin, cystatin C, GDF-15, leptin, PAI-1, TIMP-1, and smoking pack-years",
    "Composite clinical biomarker representing multi-system physiological dysregulation (albumin, creatinine, glucose, CRP, lymphocyte %, MCV, RDW, alkaline phosphatase, WBC count)",
    "Rate of biological aging (pace); trained on longitudinal within-person changes in 19 biomarkers across 20 years in the Dunedin Study birth cohort"
  ),
  `Measure Type` = c(
    "Second-generation; mortality predictor",
    "Second-generation; morbidity/mortality predictor",
    "Third-generation; pace of aging"
  ),
  `Interpretation` = c(
    "Higher values = older biological age relative to chronological age; positive EAA = accelerated aging",
    "Higher values = older biological age relative to chronological age; positive EAA = accelerated aging",
    "Higher values = faster pace of biological aging; >1.0 indicates aging faster than 1 year per calendar year"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

s1_ft <- s1_data %>%
  flextable() %>%
  set_header_labels(
    Clock = "Epigenetic Clock",
    `Original Citation` = "Original Citation",
    `Training Outcome` = "Training Outcome / Biological Basis",
    `Measure Type` = "Generation / Type",
    `Interpretation` = "Interpretation"
  ) %>%
  theme_booktabs() %>%
  autofit() %>%
  width(j = 3, width = 3) %>%
  width(j = 5, width = 2) %>%
  bold(part = "header") %>%
  fontsize(size = 9, part = "all") %>%
  set_caption("Table S1. Detailed Description of Epigenetic Clocks")

doc <- doc %>%
  body_add_par("Table S1. Detailed Description of Epigenetic Clocks", style = "heading 2") %>%
  body_add_par("") %>%
  body_add_flextable(s1_ft) %>%
  body_add_par("Note: All clocks were computed from DNA methylation data collected via Illumina 450K or EPIC arrays from saliva samples in the Fragile Families and Child Wellbeing Study (FFCWS). EAA = Epigenetic Age Acceleration.", style = "Normal") %>%
  body_add_break()

cat("S1 done.\n")

# =============================================================================
# TABLE S2: Summary of All Epigenetic Clocks in FFCWS
# =============================================================================

cat("\n=== TABLE S2 ===\n")

# Compute summary stats for each clock at each wave
clock_vars <- c("z_grimAge_9", "z_phenoAge_9", "z_poam45_9",
                "z_grimAge_15", "z_phenoAge_15", "z_poam45_15")

clock_labels <- c("DNAmGrimAge", "PhenoAge", "DunedinPACE",
                   "DNAmGrimAge", "PhenoAge", "DunedinPACE")
wave_labels <- c("Age 9", "Age 9", "Age 9",
                  "Age 15", "Age 15", "Age 15")

s2_rows <- list()
for (i in seq_along(clock_vars)) {
  v <- clock_vars[i]
  vals <- data[[v]][!is.na(data[[v]])]
  s2_rows[[i]] <- data.frame(
    Wave = wave_labels[i],
    Clock = clock_labels[i],
    N = length(vals),
    Mean = round(mean(vals), 3),
    SD = round(sd(vals), 3),
    Min = round(min(vals), 3),
    Max = round(max(vals), 3),
    stringsAsFactors = FALSE
  )
}

# Also add cell composition variables
cell_vars <- list(
  list("cell_epi_9", "Epithelial proportion", "Age 9"),
  list("cell_ic_9", "Immune cell proportion", "Age 9"),
  list("cell_epi_15", "Epithelial proportion", "Age 15"),
  list("cell_ic_15", "Immune cell proportion", "Age 15")
)

for (cv in cell_vars) {
  vals <- data[[cv[[1]]]][!is.na(data[[cv[[1]]]])]
  s2_rows[[length(s2_rows) + 1]] <- data.frame(
    Wave = cv[[3]],
    Clock = cv[[2]],
    N = length(vals),
    Mean = round(mean(vals), 3),
    SD = round(sd(vals), 3),
    Min = round(min(vals), 3),
    Max = round(max(vals), 3),
    stringsAsFactors = FALSE
  )
}

# Add array type distribution
for (wave in c("9", "15")) {
  arr_var <- paste0("array_type_", wave)
  if (arr_var %in% names(data)) {
    arr_vals <- data[[arr_var]][!is.na(data[[arr_var]])]
    n_total <- length(arr_vals)
    n_epic <- sum(arr_vals == 1, na.rm = TRUE)  # assuming 1 = EPIC
    n_450k <- n_total - n_epic
    s2_rows[[length(s2_rows) + 1]] <- data.frame(
      Wave = paste0("Age ", wave),
      Clock = paste0("Array type (N 450K / N EPIC)"),
      N = n_total,
      Mean = NA,
      SD = NA,
      Min = NA,
      Max = NA,
      stringsAsFactors = FALSE
    )
  }
}

s2_data <- bind_rows(s2_rows)

s2_ft <- s2_data %>%
  flextable() %>%
  set_header_labels(
    Wave = "Wave",
    Clock = "Variable",
    N = "N",
    Mean = "Mean",
    SD = "SD",
    Min = "Min",
    Max = "Max"
  ) %>%
  theme_booktabs() %>%
  autofit() %>%
  merge_v(j = "Wave") %>%
  bold(part = "header") %>%
  fontsize(size = 9, part = "all") %>%
  set_caption("Table S2. Summary of All Epigenetic Clocks and Cell Composition Variables in FFCWS") %>%
  colformat_double(j = c("Mean", "SD", "Min", "Max"), digits = 3, na_str = "—")

doc <- doc %>%
  body_add_par("Table S2. Summary of All Epigenetic Clocks and Cell Composition Variables in FFCWS", style = "heading 2") %>%
  body_add_par("") %>%
  body_add_flextable(s2_ft) %>%
  body_add_par("Note: All clock measures are standardized z-scores. Cell proportions estimated using EpiDISH (Middleton preferred). Array type indicates DNA methylation platform (Illumina 450K vs. EPIC).", style = "Normal") %>%
  body_add_break()

cat("S2 done.\n")

# =============================================================================
# TABLE S3: Pearson Correlations between Age at Menarche and EAA
# =============================================================================

cat("\n=== TABLE S3 ===\n")

cor_results <- read.csv("output/menarche/menarche_pearson_correlations.csv")

s3_data <- cor_results %>%
  mutate(
    r = round(r, 3),
    sig = case_when(
      p < 0.001 ~ "***",
      p < 0.01 ~ "**",
      p < 0.05 ~ "*",
      p < 0.10 ~ "\u2020",
      TRUE ~ ""
    ),
    `r` = sprintf("%.3f%s", r, sig),
    `p-value` = ifelse(p < 0.001, "< 0.001", sprintf("%.3f", p))
  ) %>%
  select(Clock, Wave, N, `r`, `p-value`)

s3_ft <- s3_data %>%
  flextable() %>%
  set_header_labels(
    Clock = "Epigenetic Clock",
    Wave = "Wave",
    N = "N",
    `r` = "Pearson r",
    `p-value` = "p-value"
  ) %>%
  theme_booktabs() %>%
  autofit() %>%
  bold(part = "header") %>%
  fontsize(size = 10, part = "all") %>%
  set_caption("Table S3. Pearson Correlations between Age at Menarche and Epigenetic Age Acceleration") %>%
  add_footer_lines(paste0("\u2020 p < 0.10, * p < 0.05, ** p < 0.01, *** p < 0.001")) %>%
  add_footer_lines("Note: Negative correlation indicates earlier menarche is associated with higher EAA. EAA adjusted for cell composition (epithelial and immune cell proportions) and array type. Sample restricted to females with menarche data.")

doc <- doc %>%
  body_add_par("Table S3. Pearson Correlations between Age at Menarche and Epigenetic Age Acceleration", style = "heading 2") %>%
  body_add_par("") %>%
  body_add_flextable(s3_ft) %>%
  body_add_break()

cat("S3 done.\n")

# =============================================================================
# TABLE S4: Mean EAA by Menarche Timing Group, Ages 9 and 15
# =============================================================================

cat("\n=== TABLE S4 ===\n")

traj_data <- read.csv("output/menarche/menarche_eaa_trajectory_summary.csv")

s4_data <- traj_data %>%
  mutate(
    `Mean EAA` = sprintf("%.3f", mean_eaa),
    SE = sprintf("%.3f", se),
    `Mean (SE)` = sprintf("%.3f (%.3f)", mean_eaa, se)
  ) %>%
  select(Clock, timing, wave, `Mean (SE)`, n) %>%
  pivot_wider(
    names_from = wave,
    values_from = c(`Mean (SE)`, n),
    names_sep = "_"
  ) %>%
  rename(
    `Menarche Timing` = timing,
    `Mean (SE) Age 9` = `Mean (SE)_Age 9`,
    `Mean (SE) Age 15` = `Mean (SE)_Age 15`,
    `N Age 9` = `n_Age 9`,
    `N Age 15` = `n_Age 15`
  ) %>%
  select(Clock, `Menarche Timing`, `N Age 9`, `Mean (SE) Age 9`, `Mean (SE) Age 15`)

s4_ft <- s4_data %>%
  flextable() %>%
  theme_booktabs() %>%
  autofit() %>%
  merge_v(j = "Clock") %>%
  bold(part = "header") %>%
  fontsize(size = 9, part = "all") %>%
  set_caption("Table S4. Mean Epigenetic Age Acceleration by Menarche Timing Group, Ages 9 and 15") %>%
  add_footer_lines("Note: Within-person sample (N = 838 girls with EAA at both waves). Early = menarche age \u226412; On Time = menarche age \u226513. EAA adjusted for cell composition and array type. Standard errors in parentheses.")

doc <- doc %>%
  body_add_par("Table S4. Mean Epigenetic Age Acceleration by Menarche Timing Group, Ages 9 and 15", style = "heading 2") %>%
  body_add_par("") %>%
  body_add_flextable(s4_ft) %>%
  body_add_break()

cat("S4 done.\n")

# =============================================================================
# TABLE S5: EAA Age 9 -> Height Growth: Full Model Specifications
# =============================================================================

cat("\n=== TABLE S5 ===\n")

growth_results <- read.csv("output/height/height_growth_results.csv")

s5_data <- growth_results %>%
  mutate(
    sig = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      p.value < 0.10 ~ "\u2020",
      TRUE ~ ""
    ),
    Coefficient = sprintf("%.3f%s", estimate, sig),
    SE = sprintf("(%.3f)", std.error),
    `p-value` = ifelse(p.value < 0.001, "< 0.001", sprintf("%.3f", p.value))
  ) %>%
  select(Clock, Model, N, Coefficient, SE, `p-value`)

# Reorder models
s5_data$Model <- factor(s5_data$Model,
                         levels = c("Unadjusted", "+ Demographics",
                                    "+ Height at 9", "+ Height 9 + Mat Height"))
s5_data <- s5_data %>% arrange(Clock, Model)

s5_ft <- s5_data %>%
  flextable() %>%
  set_header_labels(
    Clock = "Epigenetic Clock",
    Model = "Model Specification",
    N = "N",
    Coefficient = "Coefficient",
    SE = "SE",
    `p-value` = "p-value"
  ) %>%
  theme_booktabs() %>%
  autofit() %>%
  merge_v(j = "Clock") %>%
  bold(part = "header") %>%
  fontsize(size = 9, part = "all") %>%
  set_caption("Table S5. Association between EAA at Age 9 and Height Growth: Full Model Specifications") %>%
  add_footer_lines(paste0("\u2020 p < 0.10, * p < 0.05, ** p < 0.01, *** p < 0.001")) %>%
  add_footer_lines("Note: Outcome = height-for-age z-score at age 15. EAA adjusted for cell composition and array type.") %>%
  add_footer_lines("Unadjusted = EAA only; + Demographics = + sex, race/ethnicity; + Height at 9 = + baseline height; + Height 9 + Mat Height = fully adjusted (primary specification). HC2 robust standard errors.")

doc <- doc %>%
  body_add_par("Table S5. Association between EAA at Age 9 and Height Growth: Full Model Specifications", style = "heading 2") %>%
  body_add_par("") %>%
  body_add_flextable(s5_ft) %>%
  body_add_break()

cat("S5 done.\n")

# =============================================================================
# TABLE S6: EAA Age 15 -> Adult Height at 22
# =============================================================================

cat("\n=== TABLE S6 ===\n")

adult_results <- read.csv("output/height/height_adult_results.csv")

s6_data <- adult_results %>%
  mutate(
    sig = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      p.value < 0.10 ~ "\u2020",
      TRUE ~ ""
    ),
    Coefficient = sprintf("%.3f%s", estimate, sig),
    SE = sprintf("(%.3f)", std.error),
    `p-value` = sprintf("%.3f", p.value)
  ) %>%
  select(Clock, Model, N, Coefficient, SE, `p-value`)

s6_data$Model <- factor(s6_data$Model,
                         levels = c("Unadjusted", "+ Demographics",
                                    "+ Maternal Height", "+ Height at 15"))
s6_data <- s6_data %>% arrange(Clock, Model)

s6_ft <- s6_data %>%
  flextable() %>%
  set_header_labels(
    Clock = "Epigenetic Clock",
    Model = "Model Specification",
    N = "N",
    Coefficient = "Coefficient",
    SE = "SE",
    `p-value` = "p-value"
  ) %>%
  theme_booktabs() %>%
  autofit() %>%
  merge_v(j = "Clock") %>%
  bold(part = "header") %>%
  fontsize(size = 9, part = "all") %>%
  set_caption("Table S6. Association between EAA at Age 15 and Adult Height at Age 22") %>%
  add_footer_lines(paste0("\u2020 p < 0.10, * p < 0.05, ** p < 0.01, *** p < 0.001")) %>%
  add_footer_lines("Note: Outcome = height in cm at age 22 (self-reported). EAA adjusted for cell composition and array type. HC2 robust standard errors. None of the associations reach statistical significance.")

doc <- doc %>%
  body_add_par("Table S6. Association between EAA at Age 15 and Adult Height at Age 22", style = "heading 2") %>%
  body_add_par("") %>%
  body_add_flextable(s6_ft) %>%
  body_add_break()

cat("S6 done.\n")

# =============================================================================
# TABLE S7: Non-Linearity Tests (Menarche x EAA)
# =============================================================================

cat("\n=== TABLE S7 ===\n")

shape_data <- read.csv("output/menarche/menarche_shape_tests.csv")

s7_data <- shape_data %>%
  mutate(
    sig = case_when(
      p < 0.001 ~ "***",
      p < 0.01 ~ "**",
      p < 0.05 ~ "*",
      p < 0.10 ~ "\u2020",
      TRUE ~ ""
    ),
    Coefficient = sprintf("%.3f%s", Estimate, sig),
    SE = sprintf("(%.3f)", SE),
    `p-value` = ifelse(p < 0.001, "< 0.001", sprintf("%.3f", p)),
    `Turning Point` = ifelse(is.na(TurningPoint), "\u2014", sprintf("%.1f", TurningPoint)),
    Wave = gsub("Wave 4", "Age 9", Wave),
    Wave = gsub("Wave 5", "Age 15", Wave)
  ) %>%
  select(Clock, Wave, Test, Coefficient, SE, `p-value`, `Turning Point`)

s7_ft <- s7_data %>%
  flextable() %>%
  set_header_labels(
    Clock = "Epigenetic Clock",
    Wave = "Wave",
    Test = "Test",
    Coefficient = "Estimate",
    SE = "SE",
    `p-value` = "p-value",
    `Turning Point` = "Turning Point (yrs)"
  ) %>%
  theme_booktabs() %>%
  autofit() %>%
  merge_v(j = c("Clock", "Wave")) %>%
  bold(part = "header") %>%
  fontsize(size = 8, part = "all") %>%
  set_caption("Table S7. Non-Linearity Tests: Epigenetic Age Acceleration and Age at Menarche") %>%
  add_footer_lines(paste0("\u2020 p < 0.10, * p < 0.05, ** p < 0.01, *** p < 0.001")) %>%
  add_footer_lines("Note: EAA adjusted for cell composition and array type. Quadratic models use menarche age centered at the sample mean. Turning point = -\u03B2\u2081/(2\u03B2\u2082) + mean(menarche age). Categorical comparisons use Middle (11\u201313) as reference group. HC2 robust standard errors.")

doc <- doc %>%
  body_add_par("Table S7. Non-Linearity Tests: Epigenetic Age Acceleration and Age at Menarche", style = "heading 2") %>%
  body_add_par("") %>%
  body_add_flextable(s7_ft) %>%
  body_add_break()

cat("S7 done.\n")

# =============================================================================
# TABLE S8: FDR Correction Across All Outcomes
# =============================================================================

cat("\n=== TABLE S8 ===\n")

fdr_data <- read.csv("output/sensitivity_analysis/FDR_Correction_All_Outcomes.csv")

s8_data <- fdr_data %>%
  mutate(
    Estimate = sprintf("%.3f", Estimate),
    `p (raw)` = ifelse(p_value < 0.001, "< 0.001", sprintf("%.4f", p_value)),
    `p (FDR)` = ifelse(p_fdr < 0.001, "< 0.001", sprintf("%.4f", p_fdr)),
    `p (Bonferroni)` = ifelse(p_bonferroni > 1, "1.000", sprintf("%.4f", p_bonferroni)),
    `Sig (raw)` = sig_raw,
    `Sig (FDR)` = sig_fdr
  ) %>%
  select(Analysis, Clock, Test, Estimate, `p (raw)`, `Sig (raw)`, `p (FDR)`, `Sig (FDR)`, `p (Bonferroni)`)

s8_ft <- s8_data %>%
  flextable() %>%
  set_header_labels(
    Analysis = "Analysis",
    Clock = "Clock",
    Test = "Test",
    Estimate = "Estimate",
    `p (raw)` = "p (raw)",
    `Sig (raw)` = "Sig",
    `p (FDR)` = "p (FDR)",
    `Sig (FDR)` = "Sig",
    `p (Bonferroni)` = "p (Bonf.)"
  ) %>%
  theme_booktabs() %>%
  autofit() %>%
  merge_v(j = "Analysis") %>%
  bold(part = "header") %>%
  fontsize(size = 8, part = "all") %>%
  set_caption("Table S8. False Discovery Rate Correction Across All Primary Outcomes") %>%
  add_footer_lines(paste0("\u2020 p < 0.10, * p < 0.05, ** p < 0.01, *** p < 0.001")) %>%
  add_footer_lines("Note: FDR correction applied using Benjamini-Hochberg method across all 9 primary tests (3 clocks \u00D7 3 outcome families). Bonferroni correction shown for comparison.")

doc <- doc %>%
  body_add_par("Table S8. False Discovery Rate Correction Across All Primary Outcomes", style = "heading 2") %>%
  body_add_par("") %>%
  body_add_flextable(s8_ft) %>%
  body_add_break()

cat("S8 done.\n")

# =============================================================================
# TABLE S9: Gender-Stratified Height Growth Analysis
# =============================================================================

cat("\n=== TABLE S9 ===\n")

gender_data <- read.csv("output/height/height_analysis_gender_stratified.csv")

s9_data <- gender_data %>%
  mutate(
    sig = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      p.value < 0.10 ~ "\u2020",
      TRUE ~ ""
    ),
    Coefficient = sprintf("%.3f%s", estimate, sig),
    SE = sprintf("(%.3f)", std.error),
    `p-value` = sprintf("%.3f", p.value)
  ) %>%
  select(Clock, Gender, N, Coefficient, SE, `p-value`)

s9_ft <- s9_data %>%
  flextable() %>%
  set_header_labels(
    Clock = "Epigenetic Clock",
    Gender = "Sex",
    N = "N",
    Coefficient = "Coefficient",
    SE = "SE",
    `p-value` = "p-value"
  ) %>%
  theme_booktabs() %>%
  autofit() %>%
  merge_v(j = "Clock") %>%
  bold(part = "header") %>%
  fontsize(size = 9, part = "all") %>%
  set_caption("Table S9. Association between EAA at Age 9 and Height Growth, Stratified by Sex") %>%
  add_footer_lines(paste0("\u2020 p < 0.10, * p < 0.05, ** p < 0.01, *** p < 0.001")) %>%
  add_footer_lines("Note: Outcome = height-for-age z-score at age 15. Fully adjusted model: EAA + height at age 9 + maternal height. EAA adjusted for cell composition and array type. HC2 robust standard errors. No significant sex differences observed.")

doc <- doc %>%
  body_add_par("Table S9. Association between EAA at Age 9 and Height Growth, Stratified by Sex", style = "heading 2") %>%
  body_add_par("") %>%
  body_add_flextable(s9_ft) %>%
  body_add_break()

cat("S9 done.\n")

# =============================================================================
# TABLE S10: SES Models WITHOUT Race Controls
# =============================================================================

cat("\n=== TABLE S10 ===\n")

# Define controls WITHOUT race
controls_9_norace <- "female + age_9 + array_type_9 + cell_epi_9 + cell_ic_9"
controls_15_norace <- "female + age_15 + array_type_15 + cell_epi_15 + cell_ic_15"

# Outcomes
outcomes_9 <- c("z_grimAge_9", "z_poam45_9", "z_phenoAge_9")
outcomes_15 <- c("z_grimAge_15", "z_poam45_15", "z_phenoAge_15")
clock_names <- c("DNAmGrimAge", "DunedinPACE", "PhenoAge")

# --- Age 9 models without race ---
s10_age9_rows <- list()
row_idx <- 1

for (i in 1:3) {
  outcome <- outcomes_9[i]
  clock_name <- clock_names[i]

  # Model 1: Income at birth
  m1 <- run_reg(paste(outcome, "~ ln_incbirth_imp +", controls_9_norace), data)

  # Model 2: Income at birth + income at 9
  m2 <- run_reg(paste(outcome, "~ ln_incbirth_imp + ln_inc9_imp +", controls_9_norace), data)

  # Model 3: Poverty ratio
  m3 <- run_reg(paste(outcome, "~ fpl_ratio_y9 +", controls_9_norace), data)

  # Extract key coefficients
  for (model_info in list(
    list(m1, "ln_incbirth_imp", "Income at Birth (M1)"),
    list(m2, "ln_incbirth_imp", "Income at Birth (M2)"),
    list(m2, "ln_inc9_imp", "Income at Age 9 (M2)"),
    list(m3, "fpl_ratio_y9", "Poverty Ratio (M3)")
  )) {
    mod <- model_info[[1]]
    var <- model_info[[2]]
    label <- model_info[[3]]

    t_mod <- tidy(mod)
    r <- t_mod %>% filter(term == var)

    if (nrow(r) > 0) {
      fc <- format_coef(r$estimate, r$std.error, r$p.value)
      s10_age9_rows[[row_idx]] <- data.frame(
        Clock = clock_name,
        Variable = label,
        Coefficient = fc$coef,
        SE = fc$se,
        N = nobs(mod),
        stringsAsFactors = FALSE
      )
      row_idx <- row_idx + 1
    }
  }
}

# --- Age 15 models without race ---
s10_age15_rows <- list()
row_idx <- 1

for (i in 1:3) {
  outcome <- outcomes_15[i]
  clock_name <- clock_names[i]

  # Model 1: Income at birth
  m1 <- run_reg(paste(outcome, "~ ln_incbirth_imp +", controls_15_norace), data)

  # Model 2: Income at birth + income at 9
  m2 <- run_reg(paste(outcome, "~ ln_incbirth_imp + ln_inc9_imp +", controls_15_norace), data)

  # Model 3: Income at birth + income at 9 + income at 15
  m3 <- run_reg(paste(outcome, "~ ln_incbirth_imp + ln_inc9_imp + ln_inc15_imp +", controls_15_norace), data)

  # Model 4: Poverty ratio
  m4 <- run_reg(paste(outcome, "~ fpl_ratio_y15 +", controls_15_norace), data)

  for (model_info in list(
    list(m1, "ln_incbirth_imp", "Income at Birth (M1)"),
    list(m2, "ln_incbirth_imp", "Income at Birth (M2)"),
    list(m2, "ln_inc9_imp", "Income at Age 9 (M2)"),
    list(m3, "ln_incbirth_imp", "Income at Birth (M3)"),
    list(m3, "ln_inc9_imp", "Income at Age 9 (M3)"),
    list(m3, "ln_inc15_imp", "Income at Age 15 (M3)"),
    list(m4, "fpl_ratio_y15", "Poverty Ratio (M4)")
  )) {
    mod <- model_info[[1]]
    var <- model_info[[2]]
    label <- model_info[[3]]

    t_mod <- tidy(mod)
    r <- t_mod %>% filter(term == var)

    if (nrow(r) > 0) {
      fc <- format_coef(r$estimate, r$std.error, r$p.value)
      s10_age15_rows[[row_idx]] <- data.frame(
        Clock = clock_name,
        Variable = label,
        Coefficient = fc$coef,
        SE = fc$se,
        N = nobs(mod),
        stringsAsFactors = FALSE
      )
      row_idx <- row_idx + 1
    }
  }
}

# Also run education models without race
s10_edu_rows <- list()
row_idx <- 1

for (wave in c("9", "15")) {
  if (wave == "9") {
    outcomes <- outcomes_9
    controls_norace <- controls_9_norace
  } else {
    outcomes <- outcomes_15
    controls_norace <- controls_15_norace
  }

  for (i in 1:3) {
    outcome <- outcomes[i]
    clock_name <- clock_names[i]

    mod <- run_reg(paste(outcome, "~ factor(educ_cat) +", controls_norace), data)
    t_mod <- tidy(mod)

    for (educ_info in list(
      list("factor(educ_cat)2", "HS/GED vs < HS"),
      list("factor(educ_cat)3", "Some College vs < HS"),
      list("factor(educ_cat)4", "College+ vs < HS")
    )) {
      var <- educ_info[[1]]
      label <- educ_info[[2]]
      r <- t_mod %>% filter(term == var)

      if (nrow(r) > 0) {
        fc <- format_coef(r$estimate, r$std.error, r$p.value)
        s10_edu_rows[[row_idx]] <- data.frame(
          Clock = clock_name,
          Wave = paste0("Age ", wave),
          Variable = label,
          Coefficient = fc$coef,
          SE = fc$se,
          N = nobs(mod),
          stringsAsFactors = FALSE
        )
        row_idx <- row_idx + 1
      }
    }
  }
}

# Create Age 9 Income table
s10_age9_data <- bind_rows(s10_age9_rows)

s10_age9_ft <- s10_age9_data %>%
  flextable() %>%
  set_header_labels(
    Clock = "Epigenetic Clock",
    Variable = "SES Variable (Model)",
    Coefficient = "Coefficient",
    SE = "SE",
    N = "N"
  ) %>%
  theme_booktabs() %>%
  autofit() %>%
  merge_v(j = "Clock") %>%
  bold(part = "header") %>%
  fontsize(size = 9, part = "all") %>%
  set_caption("Panel A: Age 9 Income Models")

# Create Age 15 Income table
s10_age15_data <- bind_rows(s10_age15_rows)

s10_age15_ft <- s10_age15_data %>%
  flextable() %>%
  set_header_labels(
    Clock = "Epigenetic Clock",
    Variable = "SES Variable (Model)",
    Coefficient = "Coefficient",
    SE = "SE",
    N = "N"
  ) %>%
  theme_booktabs() %>%
  autofit() %>%
  merge_v(j = "Clock") %>%
  bold(part = "header") %>%
  fontsize(size = 9, part = "all") %>%
  set_caption("Panel B: Age 15 Income Models")

# Create Education table
s10_edu_data <- bind_rows(s10_edu_rows)

s10_edu_ft <- s10_edu_data %>%
  flextable() %>%
  set_header_labels(
    Clock = "Epigenetic Clock",
    Wave = "Wave",
    Variable = "Education Category",
    Coefficient = "Coefficient",
    SE = "SE",
    N = "N"
  ) %>%
  theme_booktabs() %>%
  autofit() %>%
  merge_v(j = c("Clock", "Wave")) %>%
  bold(part = "header") %>%
  fontsize(size = 9, part = "all") %>%
  set_caption("Panel C: Education Models")

doc <- doc %>%
  body_add_par("Table S10. Association between Maternal SES and Epigenetic Age Acceleration, Models Without Race Controls", style = "heading 2") %>%
  body_add_par("") %>%
  body_add_par("Panel A: Income Models, Age 9", style = "heading 3") %>%
  body_add_flextable(s10_age9_ft) %>%
  body_add_par("") %>%
  body_add_par("Panel B: Income Models, Age 15", style = "heading 3") %>%
  body_add_flextable(s10_age15_ft) %>%
  body_add_par("") %>%
  body_add_par("Panel C: Education Models", style = "heading 3") %>%
  body_add_flextable(s10_edu_ft) %>%
  body_add_par("") %>%
  body_add_par(paste0("Note: Standard errors in parentheses. \u2020 p < 0.10, * p < 0.05, ** p < 0.01, *** p < 0.001. All models include controls for sex, chronological age, array type, and cell composition (epithelial and immune cell proportions). Race/ethnicity is excluded from all models in this table. Clustered standard errors at the individual level (Stata-type)."), style = "Normal")

cat("S10 done.\n")

# =============================================================================
# SAVE DOCUMENT
# =============================================================================

print(doc, target = "output/Supplementary_Tables_S1_S10.docx")

cat("\n=== SAVED: output/Supplementary_Tables_S1_S10.docx ===\n")
cat("Contents:\n")
cat("  S1  - Detailed Description of Epigenetic Clocks\n")
cat("  S2  - Summary of All Epigenetic Clocks in FFCWS\n")
cat("  S3  - Pearson Correlations (Menarche x EAA)\n")
cat("  S4  - Mean EAA by Menarche Timing Group\n")
cat("  S5  - EAA Age 9 -> Height Growth (Full Specifications)\n")
cat("  S6  - EAA Age 15 -> Adult Height at 22\n")
cat("  S7-S9 - [Reserved]\n")
cat("  S10 - SES Models Without Race Controls\n")
cat("\n=== DONE ===\n")
