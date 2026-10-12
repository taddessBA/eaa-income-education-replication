# Replication code for: "Maternal Income, Education, and Children's Epigenetic Age Acceleration"
# Manuscript ID: DEM-2024-0306 | Journal: Demography
# Input: data_2026_01_05.dta (FFCWS restricted-use data)
# =============================================================================
# Sensitivity Analysis: Impact of Maternal Covariates on Income-Epigenetic
#                       Aging Associations
# =============================================================================
# This script tests whether income-EAA associations are robust to controlling
# for maternal age at birth, maternal relationship at birth, and maternal
# prenatal smoking — individually and jointly.
#
# Table structure (matches manuscript format):
#   Col 1: Income Variable (base controls only)
#   Col 2: + Maternal Age at Birth
#   Col 3: + Maternal Relationship at Birth
#   Col 4: + Maternal Smoking Prenatal
#   Col 5: + All Maternal Covariates
# =============================================================================

library(tidyverse)
library(haven)
library(broom)
library(estimatr)
library(flextable)
library(officer)

# NOTE: Set working directory to the folder containing this script and data_2026_01_05.dta
# setwd("path/to/replication/")

# Load data
data <- read_dta("data_2026_01_05.dta")

cat("=== DATA LOADED ===\n")
cat("Total observations:", nrow(data), "\n")
cat("Age 9 clock sample:", sum(!is.na(data$z_grimAge_9)), "\n")
cat("Age 15 clock sample:", sum(!is.na(data$z_grimAge_15)), "\n")

# =============================================================================
# Check maternal covariates
# =============================================================================

cat("\n=== MATERNAL COVARIATE SUMMARY ===\n")
cat("mat_age_birth N:", sum(!is.na(data$mat_age_birth)), "\n")
cat("mat_rel_birth N:", sum(!is.na(data$mat_rel_birth)), "\n")
cat("mat_smoke_prenatal N:", sum(!is.na(data$mat_smoke_prenatal)), "\n")

# Sample with all three maternal covariates + clock
cat("\nAge 9 sample with all maternal covariates:",
    sum(!is.na(data$z_grimAge_9) & !is.na(data$mat_age_birth) &
        !is.na(data$mat_rel_birth) & !is.na(data$mat_smoke_prenatal)), "\n")
cat("Age 15 sample with all maternal covariates:",
    sum(!is.na(data$z_grimAge_15) & !is.na(data$mat_age_birth) &
        !is.na(data$mat_rel_birth) & !is.na(data$mat_smoke_prenatal)), "\n")

# =============================================================================
# Helper Functions
# =============================================================================

run_reg <- function(formula, data) {
  lm_robust(as.formula(formula),
            data = data,
            clusters = data$idnum,
            se_type = "stata")
}

format_coef <- function(estimate, se, p_value) {
  stars <- case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    TRUE ~ ""
  )
  coef_str <- sprintf("%.3f%s", estimate, stars)
  se_str <- sprintf("(%.3f)", se)
  paste(coef_str, se_str)
}

extract_and_format <- function(model, var_name) {
  tidy_mod <- tidy(model)
  row <- tidy_mod %>% filter(term == var_name)
  if (nrow(row) == 0) return("—")
  format_coef(row$estimate, row$std.error, row$p.value)
}

extract_n <- function(model) {
  nobs(model)
}

# =============================================================================
# Define Control Sets
# =============================================================================

# Base controls (including cell composition to match Tables 2/3)
base_9 <- "female + factor(race) + age_9 + array_type_9 + cell_epi_9 + cell_ic_9"
base_15 <- "female + factor(race) + age_15 + array_type_15 + cell_epi_15 + cell_ic_15"

# Maternal covariate additions
mat_age <- "mat_age_birth"
mat_rel <- "factor(mat_rel_birth)"
mat_smoke <- "factor(mat_smoke_prenatal)"
mat_all <- paste(mat_age, mat_rel, mat_smoke, sep = " + ")

# =============================================================================
# AGE 9 SENSITIVITY ANALYSIS
# =============================================================================

cat("\n=== AGE 9: MATERNAL COVARIATE SENSITIVITY ===\n")

outcomes_9 <- c("z_phenoAge_9", "z_poam45_9", "z_grimAge_9")
outcome_labels_9 <- c("PhenoAge", "DunedinPACE", "GrimAge")

# Income variables for Age 9 models
income_specs_9 <- list(
  list(
    name = "Model 1: Income at Birth",
    vars = "ln_incbirth_imp",
    key_var = "ln_incbirth_imp"
  ),
  list(
    name = "Model 2: Income at Birth + Age 9",
    vars = "ln_incbirth_imp + ln_inc9_imp",
    key_var = "ln_incbirth_imp",
    extra_vars = list(list(var = "ln_inc9_imp", label = "Income at Age 9"))
  ),
  list(
    name = "Model 3: Poverty Ratio",
    vars = "fpl_ratio_y9",
    key_var = "fpl_ratio_y9"
  )
)

# Collect results for Age 9
results_9 <- list()
row_num <- 1

for (i in seq_along(outcomes_9)) {
  outcome <- outcomes_9[i]
  clock_label <- outcome_labels_9[i]

  for (spec in income_specs_9) {
    # Run 5 models: base, +mat_age, +mat_rel, +mat_smoke, +all
    m_base     <- run_reg(paste(outcome, "~", spec$vars, "+", base_9), data)
    m_mat_age  <- run_reg(paste(outcome, "~", spec$vars, "+", base_9, "+", mat_age), data)
    m_mat_rel  <- run_reg(paste(outcome, "~", spec$vars, "+", base_9, "+", mat_rel), data)
    m_mat_smk  <- run_reg(paste(outcome, "~", spec$vars, "+", base_9, "+", mat_smoke), data)
    m_mat_all  <- run_reg(paste(outcome, "~", spec$vars, "+", base_9, "+", mat_all), data)

    # Extract key variable results
    results_9[[row_num]] <- data.frame(
      Age = "Age 9",
      Clock = clock_label,
      Income_Variable = spec$name,
      SES_Measure = spec$key_var,
      Base = extract_and_format(m_base, spec$key_var),
      Plus_Mat_Age = extract_and_format(m_mat_age, spec$key_var),
      Plus_Mat_Rel = extract_and_format(m_mat_rel, spec$key_var),
      Plus_Mat_Smoke = extract_and_format(m_mat_smk, spec$key_var),
      Plus_All_Maternal = extract_and_format(m_mat_all, spec$key_var),
      N_Base = extract_n(m_base),
      N_All = extract_n(m_mat_all),
      stringsAsFactors = FALSE
    )
    row_num <- row_num + 1

    # If there are extra variables to report (e.g., income at 9 in Model 2)
    if (!is.null(spec$extra_vars)) {
      for (ev in spec$extra_vars) {
        results_9[[row_num]] <- data.frame(
          Age = "Age 9",
          Clock = clock_label,
          Income_Variable = paste("  ", ev$label),
          SES_Measure = ev$var,
          Base = extract_and_format(m_base, ev$var),
          Plus_Mat_Age = extract_and_format(m_mat_age, ev$var),
          Plus_Mat_Rel = extract_and_format(m_mat_rel, ev$var),
          Plus_Mat_Smoke = extract_and_format(m_mat_smk, ev$var),
          Plus_All_Maternal = extract_and_format(m_mat_all, ev$var),
          N_Base = extract_n(m_base),
          N_All = extract_n(m_mat_all),
          stringsAsFactors = FALSE
        )
        row_num <- row_num + 1
      }
    }
  }
}

# =============================================================================
# AGE 15 SENSITIVITY ANALYSIS
# =============================================================================

cat("=== AGE 15: MATERNAL COVARIATE SENSITIVITY ===\n")

outcomes_15 <- c("z_phenoAge_15", "z_poam45_15", "z_grimAge_15")
outcome_labels_15 <- c("PhenoAge", "DunedinPACE", "GrimAge")

# Income variables for Age 15 models
income_specs_15 <- list(
  list(
    name = "Model 1: Income at Birth",
    vars = "ln_incbirth_imp",
    key_var = "ln_incbirth_imp"
  ),
  list(
    name = "Model 2: Income at Birth + Age 9",
    vars = "ln_incbirth_imp + ln_inc9_imp",
    key_var = "ln_incbirth_imp",
    extra_vars = list(list(var = "ln_inc9_imp", label = "Income at Age 9"))
  ),
  list(
    name = "Model 3: Income at Birth + Age 9 + Age 15",
    vars = "ln_incbirth_imp + ln_inc9_imp + ln_inc15_imp",
    key_var = "ln_incbirth_imp",
    extra_vars = list(
      list(var = "ln_inc9_imp", label = "Income at Age 9"),
      list(var = "ln_inc15_imp", label = "Income at Age 15")
    )
  ),
  list(
    name = "Model 4: Poverty Ratio",
    vars = "fpl_ratio_y15",
    key_var = "fpl_ratio_y15"
  )
)

# Collect results for Age 15
results_15 <- list()
row_num <- 1

for (i in seq_along(outcomes_15)) {
  outcome <- outcomes_15[i]
  clock_label <- outcome_labels_15[i]

  for (spec in income_specs_15) {
    # Run 5 models: base, +mat_age, +mat_rel, +mat_smoke, +all
    m_base     <- run_reg(paste(outcome, "~", spec$vars, "+", base_15), data)
    m_mat_age  <- run_reg(paste(outcome, "~", spec$vars, "+", base_15, "+", mat_age), data)
    m_mat_rel  <- run_reg(paste(outcome, "~", spec$vars, "+", base_15, "+", mat_rel), data)
    m_mat_smk  <- run_reg(paste(outcome, "~", spec$vars, "+", base_15, "+", mat_smoke), data)
    m_mat_all  <- run_reg(paste(outcome, "~", spec$vars, "+", base_15, "+", mat_all), data)

    # Extract key variable results
    results_15[[row_num]] <- data.frame(
      Age = "Age 15",
      Clock = clock_label,
      Income_Variable = spec$name,
      SES_Measure = spec$key_var,
      Base = extract_and_format(m_base, spec$key_var),
      Plus_Mat_Age = extract_and_format(m_mat_age, spec$key_var),
      Plus_Mat_Rel = extract_and_format(m_mat_rel, spec$key_var),
      Plus_Mat_Smoke = extract_and_format(m_mat_smk, spec$key_var),
      Plus_All_Maternal = extract_and_format(m_mat_all, spec$key_var),
      N_Base = extract_n(m_base),
      N_All = extract_n(m_mat_all),
      stringsAsFactors = FALSE
    )
    row_num <- row_num + 1

    # Extra variables
    if (!is.null(spec$extra_vars)) {
      for (ev in spec$extra_vars) {
        results_15[[row_num]] <- data.frame(
          Age = "Age 15",
          Clock = clock_label,
          Income_Variable = paste("  ", ev$label),
          SES_Measure = ev$var,
          Base = extract_and_format(m_base, ev$var),
          Plus_Mat_Age = extract_and_format(m_mat_age, ev$var),
          Plus_Mat_Rel = extract_and_format(m_mat_rel, ev$var),
          Plus_Mat_Smoke = extract_and_format(m_mat_smk, ev$var),
          Plus_All_Maternal = extract_and_format(m_mat_all, ev$var),
          N_Base = extract_n(m_base),
          N_All = extract_n(m_mat_all),
          stringsAsFactors = FALSE
        )
        row_num <- row_num + 1
      }
    }
  }
}

# =============================================================================
# Combine All Results
# =============================================================================

all_results_9 <- bind_rows(results_9)
all_results_15 <- bind_rows(results_15)
all_results <- bind_rows(all_results_9, all_results_15)

# Print results
cat("\n=== FULL RESULTS ===\n")
print(as_tibble(all_results), n = Inf, width = Inf)

# =============================================================================
# Create Word Tables
# =============================================================================

if (!dir.exists("output/sensitivity_analysis")) dir.create("output/sensitivity_analysis", recursive = TRUE)

# ------- Create a clean table for Word -------

# For Word output, remove the SES_Measure column (internal use only)
table_data_9 <- all_results_9 %>%
  select(Clock, Income_Variable, Base, Plus_Mat_Age, Plus_Mat_Rel, Plus_Mat_Smoke, Plus_All_Maternal, N_Base, N_All)

table_data_15 <- all_results_15 %>%
  select(Clock, Income_Variable, Base, Plus_Mat_Age, Plus_Mat_Rel, Plus_Mat_Smoke, Plus_All_Maternal, N_Base, N_All)

# ----- AGE 9 TABLE -----

ft_9 <- table_data_9 %>%
  flextable() %>%
  set_header_labels(
    Clock = "Epigenetic\nClock",
    Income_Variable = "Income Variable",
    Base = "Base\nControls",
    Plus_Mat_Age = "+ Maternal Age\nat Birth",
    Plus_Mat_Rel = "+ Maternal\nRelationship\nat Birth",
    Plus_Mat_Smoke = "+ Maternal\nSmoking\nPrenatal",
    Plus_All_Maternal = "+ All Maternal\nCovariates",
    N_Base = "N\n(Base)",
    N_All = "N\n(All)"
  ) %>%
  add_header_lines("Table S: Sensitivity Analysis: Impact of Maternal Covariates on Income-Epigenetic Aging Associations (Age 9)") %>%
  theme_booktabs() %>%
  autofit() %>%
  merge_v(j = "Clock") %>%
  align(align = "center", part = "all") %>%
  align(j = 1:2, align = "left", part = "body") %>%
  bold(part = "header") %>%
  add_footer_lines("Note: Coefficients with clustered standard errors in parentheses. * p<0.05, ** p<0.01, *** p<0.001.") %>%
  add_footer_lines("Base controls: sex, race/ethnicity, chronological age, DNAm array type, cell composition (epithelial and immune cell proportions).") %>%
  add_footer_lines("Maternal covariates added individually then jointly: maternal age at birth (continuous), maternal relationship status at birth (categorical), maternal prenatal smoking (categorical).") %>%
  fontsize(size = 9, part = "all") %>%
  fontsize(size = 10, i = 1, part = "header")

save_as_docx(ft_9, path = "output/sensitivity_analysis/Table_S_Maternal_Covariates_Age9.docx")

# ----- AGE 15 TABLE -----

ft_15 <- table_data_15 %>%
  flextable() %>%
  set_header_labels(
    Clock = "Epigenetic\nClock",
    Income_Variable = "Income Variable",
    Base = "Base\nControls",
    Plus_Mat_Age = "+ Maternal Age\nat Birth",
    Plus_Mat_Rel = "+ Maternal\nRelationship\nat Birth",
    Plus_Mat_Smoke = "+ Maternal\nSmoking\nPrenatal",
    Plus_All_Maternal = "+ All Maternal\nCovariates",
    N_Base = "N\n(Base)",
    N_All = "N\n(All)"
  ) %>%
  add_header_lines("Table S: Sensitivity Analysis: Impact of Maternal Covariates on Income-Epigenetic Aging Associations (Age 15)") %>%
  theme_booktabs() %>%
  autofit() %>%
  merge_v(j = "Clock") %>%
  align(align = "center", part = "all") %>%
  align(j = 1:2, align = "left", part = "body") %>%
  bold(part = "header") %>%
  add_footer_lines("Note: Coefficients with clustered standard errors in parentheses. * p<0.05, ** p<0.01, *** p<0.001.") %>%
  add_footer_lines("Base controls: sex, race/ethnicity, chronological age, DNAm array type, cell composition (epithelial and immune cell proportions).") %>%
  add_footer_lines("Maternal covariates added individually then jointly: maternal age at birth (continuous), maternal relationship status at birth (categorical), maternal prenatal smoking (categorical).") %>%
  fontsize(size = 9, part = "all") %>%
  fontsize(size = 10, i = 1, part = "header")

save_as_docx(ft_15, path = "output/sensitivity_analysis/Table_S_Maternal_Covariates_Age15.docx")

# ----- COMBINED TABLE -----

table_data_all <- all_results %>%
  select(Age, Clock, Income_Variable, Base, Plus_Mat_Age, Plus_Mat_Rel, Plus_Mat_Smoke, Plus_All_Maternal, N_Base, N_All)

ft_all <- table_data_all %>%
  flextable() %>%
  set_header_labels(
    Age = "Age",
    Clock = "Epigenetic\nClock",
    Income_Variable = "Income Variable",
    Base = "Base\nControls",
    Plus_Mat_Age = "+ Maternal Age\nat Birth",
    Plus_Mat_Rel = "+ Maternal\nRelationship\nat Birth",
    Plus_Mat_Smoke = "+ Maternal\nSmoking\nPrenatal",
    Plus_All_Maternal = "+ All Maternal\nCovariates",
    N_Base = "N\n(Base)",
    N_All = "N\n(All)"
  ) %>%
  add_header_lines("Table S: Sensitivity Analysis: Impact of Maternal Covariates on Income-Epigenetic Aging Associations") %>%
  theme_booktabs() %>%
  autofit() %>%
  merge_v(j = c("Age", "Clock")) %>%
  align(align = "center", part = "all") %>%
  align(j = 1:3, align = "left", part = "body") %>%
  bold(part = "header") %>%
  hline(i = nrow(all_results_9), part = "body", border = fp_border(width = 2)) %>%
  add_footer_lines("Note: Coefficients with clustered standard errors in parentheses. * p<0.05, ** p<0.01, *** p<0.001.") %>%
  add_footer_lines("Base controls: sex, race/ethnicity, chronological age, DNAm array type, cell composition (epithelial and immune cell proportions).") %>%
  add_footer_lines("Maternal covariates added individually then jointly: maternal age at birth (continuous), maternal relationship status at birth (categorical), maternal prenatal smoking (categorical).") %>%
  fontsize(size = 8, part = "all") %>%
  fontsize(size = 9, i = 1, part = "header")

save_as_docx(ft_all, path = "output/sensitivity_analysis/Table_S_Maternal_Covariates_Combined.docx")

# Save CSV
write.csv(all_results, "output/sensitivity_analysis/maternal_covariate_sensitivity_results.csv", row.names = FALSE)

cat("\n=== FILES SAVED ===\n")
cat("- output/sensitivity_analysis/Table_S_Maternal_Covariates_Age9.docx\n")
cat("- output/sensitivity_analysis/Table_S_Maternal_Covariates_Age15.docx\n")
cat("- output/sensitivity_analysis/Table_S_Maternal_Covariates_Combined.docx\n")
cat("- output/sensitivity_analysis/maternal_covariate_sensitivity_results.csv\n")

cat("\n=== ANALYSIS COMPLETE ===\n")
