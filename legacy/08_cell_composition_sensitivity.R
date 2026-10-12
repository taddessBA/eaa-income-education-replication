# Replication code for: "Maternal Income, Education, and Children's Epigenetic Age Acceleration"
# Manuscript ID: DEM-2024-0306 | Journal: Demography
# Input: data_2026_01_05.dta (FFCWS restricted-use data)
# =============================================================================
# Table S: Sensitivity Analysis - SES-EAA Associations
#          With and Without Cell Composition Adjustment
# =============================================================================

library(tidyverse)
library(haven)
library(broom)
library(estimatr)
library(flextable)
library(officer)

# Load data
data <- read_dta("data_2026_01_05.dta")

cat("=== DATA LOADED ===\n")
cat("Age 9 combined sample:", sum(!is.na(data$z_grimAge_9)), "\n")
cat("Age 15 combined sample:", sum(!is.na(data$z_grimAge_15)), "\n")

# =============================================================================
# Helper Functions
# =============================================================================

run_reg <- function(formula, data) {
  lm_robust(as.formula(formula),
            data = data,
            clusters = data$idnum,
            se_type = "stata")
}

extract_beta_numeric <- function(model, var_name) {
  tidy_mod <- tidy(model)
  row <- tidy_mod %>% filter(term == var_name)
  if (nrow(row) == 0) return(list(est = NA, se = NA, p = NA))
  list(est = row$estimate, se = row$std.error, p = row$p.value)
}

format_coef_stars <- function(est, p) {
  if (is.na(est)) return("")
  stars <- case_when(
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ ""
  )
  sprintf("%.3f%s", est, stars)
}

# =============================================================================
# Define Controls
# =============================================================================

controls_9 <- "female + factor(race) + age_9 + array_type_9"
controls_15 <- "female + factor(race) + age_15 + array_type_15"
controls_9_cell <- "female + factor(race) + age_9 + array_type_9 + cell_epi_9 + cell_ic_9"
controls_15_cell <- "female + factor(race) + age_15 + array_type_15 + cell_epi_15 + cell_ic_15"

# =============================================================================
# Run All Models and Collect Results
# =============================================================================

results <- list()
row_num <- 1

# Clock labels
clock_labels <- c(
  "z_grimAge_9" = "DNAmGrimAge", "z_grimAge_15" = "DNAmGrimAge",
  "z_poam45_9" = "DunedinPACE", "z_poam45_15" = "DunedinPACE",
  "z_phenoAge_9" = "PhenoAge", "z_phenoAge_15" = "PhenoAge"
)

# --- INCOME MODELS: AGE 9 ---
cat("\n=== Running Age 9 Income Models ===\n")

outcomes_9 <- c("z_grimAge_9", "z_poam45_9", "z_phenoAge_9")

for (outcome in outcomes_9) {
  clock <- clock_labels[outcome]

  # Model 1: Income at birth
  m_no  <- run_reg(paste(outcome, "~ ln_incbirth_imp +", controls_9), data)
  m_yes <- run_reg(paste(outcome, "~ ln_incbirth_imp +", controls_9_cell), data)
  r_no  <- extract_beta_numeric(m_no, "ln_incbirth_imp")
  r_yes <- extract_beta_numeric(m_yes, "ln_incbirth_imp")

  results[[row_num]] <- data.frame(
    Clock = clock, SES_Measure = "Income at Birth", Wave = "Age 9",
    Without_Cell = format_coef_stars(r_no$est, r_no$p),
    With_Cell = format_coef_stars(r_yes$est, r_yes$p),
    Difference = ifelse(!is.na(r_no$est) & !is.na(r_yes$est),
                        sprintf("%.3f", r_yes$est - r_no$est), ""),
    stringsAsFactors = FALSE
  )
  row_num <- row_num + 1

  # Model 2: Income at birth (controlling for income at 9)
  m_no  <- run_reg(paste(outcome, "~ ln_incbirth_imp + ln_inc9_imp +", controls_9), data)
  m_yes <- run_reg(paste(outcome, "~ ln_incbirth_imp + ln_inc9_imp +", controls_9_cell), data)

  for (var in c("ln_incbirth_imp", "ln_inc9_imp")) {
    r_no  <- extract_beta_numeric(m_no, var)
    r_yes <- extract_beta_numeric(m_yes, var)

    var_label <- ifelse(var == "ln_incbirth_imp", "Income at Birth", "Income at Age 9")

    results[[row_num]] <- data.frame(
      Clock = clock, SES_Measure = var_label, Wave = "Age 9",
      Without_Cell = format_coef_stars(r_no$est, r_no$p),
      With_Cell = format_coef_stars(r_yes$est, r_yes$p),
      Difference = ifelse(!is.na(r_no$est) & !is.na(r_yes$est),
                          sprintf("%.3f", r_yes$est - r_no$est), ""),
      stringsAsFactors = FALSE
    )
    row_num <- row_num + 1
  }

  # Model 3: Poverty ratio
  m_no  <- run_reg(paste(outcome, "~ fpl_ratio_y9 +", controls_9), data)
  m_yes <- run_reg(paste(outcome, "~ fpl_ratio_y9 +", controls_9_cell), data)
  r_no  <- extract_beta_numeric(m_no, "fpl_ratio_y9")
  r_yes <- extract_beta_numeric(m_yes, "fpl_ratio_y9")

  results[[row_num]] <- data.frame(
    Clock = clock, SES_Measure = "Poverty Ratio", Wave = "Age 9",
    Without_Cell = format_coef_stars(r_no$est, r_no$p),
    With_Cell = format_coef_stars(r_yes$est, r_yes$p),
    Difference = ifelse(!is.na(r_no$est) & !is.na(r_yes$est),
                        sprintf("%.3f", r_yes$est - r_no$est), ""),
    stringsAsFactors = FALSE
  )
  row_num <- row_num + 1
}

# --- INCOME MODELS: AGE 15 ---
cat("=== Running Age 15 Income Models ===\n")

outcomes_15 <- c("z_grimAge_15", "z_poam45_15", "z_phenoAge_15")

for (outcome in outcomes_15) {
  clock <- clock_labels[outcome]

  # Model 1: Income at birth only
  m_no  <- run_reg(paste(outcome, "~ ln_incbirth_imp +", controls_15), data)
  m_yes <- run_reg(paste(outcome, "~ ln_incbirth_imp +", controls_15_cell), data)
  r_no  <- extract_beta_numeric(m_no, "ln_incbirth_imp")
  r_yes <- extract_beta_numeric(m_yes, "ln_incbirth_imp")

  results[[row_num]] <- data.frame(
    Clock = clock, SES_Measure = "Income at Birth", Wave = "Age 15",
    Without_Cell = format_coef_stars(r_no$est, r_no$p),
    With_Cell = format_coef_stars(r_yes$est, r_yes$p),
    Difference = ifelse(!is.na(r_no$est) & !is.na(r_yes$est),
                        sprintf("%.3f", r_yes$est - r_no$est), ""),
    stringsAsFactors = FALSE
  )
  row_num <- row_num + 1

  # Model 3: Full income model (birth + 9 + 15)
  m_no  <- run_reg(paste(outcome, "~ ln_incbirth_imp + ln_inc9_imp + ln_inc15_imp +", controls_15), data)
  m_yes <- run_reg(paste(outcome, "~ ln_incbirth_imp + ln_inc9_imp + ln_inc15_imp +", controls_15_cell), data)

  for (var in c("ln_incbirth_imp", "ln_inc9_imp", "ln_inc15_imp")) {
    r_no  <- extract_beta_numeric(m_no, var)
    r_yes <- extract_beta_numeric(m_yes, var)

    var_label <- case_when(
      var == "ln_incbirth_imp" ~ "Income at Birth",
      var == "ln_inc9_imp" ~ "Income at Age 9",
      var == "ln_inc15_imp" ~ "Income at Age 15"
    )

    results[[row_num]] <- data.frame(
      Clock = clock, SES_Measure = var_label, Wave = "Age 15",
      Without_Cell = format_coef_stars(r_no$est, r_no$p),
      With_Cell = format_coef_stars(r_yes$est, r_yes$p),
      Difference = ifelse(!is.na(r_no$est) & !is.na(r_yes$est),
                          sprintf("%.3f", r_yes$est - r_no$est), ""),
      stringsAsFactors = FALSE
    )
    row_num <- row_num + 1
  }

  # Model 4: Poverty ratio
  m_no  <- run_reg(paste(outcome, "~ fpl_ratio_y15 +", controls_15), data)
  m_yes <- run_reg(paste(outcome, "~ fpl_ratio_y15 +", controls_15_cell), data)
  r_no  <- extract_beta_numeric(m_no, "fpl_ratio_y15")
  r_yes <- extract_beta_numeric(m_yes, "fpl_ratio_y15")

  results[[row_num]] <- data.frame(
    Clock = clock, SES_Measure = "Poverty Ratio", Wave = "Age 15",
    Without_Cell = format_coef_stars(r_no$est, r_no$p),
    With_Cell = format_coef_stars(r_yes$est, r_yes$p),
    Difference = ifelse(!is.na(r_no$est) & !is.na(r_yes$est),
                        sprintf("%.3f", r_yes$est - r_no$est), ""),
    stringsAsFactors = FALSE
  )
  row_num <- row_num + 1
}

# --- EDUCATION MODELS ---
cat("=== Running Education Models ===\n")

educ_levels <- c("factor(educ_cat)2", "factor(educ_cat)3", "factor(educ_cat)4")
educ_labels <- c("High school/GED", "Some college", "College+")

# Age 9 education
for (outcome in outcomes_9) {
  clock <- clock_labels[outcome]

  m_no  <- run_reg(paste(outcome, "~ factor(educ_cat) +", controls_9), data)
  m_yes <- run_reg(paste(outcome, "~ factor(educ_cat) +", controls_9_cell), data)

  for (j in 1:length(educ_levels)) {
    r_no  <- extract_beta_numeric(m_no, educ_levels[j])
    r_yes <- extract_beta_numeric(m_yes, educ_levels[j])

    results[[row_num]] <- data.frame(
      Clock = clock, SES_Measure = educ_labels[j], Wave = "Age 9",
      Without_Cell = format_coef_stars(r_no$est, r_no$p),
      With_Cell = format_coef_stars(r_yes$est, r_yes$p),
      Difference = ifelse(!is.na(r_no$est) & !is.na(r_yes$est),
                          sprintf("%.3f", r_yes$est - r_no$est), ""),
      stringsAsFactors = FALSE
    )
    row_num <- row_num + 1
  }
}

# Age 15 education
for (outcome in outcomes_15) {
  clock <- clock_labels[outcome]

  m_no  <- run_reg(paste(outcome, "~ factor(educ_cat) +", controls_15), data)
  m_yes <- run_reg(paste(outcome, "~ factor(educ_cat) +", controls_15_cell), data)

  for (j in 1:length(educ_levels)) {
    r_no  <- extract_beta_numeric(m_no, educ_levels[j])
    r_yes <- extract_beta_numeric(m_yes, educ_levels[j])

    results[[row_num]] <- data.frame(
      Clock = clock, SES_Measure = educ_labels[j], Wave = "Age 15",
      Without_Cell = format_coef_stars(r_no$est, r_no$p),
      With_Cell = format_coef_stars(r_yes$est, r_yes$p),
      Difference = ifelse(!is.na(r_no$est) & !is.na(r_yes$est),
                          sprintf("%.3f", r_yes$est - r_no$est), ""),
      stringsAsFactors = FALSE
    )
    row_num <- row_num + 1
  }
}

# =============================================================================
# Combine Results
# =============================================================================

all_results <- bind_rows(results)

cat("\n=== FULL SENSITIVITY COMPARISON TABLE ===\n")
print(as_tibble(all_results), n = Inf)

# =============================================================================
# Create Word Table
# =============================================================================

if (!dir.exists("output/sensitivity_analysis")) dir.create("output/sensitivity_analysis", recursive = TRUE)

# Create flextable
sensitivity_ft <- all_results %>%
  flextable() %>%
  set_header_labels(
    Clock = "Clock",
    SES_Measure = "SES Measure",
    Wave = "Wave",
    Without_Cell = "Without Cell Comp (\u03B2)",
    With_Cell = "With Cell Comp (\u03B2)",
    Difference = "Difference"
  ) %>%
  add_header_lines("Table S. Sensitivity Analysis: SES-EAA Associations With and Without Cell Composition Adjustment") %>%
  theme_booktabs() %>%
  autofit() %>%
  merge_v(j = c("Clock", "Wave")) %>%
  align(align = "center", part = "all") %>%
  align(j = 1:3, align = "left", part = "body") %>%
  bold(part = "header") %>%
  bold(i = 1, part = "header") %>%
  hline(i = which(all_results$Clock != lead(all_results$Clock, default = last(all_results$Clock))), part = "body") %>%
  add_footer_lines("Note: * p<0.05, ** p<0.01, *** p<0.001. Clustered standard errors at the individual level.") %>%
  add_footer_lines("Controls without cell composition: sex, race/ethnicity, chronological age, DNAm array type.") %>%
  add_footer_lines("Controls with cell composition: above + epithelial and immune cell proportions (EpiDISH).") %>%
  add_footer_lines("Difference = With Cell Comp coefficient minus Without Cell Comp coefficient.") %>%
  add_footer_lines("Income models: Income at Birth (Model 1), Income at Birth + Age 9 (Model 2, Age 9 only), Income at Birth + Age 9 + Age 15 (Model 3, Age 15 only). Poverty Ratio entered alone. Education: reference = Less than high school.") %>%
  fontsize(size = 9, part = "all") %>%
  fontsize(size = 10, i = 1, part = "header")

save_as_docx(sensitivity_ft, path = "output/sensitivity_analysis/Table_S_Cell_Sensitivity_Comparison.docx")
write.csv(all_results, "output/sensitivity_analysis/Table_S_Cell_Sensitivity_Comparison.csv", row.names = FALSE)

cat("\n=== FILES SAVED ===\n")
cat("- output/sensitivity_analysis/Table_S_Cell_Sensitivity_Comparison.docx\n")
cat("- output/sensitivity_analysis/Table_S_Cell_Sensitivity_Comparison.csv\n")

cat("\n=== ANALYSIS COMPLETE ===\n")
