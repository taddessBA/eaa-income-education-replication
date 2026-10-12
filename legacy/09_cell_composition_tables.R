# Replication code for: "Maternal Income, Education, and Children's Epigenetic Age Acceleration"
# Manuscript ID: DEM-2024-0306 | Journal: Demography
# Input: data_2026_01_05.dta (FFCWS restricted-use data)
# =============================================================================
# Cell Composition Adjusted Tables (Matching Tables 2 & 3 Format)
# For Supplementary Materials
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

format_coef <- function(estimate, se, p_value) {
  stars <- case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    TRUE ~ ""
  )
  coef_str <- sprintf("%.3f%s", estimate, stars)
  se_str <- sprintf("(%.3f)", se)
  return(list(coef = coef_str, se = se_str))
}

extract_results <- function(model, var_name) {
  tidy_mod <- tidy(model)
  row <- tidy_mod %>% filter(term == var_name)
  if (nrow(row) == 0) return(list(coef = "", se = ""))
  format_coef(row$estimate, row$std.error, row$p.value)
}

# =============================================================================
# Define Controls WITH Cell Composition
# =============================================================================

# Age 9 controls with cell composition
controls_9_cell <- "female + factor(race) + age_9 + array_type_9 + cell_epi_9 + cell_ic_9"

# Age 15 controls with cell composition
controls_15_cell <- "female + factor(race) + age_15 + array_type_15 + cell_epi_15 + cell_ic_15"

# Outcomes
outcomes_9 <- c("z_phenoAge_9", "z_poam45_9", "z_grimAge_9")
outcomes_15 <- c("z_phenoAge_15", "z_poam45_15", "z_grimAge_15")

# =============================================================================
# TABLE S-CELL-2: Age 9 with Cell Composition (matches Table 2 format)
# =============================================================================

cat("\n=== TABLE S-CELL-2: AGE 9 WITH CELL COMPOSITION ===\n")

# Run models for each outcome
table2_cell <- list()

for (outcome in outcomes_9) {
  # Model 1: Income at birth only
  m1 <- run_reg(paste(outcome, "~ ln_incbirth_imp +", controls_9_cell), data)
  
  # Model 2: Income at birth + income at wave 4 (age 9)
  m2 <- run_reg(paste(outcome, "~ ln_incbirth_imp + ln_inc9_imp +", controls_9_cell), data)
  
  # Model 3: Income at birth + income at 9 + poverty ratio
  m3 <- run_reg(paste(outcome, "~ ln_incbirth_imp + ln_inc9_imp + fpl_ratio_y9 +", controls_9_cell), data)
  
  table2_cell[[outcome]] <- list(m1 = m1, m2 = m2, m3 = m3, n = nobs(m1))
}

# Print results
cat("\nTable S-Cell-2: Age 9 Epigenetic Clocks (Cell Composition Adjusted)\n")
cat("====================================================================\n\n")

cat(sprintf("%-20s %15s %15s %15s | %15s %15s %15s | %15s %15s %15s\n",
            "", "PhenoAge", "", "", "DunedinPACE", "", "", "GrimAge", "", ""))
cat(sprintf("%-20s %15s %15s %15s | %15s %15s %15s | %15s %15s %15s\n",
            "", "(1)", "(2)", "(3)", "(1)", "(2)", "(3)", "(1)", "(2)", "(3)"))
cat(paste(rep("-", 150), collapse = ""), "\n")

# Income at Birth row
cat(sprintf("%-20s", "Income at Birth"))
for (outcome in outcomes_9) {
  r1 <- extract_results(table2_cell[[outcome]]$m1, "ln_incbirth_imp")
  r2 <- extract_results(table2_cell[[outcome]]$m2, "ln_incbirth_imp")
  r3 <- extract_results(table2_cell[[outcome]]$m3, "ln_incbirth_imp")
  cat(sprintf(" %15s %15s %15s |", r1$coef, r2$coef, r3$coef))
}
cat("\n")
cat(sprintf("%-20s", ""))
for (outcome in outcomes_9) {
  r1 <- extract_results(table2_cell[[outcome]]$m1, "ln_incbirth_imp")
  r2 <- extract_results(table2_cell[[outcome]]$m2, "ln_incbirth_imp")
  r3 <- extract_results(table2_cell[[outcome]]$m3, "ln_incbirth_imp")
  cat(sprintf(" %15s %15s %15s |", r1$se, r2$se, r3$se))
}
cat("\n\n")

# Income at Wave 4 row
cat(sprintf("%-20s", "Income at Wave 4"))
for (outcome in outcomes_9) {
  r2 <- extract_results(table2_cell[[outcome]]$m2, "ln_inc9_imp")
  r3 <- extract_results(table2_cell[[outcome]]$m3, "ln_inc9_imp")
  cat(sprintf(" %15s %15s %15s |", "", r2$coef, r3$coef))
}
cat("\n")
cat(sprintf("%-20s", ""))
for (outcome in outcomes_9) {
  r2 <- extract_results(table2_cell[[outcome]]$m2, "ln_inc9_imp")
  r3 <- extract_results(table2_cell[[outcome]]$m3, "ln_inc9_imp")
  cat(sprintf(" %15s %15s %15s |", "", r2$se, r3$se))
}
cat("\n\n")

# Poverty Ratio row
cat(sprintf("%-20s", "Poverty Ratio"))
for (outcome in outcomes_9) {
  r3 <- extract_results(table2_cell[[outcome]]$m3, "fpl_ratio_y9")
  cat(sprintf(" %15s %15s %15s |", "", "", r3$coef))
}
cat("\n")
cat(sprintf("%-20s", ""))
for (outcome in outcomes_9) {
  r3 <- extract_results(table2_cell[[outcome]]$m3, "fpl_ratio_y9")
  cat(sprintf(" %15s %15s %15s |", "", "", r3$se))
}
cat("\n\n")

cat(sprintf("%-20s %15s %15s %15s | %15s %15s %15s | %15s %15s %15s\n",
            "Controls", "X", "X", "X", "X", "X", "X", "X", "X", "X"))
cat(sprintf("%-20s %15d\n", "N", table2_cell[[outcomes_9[1]]]$n))

# =============================================================================
# TABLE S-CELL-3: Age 15 with Cell Composition (matches Table 3 format)
# =============================================================================

cat("\n\n=== TABLE S-CELL-3: AGE 15 WITH CELL COMPOSITION ===\n")

# Run models for each outcome
table3_cell <- list()

for (outcome in outcomes_15) {
  # Model 1: Income at birth only
  m1 <- run_reg(paste(outcome, "~ ln_incbirth_imp +", controls_15_cell), data)
  
  # Model 2: Income at birth + income at wave 4 (age 9)
  m2 <- run_reg(paste(outcome, "~ ln_incbirth_imp + ln_inc9_imp +", controls_15_cell), data)
  
  # Model 3: Income at birth + income at 9 + income at 15
  m3 <- run_reg(paste(outcome, "~ ln_incbirth_imp + ln_inc9_imp + ln_inc15_imp +", controls_15_cell), data)
  
  # Model 4: Full model with poverty
  m4 <- run_reg(paste(outcome, "~ ln_incbirth_imp + ln_inc9_imp + ln_inc15_imp + fpl_ratio_y15 +", controls_15_cell), data)
  
  table3_cell[[outcome]] <- list(m1 = m1, m2 = m2, m3 = m3, m4 = m4, n = nobs(m1))
}

# Print results
cat("\nTable S-Cell-3: Age 15 Epigenetic Clocks (Cell Composition Adjusted)\n")
cat("=====================================================================\n\n")

# Header
cat(sprintf("%-20s | %12s %12s %12s %12s | %12s %12s %12s %12s | %12s %12s %12s %12s\n",
            "", "PhenoAge", "", "", "", "DunedinPACE", "", "", "", "GrimAge", "", "", ""))

# Income at Birth
cat("\nIncome at Birth\n")
for (outcome in outcomes_15) {
  r1 <- extract_results(table3_cell[[outcome]]$m1, "ln_incbirth_imp")
  r2 <- extract_results(table3_cell[[outcome]]$m2, "ln_incbirth_imp")
  r3 <- extract_results(table3_cell[[outcome]]$m3, "ln_incbirth_imp")
  r4 <- extract_results(table3_cell[[outcome]]$m4, "ln_incbirth_imp")
  cat(sprintf("  %s: %s %s | %s %s | %s %s | %s %s\n", 
              outcome, r1$coef, r1$se, r2$coef, r2$se, r3$coef, r3$se, r4$coef, r4$se))
}

# Income at Wave 4
cat("\nIncome at Wave 4\n")
for (outcome in outcomes_15) {
  r2 <- extract_results(table3_cell[[outcome]]$m2, "ln_inc9_imp")
  r3 <- extract_results(table3_cell[[outcome]]$m3, "ln_inc9_imp")
  r4 <- extract_results(table3_cell[[outcome]]$m4, "ln_inc9_imp")
  cat(sprintf("  %s: — | %s %s | %s %s | %s %s\n", 
              outcome, r2$coef, r2$se, r3$coef, r3$se, r4$coef, r4$se))
}

# Income at Wave 5
cat("\nIncome at Wave 5\n")
for (outcome in outcomes_15) {
  r3 <- extract_results(table3_cell[[outcome]]$m3, "ln_inc15_imp")
  r4 <- extract_results(table3_cell[[outcome]]$m4, "ln_inc15_imp")
  cat(sprintf("  %s: — | — | %s %s | %s %s\n", 
              outcome, r3$coef, r3$se, r4$coef, r4$se))
}

# Poverty Ratio
cat("\nPoverty Ratio\n")
for (outcome in outcomes_15) {
  r4 <- extract_results(table3_cell[[outcome]]$m4, "fpl_ratio_y15")
  cat(sprintf("  %s: — | — | — | %s %s\n", 
              outcome, r4$coef, r4$se))
}

cat("\nN:", table3_cell[[outcomes_15[1]]]$n, "\n")

# =============================================================================
# CREATE WORD TABLES
# =============================================================================

cat("\n=== CREATING WORD TABLES ===\n")

if (!dir.exists("output")) dir.create("output")

# -----------------------------------------------------------------------------
# Table S-Cell-2: Age 9 (Word format)
# -----------------------------------------------------------------------------

# Build data frame for Table 2 with cell composition
build_table2_cell_data <- function(results) {
  
  # Initialize with row names
  rows <- list()
  
  # Row 1: Income at Birth (coefficients)
  row1 <- data.frame(Variable = "Income at Birth", stringsAsFactors = FALSE)
  for (outcome in names(results)) {
    r1 <- extract_results(results[[outcome]]$m1, "ln_incbirth_imp")
    r2 <- extract_results(results[[outcome]]$m2, "ln_incbirth_imp")
    r3 <- extract_results(results[[outcome]]$m3, "ln_incbirth_imp")
    row1[[paste0(outcome, "_1")]] <- r1$coef
    row1[[paste0(outcome, "_2")]] <- r2$coef
    row1[[paste0(outcome, "_3")]] <- r3$coef
  }
  rows[[1]] <- row1
  
  # Row 2: Income at Birth (SEs)
  row2 <- data.frame(Variable = "", stringsAsFactors = FALSE)
  for (outcome in names(results)) {
    r1 <- extract_results(results[[outcome]]$m1, "ln_incbirth_imp")
    r2 <- extract_results(results[[outcome]]$m2, "ln_incbirth_imp")
    r3 <- extract_results(results[[outcome]]$m3, "ln_incbirth_imp")
    row2[[paste0(outcome, "_1")]] <- r1$se
    row2[[paste0(outcome, "_2")]] <- r2$se
    row2[[paste0(outcome, "_3")]] <- r3$se
  }
  rows[[2]] <- row2
  
  # Row 3: Income at Wave 4 (coefficients)
  row3 <- data.frame(Variable = "Income at Wave 4", stringsAsFactors = FALSE)
  for (outcome in names(results)) {
    r2 <- extract_results(results[[outcome]]$m2, "ln_inc9_imp")
    r3 <- extract_results(results[[outcome]]$m3, "ln_inc9_imp")
    row3[[paste0(outcome, "_1")]] <- ""
    row3[[paste0(outcome, "_2")]] <- r2$coef
    row3[[paste0(outcome, "_3")]] <- r3$coef
  }
  rows[[3]] <- row3
  
  # Row 4: Income at Wave 4 (SEs)
  row4 <- data.frame(Variable = "", stringsAsFactors = FALSE)
  for (outcome in names(results)) {
    r2 <- extract_results(results[[outcome]]$m2, "ln_inc9_imp")
    r3 <- extract_results(results[[outcome]]$m3, "ln_inc9_imp")
    row4[[paste0(outcome, "_1")]] <- ""
    row4[[paste0(outcome, "_2")]] <- r2$se
    row4[[paste0(outcome, "_3")]] <- r3$se
  }
  rows[[4]] <- row4
  
  # Row 5: Poverty Ratio (coefficients)
  row5 <- data.frame(Variable = "Poverty Ratio", stringsAsFactors = FALSE)
  for (outcome in names(results)) {
    r3 <- extract_results(results[[outcome]]$m3, "fpl_ratio_y9")
    row5[[paste0(outcome, "_1")]] <- ""
    row5[[paste0(outcome, "_2")]] <- ""
    row5[[paste0(outcome, "_3")]] <- r3$coef
  }
  rows[[5]] <- row5
  
  # Row 6: Poverty Ratio (SEs)
  row6 <- data.frame(Variable = "", stringsAsFactors = FALSE)
  for (outcome in names(results)) {
    r3 <- extract_results(results[[outcome]]$m3, "fpl_ratio_y9")
    row6[[paste0(outcome, "_1")]] <- ""
    row6[[paste0(outcome, "_2")]] <- ""
    row6[[paste0(outcome, "_3")]] <- r3$se
  }
  rows[[6]] <- row6
  
  # Row 7: Controls
  row7 <- data.frame(Variable = "Controls", stringsAsFactors = FALSE)
  for (outcome in names(results)) {
    row7[[paste0(outcome, "_1")]] <- "X"
    row7[[paste0(outcome, "_2")]] <- "X"
    row7[[paste0(outcome, "_3")]] <- "X"
  }
  rows[[7]] <- row7
  
  # Row 8: N
  row8 <- data.frame(Variable = "N", stringsAsFactors = FALSE)
  first_outcome <- names(results)[1]
  row8[[paste0(first_outcome, "_1")]] <- as.character(results[[first_outcome]]$n)
  for (outcome in names(results)) {
    if (outcome != first_outcome) {
      row8[[paste0(outcome, "_1")]] <- ""
    }
    row8[[paste0(outcome, "_2")]] <- ""
    row8[[paste0(outcome, "_3")]] <- ""
  }
  rows[[8]] <- row8
  
  bind_rows(rows)
}

table2_cell_data <- build_table2_cell_data(table2_cell)

# Rename columns
colnames(table2_cell_data) <- c("Variable", 
                                 "Pheno_1", "Pheno_2", "Pheno_3",
                                 "PACE_1", "PACE_2", "PACE_3",
                                 "Grim_1", "Grim_2", "Grim_3")

# Create flextable
table2_cell_ft <- table2_cell_data %>%
  flextable() %>%
  set_header_labels(
    Variable = "",
    Pheno_1 = "", Pheno_2 = "", Pheno_3 = "",
    PACE_1 = "", PACE_2 = "", PACE_3 = "",
    Grim_1 = "", Grim_2 = "", Grim_3 = ""
  ) %>%
  add_header_row(
    values = c("", "PhenoAge Age 9", "DunedinPACE Age 9", "DNAmGrimAge Age 9"),
    colwidths = c(1, 3, 3, 3)
  ) %>%
  theme_booktabs() %>%
  autofit() %>%
  align(align = "center", part = "all") %>%
  align(j = 1, align = "left", part = "body") %>%
  bold(part = "header") %>%
  add_footer_lines("Note: Standard errors in parentheses. * p<0.05 ** p<0.01 *** p<0.001.") %>%
  add_footer_lines("Controls included in all models: Self-identified race/ethnicity, sex, chronological age, DNAm array chip, and cell composition (epithelial and immune cell proportions). Income in log income.") %>%
  fontsize(size = 9, part = "all")

save_as_docx(table2_cell_ft, path = "output/Table_S_Cell_Age9.docx")

# -----------------------------------------------------------------------------
# Table S-Cell-3: Age 15 (Word format)
# -----------------------------------------------------------------------------

build_table3_cell_data <- function(results) {
  
  rows <- list()
  
  # Row 1: Income at Birth (coefficients)
  row1 <- data.frame(Variable = "Income at Birth", stringsAsFactors = FALSE)
  for (outcome in names(results)) {
    r1 <- extract_results(results[[outcome]]$m1, "ln_incbirth_imp")
    r2 <- extract_results(results[[outcome]]$m2, "ln_incbirth_imp")
    r3 <- extract_results(results[[outcome]]$m3, "ln_incbirth_imp")
    r4 <- extract_results(results[[outcome]]$m4, "ln_incbirth_imp")
    row1[[paste0(outcome, "_1")]] <- r1$coef
    row1[[paste0(outcome, "_2")]] <- r2$coef
    row1[[paste0(outcome, "_3")]] <- r3$coef
    row1[[paste0(outcome, "_4")]] <- r4$coef
  }
  rows[[1]] <- row1
  
  # Row 2: Income at Birth (SEs)
  row2 <- data.frame(Variable = "", stringsAsFactors = FALSE)
  for (outcome in names(results)) {
    r1 <- extract_results(results[[outcome]]$m1, "ln_incbirth_imp")
    r2 <- extract_results(results[[outcome]]$m2, "ln_incbirth_imp")
    r3 <- extract_results(results[[outcome]]$m3, "ln_incbirth_imp")
    r4 <- extract_results(results[[outcome]]$m4, "ln_incbirth_imp")
    row2[[paste0(outcome, "_1")]] <- r1$se
    row2[[paste0(outcome, "_2")]] <- r2$se
    row2[[paste0(outcome, "_3")]] <- r3$se
    row2[[paste0(outcome, "_4")]] <- r4$se
  }
  rows[[2]] <- row2
  
  # Row 3: Income at Wave 4 (coefficients)
  row3 <- data.frame(Variable = "Income at Wave 4", stringsAsFactors = FALSE)
  for (outcome in names(results)) {
    r2 <- extract_results(results[[outcome]]$m2, "ln_inc9_imp")
    r3 <- extract_results(results[[outcome]]$m3, "ln_inc9_imp")
    r4 <- extract_results(results[[outcome]]$m4, "ln_inc9_imp")
    row3[[paste0(outcome, "_1")]] <- ""
    row3[[paste0(outcome, "_2")]] <- r2$coef
    row3[[paste0(outcome, "_3")]] <- r3$coef
    row3[[paste0(outcome, "_4")]] <- r4$coef
  }
  rows[[3]] <- row3
  
  # Row 4: Income at Wave 4 (SEs)
  row4 <- data.frame(Variable = "", stringsAsFactors = FALSE)
  for (outcome in names(results)) {
    r2 <- extract_results(results[[outcome]]$m2, "ln_inc9_imp")
    r3 <- extract_results(results[[outcome]]$m3, "ln_inc9_imp")
    r4 <- extract_results(results[[outcome]]$m4, "ln_inc9_imp")
    row4[[paste0(outcome, "_1")]] <- ""
    row4[[paste0(outcome, "_2")]] <- r2$se
    row4[[paste0(outcome, "_3")]] <- r3$se
    row4[[paste0(outcome, "_4")]] <- r4$se
  }
  rows[[4]] <- row4
  
  # Row 5: Income at Wave 5 (coefficients)
  row5 <- data.frame(Variable = "Income at Wave 5", stringsAsFactors = FALSE)
  for (outcome in names(results)) {
    r3 <- extract_results(results[[outcome]]$m3, "ln_inc15_imp")
    r4 <- extract_results(results[[outcome]]$m4, "ln_inc15_imp")
    row5[[paste0(outcome, "_1")]] <- ""
    row5[[paste0(outcome, "_2")]] <- ""
    row5[[paste0(outcome, "_3")]] <- r3$coef
    row5[[paste0(outcome, "_4")]] <- r4$coef
  }
  rows[[5]] <- row5
  
  # Row 6: Income at Wave 5 (SEs)
  row6 <- data.frame(Variable = "", stringsAsFactors = FALSE)
  for (outcome in names(results)) {
    r3 <- extract_results(results[[outcome]]$m3, "ln_inc15_imp")
    r4 <- extract_results(results[[outcome]]$m4, "ln_inc15_imp")
    row6[[paste0(outcome, "_1")]] <- ""
    row6[[paste0(outcome, "_2")]] <- ""
    row6[[paste0(outcome, "_3")]] <- r3$se
    row6[[paste0(outcome, "_4")]] <- r4$se
  }
  rows[[6]] <- row6
  
  # Row 7: Poverty Ratio (coefficients)
  row7 <- data.frame(Variable = "Poverty Ratio", stringsAsFactors = FALSE)
  for (outcome in names(results)) {
    r4 <- extract_results(results[[outcome]]$m4, "fpl_ratio_y15")
    row7[[paste0(outcome, "_1")]] <- ""
    row7[[paste0(outcome, "_2")]] <- ""
    row7[[paste0(outcome, "_3")]] <- ""
    row7[[paste0(outcome, "_4")]] <- r4$coef
  }
  rows[[7]] <- row7
  
  # Row 8: Poverty Ratio (SEs)
  row8 <- data.frame(Variable = "", stringsAsFactors = FALSE)
  for (outcome in names(results)) {
    r4 <- extract_results(results[[outcome]]$m4, "fpl_ratio_y15")
    row8[[paste0(outcome, "_1")]] <- ""
    row8[[paste0(outcome, "_2")]] <- ""
    row8[[paste0(outcome, "_3")]] <- ""
    row8[[paste0(outcome, "_4")]] <- r4$se
  }
  rows[[8]] <- row8
  
  # Row 9: Controls
  row9 <- data.frame(Variable = "Controls", stringsAsFactors = FALSE)
  for (outcome in names(results)) {
    row9[[paste0(outcome, "_1")]] <- "X"
    row9[[paste0(outcome, "_2")]] <- "X"
    row9[[paste0(outcome, "_3")]] <- "X"
    row9[[paste0(outcome, "_4")]] <- "X"
  }
  rows[[9]] <- row9
  
  # Row 10: N
  row10 <- data.frame(Variable = "N", stringsAsFactors = FALSE)
  first_outcome <- names(results)[1]
  row10[[paste0(first_outcome, "_1")]] <- as.character(results[[first_outcome]]$n)
  for (outcome in names(results)) {
    if (outcome != first_outcome) {
      row10[[paste0(outcome, "_1")]] <- ""
    }
    row10[[paste0(outcome, "_2")]] <- ""
    row10[[paste0(outcome, "_3")]] <- ""
    row10[[paste0(outcome, "_4")]] <- ""
  }
  rows[[10]] <- row10
  
  bind_rows(rows)
}

table3_cell_data <- build_table3_cell_data(table3_cell)

# Rename columns
colnames(table3_cell_data) <- c("Variable", 
                                 "Pheno_1", "Pheno_2", "Pheno_3", "Pheno_4",
                                 "PACE_1", "PACE_2", "PACE_3", "PACE_4",
                                 "Grim_1", "Grim_2", "Grim_3", "Grim_4")

# Create flextable
table3_cell_ft <- table3_cell_data %>%
  flextable() %>%
  set_header_labels(
    Variable = "",
    Pheno_1 = "", Pheno_2 = "", Pheno_3 = "", Pheno_4 = "",
    PACE_1 = "", PACE_2 = "", PACE_3 = "", PACE_4 = "",
    Grim_1 = "", Grim_2 = "", Grim_3 = "", Grim_4 = ""
  ) %>%
  add_header_row(
    values = c("", "PhenoAge Age 15", "DunedinPACE Age 15", "DNAmGrimAge Age 15"),
    colwidths = c(1, 4, 4, 4)
  ) %>%
  theme_booktabs() %>%
  autofit() %>%
  align(align = "center", part = "all") %>%
  align(j = 1, align = "left", part = "body") %>%
  bold(part = "header") %>%
  add_footer_lines("Note: Standard errors in parentheses. * p<0.05 ** p<0.01 *** p<0.001.") %>%
  add_footer_lines("Controls included in all models: Self-identified race/ethnicity, sex, chronological age, DNAm array chip, and cell composition (epithelial and immune cell proportions). Income in log income.") %>%
  fontsize(size = 9, part = "all")

save_as_docx(table3_cell_ft, path = "output/Table_S_Cell_Age15.docx")

# -----------------------------------------------------------------------------
# Save CSV versions too
# -----------------------------------------------------------------------------

write.csv(table2_cell_data, "output/Table_S_Cell_Age9.csv", row.names = FALSE)
write.csv(table3_cell_data, "output/Table_S_Cell_Age15.csv", row.names = FALSE)

cat("\n=== FILES SAVED ===\n")
cat("- output/Table_S_Cell_Age9.docx\n")
cat("- output/Table_S_Cell_Age15.docx\n")
cat("- output/Table_S_Cell_Age9.csv\n")
cat("- output/Table_S_Cell_Age15.csv\n")

# =============================================================================
# PRINT SUMMARY FOR RESPONSE MEMO
# =============================================================================

cat("\n=== SUMMARY FOR RESPONSE MEMO ===\n")
cat("\nCell composition variables used:\n")
cat("- Epithelial cell proportion (EpiDISH)\n")
cat("- Immune cell proportion (EpiDISH)\n")
cat("\nSample sizes:\n")
cat("- Age 9:", table2_cell[[outcomes_9[1]]]$n, "\n")
cat("- Age 15:", table3_cell[[outcomes_15[1]]]$n, "\n")

cat("\n=== ANALYSIS COMPLETE ===\n")

