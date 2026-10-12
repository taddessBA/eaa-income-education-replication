# Replication code for: "Maternal Income, Education, and Children's Epigenetic Age Acceleration"
# Manuscript ID: DEM-2024-0306 | Journal: Demography
# Input: data_2026_01_05.dta (FFCWS restricted-use data)
# =============================================================================
# R&R Revision: Replicate Tables 2 & 3 with Combined Sample
# Plus Cell Composition Sensitivity Analysis
# =============================================================================

library(tidyverse)
library(haven)
library(broom)
library(estimatr)
library(flextable)
library(officer)
library(openxlsx)

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
  if (nrow(row) == 0) return(list(coef = "—", se = ""))
  format_coef(row$estimate, row$std.error, row$p.value)
}

# =============================================================================
# TABLE 2: Age 9 Epigenetic Clocks (Main Models)
# =============================================================================

cat("\n=== TABLE 2: AGE 9 MODELS ===\n")

# Define controls for Age 9
controls_9 <- "female + factor(race) + age_9 + array_type_9"
cell_controls_9 <- "female + factor(race) + age_9 + array_type_9 + cell_epi_9 + cell_ic_9"

# Outcomes
outcomes_9 <- c("z_phenoAge_9", "z_poam45_9", "z_grimAge_9")
outcome_labels_9 <- c("PhenoAge Age 9", "DunedinPACE Age 9", "GrimAge Age 9")

# Run models for Age 9
table2_results <- list()

for (i in 1:length(outcomes_9)) {
  outcome <- outcomes_9[i]
  label <- outcome_labels_9[i]
  
  # Model 1: Income at birth only
  m1 <- run_reg(paste(outcome, "~ ln_incbirth_imp +", controls_9), data)
  
  # Model 2: Income at birth + income at 9
  m2 <- run_reg(paste(outcome, "~ ln_incbirth_imp + ln_inc9_imp +", controls_9), data)
  
  # Model 3: Poverty ratio only
  m3 <- run_reg(paste(outcome, "~ fpl_ratio_y9 +", controls_9), data)
  
  table2_results[[outcome]] <- list(
    label = label,
    m1 = m1, m2 = m2, m3 = m3,
    n = nobs(m1)
  )
}

# Print Table 2 results
cat("\nTable 2: Age 9 Epigenetic Clocks\n")
cat("===============================\n\n")

for (outcome in outcomes_9) {
  res <- table2_results[[outcome]]
  cat(res$label, "\n")
  
  # Income at birth
  r1 <- extract_results(res$m1, "ln_incbirth_imp")
  r2 <- extract_results(res$m2, "ln_incbirth_imp")
  cat("  Income at Birth:  M1:", r1$coef, r1$se, " | M2:", r2$coef, r2$se, "\n")
  
  # Income at 9
  r2_inc9 <- extract_results(res$m2, "ln_inc9_imp")
  cat("  Income at Age 9:  M2:", r2_inc9$coef, r2_inc9$se, "\n")
  
  # Poverty ratio
  r3 <- extract_results(res$m3, "fpl_ratio_y9")
  cat("  Poverty Ratio:    M3:", r3$coef, r3$se, "\n")
  
  cat("  N:", res$n, "\n\n")
}

# =============================================================================
# TABLE 3: Age 15 Epigenetic Clocks (Main Models)
# =============================================================================

cat("\n=== TABLE 3: AGE 15 MODELS ===\n")

# Define controls for Age 15
controls_15 <- "female + factor(race) + age_15 + array_type_15"
cell_controls_15 <- "female + factor(race) + age_15 + array_type_15 + cell_epi_15 + cell_ic_15"

# Outcomes
outcomes_15 <- c("z_phenoAge_15", "z_poam45_15", "z_grimAge_15")
outcome_labels_15 <- c("PhenoAge Age 15", "DunedinPACE Age 15", "GrimAge Age 15")

# Run models for Age 15
table3_results <- list()

for (i in 1:length(outcomes_15)) {
  outcome <- outcomes_15[i]
  label <- outcome_labels_15[i]
  
  # Model 1: Income at birth only
  m1 <- run_reg(paste(outcome, "~ ln_incbirth_imp +", controls_15), data)
  
  # Model 2: Income at birth + income at 9
  m2 <- run_reg(paste(outcome, "~ ln_incbirth_imp + ln_inc9_imp +", controls_15), data)
  
  # Model 3: Income at birth + income at 9 + income at 15
  m3 <- run_reg(paste(outcome, "~ ln_incbirth_imp + ln_inc9_imp + ln_inc15_imp +", controls_15), data)
  
  # Model 4: Poverty ratio only
  m4 <- run_reg(paste(outcome, "~ fpl_ratio_y15 +", controls_15), data)
  
  table3_results[[outcome]] <- list(
    label = label,
    m1 = m1, m2 = m2, m3 = m3, m4 = m4,
    n = nobs(m1)
  )
}

# Print Table 3 results
cat("\nTable 3: Age 15 Epigenetic Clocks\n")
cat("=================================\n\n")

for (outcome in outcomes_15) {
  res <- table3_results[[outcome]]
  cat(res$label, "\n")
  
  # Income at birth
  r1 <- extract_results(res$m1, "ln_incbirth_imp")
  r2 <- extract_results(res$m2, "ln_incbirth_imp")
  r3 <- extract_results(res$m3, "ln_incbirth_imp")
  cat("  Income at Birth:  M1:", r1$coef, r1$se, " | M2:", r2$coef, r2$se, " | M3:", r3$coef, r3$se, "\n")
  
  # Income at 9
  r2_inc9 <- extract_results(res$m2, "ln_inc9_imp")
  r3_inc9 <- extract_results(res$m3, "ln_inc9_imp")
  cat("  Income at Age 9:  M2:", r2_inc9$coef, r2_inc9$se, " | M3:", r3_inc9$coef, r3_inc9$se, "\n")
  
  # Income at 15
  r3_inc15 <- extract_results(res$m3, "ln_inc15_imp")
  cat("  Income at Age 15: M3:", r3_inc15$coef, r3_inc15$se, "\n")
  
  # Poverty ratio
  r4 <- extract_results(res$m4, "fpl_ratio_y15")
  cat("  Poverty Ratio:    M4:", r4$coef, r4$se, "\n")
  
  cat("  N:", res$n, "\n\n")
}

# =============================================================================
# CELL COMPOSITION SENSITIVITY: Tables 2 & 3 with Cell Composition
# =============================================================================

cat("\n=== CELL COMPOSITION SENSITIVITY ANALYSIS ===\n")

# Age 9 with cell composition
table2_cell_results <- list()

for (i in 1:length(outcomes_9)) {
  outcome <- outcomes_9[i]
  label <- outcome_labels_9[i]
  
  # Model 1: Income at birth only + cell composition
  m1_cell <- run_reg(paste(outcome, "~ ln_incbirth_imp +", cell_controls_9), data)
  
  # Model 2: Income at birth + income at 9 + cell composition
  m2_cell <- run_reg(paste(outcome, "~ ln_incbirth_imp + ln_inc9_imp +", cell_controls_9), data)
  
  # Model 3: Poverty ratio + cell composition
  m3_cell <- run_reg(paste(outcome, "~ fpl_ratio_y9 +", cell_controls_9), data)
  
  table2_cell_results[[outcome]] <- list(
    label = label,
    m1 = m1_cell, m2 = m2_cell, m3 = m3_cell,
    n = nobs(m1_cell)
  )
}

# Age 15 with cell composition
table3_cell_results <- list()

for (i in 1:length(outcomes_15)) {
  outcome <- outcomes_15[i]
  label <- outcome_labels_15[i]
  
  # Model 1: Income at birth only + cell composition
  m1_cell <- run_reg(paste(outcome, "~ ln_incbirth_imp +", cell_controls_15), data)
  
  # Model 2: Income at birth + income at 9 + cell composition
  m2_cell <- run_reg(paste(outcome, "~ ln_incbirth_imp + ln_inc9_imp +", cell_controls_15), data)
  
  # Model 3: Income at birth + income at 9 + income at 15 + cell composition
  m3_cell <- run_reg(paste(outcome, "~ ln_incbirth_imp + ln_inc9_imp + ln_inc15_imp +", cell_controls_15), data)
  
  # Model 4: Poverty ratio + cell composition
  m4_cell <- run_reg(paste(outcome, "~ fpl_ratio_y15 +", cell_controls_15), data)
  
  table3_cell_results[[outcome]] <- list(
    label = label,
    m1 = m1_cell, m2 = m2_cell, m3 = m3_cell, m4 = m4_cell,
    n = nobs(m1_cell)
  )
}

# Print comparison
cat("\n=== COMPARISON: WITHOUT vs WITH CELL COMPOSITION ===\n")

cat("\n--- AGE 9 ---\n")
for (outcome in outcomes_9) {
  res_without <- table2_results[[outcome]]
  res_with <- table2_cell_results[[outcome]]
  
  cat("\n", res_without$label, "\n")
  
  # Income at birth - Model 1
  r_wo <- extract_results(res_without$m1, "ln_incbirth_imp")
  r_w <- extract_results(res_with$m1, "ln_incbirth_imp")
  cat("  Income at Birth (M1): Without:", r_wo$coef, "| With Cell:", r_w$coef, "\n")
  
  # Income at birth - Model 2
  r_wo2 <- extract_results(res_without$m2, "ln_incbirth_imp")
  r_w2 <- extract_results(res_with$m2, "ln_incbirth_imp")
  cat("  Income at Birth (M2): Without:", r_wo2$coef, "| With Cell:", r_w2$coef, "\n")
  
  # Income at 9 - Model 2
  r_wo_inc9 <- extract_results(res_without$m2, "ln_inc9_imp")
  r_w_inc9 <- extract_results(res_with$m2, "ln_inc9_imp")
  cat("  Income at Age 9 (M2): Without:", r_wo_inc9$coef, "| With Cell:", r_w_inc9$coef, "\n")
  
  # Poverty - Model 3
  r_wo_pov <- extract_results(res_without$m3, "fpl_ratio_y9")
  r_w_pov <- extract_results(res_with$m3, "fpl_ratio_y9")
  cat("  Poverty Ratio (M3):   Without:", r_wo_pov$coef, "| With Cell:", r_w_pov$coef, "\n")
}

cat("\n--- AGE 15 ---\n")
for (outcome in outcomes_15) {
  res_without <- table3_results[[outcome]]
  res_with <- table3_cell_results[[outcome]]
  
  cat("\n", res_without$label, "\n")
  
  # Income at birth - Model 1
  r_wo <- extract_results(res_without$m1, "ln_incbirth_imp")
  r_w <- extract_results(res_with$m1, "ln_incbirth_imp")
  cat("  Income at Birth (M1): Without:", r_wo$coef, "| With Cell:", r_w$coef, "\n")
  
  # Income at birth - Model 3 (full)
  r_wo3 <- extract_results(res_without$m3, "ln_incbirth_imp")
  r_w3 <- extract_results(res_with$m3, "ln_incbirth_imp")
  cat("  Income at Birth (M3): Without:", r_wo3$coef, "| With Cell:", r_w3$coef, "\n")
  
  # Income at 15 - Model 3
  r_wo_inc15 <- extract_results(res_without$m3, "ln_inc15_imp")
  r_w_inc15 <- extract_results(res_with$m3, "ln_inc15_imp")
  cat("  Income at Age 15 (M3): Without:", r_wo_inc15$coef, "| With Cell:", r_w_inc15$coef, "\n")
  
  # Poverty - Model 4
  r_wo_pov <- extract_results(res_without$m4, "fpl_ratio_y15")
  r_w_pov <- extract_results(res_with$m4, "fpl_ratio_y15")
  cat("  Poverty Ratio (M4):   Without:", r_wo_pov$coef, "| With Cell:", r_w_pov$coef, "\n")
}

# =============================================================================
# CREATE FORMATTED TABLES FOR WORD
# =============================================================================

cat("\n=== CREATING WORD TABLES ===\n")

# Create output directory
if (!dir.exists("output")) dir.create("output")

# -----------------------------------------------------------------------------
# Table 2: Age 9 (formatted for Word)
# -----------------------------------------------------------------------------

create_table2_data <- function(results) {
  
  rows <- list()
  
  # Row 1: Income at Birth
  row1 <- data.frame(
    Variable = "Income at Birth (log)",
    stringsAsFactors = FALSE
  )
  for (outcome in names(results)) {
    res <- results[[outcome]]
    r1 <- extract_results(res$m1, "ln_incbirth_imp")
    r2 <- extract_results(res$m2, "ln_incbirth_imp")
    row1[[paste0(outcome, "_M1")]] <- paste(r1$coef, r1$se)
    row1[[paste0(outcome, "_M2")]] <- paste(r2$coef, r2$se)
    row1[[paste0(outcome, "_M3")]] <- "—"
  }
  rows[[1]] <- row1
  
  # Row 2: Income at Age 9
  row2 <- data.frame(
    Variable = "Income at Age 9 (log)",
    stringsAsFactors = FALSE
  )
  for (outcome in names(results)) {
    res <- results[[outcome]]
    r2 <- extract_results(res$m2, "ln_inc9_imp")
    row2[[paste0(outcome, "_M1")]] <- "—"
    row2[[paste0(outcome, "_M2")]] <- paste(r2$coef, r2$se)
    row2[[paste0(outcome, "_M3")]] <- "—"
  }
  rows[[2]] <- row2
  
  # Row 3: Poverty Ratio
  row3 <- data.frame(
    Variable = "Poverty Ratio",
    stringsAsFactors = FALSE
  )
  for (outcome in names(results)) {
    res <- results[[outcome]]
    r3 <- extract_results(res$m3, "fpl_ratio_y9")
    row3[[paste0(outcome, "_M1")]] <- "—"
    row3[[paste0(outcome, "_M2")]] <- "—"
    row3[[paste0(outcome, "_M3")]] <- paste(r3$coef, r3$se)
  }
  rows[[3]] <- row3
  
  # Combine
  table_data <- bind_rows(rows)
  
  # Add N row
  n_row <- data.frame(Variable = "N", stringsAsFactors = FALSE)
  for (outcome in names(results)) {
    res <- results[[outcome]]
    n_row[[paste0(outcome, "_M1")]] <- as.character(res$n)
    n_row[[paste0(outcome, "_M2")]] <- ""
    n_row[[paste0(outcome, "_M3")]] <- ""
  }
  table_data <- bind_rows(table_data, n_row)
  
  return(table_data)
}

table2_data <- create_table2_data(table2_results)

# Rename columns for readability
names(table2_data) <- c("Variable", 
                        "PhenoAge_M1", "PhenoAge_M2", "PhenoAge_M3",
                        "DunedinPACE_M1", "DunedinPACE_M2", "DunedinPACE_M3",
                        "GrimAge_M1", "GrimAge_M2", "GrimAge_M3")

# Create flextable
table2_ft <- table2_data %>%
  flextable() %>%
  set_header_labels(
    Variable = "",
    PhenoAge_M1 = "(1)", PhenoAge_M2 = "(2)", PhenoAge_M3 = "(3)",
    DunedinPACE_M1 = "(1)", DunedinPACE_M2 = "(2)", DunedinPACE_M3 = "(3)",
    GrimAge_M1 = "(1)", GrimAge_M2 = "(2)", GrimAge_M3 = "(3)"
  ) %>%
  add_header_row(
    values = c("", "PhenoAge Age 9", "DunedinPACE Age 9", "GrimAge Age 9"),
    colwidths = c(1, 3, 3, 3)
  ) %>%
  theme_booktabs() %>%
  autofit() %>%
  align(align = "center", part = "all") %>%
  align(j = 1, align = "left", part = "body") %>%
  bold(part = "header") %>%
  add_footer_lines("Note: Standard errors in parentheses. * p<0.05 ** p<0.01 *** p<0.001") %>%
  add_footer_lines("Controls: Self-identified race/ethnicity, sex, chronological age, array type.") %>%
  fontsize(size = 9, part = "all")

save_as_docx(table2_ft, path = "output/Table2_Age9_Income_Models.docx")

# -----------------------------------------------------------------------------
# Table 3: Age 15 (formatted for Word)
# -----------------------------------------------------------------------------

create_table3_data <- function(results) {
  
  rows <- list()
  
  # Row 1: Income at Birth
  row1 <- data.frame(Variable = "Income at Birth (log)", stringsAsFactors = FALSE)
  for (outcome in names(results)) {
    res <- results[[outcome]]
    r1 <- extract_results(res$m1, "ln_incbirth_imp")
    r2 <- extract_results(res$m2, "ln_incbirth_imp")
    r3 <- extract_results(res$m3, "ln_incbirth_imp")
    row1[[paste0(outcome, "_M1")]] <- paste(r1$coef, r1$se)
    row1[[paste0(outcome, "_M2")]] <- paste(r2$coef, r2$se)
    row1[[paste0(outcome, "_M3")]] <- paste(r3$coef, r3$se)
    row1[[paste0(outcome, "_M4")]] <- "—"
  }
  rows[[1]] <- row1
  
  # Row 2: Income at Age 9
  row2 <- data.frame(Variable = "Income at Age 9 (log)", stringsAsFactors = FALSE)
  for (outcome in names(results)) {
    res <- results[[outcome]]
    r2 <- extract_results(res$m2, "ln_inc9_imp")
    r3 <- extract_results(res$m3, "ln_inc9_imp")
    row2[[paste0(outcome, "_M1")]] <- "—"
    row2[[paste0(outcome, "_M2")]] <- paste(r2$coef, r2$se)
    row2[[paste0(outcome, "_M3")]] <- paste(r3$coef, r3$se)
    row2[[paste0(outcome, "_M4")]] <- "—"
  }
  rows[[2]] <- row2
  
  # Row 3: Income at Age 15
  row3 <- data.frame(Variable = "Income at Age 15 (log)", stringsAsFactors = FALSE)
  for (outcome in names(results)) {
    res <- results[[outcome]]
    r3 <- extract_results(res$m3, "ln_inc15_imp")
    row3[[paste0(outcome, "_M1")]] <- "—"
    row3[[paste0(outcome, "_M2")]] <- "—"
    row3[[paste0(outcome, "_M3")]] <- paste(r3$coef, r3$se)
    row3[[paste0(outcome, "_M4")]] <- "—"
  }
  rows[[3]] <- row3
  
  # Row 4: Poverty Ratio
  row4 <- data.frame(Variable = "Poverty Ratio", stringsAsFactors = FALSE)
  for (outcome in names(results)) {
    res <- results[[outcome]]
    r4 <- extract_results(res$m4, "fpl_ratio_y15")
    row4[[paste0(outcome, "_M1")]] <- "—"
    row4[[paste0(outcome, "_M2")]] <- "—"
    row4[[paste0(outcome, "_M3")]] <- "—"
    row4[[paste0(outcome, "_M4")]] <- paste(r4$coef, r4$se)
  }
  rows[[4]] <- row4
  
  # Combine
  table_data <- bind_rows(rows)
  
  # Add N row
  n_row <- data.frame(Variable = "N", stringsAsFactors = FALSE)
  for (outcome in names(results)) {
    res <- results[[outcome]]
    n_row[[paste0(outcome, "_M1")]] <- as.character(res$n)
    n_row[[paste0(outcome, "_M2")]] <- ""
    n_row[[paste0(outcome, "_M3")]] <- ""
    n_row[[paste0(outcome, "_M4")]] <- ""
  }
  table_data <- bind_rows(table_data, n_row)
  
  return(table_data)
}

table3_data <- create_table3_data(table3_results)

# Rename columns
names(table3_data) <- c("Variable", 
                        "PhenoAge_M1", "PhenoAge_M2", "PhenoAge_M3", "PhenoAge_M4",
                        "DunedinPACE_M1", "DunedinPACE_M2", "DunedinPACE_M3", "DunedinPACE_M4",
                        "GrimAge_M1", "GrimAge_M2", "GrimAge_M3", "GrimAge_M4")

# Create flextable
table3_ft <- table3_data %>%
  flextable() %>%
  set_header_labels(
    Variable = "",
    PhenoAge_M1 = "(1)", PhenoAge_M2 = "(2)", PhenoAge_M3 = "(3)", PhenoAge_M4 = "(4)",
    DunedinPACE_M1 = "(1)", DunedinPACE_M2 = "(2)", DunedinPACE_M3 = "(3)", DunedinPACE_M4 = "(4)",
    GrimAge_M1 = "(1)", GrimAge_M2 = "(2)", GrimAge_M3 = "(3)", GrimAge_M4 = "(4)"
  ) %>%
  add_header_row(
    values = c("", "PhenoAge Age 15", "DunedinPACE Age 15", "GrimAge Age 15"),
    colwidths = c(1, 4, 4, 4)
  ) %>%
  theme_booktabs() %>%
  autofit() %>%
  align(align = "center", part = "all") %>%
  align(j = 1, align = "left", part = "body") %>%
  bold(part = "header") %>%
  add_footer_lines("Note: Standard errors in parentheses. * p<0.05 ** p<0.01 *** p<0.001") %>%
  add_footer_lines("Controls: Self-identified race/ethnicity, sex, chronological age, array type.") %>%
  fontsize(size = 9, part = "all")

save_as_docx(table3_ft, path = "output/Table3_Age15_Income_Models.docx")

# -----------------------------------------------------------------------------
# Cell Composition Sensitivity Table
# -----------------------------------------------------------------------------

# Create comparison table showing with vs without cell composition
create_cell_sensitivity_comparison <- function() {
  
  rows <- list()
  row_num <- 1
  
  # AGE 9
  for (outcome in outcomes_9) {
    res_wo <- table2_results[[outcome]]
    res_w <- table2_cell_results[[outcome]]
    label <- res_wo$label
    
    # Income at birth - M1
    r_wo <- extract_results(res_wo$m1, "ln_incbirth_imp")
    r_w <- extract_results(res_w$m1, "ln_incbirth_imp")
    
    rows[[row_num]] <- data.frame(
      Age = "Age 9",
      Clock = label,
      Variable = "Income at Birth",
      Model = "M1",
      Without_Cell = paste(r_wo$coef, r_wo$se),
      With_Cell = paste(r_w$coef, r_w$se),
      stringsAsFactors = FALSE
    )
    row_num <- row_num + 1
    
    # Income at birth - M2
    r_wo2 <- extract_results(res_wo$m2, "ln_incbirth_imp")
    r_w2 <- extract_results(res_w$m2, "ln_incbirth_imp")
    
    rows[[row_num]] <- data.frame(
      Age = "Age 9",
      Clock = label,
      Variable = "Income at Birth",
      Model = "M2",
      Without_Cell = paste(r_wo2$coef, r_wo2$se),
      With_Cell = paste(r_w2$coef, r_w2$se),
      stringsAsFactors = FALSE
    )
    row_num <- row_num + 1
    
    # Income at 9 - M2
    r_wo_inc9 <- extract_results(res_wo$m2, "ln_inc9_imp")
    r_w_inc9 <- extract_results(res_w$m2, "ln_inc9_imp")
    
    rows[[row_num]] <- data.frame(
      Age = "Age 9",
      Clock = label,
      Variable = "Income at Age 9",
      Model = "M2",
      Without_Cell = paste(r_wo_inc9$coef, r_wo_inc9$se),
      With_Cell = paste(r_w_inc9$coef, r_w_inc9$se),
      stringsAsFactors = FALSE
    )
    row_num <- row_num + 1
    
    # Poverty - M3
    r_wo_pov <- extract_results(res_wo$m3, "fpl_ratio_y9")
    r_w_pov <- extract_results(res_w$m3, "fpl_ratio_y9")
    
    rows[[row_num]] <- data.frame(
      Age = "Age 9",
      Clock = label,
      Variable = "Poverty Ratio",
      Model = "M3",
      Without_Cell = paste(r_wo_pov$coef, r_wo_pov$se),
      With_Cell = paste(r_w_pov$coef, r_w_pov$se),
      stringsAsFactors = FALSE
    )
    row_num <- row_num + 1
  }
  
  # AGE 15
  for (outcome in outcomes_15) {
    res_wo <- table3_results[[outcome]]
    res_w <- table3_cell_results[[outcome]]
    label <- res_wo$label
    
    # Income at birth - M1
    r_wo <- extract_results(res_wo$m1, "ln_incbirth_imp")
    r_w <- extract_results(res_w$m1, "ln_incbirth_imp")
    
    rows[[row_num]] <- data.frame(
      Age = "Age 15",
      Clock = label,
      Variable = "Income at Birth",
      Model = "M1",
      Without_Cell = paste(r_wo$coef, r_wo$se),
      With_Cell = paste(r_w$coef, r_w$se),
      stringsAsFactors = FALSE
    )
    row_num <- row_num + 1
    
    # Income at birth - M3 (full model)
    r_wo3 <- extract_results(res_wo$m3, "ln_incbirth_imp")
    r_w3 <- extract_results(res_w$m3, "ln_incbirth_imp")
    
    rows[[row_num]] <- data.frame(
      Age = "Age 15",
      Clock = label,
      Variable = "Income at Birth",
      Model = "M3",
      Without_Cell = paste(r_wo3$coef, r_wo3$se),
      With_Cell = paste(r_w3$coef, r_w3$se),
      stringsAsFactors = FALSE
    )
    row_num <- row_num + 1
    
    # Income at 15 - M3
    r_wo_inc15 <- extract_results(res_wo$m3, "ln_inc15_imp")
    r_w_inc15 <- extract_results(res_w$m3, "ln_inc15_imp")
    
    rows[[row_num]] <- data.frame(
      Age = "Age 15",
      Clock = label,
      Variable = "Income at Age 15",
      Model = "M3",
      Without_Cell = paste(r_wo_inc15$coef, r_wo_inc15$se),
      With_Cell = paste(r_w_inc15$coef, r_w_inc15$se),
      stringsAsFactors = FALSE
    )
    row_num <- row_num + 1
    
    # Poverty - M4
    r_wo_pov <- extract_results(res_wo$m4, "fpl_ratio_y15")
    r_w_pov <- extract_results(res_w$m4, "fpl_ratio_y15")
    
    rows[[row_num]] <- data.frame(
      Age = "Age 15",
      Clock = label,
      Variable = "Poverty Ratio",
      Model = "M4",
      Without_Cell = paste(r_wo_pov$coef, r_wo_pov$se),
      With_Cell = paste(r_w_pov$coef, r_w_pov$se),
      stringsAsFactors = FALSE
    )
    row_num <- row_num + 1
  }
  
  bind_rows(rows)
}

cell_sensitivity_data <- create_cell_sensitivity_comparison()

# Create flextable
cell_sensitivity_ft <- cell_sensitivity_data %>%
  flextable() %>%
  set_header_labels(
    Age = "Age",
    Clock = "Epigenetic Clock",
    Variable = "SES Variable",
    Model = "Model",
    Without_Cell = "Without Cell Composition",
    With_Cell = "With Cell Composition"
  ) %>%
  theme_booktabs() %>%
  autofit() %>%
  merge_v(j = c("Age", "Clock")) %>%
  align(align = "center", part = "all") %>%
  align(j = 1:4, align = "left", part = "body") %>%
  bold(part = "header") %>%
  add_footer_lines("Note: Standard errors in parentheses. * p<0.05 ** p<0.01 *** p<0.001") %>%
  add_footer_lines("Cell composition adjustment includes epithelial and immune cell proportions (EpiDISH).") %>%
  add_footer_lines("M1 = Income at birth only; M2 = + concurrent income; M3 = + all income measures; M4 = Poverty ratio only") %>%
  fontsize(size = 9, part = "all")

save_as_docx(cell_sensitivity_ft, path = "output/Table_S_Cell_Composition_Sensitivity.docx")
write.csv(cell_sensitivity_data, "output/cell_composition_sensitivity_full.csv", row.names = FALSE)

cat("\n=== TABLES SAVED ===\n")
cat("- output/Table2_Age9_Income_Models.docx\n")
cat("- output/Table3_Age15_Income_Models.docx\n")
cat("- output/Table_S_Cell_Composition_Sensitivity.docx\n")
cat("- output/cell_composition_sensitivity_full.csv\n")

cat("\n=== ANALYSIS COMPLETE ===\n")
