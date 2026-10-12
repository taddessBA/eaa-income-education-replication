# Replication code for: "Maternal Income, Education, and Children's Epigenetic Age Acceleration"
# Manuscript ID: DEM-2024-0306 | Journal: Demography
# Input: data_2026_01_05.dta (FFCWS restricted-use data)
# =============================================================================
# Table 1: Descriptive Characteristics of Study Participants
# =============================================================================
# =============================================================================

library(tidyverse)
library(haven)
library(flextable)
library(officer)

# Load data
data <- read_dta("data_2026_01_05.dta")

cat("=== DATA LOADED ===
")
cat("Full sample N:", nrow(data), "
")

# =============================================================================
# Define Samples
# =============================================================================

# Full sample: all children
full <- data
N_full <- nrow(full)

# Analytic sample: children with epigenetic clock data at either wave
analytic <- data %>%
  filter(!is.na(z_grimAge_9) | !is.na(z_grimAge_15))
N_analytic <- nrow(analytic)

cat("Analytic sample N:", N_analytic, "
")
cat("  Age 9 clock:", sum(!is.na(analytic$z_grimAge_9)), "
")
cat("  Age 15 clock:", sum(!is.na(analytic$z_grimAge_15)), "
")
cat("  Both waves:", sum(!is.na(analytic$z_grimAge_9) & !is.na(analytic$z_grimAge_15)), "
")

# =============================================================================
# Helper Functions
# =============================================================================

fmt_mean_sd <- function(x) {
  x_clean <- x[!is.na(x) & x >= 0]  # exclude negative sentinel codes
  sprintf("%.2f (%.2f)", mean(x_clean), sd(x_clean))
}

fmt_n_pct <- function(n, total) {
  sprintf("%s (%.2f%%)", format(n, big.mark = ","), n / total * 100)
}

# =============================================================================
# Compute Child Descriptives (Full Sample)
# =============================================================================

cat("
=== CHILD DESCRIPTIVES (Full Sample) ===
")

# Sex
n_male <- sum(full$female == 0, na.rm = TRUE)
n_female <- sum(full$female == 1, na.rm = TRUE)

# Race/Ethnicity (race: 1=NHWhite, 2=NHBlack, 3=Hispanic, 4=Other)
n_white <- sum(full$race == 1, na.rm = TRUE)
n_black <- sum(full$race == 2, na.rm = TRUE)
n_hispanic <- sum(full$race == 3, na.rm = TRUE)
n_other_race <- sum(full$race == 4, na.rm = TRUE)
N_race <- sum(!is.na(full$race))

# Low birth weight
n_lbw <- sum(full$birt_w) == 1, na.rm = TRUE)
N_lbw <- sum(full$birth_w %in% c(0, 1))

# Age at Wave 5 (filter positive values)
age9_vals <- full$age_9[!is.na(full$age_9) & full$age_9 > 0]
age9_str <- sprintf("%.2f (%.2f)", mean(age9_vals), sd(age9_vals))

# Age at Wave 6 (filter positive values)
age15_vals <- full$age_15[!is.na(full$age_15) & full$age_15 > 0]
age15_str <- sprintf("%.2f (%.2f)", mean(age15_vals), sd(age15_vals))

cat("Male:", n_male, "| Female:", n_female, "
")
cat("White:", n_white, "| Black:", n_black, "| Hispanic:", n_hispanic, "| Other:", n_other_race, "
")
cat("Low birth weight:", n_lbw, "/", N_lbw, "
")
cat("Age Wave 4:", age9_str, "
")
cat("Age Wave 5:", age15_str, "
")

# =============================================================================
# Compute Mother Descriptives (Analytic Sample)
# =============================================================================

cat("
=== MOTHER DESCRIPTIVES (Analytic Sample) ===
")

# Mother Race/Ethnicity
m_white <- sum(analytic$racem == 1, na.rm = TRUE)
m_black <- sum(analytic$racem == 2, na.rm = TRUE)
m_hispanic <- sum(analytic$racem == 3, na.rm = TRUE)
m_other <- sum(analytic$racem == 4, na.rm = TRUE)
M_race <- sum(!is.na(analytic$racem))

# Mother age at birth
mat_age_str <- sprintf("%.2f (%.2f)",
                        mean(analytic$mat_age_birth, na.rm = TRUE),
                        sd(analytic$mat_age_birth, na.rm = TRUE))

# Education (educ_cat: 1=<HS, 2=HS/GED, 3=Some college, 4=College+)
m_edu1 <- sum(analytic$educ_cat == 1, na.rm = TRUE)
m_edu2 <- sum(analytic$educ_cat == 2, na.rm = TRUE)
m_edu3 <- sum(analytic$educ_cat == 3, na.rm = TRUE)
m_edu4 <- sum(analytic$educ_cat == 4, na.rm = TRUE)
M_edu <- sum(!is.na(analytic$educ_cat))

# Income at child's birth (raw, positive values only)
inc_birth_vals <- analytic$incbirth[!is.na(analytic$incbirth) & analytic$incbirth >= 0]
inc_birth_str <- sprintf("%s (%s)",
                          format(round(mean(inc_birth_vals), 2), big.mark = ","),
                          format(round(sd(inc_birth_vals), 2), big.mark = ","))

# Income at Wave 5 (age 9)
inc9_vals <- analytic$inc9[!is.na(analytic$inc9) & analytic$inc9 >= 0]
inc9_str <- sprintf("%s (%s)",
                     format(round(mean(inc9_vals), 2), big.mark = ","),
                     format(round(sd(inc9_vals), 2), big.mark = ","))

# Income at Wave 6 (age 15)
inc15_vals <- analytic$inc15[!is.na(analytic$inc15) & analytic$inc15 >= 0]
inc15_str <- sprintf("%s (%s)",
                      format(round(mean(inc15_vals), 2), big.mark = ","),
                      format(round(sd(inc15_vals), 2), big.mark = ","))

# Poverty ratio
fpl_vals <- analytic$fpl_ratio[!is.na(analytic$fpl_ratio) & analytic$fpl_ratio >= 0]
fpl_str <- sprintf("%.2f (%.2f)", mean(fpl_vals), sd(fpl_vals))

cat("Mother Race: White", m_white, "| Black", m_black, "| Hispanic", m_hispanic, "| Other", m_other, "
")
cat("Mother age at birth:", mat_age_str, "
")
cat("Education: <HS", m_edu1, "| HS", m_edu2, "| Some college", m_edu3, "| College+", m_edu4, "
")
cat("Income at birth:", inc_birth_str, "
")
cat("Income at Wave 5:", inc9_str, "
")
cat("Income at Wave 6:", inc15_str, "
")
cat("Poverty ratio:", fpl_str, "
")

# =============================================================================
# Build Table 1 Data Frame
# =============================================================================

table1 <- data.frame(
  Variable = c(
    "Sample Size",
    "",
    "Sex",
    "  Male",
    "  Female",
    "Race/Ethnicity",
    "  White, non-Hispanic",
    "  Black, non-Hispanic",
    "  Hispanic",
    "  Other",
    "Low birth weight",
    "Age at Wave 5",
    "Age at Wave 6",
    "Age at Child's Birth",
    "Education",
    "  Less than high school",
    "  High school diploma/GED",
    "  Some college/Technical",
    "  College degree or higher",
    "Maternal Income",
    "  Income at child's birth",
    "  Income at Wave 5",
    "  Income at Wave 6",
    "Maternal Poverty Ratio"
  ),
  Child = c(
    format(N_full, big.mark = ","),
    "Mean (SD) or n (%)",
    "",
    fmt_n_pct(n_male, N_full),
    fmt_n_pct(n_female, N_full),
    "",
    fmt_n_pct(n_white, N_race),
    fmt_n_pct(n_black, N_race),
    fmt_n_pct(n_hispanic, N_race),
    fmt_n_pct(n_other_race, N_race),
    fmt_n_pct(n_lbw, N_lbw),
    age9_str,
    age15_str,
    "
