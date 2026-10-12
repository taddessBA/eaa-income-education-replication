# Replication code for: "Maternal Income, Education, and Children's Epigenetic Age Acceleration"
# Manuscript ID: DEM-2024-0306 | Journal: Demography
# Input: data_2026_01_05.dta (FFCWS restricted-use data)
# =============================================================================
# Height Analysis: Growth Story + Adult Height
# 1. EAA at Age 9 → Standardized Self-Reported Height at 15 (controlling for height at 9) - GROWTH
# 2. EAA at Age 15 → Height at 22 - ADULT HEIGHT ATTAINMENT
# =============================================================================

library(tidyverse)
library(haven)
library(estimatr)
library(broom)
library(ggplot2)
library(patchwork)
library(flextable)
library(officer)

# Load data
data <- read_dta("data_2026_01_05.dta")

cat("=== DATA LOADED ===\n")

# =============================================================================
# Check Height Variables
# =============================================================================

cat("\n=== HEIGHT VARIABLES ===\n")

# Check for height at age 22 (feet and inches)
cat("\nHeight at age 22 variables:\n")
cat("p7l7a (feet):\n")
print(table(data$p7l7a, useNA = "ifany"))
cat("\np7l7b (inches):\n")
print(table(data$p7l7b, useNA = "ifany"))

# Check height at age 9 (wave 4)
cat("\nHeight at age 9 variables available:\n")
height_vars <- names(data)[grep("ht|height|haz", names(data), ignore.case = TRUE)]
print(height_vars)

# =============================================================================
# Clean and Construct Height Variables
# =============================================================================

data <- data %>%
  mutate(
    # Clean maternal height (ch5mhtcm = BioMom's height)
    mat_height_clean = case_when(
      ch5mhtcm < 0 ~ NA_real_,
      ch5mhtcm < 140 ~ NA_real_,
      ch5mhtcm > 200 ~ NA_real_,
      TRUE ~ ch5mhtcm
    ),
    
    # Height at age 22 (construct from feet + inches)
    # Filter out implausible feet values (should be 4-6 feet)
    height_22_ft = case_when(
      p7l7a < 4 ~ NA_real_,
      p7l7a > 6 ~ NA_real_,  # Exclude 7, 10 feet (errors)
      TRUE ~ as.numeric(p7l7a)
    ),
    height_22_in = case_when(
      p7l7b < 0 ~ NA_real_,
      p7l7b > 11 ~ NA_real_,  # Inches should be 0-11
      TRUE ~ as.numeric(p7l7b)
    ),
    # Total height in cm at age 22
    height_22_cm = (height_22_ft * 12 + height_22_in) * 2.54,
    
    # Clean child self-reported height at age 15 (ck6chtcm)
    child_height_15 = case_when(
      ck6chtcm < 0 ~ NA_real_,    # Missing codes (-9 = not in wave, -3 = missing)
      ck6chtcm < 120 ~ NA_real_,  # Too short (implausible for age 15)
      ck6chtcm > 210 ~ NA_real_,  # Too tall (implausible)
      TRUE ~ ck6chtcm
    )
  )

# Standardize self-reported height at 15 (z-scored, mean=0, SD=1)
data$z_height_15 <- as.numeric(scale(data$child_height_15))

# Check height at age 22
cat("\n=== HEIGHT AT AGE 22 (constructed) ===\n")
cat("N:", sum(!is.na(data$height_22_cm)), "\n")
cat("Mean:", round(mean(data$height_22_cm, na.rm = TRUE), 2), "cm\n")
cat("SD:", round(sd(data$height_22_cm, na.rm = TRUE), 2), "cm\n")
cat("Range:", round(min(data$height_22_cm, na.rm = TRUE), 1), "-", 
    round(max(data$height_22_cm, na.rm = TRUE), 1), "cm\n")

# Check for height at age 9 (wave 4)
cat("\n=== HEIGHT VARIABLES AT AGE 9 ===\n")

# Check ch5haz (height-for-age z-score at wave 5/age 9)
if ("ch5haz" %in% names(data)) {
  cat("ch5haz (HAZ at age 9):\n")
  cat("  N:", sum(!is.na(data$ch5haz) & data$ch5haz > -10), "\n")
  cat("  Mean:", round(mean(data$ch5haz[data$ch5haz > -10], na.rm = TRUE), 3), "\n")
}

# Check ch5chtcm (height in cm at wave 5)
if ("ch5chtcm" %in% names(data)) {
  valid_h9 <- data$ch5chtcm[data$ch5chtcm > 0]
  cat("\nch5chtcm (height cm at age 9):\n")
  cat("  N:", length(valid_h9), "\n")
  cat("  Mean:", round(mean(valid_h9, na.rm = TRUE), 2), "cm\n")
}

# Clean height at age 9
data <- data %>%
  mutate(
    # Height at age 9 (cm)
    child_height_9 = case_when(
      ch5chtcm < 0 ~ NA_real_,
      ch5chtcm < 100 ~ NA_real_,
      ch5chtcm > 180 ~ NA_real_,
      TRUE ~ ch5chtcm
    ),
    # HAZ at age 9
    haz_9 = case_when(
      ch5haz < -10 ~ NA_real_,
      TRUE ~ ch5haz
    )
  )

cat("\n=== CLEANED HEIGHT AT AGE 9 ===\n")
cat("child_height_9 (cm):\n")
cat("  N:", sum(!is.na(data$child_height_9)), "\n")
cat("  Mean:", round(mean(data$child_height_9, na.rm = TRUE), 2), "cm\n")

# =============================================================================
# Cell Composition Adjustment
# =============================================================================

cat("\n=== CELL COMPOSITION ADJUSTMENT (+ Array Type) ===\n")

get_cell_adjusted_data <- function(data, clock_var, cell_epi_var, cell_ic_var, array_var) {
  complete_data <- data %>%
    filter(!is.na(!!sym(clock_var)) & !is.na(!!sym(cell_epi_var)) & 
             !is.na(!!sym(cell_ic_var)) & !is.na(!!sym(array_var)))
  
  # Adjust for cell composition AND array type
  formula <- paste(clock_var, "~", cell_epi_var, "+", cell_ic_var, "+ factor(", array_var, ")")
  model <- lm(as.formula(formula), data = complete_data)
  complete_data$eaa_adjusted <- residuals(model)
  
  return(complete_data)
}

# EAA at age 9 (for height at 15 analysis)
grim_adj_9 <- get_cell_adjusted_data(data, "z_grimAge_9", "cell_epi_9", "cell_ic_9", "array_type_9")
pheno_adj_9 <- get_cell_adjusted_data(data, "z_phenoAge_9", "cell_epi_9", "cell_ic_9", "array_type_9")
pace_adj_9 <- get_cell_adjusted_data(data, "z_poam45_9", "cell_epi_9", "cell_ic_9", "array_type_9")

# EAA at age 15 (for height at 22 analysis)
grim_adj_15 <- get_cell_adjusted_data(data, "z_grimAge_15", "cell_epi_15", "cell_ic_15", "array_type_15")
pheno_adj_15 <- get_cell_adjusted_data(data, "z_phenoAge_15", "cell_epi_15", "cell_ic_15", "array_type_15")
pace_adj_15 <- get_cell_adjusted_data(data, "z_poam45_15", "cell_epi_15", "cell_ic_15", "array_type_15")

cat("EAA samples at Age 9:", nrow(grim_adj_9), "\n")
cat("EAA samples at Age 15:", nrow(grim_adj_15), "\n")

# Merge height variables
merge_height_vars <- function(df, data) {
  df %>%
    left_join(data %>% select(idnum, z_height_15, child_height_15, child_height_9, haz_9,
                              mat_height_clean, height_22_cm, female, race),
              by = "idnum", suffix = c("", ".y")) %>%
    select(-ends_with(".y"))
}

grim_adj_9 <- merge_height_vars(grim_adj_9, data)
pheno_adj_9 <- merge_height_vars(pheno_adj_9, data)
pace_adj_9 <- merge_height_vars(pace_adj_9, data)

grim_adj_15 <- merge_height_vars(grim_adj_15, data)
pheno_adj_15 <- merge_height_vars(pheno_adj_15, data)
pace_adj_15 <- merge_height_vars(pace_adj_15, data)

# =============================================================================
# ANALYSIS 1: GROWTH STORY
# EAA at Age 9 → Height at 15, controlling for Height at 9
# =============================================================================

cat("\n")
cat("=============================================================================\n")
cat("ANALYSIS 1: GROWTH STORY\n")
cat("EAA at Age 9 → Height at 15, controlling for Height at 9\n")
cat("=============================================================================\n")

run_growth_analysis <- function(df, clock_label) {
  
  cat("\n=============================\n")
  cat("Clock:", clock_label, "\n")
  
  # Filter to complete cases
  d_basic <- df %>%
    filter(!is.na(eaa_adjusted) & !is.na(z_height_15) & !is.na(female) & !is.na(race))
  
  d_height9 <- d_basic %>%
    filter(!is.na(child_height_9))
  
  d_full <- d_height9 %>%
    filter(!is.na(mat_height_clean))
  
  cat("N (no height control):", nrow(d_basic), "\n")
  cat("N (with height at 9):", nrow(d_height9), "\n")
  cat("N (with height at 9 + maternal height):", nrow(d_full), "\n")
  
  # Model 1: Unadjusted
  m1 <- lm_robust(z_height_15 ~ eaa_adjusted, data = d_basic, se_type = "HC2")
  
  # Model 2: + Demographics
  m2 <- lm_robust(z_height_15 ~ eaa_adjusted + female + factor(race), 
                  data = d_basic, se_type = "HC2")
  
  # Model 3: + Height at age 9 (GROWTH MODEL)
  m3 <- lm_robust(z_height_15 ~ eaa_adjusted + female + factor(race) + child_height_9, 
                  data = d_height9, se_type = "HC2")
  
  # Model 4: + Height at 9 + Maternal height
  m4 <- lm_robust(z_height_15 ~ eaa_adjusted + female + factor(race) + child_height_9 + mat_height_clean, 
                  data = d_full, se_type = "HC2")
  
  # Print results
  cat("\nModel 1 (Unadjusted):\n")
  print(tidy(m1) %>% filter(term == "eaa_adjusted") %>% 
          select(term, estimate, std.error, p.value))
  
  cat("\nModel 2 (+ Sex, Race):\n")
  print(tidy(m2) %>% filter(term == "eaa_adjusted") %>% 
          select(term, estimate, std.error, p.value))
  
  cat("\nModel 3 (+ Height at Age 9) - GROWTH MODEL:\n")
  print(tidy(m3) %>% filter(term == "eaa_adjusted") %>% 
          select(term, estimate, std.error, p.value))
  
  cat("\nModel 4 (+ Height at 9 + Maternal Height):\n")
  print(tidy(m4) %>% filter(term == "eaa_adjusted") %>% 
          select(term, estimate, std.error, p.value))
  
  list(
    clock = clock_label,
    m1 = m1, m2 = m2, m3 = m3, m4 = m4,
    n_basic = nrow(d_basic),
    n_height9 = nrow(d_height9),
    n_full = nrow(d_full)
  )
}

growth_grim <- run_growth_analysis(grim_adj_9, "GrimAge")
growth_pheno <- run_growth_analysis(pheno_adj_9, "PhenoAge")
growth_pace <- run_growth_analysis(pace_adj_9, "DunedinPACE")

# =============================================================================
# ANALYSIS 2: ADULT HEIGHT
# EAA at Age 15 → Height at 22
# =============================================================================

cat("\n")
cat("=============================================================================\n")
cat("ANALYSIS 2: ADULT HEIGHT ATTAINMENT\n")
cat("EAA at Age 15 → Height at 22\n")
cat("=============================================================================\n")

run_adult_height_analysis <- function(df, clock_label) {
  
  cat("\n=============================\n")
  cat("Clock:", clock_label, "\n")
  
  # Filter to complete cases
  d_basic <- df %>%
    filter(!is.na(eaa_adjusted) & !is.na(height_22_cm) & !is.na(female) & !is.na(race))
  
  d_mat <- d_basic %>%
    filter(!is.na(mat_height_clean))
  
  d_height15 <- d_basic %>%
    filter(!is.na(child_height_15))
  
  cat("N (basic):", nrow(d_basic), "\n")
  cat("N (with maternal height):", nrow(d_mat), "\n")
  cat("N (with height at 15):", nrow(d_height15), "\n")
  
  if (nrow(d_basic) < 30) {
    cat("WARNING: Sample size too small for reliable analysis\n")
    return(NULL)
  }
  
  # Model 1: Unadjusted
  m1 <- lm_robust(height_22_cm ~ eaa_adjusted, data = d_basic, se_type = "HC2")
  
  # Model 2: + Demographics
  m2 <- lm_robust(height_22_cm ~ eaa_adjusted + female + factor(race), 
                  data = d_basic, se_type = "HC2")
  
  # Model 3: + Maternal height
  m3 <- NULL
  if (nrow(d_mat) >= 30) {
    m3 <- lm_robust(height_22_cm ~ eaa_adjusted + female + factor(race) + mat_height_clean, 
                    data = d_mat, se_type = "HC2")
  }
  
  # Model 4: + Height at 15 (growth from 15 to 22)
  m4 <- NULL
  if (nrow(d_height15) >= 30) {
    m4 <- lm_robust(height_22_cm ~ eaa_adjusted + female + factor(race) + child_height_15, 
                    data = d_height15, se_type = "HC2")
  }
  
  # Print results
  cat("\nModel 1 (Unadjusted):\n")
  print(tidy(m1) %>% filter(term == "eaa_adjusted") %>% 
          select(term, estimate, std.error, p.value))
  
  cat("\nModel 2 (+ Sex, Race):\n")
  print(tidy(m2) %>% filter(term == "eaa_adjusted") %>% 
          select(term, estimate, std.error, p.value))
  
  if (!is.null(m3)) {
    cat("\nModel 3 (+ Maternal Height):\n")
    print(tidy(m3) %>% filter(term == "eaa_adjusted") %>% 
            select(term, estimate, std.error, p.value))
  }
  
  if (!is.null(m4)) {
    cat("\nModel 4 (+ Height at Age 15):\n")
    print(tidy(m4) %>% filter(term == "eaa_adjusted") %>% 
            select(term, estimate, std.error, p.value))
  }
  
  list(
    clock = clock_label,
    m1 = m1, m2 = m2, m3 = m3, m4 = m4,
    n_basic = nrow(d_basic),
    n_mat = nrow(d_mat),
    n_height15 = nrow(d_height15)
  )
}

adult_grim <- run_adult_height_analysis(grim_adj_15, "GrimAge")
adult_pheno <- run_adult_height_analysis(pheno_adj_15, "PhenoAge")
adult_pace <- run_adult_height_analysis(pace_adj_15, "DunedinPACE")

# =============================================================================
# Compile Results
# =============================================================================

cat("\n")
cat("=============================================================================\n")
cat("SUMMARY: GROWTH ANALYSIS (EAA Age 9 → Height at 15, controlling for Height at 9)\n")
cat("=============================================================================\n")

extract_growth_results <- function(results, model_name, model_obj) {
  if (is.null(model_obj)) return(NULL)
  tidy(model_obj) %>% 
    filter(term == "eaa_adjusted") %>%
    mutate(Clock = results$clock, Model = model_name, N = nobs(model_obj)) %>%
    select(Clock, Model, N, estimate, std.error, p.value)
}

growth_results <- bind_rows(
  extract_growth_results(growth_grim, "Unadjusted", growth_grim$m1),
  extract_growth_results(growth_grim, "+ Demographics", growth_grim$m2),
  extract_growth_results(growth_grim, "+ Height at 9", growth_grim$m3),
  extract_growth_results(growth_grim, "+ Height 9 + Mat Height", growth_grim$m4),
  extract_growth_results(growth_pheno, "Unadjusted", growth_pheno$m1),
  extract_growth_results(growth_pheno, "+ Demographics", growth_pheno$m2),
  extract_growth_results(growth_pheno, "+ Height at 9", growth_pheno$m3),
  extract_growth_results(growth_pheno, "+ Height 9 + Mat Height", growth_pheno$m4),
  extract_growth_results(growth_pace, "Unadjusted", growth_pace$m1),
  extract_growth_results(growth_pace, "+ Demographics", growth_pace$m2),
  extract_growth_results(growth_pace, "+ Height at 9", growth_pace$m3),
  extract_growth_results(growth_pace, "+ Height 9 + Mat Height", growth_pace$m4)
)

print(as.data.frame(growth_results))

cat("\n")
cat("=============================================================================\n")
cat("SUMMARY: ADULT HEIGHT (EAA Age 15 → Height at 22)\n")
cat("=============================================================================\n")

adult_results <- bind_rows(
  extract_growth_results(adult_grim, "Unadjusted", adult_grim$m1),
  extract_growth_results(adult_grim, "+ Demographics", adult_grim$m2),
  extract_growth_results(adult_grim, "+ Maternal Height", adult_grim$m3),
  extract_growth_results(adult_grim, "+ Height at 15", adult_grim$m4),
  extract_growth_results(adult_pheno, "Unadjusted", adult_pheno$m1),
  extract_growth_results(adult_pheno, "+ Demographics", adult_pheno$m2),
  extract_growth_results(adult_pheno, "+ Maternal Height", adult_pheno$m3),
  extract_growth_results(adult_pheno, "+ Height at 15", adult_pheno$m4),
  extract_growth_results(adult_pace, "Unadjusted", adult_pace$m1),
  extract_growth_results(adult_pace, "+ Demographics", adult_pace$m2),
  extract_growth_results(adult_pace, "+ Maternal Height", adult_pace$m3),
  extract_growth_results(adult_pace, "+ Height at 15", adult_pace$m4)
)

print(as.data.frame(adult_results))

# =============================================================================
# Create Output Directory and Save
# =============================================================================

if (!dir.exists("output")) dir.create("output")

write.csv(growth_results, "output/height_growth_results.csv", row.names = FALSE)
write.csv(adult_results, "output/height_adult_results.csv", row.names = FALSE)

# =============================================================================
# Forest Plot: Growth Model + Maternal Height Results
# =============================================================================

cat("\n=== CREATING FOREST PLOT ===\n")

# Get growth model with height at 9 + maternal height
growth_full <- growth_results %>% filter(Model == "+ Height 9 + Mat Height")

# Prepare data for plotting
growth_full_plot <- growth_full %>%
  mutate(
    Clock = factor(Clock, levels = c("DunedinPACE", "PhenoAge", "GrimAge")),
    Estimate = estimate,
    Lower = estimate - 1.96 * std.error,
    Upper = estimate + 1.96 * std.error,
    sig = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      p.value < 0.10 ~ "†",
      TRUE ~ ""
    ),
    label = sprintf("%.3f%s", Estimate, sig)
  )

# Create forest plot
p_forest_growth <- ggplot(growth_full_plot, aes(x = Estimate, y = Clock)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.2, linewidth = 0.8) +
  geom_point(size = 4, shape = 18, color = "#2E86AB") +
  geom_text(aes(label = label), vjust = -1, size = 3.5) +
  labs(
    x = "Coefficient (95% CI)",
    y = "",
    title = "",
    subtitle = "",
    caption = "Note: Controlling for height at age 9 and maternal height.\nCell composition and array type adjusted.\n† p<0.10, * p<0.05, ** p<0.01, *** p<0.001"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10, color = "gray40"),
    plot.caption = element_text(size = 9, hjust = 0),
    axis.text.y = element_text(size = 11, face = "bold"),
    panel.grid.major.y = element_blank()
  ) +
  scale_x_continuous(limits = c(-0.2, 0.1))

print(p_forest_growth)
ggsave("output/Figure_Growth_Forest.png", p_forest_growth, width = 8, height = 4, dpi = 300, bg = "white")
ggsave("output/Figure_Growth_Forest.pdf", p_forest_growth, width = 8, height = 4)

# =============================================================================
# Comparison Forest Plot: All Models Side by Side
# =============================================================================

# Prepare all growth models for comparison plot
growth_compare <- growth_results %>%
  mutate(
    Clock = factor(Clock, levels = c("DunedinPACE", "PhenoAge", "GrimAge")),
    Model = factor(Model, levels = c("Unadjusted", "+ Demographics", "+ Height at 9", "+ Height 9 + Mat Height")),
    Estimate = estimate,
    Lower = estimate - 1.96 * std.error,
    Upper = estimate + 1.96 * std.error,
    sig = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      p.value < 0.10 ~ "†",
      TRUE ~ ""
    )
  )

# Create comparison forest plot
p_forest_compare <- ggplot(growth_compare, aes(x = Estimate, y = Clock, color = Model, shape = Model)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), 
                 height = 0.2, linewidth = 0.6,
                 position = position_dodge(width = 0.6)) +
  geom_point(size = 3, position = position_dodge(width = 0.6)) +
  labs(
    x = "Coefficient (95% CI)",
    y = "",
    title = "EAA at Age 9 and Height Growth: Model Comparison",
    subtitle = "Effect of EAA on standardized self-reported height at 15 across model specifications",
    caption = "Cell composition and array type adjusted. † p<0.10, * p<0.05, ** p<0.01, *** p<0.001",
    color = "Model",
    shape = "Model"
  ) +
  scale_color_manual(values = c("Unadjusted" = "#E41A1C", 
                                "+ Demographics" = "#377EB8",
                                "+ Height at 9" = "#4DAF4A",
                                "+ Height 9 + Mat Height" = "#984EA3")) +
  scale_shape_manual(values = c(16, 17, 15, 18)) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10, color = "gray40"),
    plot.caption = element_text(size = 9, hjust = 0),
    axis.text.y = element_text(size = 11, face = "bold"),
    panel.grid.major.y = element_blank(),
    legend.position = "bottom"
  ) +
  scale_x_continuous(limits = c(-0.2, 0.15))

print(p_forest_compare)
ggsave("output/Figure_Growth_Forest_Compare.png", p_forest_compare, width = 10, height = 5, dpi = 300, bg = "white")
ggsave("output/Figure_Growth_Forest_Compare.pdf", p_forest_compare, width = 10, height = 5)

# =============================================================================
# Scatter Plot: Growth Model (Height residualized on baseline + maternal height + demographics)
# =============================================================================

cat("\n=== CREATING GROWTH SCATTER PLOTS ===\n")

# Function to create growth scatter plot (adjusted for all controls to match regression)
create_growth_scatter <- function(df, clock_label, color = "#4DAF4A") {
  
  # Filter to complete cases (including all controls)
  d <- df %>% 
    filter(!is.na(eaa_adjusted) & !is.na(z_height_15) & !is.na(child_height_9) & 
             !is.na(mat_height_clean) & !is.na(female) & !is.na(race))
  
  # Residualize height at 15 on ALL controls (to match regression model)
  height_model <- lm(z_height_15 ~ child_height_9 + mat_height_clean + female + factor(race), data = d)
  d$height_growth <- residuals(height_model)
  
  # Also residualize EAA on same controls (partial regression plot approach)
  eaa_model <- lm(eaa_adjusted ~ child_height_9 + mat_height_clean + female + factor(race), data = d)
  d$eaa_resid <- residuals(eaa_model)
  
  # Fit line for residualized EAA predicting residualized height
  fit <- lm(height_growth ~ eaa_resid, data = d)
  r_sq <- summary(fit)$r.squared
  p_val <- summary(fit)$coefficients[2, 4]
  coef_val <- summary(fit)$coefficients[2, 1]
  
  p <- ggplot(d, aes(x = eaa_resid, y = height_growth)) +
    geom_point(alpha = 0.2, size = 1.5, color = color) +
    geom_smooth(method = "lm", se = TRUE, color = color, fill = color, alpha = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    labs(
      title = clock_label,
      subtitle = sprintf("b = %.3f, p = %.3f", coef_val, p_val),
      x = "EAA at Age 9",
      y = "Height Growth\n(Residualized)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9, color = "gray40")
    )
  
  return(p)
}

# Create scatter plots for each clock
p_grim_growth <- create_growth_scatter(grim_adj_9, "GrimAge", "#E41A1C")
p_pheno_growth <- create_growth_scatter(pheno_adj_9, "PhenoAge", "#377EB8")
p_pace_growth <- create_growth_scatter(pace_adj_9, "DunedinPACE", "#4DAF4A")

# Combined figure
p_growth_scatter <- p_grim_growth + p_pheno_growth + p_pace_growth +
  plot_annotation(
    title = "",
    subtitle = "",
    caption = "Note: Epigenetic Age Acceleration at Age 9 and Height Growth to Age 15. \nPartial regression plot: Both axes residualized on Height at 9, Maternal Height, Sex, and Self-reported race. \nCell composition and array type adjusted.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "gray40"),
      plot.caption = element_text(size = 9)
    )
  )

print(p_growth_scatter)
ggsave("output/Figure_Growth_Scatter.png", p_growth_scatter, 
       width = 14, height = 5, dpi = 300, bg = "white")
ggsave("output/Figure_Growth_Scatter.pdf", p_growth_scatter, 
       width = 14, height = 5)

# =============================================================================
# Key Findings
# =============================================================================

cat("\n")
cat("=============================================================================\n")
cat("KEY FINDINGS\n")
cat("=============================================================================\n")

cat("\n1. GROWTH ANALYSIS (EAA Age 9 → Height at 15, + Height at 9):\n")
growth_h9 <- growth_results %>% filter(Model == "+ Height at 9")
for (i in 1:nrow(growth_h9)) {
  r <- growth_h9[i,]
  sig_str <- if(r$p.value < 0.001) " ***" else if(r$p.value < 0.01) " **" else if(r$p.value < 0.05) " *" else if(r$p.value < 0.10) " †" else ""
  cat(sprintf("   %s: b = %.4f, SE = %.4f, p = %.4f%s\n",
              r$Clock, r$estimate, r$std.error, r$p.value, sig_str))
}

cat("\n2. GROWTH ANALYSIS (EAA Age 9 → Height at 15, + Height at 9 + Maternal Height):\n")
growth_full_results <- growth_results %>% filter(Model == "+ Height 9 + Mat Height")
for (i in 1:nrow(growth_full_results)) {
  r <- growth_full_results[i,]
  sig_str <- if(r$p.value < 0.001) " ***" else if(r$p.value < 0.01) " **" else if(r$p.value < 0.05) " *" else if(r$p.value < 0.10) " †" else ""
  cat(sprintf("   %s: b = %.4f, SE = %.4f, p = %.4f%s (N=%d)\n",
              r$Clock, r$estimate, r$std.error, r$p.value, sig_str, r$N))
}

cat("\n3. ADULT HEIGHT (EAA Age 15 → Height at 22):\n")
if (nrow(adult_results) > 0) {
  adult_demo <- adult_results %>% filter(Model == "+ Demographics")
  for (i in 1:nrow(adult_demo)) {
    r <- adult_demo[i,]
    sig_str <- if(r$p.value < 0.001) " ***" else if(r$p.value < 0.01) " **" else if(r$p.value < 0.05) " *" else if(r$p.value < 0.10) " †" else ""
    cat(sprintf("   %s: b = %.4f, SE = %.4f, p = %.4f%s (N=%d)\n",
                r$Clock, r$estimate, r$std.error, r$p.value, sig_str, r$N))
  }
} else {
  cat("   No results available\n")
}

# =============================================================================
# Save Results
# =============================================================================

cat("\n=== FILES SAVED ===\n")
cat("- output/height_growth_results.csv\n")
cat("- output/height_adult_results.csv\n")
cat("- output/Figure_Growth_Forest.png/.pdf (main forest plot)\n")
cat("- output/Figure_Growth_Forest_Compare.png/.pdf (model comparison)\n")
cat("- output/Figure_Growth_Scatter.png/.pdf (scatter plot - growth model)\n")

cat("\n=== ANALYSIS COMPLETE ===\n")
