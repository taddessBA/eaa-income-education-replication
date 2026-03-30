# Replication code for: "Maternal Income, Education, and Children's Epigenetic Age Acceleration"
# Manuscript ID: DEM-2024-0306 | Journal: Demography
# Input: data_2026_01_05.dta (FFCWS restricted-use data)
# =============================================================================
# Menarche and Epigenetic Age: Non-Linearity Tests
# Tests: Quadratic (U-shape), Extremes vs Middle, Spline
# =============================================================================

library(tidyverse)
library(haven)
library(estimatr)
library(broom)
library(ggplot2)
library(patchwork)
library(splines)

# Load data
data <- read_dta("data_2026_01_05.dta")

cat("=== DATA LOADED ===\n")

# =============================================================================
# Filter to females with menarche data
# =============================================================================

female_data <- data %>% 
  filter(female == 1) %>%
  filter(!is.na(k6d10)) %>%
  filter(!is.na(z_grimAge_9) | !is.na(z_phenoAge_9) | !is.na(z_poam45_9) |
           !is.na(z_grimAge_15) | !is.na(z_phenoAge_15) | !is.na(z_poam45_15))

cat("Female sample with menarche and epigenetic data:", nrow(female_data), "\n")
cat("\nAge at menarche distribution:\n")
print(table(female_data$k6d10))

# =============================================================================
# Cell Composition Adjustment Function
# =============================================================================

get_cell_adjusted_data <- function(data, clock_var, cell_epi_var, cell_ic_var, array_var = NULL) {

  complete_data <- data %>%
    filter(!is.na(!!sym(clock_var)) & !is.na(!!sym(cell_epi_var)) & !is.na(!!sym(cell_ic_var)))

  if (!is.null(array_var)) {
    complete_data <- complete_data %>% filter(!is.na(!!sym(array_var)))
    formula <- paste(clock_var, "~", cell_epi_var, "+", cell_ic_var, "+", array_var)
  } else {
    formula <- paste(clock_var, "~", cell_epi_var, "+", cell_ic_var)
  }
  model <- lm(as.formula(formula), data = complete_data)
  complete_data$adjusted_clock <- residuals(model)

  return(complete_data)
}

# =============================================================================
# Get Cell-Adjusted Data for Both Waves
# =============================================================================

# Wave 4 (Age 9) — including array type in residualization
grim_adj_w4 <- get_cell_adjusted_data(female_data, "z_grimAge_9", "cell_epi_9", "cell_ic_9", "array_type_9")
pheno_adj_w4 <- get_cell_adjusted_data(female_data, "z_phenoAge_9", "cell_epi_9", "cell_ic_9", "array_type_9")
pace_adj_w4 <- get_cell_adjusted_data(female_data, "z_poam45_9", "cell_epi_9", "cell_ic_9", "array_type_9")

# Wave 5 (Age 15) — including array type in residualization
grim_adj_w5 <- get_cell_adjusted_data(female_data, "z_grimAge_15", "cell_epi_15", "cell_ic_15", "array_type_15")
pheno_adj_w5 <- get_cell_adjusted_data(female_data, "z_phenoAge_15", "cell_epi_15", "cell_ic_15", "array_type_15")
pace_adj_w5 <- get_cell_adjusted_data(female_data, "z_poam45_15", "cell_epi_15", "cell_ic_15", "array_type_15")

cat("\nSample sizes (Wave 4):\n")
cat("GrimAge:", nrow(grim_adj_w4), "\n")
cat("PhenoAge:", nrow(pheno_adj_w4), "\n")
cat("DunedinPACE:", nrow(pace_adj_w4), "\n")

cat("\nSample sizes (Wave 5):\n")
cat("GrimAge:", nrow(grim_adj_w5), "\n")
cat("PhenoAge:", nrow(pheno_adj_w5), "\n")
cat("DunedinPACE:", nrow(pace_adj_w5), "\n")

# =============================================================================
# Non-Linearity Test Function
# =============================================================================

run_menarche_shape_tests <- function(df, y_var = "adjusted_clock", clock_label, wave_label) {
  
  d <- df %>%
    filter(!is.na(k6d10), !is.na(.data[[y_var]])) %>%
    mutate(
      # Center menarche age to reduce collinearity in quadratic model
      k6d10_c = k6d10 - mean(k6d10, na.rm = TRUE),
      k6d10_c2 = k6d10_c^2,
      
      # "Extremes matter" specification: tails vs middle (11-13 as reference)
      menarche_group = case_when(
        k6d10 <= 10 ~ "Early (<=10)",
        k6d10 %in% 11:13 ~ "Middle (11-13)",
        k6d10 >= 14 ~ "Late (>=14)",
        TRUE ~ NA_character_
      ),
      menarche_group = factor(menarche_group, levels = c("Middle (11-13)", "Early (<=10)", "Late (>=14)"))
    ) %>%
    filter(!is.na(menarche_group))
  
  cat("\n=============================\n")
  cat("Clock:", clock_label, "| Wave:", wave_label, "\n")
  cat("N:", nrow(d), "\n")
  cat("Menarche group counts:\n")
  print(table(d$menarche_group))
  cat("Mean menarche age:", round(mean(d$k6d10), 2), "\n")
  
  # ---- (A) Linear model (baseline)
  m_lin <- lm_robust(
    reformulate("k6d10", response = y_var),
    data = d,
    se_type = "HC2"
  )
  
  cat("\nLinear model:\n")
  print(tidy(m_lin) %>% filter(term == "k6d10"))
  
  # ---- (B) Quadratic model (U-shape test)
  m_quad <- lm_robust(
    reformulate(c("k6d10_c", "k6d10_c2"), response = y_var),
    data = d,
    se_type = "HC2"
  )
  
  # Turning point (in original k6d10 scale): -b1/(2*b2) + mean(k)
  b1 <- coef(m_quad)["k6d10_c"]
  b2 <- coef(m_quad)["k6d10_c2"]
  k_bar <- mean(d$k6d10, na.rm = TRUE)
  turning_point <- if (!is.na(b2) && b2 != 0) (-b1 / (2 * b2)) + k_bar else NA
  
  cat("\nQuadratic model:\n")
  print(tidy(m_quad))
  cat("Turning point:", round(turning_point, 2), "years\n")
  
  # ---- (C) Extremes matter model: tails vs middle
  m_ext <- lm_robust(
    reformulate("menarche_group", response = y_var),
    data = d,
    se_type = "HC2"
  )
  
  cat("\nExtremes vs Middle model:\n")
  print(tidy(m_ext))
  
  # ---- (D) Spline nonlinearity (robust check)
  m_spline <- lm_robust(
    as.formula(paste0(y_var, " ~ splines::ns(k6d10, df = 2)")),
    data = d,
    se_type = "HC2"
  )
  
  # ---- (E) Wald test for quadratic term
  wald_quad <- tidy(m_quad) %>% 
    filter(term == "k6d10_c2") %>% 
    transmute(
      test = "Quadratic term (U-shape): H0 b2=0",
      estimate, std.error, statistic, p.value
    )
  
  cat("\nWald test for quadratic term (U-shape):\n")
  print(wald_quad)
  
  # ---- (F) Spline vs Linear comparison (F-test)
  # Compare spline model to linear model using AIC/BIC
  m_lin_basic <- lm(reformulate("k6d10", response = y_var), data = d)
  m_spline_basic <- lm(as.formula(paste0(y_var, " ~ splines::ns(k6d10, df = 2)")), data = d)
  
  # Likelihood ratio / F-test
  anova_result <- tryCatch(
    anova(m_lin_basic, m_spline_basic),
    error = function(e) NULL
  )
  
  wald_spline <- list(
    anova = anova_result,
    AIC_linear = AIC(m_lin_basic),
    AIC_spline = AIC(m_spline_basic),
    BIC_linear = BIC(m_lin_basic),
    BIC_spline = BIC(m_spline_basic)
  )
  
  cat("\nSpline vs Linear comparison:\n")
  cat(sprintf("  AIC: Linear = %.2f, Spline = %.2f (diff = %.2f)\n", 
              wald_spline$AIC_linear, wald_spline$AIC_spline, 
              wald_spline$AIC_linear - wald_spline$AIC_spline))
  if (!is.null(anova_result)) {
    cat(sprintf("  F-test p-value: %.4f\n", anova_result$`Pr(>F)`[2]))
  }
  
  list(
    data = d,
    clock = clock_label,
    wave = wave_label,
    m_lin = m_lin,
    m_quad = m_quad,
    turning_point = turning_point,
    m_ext = m_ext,
    m_spline = m_spline,
    wald_quad = wald_quad,
    wald_spline = wald_spline
  )
}

# =============================================================================
# Run Tests for All Clocks, Both Waves
# =============================================================================

cat("\n")
cat("=============================================================================\n")
cat("WAVE 4 (AGE 9) NON-LINEARITY TESTS\n")
cat("=============================================================================\n")

tests_grim_w4  <- run_menarche_shape_tests(grim_adj_w4,  "adjusted_clock", "GrimAge", "Wave 4")
tests_pheno_w4 <- run_menarche_shape_tests(pheno_adj_w4, "adjusted_clock", "PhenoAge", "Wave 4")
tests_pace_w4  <- run_menarche_shape_tests(pace_adj_w4,  "adjusted_clock", "DunedinPACE", "Wave 4")

cat("\n")
cat("=============================================================================\n")
cat("WAVE 5 (AGE 15) NON-LINEARITY TESTS\n")
cat("=============================================================================\n")

tests_grim_w5  <- run_menarche_shape_tests(grim_adj_w5,  "adjusted_clock", "GrimAge", "Wave 5")
tests_pheno_w5 <- run_menarche_shape_tests(pheno_adj_w5, "adjusted_clock", "PhenoAge", "Wave 5")
tests_pace_w5  <- run_menarche_shape_tests(pace_adj_w5,  "adjusted_clock", "DunedinPACE", "Wave 5")

# =============================================================================
# Extract Key Results into Summary Table
# =============================================================================

extract_key_results <- function(tests_obj) {
  
  clock_label <- tests_obj$clock
  wave_label <- tests_obj$wave
  
  quad_tidy <- tidy(tests_obj$m_quad)
  ext_tidy  <- tidy(tests_obj$m_ext)
  lin_tidy  <- tidy(tests_obj$m_lin)
  
  # Linear term
  lin_row <- lin_tidy %>% filter(term == "k6d10") %>%
    transmute(
      Clock = clock_label,
      Wave = wave_label,
      Test = "Linear: menarche age",
      Estimate = estimate,
      SE = std.error,
      t = statistic,
      p = p.value,
      TurningPoint = as.numeric(NA)
    )
  
  # Quadratic term
  quad_b2 <- quad_tidy %>% filter(term == "k6d10_c2") %>%
    transmute(
      Clock = clock_label,
      Wave = wave_label,
      Test = "Quadratic term (U-shape): b2",
      Estimate = estimate,
      SE = std.error,
      t = statistic,
      p = p.value,
      TurningPoint = tests_obj$turning_point
    )
  
  # Extremes: Early vs Middle
  ext_early <- ext_tidy %>% filter(term == "menarche_groupEarly (<=10)") %>%
    transmute(
      Clock = clock_label,
      Wave = wave_label,
      Test = "Early (<=10) vs Middle (11-13)",
      Estimate = estimate,
      SE = std.error,
      t = statistic,
      p = p.value,
      TurningPoint = as.numeric(NA)
    )
  
  # Extremes: Late vs Middle
  ext_late <- ext_tidy %>% filter(term == "menarche_groupLate (>=14)") %>%
    transmute(
      Clock = clock_label,
      Wave = wave_label,
      Test = "Late (>=14) vs Middle (11-13)",
      Estimate = estimate,
      SE = std.error,
      t = statistic,
      p = p.value,
      TurningPoint = as.numeric(NA)
    )
  
  bind_rows(lin_row, quad_b2, ext_early, ext_late)
}

# Combine all results
shape_results <- bind_rows(
  extract_key_results(tests_grim_w4),
  extract_key_results(tests_pheno_w4),
  extract_key_results(tests_pace_w4),
  extract_key_results(tests_grim_w5),
  extract_key_results(tests_pheno_w5),
  extract_key_results(tests_pace_w5)
)

cat("\n")
cat("=============================================================================\n")
cat("SUMMARY: SHAPE TEST RESULTS\n")
cat("=============================================================================\n")

# Convert to data.frame for cleaner printing
shape_results_df <- as.data.frame(shape_results)
print(shape_results_df)

# =============================================================================
# Create Output Directory and Save
# =============================================================================

if (!dir.exists("output")) dir.create("output")

write.csv(shape_results, "output/menarche_shape_tests.csv", row.names = FALSE)

# =============================================================================
# Create Summary Table for Manuscript
# =============================================================================

# Focus on key tests: Linear, Quadratic, Early vs Middle
summary_table <- shape_results %>%
  mutate(
    sig = case_when(
      p < 0.001 ~ "***",
      p < 0.01 ~ "**",
      p < 0.05 ~ "*",
      p < 0.10 ~ "†",
      TRUE ~ ""
    ),
    result = sprintf("%.3f%s (%.3f)", Estimate, sig, SE)
  ) %>%
  select(Clock, Wave, Test, result, p, TurningPoint)

cat("\n=== FORMATTED SUMMARY TABLE ===\n")
summary_table_df <- as.data.frame(summary_table)
print(summary_table_df)

# =============================================================================
# Quadratic Fit Plots
# =============================================================================

cat("\n=== CREATING QUADRATIC FIT PLOTS ===\n")

plot_quad_fit <- function(tests_obj) {
  
  d <- tests_obj$data
  clock_label <- tests_obj$clock
  wave_label <- tests_obj$wave
  
  # Create prediction grid
  grid <- tibble(k6d10 = seq(min(d$k6d10), max(d$k6d10), by = 0.1)) %>%
    mutate(
      k6d10_c  = k6d10 - mean(d$k6d10, na.rm = TRUE),
      k6d10_c2 = k6d10_c^2
    )
  
  # Get predictions
  grid$yhat <- predict(tests_obj$m_quad, newdata = grid)
  
  # Get turning point
  tp <- tests_obj$turning_point
  
  p <- ggplot(d, aes(x = k6d10, y = adjusted_clock)) +
    geom_jitter(alpha = 0.15, width = 0.1, height = 0, size = 1) +
    geom_line(data = grid, aes(x = k6d10, y = yhat), 
              color = "#E41A1C", linewidth = 1.2) +
    geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed") +
    {if (!is.na(tp) && tp >= min(d$k6d10) && tp <= max(d$k6d10)) 
      geom_vline(xintercept = tp, linewidth = 0.5, linetype = "dotted", color = "blue")} +
    scale_x_continuous(breaks = 10:15) +
    labs(
      x = "Age at First Menstruation",
      y = paste0("Cell-Adjusted EAA (", clock_label, ")"),
      title = paste0(clock_label, " - ", wave_label),
      subtitle = if (!is.na(tp)) paste0("Turning point: ", round(tp, 1), " years") else NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9, color = "gray40")
    )
  
  return(p)
}

# Create plots for all clocks and waves
p_grim_w4  <- plot_quad_fit(tests_grim_w4)
p_pheno_w4 <- plot_quad_fit(tests_pheno_w4)
p_pace_w4  <- plot_quad_fit(tests_pace_w4)

p_grim_w5  <- plot_quad_fit(tests_grim_w5)
p_pheno_w5 <- plot_quad_fit(tests_pheno_w5)
p_pace_w5  <- plot_quad_fit(tests_pace_w5)

# Save individual plots
ggsave("output/Figure_QuadFit_GrimAge_W4.png", p_grim_w4, width = 6, height = 5, dpi = 300, bg = "white")
ggsave("output/Figure_QuadFit_PhenoAge_W4.png", p_pheno_w4, width = 6, height = 5, dpi = 300, bg = "white")
ggsave("output/Figure_QuadFit_DunedinPACE_W4.png", p_pace_w4, width = 6, height = 5, dpi = 300, bg = "white")

ggsave("output/Figure_QuadFit_GrimAge_W5.png", p_grim_w5, width = 6, height = 5, dpi = 300, bg = "white")
ggsave("output/Figure_QuadFit_PhenoAge_W5.png", p_pheno_w5, width = 6, height = 5, dpi = 300, bg = "white")
ggsave("output/Figure_QuadFit_DunedinPACE_W5.png", p_pace_w5, width = 6, height = 5, dpi = 300, bg = "white")

# Combined panel figure - ALL (both waves)
p_all_quad <- (p_grim_w4 | p_grim_w5) /
  (p_pheno_w4 | p_pheno_w5) /
  (p_pace_w4 | p_pace_w5) +
  plot_annotation(
    title = "Quadratic Fit: Epigenetic Age Acceleration by Age at Menarche",
    subtitle = "Red line shows quadratic fit; blue dotted line shows turning point (if within range)",
    caption = "Cell composition adjusted. Points jittered for visibility.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "gray40"),
      plot.caption = element_text(size = 9)
    )
  )

ggsave("output/Figure_QuadFit_All.png", p_all_quad, width = 12, height = 14, dpi = 300, bg = "white")
ggsave("output/Figure_QuadFit_All.pdf", p_all_quad, width = 12, height = 14)

# =============================================================================
# MANUSCRIPT FIGURE: Wave 4 Only (3 panels)
# Most informative since EAA measured BEFORE menarche
# =============================================================================

cat("\n=== CREATING WAVE 4 MANUSCRIPT FIGURE ===\n")

# Wave 4 combined panel (horizontal layout)
p_wave4_quad <- p_grim_w4 + p_pace_w4 + p_pheno_w4 +
  plot_annotation(
    title = "",
    subtitle = "",
    caption = "Note: Quadratic fit shown in red; vertical blue line indicates turning point. Cell composition adjusted. Points jittered for visibility. N = 865 females.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "gray40"),
      plot.caption = element_text(size = 9, hjust = 0)
    )
  )

print(p_wave4_quad)

ggsave("output/Figure_Manuscript_QuadFit_Wave4.png", p_wave4_quad, 
       width = 14, height = 5, dpi = 300, bg = "white")
ggsave("output/Figure_Manuscript_QuadFit_Wave4.pdf", p_wave4_quad, 
       width = 14, height = 5)

# Alternative: Vertical layout
p_wave4_quad_vert <- p_grim_w4 / p_pace_w4 / p_pheno_w4 +
  plot_annotation(
    title = "Epigenetic Age Acceleration at Age 9 Predicts Age at Menarche",
    subtitle = "Quadratic fit shown in red; vertical blue line indicates turning point",
    caption = "Note: Cell composition adjusted. Points jittered for visibility. N = 865 females.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "gray40"),
      plot.caption = element_text(size = 9, hjust = 0)
    )
  )

ggsave("output/Figure_Manuscript_QuadFit_Wave4_Vert.png", p_wave4_quad_vert, 
       width = 7, height = 12, dpi = 300, bg = "white")
ggsave("output/Figure_Manuscript_QuadFit_Wave4_Vert.pdf", p_wave4_quad_vert, 
       width = 7, height = 12)

cat("Saved: Figure_Manuscript_QuadFit_Wave4.png/.pdf (horizontal)\n")
cat("Saved: Figure_Manuscript_QuadFit_Wave4_Vert.png/.pdf (vertical)\n")

# =============================================================================
# Extremes Bar Plot
# =============================================================================

cat("\n=== CREATING EXTREMES VS MIDDLE PLOTS ===\n")

# Calculate means by group
calc_group_means <- function(tests_obj) {
  tests_obj$data %>%
    group_by(menarche_group) %>%
    summarise(
      mean = mean(adjusted_clock, na.rm = TRUE),
      se = sd(adjusted_clock, na.rm = TRUE) / sqrt(n()),
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(
      Clock = tests_obj$clock,
      Wave = tests_obj$wave
    )
}

# Combine all group means
group_means <- bind_rows(
  calc_group_means(tests_grim_w4),
  calc_group_means(tests_pheno_w4),
  calc_group_means(tests_pace_w4),
  calc_group_means(tests_grim_w5),
  calc_group_means(tests_pheno_w5),
  calc_group_means(tests_pace_w5)
)

cat("\n=== GROUP MEANS ===\n")
print(as.data.frame(group_means))

write.csv(group_means, "output/menarche_group_means.csv", row.names = FALSE)

# Create extremes plot for GrimAge (significant result)
grim_groups <- group_means %>% filter(Clock == "GrimAge")

p_extremes_grim <- ggplot(grim_groups, aes(x = menarche_group, y = mean, fill = Wave)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(width = 0.8), width = 0.25) +
  scale_fill_manual(values = c("Wave 4" = "#A8D5BA", "Wave 5" = "#5B8C5A")) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  labs(
    x = "Menarche Timing Group",
    y = "Mean Epigenetic Age Acceleration\n(GrimAge)",
    fill = NULL,
    title = "",
    caption = "Note: GrimAge by Menarche Timing: Early vs Middle vs Late."
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

print(p_extremes_grim)

ggsave("output/Figure_Extremes_GrimAge.png", p_extremes_grim, width = 8, height = 6, dpi = 300, bg = "white")

# =============================================================================
# Spline vs Linear Comparison
# =============================================================================

cat("\n=== SPLINE VS LINEAR COMPARISON (AIC) ===\n")

cat("\nWave 4:\n")
cat(sprintf("GrimAge:     AIC diff = %.2f (positive = spline better)\n", 
            tests_grim_w4$wald_spline$AIC_linear - tests_grim_w4$wald_spline$AIC_spline))
cat(sprintf("PhenoAge:    AIC diff = %.2f\n", 
            tests_pheno_w4$wald_spline$AIC_linear - tests_pheno_w4$wald_spline$AIC_spline))
cat(sprintf("DunedinPACE: AIC diff = %.2f\n", 
            tests_pace_w4$wald_spline$AIC_linear - tests_pace_w4$wald_spline$AIC_spline))

cat("\nWave 5:\n")
cat(sprintf("GrimAge:     AIC diff = %.2f\n", 
            tests_grim_w5$wald_spline$AIC_linear - tests_grim_w5$wald_spline$AIC_spline))
cat(sprintf("PhenoAge:    AIC diff = %.2f\n", 
            tests_pheno_w5$wald_spline$AIC_linear - tests_pheno_w5$wald_spline$AIC_spline))
cat(sprintf("DunedinPACE: AIC diff = %.2f\n", 
            tests_pace_w5$wald_spline$AIC_linear - tests_pace_w5$wald_spline$AIC_spline))

# =============================================================================
# Key Findings Summary
# =============================================================================

cat("\n")
cat("=============================================================================\n")
cat("KEY FINDINGS SUMMARY\n")
cat("=============================================================================\n")

# Filter to significant results
sig_results <- shape_results %>% filter(p < 0.10)

cat("\nSignificant results (p < 0.10):\n")
print(as.data.frame(sig_results))

cat("\n=== INTERPRETATION ===\n")
cat("\n1. LINEAR EFFECTS:\n")
lin_results <- shape_results %>% filter(Test == "Linear: menarche age")
for (i in 1:nrow(lin_results)) {
  r <- lin_results[i,]
  sig_str <- if(r$p < 0.05) " *SIGNIFICANT*" else ""
  cat(sprintf("   %s %s: b = %.3f, p = %.4f%s\n", 
              r$Clock, r$Wave, r$Estimate, r$p, sig_str))
}

cat("\n2. QUADRATIC (U-SHAPE) EFFECTS:\n")
quad_results <- shape_results %>% filter(Test == "Quadratic term (U-shape): b2")
for (i in 1:nrow(quad_results)) {
  r <- quad_results[i,]
  sig_str <- if(r$p < 0.05) " *SIGNIFICANT*" else ""
  tp_str <- if(!is.na(r$TurningPoint)) sprintf(", turning point = %.1f", r$TurningPoint) else ""
  cat(sprintf("   %s %s: b2 = %.3f, p = %.4f%s%s\n", 
              r$Clock, r$Wave, r$Estimate, r$p, sig_str, tp_str))
}

cat("\n3. EARLY (<=10) VS MIDDLE (11-13):\n")
early_results <- shape_results %>% filter(Test == "Early (<=10) vs Middle (11-13)")
for (i in 1:nrow(early_results)) {
  r <- early_results[i,]
  sig_str <- if(r$p < 0.05) " *SIGNIFICANT*" else ""
  cat(sprintf("   %s %s: b = %.3f, p = %.4f%s\n", 
              r$Clock, r$Wave, r$Estimate, r$p, sig_str))
}

cat("\n4. LATE (>=14) VS MIDDLE (11-13):\n")
late_results <- shape_results %>% filter(Test == "Late (>=14) vs Middle (11-13)")
for (i in 1:nrow(late_results)) {
  r <- late_results[i,]
  sig_str <- if(r$p < 0.05) " *SIGNIFICANT*" else ""
  cat(sprintf("   %s %s: b = %.3f, p = %.4f%s\n", 
              r$Clock, r$Wave, r$Estimate, r$p, sig_str))
}

# =============================================================================
# Output Files
# =============================================================================

cat("\n=== FILES SAVED ===\n")
cat("\nData:\n")
cat("- output/menarche_shape_tests.csv\n")
cat("- output/menarche_group_means.csv\n")

cat("\nManuscript Figures:\n")
cat("- output/Figure_Manuscript_QuadFit_Wave4.png/.pdf (3 panels, horizontal)\n")
cat("- output/Figure_Manuscript_QuadFit_Wave4_Vert.png/.pdf (3 panels, vertical)\n")

cat("\nSupplementary Figures:\n")
cat("- output/Figure_QuadFit_All.png/.pdf (6-panel, both waves)\n")
cat("- output/Figure_QuadFit_[Clock]_[Wave].png (individual)\n")
cat("- output/Figure_Extremes_GrimAge.png\n")

cat("\n=== ANALYSIS COMPLETE ===\n")

