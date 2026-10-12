# Replication code for: "Maternal Income, Education, and Children's Epigenetic Age Acceleration"
# Manuscript ID: DEM-2024-0306 | Journal: Demography
# Input: data_2026_01_05.dta (FFCWS restricted-use data)
# =============================================================================
# Menarche & EAA: Pearson Correlation
# Individual-level correlation between age at menarche and EAA
# All three clocks, both waves, cell-composition adjusted
# =============================================================================

library(tidyverse)
library(haven)
library(ggplot2)
library(patchwork)

# Load data
data <- read_dta("data_2026_01_05.dta")

cat("=== DATA LOADED ===\n")

# =============================================================================
# Filter to females with menarche data
# =============================================================================

female_data <- data %>%
  filter(female == 1) %>%
  filter(!is.na(k6d10))

cat("Female sample with menarche:", nrow(female_data), "\n")

# =============================================================================
# Cell Composition Adjustment
# =============================================================================

get_cell_adjusted_data <- function(data, clock_var, cell_epi_var, cell_ic_var, array_var) {
  complete_data <- data %>%
    filter(!is.na(!!sym(clock_var)) & !is.na(!!sym(cell_epi_var)) &
           !is.na(!!sym(cell_ic_var)) & !is.na(!!sym(array_var)))
  formula <- paste(clock_var, "~", cell_epi_var, "+", cell_ic_var, "+ factor(", array_var, ")")
  model <- lm(as.formula(formula), data = complete_data)
  complete_data$adjusted_clock <- residuals(model)
  return(complete_data)
}

# Wave 4 (Age 9)
grim_w4  <- get_cell_adjusted_data(female_data, "z_grimAge_9", "cell_epi_9", "cell_ic_9", "array_type_9")
pheno_w4 <- get_cell_adjusted_data(female_data, "z_phenoAge_9", "cell_epi_9", "cell_ic_9", "array_type_9")
pace_w4  <- get_cell_adjusted_data(female_data, "z_poam45_9", "cell_epi_9", "cell_ic_9", "array_type_9")

# Wave 5 (Age 15)
grim_w5  <- get_cell_adjusted_data(female_data, "z_grimAge_15", "cell_epi_15", "cell_ic_15", "array_type_15")
pheno_w5 <- get_cell_adjusted_data(female_data, "z_phenoAge_15", "cell_epi_15", "cell_ic_15", "array_type_15")
pace_w5  <- get_cell_adjusted_data(female_data, "z_poam45_15", "cell_epi_15", "cell_ic_15", "array_type_15")

# =============================================================================
# Run Pearson correlations
# =============================================================================

run_cor <- function(df, clock_label, wave_label) {
  ct <- cor.test(df$k6d10, df$adjusted_clock)
  tibble(
    Clock = clock_label,
    Wave = wave_label,
    N = nrow(df),
    r = ct$estimate,
    p = ct$p.value
  )
}

cor_results <- bind_rows(
  run_cor(grim_w4,  "GrimAge",     "Age 9"),
  run_cor(pace_w4,  "DunedinPACE", "Age 9"),
  run_cor(pheno_w4, "PhenoAge",    "Age 9"),
  run_cor(grim_w5,  "GrimAge",     "Age 15"),
  run_cor(pace_w5,  "DunedinPACE", "Age 15"),
  run_cor(pheno_w5, "PhenoAge",    "Age 15")
)

cat("\n=== PEARSON CORRELATIONS ===\n")
print(as.data.frame(cor_results))

# =============================================================================
# Plotting function — jittered scatter with correlation annotation
# =============================================================================

plot_cor <- function(df, clock_label, color) {

  ct <- cor.test(df$k6d10, df$adjusted_clock)
  r_val <- ct$estimate
  p_val <- ct$p.value

  # Format p-value
  if (p_val < 0.001) {
    p_label <- "p < 0.001"
  } else {
    p_label <- sprintf("p = %.3f", p_val)
  }
  cor_label <- sprintf("r = %.3f, %s", r_val, p_label)

  p <- ggplot(df, aes(x = k6d10, y = adjusted_clock)) +
    geom_jitter(width = 0.15, height = 0, alpha = 0.2, size = 1.2, color = color) +
    geom_smooth(method = "lm", se = TRUE, color = color, fill = color,
                alpha = 0.15, linewidth = 1) +
    geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed", color = "gray50") +
    annotate("text", x = 10.2, y = max(df$adjusted_clock) * 0.85,
             label = cor_label, hjust = 0, size = 3.5, fontface = "bold") +
    scale_x_continuous(breaks = 10:14) +
    labs(
      x = "Age at First Menstruation",
      y = "EAA (Cell-Adjusted)",
      title = clock_label
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      plot.title = element_text(face = "bold", size = 13)
    )

  return(p)
}

if (!dir.exists("output")) dir.create("output")

# =============================================================================
# Age 15 (Wave 5) — main figure
# =============================================================================

cat("\n=== CREATING AGE 15 FIGURE ===\n")

p_grim_15  <- plot_cor(grim_w5,  "GrimAge",     "#5B8C5A")
p_pace_15  <- plot_cor(pace_w5,  "DunedinPACE", "#B5874A")
p_pheno_15 <- plot_cor(pheno_w5, "PhenoAge",    "#4A7FB5")

p_age15 <- p_grim_15 + p_pace_15 + p_pheno_15 +
  plot_layout(ncol = 3) +
  plot_annotation(
    title = "Epigenetic Age Acceleration at Age 15 and Age at Menarche",
    caption = "Note: EAA adjusted for cell composition and array type. Points jittered for visibility.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.caption = element_text(size = 9, hjust = 0)
    )
  )

ggsave("output/Figure_Menarche_Correlation_Age15.png", p_age15,
       width = 15, height = 5.5, dpi = 300, bg = "white")
ggsave("output/Figure_Menarche_Correlation_Age15.pdf", p_age15,
       width = 15, height = 5.5)

cat("Saved Age 15 figure.\n")

# =============================================================================
# Age 9 (Wave 4)
# =============================================================================

cat("\n=== CREATING AGE 9 FIGURE ===\n")

p_grim_9  <- plot_cor(grim_w4,  "GrimAge",     "#5B8C5A")
p_pace_9  <- plot_cor(pace_w4,  "DunedinPACE", "#B5874A")
p_pheno_9 <- plot_cor(pheno_w4, "PhenoAge",    "#4A7FB5")

p_age9 <- p_grim_9 + p_pace_9 + p_pheno_9 +
  plot_layout(ncol = 3) +
  plot_annotation(
    title = "Epigenetic Age Acceleration at Age 9 and Age at Menarche",
    caption = "Note: EAA adjusted for cell composition and array type. Points jittered for visibility.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.caption = element_text(size = 9, hjust = 0)
    )
  )

ggsave("output/Figure_Menarche_Correlation_Age9.png", p_age9,
       width = 15, height = 5.5, dpi = 300, bg = "white")
ggsave("output/Figure_Menarche_Correlation_Age9.pdf", p_age9,
       width = 15, height = 5.5)

cat("Saved Age 9 figure.\n")

# =============================================================================
# Save correlation results
# =============================================================================

write.csv(cor_results, "output/menarche_pearson_correlations.csv", row.names = FALSE)

cat("\n=== CORRELATION RESULTS ===\n")
cor_results %>%
  mutate(
    r = round(r, 3),
    p = round(p, 4),
    sig = case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", p < 0.10 ~ "\u2020", TRUE ~ "")
  ) %>%
  print()

cat("\n=== DONE ===\n")
