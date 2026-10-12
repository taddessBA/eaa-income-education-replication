# Replication code for: "Maternal Income, Education, and Children's Epigenetic Age Acceleration"
# Manuscript ID: DEM-2024-0306 | Journal: Demography
# Input: data_2026_01_05.dta (FFCWS restricted-use data)
# =============================================================================
# Menarche & EAA Figure — Refined Option 1, Wave 4 (Age 9)
# Mean EAA at Age 9 by age at menarche (10-14)
# Order: GrimAge, DunedinPACE, PhenoAge
# Error bars = ±1 SE; point size scaled by sample size
# =============================================================================

library(tidyverse)
library(haven)
library(estimatr)
library(ggplot2)
library(patchwork)

# Load data
data <- read_dta("data_2026_01_05.dta")

cat("=== DATA LOADED ===\n")

# =============================================================================
# Filter to females with menarche data + Wave 4 epigenetic data
# =============================================================================

female_data <- data %>%
  filter(female == 1) %>%
  filter(!is.na(k6d10)) %>%
  filter(!is.na(z_grimAge_9) | !is.na(z_phenoAge_9) | !is.na(z_poam45_9))

cat("Female sample:", nrow(female_data), "\n")
cat("\nAge at menarche distribution:\n")
print(table(female_data$k6d10))

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

grim_adj  <- get_cell_adjusted_data(female_data, "z_grimAge_9", "cell_epi_9", "cell_ic_9", "array_type_9")
pheno_adj <- get_cell_adjusted_data(female_data, "z_phenoAge_9", "cell_epi_9", "cell_ic_9", "array_type_9")
pace_adj  <- get_cell_adjusted_data(female_data, "z_poam45_9", "cell_epi_9", "cell_ic_9", "array_type_9")

# =============================================================================
# Calculate per-age means for all clocks
# =============================================================================

calc_per_age_means <- function(df, clock_label) {
  df %>%
    filter(k6d10 >= 10 & k6d10 <= 14) %>%
    group_by(k6d10) %>%
    summarise(
      mean_eaa = mean(adjusted_clock, na.rm = TRUE),
      se = sd(adjusted_clock, na.rm = TRUE) / sqrt(n()),
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(Clock = clock_label)
}

per_age <- bind_rows(
  calc_per_age_means(grim_adj,  "GrimAge"),
  calc_per_age_means(pheno_adj, "PhenoAge"),
  calc_per_age_means(pace_adj,  "DunedinPACE")
)

cat("\n=== PER-AGE MEANS ===\n")
print(as.data.frame(per_age))

if (!dir.exists("output")) dir.create("output")

# =============================================================================
# Plotting function
# =============================================================================

plot_menarche_eaa <- function(per_age_df, clock_name, color) {

  df <- per_age_df %>% filter(Clock == clock_name)

  p <- ggplot(df, aes(x = k6d10, y = mean_eaa)) +
    geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed", color = "gray50") +
    geom_errorbar(aes(ymin = mean_eaa - se, ymax = mean_eaa + se),
                  width = 0.12, color = color, linewidth = 0.6) +
    geom_point(aes(size = n), color = color) +
    scale_size_continuous(range = c(3, 7), guide = "none") +
    scale_x_continuous(breaks = 10:14) +
    labs(
      x = "Age at First Menstruation",
      y = "Mean EAA (Cell-Adjusted)",
      title = clock_name
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      plot.title = element_text(face = "bold", size = 13)
    )

  return(p)
}

# =============================================================================
# Create panels: GrimAge, DunedinPACE, PhenoAge
# =============================================================================

p_grim  <- plot_menarche_eaa(per_age, "GrimAge",     "#5B8C5A")
p_pace  <- plot_menarche_eaa(per_age, "DunedinPACE", "#B5874A")
p_pheno <- plot_menarche_eaa(per_age, "PhenoAge",    "#4A7FB5")

p_combined <- p_grim + p_pace + p_pheno +
  plot_layout(ncol = 3) +
  plot_annotation(
    title = "Mean Epigenetic Age Acceleration at Age 9 by Age at Menarche",
    subtitle = "Error bars = \u00B11 SE | Point size reflects sample size",
    caption = "Note: EAA adjusted for cell composition and array type.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "gray40"),
      plot.caption = element_text(size = 9, hjust = 0)
    )
  )

ggsave("output/Figure_Menarche_Option1_Refined_W4.png", p_combined,
       width = 15, height = 5.5, dpi = 300, bg = "white")
ggsave("output/Figure_Menarche_Option1_Refined_W4.pdf", p_combined,
       width = 15, height = 5.5)

cat("\nSaved: output/Figure_Menarche_Option1_Refined_W4.png/.pdf\n")

# =============================================================================
# Print summary
# =============================================================================

cat("\n=== SUMMARY TABLE ===\n")
summary_wide <- per_age %>%
  select(k6d10, Clock, mean_eaa, se, n) %>%
  mutate(mean_eaa = round(mean_eaa, 3), se = round(se, 3)) %>%
  pivot_wider(names_from = Clock, values_from = c(mean_eaa, se, n))
print(as.data.frame(summary_wide))

write.csv(per_age, "output/menarche_per_age_means_age9_refined.csv", row.names = FALSE)

cat("\n=== DONE ===\n")
