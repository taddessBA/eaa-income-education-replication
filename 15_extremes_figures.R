# Replication code for: "Maternal Income, Education, and Children's Epigenetic Age Acceleration"
# Manuscript ID: DEM-2024-0306 | Journal: Demography
# Input: data_2026_01_05.dta (FFCWS restricted-use data)
# =============================================================================
# Extremes vs On Time Bar Plots: All Three Clocks — Age 9 (Wave 4)
# Categories: Early (<=12), On Time (>=13)
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
# Filter to females with menarche data
# =============================================================================

female_data <- data %>%
  filter(female == 1) %>%
  filter(!is.na(k6d10)) %>%
  filter(!is.na(z_grimAge_9) | !is.na(z_phenoAge_9) | !is.na(z_poam45_9))

cat("Female sample with menarche and epigenetic data:", nrow(female_data), "\n")
cat("\nAge at menarche distribution:\n")
print(table(female_data$k6d10))

# =============================================================================
# Cell Composition Adjustment Function
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

# =============================================================================
# Get Cell + Array Adjusted Data — Wave 4 (Age 9) Only
# =============================================================================

grim_adj_w4 <- get_cell_adjusted_data(female_data, "z_grimAge_9", "cell_epi_9", "cell_ic_9", "array_type_9")
pheno_adj_w4 <- get_cell_adjusted_data(female_data, "z_phenoAge_9", "cell_epi_9", "cell_ic_9", "array_type_9")
pace_adj_w4 <- get_cell_adjusted_data(female_data, "z_poam45_9", "cell_epi_9", "cell_ic_9", "array_type_9")

cat("\nSample sizes (Wave 4 / Age 9):\n")
cat("GrimAge:", nrow(grim_adj_w4), "\n")
cat("PhenoAge:", nrow(pheno_adj_w4), "\n")
cat("DunedinPACE:", nrow(pace_adj_w4), "\n")

# =============================================================================
# Create menarche groups with updated labels
# =============================================================================

add_menarche_groups <- function(df) {
  df %>%
    filter(!is.na(k6d10)) %>%
    mutate(
      menarche_group = case_when(
        k6d10 <= 12 ~ "Early (\u226412)",
        k6d10 >= 13 ~ "On Time (\u226513)",
        TRUE ~ NA_character_
      ),
      menarche_group = factor(menarche_group,
                              levels = c("Early (\u226412)", "On Time (\u226513)"))
    ) %>%
    filter(!is.na(menarche_group))
}

# Apply to Wave 4 datasets
grim_adj_w4  <- add_menarche_groups(grim_adj_w4)
pheno_adj_w4 <- add_menarche_groups(pheno_adj_w4)
pace_adj_w4  <- add_menarche_groups(pace_adj_w4)

# =============================================================================
# Calculate group means and SEs
# =============================================================================

calc_group_means <- function(df, clock_label, wave_label) {
  df %>%
    group_by(menarche_group) %>%
    summarise(
      mean = mean(adjusted_clock, na.rm = TRUE),
      se = sd(adjusted_clock, na.rm = TRUE) / sqrt(n()),
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(
      Clock = clock_label,
      Wave = wave_label
    )
}

group_means <- bind_rows(
  calc_group_means(grim_adj_w4,  "GrimAge",     "Age 9"),
  calc_group_means(pheno_adj_w4, "PhenoAge",    "Age 9"),
  calc_group_means(pace_adj_w4,  "DunedinPACE", "Age 9")
)

cat("\n=== GROUP MEANS (ALL CLOCKS) ===\n")
print(as.data.frame(group_means))

# =============================================================================
# Extremes model: significance tests (Early as reference)
# =============================================================================

run_extremes_test <- function(df, clock_label, wave_label) {
  m <- lm_robust(adjusted_clock ~ menarche_group, data = df, se_type = "HC2")
  tidy_m <- broom::tidy(m)
  tidy_m$Clock <- clock_label
  tidy_m$Wave <- wave_label
  cat(sprintf("\n%s - %s:\n", clock_label, wave_label))
  print(tidy_m)
  return(tidy_m)
}

cat("\n=== EARLY VS ON TIME TESTS (AGE 9) ===\n")
ext_tests <- bind_rows(
  run_extremes_test(grim_adj_w4,  "GrimAge",     "Age 9"),
  run_extremes_test(pheno_adj_w4, "PhenoAge",    "Age 9"),
  run_extremes_test(pace_adj_w4,  "DunedinPACE", "Age 9")
)

# =============================================================================
# Plotting function for individual clock extremes figures
# =============================================================================

plot_extremes <- function(group_means_df, clock_name, bar_color) {

  df <- group_means_df %>% filter(Clock == clock_name)

  p <- ggplot(df, aes(x = menarche_group, y = mean)) +
    geom_bar(stat = "identity", fill = bar_color, width = 0.6) +
    geom_errorbar(aes(ymin = mean - 1.96 * se, ymax = mean + 1.96 * se),
                  width = 0.2) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    labs(
      x = "Menarche Timing Group",
      y = paste0("Mean Epigenetic Age Acceleration\n(", clock_name, ", Cell-Adjusted)"),
      title = clock_name
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.x = element_blank(),
      plot.title = element_text(face = "bold", size = 13)
    )

  return(p)
}

# Colors for each clock
color_grim  <- "#A8D5BA"
color_pheno <- "#A8C8E8"
color_pace  <- "#E8C8A8"

# =============================================================================
# Create individual figures for each clock
# =============================================================================

if (!dir.exists("output")) dir.create("output")

p_grim  <- plot_extremes(group_means, "GrimAge",     color_grim)
p_pheno <- plot_extremes(group_means, "PhenoAge",    color_pheno)
p_pace  <- plot_extremes(group_means, "DunedinPACE", color_pace)

ggsave("output/Figure_Extremes_GrimAge_W4.png",     p_grim,  width = 8, height = 6, dpi = 300, bg = "white")
ggsave("output/Figure_Extremes_GrimAge_W4.pdf",     p_grim,  width = 8, height = 6)
ggsave("output/Figure_Extremes_PhenoAge_W4.png",    p_pheno, width = 8, height = 6, dpi = 300, bg = "white")
ggsave("output/Figure_Extremes_PhenoAge_W4.pdf",    p_pheno, width = 8, height = 6)
ggsave("output/Figure_Extremes_DunedinPACE_W4.png", p_pace,  width = 8, height = 6, dpi = 300, bg = "white")
ggsave("output/Figure_Extremes_DunedinPACE_W4.pdf", p_pace,  width = 8, height = 6)

cat("\nSaved individual figures for all three clocks.\n")

# =============================================================================
# Combined 3-panel figure (all clocks together)
# =============================================================================

p_combined <- p_grim + p_pheno + p_pace +
  plot_layout(ncol = 3) +
  plot_annotation(
    title = "Epigenetic Age Acceleration at Age 9 by Menarche Timing",
    subtitle = "Early (\u226412) vs On Time (\u226513) menarche | Error bars = 95% CI",
    caption = "Note: Cell composition and array type adjusted. EAA measured at age 9.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "gray40"),
      plot.caption = element_text(size = 9, hjust = 0)
    )
  )

ggsave("output/Figure_Extremes_All_Clocks_W4.png", p_combined, width = 16, height = 6, dpi = 300, bg = "white")
ggsave("output/Figure_Extremes_All_Clocks_W4.pdf", p_combined, width = 16, height = 6)

cat("Saved combined 3-panel figure.\n")

# =============================================================================
# Also save vertical layout
# =============================================================================

p_combined_vert <- p_grim / p_pheno / p_pace +
  plot_annotation(
    title = "Epigenetic Age Acceleration at Age 9 by Menarche Timing",
    subtitle = "Early (\u226412) vs On Time (\u226513) menarche | Error bars = 95% CI",
    caption = "Note: Cell composition and array type adjusted. EAA measured at age 9.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "gray40"),
      plot.caption = element_text(size = 9, hjust = 0)
    )
  )

ggsave("output/Figure_Extremes_All_Clocks_W4_Vert.png", p_combined_vert, width = 8, height = 16, dpi = 300, bg = "white")
ggsave("output/Figure_Extremes_All_Clocks_W4_Vert.pdf", p_combined_vert, width = 8, height = 16)

cat("Saved vertical layout.\n")

# =============================================================================
# Save group means data
# =============================================================================

write.csv(group_means, "output/menarche_extremes_group_means_all_clocks_w4.csv", row.names = FALSE)

cat("\n=== ALL FILES SAVED ===\n")
cat("Individual figures:\n")
cat("  - output/Figure_Extremes_GrimAge_W4.png/.pdf\n")
cat("  - output/Figure_Extremes_PhenoAge_W4.png/.pdf\n")
cat("  - output/Figure_Extremes_DunedinPACE_W4.png/.pdf\n")
cat("Combined figures:\n")
cat("  - output/Figure_Extremes_All_Clocks_W4.png/.pdf (horizontal)\n")
cat("  - output/Figure_Extremes_All_Clocks_W4_Vert.png/.pdf (vertical)\n")
cat("Data:\n")
cat("  - output/menarche_extremes_group_means_all_clocks_w4.csv\n")
cat("\n=== DONE ===\n")
