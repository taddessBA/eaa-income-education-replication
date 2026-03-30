# Replication code for: "Maternal Income, Education, and Children's Epigenetic Age Acceleration"
# Manuscript ID: DEM-2024-0306 | Journal: Demography
# Input: data_2026_01_05.dta (FFCWS restricted-use data)
# =============================================================================
# Menarche & EAA Trajectory: Within-Person Change from Age 9 to Age 15
# Girls with EAA at both waves, split by Early (<=12) vs On Time (>=13)
# Shows: high EAA at 9 → early menarche → EAA continues rising by 15
# =============================================================================

library(tidyverse)
library(haven)
library(ggplot2)
library(patchwork)

data <- read_dta("data_2026_01_05.dta")

cat("=== DATA LOADED ===\n")

# =============================================================================
# Get females with complete EAA at BOTH waves + menarche
# =============================================================================

females <- data %>%
  filter(female == 1, !is.na(k6d10)) %>%
  filter(!is.na(z_grimAge_9) & !is.na(z_grimAge_15) &
         !is.na(z_phenoAge_9) & !is.na(z_phenoAge_15) &
         !is.na(z_poam45_9) & !is.na(z_poam45_15) &
         !is.na(cell_epi_9) & !is.na(cell_ic_9) & !is.na(array_type_9) &
         !is.na(cell_epi_15) & !is.na(cell_ic_15) & !is.na(array_type_15))

cat("Girls with EAA at both waves + menarche:", nrow(females), "\n")

# =============================================================================
# Cell-adjust each clock at each wave
# =============================================================================

adjust <- function(df, clock_var, cell_epi, cell_ic, array_var) {
  formula <- paste(clock_var, "~", cell_epi, "+", cell_ic, "+ factor(", array_var, ")")
  residuals(lm(as.formula(formula), data = df))
}

females$eaa_grim_9   <- adjust(females, "z_grimAge_9",   "cell_epi_9",  "cell_ic_9",  "array_type_9")
females$eaa_grim_15  <- adjust(females, "z_grimAge_15",  "cell_epi_15", "cell_ic_15", "array_type_15")
females$eaa_pheno_9  <- adjust(females, "z_phenoAge_9",  "cell_epi_9",  "cell_ic_9",  "array_type_9")
females$eaa_pheno_15 <- adjust(females, "z_phenoAge_15", "cell_epi_15", "cell_ic_15", "array_type_15")
females$eaa_pace_9   <- adjust(females, "z_poam45_9",    "cell_epi_9",  "cell_ic_9",  "array_type_9")
females$eaa_pace_15  <- adjust(females, "z_poam45_15",   "cell_epi_15", "cell_ic_15", "array_type_15")

# =============================================================================
# Create long-format data for plotting
# =============================================================================

make_long <- function(df, var9, var15, clock_label) {
  bind_rows(
    df %>% transmute(idnum, k6d10, timing = ifelse(k6d10 <= 12, "Early (\u226412)", "On Time (\u226513)"),
                     wave = "Age 9", eaa = .data[[var9]], Clock = clock_label),
    df %>% transmute(idnum, k6d10, timing = ifelse(k6d10 <= 12, "Early (\u226412)", "On Time (\u226513)"),
                     wave = "Age 15", eaa = .data[[var15]], Clock = clock_label)
  ) %>%
    mutate(
      timing = factor(timing, levels = c("Early (\u226412)", "On Time (\u226513)")),
      wave = factor(wave, levels = c("Age 9", "Age 15"))
    )
}

long_all <- bind_rows(
  make_long(females, "eaa_grim_9",  "eaa_grim_15",  "GrimAge"),
  make_long(females, "eaa_pace_9",  "eaa_pace_15",  "DunedinPACE"),
  make_long(females, "eaa_pheno_9", "eaa_pheno_15", "PhenoAge")
)

# =============================================================================
# Calculate group means for the trajectory plot
# =============================================================================

group_summary <- long_all %>%
  group_by(Clock, timing, wave) %>%
  summarise(
    mean_eaa = mean(eaa, na.rm = TRUE),
    se = sd(eaa, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

cat("\n=== GROUP SUMMARY ===\n")
print(as.data.frame(group_summary))

# For the arrow segments, create wide format
arrow_data <- group_summary %>%
  pivot_wider(names_from = wave, values_from = c(mean_eaa, se, n),
              names_sep = "_") %>%
  rename(mean_9 = `mean_eaa_Age 9`, mean_15 = `mean_eaa_Age 15`,
         se_9 = `se_Age 9`, se_15 = `se_Age 15`)

# =============================================================================
# MAIN FIGURE: Trajectory plot — connected means with arrows
# =============================================================================

if (!dir.exists("output")) dir.create("output")

plot_trajectory <- function(group_summary_df, arrow_df, clock_name, color_early, color_ontime) {

  gs <- group_summary_df %>% filter(Clock == clock_name)
  ad <- arrow_df %>% filter(Clock == clock_name)

  p <- ggplot(gs, aes(x = wave, y = mean_eaa, group = timing, color = timing)) +
    geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed", color = "gray50") +
    # Error bars (±1 SE)
    geom_errorbar(aes(ymin = mean_eaa - se, ymax = mean_eaa + se),
                  width = 0.08, linewidth = 0.6) +
    # Connecting line
    geom_line(linewidth = 1.2) +
    # Points
    geom_point(size = 4) +
    scale_color_manual(values = c(setNames(color_early, "Early (\u226412)"),
                                  setNames(color_ontime, "On Time (\u226513)"))) +
    labs(
      x = NULL,
      y = "Mean EAA (Cell-Adjusted)",
      title = clock_name,
      color = "Menarche Timing"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 13),
      legend.position = "bottom"
    )

  return(p)
}

p_grim  <- plot_trajectory(group_summary, arrow_data, "GrimAge",     "#C0392B", "#2E86C1")
p_pace  <- plot_trajectory(group_summary, arrow_data, "DunedinPACE", "#C0392B", "#2E86C1")
p_pheno <- plot_trajectory(group_summary, arrow_data, "PhenoAge",    "#C0392B", "#2E86C1")

p_combined <- p_grim + p_pace + p_pheno +
  plot_layout(ncol = 3, guides = "collect") +
  plot_annotation(
    caption = "Within-person EAA trajectories from age 9 to 15 by menarche timing. Girls with early menarche start with higher EAA and continue to accelerate.\nSame individuals measured at both waves (N = 838). EAA adjusted for cell composition and array type. Error bars = \u00B11 SE.",
    theme = theme(
      plot.caption = element_text(size = 9, hjust = 0, lineheight = 1.3),
      legend.position = "bottom"
    )
  )

ggsave("output/Figure_Menarche_EAA_Trajectory.png", p_combined,
       width = 16, height = 7, dpi = 300, bg = "white")
ggsave("output/Figure_Menarche_EAA_Trajectory.pdf", p_combined,
       width = 16, height = 7)

cat("\nSaved: output/Figure_Menarche_EAA_Trajectory.png/.pdf\n")

# =============================================================================
# ALSO: Per menarche age (10-14) trajectory — shows the full gradient
# =============================================================================

group_by_age <- long_all %>%
  filter(k6d10 >= 10, k6d10 <= 14) %>%
  group_by(Clock, k6d10, wave) %>%
  summarise(
    mean_eaa = mean(eaa, na.rm = TRUE),
    se = sd(eaa, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(k6d10 = factor(k6d10))

plot_trajectory_by_age <- function(df, clock_name) {

  d <- df %>% filter(Clock == clock_name)

  p <- ggplot(d, aes(x = wave, y = mean_eaa, group = k6d10, color = k6d10)) +
    geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed", color = "gray50") +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = mean_eaa - se, ymax = mean_eaa + se),
                  width = 0.06, linewidth = 0.5) +
    scale_color_manual(
      values = c("10" = "#C0392B", "11" = "#E67E22", "12" = "#F4D03F",
                 "13" = "#27AE60", "14" = "#2E86C1"),
      labels = c("10" = "Age 10", "11" = "Age 11", "12" = "Age 12",
                 "13" = "Age 13", "14" = "Age 14")
    ) +
    labs(
      x = NULL,
      y = "Mean EAA (Cell-Adjusted)",
      title = clock_name,
      color = "Menarche Age"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 13),
      legend.position = "bottom"
    )

  return(p)
}

p2_grim  <- plot_trajectory_by_age(group_by_age, "GrimAge")
p2_pace  <- plot_trajectory_by_age(group_by_age, "DunedinPACE")
p2_pheno <- plot_trajectory_by_age(group_by_age, "PhenoAge")

p2_combined <- p2_grim + p2_pace + p2_pheno +
  plot_layout(ncol = 3, guides = "collect") +
  plot_annotation(
    title = "Within-Person EAA Trajectories from Age 9 to 15 by Age at Menarche",
    subtitle = "Each line tracks the same girls from age 9 to 15, grouped by when they started menstruation",
    caption = "Note: Same individuals at both waves (N = 838). EAA adjusted for cell composition and array type. Error bars = \u00B11 SE.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "gray40"),
      plot.caption = element_text(size = 9, hjust = 0),
      legend.position = "bottom"
    )
  )

ggsave("output/Figure_Menarche_EAA_Trajectory_ByAge.png", p2_combined,
       width = 15, height = 6, dpi = 300, bg = "white")
ggsave("output/Figure_Menarche_EAA_Trajectory_ByAge.pdf", p2_combined,
       width = 15, height = 6)

cat("Saved: output/Figure_Menarche_EAA_Trajectory_ByAge.png/.pdf\n")

# =============================================================================
# Save data
# =============================================================================

write.csv(group_summary, "output/menarche_eaa_trajectory_summary.csv", row.names = FALSE)
write.csv(as.data.frame(group_by_age), "output/menarche_eaa_trajectory_by_age.csv", row.names = FALSE)

cat("\n=== DONE ===\n")
