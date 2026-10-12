# Replication code for: "Maternal Income, Education, and Children's Epigenetic Age Acceleration"
# Manuscript ID: DEM-2024-0306 | Journal: Demography
# Input: data_2026_01_05.dta (FFCWS restricted-use data)
# =============================================================================
# Menarche Figures: Manuscript and Supplementary Versions
# Excludes menarche age 15 (n=2)
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
# Filter to females with menarche data (EXCLUDING age 15 menarche)
# =============================================================================

female_data <- data %>% 
  filter(female == 1) %>%
  filter(!is.na(k6d10)) %>%
  filter(k6d10 >= 10 & k6d10 <= 14)  # Exclude age 15 (n=2)

cat("Female sample (menarche ages 10-14):", nrow(female_data), "\n")
cat("\nAge at menarche distribution:\n")
print(table(female_data$k6d10))

# =============================================================================
# Cell Composition Adjustment Functions
# =============================================================================

get_cell_adjusted_residuals <- function(data, clock_var, cell_epi_var, cell_ic_var) {
  complete_data <- data %>%
    filter(!is.na(!!sym(clock_var)) & !is.na(!!sym(cell_epi_var)) & !is.na(!!sym(cell_ic_var)))
  
  formula <- paste(clock_var, "~", cell_epi_var, "+", cell_ic_var)
  model <- lm(as.formula(formula), data = complete_data)
  complete_data$adjusted_clock <- residuals(model)
  
  return(complete_data %>% select(idnum, k6d10, adjusted_clock, !!sym(clock_var)))
}

calc_mean_se <- function(adj_data, adj_col_name) {
  adj_data %>%
    group_by(k6d10) %>%
    summarise(
      mean_adj = mean(adjusted_clock, na.rm = TRUE),
      se_adj = sd(adjusted_clock, na.rm = TRUE) / sqrt(n()),
      n = n(),
      .groups = "drop"
    )
}

# =============================================================================
# Get Adjusted Data for Both Waves
# =============================================================================

# Wave 4 (Age 9)
grim_adj_w4 <- get_cell_adjusted_residuals(female_data, "z_grimAge_9", "cell_epi_9", "cell_ic_9")
pheno_adj_w4 <- get_cell_adjusted_residuals(female_data, "z_phenoAge_9", "cell_epi_9", "cell_ic_9")
pace_adj_w4 <- get_cell_adjusted_residuals(female_data, "z_poam45_9", "cell_epi_9", "cell_ic_9")

# Wave 5 (Age 15)
grim_adj_w5 <- get_cell_adjusted_residuals(female_data, "z_grimAge_15", "cell_epi_15", "cell_ic_15")
pheno_adj_w5 <- get_cell_adjusted_residuals(female_data, "z_phenoAge_15", "cell_epi_15", "cell_ic_15")
pace_adj_w5 <- get_cell_adjusted_residuals(female_data, "z_poam45_15", "cell_epi_15", "cell_ic_15")

# Calculate summaries
grim_sum_w4 <- calc_mean_se(grim_adj_w4) %>% mutate(Clock = "GrimAge", Wave = "Age 9")
pheno_sum_w4 <- calc_mean_se(pheno_adj_w4) %>% mutate(Clock = "PhenoAge", Wave = "Age 9")
pace_sum_w4 <- calc_mean_se(pace_adj_w4) %>% mutate(Clock = "DunedinPACE", Wave = "Age 9")

grim_sum_w5 <- calc_mean_se(grim_adj_w5) %>% mutate(Clock = "GrimAge", Wave = "Age 15")
pheno_sum_w5 <- calc_mean_se(pheno_adj_w5) %>% mutate(Clock = "PhenoAge", Wave = "Age 15")
pace_sum_w5 <- calc_mean_se(pace_adj_w5) %>% mutate(Clock = "DunedinPACE", Wave = "Age 15")

# Combine all
all_summary <- bind_rows(
  grim_sum_w4, pheno_sum_w4, pace_sum_w4,
  grim_sum_w5, pheno_sum_w5, pace_sum_w5
)

cat("\n=== SUMMARY DATA ===\n")
print(all_summary)

# =============================================================================
# Correlations (for reporting)
# =============================================================================

cat("\n=== CORRELATIONS ===\n")

# Wave 4
cor_grim_w4 <- cor.test(grim_adj_w4$k6d10, grim_adj_w4$adjusted_clock)
cor_pheno_w4 <- cor.test(pheno_adj_w4$k6d10, pheno_adj_w4$adjusted_clock)
cor_pace_w4 <- cor.test(pace_adj_w4$k6d10, pace_adj_w4$adjusted_clock)

# Wave 5
cor_grim_w5 <- cor.test(grim_adj_w5$k6d10, grim_adj_w5$adjusted_clock)
cor_pheno_w5 <- cor.test(pheno_adj_w5$k6d10, pheno_adj_w5$adjusted_clock)
cor_pace_w5 <- cor.test(pace_adj_w5$k6d10, pace_adj_w5$adjusted_clock)

cat("\n                  Wave 4 (Age 9)          Wave 5 (Age 15)\n")
cat("Clock           r        p              r        p\n")
cat(paste(rep("-", 60), collapse = ""), "\n")
cat(sprintf("GrimAge      %6.3f   %6.4f         %6.3f   %6.4f\n", 
            cor_grim_w4$estimate, cor_grim_w4$p.value, cor_grim_w5$estimate, cor_grim_w5$p.value))
cat(sprintf("PhenoAge     %6.3f   %6.4f         %6.3f   %6.4f\n", 
            cor_pheno_w4$estimate, cor_pheno_w4$p.value, cor_pheno_w5$estimate, cor_pheno_w5$p.value))
cat(sprintf("DunedinPACE  %6.3f   %6.4f         %6.3f   %6.4f\n", 
            cor_pace_w4$estimate, cor_pace_w4$p.value, cor_pace_w5$estimate, cor_pace_w5$p.value))

# =============================================================================
# Create output directory
# =============================================================================

if (!dir.exists("output")) dir.create("output")

# =============================================================================
# MANUSCRIPT FIGURE: GrimAge at Both Waves (the significant result)
# =============================================================================

cat("\n=== CREATING MANUSCRIPT FIGURE ===\n")

# Filter to GrimAge only
grim_data <- all_summary %>% filter(Clock == "GrimAge")

# Create a clean, publication-ready figure
p_manuscript <- ggplot(grim_data, aes(x = factor(k6d10), y = mean_adj, fill = Wave)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = mean_adj - se_adj, ymax = mean_adj + se_adj),
                position = position_dodge(width = 0.8), width = 0.25, linewidth = 0.5) +
  scale_fill_manual(values = c("Age 9" = "#A8D5BA", "Age 15" = "#5B8C5A"),
                    labels = c("Age 9" = "Wave 4 (Age 9)", "Age 15" = "Wave 5 (Age 15)")) +
  labs(
    x = "Age at First Menstruation (years)",
    y = "Mean Epigenetic Age Acceleration\n(GrimAge, Cell-Adjusted)",
    fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  geom_hline(yintercept = 0, linetype = "solid", color = "gray40", linewidth = 0.3) 
  # Add sample sizes as text below bars
 # geom_text(aes(y = -0.25, label = paste0("n=", n)), 
            #position = position_dodge(width = 0.8), size = 2.5, color = "gray50")

print(p_manuscript)

ggsave("output/Figure_Manuscript_GrimAge_Menarche.png", p_manuscript, 
       width = 8, height = 6, dpi = 300, bg = "white")
ggsave("output/Figure_Manuscript_GrimAge_Menarche.pdf", p_manuscript, 
       width = 8, height = 6)

# Alternative: Faceted version
p_manuscript_facet <- ggplot(grim_data, aes(x = factor(k6d10), y = mean_adj)) +
  geom_bar(stat = "identity", fill = "#6B8E6B", width = 0.7) +
  geom_errorbar(aes(ymin = mean_adj - se_adj, ymax = mean_adj + se_adj),
                width = 0.25, linewidth = 0.5) +
  facet_wrap(~Wave, ncol = 2) +
  labs(
    x = "Age at First Menstruation (years)",
    y = "Mean Epigenetic Age Acceleration\n(GrimAge, Cell-Adjusted)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10),
    strip.text = element_text(size = 11, face = "bold"),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  geom_hline(yintercept = 0, linetype = "solid", color = "gray40", linewidth = 0.3)

print(p_manuscript_facet)

ggsave("output/Figure_Manuscript_GrimAge_Menarche_Facet.png", p_manuscript_facet, 
       width = 10, height = 5, dpi = 300)
ggsave("output/Figure_Manuscript_GrimAge_Menarche_Facet.pdf", p_manuscript_facet, 
       width = 10, height = 5)

# =============================================================================
# SUPPLEMENTARY FIGURE: All Clocks, Both Waves
# =============================================================================

cat("\n=== CREATING SUPPLEMENTARY FIGURES ===\n")

# Function for individual bar charts
create_bar_chart <- function(summary_data, clock_name, wave_name, fill_color) {
  
  plot_data <- summary_data %>% filter(Clock == clock_name, Wave == wave_name)
  wave_label <- ifelse(wave_name == "Age 9", "Wave 4", "Wave 5")
  
  p <- ggplot(plot_data, aes(x = factor(k6d10), y = mean_adj)) +
    geom_bar(stat = "identity", fill = fill_color, width = 0.7) +
    geom_errorbar(aes(ymin = mean_adj - se_adj, ymax = mean_adj + se_adj),
                  width = 0.25, linewidth = 0.5) +
    labs(
      title = paste(clock_name, "-", wave_label),
      x = "Age at First Menstruation",
      y = "Mean EAA (Cell-Adjusted)"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8)
    ) +
    geom_hline(yintercept = 0, linetype = "solid", color = "gray40", linewidth = 0.3)
  
  return(p)
}

# Colors
colors <- list(
  GrimAge = "#E8B8B8",
  PhenoAge = "#B8B8D1", 
  DunedinPACE = "#B8D1B8"
)

# Create all 6 panels
p_grim_w4 <- create_bar_chart(all_summary, "GrimAge", "Age 9", colors$GrimAge)
p_grim_w5 <- create_bar_chart(all_summary, "GrimAge", "Age 15", colors$GrimAge)
p_pheno_w4 <- create_bar_chart(all_summary, "PhenoAge", "Age 9", colors$PhenoAge)
p_pheno_w5 <- create_bar_chart(all_summary, "PhenoAge", "Age 15", colors$PhenoAge)
p_pace_w4 <- create_bar_chart(all_summary, "DunedinPACE", "Age 9", colors$DunedinPACE)
p_pace_w5 <- create_bar_chart(all_summary, "DunedinPACE", "Age 15", colors$DunedinPACE)

# Combine into supplementary figure (3 rows x 2 columns)
p_supplement <- (p_grim_w4 | p_grim_w5) /
                (p_pheno_w4 | p_pheno_w5) /
                (p_pace_w4 | p_pace_w5) +
  plot_annotation(
    title = "Epigenetic Age Acceleration by Age at Menarche",
    subtitle = "Cell composition adjusted. Menarche age 15 excluded (n=2).",
    caption = "Note: Sample sizes vary by age group. EAA = Epigenetic Age Acceleration.",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40"),
      plot.caption = element_text(size = 9, hjust = 0)
    )
  )

print(p_supplement)

ggsave("output/Figure_Supplement_All_Clocks_Menarche.png", p_supplement, 
       width = 10, height = 12, dpi = 300, bg = "white")
ggsave("output/Figure_Supplement_All_Clocks_Menarche.pdf", p_supplement, 
       width = 10, height = 12)

# =============================================================================
# Individual Supplementary Figures (one per clock)
# =============================================================================

# GrimAge both waves
p_grim_both <- (p_grim_w4 | p_grim_w5) +
  plot_annotation(caption = "Cell composition adjusted. Menarche age 15 excluded (n=2).")
ggsave("output/Figure_S_GrimAge_Menarche.png", p_grim_both, width = 10, height = 4, dpi = 300, bg = "white")

# PhenoAge both waves
p_pheno_both <- (p_pheno_w4 | p_pheno_w5) +
  plot_annotation(caption = "Cell composition adjusted. Menarche age 15 excluded (n=2).")
ggsave("output/Figure_S_PhenoAge_Menarche.png", p_pheno_both, width = 10, height = 4, dpi = 300, bg = "white")

# DunedinPACE both waves
p_pace_both <- (p_pace_w4 | p_pace_w5) +
  plot_annotation(caption = "Cell composition adjusted. Menarche age 15 excluded (n=2).")
ggsave("output/Figure_S_DunedinPACE_Menarche.png", p_pace_both, width = 10, height = 4, dpi = 300, bg = "white")

# =============================================================================
# Save Summary Data
# =============================================================================

write.csv(all_summary, "output/menarche_bioage_summary_all.csv", row.names = FALSE)

# Create correlation summary table
cor_summary <- data.frame(
  Clock = c("GrimAge", "PhenoAge", "DunedinPACE"),
  Wave4_r = c(cor_grim_w4$estimate, cor_pheno_w4$estimate, cor_pace_w4$estimate),
  Wave4_p = c(cor_grim_w4$p.value, cor_pheno_w4$p.value, cor_pace_w4$p.value),
  Wave4_n = c(nrow(grim_adj_w4), nrow(pheno_adj_w4), nrow(pace_adj_w4)),
  Wave5_r = c(cor_grim_w5$estimate, cor_pheno_w5$estimate, cor_pace_w5$estimate),
  Wave5_p = c(cor_grim_w5$p.value, cor_pheno_w5$p.value, cor_pace_w5$p.value),
  Wave5_n = c(nrow(grim_adj_w5), nrow(pheno_adj_w5), nrow(pace_adj_w5))
)

write.csv(cor_summary, "output/menarche_correlations.csv", row.names = FALSE)

cat("\n=== CORRELATION SUMMARY TABLE ===\n")
print(cor_summary)

# =============================================================================
# Output Summary
# =============================================================================

cat("\n=== FILES SAVED ===\n")
cat("\nManuscript Figures:\n")
cat("- output/Figure_Manuscript_GrimAge_Menarche.png/.pdf (grouped bars)\n")
cat("- output/Figure_Manuscript_GrimAge_Menarche_Facet.png/.pdf (faceted)\n")

cat("\nSupplementary Figures:\n")
cat("- output/Figure_Supplement_All_Clocks_Menarche.png/.pdf (all 6 panels)\n")
cat("- output/Figure_S_GrimAge_Menarche.png\n")
cat("- output/Figure_S_PhenoAge_Menarche.png\n")
cat("- output/Figure_S_DunedinPACE_Menarche.png\n")

cat("\nData Files:\n")
cat("- output/menarche_bioage_summary_all.csv\n")
cat("- output/menarche_correlations.csv\n")

# =============================================================================
# Figure Caption Text
# =============================================================================

cat("\n=== SUGGESTED FIGURE CAPTIONS ===\n")

cat("\nMANUSCRIPT FIGURE:\n")
cat("Figure X. Mean epigenetic age acceleration (GrimAge) by age at first menstruation.\n")
cat("Bars represent mean cell composition-adjusted GrimAge acceleration at Wave 4 (age 9)\n")
cat("and Wave 5 (age 15). Error bars indicate standard errors. Girls with earlier menarche\n")
cat("(age 10) show higher biological age acceleration compared to those with later menarche.\n")
cat(sprintf("Wave 4: r = %.3f, p = %.3f, n = %d; Wave 5: r = %.3f, p = %.3f, n = %d.\n",
            cor_grim_w4$estimate, cor_grim_w4$p.value, nrow(grim_adj_w4),
            cor_grim_w5$estimate, cor_grim_w5$p.value, nrow(grim_adj_w5)))
cat("Participants with menarche at age 15 excluded due to small sample size (n=2).\n")

cat("\nSUPPLEMENTARY FIGURE:\n")
cat("Figure SX. Mean epigenetic age acceleration by age at first menstruation across\n")
cat("three epigenetic clocks (GrimAge, PhenoAge, DunedinPACE) at Wave 4 (age 9) and\n")
cat("Wave 5 (age 15). Bars represent mean cell composition-adjusted values; error bars\n")
cat("indicate standard errors. Participants with menarche at age 15 excluded (n=2).\n")

cat("\n=== ANALYSIS COMPLETE ===\n")
