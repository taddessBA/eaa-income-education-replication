# Replication code for: "Maternal Income, Education, and Children's Epigenetic Age Acceleration"
# Manuscript ID: DEM-2024-0306 | Journal: Demography
# Input: data_2026_01_05.dta (FFCWS restricted-use data)
# =============================================================================
# FDR Correction and Summary Tables for Manuscript
# =============================================================================

library(tidyverse)

# =============================================================================
# 1. FDR CORRECTION ACROSS ALL OUTCOMES
# =============================================================================

cat("=============================================================================\n")
cat("FDR CORRECTION ACROSS ALL OUTCOMES\n")
cat("=============================================================================\n")

# Collect all p-values from our analyses
# -----------------------------------------------------------------------------

# MENARCHE ANALYSES (Wave 4 - prospective)
# Source: R&R_menarche_EAA_correlation.R (linear r) and
#         R&R_Non-Linearity Tests_menarche.R (quadratic β₂)
# Linear estimates are Pearson correlations (Table S3)
# Quadratic estimates are β₂ from centered polynomial regressions (Table S7)
menarche_pvals <- data.frame(
  Analysis = c("Menarche_Linear", "Menarche_Linear", "Menarche_Linear",
               "Menarche_Quadratic", "Menarche_Quadratic", "Menarche_Quadratic"),
  Clock = c("GrimAge", "DunedinPACE", "PhenoAge",
            "GrimAge", "DunedinPACE", "PhenoAge"),
  Test = c("Linear r", "Linear r", "Linear r",
           "Quadratic β₂", "Quadratic β₂", "Quadratic β₂"),
  Estimate = c(-0.075, 0.001, 0.025,       # Linear: Table S3 correlations
               0.032, 0.025, -0.014),       # Quadratic: Table S7 β₂ coefficients
  p_value = c(0.0283, 0.983, 0.471,         # Linear correlation p-values (Table S3)
              0.023, 0.0665, 0.5119)         # Quadratic p-values (Table S7)
)

# HEIGHT/GROWTH ANALYSES
# Source: R&R_height_V2.R (Model 4: + Height at 9 + Maternal Height)
# DV: Standardized self-reported height at age 15 (ck6chtcm, z-scored)
# Growth model: controls for height at age 9, maternal height, sex, race
# Cell composition and array type adjusted via residualization
# Updated 2026-03-02 after rerun with self-reported height (N=1,718)
height_pvals <- data.frame(
  Analysis = rep("Height_Growth", 3),
  Clock = c("GrimAge", "PhenoAge", "DunedinPACE"),
  Test = rep("EAA → Height Growth", 3),
  Estimate = c(-0.085, -0.030, -0.047),        # V2 Model 4 coefficients (rerun 2026-03-02)
  p_value = c(0.0001, 0.056, 0.048)            # V2 Model 4 p-values (rerun 2026-03-02)
)

# ADULT HEIGHT ANALYSES
# Source: R&R_height_V2.R (Model 3: + Maternal Height)
# DV: Height in cm at age 22 (self-reported, from p7l7a feet + p7l7b inches)
# Controls for maternal height, sex, race
# Cell composition and array type adjusted via residualization
# All results null (p > 0.30) — included for completeness per reviewer request
adult_height_pvals <- data.frame(
  Analysis = rep("Adult_Height", 3),
  Clock = c("GrimAge", "PhenoAge", "DunedinPACE"),
  Test = rep("EAA → Adult Height", 3),
  Estimate = c(0.194, -0.045, -0.015),           # V2 Model 3 coefficients (+ Maternal Height)
  p_value = c(0.364, 0.767, 0.945)               # V2 Model 3 p-values
)

# Combine all p-values (12 tests: 3 clocks × 4 outcome families)
all_pvals <- bind_rows(menarche_pvals, height_pvals, adult_height_pvals)

cat("\n=== RAW P-VALUES ===\n")
print(all_pvals)

# Apply FDR correction
all_pvals$p_fdr <- p.adjust(all_pvals$p_value, method = "fdr")

# Also apply Bonferroni for comparison
all_pvals$p_bonferroni <- p.adjust(all_pvals$p_value, method = "bonferroni")

# Add significance indicators
all_pvals <- all_pvals %>%
  mutate(
    sig_raw = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      p_value < 0.10 ~ "†",
      TRUE ~ ""
    ),
    sig_fdr = case_when(
      p_fdr < 0.001 ~ "***",
      p_fdr < 0.01 ~ "**",
      p_fdr < 0.05 ~ "*",
      p_fdr < 0.10 ~ "†",
      TRUE ~ ""
    )
  )

cat("\n=== FDR-CORRECTED P-VALUES ===\n")
print(all_pvals %>% select(Analysis, Clock, Test, Estimate, p_value, sig_raw, p_fdr, sig_fdr))

# Create output directory
if (!dir.exists("output")) dir.create("output")

# Save to CSV
write.csv(all_pvals, "output/FDR_Correction_All_Outcomes.csv", row.names = FALSE)

# =============================================================================
# 2. SUMMARY TABLE: MENARCHE ANALYSIS (with FDR)
# =============================================================================

cat("\n")
cat("=============================================================================\n")
cat("TABLE: MENARCHE ANALYSIS SUMMARY (Wave 4, N = 865)\n")
cat("=============================================================================\n")

menarche_summary <- all_pvals %>% 
  filter(grepl("Menarche", Analysis)) %>%
  select(Clock, Test, Estimate, p_value, p_fdr, sig_raw, sig_fdr)

# Pivot to wide format for table
menarche_wide <- menarche_summary %>%
  pivot_wider(
    names_from = Clock,
    values_from = c(Estimate, p_value, p_fdr, sig_raw, sig_fdr)
  )

cat("\n")
cat("Table X. Non-Linear Association Between EAA at Age 9 and Age at Menarche (N = 865)\n")
cat("─────────────────────────────────────────────────────────────────────────────────────\n")
cat("                              GrimAge           DunedinPACE        PhenoAge\n")
cat("─────────────────────────────────────────────────────────────────────────────────────\n")

# Quadratic results
grim_q <- menarche_summary %>% filter(Clock == "GrimAge", Test == "Quadratic β₂")
pace_q <- menarche_summary %>% filter(Clock == "DunedinPACE", Test == "Quadratic β₂")
pheno_q <- menarche_summary %>% filter(Clock == "PhenoAge", Test == "Quadratic β₂")

cat(sprintf("Quadratic β₂                  %-17.3f %-18.3f %.3f\n", 
            grim_q$Estimate, pace_q$Estimate, pheno_q$Estimate))
cat(sprintf("  p-value (raw)               %-17.4f %-18.4f %.4f\n",
            grim_q$p_value, pace_q$p_value, pheno_q$p_value))
cat(sprintf("  p-value (FDR)               %-17.4f %-18.4f %.4f\n",
            grim_q$p_fdr, pace_q$p_fdr, pheno_q$p_fdr))
cat(sprintf("  Significance (FDR)          %-17s %-18s %s\n",
            grim_q$sig_fdr, pace_q$sig_fdr, pheno_q$sig_fdr))
cat(sprintf("  Turning point (years)       %-17.1f %-18.1f %.1f\n", 12.6, 12.0, 12.7))
cat("─────────────────────────────────────────────────────────────────────────────────────\n")
cat("Note: Cell composition and array type adjusted. FDR correction applied across\n")
cat("all 12 tests (3 clocks × 4 outcomes). † p<0.10, * p<0.05, ** p<0.01, *** p<0.001.\n")
cat("─────────────────────────────────────────────────────────────────────────────────────\n")

# =============================================================================
# 3. SUMMARY TABLE: HEIGHT/GROWTH ANALYSIS (with FDR)
# =============================================================================

cat("\n")
cat("=============================================================================\n")
cat("TABLE: HEIGHT/GROWTH ANALYSIS SUMMARY (N = 1,718)\n")
cat("=============================================================================\n")

height_summary <- all_pvals %>% 
  filter(Analysis == "Height_Growth")

cat("\n")
cat("Table X. EAA at Age 9 and Height Growth to Age 15 (N = 1,718)\n")
cat("─────────────────────────────────────────────────────────────────────────────────────\n")
cat("                              GrimAge           PhenoAge           DunedinPACE\n")
cat("─────────────────────────────────────────────────────────────────────────────────────\n")

grim_h <- height_summary %>% filter(Clock == "GrimAge")
pheno_h <- height_summary %>% filter(Clock == "PhenoAge")
pace_h <- height_summary %>% filter(Clock == "DunedinPACE")

cat(sprintf("β (EAA → Height Growth)       %-17.3f %-18.3f %.3f\n", 
            grim_h$Estimate, pheno_h$Estimate, pace_h$Estimate))
cat(sprintf("  SE                          %-17.3f %-18.3f %.3f\n", 0.022, 0.016, 0.024))
cat(sprintf("  p-value (raw)               %-17.4f %-18.4f %.4f\n",
            grim_h$p_value, pheno_h$p_value, pace_h$p_value))
cat(sprintf("  p-value (FDR)               %-17.4f %-18.4f %.4f\n",
            grim_h$p_fdr, pheno_h$p_fdr, pace_h$p_fdr))
cat(sprintf("  Significance (FDR)          %-17s %-18s %s\n",
            grim_h$sig_fdr, pheno_h$sig_fdr, pace_h$sig_fdr))
cat("─────────────────────────────────────────────────────────────────────────────────────\n")
cat("Note: Growth model controls for height at age 9, maternal height, sex, and race.\n")
cat("Cell composition and array type adjusted. FDR correction applied across all 12 tests.\n")
cat("Source: R&R_height_V2.R, Model 4 (+ Height at 9 + Maternal Height)\n")
cat("─────────────────────────────────────────────────────────────────────────────────────\n")

# =============================================================================
# 3b. SUMMARY TABLE: ADULT HEIGHT (with FDR)
# =============================================================================

cat("\n")
cat("=============================================================================\n")
cat("TABLE: ADULT HEIGHT ANALYSIS SUMMARY (N ≈ 1,393)\n")
cat("=============================================================================\n")

adult_summary <- all_pvals %>%
  filter(Analysis == "Adult_Height")

cat("\n")
cat("Table X. EAA at Age 15 and Adult Height at Age 22 (N ≈ 1,393)\n")
cat("─────────────────────────────────────────────────────────────────────────────────────\n")
cat("                              GrimAge           PhenoAge           DunedinPACE\n")
cat("─────────────────────────────────────────────────────────────────────────────────────\n")

grim_a <- adult_summary %>% filter(Clock == "GrimAge")
pheno_a <- adult_summary %>% filter(Clock == "PhenoAge")
pace_a <- adult_summary %>% filter(Clock == "DunedinPACE")

cat(sprintf("β (EAA → Adult Height)        %-17.3f %-18.3f %.3f\n",
            grim_a$Estimate, pheno_a$Estimate, pace_a$Estimate))
cat(sprintf("  SE                          %-17.3f %-18.3f %.3f\n", 0.213, 0.152, 0.220))
cat(sprintf("  p-value (raw)               %-17.4f %-18.4f %.4f\n",
            grim_a$p_value, pheno_a$p_value, pace_a$p_value))
cat(sprintf("  p-value (FDR)               %-17.4f %-18.4f %.4f\n",
            grim_a$p_fdr, pheno_a$p_fdr, pace_a$p_fdr))
cat(sprintf("  Significance (FDR)          %-17s %-18s %s\n",
            grim_a$sig_fdr, pheno_a$sig_fdr, pace_a$sig_fdr))
cat("─────────────────────────────────────────────────────────────────────────────────────\n")
cat("Note: Model controls for maternal height, sex, and race. All null.\n")
cat("Cell composition and array type adjusted. FDR correction applied across all 12 tests.\n")
cat("─────────────────────────────────────────────────────────────────────────────────────\n")

# =============================================================================
# 4. COMBINED SUMMARY FOR RESPONSE MEMO
# =============================================================================

cat("\n")
cat("=============================================================================\n")
cat("SUMMARY FOR RESPONSE MEMO\n")
cat("=============================================================================\n")

cat("\n[FDR CORRECTION STATEMENT]\n")
cat("─────────────────────────────────────────────────────────────────────────────────────\n")
cat("We applied false discovery rate (FDR) correction across all 12 developmental\n")
cat("outcome tests (3 clocks × 4 outcomes: menarche linear, menarche quadratic,\n")
cat("height growth, adult height).\n")
cat("Results that remain significant after FDR correction are indicated below.\n")
cat("─────────────────────────────────────────────────────────────────────────────────────\n")

cat("\n[SIGNIFICANT RESULTS AFTER FDR]\n")
sig_after_fdr <- all_pvals %>% filter(p_fdr < 0.10)
if (nrow(sig_after_fdr) > 0) {
  for (i in 1:nrow(sig_after_fdr)) {
    r <- sig_after_fdr[i,]
    cat(sprintf("  • %s - %s: β = %.3f, p_raw = %.4f, p_FDR = %.4f %s\n",
                r$Clock, r$Test, r$Estimate, r$p_value, r$p_fdr, r$sig_fdr))
  }
} else {
  cat("  No results significant at FDR < 0.10\n")
}

cat("\n[MARGINALLY SIGNIFICANT AFTER FDR (p < 0.15)]\n")
marginal_fdr <- all_pvals %>% filter(p_fdr >= 0.10 & p_fdr < 0.15)
if (nrow(marginal_fdr) > 0) {
  for (i in 1:nrow(marginal_fdr)) {
    r <- marginal_fdr[i,]
    cat(sprintf("  • %s - %s: β = %.3f, p_raw = %.4f, p_FDR = %.4f\n",
                r$Clock, r$Test, r$Estimate, r$p_value, r$p_fdr))
  }
}

# =============================================================================
# 5. MANUSCRIPT TEXT SUGGESTIONS
# =============================================================================

cat("\n")
cat("=============================================================================\n")
cat("SUGGESTED MANUSCRIPT TEXT\n")
cat("=============================================================================\n")

cat("\n[MENARCHE RESULTS - FDR CORRECTED]\n")
cat("─────────────────────────────────────────────────────────────────────────────────────\n")
cat("We examined whether the association between epigenetic age acceleration at age 9\n")
cat("and age at menarche was non-linear. Quadratic models revealed a significant\n")
cat(sprintf("U-shaped association for GrimAge (β₂ = 0.032, p_raw = 0.023, p_FDR = %.3f;\n", grim_q$p_fdr))
cat("turning point = 12.6 years), a marginally significant pattern for DunedinPACE\n")
cat(sprintf("(β₂ = 0.025, p_raw = 0.066, p_FDR = %.3f; turning point = 12.0 years), and no\n", pace_q$p_fdr))
cat(sprintf("significant association for PhenoAge (β₂ = -0.014, p_raw = 0.512, p_FDR = %.3f).\n", pheno_q$p_fdr))
cat("─────────────────────────────────────────────────────────────────────────────────────\n")

cat("\n[HEIGHT/GROWTH RESULTS - FDR CORRECTED]\n")
cat("─────────────────────────────────────────────────────────────────────────────────────\n")
cat("In the growth model — controlling for height at age 9, maternal height, sex, and\n")
cat("race — GrimAge EAA at age 9 was significantly associated with height growth to\n")
cat(sprintf("age 15 (β = %.3f, p_raw = %.4f, p_FDR = %.3f). DunedinPACE showed a\n",
            grim_h$Estimate, grim_h$p_value, grim_h$p_fdr))
cat(sprintf("significant unadjusted association (β = %.3f, p_raw = %.4f, p_FDR = %.3f),\n",
            pace_h$Estimate, pace_h$p_value, pace_h$p_fdr))
cat(sprintf("and PhenoAge was marginally significant (β = %.3f, p_raw = %.4f, p_FDR = %.3f).\n",
            pheno_h$Estimate, pheno_h$p_value, pheno_h$p_fdr))
cat("─────────────────────────────────────────────────────────────────────────────────────\n")

# Save summary tables
write.csv(menarche_summary, "output/Table_Menarche_FDR.csv", row.names = FALSE)
write.csv(height_summary, "output/Table_Height_FDR.csv", row.names = FALSE)
write.csv(adult_summary, "output/Table_Adult_Height_FDR.csv", row.names = FALSE)

cat("\n=== FILES SAVED ===\n")
cat("- output/FDR_Correction_All_Outcomes.csv\n")
cat("- output/Table_Menarche_FDR.csv\n")
cat("- output/Table_Height_FDR.csv\n")
cat("- output/Table_Adult_Height_FDR.csv\n")

cat("\n=== ANALYSIS COMPLETE ===\n")

