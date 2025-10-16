#!/usr/bin/env Rscript
# =============================================================================
# Create Publication-Ready Summary Tables
# =============================================================================

library(readr)
library(dplyr)
library(knitr)

# Paths
WORKSPACE_ROOT <- "/Users/kellywang/Documents/compSci/p1/greenLab/GreenLab"
OUT_DIR <- file.path(WORKSPACE_ROOT, "analysis", "outputs")

# Read the cleaned data
data <- read_csv(file.path(OUT_DIR, "cleaned_experiment_data.csv"), show_col_types = FALSE)

# Filter out CLEAN and ORIGINAL_FIXED variants
message("Filtering out CLEAN and ORIGINAL_FIXED variants...")
n_before <- nrow(data)
data <- data %>%
  filter(!variant %in% c("CLEAN", "ORIGINAL_FIXED"))
n_after <- nrow(data)
message("Removed ", n_before - n_after, " rows")
message("Remaining: ", n_after, " measurements")

message("\nCreating publication-ready summary tables...")

# =============================================================================
# Table 1: Overall Descriptive Statistics (Main Metrics)
# =============================================================================

table1 <- data %>%
  summarise(
    Metric = "All Data",
    N = n(),
    `Energy Mean (J)` = sprintf("%.2f", mean(energy_joules, na.rm = TRUE)),
    `Energy Median (J)` = sprintf("%.2f", median(energy_joules, na.rm = TRUE)),
    `Energy SD (J)` = sprintf("%.2f", sd(energy_joules, na.rm = TRUE)),
    `Energy IQR (J)` = sprintf("%.2f", IQR(energy_joules, na.rm = TRUE)),
    `Power Mean (W)` = sprintf("%.2f", mean(power_watts, na.rm = TRUE)),
    `Power Median (W)` = sprintf("%.2f", median(power_watts, na.rm = TRUE)),
    `CPU Mean (%)` = sprintf("%.2f", mean(cpu_usage_avg, na.rm = TRUE)),
    `Memory Mean (MB)` = sprintf("%.2f", mean(memory_usage_avg, na.rm = TRUE))
  )

write_csv(table1, file.path(OUT_DIR, "table1_overall_summary.csv"))
message("✓ Table 1: Overall Summary")

# =============================================================================
# Table 2: Summary by Variant (for RQ analysis)
# =============================================================================

table2 <- data %>%
  group_by(Variant = variant) %>%
  summarise(
    N = n(),
    `Energy Mean (J)` = sprintf("%.2f", mean(energy_joules, na.rm = TRUE)),
    `Energy Median (J)` = sprintf("%.2f", median(energy_joules, na.rm = TRUE)),
    `Energy SD (J)` = sprintf("%.2f", sd(energy_joules, na.rm = TRUE)),
    `Energy Min (J)` = sprintf("%.2f", min(energy_joules, na.rm = TRUE)),
    `Energy Max (J)` = sprintf("%.2f", max(energy_joules, na.rm = TRUE)),
    `Energy IQR (J)` = sprintf("%.2f", IQR(energy_joules, na.rm = TRUE)),
    `CV (%)` = sprintf("%.2f", (sd(energy_joules, na.rm = TRUE) / mean(energy_joules, na.rm = TRUE)) * 100),
    `Power Mean (W)` = sprintf("%.2f", mean(power_watts, na.rm = TRUE)),
    `CPU Mean (%)` = sprintf("%.2f", mean(cpu_usage_avg, na.rm = TRUE)),
    `Memory Mean (MB)` = sprintf("%.2f", mean(memory_usage_avg, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  arrange(as.numeric(`Energy Median (J)`))

write_csv(table2, file.path(OUT_DIR, "table2_summary_by_variant.csv"))
message("✓ Table 2: Summary by Variant")

# =============================================================================
# Table 3: Summary by Benchmark
# =============================================================================

table3 <- data %>%
  group_by(Benchmark = benchmark) %>%
  summarise(
    N = n(),
    `Variants Tested` = n_distinct(variant),
    `Energy Mean (J)` = sprintf("%.2f", mean(energy_joules, na.rm = TRUE)),
    `Energy Median (J)` = sprintf("%.2f", median(energy_joules, na.rm = TRUE)),
    `Energy SD (J)` = sprintf("%.2f", sd(energy_joules, na.rm = TRUE)),
    `Energy Range (J)` = sprintf("%.2f - %.2f", min(energy_joules, na.rm = TRUE), max(energy_joules, na.rm = TRUE)),
    `CV (%)` = sprintf("%.2f", (sd(energy_joules, na.rm = TRUE) / mean(energy_joules, na.rm = TRUE)) * 100),
    `Power Mean (W)` = sprintf("%.2f", mean(power_watts, na.rm = TRUE)),
    `CPU Mean (%)` = sprintf("%.2f", mean(cpu_usage_avg, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  arrange(desc(as.numeric(`Energy Median (J)`)))

write_csv(table3, file.path(OUT_DIR, "table3_summary_by_benchmark.csv"))
message("✓ Table 3: Summary by Benchmark")

# =============================================================================
# Table 4: Best and Worst Performers by Benchmark
# =============================================================================

baseline_energy <- data %>%
  filter(variant == "BASELINE") %>%
  group_by(benchmark) %>%
  summarise(baseline_median = median(energy_joules, na.rm = TRUE), .groups = "drop")

variant_performance <- data %>%
  filter(variant != "BASELINE") %>%
  group_by(benchmark, variant) %>%
  summarise(
    median_energy = median(energy_joules, na.rm = TRUE),
    mean_energy = mean(energy_joules, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(baseline_energy, by = "benchmark") %>%
  mutate(
    pct_change = ((median_energy - baseline_median) / baseline_median) * 100
  )

best_per_benchmark <- variant_performance %>%
  group_by(benchmark) %>%
  slice_min(median_energy, n = 1) %>%
  ungroup() %>%
  transmute(
    Benchmark = benchmark,
    `Best Variant` = variant,
    `Energy (J)` = sprintf("%.2f", median_energy),
    `vs Baseline (%)` = sprintf("%+.2f", pct_change)
  )

worst_per_benchmark <- variant_performance %>%
  group_by(benchmark) %>%
  slice_max(median_energy, n = 1) %>%
  ungroup() %>%
  transmute(
    Benchmark = benchmark,
    `Worst Variant` = variant,
    `Energy (J)` = sprintf("%.2f", median_energy),
    `vs Baseline (%)` = sprintf("%+.2f", pct_change)
  )

table4 <- best_per_benchmark %>%
  left_join(worst_per_benchmark, by = "Benchmark") %>%
  arrange(Benchmark)

write_csv(table4, file.path(OUT_DIR, "table4_best_worst_by_benchmark.csv"))
message("✓ Table 4: Best and Worst Performers")

# =============================================================================
# Table 5: Benchmark × Variant Matrix (Median Energy)
# =============================================================================

table5 <- data %>%
  group_by(benchmark, variant) %>%
  summarise(median_energy = median(energy_joules, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = variant, values_from = median_energy) %>%
  mutate(across(where(is.numeric), ~sprintf("%.2f", .)))

write_csv(table5, file.path(OUT_DIR, "table5_energy_matrix_median.csv"))
message("✓ Table 5: Energy Matrix (Median)")

# =============================================================================
# Table 6: Benchmark × Variant Matrix (% Change vs Baseline)
# =============================================================================

baseline_by_benchmark <- data %>%
  filter(variant == "BASELINE") %>%
  group_by(benchmark) %>%
  summarise(baseline_median = median(energy_joules, na.rm = TRUE), .groups = "drop")

table6 <- data %>%
  group_by(benchmark, variant) %>%
  summarise(median_energy = median(energy_joules, na.rm = TRUE), .groups = "drop") %>%
  left_join(baseline_by_benchmark, by = "benchmark") %>%
  mutate(pct_change = ((median_energy - baseline_median) / baseline_median) * 100) %>%
  select(benchmark, variant, pct_change) %>%
  tidyr::pivot_wider(names_from = variant, values_from = pct_change) %>%
  mutate(across(where(is.numeric), ~sprintf("%+.1f", .)))

write_csv(table6, file.path(OUT_DIR, "table6_pct_change_matrix.csv"))
message("✓ Table 6: % Change Matrix")

# =============================================================================
# Table 7: Data Quality Summary
# =============================================================================

table7 <- data %>%
  group_by(Benchmark = benchmark, Variant = variant) %>%
  summarise(
    N = n(),
    `Outliers (Moderate)` = sum(is_outlier_moderate, na.rm = TRUE),
    `Outliers (Extreme)` = sum(is_outlier_extreme, na.rm = TRUE),
    `Outlier Rate (%)` = sprintf("%.1f", (sum(is_outlier_moderate, na.rm = TRUE) / n()) * 100),
    `CV (%)` = sprintf("%.2f", (sd(energy_joules, na.rm = TRUE) / mean(energy_joules, na.rm = TRUE)) * 100),
    Stability = case_when(
      (sd(energy_joules, na.rm = TRUE) / mean(energy_joules, na.rm = TRUE)) * 100 < 10 ~ "Stable",
      (sd(energy_joules, na.rm = TRUE) / mean(energy_joules, na.rm = TRUE)) * 100 < 20 ~ "Moderate",
      (sd(energy_joules, na.rm = TRUE) / mean(energy_joules, na.rm = TRUE)) * 100 < 30 ~ "Variable",
      TRUE ~ "High Var"
    ),
    .groups = "drop"
  )

write_csv(table7, file.path(OUT_DIR, "table7_data_quality_matrix.csv"))
message("✓ Table 7: Data Quality Matrix")

# =============================================================================
# Table 8: Top 10 Most Energy-Efficient Combinations
# =============================================================================

table8 <- data %>%
  group_by(Benchmark = benchmark, Variant = variant) %>%
  summarise(
    N = n(),
    `Energy Median (J)` = sprintf("%.3f", median(energy_joules, na.rm = TRUE)),
    `Energy Mean (J)` = sprintf("%.3f", mean(energy_joules, na.rm = TRUE)),
    `SD (J)` = sprintf("%.3f", sd(energy_joules, na.rm = TRUE)),
    `Power (W)` = sprintf("%.2f", mean(power_watts, na.rm = TRUE)),
    `CPU (%)` = sprintf("%.1f", mean(cpu_usage_avg, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  arrange(as.numeric(`Energy Median (J)`)) %>%
  head(10) %>%
  mutate(Rank = row_number()) %>%
  select(Rank, everything())

write_csv(table8, file.path(OUT_DIR, "table8_top10_efficient.csv"))
message("✓ Table 8: Top 10 Most Efficient")

# =============================================================================
# Table 9: Top 10 Least Energy-Efficient Combinations
# =============================================================================

table9 <- data %>%
  group_by(Benchmark = benchmark, Variant = variant) %>%
  summarise(
    N = n(),
    `Energy Median (J)` = sprintf("%.2f", median(energy_joules, na.rm = TRUE)),
    `Energy Mean (J)` = sprintf("%.2f", mean(energy_joules, na.rm = TRUE)),
    `SD (J)` = sprintf("%.2f", sd(energy_joules, na.rm = TRUE)),
    `Power (W)` = sprintf("%.2f", mean(power_watts, na.rm = TRUE)),
    `CPU (%)` = sprintf("%.1f", mean(cpu_usage_avg, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  arrange(desc(as.numeric(`Energy Median (J)`))) %>%
  head(10) %>%
  mutate(Rank = row_number()) %>%
  select(Rank, everything())

write_csv(table9, file.path(OUT_DIR, "table9_top10_inefficient.csv"))
message("✓ Table 9: Top 10 Most Inefficient")

# =============================================================================
# Table 10: Guideline Effectiveness Summary
# =============================================================================

baseline_overall <- data %>%
  filter(variant == "BASELINE") %>%
  summarise(
    baseline_energy_median = median(energy_joules, na.rm = TRUE),
    baseline_power_median = median(power_watts, na.rm = TRUE)
  )

table10 <- data %>%
  filter(variant %in% c("G1", "G3", "G4", "G6", "G7", "G9", "G12", "G14")) %>%
  group_by(Guideline = variant) %>%
  summarise(
    `N Runs` = n(),
    `N Benchmarks` = n_distinct(benchmark),
    `Energy Median (J)` = sprintf("%.2f", median(energy_joules, na.rm = TRUE)),
    `Energy Mean (J)` = sprintf("%.2f", mean(energy_joules, na.rm = TRUE)),
    `vs Baseline (%)` = sprintf("%+.1f", 
      ((median(energy_joules, na.rm = TRUE) - baseline_overall$baseline_energy_median) / 
       baseline_overall$baseline_energy_median) * 100),
    `Power (W)` = sprintf("%.2f", median(power_watts, na.rm = TRUE)),
    `CPU (%)` = sprintf("%.1f", mean(cpu_usage_avg, na.rm = TRUE)),
    `CV (%)` = sprintf("%.1f", (sd(energy_joules, na.rm = TRUE) / mean(energy_joules, na.rm = TRUE)) * 100),
    .groups = "drop"
  ) %>%
  arrange(as.numeric(`Energy Median (J)`))

write_csv(table10, file.path(OUT_DIR, "table10_guideline_effectiveness.csv"))
message("✓ Table 10: Guideline Effectiveness")

# =============================================================================
# Summary
# =============================================================================

message("\n=============================================================================")
message("All summary tables created successfully!")
message("=============================================================================")
message("\nTables saved to: ", OUT_DIR)
message("\nTable descriptions:")
message("  1. Overall Summary - High-level statistics across all data")
message("  2. Summary by Variant - Performance metrics for each variant")
message("  3. Summary by Benchmark - Performance metrics for each benchmark")
message("  4. Best/Worst by Benchmark - Top and bottom performers per benchmark")
message("  5. Energy Matrix (Median) - Pivot table of median energy values")
message("  6. % Change Matrix - Pivot table of % change vs baseline")
message("  7. Data Quality Matrix - Outliers and stability by combination")
message("  8. Top 10 Efficient - Most energy-efficient combinations")
message("  9. Top 10 Inefficient - Least energy-efficient combinations")
message(" 10. Guideline Effectiveness - Summary of optimization guidelines")
message("\n=============================================================================\n")

