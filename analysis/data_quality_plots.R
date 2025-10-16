#!/usr/bin/env Rscript

# Data Quality Visualization Plots
# Load required packages
packages <- c("ggplot2", "readr", "dplyr", "stringr", "tidyr", "scales", "forcats")
for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}
suppressPackageStartupMessages({
  invisible(lapply(packages, function(p) library(p, character.only = TRUE)))
})

# Read data
OUTPUT_DIR <- "/Users/kellywang/Desktop/GreenLab/analysis"
PLOTS_DIR <- file.path(OUTPUT_DIR, "plots")
OUT_DIR <- file.path(OUTPUT_DIR, "outputs")

combined_cleaned <- read_csv(file.path(OUT_DIR, "combined_runs_cleaned.csv"), show_col_types = FALSE)
data_quality <- read_csv(file.path(OUT_DIR, "data_quality_summary.csv"), show_col_types = FALSE)
outliers_data <- read_csv(file.path(OUT_DIR, "outliers_detected.csv"), show_col_types = FALSE)

message("Generating data quality plots...")

# Helper function
safe_save <- function(plot, filename, width = 10, height = 6) {
  out_path <- file.path(PLOTS_DIR, filename)
  tryCatch({
    ggsave(out_path, plot = plot, width = width, height = height, dpi = 150)
    message("Saved ", out_path)
  }, error = function(e) warning("Failed to save ", out_path, ": ", conditionMessage(e)))
}

# 1) Outlier counts by benchmark-variant
p_outliers <- ggplot(data_quality, aes(x = reorder(paste(benchmark, variant, sep = " - "), n_outliers), 
                                        y = n_outliers)) +
  geom_col(aes(fill = outlier_rate), width = 0.7) +
  scale_fill_gradient(low = "#fee5d9", high = "#a50f15", name = "Outlier\nRate (%)") +
  coord_flip() +
  labs(
    title = "Outlier Detection by Benchmark-Variant",
    subtitle = "Using 3×IQR method (extreme outliers only)",
    x = "Benchmark - Variant",
    y = "Number of Outliers"
  ) +
  theme_minimal(base_size = 10) +
  theme(axis.text.y = element_text(size = 7))

safe_save(p_outliers, "quality_outliers_by_variant.png", width = 11, height = 12)

# 2) Coefficient of Variation heatmap (consistency metric)
p_cv_heat <- ggplot(data_quality, aes(x = variant, y = benchmark, fill = energy_cv)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.1f", energy_cv)), size = 2.5, color = "black") +
  scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b", 
                       midpoint = median(data_quality$energy_cv, na.rm = TRUE),
                       name = "CV (%)") +
  labs(
    title = "Measurement Consistency: Coefficient of Variation",
    subtitle = "Lower CV = more consistent measurements across runs",
    x = "Variant",
    y = "Benchmark"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

safe_save(p_cv_heat, "quality_cv_heatmap.png", width = 12, height = 8)

# 3) Energy distribution with outliers flagged
p_energy_dist <- ggplot(combined_cleaned, aes(x = energy_joules, fill = is_outlier)) +
  geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
  scale_fill_manual(values = c("FALSE" = "#4daf4a", "TRUE" = "#e41a1c"),
                    labels = c("Normal", "Outlier"), name = NULL) +
  scale_x_log10(labels = label_number(scale_cut = cut_short_scale(), suffix = "J")) +
  labs(
    title = "Energy Distribution with Outliers Flagged",
    subtitle = "Log scale for better visualization",
    x = "Energy (J, log scale)",
    y = "Count"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

safe_save(p_energy_dist, "quality_energy_distribution.png", width = 10, height = 6)

# 4) Z-score distribution (check for normality)
p_zscore <- ggplot(combined_cleaned, aes(x = energy_z_score)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40, fill = "#4daf4a", alpha = 0.7) +
  geom_density(color = "#e41a1c", linewidth = 1) +
  geom_vline(xintercept = c(-3, 3), linetype = "dashed", color = "blue") +
  annotate("text", x = -3, y = 0.3, label = "-3σ", hjust = 1.2, color = "blue") +
  annotate("text", x = 3, y = 0.3, label = "+3σ", hjust = -0.2, color = "blue") +
  labs(
    title = "Z-Score Distribution of Energy Values",
    subtitle = "Should approximate normal distribution if data is well-behaved",
    x = "Z-Score",
    y = "Density"
  ) +
  theme_minimal(base_size = 12)

safe_save(p_zscore, "quality_zscore_distribution.png", width = 10, height = 6)

# 5) Normalized energy comparison (BASELINE vs OPTIMIZED)
norm_comparison <- combined_cleaned %>%
  filter(variant %in% c("BASELINE", "OPTIMIZED")) %>%
  select(benchmark, variant, energy_normalized, results_set)

p_norm <- ggplot(norm_comparison, aes(x = benchmark, y = energy_normalized, fill = variant)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = c("BASELINE" = "#377eb8", "OPTIMIZED" = "#e41a1c")) +
  coord_flip() +
  labs(
    title = "Normalized Energy: BASELINE vs OPTIMIZED",
    subtitle = "0-1 scale normalized per benchmark (0=best, 1=worst)",
    x = "Benchmark",
    y = "Normalized Energy (0-1)",
    fill = "Variant"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

safe_save(p_norm, "quality_normalized_comparison.png", width = 10, height = 8)

# 6) Efficiency ratio distribution
eff_data <- combined_cleaned %>%
  filter(!is.na(efficiency_ratio), is.finite(efficiency_ratio), variant != "BASELINE")

p_efficiency <- ggplot(eff_data, aes(x = variant, y = efficiency_ratio, fill = variant)) +
  geom_violin(alpha = 0.6, trim = TRUE) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 1) +
  annotate("text", x = 1, y = 1.05, label = "Baseline", color = "red", hjust = 0) +
  scale_fill_brewer(palette = "Set3", guide = "none") +
  labs(
    title = "Efficiency Ratio Distribution by Variant",
    subtitle = "Ratio = Energy / Baseline Energy (< 1 = improvement, > 1 = regression)",
    x = "Variant",
    y = "Efficiency Ratio"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

safe_save(p_efficiency, "quality_efficiency_ratio.png", width = 11, height = 7)

# 7) Benchmark stability (average CV per benchmark)
benchmark_stability <- data_quality %>%
  group_by(benchmark) %>%
  summarise(
    avg_cv = mean(energy_cv, na.rm = TRUE),
    max_cv = max(energy_cv, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(avg_cv))

p_stability <- ggplot(benchmark_stability, aes(x = reorder(benchmark, avg_cv), y = avg_cv)) +
  geom_col(aes(fill = avg_cv), width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", avg_cv)), hjust = -0.2, size = 3.5) +
  scale_fill_gradient(low = "#4daf4a", high = "#e41a1c", guide = "none") +
  coord_flip() +
  labs(
    title = "Benchmark Stability: Average CV Across Variants",
    subtitle = "Lower = more stable/reproducible measurements",
    x = "Benchmark",
    y = "Average Coefficient of Variation (%)"
  ) +
  theme_minimal(base_size = 12) +
  ylim(0, max(benchmark_stability$avg_cv) * 1.15)

safe_save(p_stability, "quality_benchmark_stability.png", width = 10, height = 7)

# 8) Outlier impact on mean values
mean_comparison <- combined_cleaned %>%
  group_by(benchmark, variant) %>%
  summarise(
    mean_with_outliers = mean(energy_joules, na.rm = TRUE),
    mean_without_outliers = mean(energy_joules[!is_outlier], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    pct_diff = (mean_with_outliers - mean_without_outliers) / mean_without_outliers * 100,
    has_impact = abs(pct_diff) > 5  # Flag if outliers change mean by >5%
  ) %>%
  arrange(desc(abs(pct_diff))) %>%
  head(30)

p_outlier_impact <- ggplot(mean_comparison, 
                            aes(x = reorder(paste(benchmark, variant, sep = " - "), abs(pct_diff)), 
                                y = pct_diff)) +
  geom_col(aes(fill = has_impact), width = 0.7) +
  scale_fill_manual(values = c("TRUE" = "#e41a1c", "FALSE" = "#4daf4a"),
                    labels = c("< 5% impact", "> 5% impact"), name = "Outlier Impact") +
  coord_flip() +
  labs(
    title = "Top 30: Impact of Outliers on Mean Energy",
    subtitle = "% difference in mean when outliers are included vs excluded",
    x = "Benchmark - Variant",
    y = "% Change in Mean"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.y = element_text(size = 7),
    legend.position = "top"
  )

safe_save(p_outlier_impact, "quality_outlier_impact.png", width = 11, height = 10)

message("Done! All quality plots generated in: ", PLOTS_DIR)
message("Summary:")
message("  - Total outliers detected: ", nrow(outliers_data))
message("  - Benchmarks analyzed: ", length(unique(data_quality$benchmark)))
message("  - Variants analyzed: ", length(unique(data_quality$variant)))
message("  - Average CV across all: ", sprintf("%.2f%%", mean(data_quality$energy_cv, na.rm = TRUE)))

