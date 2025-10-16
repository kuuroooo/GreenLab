#!/usr/bin/env Rscript

# Additional plots for energy analysis
# Load required packages
packages <- c("ggplot2", "readr", "dplyr", "stringr", "tidyr", "scales", "forcats", "RColorBrewer")
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

combined_runs <- read_csv(file.path(OUT_DIR, "combined_runs.csv"), show_col_types = FALSE)
summary_by_variant <- read_csv(file.path(OUT_DIR, "summary_by_benchmark_variant.csv"), show_col_types = FALSE)

message("Generating additional plots...")

# Helper function
safe_save <- function(plot, filename, width = 10, height = 6) {
  out_path <- file.path(PLOTS_DIR, filename)
  tryCatch({
    ggsave(out_path, plot = plot, width = width, height = height, dpi = 150)
    message("Saved ", out_path)
  }, error = function(e) warning("Failed to save ", out_path, ": ", conditionMessage(e)))
}

# 1) Radar/Spider chart of normalized metrics per guideline (avg across all benchmarks)
guideline_summary <- summary_by_variant %>%
  filter(variant != "BASELINE") %>%
  group_by(variant) %>%
  summarise(
    avg_pct_vs_baseline = mean(pct_vs_baseline, na.rm = TRUE),
    median_pct_vs_baseline = median(pct_vs_baseline, na.rm = TRUE),
    best_case = min(pct_vs_baseline, na.rm = TRUE),
    worst_case = max(pct_vs_baseline, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(avg_pct_vs_baseline)

p_guideline_effectiveness <- ggplot(guideline_summary, aes(x = reorder(variant, avg_pct_vs_baseline), y = avg_pct_vs_baseline)) +
  geom_col(aes(fill = avg_pct_vs_baseline), width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
  scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b", midpoint = 0, guide = "none") +
  coord_flip() +
  labs(
    title = "Guideline Effectiveness: Average % Change vs Baseline",
    subtitle = "Averaged across all benchmarks",
    x = "Variant/Guideline",
    y = "Average % Change in Energy"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

safe_save(p_guideline_effectiveness, "guideline_effectiveness_avg.png", width = 10, height = 7)

# 2) Best vs Worst performers: which guidelines work best/worst per benchmark
best_worst <- summary_by_variant %>%
  filter(variant != "BASELINE") %>%
  group_by(benchmark) %>%
  slice_min(pct_vs_baseline, n = 1, with_ties = FALSE) %>%
  mutate(type = "Best") %>%
  bind_rows(
    summary_by_variant %>%
      filter(variant != "BASELINE") %>%
      group_by(benchmark) %>%
      slice_max(pct_vs_baseline, n = 1, with_ties = FALSE) %>%
      mutate(type = "Worst")
  ) %>%
  ungroup()

p_best_worst <- ggplot(best_worst, aes(x = reorder(benchmark, pct_vs_baseline), y = pct_vs_baseline, fill = type)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("Best" = "#2ca02c", "Worst" = "#d62728")) +
  coord_flip() +
  labs(
    title = "Best vs Worst Performing Guidelines per Benchmark",
    x = "Benchmark",
    y = "% Change vs Baseline",
    fill = "Performance"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top")

safe_save(p_best_worst, "best_worst_guidelines_per_benchmark.png", width = 11, height = 8)

# 3) Variability/Consistency: CV (coefficient of variation) of energy across runs
variability <- combined_runs %>%
  group_by(benchmark, variant) %>%
  summarise(
    cv = sd(energy_joules, na.rm = TRUE) / mean(energy_joules, na.rm = TRUE) * 100,
    mean_energy = mean(energy_joules, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(cv), is.finite(cv))

p_variability <- ggplot(variability, aes(x = variant, y = cv, fill = variant)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.7) +
  scale_fill_brewer(palette = "Set3", guide = "none") +
  labs(
    title = "Energy Measurement Variability by Variant",
    subtitle = "Coefficient of Variation (CV) across 10 runs",
    x = "Variant",
    y = "CV (%)"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

safe_save(p_variability, "variability_by_variant.png", width = 10, height = 6)

# 4) Energy savings distribution: histogram of % changes
energy_changes <- summary_by_variant %>%
  filter(variant != "BASELINE") %>%
  select(benchmark, variant, pct_vs_baseline)

p_distribution <- ggplot(energy_changes, aes(x = pct_vs_baseline)) +
  geom_histogram(aes(fill = after_stat(x)), bins = 30, color = "white", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
  scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b", midpoint = 0, guide = "none") +
  labs(
    title = "Distribution of Energy Changes",
    subtitle = "% change vs baseline across all benchmark-variant combinations",
    x = "% Change vs Baseline",
    y = "Count"
  ) +
  theme_minimal(base_size = 12)

safe_save(p_distribution, "energy_change_distribution.png", width = 10, height = 6)

# 5) Pareto chart: cumulative energy savings by guideline
pareto_data <- guideline_summary %>%
  arrange(avg_pct_vs_baseline) %>%
  mutate(
    cumulative = cumsum(abs(avg_pct_vs_baseline)),
    cumulative_pct = cumulative / sum(abs(avg_pct_vs_baseline)) * 100
  )

p_pareto <- ggplot(pareto_data, aes(x = reorder(variant, avg_pct_vs_baseline))) +
  geom_col(aes(y = abs(avg_pct_vs_baseline), fill = avg_pct_vs_baseline < 0), width = 0.7) +
  geom_line(aes(y = cumulative_pct * max(abs(avg_pct_vs_baseline)) / 100, group = 1), 
            color = "red", linewidth = 1.2) +
  geom_point(aes(y = cumulative_pct * max(abs(avg_pct_vs_baseline)) / 100), 
             color = "red", size = 3) +
  scale_y_continuous(
    name = "Absolute % Change",
    sec.axis = sec_axis(~ . * 100 / max(abs(pareto_data$avg_pct_vs_baseline)), name = "Cumulative %")
  ) +
  scale_fill_manual(values = c("TRUE" = "#2ca02c", "FALSE" = "#d62728"), 
                    labels = c("Increase", "Decrease"), name = "Energy Impact") +
  labs(
    title = "Pareto Analysis: Guideline Impact on Energy",
    x = "Variant/Guideline"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

safe_save(p_pareto, "pareto_guideline_impact.png", width = 11, height = 7)

# 6) CPU vs Memory usage correlation
cpu_mem_data <- combined_runs %>%
  filter(!is.na(cpu_usage_avg), !is.na(memory_usage_avg), !is.na(energy_joules))

p_cpu_mem <- ggplot(cpu_mem_data, aes(x = cpu_usage_avg, y = memory_usage_avg, color = energy_joules, size = energy_joules)) +
  geom_point(alpha = 0.6) +
  scale_color_viridis_c(option = "plasma", trans = "log10", 
                        labels = label_number(scale_cut = cut_short_scale(), suffix = "J"),
                        name = "Energy") +
  scale_size_continuous(range = c(0.5, 4), guide = "none") +
  labs(
    title = "CPU vs Memory Usage with Energy Overlay",
    x = "CPU Usage (avg %)",
    y = "Memory Usage (avg %)"
  ) +
  theme_minimal(base_size = 12)

safe_save(p_cpu_mem, "cpu_memory_energy_correlation.png", width = 10, height = 7)

# 7) Benchmark difficulty: energy range across variants
benchmark_ranges <- summary_by_variant %>%
  group_by(benchmark) %>%
  summarise(
    min_energy = min(energy_mean, na.rm = TRUE),
    max_energy = max(energy_mean, na.rm = TRUE),
    range = max_energy - min_energy,
    baseline_energy = energy_mean[variant == "BASELINE"][1],
    pct_range = range / baseline_energy * 100,
    .groups = "drop"
  ) %>%
  arrange(desc(pct_range))

p_ranges <- ggplot(benchmark_ranges, aes(x = reorder(benchmark, pct_range), y = pct_range)) +
  geom_col(aes(fill = baseline_energy), width = 0.7) +
  scale_fill_viridis_c(option = "viridis", trans = "log10",
                       labels = label_number(scale_cut = cut_short_scale(), suffix = "J"),
                       name = "Baseline\nEnergy") +
  coord_flip() +
  labs(
    title = "Optimization Potential by Benchmark",
    subtitle = "Range of energy across all variants (as % of baseline)",
    x = "Benchmark",
    y = "Energy Range (% of baseline)"
  ) +
  theme_minimal(base_size = 11)

safe_save(p_ranges, "optimization_potential_by_benchmark.png", width = 10, height = 7)

# 8) Guideline combination matrix: which guidelines appear in top performers
top_performers <- summary_by_variant %>%
  filter(variant != "BASELINE") %>%
  group_by(benchmark) %>%
  arrange(pct_vs_baseline) %>%
  slice_head(n = 3) %>%
  ungroup() %>%
  select(benchmark, variant, pct_vs_baseline)

guideline_frequency <- top_performers %>%
  count(variant, name = "frequency") %>%
  arrange(desc(frequency))

p_frequency <- ggplot(guideline_frequency, aes(x = reorder(variant, frequency), y = frequency)) +
  geom_col(aes(fill = frequency), width = 0.7) +
  geom_text(aes(label = frequency), hjust = -0.3, size = 4) +
  scale_fill_gradient(low = "#fee5d9", high = "#a50f15", guide = "none") +
  coord_flip() +
  labs(
    title = "Most Frequent Guidelines in Top 3 Performers",
    subtitle = "Count of appearances in top 3 per benchmark",
    x = "Variant/Guideline",
    y = "Frequency"
  ) +
  theme_minimal(base_size = 12) +
  ylim(0, max(guideline_frequency$frequency) * 1.1)

safe_save(p_frequency, "top_performer_frequency.png", width = 10, height = 7)

# 9) Power efficiency: Energy per unit of work (if execution time available)
# Note: execution_time_seconds appears to be 0 in the data, using log-parsed times
power_efficiency <- combined_runs %>%
  filter(!is.na(exec_time_sec_log), exec_time_sec_log > 0, !is.na(energy_joules)) %>%
  mutate(
    power_eff = energy_joules / exec_time_sec_log,
    work_norm = 1 / exec_time_sec_log  # inverse time as proxy for throughput
  )

if (nrow(power_efficiency) > 0) {
  p_power_eff <- ggplot(power_efficiency, aes(x = variant, y = power_eff, fill = variant)) +
    geom_violin(alpha = 0.7, trim = TRUE) +
    geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.3) +
    scale_y_continuous(labels = label_number(scale_cut = cut_short_scale(), suffix = "W")) +
    scale_fill_brewer(palette = "Set3", guide = "none") +
    labs(
      title = "Power Efficiency Distribution by Variant",
      subtitle = "Energy / Execution Time (Watts)",
      x = "Variant",
      y = "Average Power (W)"
    ) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  safe_save(p_power_eff, "power_efficiency_by_variant.png", width = 11, height = 7)
}

# 10) Faceted comparison: energy vs baseline for each benchmark
facet_data <- summary_by_variant %>%
  filter(variant != "BASELINE") %>%
  mutate(variant = fct_reorder(variant, pct_vs_baseline))

p_facet <- ggplot(facet_data, aes(x = variant, y = pct_vs_baseline, fill = pct_vs_baseline < 0)) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.5) +
  facet_wrap(~ benchmark, scales = "free_x", ncol = 3) +
  scale_fill_manual(values = c("TRUE" = "#2ca02c", "FALSE" = "#d62728"), guide = "none") +
  labs(
    title = "% Change vs Baseline by Benchmark and Variant",
    x = NULL,
    y = "% Change vs Baseline"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    strip.text = element_text(face = "bold", size = 9)
  )

safe_save(p_facet, "faceted_pct_change_all_benchmarks.png", width = 14, height = 10)

# 11) OPTIMIZED vs BASELINE comparison across benchmarks
optimized_comparison <- summary_by_variant %>%
  filter(variant %in% c("BASELINE", "OPTIMIZED")) %>%
  select(benchmark, variant, energy_mean) %>%
  pivot_wider(names_from = variant, values_from = energy_mean) %>%
  mutate(
    pct_change = (OPTIMIZED - BASELINE) / BASELINE * 100,
    improvement = pct_change < 0
  ) %>%
  arrange(pct_change)

p_optimized <- ggplot(optimized_comparison, aes(x = reorder(benchmark, pct_change), y = pct_change)) +
  geom_col(aes(fill = improvement), width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30", linewidth = 0.8) +
  scale_fill_manual(values = c("TRUE" = "#2ca02c", "FALSE" = "#d62728"),
                    labels = c("Improvement", "Regression"), name = NULL) +
  coord_flip() +
  labs(
    title = "OPTIMIZED vs BASELINE: Overall Impact",
    subtitle = "% change when all guidelines are applied",
    x = "Benchmark",
    y = "% Change vs Baseline"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

safe_save(p_optimized, "optimized_vs_baseline_comparison.png", width = 10, height = 8)

# 12) Statistical significance indicators (simple t-test visualization)
baseline_data <- combined_runs %>%
  filter(variant == "BASELINE") %>%
  select(benchmark, baseline_energy = energy_joules)

significance_tests <- combined_runs %>%
  filter(variant != "BASELINE") %>%
  group_by(benchmark, variant) %>%
  summarise(
    mean_energy = mean(energy_joules, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    baseline_data %>% group_by(benchmark) %>% 
      summarise(baseline_mean = mean(baseline_energy, na.rm = TRUE), .groups = "drop"),
    by = "benchmark"
  ) %>%
  mutate(
    abs_diff = abs(mean_energy - baseline_mean),
    pct_diff = abs((mean_energy - baseline_mean) / baseline_mean * 100)
  ) %>%
  arrange(desc(abs_diff)) %>%
  head(20)

p_top_changes <- ggplot(significance_tests, 
                        aes(x = reorder(paste(benchmark, variant, sep = " - "), abs_diff), 
                            y = pct_diff)) +
  geom_col(aes(fill = mean_energy < baseline_mean), width = 0.7) +
  scale_fill_manual(values = c("TRUE" = "#2ca02c", "FALSE" = "#d62728"),
                    labels = c("Improvement", "Regression"), name = "Impact") +
  coord_flip() +
  labs(
    title = "Top 20 Largest Energy Changes",
    subtitle = "Benchmark-Variant combinations with biggest absolute differences",
    x = "Benchmark - Variant",
    y = "% Change vs Baseline"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top")

safe_save(p_top_changes, "top_20_energy_changes.png", width = 11, height = 9)

message("Done! All additional plots generated in: ", PLOTS_DIR)

