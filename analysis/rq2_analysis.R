#!/usr/bin/env Rscript
# =============================================================================
# RQ2 Analysis: Multithreading Guidelines vs Baseline
# =============================================================================
# Research Question 2: Do multithreading optimization guidelines significantly
# impact energy consumption and related metrics compared to baseline?
# Multithreading guidelines: G9, G12, G14

# Load required packages
packages <- c("ggplot2", "readr", "dplyr", "tidyr", "scales", "forcats", "ggrepel")
for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}
suppressPackageStartupMessages({
  invisible(lapply(packages, function(p) library(p, character.only = TRUE)))
})

# =============================================================================
# CONFIGURATION
# =============================================================================
WORKSPACE_ROOT <- "/Users/kellywang/Documents/compSci/p1/greenLab/GreenLab"
OUT_DIR <- file.path(WORKSPACE_ROOT, "analysis", "outputs")
PLOTS_DIR <- file.path(WORKSPACE_ROOT, "analysis", "plots", "rq2")
dir.create(PLOTS_DIR, recursive = TRUE, showWarnings = FALSE)

message("=============================================================================")
message("RQ2 Analysis: Multithreading Guidelines vs Baseline")
message("=============================================================================")
message("Output directory: ", PLOTS_DIR)
message("")

# =============================================================================
# LOAD DATA
# =============================================================================
data <- read_csv(file.path(OUT_DIR, "cleaned_experiment_data.csv"), show_col_types = FALSE)

# Define multithreading guidelines (G8-G14, but only G9, G12, G14 exist in data)
multithreading_guidelines <- c("G9", "G12", "G14")

message("Data loaded: ", nrow(data), " measurements")
message("Multithreading guidelines: ", paste(multithreading_guidelines, collapse = ", "))

# Filter to only BASELINE and multithreading guidelines
rq2_data <- data %>%
  filter(variant %in% c("BASELINE", multithreading_guidelines))

message("RQ2 data: ", nrow(rq2_data), " measurements")
message("Benchmarks: ", n_distinct(rq2_data$benchmark))

# =============================================================================
# CALCULATE Δ% CHANGES PER BENCHMARK AND GUIDELINE
# =============================================================================
message("\nCalculating Δ% changes per benchmark and guideline...")

# Get baseline values per benchmark
baseline_values <- rq2_data %>%
  filter(variant == "BASELINE") %>%
  group_by(benchmark) %>%
  summarise(
    baseline_energy = mean(energy_joules, na.rm = TRUE),
    baseline_time = mean(exec_time_sec, na.rm = TRUE),
    baseline_cpu = mean(cpu_usage_avg, na.rm = TRUE),
    baseline_memory = mean(memory_usage_avg, na.rm = TRUE),
    baseline_power = mean(power_watts, na.rm = TRUE),
    .groups = "drop"
  )

# Calculate Δ% for each guideline-benchmark combination
guideline_deltas <- rq2_data %>%
  filter(variant %in% multithreading_guidelines) %>%
  left_join(baseline_values, by = "benchmark") %>%
  mutate(
    delta_energy_pct = ((energy_joules - baseline_energy) / baseline_energy) * 100,
    delta_time_pct = ((exec_time_sec - baseline_time) / baseline_time) * 100,
    delta_cpu_pct = ((cpu_usage_avg - baseline_cpu) / baseline_cpu) * 100,
    delta_memory_pct = ((memory_usage_avg - baseline_memory) / baseline_memory) * 100,
    delta_power_pct = ((power_watts - baseline_power) / baseline_power) * 100
  ) %>%
  group_by(benchmark, variant) %>%
  summarise(
    # Energy
    energy_mean = mean(delta_energy_pct, na.rm = TRUE),
    energy_se = sd(delta_energy_pct, na.rm = TRUE) / sqrt(n()),
    energy_ci_lower = energy_mean - 1.96 * energy_se,
    energy_ci_upper = energy_mean + 1.96 * energy_se,
    # Time
    time_mean = mean(delta_time_pct, na.rm = TRUE),
    time_se = sd(delta_time_pct, na.rm = TRUE) / sqrt(n()),
    time_ci_lower = time_mean - 1.96 * time_se,
    time_ci_upper = time_mean + 1.96 * time_se,
    # CPU
    cpu_mean = mean(delta_cpu_pct, na.rm = TRUE),
    cpu_se = sd(delta_cpu_pct, na.rm = TRUE) / sqrt(n()),
    cpu_ci_lower = cpu_mean - 1.96 * cpu_se,
    cpu_ci_upper = cpu_mean + 1.96 * cpu_se,
    # Memory
    memory_mean = mean(delta_memory_pct, na.rm = TRUE),
    memory_se = sd(delta_memory_pct, na.rm = TRUE) / sqrt(n()),
    memory_ci_lower = memory_mean - 1.96 * memory_se,
    memory_ci_upper = memory_mean + 1.96 * memory_se,
    # Power
    power_mean = mean(delta_power_pct, na.rm = TRUE),
    power_se = sd(delta_power_pct, na.rm = TRUE) / sqrt(n()),
    power_ci_lower = power_mean - 1.96 * power_se,
    power_ci_upper = power_mean + 1.96 * power_se,
    .groups = "drop"
  )

# Save results
write_csv(guideline_deltas, file.path(OUT_DIR, "rq2_multithreading_deltas.csv"))
message("Saved: rq2_multithreading_deltas.csv")

# =============================================================================
# PLOT 1: Grouped Bar Chart - All Metrics per Benchmark (Averaged Multithreading)
# =============================================================================
message("\nCreating Plot 1: All Metrics per Benchmark (averaged multithreading)...")

# Calculate average multithreading effect (excluding basic_string + G12 outlier)
plot1_data_raw <- rq2_data %>%
  filter(variant %in% multithreading_guidelines) %>%
  filter(!(benchmark == "basic_string" & variant == "G12")) %>%
  select(benchmark, variant, results_set, energy_joules, cpu_usage_avg, memory_usage_avg, power_watts)

# Calculate baseline values for each benchmark
baseline_values_plot1 <- rq2_data %>%
  filter(variant == "BASELINE") %>%
  group_by(benchmark) %>%
  summarise(
    baseline_energy = mean(energy_joules, na.rm = TRUE),
    baseline_cpu = mean(cpu_usage_avg, na.rm = TRUE),
    baseline_memory = mean(memory_usage_avg, na.rm = TRUE),
    baseline_power = mean(power_watts, na.rm = TRUE),
    .groups = "drop"
  )

# Calculate Δ% and statistics for averaged multithreading
plot1_stats <- plot1_data_raw %>%
  left_join(baseline_values_plot1, by = "benchmark") %>%
  mutate(
    delta_energy_pct = ((energy_joules - baseline_energy) / baseline_energy) * 100,
    delta_cpu_pct = ((cpu_usage_avg - baseline_cpu) / baseline_cpu) * 100,
    delta_memory_pct = ((memory_usage_avg - baseline_memory) / baseline_memory) * 100,
    delta_power_pct = ((power_watts - baseline_power) / baseline_power) * 100
  ) %>%
  group_by(benchmark) %>%
  summarise(
    # Energy
    energy_mean = mean(delta_energy_pct, na.rm = TRUE),
    energy_se = sd(delta_energy_pct, na.rm = TRUE) / sqrt(n()),
    energy_ci_lower = energy_mean - 1.96 * energy_se,
    energy_ci_upper = energy_mean + 1.96 * energy_se,
    # CPU
    cpu_mean = mean(delta_cpu_pct, na.rm = TRUE),
    cpu_se = sd(delta_cpu_pct, na.rm = TRUE) / sqrt(n()),
    cpu_ci_lower = cpu_mean - 1.96 * cpu_se,
    cpu_ci_upper = cpu_mean + 1.96 * cpu_se,
    # Memory
    memory_mean = mean(delta_memory_pct, na.rm = TRUE),
    memory_se = sd(delta_memory_pct, na.rm = TRUE) / sqrt(n()),
    memory_ci_lower = memory_mean - 1.96 * memory_se,
    memory_ci_upper = memory_mean + 1.96 * memory_se,
    # Power
    power_mean = mean(delta_power_pct, na.rm = TRUE),
    power_se = sd(delta_power_pct, na.rm = TRUE) / sqrt(n()),
    power_ci_lower = power_mean - 1.96 * power_se,
    power_ci_upper = power_mean + 1.96 * power_se,
    .groups = "drop"
  )

# Reshape for plotting
plot1_long <- plot1_stats %>%
  pivot_longer(
    cols = c(energy_mean, cpu_mean, memory_mean, power_mean),
    names_to = "metric",
    values_to = "delta_pct"
  ) %>%
  mutate(
    metric_clean = case_when(
      metric == "energy_mean" ~ "Energy %",
      metric == "cpu_mean" ~ "CPU %",
      metric == "memory_mean" ~ "Memory %",
      metric == "power_mean" ~ "Power %",
      TRUE ~ metric
    ),
    ci_lower = case_when(
      metric == "energy_mean" ~ energy_ci_lower,
      metric == "cpu_mean" ~ cpu_ci_lower,
      metric == "memory_mean" ~ memory_ci_lower,
      metric == "power_mean" ~ power_ci_lower,
      TRUE ~ NA_real_
    ),
    ci_upper = case_when(
      metric == "energy_mean" ~ energy_ci_upper,
      metric == "cpu_mean" ~ cpu_ci_upper,
      metric == "memory_mean" ~ memory_ci_upper,
      metric == "power_mean" ~ power_ci_upper,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(!is.na(delta_pct), is.finite(delta_pct)) %>%
  mutate(metric_clean = factor(metric_clean, levels = c("Energy %", "CPU %", "Memory %", "Power %")))

# Sort benchmarks by energy change
benchmark_order <- plot1_long %>%
  filter(metric_clean == "Energy %") %>%
  arrange(delta_pct) %>%
  pull(benchmark)

plot1_long <- plot1_long %>%
  mutate(benchmark = factor(benchmark, levels = benchmark_order))

# Calculate median energy change for annotation
median_energy <- plot1_long %>%
  filter(metric_clean == "Energy %") %>%
  summarise(median = median(delta_pct, na.rm = TRUE)) %>%
  pull(median)

# Check for extreme values and add truncation
y_limit <- 100  # Cap at ±100%

# Identify truncated values
plot1_data_labels <- plot1_long %>%
  mutate(
    is_truncated = abs(delta_pct) > y_limit,
    label_text = ifelse(is_truncated, sprintf("%+.0f%%", round(delta_pct)), ""),
    label_y = ifelse(delta_pct > y_limit, y_limit - 8,
                    ifelse(delta_pct < -y_limit, -y_limit + 8, NA)),
    delta_pct_display = case_when(
      delta_pct > y_limit ~ y_limit,
      delta_pct < -y_limit ~ -y_limit,
      TRUE ~ delta_pct
    )
  )

# Create diverging color scale (red = worse, blue = better)
p1 <- ggplot(plot1_data_labels, aes(x = benchmark, y = delta_pct_display, fill = delta_pct)) +
  geom_col(position = position_dodge(width = 0.8), alpha = 0.85, width = 0.75) +
  # Add text labels for truncated values
  geom_text(data = plot1_data_labels %>% filter(is_truncated),
            aes(x = benchmark, y = label_y, label = label_text),
            size = 2.5, fontface = "bold", color = "black", angle = 0) +
  geom_errorbar(aes(ymin = pmax(ci_lower, -y_limit), ymax = pmin(ci_upper, y_limit)),
                position = position_dodge(width = 0.8),
                width = 0.25, linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  facet_wrap(~metric_clean, ncol = 2, scales = "fixed") +
  coord_cartesian(ylim = c(-y_limit, y_limit)) +
  scale_fill_gradient2(
    low = "#2166ac", mid = "#f7f7f7", high = "#b2182b",
    midpoint = 0,
    limits = c(-100, 100),
    name = "Δ%",
    labels = function(x) sprintf("%+.0f%%", x),
    na.value = "gray50"
  ) +
  scale_y_continuous(labels = function(x) sprintf("%+.0f%%", x)) +
  labs(
    title = "Per-Benchmark Impact of Multithreading vs Baseline",
    subtitle = sprintf("Median ΔEnergy: %+.1f%% (positive = higher energy consumption)", median_energy),
    x = "Benchmark",
    y = "Δ% = (Multithreading − Baseline) / Baseline × 100",
    caption = sprintf("Bars show mean percentage change; error bars show 95%% CI. Y-axis capped at ±%d%% for readability.\nDashed line = no change. Red = increase/worse, Blue = decrease/better. Actual values shown for truncated bars.\nMultithreading = average of G9, G12, G14 (basic_string+G12 outlier excluded).", y_limit)
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10, color = "gray30"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "right",
    panel.grid.major.x = element_blank(),
    plot.caption = element_text(hjust = 0, size = 8, color = "gray40")
  )

# Save plot
ggsave(file.path(PLOTS_DIR, "plot1_multithreading_all_metrics_per_benchmark.png"), 
       plot = p1, width = 14, height = 10, dpi = 300)
message("Saved: plot1_multithreading_all_metrics_per_benchmark.png")

# =============================================================================
# PLOT 2: Paired Slope Chart - Baseline vs Multithreading (Averaged)
# =============================================================================
message("\nCreating Plot 2: Paired Per-Benchmark Energy Comparison...")

# Calculate average multithreading effect for each benchmark
# (Average across G9, G12, G14, excluding the basic_string + G12 outlier)
multithreading_avg <- rq2_data %>%
  filter(variant %in% multithreading_guidelines) %>%
  filter(!(benchmark == "basic_string" & variant == "G12")) %>%  # Exclude outlier
  group_by(benchmark, results_set) %>%
  summarise(
    energy_joules = mean(energy_joules, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(variant = "MULTITHREADING")

# Combine baseline and averaged multithreading
slope_data_raw <- bind_rows(
  rq2_data %>% filter(variant == "BASELINE") %>% select(benchmark, variant, results_set, energy_joules),
  multithreading_avg
)

# Calculate means for slope chart
slope_data <- slope_data_raw %>%
  group_by(benchmark, variant) %>%
  summarise(
    energy_mean = mean(energy_joules, na.rm = TRUE),
    energy_se = sd(energy_joules, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(variant = factor(variant, levels = c("BASELINE", "MULTITHREADING")))

# Determine if multithreading improved or worsened energy for each benchmark
benchmark_improvement <- slope_data %>%
  pivot_wider(names_from = variant, values_from = c(energy_mean, energy_se)) %>%
  mutate(
    improved = energy_mean_MULTITHREADING < energy_mean_BASELINE,
    improvement_label = ifelse(improved, "Improved (↓)", "Worsened (↑)")
  ) %>%
  select(benchmark, improved, improvement_label)

# Add improvement info back to slope_data
slope_data <- slope_data %>%
  left_join(benchmark_improvement, by = "benchmark") %>%
  mutate(improvement_label = factor(improvement_label, levels = c("Improved (↓)", "Worsened (↑)")))

# Calculate overall means
overall_means <- slope_data %>%
  group_by(variant) %>%
  summarise(
    overall_mean = mean(energy_mean, na.rm = TRUE),
    overall_se = sd(energy_mean, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# Calculate median delta
baseline_avg <- slope_data %>% filter(variant == "BASELINE") %>% pull(energy_mean)
multithreading_avg_values <- slope_data %>% filter(variant == "MULTITHREADING") %>% pull(energy_mean)
delta_values <- ((multithreading_avg_values - baseline_avg) / baseline_avg) * 100
median_delta_energy <- median(delta_values, na.rm = TRUE)

# Define colors for benchmarks (using a colorblind-friendly palette)
benchmark_colors <- setNames(
  scales::hue_pal()(length(unique(slope_data$benchmark))),
  unique(slope_data$benchmark)
)

# Create the slope plot with linetype based on improvement
p2 <- ggplot(slope_data, aes(x = variant, y = energy_mean, group = benchmark, 
                              color = benchmark, linetype = improvement_label)) +
  # Individual benchmark lines (with color and linetype)
  geom_line(linewidth = 1.2, alpha = 0.8) +
  geom_point(size = 3, alpha = 0.9) +
  # Overall mean lines (bold black)
  geom_line(data = overall_means, 
            aes(x = variant, y = overall_mean, group = 1),
            color = "black", linewidth = 2, inherit.aes = FALSE) +
  geom_point(data = overall_means,
             aes(x = variant, y = overall_mean),
             color = "black", size = 5, shape = 18, inherit.aes = FALSE) +
  # Error bars for overall means
  geom_errorbar(data = overall_means,
                aes(x = variant, y = overall_mean, 
                    ymin = overall_mean - overall_se, 
                    ymax = overall_mean + overall_se),
                color = "black", width = 0.1, linewidth = 1.2,
                inherit.aes = FALSE) +
  # Annotation for median delta
  annotate("text", x = 1.5, y = max(slope_data$energy_mean) * 0.95,
           label = sprintf("Median ΔEnergy: %+.1f%%\nOverall mean shown in black", median_delta_energy),
           size = 4, fontface = "bold", color = "black") +
  # Styling
  scale_y_log10(labels = label_number(scale_cut = cut_short_scale(), suffix = " J")) +
  scale_color_manual(values = benchmark_colors, name = "Benchmark") +
  scale_linetype_manual(
    values = c("Improved (↓)" = "dashed", "Worsened (↑)" = "solid"),
    name = "Effect",
    guide = guide_legend(order = 1)
  ) +
  labs(
    title = "Per-Benchmark Energy Consumption: Baseline vs Multithreading",
    subtitle = "Each colored line connects paired runs from the same benchmark",
    x = NULL,
    y = "Energy Consumption (J, log scale)",
    caption = "Dashed lines = energy improved (↓); Solid lines = energy worsened (↑).\nBlack line = overall mean ± SE. Upward slope indicates energy increase.\nMultithreading = average of G9, G12, G14 (basic_string+G12 outlier excluded)."
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "gray30"),
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    plot.caption = element_text(hjust = 0, size = 8, color = "gray40"),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 9),
    legend.box = "vertical"
  ) +
  guides(
    linetype = guide_legend(order = 1, override.aes = list(linewidth = 1.5, alpha = 1, color = "black")),
    color = guide_legend(order = 2, override.aes = list(linewidth = 1.2, alpha = 1))
  )

# Save plot
ggsave(file.path(PLOTS_DIR, "plot2_paired_energy_comparison.png"), 
       plot = p2, width = 10, height = 8, dpi = 300)
message("Saved: plot2_paired_energy_comparison.png")

# =============================================================================
# SUMMARY STATISTICS
# =============================================================================
message("\nCalculating summary statistics...")

summary_stats <- guideline_deltas %>%
  group_by(variant) %>%
  summarise(
    n_benchmarks = n(),
    energy_median = median(energy_mean, na.rm = TRUE),
    energy_mean_overall = mean(energy_mean, na.rm = TRUE),
    time_median = median(time_mean, na.rm = TRUE),
    cpu_median = median(cpu_mean, na.rm = TRUE),
    memory_median = median(memory_mean, na.rm = TRUE),
    power_median = median(power_mean, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(summary_stats, file.path(OUT_DIR, "rq2_summary_statistics.csv"))
message("Saved: rq2_summary_statistics.csv")

# =============================================================================
# FINAL REPORT
# =============================================================================
message("\n=============================================================================")
message("RQ2 ANALYSIS SUMMARY")
message("=============================================================================")

message("\nResearch Question 2:")
message("Do multithreading optimization guidelines significantly impact energy")
message("consumption and related metrics compared to baseline?")

message("\n--- SUMMARY STATISTICS ---")
for (i in 1:nrow(summary_stats)) {
  row <- summary_stats[i, ]
  message(sprintf("\n%s:", row$variant))
  message(sprintf("  Benchmarks analyzed: %d", row$n_benchmarks))
  message(sprintf("  Median ΔEnergy: %+.2f%%", row$energy_median))
  message(sprintf("  Median ΔTime: %+.2f%%", row$time_median))
  message(sprintf("  Median ΔCPU: %+.2f%%", row$cpu_median))
  message(sprintf("  Median ΔMemory: %+.2f%%", row$memory_median))
  message(sprintf("  Median ΔPower: %+.2f%%", row$power_median))
}

message("\n--- OUTPUT FILES ---")
message("Plots:")
message("  - ", file.path(PLOTS_DIR, "plot1_multithreading_all_metrics_per_benchmark.png"))
message("  - ", file.path(PLOTS_DIR, "plot2_paired_energy_comparison.png"))
message("\nData tables:")
message("  - ", file.path(OUT_DIR, "rq2_multithreading_deltas.csv"))
message("  - ", file.path(OUT_DIR, "rq2_summary_statistics.csv"))

message("\n=============================================================================")
message("RQ2 Analysis Complete!")
message("=============================================================================\n")

