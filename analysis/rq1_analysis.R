#!/usr/bin/env Rscript
# =============================================================================
# RQ1 Analysis: Code Optimization vs Baseline
# =============================================================================
# Research Question 1: Does code optimization significantly impact energy 
# consumption and related metrics compared to baseline?

# Load required packages
packages <- c("ggplot2", "readr", "dplyr", "tidyr", "scales", "forcats")
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
PLOTS_DIR <- file.path(WORKSPACE_ROOT, "analysis", "plots", "rq1")
dir.create(PLOTS_DIR, recursive = TRUE, showWarnings = FALSE)

message("=============================================================================")
message("RQ1 Analysis: Code Optimization vs Baseline")
message("=============================================================================")
message("Output directory: ", PLOTS_DIR)
message("")

# =============================================================================
# LOAD DATA
# =============================================================================
data <- read_csv(file.path(OUT_DIR, "cleaned_experiment_data.csv"), show_col_types = FALSE)

message("Data loaded: ", nrow(data), " measurements")

# Filter to only BASELINE and OPTIMIZED variants
rq1_data <- data %>%
  filter(variant %in% c("BASELINE", "OPTIMIZED"))

message("RQ1 data: ", nrow(rq1_data), " measurements (BASELINE + OPTIMIZED)")
message("Benchmarks: ", n_distinct(rq1_data$benchmark))

# =============================================================================
# CALCULATE Δ% CHANGES PER BENCHMARK
# =============================================================================
message("\nCalculating Δ% changes per benchmark...")

# Get baseline values per benchmark (averaged across runs)
baseline_values <- rq1_data %>%
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

# Get optimized values per benchmark
optimized_values <- rq1_data %>%
  filter(variant == "OPTIMIZED") %>%
  group_by(benchmark) %>%
  summarise(
    opt_energy = mean(energy_joules, na.rm = TRUE),
    opt_time = mean(exec_time_sec, na.rm = TRUE),
    opt_cpu = mean(cpu_usage_avg, na.rm = TRUE),
    opt_memory = mean(memory_usage_avg, na.rm = TRUE),
    opt_power = mean(power_watts, na.rm = TRUE),
    .groups = "drop"
  )

# Calculate Δ% for each metric
delta_pct <- baseline_values %>%
  left_join(optimized_values, by = "benchmark") %>%
  mutate(
    delta_energy_pct = ((opt_energy - baseline_energy) / baseline_energy) * 100,
    delta_time_pct = ((opt_time - baseline_time) / baseline_time) * 100,
    delta_cpu_pct = ((opt_cpu - baseline_cpu) / baseline_cpu) * 100,
    delta_memory_pct = ((opt_memory - baseline_memory) / baseline_memory) * 100,
    delta_power_pct = ((opt_power - baseline_power) / baseline_power) * 100
  ) %>%
  select(benchmark, starts_with("delta_"))

message("Calculated Δ% changes for ", nrow(delta_pct), " benchmarks")

# Save delta calculations
write_csv(delta_pct, file.path(OUT_DIR, "rq1_delta_pct_by_benchmark.csv"))
message("Saved: rq1_delta_pct_by_benchmark.csv")

# =============================================================================
# STATISTICAL TESTS
# =============================================================================
message("\nPerforming statistical tests...")

# Reshape data for testing
delta_long <- delta_pct %>%
  pivot_longer(
    cols = starts_with("delta_"),
    names_to = "metric",
    values_to = "delta_pct"
  ) %>%
  mutate(
    metric = case_when(
      metric == "delta_energy_pct" ~ "Energy %",
      metric == "delta_time_pct" ~ "Time %",
      metric == "delta_cpu_pct" ~ "CPU %",
      metric == "delta_memory_pct" ~ "Memory %",
      metric == "delta_power_pct" ~ "Power %",
      TRUE ~ metric
    )
  )

# Perform tests for each metric
test_results <- delta_long %>%
  group_by(metric) %>%
  summarise(
    n = n(),
    median = median(delta_pct, na.rm = TRUE),
    mean = mean(delta_pct, na.rm = TRUE),
    sd = sd(delta_pct, na.rm = TRUE),
    # Shapiro-Wilk test for normality
    shapiro_p = tryCatch(
      shapiro.test(delta_pct)$p.value,
      error = function(e) NA
    ),
    # Wilcoxon signed-rank test (H0: median = 0)
    wilcox_p = tryCatch(
      wilcox.test(delta_pct, mu = 0, alternative = "two.sided")$p.value,
      error = function(e) NA
    ),
    # Paired t-test (H0: mean = 0)
    ttest_p = tryCatch(
      t.test(delta_pct, mu = 0, alternative = "two.sided")$p.value,
      error = function(e) NA
    ),
    .groups = "drop"
  ) %>%
  mutate(
    # Use Wilcoxon if non-normal, t-test if normal
    test_used = ifelse(shapiro_p < 0.05, "Wilcoxon", "t-test"),
    p_value = ifelse(shapiro_p < 0.05, wilcox_p, ttest_p),
    significant = p_value < 0.05,
    significance_label = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

print(test_results)

# Save test results
write_csv(test_results, file.path(OUT_DIR, "rq1_statistical_tests.csv"))
message("\nSaved: rq1_statistical_tests.csv")

# =============================================================================
# PLOT 1a: Grouped Bar Chart - All Metrics per Benchmark (Averaged Code Optimization)
# =============================================================================
message("\nCreating Plot 1a: All Metrics per Benchmark (averaged code optimization)...")

# Get all code optimization guidelines (excluding G12 for basic_string due to outlier)
code_opt_guidelines <- data %>%
  filter(!variant %in% c("BASELINE", "OPTIMIZED", "CLEAN", "ORIGINAL_FIXED")) %>%
  pull(variant) %>%
  unique()

message("Code optimization guidelines: ", paste(code_opt_guidelines, collapse = ", "))

# Calculate average code optimization effect (excluding basic_string + G12 outlier)
plot1a_data_raw <- data %>%
  filter(variant %in% code_opt_guidelines) %>%
  filter(!(benchmark == "basic_string" & variant == "G12")) %>%
  select(benchmark, variant, results_set, energy_joules, cpu_usage_avg, memory_usage_avg, power_watts)

# Calculate baseline values for each benchmark
baseline_values_plot1a <- data %>%
  filter(variant == "BASELINE") %>%
  group_by(benchmark) %>%
  summarise(
    baseline_energy = mean(energy_joules, na.rm = TRUE),
    baseline_cpu = mean(cpu_usage_avg, na.rm = TRUE),
    baseline_memory = mean(memory_usage_avg, na.rm = TRUE),
    baseline_power = mean(power_watts, na.rm = TRUE),
    .groups = "drop"
  )

# Calculate Δ% and statistics for averaged code optimization
plot1a_stats <- plot1a_data_raw %>%
  left_join(baseline_values_plot1a, by = "benchmark") %>%
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
plot1a_long <- plot1a_stats %>%
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
benchmark_order <- plot1a_long %>%
  filter(metric_clean == "Energy %") %>%
  arrange(delta_pct) %>%
  pull(benchmark)

plot1a_long <- plot1a_long %>%
  mutate(benchmark = factor(benchmark, levels = benchmark_order))

# Calculate median energy change for annotation
median_energy <- plot1a_long %>%
  filter(metric_clean == "Energy %") %>%
  summarise(median = median(delta_pct, na.rm = TRUE)) %>%
  pull(median)

# Check for extreme values and add truncation
y_limit <- 100  # Cap at ±100%

# Identify truncated values
plot1a_data_labels <- plot1a_long %>%
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
p1a <- ggplot(plot1a_data_labels, aes(x = benchmark, y = delta_pct_display, fill = delta_pct)) +
  geom_col(position = position_dodge(width = 0.8), alpha = 0.85, width = 0.75) +
  # Add text labels for truncated values
  geom_text(data = plot1a_data_labels %>% filter(is_truncated),
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
    title = "Per-Benchmark Impact of Code Optimization vs Baseline",
    subtitle = sprintf("Median ΔEnergy: %+.1f%% (negative = lower energy consumption)", median_energy),
    x = "Benchmark",
    y = "Δ% = (Code Optimization − Baseline) / Baseline × 100",
    caption = sprintf("Bars show mean percentage change; error bars show 95%% CI. Y-axis capped at ±%d%% for readability.\nDashed line = no change. Red = increase/worse, Blue = decrease/better. Actual values shown for truncated bars.\nCode Optimization = average of G1, G3, G4, G6, G7, G9, G14 (basic_string+G12 outlier excluded).", y_limit)
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
ggsave(file.path(PLOTS_DIR, "plot1a_code_optimization_all_metrics_per_benchmark.png"), 
       plot = p1a, width = 14, height = 10, dpi = 300)
message("Saved: plot1a_code_optimization_all_metrics_per_benchmark.png")

# =============================================================================
# PLOT 1b: Per-Guideline Impact on Each Metric for Representative Benchmarks
# =============================================================================
message("\nCreating Plot 1b: Per-Guideline Impact on Metrics...")

# Select representative benchmarks
target_benchmarks <- c("sequential_minimum_optimization", "basic_string")
message("Using benchmarks: ", paste(target_benchmarks, collapse = ", "))

# Get data for all variants (including guidelines) for these benchmarks
plot1_data_raw <- data %>%
  filter(benchmark %in% target_benchmarks,
         variant != "OPTIMIZED") %>%  # Exclude OPTIMIZED, keep BASELINE and guidelines
  select(benchmark, variant, results_set, energy_joules, exec_time_sec, 
         cpu_usage_avg, memory_usage_avg, power_watts)

# Calculate baseline values for each benchmark
baseline_values_plot1 <- plot1_data_raw %>%
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

# Calculate Δ% and statistics for each guideline
plot1_stats <- plot1_data_raw %>%
  filter(variant != "BASELINE") %>%  # Only guidelines
  left_join(baseline_values_plot1, by = "benchmark") %>%
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
    # Time (will be NA)
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
  mutate(
    metric_clean = factor(metric_clean, levels = c("Energy %", "CPU %", "Memory %", "Power %")),
    benchmark_clean = case_when(
      benchmark == "sequential_minimum_optimization" ~ "Sequential Min. Opt.",
      benchmark == "basic_string" ~ "Basic String",
      TRUE ~ benchmark
    )
  )

# Create separate plots for each benchmark
for (bench in target_benchmarks) {
  bench_clean <- ifelse(bench == "sequential_minimum_optimization", 
                        "Sequential Min. Opt.", "Basic String")
  
  plot_data <- plot1_long %>% filter(benchmark == bench)
  
  # Check for extreme outliers (>1000% or <-1000%)
  has_extreme <- any(abs(plot_data$delta_pct) > 1000, na.rm = TRUE)
  extreme_variants <- plot_data %>% 
    filter(abs(delta_pct) > 1000) %>% 
    pull(variant) %>% 
    unique()
  
  # For basic_string, filter out extreme outliers for readability
  if (bench == "basic_string" && has_extreme) {
    message("  Note: Removing extreme outliers for ", bench, ": ", 
            paste(extreme_variants, collapse = ", "))
    plot_data_filtered <- plot_data %>% filter(abs(delta_pct) <= 1000)
    
    # Create note about excluded data
    excluded_text <- sprintf("Note: G12 excluded (Energy: +%s%%)", 
                            format(round(max(plot_data$delta_pct[plot_data$variant == "G12"], 
                                           na.rm = TRUE)), 
                                  big.mark = ","))
    subtitle_text <- sprintf("Percentage change relative to baseline (extreme outliers >±1000%% excluded)\n%s", 
                            excluded_text)
  } else {
    plot_data_filtered <- plot_data
    subtitle_text <- "Percentage change relative to baseline for each optimization guideline"
  }
  
  # Set reasonable y-axis limits for readability
  # Different limits for different benchmarks
  y_limit <- ifelse(bench == "basic_string", 100, 50)
  
  # Identify truncated values (beyond ±50%)
  plot_data_with_labels <- plot_data_filtered %>%
    mutate(
      is_truncated = abs(delta_pct) > y_limit,
      display_value = ifelse(delta_pct > y_limit, y_limit - 5,
                            ifelse(delta_pct < -y_limit, -y_limit + 5, delta_pct)),
      label_text = ifelse(is_truncated, sprintf("%+.0f%%", round(delta_pct)), ""),
      label_y = ifelse(delta_pct > y_limit, y_limit - 2,
                      ifelse(delta_pct < -y_limit, -y_limit + 2, NA))
    )
  
  # Truncate values for display
  plot_data_display <- plot_data_with_labels %>%
    mutate(
      delta_pct_display = case_when(
        delta_pct > y_limit ~ y_limit,
        delta_pct < -y_limit ~ -y_limit,
        TRUE ~ delta_pct
      )
    )
  
  p <- ggplot(plot_data_display, aes(x = metric_clean, y = delta_pct_display, fill = variant)) +
    geom_col(position = position_dodge(width = 0.8), alpha = 0.8, width = 0.75) +
    # Add text labels for truncated values
    geom_text(data = plot_data_with_labels %>% filter(is_truncated),
              aes(x = metric_clean, y = label_y, label = label_text),
              position = position_dodge(width = 0.8),
              size = 3, fontface = "bold", color = "black") +
    geom_errorbar(aes(ymin = pmax(ci_lower, -y_limit), ymax = pmin(ci_upper, y_limit)),
                  position = position_dodge(width = 0.8),
                  width = 0.25, linewidth = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
    scale_fill_brewer(palette = "Set3", name = "Guideline") +
    scale_y_continuous(
      labels = function(x) sprintf("%+.0f%%", x),
      limits = c(-y_limit, y_limit),
      breaks = seq(-y_limit, y_limit, by = 10)
    ) +
    labs(
      title = sprintf("Per-Guideline Impact on Metrics: %s", bench_clean),
      subtitle = subtitle_text,
      x = "Metric",
      y = "Δ% = (Guideline − Baseline) / Baseline × 100",
      caption = sprintf("Bars show mean percentage change; error bars show 95%% CI over repeated runs.\nDashed red line = no change. Y-axis capped at ±%d%% for readability; actual values shown for truncated bars.", y_limit)
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "gray30"),
      axis.text.x = element_text(size = 11, face = "bold"),
      axis.text.y = element_text(size = 10),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      panel.grid.major.x = element_blank(),
      plot.caption = element_text(hjust = 0, size = 8, color = "gray40")
    )
  
  # Save plot
  filename <- sprintf("plot1b_%s_guideline_impact.png", bench)
  ggsave(file.path(PLOTS_DIR, filename), 
         plot = p, width = 14, height = 8, dpi = 300)
  message("Saved: ", filename)
}

# =============================================================================
# PLOT 2: Paired Slope Chart for Energy
# =============================================================================
message("\nCreating Plot 2: Paired Per-Benchmark Energy Comparison...")

# Prepare data for slope chart
slope_data <- rq1_data %>%
  filter(variant %in% c("BASELINE", "OPTIMIZED")) %>%
  group_by(benchmark, variant) %>%
  summarise(
    energy_mean = mean(energy_joules, na.rm = TRUE),
    energy_se = sd(energy_joules, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(variant = factor(variant, levels = c("BASELINE", "OPTIMIZED")))

# Determine if optimization improved or worsened energy for each benchmark
benchmark_improvement <- slope_data %>%
  pivot_wider(names_from = variant, values_from = c(energy_mean, energy_se)) %>%
  mutate(
    improved = energy_mean_OPTIMIZED < energy_mean_BASELINE,
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
median_delta_energy <- delta_pct %>%
  summarise(median_delta = median(delta_energy_pct, na.rm = TRUE)) %>%
  pull(median_delta)

# Define colors for benchmarks (using a colorblind-friendly palette)
benchmark_colors <- setNames(
  scales::hue_pal()(10),
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
    title = "Per-Benchmark Energy Consumption: Baseline vs Code Optimization",
    subtitle = "Each colored line connects paired runs from the same benchmark",
    x = NULL,
    y = "Energy Consumption (J, log scale)",
    caption = "Dashed lines = energy improved (↓); Solid lines = energy worsened (↑).\nBlack line = overall mean ± SE. Downward slope indicates energy reduction."
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
# SUMMARY STATISTICS TABLE
# =============================================================================
message("\nCreating summary statistics table...")

summary_stats <- delta_pct %>%
  summarise(
    across(starts_with("delta_"),
           list(
             median = ~median(., na.rm = TRUE),
             mean = ~mean(., na.rm = TRUE),
             sd = ~sd(., na.rm = TRUE),
             min = ~min(., na.rm = TRUE),
             max = ~max(., na.rm = TRUE),
             q25 = ~quantile(., 0.25, na.rm = TRUE),
             q75 = ~quantile(., 0.75, na.rm = TRUE)
           ),
           .names = "{.col}_{.fn}")
  ) %>%
  pivot_longer(everything(), names_to = "stat", values_to = "value") %>%
  separate(stat, into = c("metric", "stat_name"), sep = "_(?=[^_]+$)") %>%
  pivot_wider(names_from = stat_name, values_from = value) %>%
  mutate(
    metric = case_when(
      metric == "delta_energy_pct" ~ "Energy %",
      metric == "delta_time_pct" ~ "Time %",
      metric == "delta_cpu_pct" ~ "CPU %",
      metric == "delta_memory_pct" ~ "Memory %",
      metric == "delta_power_pct" ~ "Power %",
      TRUE ~ metric
    )
  )

write_csv(summary_stats, file.path(OUT_DIR, "rq1_summary_statistics.csv"))
message("Saved: rq1_summary_statistics.csv")

# =============================================================================
# FINAL REPORT
# =============================================================================
message("\n=============================================================================")
message("RQ1 ANALYSIS SUMMARY")
message("=============================================================================")

message("\nResearch Question 1:")
message("Does code optimization significantly impact energy consumption and")
message("related metrics compared to baseline?")

message("\n--- STATISTICAL TEST RESULTS ---")
for (i in 1:nrow(test_results)) {
  row <- test_results[i, ]
  message(sprintf("\n%s:", row$metric))
  message(sprintf("  Median Δ%%: %+.2f%%", row$median))
  message(sprintf("  Mean Δ%%: %+.2f%% (±%.2f)", row$mean, row$sd))
  message(sprintf("  Test: %s", row$test_used))
  message(sprintf("  p-value: %.4f %s", row$p_value, 
                  ifelse(row$significant, "(SIGNIFICANT)", "(not significant)")))
}

message("\n--- CONCLUSION ---")
sig_count <- sum(test_results$significant)
total_count <- nrow(test_results)
message(sprintf("%d out of %d metrics show significant deviation from zero (α = 0.05)", 
                sig_count, total_count))

if (test_results %>% filter(metric == "Energy %") %>% pull(significant)) {
  energy_median <- test_results %>% filter(metric == "Energy %") %>% pull(median)
  direction <- ifelse(energy_median < 0, "REDUCTION", "INCREASE")
  message(sprintf("\n✓ Energy consumption shows significant %s (median: %+.2f%%)", 
                  direction, energy_median))
} else {
  message("\n✗ Energy consumption does NOT show significant change")
}

message("\n--- OUTPUT FILES ---")
message("Plots:")
message("  - ", file.path(PLOTS_DIR, "plot1a_code_optimization_all_metrics_per_benchmark.png"))
message("  - ", file.path(PLOTS_DIR, "plot1b_sequential_minimum_optimization_guideline_impact.png"))
message("  - ", file.path(PLOTS_DIR, "plot1b_basic_string_guideline_impact.png"))
message("  - ", file.path(PLOTS_DIR, "plot2_paired_energy_comparison.png"))
message("\nData tables:")
message("  - ", file.path(OUT_DIR, "rq1_delta_pct_by_benchmark.csv"))
message("  - ", file.path(OUT_DIR, "rq1_statistical_tests.csv"))
message("  - ", file.path(OUT_DIR, "rq1_summary_statistics.csv"))

message("\n=============================================================================")
message("RQ1 Analysis Complete!")
message("=============================================================================\n")

