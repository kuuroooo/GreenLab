#!/usr/bin/env Rscript
# =============================================================================
# Comprehensive Data Exploration for Energy Benchmark Experiments
# =============================================================================
# This script performs exploratory data analysis on energy consumption data
# from multiple experimental runs, including descriptive statistics, outlier
# detection, and visualizations.

# Install and load required packages
packages <- c("ggplot2", "readr", "dplyr", "stringr", "purrr", "tidyr", 
              "scales", "forcats", "ggridges", "gridExtra", "knitr")
for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}
suppressPackageStartupMessages({
  invisible(lapply(packages, function(p) library(p, character.only = TRUE)))
})

# =============================================================================
# CONFIGURATION - Updated paths to use actual workspace
# =============================================================================
WORKSPACE_ROOT <- "/Users/kellywang/Documents/compSci/p1/greenLab/GreenLab"
RESULTS_ROOT <- file.path(WORKSPACE_ROOT, "RESULTS")
EXPERIMENTS_ROOT <- file.path(WORKSPACE_ROOT, "EXPERIMENTS")
SUBDIR <- "benchmarks_energy_analysis"

OUTPUT_DIR <- file.path(WORKSPACE_ROOT, "analysis")
PLOTS_DIR <- file.path(OUTPUT_DIR, "plots")
OUT_DIR <- file.path(OUTPUT_DIR, "outputs")
dir.create(PLOTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

message("=============================================================================")
message("Data Exploration Script")
message("=============================================================================")
message("Using RESULTS_ROOT: ", RESULTS_ROOT)
message("Using EXPERIMENTS_ROOT: ", EXPERIMENTS_ROOT)
message("Output directory: ", OUTPUT_DIR)
message("")

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

is_readable_file <- function(path) {
  isTRUE(file.exists(path)) && isTRUE(file.access(path, 4) == 0)
}

read_run_table_safe <- function(csv_path, results_set) {
  tryCatch(
    readr::read_csv(csv_path, show_col_types = FALSE) %>% 
      mutate(results_set = results_set),
    error = function(e) {
      warning("Failed to read ", csv_path, ": ", conditionMessage(e))
      tibble()
    }
  )
}

parse_variant_and_benchmark <- function(python_file) {
  path_no_ext <- str_replace(python_file, "\\.py$", "")
  parts <- str_split(path_no_ext, "/")[[1]]
  num_parts <- length(parts)
  file_base <- parts[num_parts]
  benchmark <- if (num_parts >= 2) parts[num_parts - 1] else file_base

  if (file_base == benchmark) {
    variant <- "BASELINE"
  } else if (str_detect(file_base, paste0("^", benchmark, "_"))) {
    variant <- str_replace(file_base, paste0("^", benchmark, "_"), "")
  } else {
    variant <- "BASELINE"
  }
  c(benchmark = benchmark, variant = variant)
}

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(x))
}

parse_energy_from_log <- function(log_path) {
  if (!is_readable_file(log_path)) return(tibble(energy_joules_log = NA_real_, exec_time_sec_log = NA_real_))
  lines <- tryCatch(readLines(log_path, warn = FALSE), error = function(e) character())
  if (length(lines) == 0) return(tibble(energy_joules_log = NA_real_, exec_time_sec_log = NA_real_))
  
  re <- "Energy consumption in joules: ([0-9eE+.-]+) for ([0-9eE+.-]+) sec"
  m <- stringr::str_match(lines, re)
  m <- m[!is.na(m[,1]), , drop = FALSE]
  if (nrow(m) == 0) return(tibble(energy_joules_log = NA_real_, exec_time_sec_log = NA_real_))
  last <- m[nrow(m), , drop = FALSE]
  tibble(
    energy_joules_log = safe_numeric(last[,2]),
    exec_time_sec_log = safe_numeric(last[,3])
  )
}

# Calculate comprehensive descriptive statistics
calc_descriptive_stats <- function(data, var_name) {
  x <- data[[var_name]]
  x <- x[!is.na(x) & is.finite(x)]
  
  if (length(x) == 0) {
    return(tibble(
      n = 0, mean = NA, median = NA, sd = NA, variance = NA, 
      min = NA, max = NA, q25 = NA, q75 = NA, iqr = NA, 
      cv = NA, skewness = NA, kurtosis = NA
    ))
  }
  
  mean_val <- mean(x)
  sd_val <- sd(x)
  
  # Skewness and kurtosis
  n <- length(x)
  skewness <- if (n > 2 && sd_val > 0) {
    sum((x - mean_val)^3) / (n * sd_val^3)
  } else NA
  
  kurtosis <- if (n > 3 && sd_val > 0) {
    sum((x - mean_val)^4) / (n * sd_val^4) - 3
  } else NA
  
  tibble(
    n = n,
    mean = mean_val,
    median = median(x),
    sd = sd_val,
    variance = var(x),
    min = min(x),
    max = max(x),
    q25 = quantile(x, 0.25),
    q75 = quantile(x, 0.75),
    iqr = IQR(x),
    cv = if (mean_val != 0) (sd_val / mean_val * 100) else NA,
    skewness = skewness,
    kurtosis = kurtosis
  )
}

# =============================================================================
# DATA INGESTION
# =============================================================================
message("Step 1: Ingesting data from RESULTS folder...")

# Find all run_table.csv files across results_* directories
run_table_paths <- Sys.glob(file.path(RESULTS_ROOT, "results_*", SUBDIR, "run_table.csv"))
message("Found ", length(run_table_paths), " run_table.csv files")

if (length(run_table_paths) == 0) {
  stop("No run_table.csv files found under ", file.path(RESULTS_ROOT, "results_*/", SUBDIR))
}

# Read and combine all run tables
run_tables <- purrr::map_dfr(run_table_paths, function(p) {
  results_set <- basename(dirname(dirname(p)))
  read_run_table_safe(p, results_set)
})

message("Total records loaded: ", nrow(run_tables))

# Check required columns
if (!all(c("__run_id", "python_file") %in% names(run_tables))) {
  stop("run_table.csv missing required columns __run_id and python_file")
}

# =============================================================================
# DATA PREPARATION
# =============================================================================
message("\nStep 2: Preparing and cleaning data...")

# Normalize columns and convert to numeric
run_tables <- run_tables %>%
  rename(run_id = `__run_id`) %>%
  mutate(
    total_energy_joules = safe_numeric(total_energy_joules),
    execution_time_seconds = safe_numeric(execution_time_seconds),
    avg_power_watts = safe_numeric(avg_power_watts),
    max_power_watts = safe_numeric(max_power_watts),
    min_power_watts = safe_numeric(min_power_watts),
    cpu_usage_avg = safe_numeric(cpu_usage_avg),
    memory_usage_avg = safe_numeric(memory_usage_avg)
  )

# Parse benchmark and variant from file path
bench_variant <- purrr::map_dfr(run_tables$python_file, function(pf) {
  res <- parse_variant_and_benchmark(pf)
  tibble(benchmark = res["benchmark"], variant = res["variant"])
})
run_tables <- bind_cols(run_tables, bench_variant)

# Filter out CLEAN and ORIGINAL_FIXED variants
message("Filtering out CLEAN and ORIGINAL_FIXED variants...")
n_before <- nrow(run_tables)
run_tables <- run_tables %>%
  filter(!variant %in% c("CLEAN", "ORIGINAL_FIXED"))
n_after <- nrow(run_tables)
message("Removed ", n_before - n_after, " rows (", n_before - n_after, " measurements)")
message("Remaining: ", n_after, " measurements")

# Parse energy from logs as fallback
run_tables <- run_tables %>%
  mutate(log_path = file.path(RESULTS_ROOT, results_set, SUBDIR, run_id, "energibridge.log"))

log_info <- purrr::map_dfr(run_tables$log_path, parse_energy_from_log)
run_tables <- bind_cols(run_tables, log_info)

# Coalesce energy and execution time from multiple sources
run_tables <- run_tables %>%
  mutate(
    energy_joules = coalesce(total_energy_joules, energy_joules_log),
    exec_time_sec = coalesce(execution_time_seconds, exec_time_sec_log),
    # Calculate power if missing
    power_watts = if_else(
      is.na(avg_power_watts) & !is.na(energy_joules) & !is.na(exec_time_sec) & exec_time_sec > 0,
      energy_joules / exec_time_sec,
      avg_power_watts
    )
  )

message("Data prepared successfully")
message("Unique benchmarks: ", n_distinct(run_tables$benchmark))
message("Unique variants: ", n_distinct(run_tables$variant))
message("Total experimental runs: ", n_distinct(run_tables$results_set))

# =============================================================================
# OUTLIER DETECTION
# =============================================================================
message("\nStep 3: Detecting outliers...")

# Detect outliers using IQR method (per benchmark-variant group)
run_tables <- run_tables %>%
  group_by(benchmark, variant) %>%
  mutate(
    q1_energy = quantile(energy_joules, 0.25, na.rm = TRUE),
    q3_energy = quantile(energy_joules, 0.75, na.rm = TRUE),
    iqr_energy = q3_energy - q1_energy,
    lower_bound_energy = q1_energy - 1.5 * iqr_energy,
    upper_bound_energy = q3_energy + 1.5 * iqr_energy,
    is_outlier_moderate = energy_joules < lower_bound_energy | energy_joules > upper_bound_energy,
    # Extreme outliers (3×IQR)
    lower_bound_extreme = q1_energy - 3 * iqr_energy,
    upper_bound_extreme = q3_energy + 3 * iqr_energy,
    is_outlier_extreme = energy_joules < lower_bound_extreme | energy_joules > upper_bound_extreme
  ) %>%
  ungroup() %>%
  select(-q1_energy, -q3_energy, -iqr_energy, -lower_bound_energy, -upper_bound_energy,
         -lower_bound_extreme, -upper_bound_extreme)

# Calculate Z-scores per benchmark (to compare across variants)
run_tables <- run_tables %>%
  group_by(benchmark) %>%
  mutate(
    energy_z_score = (energy_joules - mean(energy_joules, na.rm = TRUE)) / 
                     sd(energy_joules, na.rm = TRUE)
  ) %>%
  ungroup()

# Mark outliers based on Z-score (|z| > 3 is extreme)
run_tables <- run_tables %>%
  mutate(
    is_outlier_zscore = abs(energy_z_score) > 2.5,
    is_outlier_extreme_zscore = abs(energy_z_score) > 3
  )

n_outliers_moderate <- sum(run_tables$is_outlier_moderate, na.rm = TRUE)
n_outliers_extreme <- sum(run_tables$is_outlier_extreme, na.rm = TRUE)
message("Moderate outliers detected (1.5×IQR): ", n_outliers_moderate)
message("Extreme outliers detected (3×IQR): ", n_outliers_extreme)

# =============================================================================
# COMPREHENSIVE DESCRIPTIVE STATISTICS
# =============================================================================
message("\nStep 4: Calculating descriptive statistics...")

# Overall statistics across all data
overall_stats <- tibble(
  metric = c("Energy (J)", "Execution Time (s)", "Avg Power (W)", 
             "CPU Usage (%)", "Memory Usage (MB)"),
  variable = c("energy_joules", "exec_time_sec", "power_watts", 
               "cpu_usage_avg", "memory_usage_avg")
) %>%
  rowwise() %>%
  mutate(
    stats = list(calc_descriptive_stats(run_tables, variable))
  ) %>%
  ungroup() %>%
  unnest(stats) %>%
  select(-variable)

overall_stats_out <- file.path(OUT_DIR, "overall_descriptive_statistics.csv")
readr::write_csv(overall_stats, overall_stats_out)
message("Saved: ", overall_stats_out)

# Statistics by benchmark
stats_by_benchmark <- run_tables %>%
  group_by(benchmark) %>%
  summarise(
    n_runs = n(),
    n_variants = n_distinct(variant),
    # Energy statistics
    energy_mean = mean(energy_joules, na.rm = TRUE),
    energy_median = median(energy_joules, na.rm = TRUE),
    energy_sd = sd(energy_joules, na.rm = TRUE),
    energy_variance = var(energy_joules, na.rm = TRUE),
    energy_min = min(energy_joules, na.rm = TRUE),
    energy_max = max(energy_joules, na.rm = TRUE),
    energy_q25 = quantile(energy_joules, 0.25, na.rm = TRUE),
    energy_q75 = quantile(energy_joules, 0.75, na.rm = TRUE),
    energy_iqr = IQR(energy_joules, na.rm = TRUE),
    energy_cv = (energy_sd / energy_mean * 100),
    # Execution time statistics
    time_mean = mean(exec_time_sec, na.rm = TRUE),
    time_median = median(exec_time_sec, na.rm = TRUE),
    time_sd = sd(exec_time_sec, na.rm = TRUE),
    time_min = min(exec_time_sec, na.rm = TRUE),
    time_max = max(exec_time_sec, na.rm = TRUE),
    # Power statistics
    power_mean = mean(power_watts, na.rm = TRUE),
    power_median = median(power_watts, na.rm = TRUE),
    power_sd = sd(power_watts, na.rm = TRUE),
    power_min = min(power_watts, na.rm = TRUE),
    power_max = max(power_watts, na.rm = TRUE),
    # CPU statistics
    cpu_mean = mean(cpu_usage_avg, na.rm = TRUE),
    cpu_median = median(cpu_usage_avg, na.rm = TRUE),
    cpu_sd = sd(cpu_usage_avg, na.rm = TRUE),
    # Memory statistics
    memory_mean = mean(memory_usage_avg, na.rm = TRUE),
    memory_median = median(memory_usage_avg, na.rm = TRUE),
    memory_sd = sd(memory_usage_avg, na.rm = TRUE),
    .groups = "drop"
  )

stats_by_benchmark_out <- file.path(OUT_DIR, "stats_by_benchmark.csv")
readr::write_csv(stats_by_benchmark, stats_by_benchmark_out)
message("Saved: ", stats_by_benchmark_out)

# Statistics by variant
stats_by_variant <- run_tables %>%
  group_by(variant) %>%
  summarise(
    n_runs = n(),
    n_benchmarks = n_distinct(benchmark),
    # Energy statistics
    energy_mean = mean(energy_joules, na.rm = TRUE),
    energy_median = median(energy_joules, na.rm = TRUE),
    energy_sd = sd(energy_joules, na.rm = TRUE),
    energy_variance = var(energy_joules, na.rm = TRUE),
    energy_min = min(energy_joules, na.rm = TRUE),
    energy_max = max(energy_joules, na.rm = TRUE),
    energy_q25 = quantile(energy_joules, 0.25, na.rm = TRUE),
    energy_q75 = quantile(energy_joules, 0.75, na.rm = TRUE),
    energy_iqr = IQR(energy_joules, na.rm = TRUE),
    energy_cv = (energy_sd / energy_mean * 100),
    # Power statistics
    power_mean = mean(power_watts, na.rm = TRUE),
    power_median = median(power_watts, na.rm = TRUE),
    power_sd = sd(power_watts, na.rm = TRUE),
    # CPU statistics
    cpu_mean = mean(cpu_usage_avg, na.rm = TRUE),
    cpu_median = median(cpu_usage_avg, na.rm = TRUE),
    cpu_sd = sd(cpu_usage_avg, na.rm = TRUE),
    # Memory statistics
    memory_mean = mean(memory_usage_avg, na.rm = TRUE),
    memory_median = median(memory_usage_avg, na.rm = TRUE),
    memory_sd = sd(memory_usage_avg, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(energy_mean))

stats_by_variant_out <- file.path(OUT_DIR, "stats_by_variant.csv")
readr::write_csv(stats_by_variant, stats_by_variant_out)
message("Saved: ", stats_by_variant_out)

# Detailed statistics by benchmark and variant (most important!)
stats_by_benchmark_variant <- run_tables %>%
  group_by(benchmark, variant) %>%
  summarise(
    n_runs = n(),
    # Energy statistics
    energy_mean = mean(energy_joules, na.rm = TRUE),
    energy_median = median(energy_joules, na.rm = TRUE),
    energy_sd = sd(energy_joules, na.rm = TRUE),
    energy_variance = var(energy_joules, na.rm = TRUE),
    energy_min = min(energy_joules, na.rm = TRUE),
    energy_max = max(energy_joules, na.rm = TRUE),
    energy_q25 = quantile(energy_joules, 0.25, na.rm = TRUE),
    energy_q75 = quantile(energy_joules, 0.75, na.rm = TRUE),
    energy_iqr = IQR(energy_joules, na.rm = TRUE),
    energy_cv = (energy_sd / energy_mean * 100),
    # Execution time
    time_mean = mean(exec_time_sec, na.rm = TRUE),
    time_median = median(exec_time_sec, na.rm = TRUE),
    time_sd = sd(exec_time_sec, na.rm = TRUE),
    time_variance = var(exec_time_sec, na.rm = TRUE),
    # Power
    power_mean = mean(power_watts, na.rm = TRUE),
    power_median = median(power_watts, na.rm = TRUE),
    power_sd = sd(power_watts, na.rm = TRUE),
    power_variance = var(power_watts, na.rm = TRUE),
    # CPU
    cpu_mean = mean(cpu_usage_avg, na.rm = TRUE),
    cpu_median = median(cpu_usage_avg, na.rm = TRUE),
    cpu_sd = sd(cpu_usage_avg, na.rm = TRUE),
    cpu_variance = var(cpu_usage_avg, na.rm = TRUE),
    # Memory
    memory_mean = mean(memory_usage_avg, na.rm = TRUE),
    memory_median = median(memory_usage_avg, na.rm = TRUE),
    memory_sd = sd(memory_usage_avg, na.rm = TRUE),
    memory_variance = var(memory_usage_avg, na.rm = TRUE),
    # Outliers
    n_outliers_moderate = sum(is_outlier_moderate, na.rm = TRUE),
    n_outliers_extreme = sum(is_outlier_extreme, na.rm = TRUE),
    outlier_rate = (n_outliers_moderate / n()) * 100,
    .groups = "drop"
  )

stats_by_benchmark_variant_out <- file.path(OUT_DIR, "stats_by_benchmark_variant.csv")
readr::write_csv(stats_by_benchmark_variant, stats_by_benchmark_variant_out)
message("Saved: ", stats_by_benchmark_variant_out)

# Comparison with baseline
baseline_stats <- run_tables %>%
  filter(variant == "BASELINE") %>%
  group_by(benchmark) %>%
  summarise(
    baseline_energy_mean = mean(energy_joules, na.rm = TRUE),
    baseline_energy_median = median(energy_joules, na.rm = TRUE),
    baseline_power_mean = mean(power_watts, na.rm = TRUE),
    baseline_time_mean = mean(exec_time_sec, na.rm = TRUE),
    .groups = "drop"
  )

stats_vs_baseline <- stats_by_benchmark_variant %>%
  left_join(baseline_stats, by = "benchmark") %>%
  mutate(
    energy_pct_vs_baseline = ((energy_mean - baseline_energy_mean) / baseline_energy_mean) * 100,
    energy_median_pct_vs_baseline = ((energy_median - baseline_energy_median) / baseline_energy_median) * 100,
    power_pct_vs_baseline = ((power_mean - baseline_power_mean) / baseline_power_mean) * 100,
    time_pct_vs_baseline = ((time_mean - baseline_time_mean) / baseline_time_mean) * 100
  ) %>%
  select(benchmark, variant, n_runs, 
         energy_mean, energy_median, energy_sd, energy_iqr, energy_cv,
         baseline_energy_mean, energy_pct_vs_baseline, energy_median_pct_vs_baseline,
         power_mean, power_pct_vs_baseline,
         time_mean, time_pct_vs_baseline,
         cpu_mean, memory_mean,
         outlier_rate)

stats_vs_baseline_out <- file.path(OUT_DIR, "stats_comparison_vs_baseline.csv")
readr::write_csv(stats_vs_baseline, stats_vs_baseline_out)
message("Saved: ", stats_vs_baseline_out)

# Data quality summary
data_quality <- run_tables %>%
  group_by(benchmark, variant) %>%
  summarise(
    n_runs = n(),
    missing_energy = sum(is.na(energy_joules)),
    missing_power = sum(is.na(power_watts)),
    missing_cpu = sum(is.na(cpu_usage_avg)),
    missing_memory = sum(is.na(memory_usage_avg)),
    n_outliers_moderate = sum(is_outlier_moderate, na.rm = TRUE),
    n_outliers_extreme = sum(is_outlier_extreme, na.rm = TRUE),
    outlier_rate_pct = (n_outliers_moderate / n_runs) * 100,
    energy_cv = (sd(energy_joules, na.rm = TRUE) / mean(energy_joules, na.rm = TRUE)) * 100,
    stability = case_when(
      energy_cv < 5 ~ "Very Stable",
      energy_cv < 10 ~ "Stable",
      energy_cv < 20 ~ "Moderate",
      energy_cv < 30 ~ "Variable",
      TRUE ~ "Highly Variable"
    ),
    .groups = "drop"
  )

data_quality_out <- file.path(OUT_DIR, "data_quality_report.csv")
readr::write_csv(data_quality, data_quality_out)
message("Saved: ", data_quality_out)

# Outlier details
outliers_detailed <- run_tables %>%
  filter(is_outlier_moderate == TRUE) %>%
  select(results_set, run_id, benchmark, variant,
         energy_joules, power_watts, cpu_usage_avg, memory_usage_avg,
         energy_z_score, is_outlier_extreme) %>%
  arrange(benchmark, variant, desc(abs(energy_z_score)))

outliers_out <- file.path(OUT_DIR, "outliers_detailed.csv")
readr::write_csv(outliers_detailed, outliers_out)
message("Saved: ", outliers_out)

# Save cleaned dataset
cleaned_data_out <- file.path(OUT_DIR, "cleaned_experiment_data.csv")
readr::write_csv(run_tables, cleaned_data_out)
message("Saved: ", cleaned_data_out)

# =============================================================================
# VISUALIZATIONS
# =============================================================================
message("\nStep 5: Creating visualizations...")

safe_save <- function(plot, filename, width = 10, height = 6) {
  out_path <- file.path(PLOTS_DIR, filename)
  tryCatch({
    ggsave(out_path, plot = plot, width = width, height = height, dpi = 150)
    message("  Saved: ", filename)
  }, error = function(e) warning("Failed to save ", filename, ": ", conditionMessage(e)))
}

# 3. Boxplot of energy by variant
p3 <- run_tables %>%
  filter(!is.na(energy_joules)) %>%
  mutate(variant = fct_reorder(variant, energy_joules, .fun = median)) %>%
  ggplot(aes(x = variant, y = energy_joules, fill = variant)) +
  geom_boxplot(outlier.alpha = 0.5, outlier.size = 1) +
  scale_y_log10(labels = label_number(scale_cut = cut_short_scale(), suffix = "J")) +
  scale_fill_brewer(palette = "Set3") +
  labs(title = "Energy Consumption by Variant (All Benchmarks)",
       x = "Variant", y = "Energy (J, log scale)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
safe_save(p3, "03_energy_boxplot_by_variant.png", width = 12, height = 7)

# 4. Boxplot of energy by benchmark
p4 <- run_tables %>%
  filter(!is.na(energy_joules)) %>%
  mutate(benchmark = fct_reorder(benchmark, energy_joules, .fun = median)) %>%
  ggplot(aes(x = benchmark, y = energy_joules, fill = benchmark)) +
  geom_boxplot(outlier.alpha = 0.5, outlier.size = 1) +
  scale_y_log10(labels = label_number(scale_cut = cut_short_scale(), suffix = "J")) +
  scale_fill_brewer(palette = "Set3") +
  labs(title = "Energy Consumption by Benchmark (All Variants)",
       x = "Benchmark", y = "Energy (J, log scale)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
safe_save(p4, "04_energy_boxplot_by_benchmark.png", width = 14, height = 7)

# 5. Faceted boxplot: energy by variant for each benchmark
p5 <- run_tables %>%
  filter(!is.na(energy_joules)) %>%
  ggplot(aes(x = variant, y = energy_joules, fill = variant)) +
  geom_boxplot(outlier.size = 0.8) +
  facet_wrap(~benchmark, scales = "free_y", ncol = 3) +
  scale_y_log10(labels = label_number(scale_cut = cut_short_scale(), suffix = "J")) +
  scale_fill_brewer(palette = "Set3") +
  labs(title = "Energy Consumption by Variant for Each Benchmark",
       x = "Variant", y = "Energy (J, log scale)") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none",
        strip.text = element_text(size = 9))
safe_save(p5, "05_energy_boxplot_faceted.png", width = 16, height = 12)

# 6. Violin plot with quartiles
p6 <- run_tables %>%
  filter(!is.na(energy_joules)) %>%
  mutate(variant = fct_reorder(variant, energy_joules, .fun = median)) %>%
  ggplot(aes(x = variant, y = energy_joules, fill = variant)) +
  geom_violin(alpha = 0.6, draw_quantiles = c(0.25, 0.5, 0.75)) +
  scale_y_log10(labels = label_number(scale_cut = cut_short_scale(), suffix = "J")) +
  scale_fill_brewer(palette = "Set3") +
  labs(title = "Energy Distribution by Variant (Violin Plot with Quartiles)",
       x = "Variant", y = "Energy (J, log scale)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
safe_save(p6, "06_energy_violin_by_variant.png", width = 12, height = 7)

# 10. Heatmap: mean energy by benchmark and variant
heatmap_data <- stats_by_benchmark_variant %>%
  select(benchmark, variant, energy_mean)

p10 <- ggplot(heatmap_data, aes(x = variant, y = benchmark, fill = energy_mean)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "plasma", trans = "log10",
                       labels = label_number(scale_cut = cut_short_scale(), suffix = "J")) +
  labs(title = "Mean Energy Consumption Heatmap",
       x = "Variant", y = "Benchmark", fill = "Mean Energy (J)") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
safe_save(p10, "10_energy_heatmap_mean.png", width = 12, height = 10)

# 15. Summary comparison: baseline vs optimized vs guidelines
summary_comparison <- stats_by_benchmark_variant %>%
  mutate(
    variant_type = case_when(
      variant == "BASELINE" ~ "Baseline",
      variant == "OPTIMIZED" ~ "Optimized",
      TRUE ~ "Guideline"
    )
  ) %>%
  group_by(benchmark, variant_type) %>%
  summarise(
    energy_mean = mean(energy_mean, na.rm = TRUE),
    energy_median = median(energy_median, na.rm = TRUE),
    .groups = "drop"
  )

p15 <- ggplot(summary_comparison, aes(x = benchmark, y = energy_mean, fill = variant_type)) +
  geom_col(position = "dodge") +
  scale_y_log10(labels = label_number(scale_cut = cut_short_scale(), suffix = "J")) +
  scale_fill_manual(values = c("Baseline" = "#e74c3c", "Optimized" = "#27ae60", "Guideline" = "#3498db")) +
  labs(title = "Mean Energy: Baseline vs Optimized vs Guidelines",
       x = "Benchmark", y = "Mean Energy (J, log scale)", fill = "Type") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
safe_save(p15, "15_energy_comparison_by_type.png", width = 14, height = 7)

# =============================================================================
# SUMMARY REPORT
# =============================================================================
message("\n=============================================================================")
message("SUMMARY REPORT")
message("=============================================================================")
message("\nData Overview:")
message("  Total runs: ", nrow(run_tables))
message("  Unique benchmarks: ", n_distinct(run_tables$benchmark))
message("  Unique variants: ", n_distinct(run_tables$variant))
message("  Experimental runs (results_*): ", n_distinct(run_tables$results_set))

message("\nEnergy Statistics (Overall):")
overall_energy <- overall_stats %>% filter(metric == "Energy (J)")
message("  Mean: ", sprintf("%.3f J", overall_energy$mean))
message("  Median: ", sprintf("%.3f J", overall_energy$median))
message("  SD: ", sprintf("%.3f J", overall_energy$sd))
message("  Range: [", sprintf("%.3f", overall_energy$min), ", ", 
        sprintf("%.3f", overall_energy$max), "] J")
message("  IQR: ", sprintf("%.3f J", overall_energy$iqr))
message("  CV: ", sprintf("%.2f%%", overall_energy$cv))

message("\nOutliers:")
message("  Moderate outliers (1.5×IQR): ", n_outliers_moderate, 
        " (", sprintf("%.2f%%", n_outliers_moderate/nrow(run_tables)*100), ")")
message("  Extreme outliers (3×IQR): ", n_outliers_extreme,
        " (", sprintf("%.2f%%", n_outliers_extreme/nrow(run_tables)*100), ")")

message("\nData Quality:")
high_var <- data_quality %>% filter(energy_cv > 20) %>% nrow()
message("  High variability combinations (CV > 20%): ", high_var)
message("  Most stable variant: ", 
        (stats_by_variant %>% arrange(energy_cv) %>% slice(1))$variant)
message("  Most variable variant: ", 
        (stats_by_variant %>% arrange(desc(energy_cv)) %>% slice(1))$variant)

message("\nOutput Files:")
message("  - ", OUT_DIR)
message("    • overall_descriptive_statistics.csv")
message("    • stats_by_benchmark.csv")
message("    • stats_by_variant.csv")
message("    • stats_by_benchmark_variant.csv")
message("    • stats_comparison_vs_baseline.csv")
message("    • data_quality_report.csv")
message("    • outliers_detailed.csv")
message("    • cleaned_experiment_data.csv")
message("\n  - ", PLOTS_DIR)
message("    • 6 visualization files (03-06, 10, 15)")

message("\n=============================================================================")
message("Data exploration complete!")
message("=============================================================================\n")

