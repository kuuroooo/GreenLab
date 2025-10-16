#!/usr/bin/env Rscript

 # Ensure required packages are installed and loaded before use
 # Use specific packages to avoid heavy meta-package deps (e.g., textshaping/ragg via tidyverse)
 packages <- c("ggplot2", "readr", "dplyr", "stringr", "purrr", "tidyr", "scales", "forcats")
 for (p in packages) {
   if (!requireNamespace(p, quietly = TRUE)) {
     install.packages(p, repos = "https://cloud.r-project.org")
   }
 }
 suppressPackageStartupMessages({
   invisible(lapply(packages, function(p) library(p, character.only = TRUE)))
 })

# Configuration: absolute paths by default (override via env vars)
RESULTS_ROOT <- Sys.getenv("RESULTS_ROOT", unset = "/Users/kellywang/Desktop/GreenLab/RESULTS")
EXPERIMENTS_ROOT <- Sys.getenv("EXPERIMENTS_ROOT", unset = "/Users/kellywang/Desktop/GreenLab/EXPERIMENTS")
SUBDIR <- "benchmarks_energy_analysis"

OUTPUT_DIR <- file.path("/Users/kellywang/Desktop/GreenLab/analysis")
PLOTS_DIR <- file.path(OUTPUT_DIR, "plots")
OUT_DIR <- file.path(OUTPUT_DIR, "outputs")
dir.create(PLOTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

message("Using RESULTS_ROOT=", RESULTS_ROOT)
message("Using EXPERIMENTS_ROOT=", EXPERIMENTS_ROOT)

# Utilities ---------------------------------------------------------------

is_readable_file <- function(path) {
  isTRUE(file.exists(path)) && isTRUE(file.access(path, 4) == 0)
}

read_run_table_safe <- function(csv_path, results_set) {
  tryCatch(
    readr::read_csv(csv_path, show_col_types = FALSE) %>% mutate(results_set = results_set),
    error = function(e) {
      warning("Failed to read ", csv_path, ": ", conditionMessage(e))
      tibble()
    }
  )
}

parse_variant_and_benchmark_single <- function(python_file) {
  # Example inputs:
  #  - "kmeans/kmeans.py" -> benchmark=kmeans, variant=BASELINE
  #  - "kmeans/kmeans_G1.py" -> benchmark=kmeans, variant=G1
  #  - "two_hidden_layers_neural_network/two_hidden_layers_neural_network_ORIGINAL_FIXED.py" -> variant=ORIGINAL_FIXED
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
    # Fallback when naming deviates; treat as BASELINE
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
  # Look for: "Energy consumption in joules: <num> for <num> sec of execution."
  # Hyphen placed at end of class to avoid escape issues in R strings
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

# Ingest all run_table.csv across results_* --------------------------------

run_table_paths <- Sys.glob(file.path(RESULTS_ROOT, "results_*", SUBDIR, "run_table.csv"))
if (length(run_table_paths) == 0) {
  stop("No run_table.csv files found under ", file.path(RESULTS_ROOT, "results_*/", SUBDIR))
}

run_tables <- purrr::map_dfr(run_table_paths, function(p) {
  results_set <- basename(dirname(dirname(p))) # results_X
  read_run_table_safe(p, results_set)
})

if (!all(c("__run_id", "python_file") %in% names(run_tables))) {
  stop("run_table.csv missing required columns __run_id and python_file")
}

# Normalize columns and coerce numeric
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

# Derive benchmark and variant (vectorized)
bench_variant <- purrr::map_dfr(run_tables$python_file, function(pf) {
  res <- parse_variant_and_benchmark_single(pf)
  tibble(benchmark = res["benchmark"], variant = res["variant"])
})
run_tables <- bind_cols(run_tables, bench_variant)

# Parse logs as fallback for energy/time
run_tables <- run_tables %>%
  mutate(log_path = file.path(RESULTS_ROOT, results_set, SUBDIR, run_id, "energibridge.log"))

log_info <- purrr::map_dfr(run_tables$log_path, parse_energy_from_log)
run_tables <- bind_cols(run_tables, log_info)

run_tables <- run_tables %>%
  mutate(
    energy_joules = coalesce(total_energy_joules, energy_joules_log),
    exec_time_sec = coalesce(execution_time_seconds, exec_time_sec_log)
  )

# Persist raw combined runs
combined_raw_out <- file.path(OUT_DIR, "combined_runs_raw.csv")
readr::write_csv(run_tables, combined_raw_out)
message("Wrote ", combined_raw_out)

# Data Cleaning and Normalization -----------------------------------------
message("Cleaning and normalizing data...")

# 1. Detect outliers using IQR method (per benchmark-variant group)
outlier_detection <- run_tables %>%
  group_by(benchmark, variant) %>%
  mutate(
    q1 = quantile(energy_joules, 0.25, na.rm = TRUE),
    q3 = quantile(energy_joules, 0.75, na.rm = TRUE),
    iqr = q3 - q1,
    lower_bound = q1 - 3 * iqr,  # 3×IQR for extreme outliers
    upper_bound = q3 + 3 * iqr,
    is_outlier = energy_joules < lower_bound | energy_joules > upper_bound
  ) %>%
  ungroup() %>%
  select(-q1, -q3, -iqr, -lower_bound, -upper_bound)

# 2. Calculate Z-scores per benchmark (to compare across variants)
run_tables_with_zscore <- outlier_detection %>%
  group_by(benchmark) %>%
  mutate(
    energy_z_score = (energy_joules - mean(energy_joules, na.rm = TRUE)) / sd(energy_joules, na.rm = TRUE)
  ) %>%
  ungroup()

# 3. Min-max normalization per benchmark (0-1 scale)
run_tables_normalized <- run_tables_with_zscore %>%
  group_by(benchmark) %>%
  mutate(
    energy_min_bench = min(energy_joules, na.rm = TRUE),
    energy_max_bench = max(energy_joules, na.rm = TRUE),
    energy_range = energy_max_bench - energy_min_bench,
    energy_normalized = if_else(
      energy_range > 0,
      (energy_joules - energy_min_bench) / energy_range,
      0.5  # If all values are same, center at 0.5
    )
  ) %>%
  ungroup() %>%
  select(-energy_min_bench, -energy_max_bench, -energy_range)

# 4. Calculate % change vs baseline for each run
baseline_per_run <- run_tables_normalized %>%
  filter(variant == "BASELINE") %>%
  select(benchmark, results_set, baseline_energy_run = energy_joules)

run_tables_with_baseline <- run_tables_normalized %>%
  left_join(baseline_per_run, by = c("benchmark", "results_set")) %>%
  mutate(
    pct_vs_baseline_run = if_else(
      !is.na(baseline_energy_run) & baseline_energy_run > 0,
      (energy_joules - baseline_energy_run) / baseline_energy_run * 100,
      NA_real_
    ),
    efficiency_ratio = if_else(
      !is.na(baseline_energy_run) & baseline_energy_run > 0,
      energy_joules / baseline_energy_run,
      NA_real_
    )
  )

# 5. Calculate power if missing (using log-parsed execution time)
run_tables_cleaned <- run_tables_with_baseline %>%
  mutate(
    calculated_power = if_else(
      !is.na(energy_joules) & !is.na(exec_time_sec_log) & exec_time_sec_log > 0,
      energy_joules / exec_time_sec_log,
      NA_real_
    ),
    # Use calculated power if avg_power_watts is missing
    power_final = coalesce(avg_power_watts, calculated_power)
  )

# 6. Data quality metrics per variant
data_quality <- run_tables_cleaned %>%
  group_by(benchmark, variant) %>%
  summarise(
    n_runs = n(),
    n_outliers = sum(is_outlier, na.rm = TRUE),
    outlier_rate = n_outliers / n_runs * 100,
    energy_mean = mean(energy_joules, na.rm = TRUE),
    energy_sd = sd(energy_joules, na.rm = TRUE),
    energy_cv = energy_sd / energy_mean * 100,  # Coefficient of variation
    energy_range = max(energy_joules, na.rm = TRUE) - min(energy_joules, na.rm = TRUE),
    has_missing = any(is.na(energy_joules)),
    .groups = "drop"
  )

# 7. Extract outliers for review
outliers_data <- run_tables_cleaned %>%
  filter(is_outlier == TRUE) %>%
  select(
    results_set, run_id, benchmark, variant,
    energy_joules, avg_power_watts, cpu_usage_avg, memory_usage_avg,
    energy_z_score, pct_vs_baseline_run
  ) %>%
  arrange(benchmark, variant, desc(abs(energy_z_score)))

# Save cleaned and quality data
combined_clean_out <- file.path(OUT_DIR, "combined_runs_cleaned.csv")
quality_out <- file.path(OUT_DIR, "data_quality_summary.csv")
outliers_out <- file.path(OUT_DIR, "outliers_detected.csv")

readr::write_csv(run_tables_cleaned, combined_clean_out)
readr::write_csv(data_quality, quality_out)
readr::write_csv(outliers_data, outliers_out)

message("Wrote ", combined_clean_out)
message("Wrote ", quality_out, " (", nrow(data_quality), " benchmark-variant combinations)")
message("Wrote ", outliers_out, " (", nrow(outliers_data), " outliers detected)")

# Use cleaned data for further analysis
run_tables <- run_tables_cleaned

# Summaries ---------------------------------------------------------------

# Compute per-benchmark baseline means for deltas
baseline_means <- run_tables %>%
  filter(variant == "BASELINE") %>%
  group_by(benchmark) %>%
  summarise(
    baseline_energy_mean = mean(energy_joules, na.rm = TRUE),
    baseline_power_mean = mean(avg_power_watts, na.rm = TRUE),
    .groups = "drop"
  )

optimized_means <- run_tables %>%
  filter(variant == "OPTIMIZED") %>%
  group_by(benchmark) %>%
  summarise(
    optimized_energy_mean = mean(energy_joules, na.rm = TRUE),
    optimized_power_mean = mean(avg_power_watts, na.rm = TRUE),
    .groups = "drop"
  )

summary_by_variant <- run_tables %>%
  group_by(benchmark, variant) %>%
  summarise(
    n = n(),
    energy_mean = mean(energy_joules, na.rm = TRUE),
    energy_sd = sd(energy_joules, na.rm = TRUE),
    energy_median = median(energy_joules, na.rm = TRUE),
    energy_min = min(energy_joules, na.rm = TRUE),
    energy_max = max(energy_joules, na.rm = TRUE),
    power_mean = mean(avg_power_watts, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(baseline_means, by = "benchmark") %>%
  left_join(optimized_means, by = "benchmark") %>%
  mutate(
    delta_vs_baseline = energy_mean - baseline_energy_mean,
    pct_vs_baseline = (energy_mean - baseline_energy_mean) / baseline_energy_mean * 100,
    delta_vs_optimized = energy_mean - optimized_energy_mean,
    pct_vs_optimized = (energy_mean - optimized_energy_mean) / optimized_energy_mean * 100
  )

summary_out <- file.path(OUT_DIR, "summary_by_benchmark_variant.csv")
readr::write_csv(summary_by_variant, summary_out)
message("Wrote ", summary_out)

# Plots ------------------------------------------------------------------

safe_save <- function(plot, filename, width = 10, height = 6) {
  out_path <- file.path(PLOTS_DIR, filename)
  tryCatch({
    ggsave(out_path, plot = plot, width = width, height = height, dpi = 150)
    message("Saved ", out_path)
  }, error = function(e) warning("Failed to save ", out_path, ": ", conditionMessage(e)))
}

# 1) Per-benchmark energy bars by variant
benchmarks <- sort(unique(summary_by_variant$benchmark))

for (b in benchmarks) {
  df <- summary_by_variant %>% filter(benchmark == b) %>%
    mutate(variant = fct_reorder(variant, energy_mean))

  p <- ggplot(df, aes(x = variant, y = energy_mean, fill = variant)) +
    geom_col() +
    geom_errorbar(aes(ymin = pmax(0, energy_mean - energy_sd), ymax = energy_mean + energy_sd), width = 0.2) +
    scale_y_continuous(labels = label_number(scale_cut = cut_short_scale(), suffix = "J")) +
    scale_fill_brewer(palette = "Set3", guide = "none") +
    labs(title = paste0("Energy by Variant — ", b), x = "Variant", y = "Mean Energy (J)") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  safe_save(p, paste0("energy_by_variant_", b, ".png"))
}

# 2) Heatmap of % change vs baseline (rows: benchmark, cols: variant)
heat_df <- summary_by_variant %>%
  mutate(variant = factor(variant)) %>%
  select(benchmark, variant, pct_vs_baseline)

p_heat <- ggplot(heat_df, aes(x = variant, y = benchmark, fill = pct_vs_baseline)) +
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdYlGn", direction = -1, na.value = "grey90",
                       labels = function(x) sprintf("%+.0f%%", x)) +
  labs(title = "% Change vs Baseline (Mean Energy)", x = "Variant", y = "Benchmark", fill = "% vs Baseline") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

safe_save(p_heat, "heatmap_pct_vs_baseline.png", width = 12, height = 8)

# 3) Boxplot of per-run % change vs baseline across all benchmarks
# (pct_vs_baseline_run already calculated during data cleaning)
p_box <- run_tables %>%
  filter(!is.na(pct_vs_baseline_run), is.finite(pct_vs_baseline_run)) %>%
  mutate(variant = fct_reorder(variant, pct_vs_baseline_run, .fun = median, na.rm = TRUE)) %>%
  ggplot(aes(x = variant, y = pct_vs_baseline_run, fill = variant)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_y_continuous(labels = function(x) sprintf("%+.0f%%", x)) +
  scale_fill_brewer(palette = "Set3", guide = "none") +
  labs(title = "% Change vs Baseline by Variant (Per-Run)", x = "Variant", y = "% vs Baseline (per run)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

safe_save(p_box, "boxplot_pct_vs_baseline_per_run.png", width = 12, height = 7)

# 4) Scatter: Energy vs Avg Power (colored by variant)
p_scatter <- run_tables %>%
  ggplot(aes(x = avg_power_watts, y = energy_joules, color = variant)) +
  geom_point(alpha = 0.8) +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale(), suffix = "J")) +
  scale_x_continuous(labels = label_number(scale_cut = cut_short_scale(), suffix = "W")) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  labs(title = "Energy vs Average Power", x = "Average Power (W)", y = "Energy (J)") +
  theme_minimal(base_size = 12)

safe_save(p_scatter, "scatter_energy_vs_avg_power.png", width = 9, height = 7)

# 5) Optional: per-benchmark line of pct_vs_baseline across variants
line_df <- summary_by_variant %>%
  mutate(variant = factor(variant))

p_line <- ggplot(line_df, aes(x = variant, y = pct_vs_baseline, group = benchmark, color = benchmark)) +
  geom_line(alpha = 0.6) +
  geom_point(size = 1.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_y_continuous(labels = function(x) sprintf("%+.0f%%", x)) +
  labs(title = "% Change vs Baseline by Variant (Benchmark Lines)", x = "Variant", y = "% vs Baseline") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

safe_save(p_line, "lines_pct_vs_baseline_by_benchmark.png", width = 12, height = 7)

message("Done. Outputs under:")
message(" - ", PLOTS_DIR)
message(" - ", OUT_DIR)


