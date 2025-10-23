#!/usr/bin/env Rscript
# =============================================================================
# =============================================================================
# =============================================================================

packages <- c("ggplot2", "readr", "dplyr", "stringr", "purrr", "tidyr", 
              "scales", "forcats", "ggridges", "gridExtra", "knitr", "ggrepel")
for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}
suppressPackageStartupMessages({
  invisible(lapply(packages, function(p) library(p, character.only = TRUE)))
})


WORKSPACE_ROOT <- "/Users/kellywang/Documents/compSci/p1/greenLab/GreenLab"
RESULTS_ROOT <- file.path(WORKSPACE_ROOT, "RESULTS")
EXPERIMENTS_ROOT <- file.path(WORKSPACE_ROOT, "EXPERIMENTS")
SUBDIR <- "benchmarks_energy_analysis"
OUTPUT_DIR <- file.path(WORKSPACE_ROOT, "analysis")
PLOTS_DIR <- file.path(OUTPUT_DIR, "plots")
OUT_DIR <- file.path(OUTPUT_DIR, "outputs")

dir.create(PLOTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# =============================================================================

# =============================================================================
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
# =============================================================================

run_table_paths <- Sys.glob(file.path(RESULTS_ROOT, "results_*", SUBDIR, "run_table.csv"))

if (length(run_table_paths) == 0) {
  stop("No run_table.csv files found under ", file.path(RESULTS_ROOT, "results_*/", SUBDIR))
}

run_tables <- purrr::map_dfr(run_table_paths, function(p) {
  results_set <- basename(dirname(dirname(p)))
  read_run_table_safe(p, results_set)
})


if (!all(c("__run_id", "python_file") %in% names(run_tables))) {
  stop("run_table.csv missing required columns __run_id and python_file")
}

# =============================================================================
# =============================================================================

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

bench_variant <- purrr::map_dfr(run_tables$python_file, function(pf) {
  res <- parse_variant_and_benchmark(pf)
  tibble(benchmark = res["benchmark"], variant = res["variant"])
})
run_tables <- bind_cols(run_tables, bench_variant)

n_before <- nrow(run_tables)
run_tables <- run_tables %>%
  filter(!variant %in% c("CLEAN", "ORIGINAL_FIXED"))
n_after <- nrow(run_tables)

run_tables <- run_tables %>%
  mutate(log_path = file.path(RESULTS_ROOT, results_set, SUBDIR, run_id, "energibridge.log"))

log_info <- purrr::map_dfr(run_tables$log_path, parse_energy_from_log)
run_tables <- bind_cols(run_tables, log_info)

run_tables <- run_tables %>%
  mutate(
    energy_joules = coalesce(total_energy_joules, energy_joules_log),
    exec_time_sec = coalesce(execution_time_seconds, exec_time_sec_log),
    power_watts = if_else(
      is.na(avg_power_watts) & !is.na(energy_joules) & !is.na(exec_time_sec) & exec_time_sec > 0,
      energy_joules / exec_time_sec,
      avg_power_watts
    )
  )


# =============================================================================
# =============================================================================

run_tables <- run_tables %>%
  group_by(benchmark, variant) %>%
  mutate(
    q1_energy = quantile(energy_joules, 0.25, na.rm = TRUE),
    q3_energy = quantile(energy_joules, 0.75, na.rm = TRUE),
    iqr_energy = q3_energy - q1_energy,
    lower_bound_energy = q1_energy - 1.5 * iqr_energy,
    upper_bound_energy = q3_energy + 1.5 * iqr_energy,
    is_outlier_moderate = energy_joules < lower_bound_energy | energy_joules > upper_bound_energy,
    lower_bound_extreme = q1_energy - 3 * iqr_energy,
    upper_bound_extreme = q3_energy + 3 * iqr_energy,
    is_outlier_extreme = energy_joules < lower_bound_extreme | energy_joules > upper_bound_extreme
  ) %>%
  ungroup() %>%
  select(-q1_energy, -q3_energy, -iqr_energy, -lower_bound_energy, -upper_bound_energy,
         -lower_bound_extreme, -upper_bound_extreme)

run_tables <- run_tables %>%
  group_by(benchmark) %>%
  mutate(
    energy_z_score = (energy_joules - mean(energy_joules, na.rm = TRUE)) / 
                     sd(energy_joules, na.rm = TRUE)
  ) %>%
  ungroup()

run_tables <- run_tables %>%
  mutate(
    is_outlier_zscore = abs(energy_z_score) > 2.5,
    is_outlier_extreme_zscore = abs(energy_z_score) > 3
  )

n_outliers_moderate <- sum(run_tables$is_outlier_moderate, na.rm = TRUE)
n_outliers_extreme <- sum(run_tables$is_outlier_extreme, na.rm = TRUE)

# =============================================================================
# =============================================================================

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

stats_by_benchmark <- run_tables %>%
  group_by(benchmark) %>%
  summarise(
    n_runs = n(),
    n_variants = n_distinct(variant),
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
    time_mean = mean(exec_time_sec, na.rm = TRUE),
    time_median = median(exec_time_sec, na.rm = TRUE),
    time_sd = sd(exec_time_sec, na.rm = TRUE),
    time_min = min(exec_time_sec, na.rm = TRUE),
    time_max = max(exec_time_sec, na.rm = TRUE),
    power_mean = mean(power_watts, na.rm = TRUE),
    power_median = median(power_watts, na.rm = TRUE),
    power_sd = sd(power_watts, na.rm = TRUE),
    power_min = min(power_watts, na.rm = TRUE),
    power_max = max(power_watts, na.rm = TRUE),
    cpu_mean = mean(cpu_usage_avg, na.rm = TRUE),
    cpu_median = median(cpu_usage_avg, na.rm = TRUE),
    cpu_sd = sd(cpu_usage_avg, na.rm = TRUE),
    memory_mean = mean(memory_usage_avg, na.rm = TRUE),
    memory_median = median(memory_usage_avg, na.rm = TRUE),
    memory_sd = sd(memory_usage_avg, na.rm = TRUE),
    .groups = "drop"
  )

stats_by_benchmark_out <- file.path(OUT_DIR, "stats_by_benchmark.csv")
readr::write_csv(stats_by_benchmark, stats_by_benchmark_out)

stats_by_variant <- run_tables %>%
  group_by(variant) %>%
  summarise(
    n_runs = n(),
    n_benchmarks = n_distinct(benchmark),
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
    power_mean = mean(power_watts, na.rm = TRUE),
    power_median = median(power_watts, na.rm = TRUE),
    power_sd = sd(power_watts, na.rm = TRUE),
    cpu_mean = mean(cpu_usage_avg, na.rm = TRUE),
    cpu_median = median(cpu_usage_avg, na.rm = TRUE),
    cpu_sd = sd(cpu_usage_avg, na.rm = TRUE),
    memory_mean = mean(memory_usage_avg, na.rm = TRUE),
    memory_median = median(memory_usage_avg, na.rm = TRUE),
    memory_sd = sd(memory_usage_avg, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(energy_mean))

stats_by_variant_out <- file.path(OUT_DIR, "stats_by_variant.csv")
readr::write_csv(stats_by_variant, stats_by_variant_out)

stats_by_benchmark_variant <- run_tables %>%
  group_by(benchmark, variant) %>%
  summarise(
    n_runs = n(),
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
    time_mean = mean(exec_time_sec, na.rm = TRUE),
    time_median = median(exec_time_sec, na.rm = TRUE),
    time_sd = sd(exec_time_sec, na.rm = TRUE),
    time_variance = var(exec_time_sec, na.rm = TRUE),
    power_mean = mean(power_watts, na.rm = TRUE),
    power_median = median(power_watts, na.rm = TRUE),
    power_sd = sd(power_watts, na.rm = TRUE),
    power_variance = var(power_watts, na.rm = TRUE),
    cpu_mean = mean(cpu_usage_avg, na.rm = TRUE),
    cpu_median = median(cpu_usage_avg, na.rm = TRUE),
    cpu_sd = sd(cpu_usage_avg, na.rm = TRUE),
    cpu_variance = var(cpu_usage_avg, na.rm = TRUE),
    memory_mean = mean(memory_usage_avg, na.rm = TRUE),
    memory_median = median(memory_usage_avg, na.rm = TRUE),
    memory_sd = sd(memory_usage_avg, na.rm = TRUE),
    memory_variance = var(memory_usage_avg, na.rm = TRUE),
    n_outliers_moderate = sum(is_outlier_moderate, na.rm = TRUE),
    n_outliers_extreme = sum(is_outlier_extreme, na.rm = TRUE),
    outlier_rate = (n_outliers_moderate / n()) * 100,
    .groups = "drop"
  )

stats_by_benchmark_variant_out <- file.path(OUT_DIR, "stats_by_benchmark_variant.csv")
readr::write_csv(stats_by_benchmark_variant, stats_by_benchmark_variant_out)

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

outliers_detailed <- run_tables %>%
  filter(is_outlier_moderate == TRUE) %>%
  select(results_set, run_id, benchmark, variant,
         energy_joules, power_watts, cpu_usage_avg, memory_usage_avg,
         energy_z_score, is_outlier_extreme) %>%
  arrange(benchmark, variant, desc(abs(energy_z_score)))

outliers_out <- file.path(OUT_DIR, "outliers_detailed.csv")
readr::write_csv(outliers_detailed, outliers_out)

cleaned_data_out <- file.path(OUT_DIR, "cleaned_experiment_data.csv")
readr::write_csv(run_tables, cleaned_data_out)

# =============================================================================
# =============================================================================

data <- read_csv(file.path(OUT_DIR, "cleaned_experiment_data.csv"), show_col_types = FALSE)

data <- data %>% filter(!variant %in% c("CLEAN", "ORIGINAL_FIXED"))

n_before <- nrow(data)
data <- data %>%
  filter(!variant %in% c("CLEAN", "ORIGINAL_FIXED"))
n_after <- nrow(data)


# =============================================================================
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

# =============================================================================
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

# =============================================================================
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

# =============================================================================
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

# =============================================================================
# =============================================================================

table5 <- data %>%
  group_by(benchmark, variant) %>%
  summarise(median_energy = median(energy_joules, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = variant, values_from = median_energy) %>%
  mutate(across(where(is.numeric), ~sprintf("%.2f", .)))

write_csv(table5, file.path(OUT_DIR, "table5_energy_matrix_median.csv"))

# =============================================================================
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

# =============================================================================
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

# =============================================================================
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

# =============================================================================
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

# =============================================================================
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

# =============================================================================
# =============================================================================


# =============================================================================
# =============================================================================

run_tables <- read_csv(file.path(OUT_DIR, "cleaned_experiment_data.csv"), show_col_types = FALSE)
stats_by_benchmark_variant <- read_csv(file.path(OUT_DIR, "stats_by_benchmark_variant.csv"), show_col_types = FALSE)

safe_save <- function(plot, filename, width = 10, height = 6) {
  out_path <- file.path(PLOTS_DIR, filename)
  tryCatch({
    ggsave(out_path, plot = plot, width = width, height = height, dpi = 150)
  }, error = function(e) warning("Failed to save ", filename, ": ", conditionMessage(e)))
}

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
# plot: 03_energy_boxplot_by_variant.png
safe_save(p3, "03_energy_boxplot_by_variant.png", width = 12, height = 7)

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
# plot: 04_energy_boxplot_by_benchmark.png
safe_save(p4, "04_energy_boxplot_by_benchmark.png", width = 14, height = 7)

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
# plot: 05_energy_boxplot_faceted.png
safe_save(p5, "05_energy_boxplot_faceted.png", width = 16, height = 12)

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
# plot: 06_energy_violin_by_variant.png
safe_save(p6, "06_energy_violin_by_variant.png", width = 12, height = 7)

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
# plot: 10_energy_heatmap_mean.png
safe_save(p10, "10_energy_heatmap_mean.png", width = 12, height = 10)

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
# plot: 15_energy_comparison_by_type.png
safe_save(p15, "15_energy_comparison_by_type.png", width = 14, height = 7)

# =============================================================================
# =============================================================================

PLOTS_DIR <- file.path(WORKSPACE_ROOT, "analysis", "plots", "rq1")
dir.create(PLOTS_DIR, recursive = TRUE, showWarnings = FALSE)

data <- read_csv(file.path(OUT_DIR, "cleaned_experiment_data.csv"), show_col_types = FALSE)


rq1_data <- data %>%
  filter(variant %in% c("BASELINE", "OPTIMIZED"))


# =============================================================================
# =============================================================================

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


write_csv(delta_pct, file.path(OUT_DIR, "rq1_delta_pct_by_benchmark.csv"))

# =============================================================================
# =============================================================================

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

test_results <- delta_long %>%
  group_by(metric) %>%
  summarise(
    n = n(),
    median = median(delta_pct, na.rm = TRUE),
    mean = mean(delta_pct, na.rm = TRUE),
    sd = sd(delta_pct, na.rm = TRUE),
    shapiro_p = tryCatch(
      shapiro.test(delta_pct)$p.value,
      error = function(e) NA
    ),
    wilcox_p = tryCatch(
      wilcox.test(delta_pct, mu = 0, alternative = "two.sided")$p.value,
      error = function(e) NA
    ),
    ttest_p = tryCatch(
      t.test(delta_pct, mu = 0, alternative = "two.sided")$p.value,
      error = function(e) NA
    ),
    .groups = "drop"
  ) %>%
  mutate(
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

write_csv(test_results, file.path(OUT_DIR, "rq1_statistical_tests.csv"))

# =============================================================================
# =============================================================================

code_opt_guidelines <- data %>%
  filter(!variant %in% c("BASELINE", "OPTIMIZED", "CLEAN", "ORIGINAL_FIXED")) %>%
  pull(variant) %>%
  unique()


plot1a_data_raw <- data %>%
  filter(variant %in% code_opt_guidelines) %>%
  filter(!(benchmark == "basic_string" & variant == "G12")) %>%
  select(benchmark, variant, results_set, energy_joules, cpu_usage_avg, memory_usage_avg, power_watts)

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
    energy_mean = mean(delta_energy_pct, na.rm = TRUE),
    energy_se = sd(delta_energy_pct, na.rm = TRUE) / sqrt(n()),
    energy_ci_lower = energy_mean - 1.96 * energy_se,
    energy_ci_upper = energy_mean + 1.96 * energy_se,
    cpu_mean = mean(delta_cpu_pct, na.rm = TRUE),
    cpu_se = sd(delta_cpu_pct, na.rm = TRUE) / sqrt(n()),
    cpu_ci_lower = cpu_mean - 1.96 * cpu_se,
    cpu_ci_upper = cpu_mean + 1.96 * cpu_se,
    memory_mean = mean(delta_memory_pct, na.rm = TRUE),
    memory_se = sd(delta_memory_pct, na.rm = TRUE) / sqrt(n()),
    memory_ci_lower = memory_mean - 1.96 * memory_se,
    memory_ci_upper = memory_mean + 1.96 * memory_se,
    power_mean = mean(delta_power_pct, na.rm = TRUE),
    power_se = sd(delta_power_pct, na.rm = TRUE) / sqrt(n()),
    power_ci_lower = power_mean - 1.96 * power_se,
    power_ci_upper = power_mean + 1.96 * power_se,
    .groups = "drop"
  )

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

benchmark_order <- plot1a_long %>%
  filter(metric_clean == "Energy %") %>%
  arrange(delta_pct) %>%
  pull(benchmark)

plot1a_long <- plot1a_long %>%
  mutate(benchmark = factor(benchmark, levels = benchmark_order))

median_energy <- plot1a_long %>%
  filter(metric_clean == "Energy %") %>%
  summarise(median = median(delta_pct, na.rm = TRUE)) %>%
  pull(median)

y_limit <- 100  # Cap at ±100%

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

p1a <- ggplot(plot1a_data_labels, aes(x = benchmark, y = delta_pct_display, fill = delta_pct)) +
  geom_col(position = position_dodge(width = 0.8), alpha = 0.85, width = 0.75) +
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

# plot: plot1a_code_optimization_all_metrics_per_benchmark.png
ggsave(file.path(PLOTS_DIR, "plot1a_code_optimization_all_metrics_per_benchmark.png"), 
       plot = p1a, width = 14, height = 10, dpi = 300)
# plot: plot1a_code_optimization_all_metrics_per_benchmark.png

# =============================================================================
# =============================================================================

target_benchmarks <- c("sequential_minimum_optimization", "basic_string")

plot1_data_raw <- data %>%
  filter(benchmark %in% target_benchmarks,
         variant != "OPTIMIZED") %>%  # Exclude OPTIMIZED, keep BASELINE and guidelines
  select(benchmark, variant, results_set, energy_joules, exec_time_sec, 
         cpu_usage_avg, memory_usage_avg, power_watts)

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
    energy_mean = mean(delta_energy_pct, na.rm = TRUE),
    energy_se = sd(delta_energy_pct, na.rm = TRUE) / sqrt(n()),
    energy_ci_lower = energy_mean - 1.96 * energy_se,
    energy_ci_upper = energy_mean + 1.96 * energy_se,
    time_mean = mean(delta_time_pct, na.rm = TRUE),
    time_se = sd(delta_time_pct, na.rm = TRUE) / sqrt(n()),
    time_ci_lower = time_mean - 1.96 * time_se,
    time_ci_upper = time_mean + 1.96 * time_se,
    cpu_mean = mean(delta_cpu_pct, na.rm = TRUE),
    cpu_se = sd(delta_cpu_pct, na.rm = TRUE) / sqrt(n()),
    cpu_ci_lower = cpu_mean - 1.96 * cpu_se,
    cpu_ci_upper = cpu_mean + 1.96 * cpu_se,
    memory_mean = mean(delta_memory_pct, na.rm = TRUE),
    memory_se = sd(delta_memory_pct, na.rm = TRUE) / sqrt(n()),
    memory_ci_lower = memory_mean - 1.96 * memory_se,
    memory_ci_upper = memory_mean + 1.96 * memory_se,
    power_mean = mean(delta_power_pct, na.rm = TRUE),
    power_se = sd(delta_power_pct, na.rm = TRUE) / sqrt(n()),
    power_ci_lower = power_mean - 1.96 * power_se,
    power_ci_upper = power_mean + 1.96 * power_se,
    .groups = "drop"
  )

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

for (bench in target_benchmarks) {
  bench_clean <- ifelse(bench == "sequential_minimum_optimization", 
                        "Sequential Min. Opt.", "Basic String")
  
  plot_data <- plot1_long %>% filter(benchmark == bench)
  
  has_extreme <- any(abs(plot_data$delta_pct) > 1000, na.rm = TRUE)
  extreme_variants <- plot_data %>% 
    filter(abs(delta_pct) > 1000) %>% 
    pull(variant) %>% 
    unique()
  
  if (bench == "basic_string" && has_extreme) {
    plot_data_filtered <- plot_data %>% filter(abs(delta_pct) <= 1000)
    
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
  
  y_limit <- ifelse(bench == "basic_string", 100, 50)
  
  plot_data_with_labels <- plot_data_filtered %>%
    mutate(
      is_truncated = abs(delta_pct) > y_limit,
      display_value = ifelse(delta_pct > y_limit, y_limit - 5,
                            ifelse(delta_pct < -y_limit, -y_limit + 5, delta_pct)),
      label_text = ifelse(is_truncated, sprintf("%+.0f%%", round(delta_pct)), ""),
      label_y = ifelse(delta_pct > y_limit, y_limit - 2,
                      ifelse(delta_pct < -y_limit, -y_limit + 2, NA))
    )
  
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
  
  filename <- sprintf("plot1b_%s_guideline_impact.png", bench)
  ggsave(file.path(PLOTS_DIR, filename), 
         plot = p, width = 14, height = 8, dpi = 300)
}

# =============================================================================
# =============================================================================

slope_data <- rq1_data %>%
  filter(variant %in% c("BASELINE", "OPTIMIZED")) %>%
  group_by(benchmark, variant) %>%
  summarise(
    energy_mean = mean(energy_joules, na.rm = TRUE),
    energy_se = sd(energy_joules, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(variant = factor(variant, levels = c("BASELINE", "OPTIMIZED")))

benchmark_improvement <- slope_data %>%
  pivot_wider(names_from = variant, values_from = c(energy_mean, energy_se)) %>%
  mutate(
    improved = energy_mean_OPTIMIZED < energy_mean_BASELINE,
    improvement_label = ifelse(improved, "Improved (↓)", "Worsened (↑)")
  ) %>%
  select(benchmark, improved, improvement_label)

slope_data <- slope_data %>%
  left_join(benchmark_improvement, by = "benchmark") %>%
  mutate(improvement_label = factor(improvement_label, levels = c("Improved (↓)", "Worsened (↑)")))

overall_means <- slope_data %>%
  group_by(variant) %>%
  summarise(
    overall_mean = mean(energy_mean, na.rm = TRUE),
    overall_se = sd(energy_mean, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

median_delta_energy <- delta_pct %>%
  summarise(median_delta = median(delta_energy_pct, na.rm = TRUE)) %>%
  pull(median_delta)

benchmark_colors <- setNames(
  scales::hue_pal()(10),
  unique(slope_data$benchmark)
)

p2 <- ggplot(slope_data, aes(x = variant, y = energy_mean, group = benchmark, 
                              color = benchmark, linetype = improvement_label)) +
  geom_line(linewidth = 1.2, alpha = 0.8) +
  geom_point(size = 3, alpha = 0.9) +
  geom_line(data = overall_means, 
            aes(x = variant, y = overall_mean, group = 1),
            color = "black", linewidth = 2, inherit.aes = FALSE) +
  geom_point(data = overall_means,
             aes(x = variant, y = overall_mean),
             color = "black", size = 5, shape = 18, inherit.aes = FALSE) +
  geom_errorbar(data = overall_means,
                aes(x = variant, y = overall_mean, 
                    ymin = overall_mean - overall_se, 
                    ymax = overall_mean + overall_se),
                color = "black", width = 0.1, linewidth = 1.2,
                inherit.aes = FALSE) +
  annotate("text", x = 1.5, y = max(slope_data$energy_mean) * 0.95,
           label = sprintf("Median ΔEnergy: %+.1f%%\nOverall mean shown in black", median_delta_energy),
           size = 4, fontface = "bold", color = "black") +
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

# plot: plot2_paired_energy_comparison.png
ggsave(file.path(PLOTS_DIR, "plot2_paired_energy_comparison.png"), 
       plot = p2, width = 10, height = 8, dpi = 300)

# =============================================================================
# =============================================================================

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

# =============================================================================
# =============================================================================



sig_count <- sum(test_results$significant)
total_count <- nrow(test_results)

if (test_results %>% filter(metric == "Energy %") %>% pull(significant)) {
  energy_median <- test_results %>% filter(metric == "Energy %") %>% pull(median)
  direction <- ifelse(energy_median < 0, "REDUCTION", "INCREASE")
} else {
}

# plot: plot1a_code_optimization_all_metrics_per_benchmark.png
# plot: plot2_paired_energy_comparison.png



# =============================================================================
# =============================================================================

PLOTS_DIR <- file.path(WORKSPACE_ROOT, "analysis", "plots", "rq2")
dir.create(PLOTS_DIR, recursive = TRUE, showWarnings = FALSE)

data <- read_csv(file.path(OUT_DIR, "cleaned_experiment_data.csv"), show_col_types = FALSE)

multithreading_guidelines <- c("G9", "G12", "G14")


rq2_data <- data %>%
  filter(variant %in% c("BASELINE", multithreading_guidelines))


# =============================================================================
# =============================================================================

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
    energy_mean = mean(delta_energy_pct, na.rm = TRUE),
    energy_se = sd(delta_energy_pct, na.rm = TRUE) / sqrt(n()),
    energy_ci_lower = energy_mean - 1.96 * energy_se,
    energy_ci_upper = energy_mean + 1.96 * energy_se,
    time_mean = mean(delta_time_pct, na.rm = TRUE),
    time_se = sd(delta_time_pct, na.rm = TRUE) / sqrt(n()),
    time_ci_lower = time_mean - 1.96 * time_se,
    time_ci_upper = time_mean + 1.96 * time_se,
    cpu_mean = mean(delta_cpu_pct, na.rm = TRUE),
    cpu_se = sd(delta_cpu_pct, na.rm = TRUE) / sqrt(n()),
    cpu_ci_lower = cpu_mean - 1.96 * cpu_se,
    cpu_ci_upper = cpu_mean + 1.96 * cpu_se,
    memory_mean = mean(delta_memory_pct, na.rm = TRUE),
    memory_se = sd(delta_memory_pct, na.rm = TRUE) / sqrt(n()),
    memory_ci_lower = memory_mean - 1.96 * memory_se,
    memory_ci_upper = memory_mean + 1.96 * memory_se,
    power_mean = mean(delta_power_pct, na.rm = TRUE),
    power_se = sd(delta_power_pct, na.rm = TRUE) / sqrt(n()),
    power_ci_lower = power_mean - 1.96 * power_se,
    power_ci_upper = power_mean + 1.96 * power_se,
    .groups = "drop"
  )

write_csv(guideline_deltas, file.path(OUT_DIR, "rq2_multithreading_deltas.csv"))

# =============================================================================
# =============================================================================

plot1_data_raw <- rq2_data %>%
  filter(variant %in% multithreading_guidelines) %>%
  filter(!(benchmark == "basic_string" & variant == "G12")) %>%
  select(benchmark, variant, results_set, energy_joules, cpu_usage_avg, memory_usage_avg, power_watts)

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
    energy_mean = mean(delta_energy_pct, na.rm = TRUE),
    energy_se = sd(delta_energy_pct, na.rm = TRUE) / sqrt(n()),
    energy_ci_lower = energy_mean - 1.96 * energy_se,
    energy_ci_upper = energy_mean + 1.96 * energy_se,
    cpu_mean = mean(delta_cpu_pct, na.rm = TRUE),
    cpu_se = sd(delta_cpu_pct, na.rm = TRUE) / sqrt(n()),
    cpu_ci_lower = cpu_mean - 1.96 * cpu_se,
    cpu_ci_upper = cpu_mean + 1.96 * cpu_se,
    memory_mean = mean(delta_memory_pct, na.rm = TRUE),
    memory_se = sd(delta_memory_pct, na.rm = TRUE) / sqrt(n()),
    memory_ci_lower = memory_mean - 1.96 * memory_se,
    memory_ci_upper = memory_mean + 1.96 * memory_se,
    power_mean = mean(delta_power_pct, na.rm = TRUE),
    power_se = sd(delta_power_pct, na.rm = TRUE) / sqrt(n()),
    power_ci_lower = power_mean - 1.96 * power_se,
    power_ci_upper = power_mean + 1.96 * power_se,
    .groups = "drop"
  )

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

benchmark_order <- plot1_long %>%
  filter(metric_clean == "Energy %") %>%
  arrange(delta_pct) %>%
  pull(benchmark)

plot1_long <- plot1_long %>%
  mutate(benchmark = factor(benchmark, levels = benchmark_order))

median_energy <- plot1_long %>%
  filter(metric_clean == "Energy %") %>%
  summarise(median = median(delta_pct, na.rm = TRUE)) %>%
  pull(median)

y_limit <- 100  # Cap at ±100%

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

p1 <- ggplot(plot1_data_labels, aes(x = benchmark, y = delta_pct_display, fill = delta_pct)) +
  geom_col(position = position_dodge(width = 0.8), alpha = 0.85, width = 0.75) +
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

# plot: plot1_multithreading_all_metrics_per_benchmark.png
ggsave(file.path(PLOTS_DIR, "plot1_multithreading_all_metrics_per_benchmark.png"), 
       plot = p1, width = 14, height = 10, dpi = 300)

# =============================================================================
# =============================================================================

multithreading_avg <- rq2_data %>%
  filter(variant %in% multithreading_guidelines) %>%
  filter(!(benchmark == "basic_string" & variant == "G12")) %>%  # Exclude outlier
  group_by(benchmark, results_set) %>%
  summarise(
    energy_joules = mean(energy_joules, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(variant = "MULTITHREADING")

slope_data_raw <- bind_rows(
  rq2_data %>% filter(variant == "BASELINE") %>% select(benchmark, variant, results_set, energy_joules),
  multithreading_avg
)

slope_data <- slope_data_raw %>%
  group_by(benchmark, variant) %>%
  summarise(
    energy_mean = mean(energy_joules, na.rm = TRUE),
    energy_se = sd(energy_joules, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(variant = factor(variant, levels = c("BASELINE", "MULTITHREADING")))

benchmark_improvement <- slope_data %>%
  pivot_wider(names_from = variant, values_from = c(energy_mean, energy_se)) %>%
  mutate(
    improved = energy_mean_MULTITHREADING < energy_mean_BASELINE,
    improvement_label = ifelse(improved, "Improved (↓)", "Worsened (↑)")
  ) %>%
  select(benchmark, improved, improvement_label)

slope_data <- slope_data %>%
  left_join(benchmark_improvement, by = "benchmark") %>%
  mutate(improvement_label = factor(improvement_label, levels = c("Improved (↓)", "Worsened (↑)")))

overall_means <- slope_data %>%
  group_by(variant) %>%
  summarise(
    overall_mean = mean(energy_mean, na.rm = TRUE),
    overall_se = sd(energy_mean, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

baseline_avg <- slope_data %>% filter(variant == "BASELINE") %>% pull(energy_mean)
multithreading_avg_values <- slope_data %>% filter(variant == "MULTITHREADING") %>% pull(energy_mean)
delta_values <- ((multithreading_avg_values - baseline_avg) / baseline_avg) * 100
median_delta_energy <- median(delta_values, na.rm = TRUE)

benchmark_colors <- setNames(
  scales::hue_pal()(length(unique(slope_data$benchmark))),
  unique(slope_data$benchmark)
)

p2 <- ggplot(slope_data, aes(x = variant, y = energy_mean, group = benchmark, 
                              color = benchmark, linetype = improvement_label)) +
  geom_line(linewidth = 1.2, alpha = 0.8) +
  geom_point(size = 3, alpha = 0.9) +
  geom_line(data = overall_means, 
            aes(x = variant, y = overall_mean, group = 1),
            color = "black", linewidth = 2, inherit.aes = FALSE) +
  geom_point(data = overall_means,
             aes(x = variant, y = overall_mean),
             color = "black", size = 5, shape = 18, inherit.aes = FALSE) +
  geom_errorbar(data = overall_means,
                aes(x = variant, y = overall_mean, 
                    ymin = overall_mean - overall_se, 
                    ymax = overall_mean + overall_se),
                color = "black", width = 0.1, linewidth = 1.2,
                inherit.aes = FALSE) +
  annotate("text", x = 1.5, y = max(slope_data$energy_mean) * 0.95,
           label = sprintf("Median ΔEnergy: %+.1f%%\nOverall mean shown in black", median_delta_energy),
           size = 4, fontface = "bold", color = "black") +
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

# plot: plot2_paired_energy_comparison.png
ggsave(file.path(PLOTS_DIR, "plot2_paired_energy_comparison.png"), 
       plot = p2, width = 10, height = 8, dpi = 300)

# =============================================================================
# =============================================================================

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

# =============================================================================
# =============================================================================



# plot: plot1_multithreading_all_metrics_per_benchmark.png
# plot: plot2_paired_energy_comparison.png



# =============================================================================
# =============================================================================

PLOTS_DIR <- file.path(WORKSPACE_ROOT, "analysis", "plots", "rq3")
dir.create(PLOTS_DIR, recursive = TRUE, showWarnings = FALSE)

data <- read_csv(file.path(OUT_DIR, "cleaned_experiment_data.csv"), show_col_types = FALSE)


code_opt_guidelines <- c("G1", "G3", "G4", "G6", "G7")  # Exclude G12 (outlier), G9, G14 (multithreading)
multithreading_guidelines <- c("G9", "G12", "G14")


rq3_data <- data %>%
  filter(variant %in% c(code_opt_guidelines, multithreading_guidelines))


# =============================================================================
# =============================================================================

code_opt_avg <- rq3_data %>%
  filter(variant %in% code_opt_guidelines) %>%
  group_by(benchmark, results_set) %>%
  summarise(
    energy_joules = mean(energy_joules, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(treatment = "CODE_OPTIMIZATION")

multithreading_avg <- rq3_data %>%
  filter(variant %in% multithreading_guidelines) %>%
  filter(!(benchmark == "basic_string" & variant == "G12")) %>%  # Exclude outlier
  group_by(benchmark, results_set) %>%
  summarise(
    energy_joules = mean(energy_joules, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(treatment = "MULTITHREADING")

slope_data_raw <- bind_rows(code_opt_avg, multithreading_avg)

slope_data <- slope_data_raw %>%
  group_by(benchmark, treatment) %>%
  summarise(
    energy_mean = mean(energy_joules, na.rm = TRUE),
    energy_se = sd(energy_joules, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(treatment = factor(treatment, levels = c("CODE_OPTIMIZATION", "MULTITHREADING")))

benchmark_improvement <- slope_data %>%
  pivot_wider(names_from = treatment, values_from = c(energy_mean, energy_se)) %>%
  mutate(
    improved = energy_mean_MULTITHREADING < energy_mean_CODE_OPTIMIZATION,
    improvement_label = ifelse(improved, "Improved (↓)", "Worsened (↑)")
  ) %>%
  select(benchmark, improved, improvement_label)

slope_data <- slope_data %>%
  left_join(benchmark_improvement, by = "benchmark") %>%
  mutate(improvement_label = factor(improvement_label, levels = c("Improved (↓)", "Worsened (↑)")))

overall_means <- slope_data %>%
  group_by(treatment) %>%
  summarise(
    overall_mean = mean(energy_mean, na.rm = TRUE),
    overall_se = sd(energy_mean, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

code_opt_values <- slope_data %>% filter(treatment == "CODE_OPTIMIZATION") %>% pull(energy_mean)
multithreading_values <- slope_data %>% filter(treatment == "MULTITHREADING") %>% pull(energy_mean)
delta_values <- ((multithreading_values - code_opt_values) / code_opt_values) * 100
median_delta_energy <- median(delta_values, na.rm = TRUE)

benchmark_colors <- setNames(
  scales::hue_pal()(length(unique(slope_data$benchmark))),
  unique(slope_data$benchmark)
)

p <- ggplot(slope_data, aes(x = treatment, y = energy_mean, group = benchmark, 
                             color = benchmark, linetype = improvement_label)) +
  geom_line(linewidth = 1.2, alpha = 0.8) +
  geom_point(size = 3, alpha = 0.9) +
  geom_line(data = overall_means, 
            aes(x = treatment, y = overall_mean, group = 1),
            color = "black", linewidth = 2, inherit.aes = FALSE) +
  geom_point(data = overall_means,
             aes(x = treatment, y = overall_mean),
             color = "black", size = 5, shape = 18, inherit.aes = FALSE) +
  geom_errorbar(data = overall_means,
                aes(x = treatment, y = overall_mean, 
                    ymin = overall_mean - overall_se, 
                    ymax = overall_mean + overall_se),
                color = "black", width = 0.1, linewidth = 1.2,
                inherit.aes = FALSE) +
  annotate("text", x = 1.5, y = max(slope_data$energy_mean) * 0.95,
           label = sprintf("Median Δ: %+.1f%%\nOverall mean shown in black", median_delta_energy),
           size = 4, fontface = "bold", color = "black") +
  scale_y_log10(labels = label_number(scale_cut = cut_short_scale(), suffix = " J")) +
  scale_color_manual(values = benchmark_colors, name = "Benchmark") +
  scale_linetype_manual(
    values = c("Improved (↓)" = "dashed", "Worsened (↑)" = "solid"),
    name = "Effect",
    guide = guide_legend(order = 1)
  ) +
  scale_x_discrete(labels = c("CODE_OPTIMIZATION" = "Code Optimization", 
                              "MULTITHREADING" = "Multithreading")) +
  labs(
    title = "Per-Benchmark Energy Consumption: Code Optimization vs Multithreading",
    subtitle = "Each colored line connects paired runs from the same benchmark",
    x = NULL,
    y = "Energy Consumption (J, log scale)",
    caption = "Dashed lines = multithreading improved (↓ energy) vs code opt; Solid lines = worsened (↑).\nBlack line = overall mean ± SE. Upward slope indicates multithreading uses more energy.\nCode Opt = avg of G1, G3, G4, G6, G7; Multithreading = avg of G9, G12, G14 (basic_string+G12 excluded)."
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

# plot: plot_code_opt_vs_multithreading_comparison.png
ggsave(file.path(PLOTS_DIR, "plot_code_opt_vs_multithreading_comparison.png"), 
       plot = p, width = 10, height = 8, dpi = 300)

# =============================================================================
# =============================================================================

summary_stats <- slope_data %>%
  group_by(treatment) %>%
  summarise(
    n_benchmarks = n(),
    energy_median = median(energy_mean, na.rm = TRUE),
    energy_mean = mean(energy_mean, na.rm = TRUE),
    energy_sd = sd(energy_mean, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(summary_stats, file.path(OUT_DIR, "rq3_summary_statistics.csv"))

delta_by_benchmark <- slope_data %>%
  pivot_wider(names_from = treatment, values_from = c(energy_mean, energy_se)) %>%
  mutate(
    delta_energy_pct = ((energy_mean_MULTITHREADING - energy_mean_CODE_OPTIMIZATION) / 
                        energy_mean_CODE_OPTIMIZATION) * 100
  ) %>%
  select(benchmark, energy_mean_CODE_OPTIMIZATION, energy_mean_MULTITHREADING, delta_energy_pct)

write_csv(delta_by_benchmark, file.path(OUT_DIR, "rq3_delta_by_benchmark.csv"))

# =============================================================================
# =============================================================================



if (median_delta_energy > 0) {
} else {
}

n_improved <- sum(benchmark_improvement$improved)
n_worsened <- sum(!benchmark_improvement$improved)

# plot: plot_code_opt_vs_multithreading_comparison.png






# =============================================================================
# =============================================================================
