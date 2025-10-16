#!/usr/bin/env Rscript
# =============================================================================
# RQ3 Analysis: Code Optimization vs Multithreading
# =============================================================================
# Research Question 3: How do code optimization guidelines compare to 
# multithreading guidelines in terms of energy consumption?

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
PLOTS_DIR <- file.path(WORKSPACE_ROOT, "analysis", "plots", "rq3")
dir.create(PLOTS_DIR, recursive = TRUE, showWarnings = FALSE)

message("=============================================================================")
message("RQ3 Analysis: Code Optimization vs Multithreading")
message("=============================================================================")
message("Output directory: ", PLOTS_DIR)
message("")

# =============================================================================
# LOAD DATA
# =============================================================================
data <- read_csv(file.path(OUT_DIR, "cleaned_experiment_data.csv"), show_col_types = FALSE)

message("Data loaded: ", nrow(data), " measurements")

# Define guideline groups
code_opt_guidelines <- c("G1", "G3", "G4", "G6", "G7")  # Exclude G12 (outlier), G9, G14 (multithreading)
multithreading_guidelines <- c("G9", "G12", "G14")

message("Code optimization guidelines: ", paste(code_opt_guidelines, collapse = ", "))
message("Multithreading guidelines: ", paste(multithreading_guidelines, collapse = ", "))

# Filter data
rq3_data <- data %>%
  filter(variant %in% c(code_opt_guidelines, multithreading_guidelines))

message("RQ3 data: ", nrow(rq3_data), " measurements")
message("Benchmarks: ", n_distinct(rq3_data$benchmark))

# =============================================================================
# PAIRED ENERGY COMPARISON: Code Optimization vs Multithreading
# =============================================================================
message("\nCreating paired energy comparison plot...")

# Calculate average code optimization effect for each benchmark (excluding basic_string + G12)
code_opt_avg <- rq3_data %>%
  filter(variant %in% code_opt_guidelines) %>%
  group_by(benchmark, results_set) %>%
  summarise(
    energy_joules = mean(energy_joules, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(treatment = "CODE_OPTIMIZATION")

# Calculate average multithreading effect for each benchmark (excluding basic_string + G12)
multithreading_avg <- rq3_data %>%
  filter(variant %in% multithreading_guidelines) %>%
  filter(!(benchmark == "basic_string" & variant == "G12")) %>%  # Exclude outlier
  group_by(benchmark, results_set) %>%
  summarise(
    energy_joules = mean(energy_joules, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(treatment = "MULTITHREADING")

# Combine both treatments
slope_data_raw <- bind_rows(code_opt_avg, multithreading_avg)

# Calculate means for slope chart
slope_data <- slope_data_raw %>%
  group_by(benchmark, treatment) %>%
  summarise(
    energy_mean = mean(energy_joules, na.rm = TRUE),
    energy_se = sd(energy_joules, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(treatment = factor(treatment, levels = c("CODE_OPTIMIZATION", "MULTITHREADING")))

# Determine if multithreading improved or worsened compared to code optimization
benchmark_improvement <- slope_data %>%
  pivot_wider(names_from = treatment, values_from = c(energy_mean, energy_se)) %>%
  mutate(
    improved = energy_mean_MULTITHREADING < energy_mean_CODE_OPTIMIZATION,
    improvement_label = ifelse(improved, "Improved (↓)", "Worsened (↑)")
  ) %>%
  select(benchmark, improved, improvement_label)

# Add improvement info back to slope_data
slope_data <- slope_data %>%
  left_join(benchmark_improvement, by = "benchmark") %>%
  mutate(improvement_label = factor(improvement_label, levels = c("Improved (↓)", "Worsened (↑)")))

# Calculate overall means
overall_means <- slope_data %>%
  group_by(treatment) %>%
  summarise(
    overall_mean = mean(energy_mean, na.rm = TRUE),
    overall_se = sd(energy_mean, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# Calculate median delta
code_opt_values <- slope_data %>% filter(treatment == "CODE_OPTIMIZATION") %>% pull(energy_mean)
multithreading_values <- slope_data %>% filter(treatment == "MULTITHREADING") %>% pull(energy_mean)
delta_values <- ((multithreading_values - code_opt_values) / code_opt_values) * 100
median_delta_energy <- median(delta_values, na.rm = TRUE)

# Define colors for benchmarks
benchmark_colors <- setNames(
  scales::hue_pal()(length(unique(slope_data$benchmark))),
  unique(slope_data$benchmark)
)

# Create the slope plot
p <- ggplot(slope_data, aes(x = treatment, y = energy_mean, group = benchmark, 
                             color = benchmark, linetype = improvement_label)) +
  # Individual benchmark lines
  geom_line(linewidth = 1.2, alpha = 0.8) +
  geom_point(size = 3, alpha = 0.9) +
  # Overall mean lines (bold black)
  geom_line(data = overall_means, 
            aes(x = treatment, y = overall_mean, group = 1),
            color = "black", linewidth = 2, inherit.aes = FALSE) +
  geom_point(data = overall_means,
             aes(x = treatment, y = overall_mean),
             color = "black", size = 5, shape = 18, inherit.aes = FALSE) +
  # Error bars for overall means
  geom_errorbar(data = overall_means,
                aes(x = treatment, y = overall_mean, 
                    ymin = overall_mean - overall_se, 
                    ymax = overall_mean + overall_se),
                color = "black", width = 0.1, linewidth = 1.2,
                inherit.aes = FALSE) +
  # Annotation for median delta
  annotate("text", x = 1.5, y = max(slope_data$energy_mean) * 0.95,
           label = sprintf("Median Δ: %+.1f%%\nOverall mean shown in black", median_delta_energy),
           size = 4, fontface = "bold", color = "black") +
  # Styling
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

# Save plot
ggsave(file.path(PLOTS_DIR, "plot_code_opt_vs_multithreading_comparison.png"), 
       plot = p, width = 10, height = 8, dpi = 300)
message("Saved: plot_code_opt_vs_multithreading_comparison.png")

# =============================================================================
# SUMMARY STATISTICS
# =============================================================================
message("\nCalculating summary statistics...")

# Calculate statistics by treatment
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
message("Saved: rq3_summary_statistics.csv")

# Save delta values
delta_by_benchmark <- slope_data %>%
  pivot_wider(names_from = treatment, values_from = c(energy_mean, energy_se)) %>%
  mutate(
    delta_energy_pct = ((energy_mean_MULTITHREADING - energy_mean_CODE_OPTIMIZATION) / 
                        energy_mean_CODE_OPTIMIZATION) * 100
  ) %>%
  select(benchmark, energy_mean_CODE_OPTIMIZATION, energy_mean_MULTITHREADING, delta_energy_pct)

write_csv(delta_by_benchmark, file.path(OUT_DIR, "rq3_delta_by_benchmark.csv"))
message("Saved: rq3_delta_by_benchmark.csv")

# =============================================================================
# FINAL REPORT
# =============================================================================
message("\n=============================================================================")
message("RQ3 ANALYSIS SUMMARY")
message("=============================================================================")

message("\nResearch Question 3:")
message("How do code optimization guidelines compare to multithreading guidelines")
message("in terms of energy consumption?")

message("\n--- SUMMARY STATISTICS ---")
for (i in 1:nrow(summary_stats)) {
  row <- summary_stats[i, ]
  message(sprintf("\n%s:", row$treatment))
  message(sprintf("  Benchmarks analyzed: %d", row$n_benchmarks))
  message(sprintf("  Median energy: %.2f J", row$energy_median))
  message(sprintf("  Mean energy: %.2f J (±%.2f)", row$energy_mean, row$energy_sd))
}

message(sprintf("\n--- OVERALL COMPARISON ---"))
message(sprintf("Median Δ (Multithreading vs Code Opt): %+.1f%%", median_delta_energy))
if (median_delta_energy > 0) {
  message("→ Multithreading uses MORE energy than code optimization on average")
} else {
  message("→ Multithreading uses LESS energy than code optimization on average")
}

message("\n--- PER-BENCHMARK RESULTS ---")
n_improved <- sum(benchmark_improvement$improved)
n_worsened <- sum(!benchmark_improvement$improved)
message(sprintf("Benchmarks where multithreading improved (↓): %d", n_improved))
message(sprintf("Benchmarks where multithreading worsened (↑): %d", n_worsened))

message("\n--- OUTPUT FILES ---")
message("Plots:")
message("  - ", file.path(PLOTS_DIR, "plot_code_opt_vs_multithreading_comparison.png"))
message("\nData tables:")
message("  - ", file.path(OUT_DIR, "rq3_summary_statistics.csv"))
message("  - ", file.path(OUT_DIR, "rq3_delta_by_benchmark.csv"))

message("\n=============================================================================")
message("RQ3 Analysis Complete!")
message("=============================================================================\n")

