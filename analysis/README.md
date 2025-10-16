## Energy analysis in R

This folder contains an R script to ingest all `RESULTS/results_*/benchmarks_energy_analysis/run_table.csv` files, join with per-run `energibridge.log` as a fallback, compute summaries, and produce plots.

### Requirements

- R (>= 4.1)
- Packages: `tidyverse`, `readr`, `dplyr`, `stringr`, `purrr`, `tidyr`, `scales`

Install packages:

```r
install.packages(c("tidyverse", "readr", "dplyr", "stringr", "purrr", "tidyr", "scales"))
```

### Run

```bash
Rscript /Users/kellywang/Desktop/GreenLab/analysis/analyze_energy.R
```

Optional environment overrides (use absolute paths):

```bash
RESULTS_ROOT=/Users/kellywang/Desktop/GreenLab/RESULTS \
EXPERIMENTS_ROOT=/Users/kellywang/Desktop/GreenLab/EXPERIMENTS \
Rscript /Users/kellywang/Desktop/GreenLab/analysis/analyze_energy.R
```

### Outputs

**CSV Files:**
- `combined_runs_raw.csv` – All runs before cleaning (original data)
- `combined_runs_cleaned.csv` – Cleaned data with normalized metrics and quality flags
- `summary_by_benchmark_variant.csv` – Aggregated statistics with % changes
- `outliers_detected.csv` – Flagged outliers (>3×IQR method)
- `data_quality_summary.csv` – Per-variant quality metrics (CV, outlier counts)

**New Columns in Cleaned Data:**
- `is_outlier` – Boolean flag for extreme outliers
- `energy_z_score` – Z-score normalized energy per benchmark
- `energy_normalized` – Min-max normalized (0-1) per benchmark
- `pct_vs_baseline_run` – % change vs baseline for this run
- `efficiency_ratio` – Energy / Baseline (< 1 = improvement)
- `calculated_power` – Derived from energy/time if missing
- `power_final` – Best available power value

### Plots Generated

Run the main analysis:
```bash
Rscript /Users/kellywang/Desktop/GreenLab/analysis/analyze_energy.R
```

Run additional plots:
```bash
Rscript /Users/kellywang/Desktop/GreenLab/analysis/additional_plots.R
```

Run data quality plots:
```bash
Rscript /Users/kellywang/Desktop/GreenLab/analysis/data_quality_plots.R
```

**Basic Plots** (`analyze_energy.R`):
- `energy_by_variant_<benchmark>.png` – Per-benchmark energy bars by variant (10 plots)
- `heatmap_pct_vs_baseline.png` – Heatmap of % change vs baseline
- `boxplot_pct_vs_baseline_per_run.png` – Per-run % change distribution
- `scatter_energy_vs_avg_power.png` – Energy vs power scatter
- `lines_pct_vs_baseline_by_benchmark.png` – Benchmark lines across variants

**Advanced Plots** (`additional_plots.R`):
1. `guideline_effectiveness_avg.png` – Average guideline impact across benchmarks
2. `best_worst_guidelines_per_benchmark.png` – Best/worst performing guidelines per benchmark
3. `variability_by_variant.png` – Measurement consistency (CV) by variant
4. `energy_change_distribution.png` – Histogram of all energy changes
5. `pareto_guideline_impact.png` – Pareto analysis of guideline contributions
6. `cpu_memory_energy_correlation.png` – CPU vs memory usage with energy overlay
7. `optimization_potential_by_benchmark.png` – Energy range as % of baseline
8. `top_performer_frequency.png` – Most common guidelines in top performers
9. `power_efficiency_by_variant.png` – Power (Watts) distribution by variant
10. `faceted_pct_change_all_benchmarks.png` – All benchmarks in faceted view
11. `optimized_vs_baseline_comparison.png` – Impact of applying all guidelines
12. `top_20_energy_changes.png` – Largest absolute energy differences

**Data Quality Plots** (`data_quality_plots.R`):
1. `quality_outliers_by_variant.png` – Outlier counts per benchmark-variant
2. `quality_cv_heatmap.png` – Coefficient of variation heatmap (consistency)
3. `quality_energy_distribution.png` – Energy histogram with outliers flagged
4. `quality_zscore_distribution.png` – Z-score distribution (normality check)
5. `quality_normalized_comparison.png` – BASELINE vs OPTIMIZED on 0-1 scale
6. `quality_efficiency_ratio.png` – Efficiency ratio distributions
7. `quality_benchmark_stability.png` – Average CV per benchmark
8. `quality_outlier_impact.png` – Impact of outliers on mean values


