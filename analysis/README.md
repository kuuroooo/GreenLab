# Energy Benchmark Analysis - Complete Documentation

This directory contains comprehensive data exploration and statistical analysis of energy consumption experiments testing optimization guidelines across multiple Python benchmarks.

## 📁 Directory Structure

```
analysis/
├── data_exploration.R              # Main exploratory analysis script
├── create_summary_tables.R         # Publication-ready summary tables
├── analyze_energy.R                # Original analysis script
├── DATA_EXPLORATION_SUMMARY.md     # Detailed findings report
├── README.md                       # This file
├── outputs/                        # CSV tables and statistics
│   ├── [20+ summary CSV files]
└── plots/                          # Visualization files
    └── [15 visualization PNG files]
```

---

## 🚀 Quick Start

### Running the Analysis

1. **Complete Data Exploration** (recommended):
   ```bash
   Rscript analysis/data_exploration.R
   ```
   This generates all statistics, quality metrics, outlier detection, and visualizations.

2. **Create Summary Tables**:
   ```bash
   Rscript analysis/create_summary_tables.R
   ```
   This creates 10 publication-ready summary tables.

3. **View the Report**:
   Open `DATA_EXPLORATION_SUMMARY.md` for a comprehensive written report.

---

## 📊 Output Files Guide

### Key Statistical Summaries

| File | Description | Use Case |
|------|-------------|----------|
| **table1_overall_summary.csv** | Overall statistics across all data | Quick overview |
| **table2_summary_by_variant.csv** | Performance by variant (G1, G7, BASELINE, etc.) | Compare optimization guidelines |
| **table3_summary_by_benchmark.csv** | Performance by benchmark | Identify which programs vary most |
| **table10_guideline_effectiveness.csv** | Guideline performance vs baseline | **Key table for RQ answers** |
| **table4_best_worst_by_benchmark.csv** | Top/bottom performers per benchmark | Identify best optimization per program |

### Detailed Analysis Files

| File | Description |
|------|-------------|
| **stats_by_benchmark_variant.csv** | Full breakdown of every benchmark×variant combination |
| **stats_comparison_vs_baseline.csv** | % change calculations vs baseline for all combinations |
| **overall_descriptive_statistics.csv** | Comprehensive stats (mean, median, SD, variance, IQR, CV, skewness, kurtosis) |
| **data_quality_report.csv** | Outlier rates, stability metrics, CV by combination |
| **outliers_detailed.csv** | List of all detected outliers with context |
| **cleaned_experiment_data.csv** | Full dataset with calculated metrics (930 rows) |

### Summary Tables for Publications

| Table # | File | Purpose |
|---------|------|---------|
| 1 | table1_overall_summary.csv | Overall descriptive statistics |
| 2 | table2_summary_by_variant.csv | Variant comparison (sorted by efficiency) |
| 3 | table3_summary_by_benchmark.csv | Benchmark characteristics |
| 4 | table4_best_worst_by_benchmark.csv | Best/worst variant per benchmark |
| 5 | table5_energy_matrix_median.csv | Pivot: Benchmark × Variant (median energy) |
| 6 | table6_pct_change_matrix.csv | Pivot: % change vs baseline |
| 7 | table7_data_quality_matrix.csv | Outliers and stability assessment |
| 8 | table8_top10_efficient.csv | Most energy-efficient combinations |
| 9 | table9_top10_inefficient.csv | Least energy-efficient combinations |
| 10 | table10_guideline_effectiveness.csv | **Guideline summary for RQ** |

---

## 📈 Visualizations

### Distribution Plots (Understanding Energy Patterns)

| Plot | File | Purpose |
|------|------|---------|
| 1 | 01_energy_distribution_histogram.png | Overall energy histogram (log scale) |
| 2 | 02_energy_density_by_variant.png | Overlaid density curves per variant |
| 3 | 14_energy_ridgeplot_by_variant.png | Stacked density plots (ridge plot) |

### Comparison Plots (Variant vs Benchmark Performance)

| Plot | File | Purpose |
|------|------|---------|
| 4 | 03_energy_boxplot_by_variant.png | Energy by variant (all benchmarks combined) |
| 5 | 04_energy_boxplot_by_benchmark.png | Energy by benchmark (all variants combined) |
| 6 | 05_energy_boxplot_faceted.png | **Each benchmark's variants side-by-side** |
| 7 | 06_energy_violin_by_variant.png | Violin plots showing quartiles |
| 8 | 15_energy_comparison_by_type.png | Baseline vs Optimized vs Guidelines |

### Correlation & Relationship Plots

| Plot | File | Purpose |
|------|------|---------|
| 9 | 07_energy_vs_power_scatter.png | Relationship between energy and power |
| 10 | 08_energy_vs_cpu_scatter.png | CPU usage impact on energy |
| 11 | 09_energy_vs_memory_scatter.png | Memory usage impact on energy |

### Heatmaps (Matrix Visualizations)

| Plot | File | Purpose |
|------|------|---------|
| 12 | 10_energy_heatmap_mean.png | Mean energy: Benchmark × Variant |
| 13 | 11_energy_heatmap_pct_change.png | **% change vs baseline** (green=better, red=worse) |

### Data Quality Plots

| Plot | File | Purpose |
|------|------|---------|
| 14 | 12_stability_cv_by_variant.png | Coefficient of variation (measurement stability) |
| 15 | 13_outlier_rate_by_combination.png | Outlier frequency by combination |

---

## 🎯 How to Use These Results

### For Research Questions (RQ)

1. **Which guideline is most effective?**
   - See `table10_guideline_effectiveness.csv`
   - G4, G9, G3, G7 show ~6-9% median energy reduction vs baseline
   - G12 shows +294.6% increase (avoid!)

2. **Does effectiveness vary by benchmark?**
   - See `table4_best_worst_by_benchmark.csv`
   - Yes! Different guidelines work best for different programs
   - Example: G7 best for kmeans, G4 best for convolution_neural_network

3. **Are there interaction effects?**
   - See `11_energy_heatmap_pct_change.png` (visual)
   - See `table6_pct_change_matrix.csv` (numerical)
   - Clear evidence: guideline effectiveness is benchmark-dependent

### For Publication/Thesis

**Tables to Include:**

1. **Table 2** (summary_by_variant.csv) - Shows variant characteristics
2. **Table 4** (best_worst_by_benchmark.csv) - Shows benchmark-specific results
3. **Table 10** (guideline_effectiveness.csv) - Main RQ results
4. **Table 7** (data_quality_matrix.csv) - Supports methodology rigor

**Figures to Include:**

1. **Figure 1**: `05_energy_boxplot_faceted.png` - Shows all data compactly
2. **Figure 2**: `11_energy_heatmap_pct_change.png` - Shows % improvements/regressions
3. **Figure 3**: `12_stability_cv_by_variant.png` - Shows measurement reliability
4. **Figure 4**: `15_energy_comparison_by_type.png` - Shows overall trends

### For Statistical Testing (Next Steps)

The cleaned dataset (`cleaned_experiment_data.csv`) is ready for:

1. **ANOVA / Kruskal-Wallis Test**
   - Test if variant differences are statistically significant
   - Use energy_joules as dependent variable
   - Use variant as independent variable
   - Consider benchmark as blocking factor

2. **Post-hoc Tests**
   - Dunn's test with Bonferroni correction
   - Compare each guideline to baseline
   - Report adjusted p-values

3. **Effect Size Calculations**
   - Cohen's d for each guideline vs baseline
   - Helps determine practical significance

4. **Mixed-Effects Model**
   - Account for repeated measures (10 runs)
   - Model: energy ~ variant + (1|benchmark) + (1|results_set)

---

## 📊 Key Findings Summary

### Overall Statistics

- **930 measurements** from 10 benchmarks × 12 variants × 10 runs
- **Median energy: 0.98 J**, **Mean energy: 87.58 J** (highly skewed distribution)
- **13.98% outliers** detected (moderate: 1.5×IQR method)
- **Memory usage very stable** (CV = 1.92%) across all runs

### Best Performers

| Rank | Variant | Median Energy | vs Baseline |
|------|---------|---------------|-------------|
| 1 | CLEAN | 0.46 J | - |
| 2 | G4 | 0.94 J | -8.8% |
| 3 | G9 | 0.95 J | -7.7% |
| 4 | G3 | 0.96 J | -7.0% |
| 5 | G7 | 0.97 J | -6.3% |

### Worst Performers

- **OPTIMIZED**: 452.53 J mean (needs investigation!)
- **G12**: 210.95 J mean (+294.6% vs baseline)
- **G14**: 39.44 J mean (+32.0% vs baseline)

### Data Quality

- **No missing data** ✓
- **High stability**: CLEAN variant (CV = 14.37%)
- **High variability**: 28 combinations with CV > 20%
- Most combinations show **moderate to high measurement variability**

---

## 🔧 Technical Details

### Statistical Methods

1. **Descriptive Statistics**: Mean, median, SD, variance, min, max, Q25, Q75, IQR, CV
2. **Outlier Detection**: 
   - IQR method (1.5×IQR for moderate, 3×IQR for extreme)
   - Z-score method (|z| > 2.5 flagged)
3. **Normalization**: 
   - Min-max scaling per benchmark (0-1 scale)
   - Z-score standardization for cross-benchmark comparison
4. **Stability Metrics**: Coefficient of Variation (CV = SD/mean × 100%)

### Data Sources

- **RESULTS folder**: run_table.csv files (aggregated metrics)
- **EXPERIMENTS folder**: Individual CSV files (time-series data)
- **Log files**: energibridge.log files (backup energy values)

### Software

- **R version**: 4.5+
- **Packages**: ggplot2, dplyr, tidyr, readr, purrr, stringr, scales, forcats, ggridges

---

## 📝 Interpretation Guidelines

### Coefficient of Variation (CV)

- **< 10%**: Stable measurements
- **10-20%**: Moderate variability
- **20-30%**: High variability
- **> 30%**: Very high variability (consider more runs)

### % Change vs Baseline

- **Negative values**: Energy reduction (good! 🎉)
- **Positive values**: Energy increase (bad 😞)
- **Near 0%**: No meaningful change
- **|% change| < 5%**: Might not be practically significant

### Statistical vs Practical Significance

- **Statistical significance** (p < 0.05): Difference is real, not due to chance
- **Practical significance** (effect size): Difference is large enough to matter
- Always consider **both** when interpreting results

---

## 🐛 Known Issues & Warnings

1. **Log-scale warnings**: Some plots use log scale, which produces warnings for zero values (normal)
2. **Ridge plot error**: One ridge plot failed due to insufficient data in some groups (non-critical)
3. **Skewed distribution**: Use median instead of mean for comparisons
4. **High variability**: Some combinations need more runs for stable estimates

---

## 📚 Additional Resources

- **DATA_EXPLORATION_SUMMARY.md**: Detailed written report with interpretation
- **Original analysis script**: `analyze_energy.R` (legacy, kept for reference)
- **Raw data**: See `RESULTS/` and `EXPERIMENTS/` folders

---

## 💡 Tips for Further Analysis

1. **Focus on median values** - more robust to outliers than mean
2. **Use IQR instead of SD** when reporting variability (for skewed data)
3. **Check data quality table** before making conclusions
4. **Consider benchmark-specific recommendations** - no one-size-fits-all
5. **Report effect sizes** along with p-values

---

## ✅ Checklist for Paper/Thesis

- [ ] Read `DATA_EXPLORATION_SUMMARY.md`
- [ ] Review all summary tables (especially Tables 2, 4, 10)
- [ ] Examine key plots (especially plots 5, 11, 12, 15)
- [ ] Check outliers in `outliers_detailed.csv`
- [ ] Review data quality in `data_quality_report.csv`
- [ ] Run statistical tests (ANOVA, post-hoc)
- [ ] Calculate effect sizes
- [ ] Verify assumptions (normality, homoscedasticity)
- [ ] Prepare figures with proper captions
- [ ] Write methods section describing analysis
- [ ] Interpret results in context of RQ

---

## 🤝 Questions or Issues?

If you need:
- Different visualizations
- Additional statistical tests
- Custom summary tables
- Help interpreting results

Modify the R scripts or create new ones based on the existing templates!

---

**Last Updated:** October 16, 2025  
**Analysis Version:** 1.1  
**Data:** 910 measurements, 10 benchmarks, 10 variants, 10 runs  
**Note:** CLEAN and ORIGINAL_FIXED variants excluded from analysis

