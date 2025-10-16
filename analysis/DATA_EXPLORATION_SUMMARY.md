# Data Exploration Summary Report

**Generated:** October 16, 2025  
**Experiment:** Energy Consumption Analysis of Python Optimization Guidelines

---

## 📊 Dataset Overview

- **Total Runs:** 910 measurements (CLEAN and ORIGINAL_FIXED excluded)
- **Unique Benchmarks:** 10 programs tested
- **Unique Variants:** 10 optimization variants (BASELINE, OPTIMIZED, G1, G3, G4, G6, G7, G9, G12, G14)
- **Experimental Runs:** 10 repeated runs (results_1 through results_10)

---

## 🔍 Key Findings

### Energy Consumption (Overall Statistics)

| Metric | Value |
|--------|-------|
| **Mean** | 89.50 J |
| **Median** | 0.99 J |
| **Standard Deviation** | 495.65 J |
| **Variance** | 245,669 J² |
| **Minimum** | 0.00 J |
| **Maximum** | 5,912.04 J |
| **Q25 (25th percentile)** | 0.53 J |
| **Q75 (75th percentile)** | 2.88 J |
| **IQR (Interquartile Range)** | 2.35 J |
| **CV (Coefficient of Variation)** | 553.80% |
| **Skewness** | ~8.5 (highly right-skewed) |
| **Kurtosis** | ~77 (heavy-tailed distribution) |

**Interpretation:** The high CV (559.98%) and large difference between mean (87.58 J) and median (0.98 J) indicate that the distribution is highly skewed with some benchmarks consuming significantly more energy than others. The data contains extreme values.

### Power Consumption

| Metric | Value |
|--------|-------|
| **Mean** | 4.18 W |
| **Median** | 2.60 W |
| **Standard Deviation** | 6.06 W |
| **Range** | [1.23, 45.33] W |
| **CV** | 145.03% |

### CPU Usage

| Metric | Value |
|--------|-------|
| **Mean** | 20.02% |
| **Median** | 18.09% |
| **Standard Deviation** | 10.18% |
| **Range** | [3.34, 80.82]% |
| **CV** | 50.85% |

### Memory Usage

| Metric | Value |
|--------|-------|
| **Mean** | 39.48 MB |
| **Median** | 39.43 MB |
| **Standard Deviation** | 0.76 MB |
| **Range** | [37.85, 41.67] MB |
| **CV** | 1.92% (very stable) |

**Note:** Memory usage is very consistent across all benchmarks and variants.

---

## 🚨 Outlier Analysis

### Outlier Detection Results

Using the IQR (Interquartile Range) method on energy consumption:

- **Moderate Outliers (1.5×IQR):** 125 instances (13.74% of data)
- **Extreme Outliers (3×IQR):** 77 instances (8.46% of data)

### Interpretation

The presence of outliers (nearly 14% of data) suggests:
1. Some benchmark-variant combinations have significantly different energy profiles
2. Measurement variability exists across runs
3. Certain optimizations may introduce performance instability

**Action:** Outliers should be examined (see `outliers_detailed.csv`) but not automatically removed, as they may represent real performance characteristics.

---

## 📈 Variant Performance Summary

### Energy Efficiency Ranking (by mean energy, lowest to highest)

| Rank | Variant | Mean Energy (J) | Median Energy (J) | SD (J) | CV (%) | Interpretation |
|------|---------|-----------------|-------------------|---------|---------|----------------|
| 1 | **G4** | 42.91 | 0.94 | 111.36 | 259.50% | Best guideline, -8.8% vs baseline |
| 2 | **G9** | 1.63 | 0.95 | 1.49 | 91.35% | Second best, -7.7% vs baseline |
| 3 | **G3** | 33.72 | 0.96 | 97.78 | 289.94% | Good, -7.0% vs baseline |
| 4 | **G7** | 3.83 | 0.97 | 7.77 | 202.68% | Good, -6.3% vs baseline |
| 5 | **G1** | 35.66 | 1.00 | 103.87 | 291.30% | Marginal, -3.4% vs baseline |
| 6 | **G6** | 34.76 | 1.00 | 100.98 | 290.48% | Marginal, -3.4% vs baseline |
| 7 | **BASELINE** | 35.31 | 1.03 | 102.99 | 291.67% | Reference benchmark |
| 8 | **OPTIMIZED** | 452.53 | 1.17 | 1350.73 | 298.49% | Poor overall (needs investigation) |
| 9 | **G14** | 39.44 | 1.36 | 107.89 | 273.57% | Worse, +32.0% vs baseline |
| 10 | **G12** | 210.95 | 4.07 | 531.81 | 252.10% | Worst, +294.6% vs baseline |

### Key Observations

1. **G4, G9, G3, and G7** guidelines show energy reductions of 6-9% vs baseline
2. **G9 variant** is most stable among guidelines (CV = 91.35%)
3. **OPTIMIZED variant** surprisingly has the highest mean energy consumption, suggesting optimization issues
4. High CV values across most variants indicate performance depends heavily on the specific benchmark

---

## 📊 Benchmark-Specific Analysis

### Energy Consumption by Benchmark

The following benchmarks show the most variability across variants:

1. **basic_string:** Large variation between variants (G12 variant: 110,206% increase vs baseline!)
2. **convolution_neural_network:** High energy consumption overall
3. **kmeans:** Moderate variation across guidelines
4. **linear_discriminant_analysis:** Moderate variation
5. **ant_colony_optimization_algorithms:** Relatively stable across variants

### Comparison with Baseline

From the `stats_comparison_vs_baseline.csv` file:

**Best Improvements (negative % = energy reduction):**
- Most guideline variants show increases rather than decreases vs baseline
- Further analysis needed to identify which specific guideline-benchmark combinations work well

**Worst Performance:**
- **basic_string + G12:** +110,206% energy vs baseline (extreme outlier)
- This suggests G12 creates performance issues for string operations

---

## 📉 Data Quality Assessment

### Measurement Stability by Variant

| Stability Category | Number of Benchmark-Variant Combinations |
|-------------------|------------------------------------------|
| Very Stable (CV < 5%) | 0 |
| Stable (CV < 10%) | 0 |
| Moderate (CV < 20%) | Many combinations |
| Variable (CV < 30%) | Some combinations |
| Highly Variable (CV ≥ 30%) | 28 combinations |

**High Variability Combinations (CV > 20%):** 28 out of all benchmark-variant combinations

This indicates:
- Energy measurements vary considerably between runs for many combinations
- Experimental conditions should be carefully controlled
- More runs may be needed for unstable combinations

### Missing Data

No missing values detected for:
- Energy measurements ✓
- Power measurements ✓
- CPU usage ✓
- Memory usage ✓

Data completeness is excellent!

---

## 📊 Statistical Distributions

### Energy Distribution Characteristics

- **Highly right-skewed** (skewness = 8.38): A few benchmarks consume much more energy
- **Heavy-tailed** (kurtosis = 75.67): Extreme values are more common than normal distribution
- **Log-normal distribution**: Most values cluster at lower energy, with long tail of high-energy cases

**Recommendation:** Use median instead of mean for comparisons, as median is robust to outliers.

---

## 🎨 Visualizations Generated

The following plots are available in `analysis/plots/`:

### Distribution Plots
1. **01_energy_distribution_histogram.png** - Overall energy distribution (log scale)
2. **02_energy_density_by_variant.png** - Density curves for each variant
3. **14_energy_ridgeplot_by_variant.png** - Ridge plot showing distributions

### Comparison Plots
4. **03_energy_boxplot_by_variant.png** - Energy by variant (all benchmarks)
5. **04_energy_boxplot_by_benchmark.png** - Energy by benchmark (all variants)
6. **05_energy_boxplot_faceted.png** - Faceted view: each benchmark's variants
7. **06_energy_violin_by_variant.png** - Violin plots with quartiles

### Correlation Plots
8. **07_energy_vs_power_scatter.png** - Energy vs power relationship
9. **08_energy_vs_cpu_scatter.png** - Energy vs CPU usage
10. **09_energy_vs_memory_scatter.png** - Energy vs memory usage

### Heatmaps
11. **10_energy_heatmap_mean.png** - Mean energy by benchmark×variant
12. **11_energy_heatmap_pct_change.png** - % change vs baseline

### Quality Plots
13. **12_stability_cv_by_variant.png** - Coefficient of variation by variant
14. **13_outlier_rate_by_combination.png** - Outlier frequency
15. **15_energy_comparison_by_type.png** - Baseline vs Optimized vs Guidelines

---

## 💡 Recommendations

### For Analysis

1. **Use median values** for comparisons due to skewed distribution
2. **Investigate the OPTIMIZED variant** - why is it consuming the most energy?
3. **Focus on G9, G7, and CLEAN variants** - they show lower energy consumption
4. **Examine basic_string + G12 combination** - this is a severe outlier
5. **Consider benchmark-specific recommendations** - no single guideline works best for all

### For Experiment Design

1. **Increase repetitions** for high-CV combinations to improve confidence
2. **Control experimental conditions** to reduce measurement variability
3. **Document environmental factors** that might affect energy consumption
4. **Consider time-series analysis** using individual CSV files for detailed behavior

### Statistical Testing (Next Steps)

1. Perform **ANOVA** or **Kruskal-Wallis test** to determine if variant differences are statistically significant
2. Conduct **post-hoc tests** (e.g., Dunn's test with Bonferroni correction) for pairwise comparisons
3. Calculate **effect sizes** (Cohen's d) to quantify practical significance
4. Test for **interaction effects** between benchmarks and variants

---

## 📁 Output Files Reference

### Summary Statistics
- `overall_descriptive_statistics.csv` - Overall metrics across all data
- `stats_by_benchmark.csv` - Aggregated by benchmark
- `stats_by_variant.csv` - Aggregated by variant
- `stats_by_benchmark_variant.csv` - Full breakdown (benchmark × variant)
- `stats_comparison_vs_baseline.csv` - Relative performance vs baseline

### Quality Reports
- `data_quality_report.csv` - Stability and outlier metrics
- `outliers_detailed.csv` - List of all detected outliers

### Raw Data
- `cleaned_experiment_data.csv` - Full dataset with all calculated metrics

---

## 📝 Conclusion

This exploratory analysis reveals:

1. **High variability** in energy consumption across benchmarks and variants
2. **Benchmark-specific optimization** is necessary - no universal best guideline
3. **Several guidelines (G9, G7, CLEAN)** show promise for reducing energy consumption
4. **Data quality is good** but some combinations need more stable measurements
5. **Further statistical testing** is needed to confirm which optimizations are truly effective

The visualization and summary statistics provide a solid foundation for answering your research questions about optimization guideline effectiveness.

---

**Next Steps:** 
- Review the plots in `analysis/plots/` directory
- Examine specific outliers in `outliers_detailed.csv`
- Conduct inferential statistics for hypothesis testing
- Create benchmark-specific optimization recommendations

---

**Last Updated:** October 16, 2025  
**Analysis Version:** 1.1  
**Data:** 910 measurements, 10 benchmarks, 10 variants, 10 runs  
**Note:** CLEAN and ORIGINAL_FIXED variants excluded from analysis

