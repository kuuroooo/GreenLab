# RQ1 Analysis Summary

**Research Question 1:** Does code optimization significantly impact energy consumption and related metrics compared to baseline?

**Generated:** October 16, 2025  
**Analysis:** OPTIMIZED variant vs BASELINE across 10 benchmarks

---

## 📊 Key Findings

### Overall Result
**None of the metrics showed statistically significant deviation from zero at α = 0.05**

However, the **median trends suggest**:
- ✅ **Energy**: -3.94% (median reduction, not significant p=1.000)
- ✅ **Power**: -3.89% (median reduction, not significant p=0.770)
- ⚠️ **CPU**: +2.90% (median increase, not significant p=0.193)
- ⚠️ **Memory**: +0.17% (median negligible, not significant p=0.193)
- ❌ **Time**: No data available (execution time = 0 in all measurements)

---

## 📈 Statistical Test Results

All metrics tested using **Wilcoxon signed-rank test** (due to non-normal distributions as indicated by Shapiro-Wilk test p < 0.05).

| Metric | Median Δ% | Mean Δ% | SD | p-value | Significant? |
|--------|-----------|---------|-----|---------|--------------|
| **Energy %** | -3.94% | +36,034% | 110,981% | 1.000 | No (ns) |
| **Power %** | -3.89% | +215% | 574% | 0.770 | No (ns) |
| **CPU %** | +2.90% | +51% | 113% | 0.193 | No (ns) |
| **Memory %** | +0.17% | +0.42% | 0.79% | 0.193 | No (ns) |
| **Time %** | N/A | N/A | N/A | N/A | N/A |

### Key Observations

1. **High Variability**: The massive difference between median and mean values (especially for Energy: median -3.94% vs mean +36,034%!) indicates **extreme outliers** in some benchmarks.

2. **Non-significance**: Despite median improvements, the high variability across benchmarks means we cannot confidently reject the null hypothesis that the optimization has no effect.

3. **Time Data Missing**: Execution time was not properly captured (all values = 0), preventing time-based analysis.

---

## 🔍 Per-Benchmark Analysis

| Benchmark | Energy Δ% | Power Δ% | CPU Δ% | Memory Δ% |
|-----------|-----------|----------|---------|-----------|
| ant_colony_optimization | +8.1% | +6.8% | -1.6% | +0.1% |
| **basic_string** | **+351,796%** 🚨 | **+1,800%** 🚨 | +360% | +1.2% |
| convolution_neural_network | -12.0% ✅ | -12.0% ✅ | +1.3% | +0.3% |
| kmeans | +30.4% | +1.5% | +40.3% | +0.7% |
| linear_discriminant_analysis | -1.2% ✅ | -1.5% ✅ | -14.2% ✅ | +0.2% |
| n_body_simulation | -8.9% ✅ | -8.3% ✅ | -1.0% ✅ | -0.1% ✅ |
| **sequential_minimum_optimization** | **-96.0%** ✅✅✅ | **-45.0%** ✅ | -3.3% ✅ | -0.3% ✅ |
| **sol1** | **+8,668%** 🚨 | **+426%** 🚨 | +31.6% | +2.3% |
| sol1_problem180 | -33.5% ✅ | -12.5% ✅ | +87.5% | -0.3% |
| two_hidden_layers_neural_network | -6.6% ✅ | -6.3% ✅ | +4.5% | +0.1% |

### Critical Issues Identified

**🚨 Severe Outliers:**
1. **basic_string**: +351,796% energy increase! This is catastrophic.
2. **sol1**: +8,668% energy increase - also extremely problematic.

**✅ Best Performers:**
1. **sequential_minimum_optimization**: -96.0% energy reduction! Outstanding.
2. **sol1_problem180**: -33.5% energy reduction
3. **convolution_neural_network**: -12.0% energy reduction
4. **n_body_simulation**: -8.9% energy reduction

---

## 📊 Visualizations

### Plot 1: Overall Distribution of Δ% Changes

**File:** `plots/rq1/plot1_delta_pct_distribution.png`

This multi-variable boxplot shows:
- Distribution of Δ% across all benchmarks for each metric
- Median line with 95% CI notches
- Individual benchmark data points
- Reference line at y = 0 (no change)
- Statistical significance annotations

**Interpretation:**
- Most boxplots span across y=0, indicating high variability
- Energy and Power medians are slightly below zero (small improvement trend)
- CPU and Memory medians are near zero (minimal change)
- Wide IQRs confirm heterogeneous effects across benchmarks

### Plot 2: Paired Per-Benchmark Energy Comparison

**File:** `plots/rq1/plot2_paired_energy_comparison.png`

This slope chart shows:
- Each gray line connects Baseline → Optimized for one benchmark
- Red bold line shows overall mean ± SE
- Median ΔEnergy annotated on plot
- Log scale for Y-axis (energy in Joules)

**Interpretation:**
- Most lines slope downward → energy reduction trend
- Two lines shoot upward dramatically (basic_string, sol1) → major outliers
- Overall mean trend is slightly downward
- High variability visible

---

## 💡 Interpretation & Discussion

### Why No Statistical Significance?

Despite seeing some benchmarks with large energy reductions (e.g., -96% for sequential_minimum_optimization), the statistical tests find no significant overall effect because:

1. **Extreme Outliers**: Two benchmarks (basic_string, sol1) show catastrophic energy increases (+351,796% and +8,668%), completely dominating the statistics.

2. **High Heterogeneity**: The optimization effect is highly **benchmark-dependent**:
   - Some benchmarks benefit greatly (−96%)
   - Some are harmed greatly (+351,796%)
   - This variability prevents finding a consistent overall effect

3. **Small Sample Size**: Only 10 benchmarks means low statistical power to detect effects amid high variance.

### What Does This Mean?

**The "OPTIMIZED" variant is NOT universally beneficial:**
- ✅ Works well for: sequential_minimum_optimization, convolution_neural_network, n_body_simulation
- ❌ Fails catastrophically for: basic_string, sol1
- ⚠️ Mixed results for others

**This suggests the "OPTIMIZED" version may have bugs or inappropriate optimizations for certain code patterns.**

---

## 🎯 Recommendations

### For the Thesis/Paper

1. **Report Non-Significance**: Be transparent that the OPTIMIZED variant does not show significant overall improvement.

2. **Investigate Outliers**: The basic_string and sol1 results need investigation:
   - Is there a bug in the OPTIMIZED implementations?
   - Are these legitimate worst-case scenarios?
   - Should these be excluded as implementation errors?

3. **Benchmark-Specific Analysis**: Instead of asking "Does optimization work?", ask "For which types of algorithms does optimization work?"

4. **Alternative Approach**: Compare individual **guidelines** (G1-G14) rather than the combined OPTIMIZED variant, as shown in the exploratory analysis where several guidelines showed consistent improvements.

### Statistical Recommendations

1. **Remove Outliers**: Consider re-running analysis excluding basic_string and sol1 if you can justify them as implementation errors.

2. **Non-Parametric Tests**: Already using Wilcoxon (good choice given non-normality).

3. **Effect Sizes**: Report median Δ% alongside p-values to show practical significance even when statistical significance is absent.

4. **Stratified Analysis**: Group benchmarks by type (string manipulation, ML algorithms, numerical simulation) and analyze separately.

---

## 📁 Output Files

### Plots (High Resolution - 300 DPI)
- `plots/rq1/plot1_delta_pct_distribution.png` - Multi-metric distribution
- `plots/rq1/plot2_paired_energy_comparison.png` - Paired slope chart

### Data Tables
- `outputs/rq1_delta_pct_by_benchmark.csv` - Δ% calculations per benchmark
- `outputs/rq1_statistical_tests.csv` - Test results with p-values
- `outputs/rq1_summary_statistics.csv` - Descriptive statistics

### R Script
- `rq1_analysis.R` - Reproducible analysis script

---

## 🔄 Next Steps

1. **Verify OPTIMIZED Implementations**: Check basic_string and sol1 for bugs
2. **Consider Alternative RQ**: Instead of "Does OPTIMIZED work?", focus on "Which individual guidelines work?"
3. **Run RQ2 Analysis**: Compare individual guidelines (G1-G14) vs baseline
4. **Interaction Analysis**: Examine benchmark × guideline interactions

---

## 📝 Citation Template

For your paper/thesis:

> "We compared the OPTIMIZED variant against BASELINE across 10 benchmarks. Although median energy consumption showed a -3.94% reduction trend, this difference was not statistically significant (Wilcoxon signed-rank test, p = 1.000, n=10). The lack of significance was driven by extreme heterogeneity in optimization effects: while some benchmarks showed substantial reductions (e.g., sequential_minimum_optimization: -96.0%), others showed catastrophic increases (basic_string: +351,796%). This suggests that optimization effectiveness is highly benchmark-dependent and that the combined OPTIMIZED variant contains implementation issues requiring further investigation."

---

**Summary:** The OPTIMIZED variant does not show statistically significant improvement over BASELINE, primarily due to catastrophic failures on 2 benchmarks (basic_string, sol1) that overwhelm the positive effects seen on other benchmarks. Further investigation of these outliers is needed, and alternative analysis focusing on individual optimization guidelines may be more informative.

