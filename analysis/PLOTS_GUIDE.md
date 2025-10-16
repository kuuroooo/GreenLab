# Energy Analysis Plots Guide

This document describes all the plots generated from the energy analysis and what insights they provide.

## Quick Start

```bash
# Generate all plots
cd /Users/kellywang/Desktop/GreenLab
Rscript analysis/analyze_energy.R
Rscript analysis/additional_plots.R
```

All plots are saved in: `analysis/plots/`

---

## Basic Analysis Plots

### 1. Per-Benchmark Energy Bars
**Files:** `energy_by_variant_<benchmark>.png` (10 files)
- **What it shows:** Mean energy consumption with error bars for each variant within a single benchmark
- **Use for:** Comparing guideline performance within a specific benchmark
- **Key insight:** Identify which optimization guidelines work best for each algorithm

### 2. Heatmap: % Change vs Baseline
**File:** `heatmap_pct_vs_baseline.png`
- **What it shows:** Color-coded matrix of % energy change (rows=benchmarks, cols=variants)
- **Use for:** Quick visual overview of which guidelines help/hurt which benchmarks
- **Key insight:** Green = improvement, Red = regression, patterns show guideline consistency

### 3. Boxplot: % Change Distribution
**File:** `boxplot_pct_vs_baseline_per_run.png`
- **What it shows:** Distribution of % changes across all runs and benchmarks per variant
- **Use for:** Understanding overall guideline effectiveness and variance
- **Key insight:** Medians show typical impact, boxes show consistency

### 4. Energy vs Average Power Scatter
**File:** `scatter_energy_vs_avg_power.png`
- **What it shows:** Relationship between total energy and average power draw
- **Use for:** Identifying power-efficient vs energy-intensive variants
- **Key insight:** Points cluster by variant; outliers indicate anomalies

### 5. Benchmark Lines Across Variants
**File:** `lines_pct_vs_baseline_by_benchmark.png`
- **What it shows:** Each benchmark as a line showing % change trajectory across variants
- **Use for:** Comparing how different benchmarks respond to the same guidelines
- **Key insight:** Parallel lines = similar behavior; diverging = benchmark-specific effects

---

## Advanced Analysis Plots

### 6. Guideline Effectiveness Average
**File:** `guideline_effectiveness_avg.png`
- **What it shows:** Horizontal bar chart of average % change per guideline (averaged across all benchmarks)
- **Use for:** Ranking guidelines from most beneficial to most harmful overall
- **Key insight:** Identifies universally good/bad optimization strategies

### 7. Best vs Worst Guidelines per Benchmark
**File:** `best_worst_guidelines_per_benchmark.png`
- **What it shows:** Side-by-side comparison of best and worst performing guideline for each benchmark
- **Use for:** Understanding optimization range and benchmark-specific sensitivities
- **Key insight:** Large gaps = high optimization potential

### 8. Variability by Variant
**File:** `variability_by_variant.png`
- **What it shows:** Coefficient of variation (CV) showing measurement consistency
- **Use for:** Assessing reliability of measurements across 10 experimental runs
- **Key insight:** Low CV = stable/reproducible, High CV = noisy measurements

### 9. Energy Change Distribution
**File:** `energy_change_distribution.png`
- **What it shows:** Histogram of all % changes across all benchmark-variant combinations
- **Use for:** Understanding overall distribution (normal? bimodal? skewed?)
- **Key insight:** Shape reveals whether optimizations are typically helpful or harmful

### 10. Pareto Analysis: Guideline Impact
**File:** `pareto_guideline_impact.png`
- **What it shows:** Bars = absolute impact, line = cumulative contribution
- **Use for:** Identifying which guidelines contribute most to total energy changes (80/20 rule)
- **Key insight:** Focus optimization efforts on high-impact guidelines

### 11. CPU vs Memory with Energy Overlay
**File:** `cpu_memory_energy_correlation.png`
- **What it shows:** Scatter plot of CPU usage vs memory usage, colored/sized by energy
- **Use for:** Identifying resource usage patterns and their relationship to energy
- **Key insight:** Clusters reveal computational profiles; hot colors = high energy

### 12. Optimization Potential by Benchmark
**File:** `optimization_potential_by_benchmark.png`
- **What it shows:** Energy range across all variants as % of baseline energy
- **Use for:** Identifying which benchmarks have most room for optimization
- **Key insight:** Large ranges = sensitive to guidelines, small ranges = stable

### 13. Top Performer Frequency
**File:** `top_performer_frequency.png`
- **What it shows:** Count of how often each guideline appears in top 3 performers per benchmark
- **Use for:** Finding consistently effective guidelines
- **Key insight:** High frequency = reliable optimization strategy

### 14. Power Efficiency by Variant
**File:** `power_efficiency_by_variant.png`
- **What it shows:** Violin + box plots of average power (Watts) distribution
- **Use for:** Comparing power draw patterns across guidelines
- **Key insight:** Lower power = more efficient execution

### 15. Faceted % Change: All Benchmarks
**File:** `faceted_pct_change_all_benchmarks.png`
- **What it shows:** Small multiples showing all variants for all benchmarks in one view
- **Use for:** Comprehensive overview of entire experiment at a glance
- **Key insight:** Patterns across facets reveal global vs local effects

### 16. OPTIMIZED vs BASELINE Comparison
**File:** `optimized_vs_baseline_comparison.png`
- **What it shows:** % change when ALL guidelines are applied together
- **Use for:** Assessing cumulative impact of combining all optimization strategies
- **Key insight:** Reveals whether guidelines are additive, synergistic, or conflicting

### 17. Top 20 Energy Changes
**File:** `top_20_energy_changes.png`
- **What it shows:** Largest absolute energy differences (both improvements and regressions)
- **Use for:** Identifying most dramatic optimization successes and failures
- **Key insight:** Extreme cases warrant deeper investigation

---

## How to Interpret the Plots

### Color Schemes
- **Green/Blue:** Energy reduction (good)
- **Red/Orange:** Energy increase (bad)
- **Grey:** Neutral or baseline

### Statistical Notes
- **Error bars:** Standard deviation across 10 runs
- **Box plots:** Box = IQR (25th-75th percentile), whiskers = 1.5×IQR
- **Violin plots:** Density distribution shape

### Key Metrics
- **% vs Baseline:** Positive = increase, Negative = decrease
- **Energy (J):** Total joules consumed during execution
- **Power (W):** Average watts = Energy / Time
- **CV (%):** Lower = more consistent measurements

---

## Recommended Analysis Workflow

1. **Start with:** `guideline_effectiveness_avg.png` – Get overall picture
2. **Then check:** `optimized_vs_baseline_comparison.png` – See combined effect
3. **Explore:** `heatmap_pct_vs_baseline.png` – Find patterns
4. **Deep dive:** Individual benchmark plots for interesting cases
5. **Validate:** `variability_by_variant.png` – Check data quality
6. **Investigate:** `top_20_energy_changes.png` – Understand extremes

---

## Questions These Plots Answer

- Which optimization guideline is most effective overall? → Plot #6
- Which benchmark benefits most from optimization? → Plot #12
- Are the measurements reliable? → Plot #8
- Do guidelines combine well or conflict? → Plot #16
- Which guideline is most consistent across benchmarks? → Plot #13
- What's the typical energy impact distribution? → Plot #9
- How do CPU and memory relate to energy? → Plot #11
- Which specific cases show the biggest changes? → Plot #17

---

## Next Steps

After reviewing the plots:
1. Identify best-performing guidelines for your use case
2. Check for benchmark-specific recommendations
3. Investigate any unexpected results (outliers, regressions)
4. Consider interaction effects between guidelines
5. Validate findings with additional runs if CV is high
6. Document insights for reporting or publication

For questions or issues, refer to the main README.md in the analysis directory.

