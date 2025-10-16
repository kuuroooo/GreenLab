# Analysis Changelog

## Version 1.1 - October 16, 2025

### Changes Made

**Excluded Variants:**
- ✅ Removed **CLEAN** variant from all analyses
- ✅ Removed **ORIGINAL_FIXED** variant from all analyses

**Rationale:** Focus on the core experiment comparing BASELINE, OPTIMIZED, and optimization guidelines (G1, G3, G4, G6, G7, G9, G12, G14).

### Updated Statistics

| Metric | Before (v1.0) | After (v1.1) | Change |
|--------|---------------|--------------|--------|
| **Total Measurements** | 930 | 910 | -20 |
| **Unique Variants** | 12 | 10 | -2 |
| **Median Energy (J)** | 0.98 | 0.99 | +0.01 |
| **Mean Energy (J)** | 87.58 | 89.50 | +1.92 |
| **Outliers (moderate)** | 130 (13.98%) | 125 (13.74%) | -5 |
| **Outliers (extreme)** | 81 (8.71%) | 77 (8.46%) | -4 |

### Impact on Results

**Guideline Effectiveness (vs Baseline):**
- **G4**: -8.8% ✅ (Best)
- **G9**: -7.7% ✅ (Second best)
- **G3**: -7.0% ✅ (Good)
- **G7**: -6.3% ✅ (Good)
- **G1**: -3.4% ⚠️ (Marginal)
- **G6**: -3.4% ⚠️ (Marginal)
- **G14**: +32.0% ❌ (Worse)
- **G12**: +294.6% ❌ (Much worse)

**Most Stable Variant:** G9 (CV = 91.35%)

### Files Updated

**R Scripts:**
- ✅ `data_exploration.R` - Added filtering after variant parsing
- ✅ `create_summary_tables.R` - Added filtering after data load

**Documentation:**
- ✅ `QUICK_START_GUIDE.md` - Updated all statistics
- ✅ `DATA_EXPLORATION_SUMMARY.md` - Updated dataset info and variant rankings
- ✅ `README.md` - Updated version info

**Data Files (all regenerated):**
- ✅ All 18 CSV files in `outputs/`
- ✅ All 15 PNG files in `plots/`

### How to Verify

```bash
# Check the filtered data
wc -l analysis/outputs/cleaned_experiment_data.csv
# Should show 911 lines (910 data + 1 header)

# Check unique variants
cut -d',' -f'variant' analysis/outputs/cleaned_experiment_data.csv | sort -u
# Should show 10 variants (no CLEAN, no ORIGINAL_FIXED)

# Verify summary tables
head analysis/outputs/table2_summary_by_variant.csv
# Should show only 10 variants
```

### Scripts Are Idempotent

Running the scripts multiple times will produce the same results. The filtering happens early in the pipeline, so all downstream analyses automatically exclude CLEAN and ORIGINAL_FIXED.

---

## Version 1.0 - October 16, 2025

### Initial Release

**Features:**
- Complete data exploration of 930 measurements
- Analysis of 12 variants (including CLEAN and ORIGINAL_FIXED)
- 18 summary CSV files
- 15 visualization plots
- Comprehensive documentation

---

**How to Roll Back (if needed):**

If you need to include CLEAN and ORIGINAL_FIXED again, simply comment out or remove these lines:

In `data_exploration.R` (lines ~193-200):
```r
# # Filter out CLEAN and ORIGINAL_FIXED variants
# message("Filtering out CLEAN and ORIGINAL_FIXED variants...")
# n_before <- nrow(run_tables)
# run_tables <- run_tables %>%
#   filter(!variant %in% c("CLEAN", "ORIGINAL_FIXED"))
# n_after <- nrow(run_tables)
# message("Removed ", n_before - n_after, " rows")
```

In `create_summary_tables.R` (lines ~17-24):
```r
# # Filter out CLEAN and ORIGINAL_FIXED variants
# message("Filtering out CLEAN and ORIGINAL_FIXED variants...")
# n_before <- nrow(data)
# data <- data %>%
#   filter(!variant %in% c("CLEAN", "ORIGINAL_FIXED"))
# n_after <- nrow(data)
# message("Removed ", n_before - n_after, " rows")
```

Then re-run both scripts.

