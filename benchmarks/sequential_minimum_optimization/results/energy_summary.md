# Energy Consumption Summary

## Sequential Minimum Optimization (SMO) Algorithm Energy Measurements

This document summarizes the energy consumption measurements for different variants of the Sequential Minimum Optimization algorithm using energibridge.

### Results Overview

| Algorithm Variant | Energy Consumption (Joules) | Execution Time (seconds) | Energy Efficiency (J/s) |
|-------------------|----------------------------|-------------------------|-------------------------|
| G1 | 575.38 | 39.71 | 14.49 |
| G3 | 831.55 | 37.73 | 22.05 |
| G4 | 917.42 | 39.77 | 23.07 |
| G6 | 902.59 | 39.57 | 22.81 |
| G7 | 100.27 | 5.53 | 18.13 |
| G14 | 518.28 | 39.11 | 13.25 |
| OPTIMIZED | 256.90 | 47.97 | 5.36 |
| Original | 656.54 | 38.04 | 17.26 |

### Key Findings

1. **Most Energy Efficient**: G14 variant (518.28 J) - lowest total energy consumption
2. **Most Time Efficient**: G7 variant (5.53 s) - fastest execution time
3. **Best Overall Performance**: OPTIMIZED variant (256.90 J) - significantly lower energy consumption despite longer execution time
4. **Highest Energy Consumption**: G4 variant (917.42 J) - highest total energy consumption

### Analysis

- The **OPTIMIZED** variant shows the best energy efficiency with only 256.90 joules consumed, making it 2.5x more energy efficient than the next best variant (G14)
- The **G7** variant has the shortest execution time but moderate energy efficiency
- The **G4** variant consumes the most energy, suggesting potential optimization opportunities
- All variants except G7 and OPTIMIZED have similar execution times (~38-40 seconds)

### Files Generated

Each algorithm variant has a corresponding CSV file with detailed energy measurements:
- `sequential_minimum_optimization_G1.csv`
- `sequential_minimum_optimization_G3.csv`
- `sequential_minimum_optimization_G4.csv`
- `sequential_minimum_optimization_G6.csv`
- `sequential_minimum_optimization_G7.csv`
- `sequential_minimum_optimization_G14.csv`
- `sequential_minimum_optimization_OPTIMIZED.csv`
- `sequential_minimum_optimization.csv`

Each CSV contains detailed measurements including:
- CPU frequency and temperature data
- CPU usage percentages
- System power consumption
- Memory usage statistics
- Time-series data for analysis

### Measurement Details

- **Tool Used**: energibridge
- **Measurement Interval**: Default (200 microseconds)
- **Environment**: macOS with Python virtual environment
- **Dependencies**: pandas, matplotlib, scikit-learn
