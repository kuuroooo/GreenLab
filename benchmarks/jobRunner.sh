#!/bin/bash

echo "[*] Setting Macbook into headless mode..."
ssh macbook "nohup ./headless.sh start > /dev/null 2>&1 &"

# Run experiment 10 times with 30 seconds between each run
for i in {1..10}; do
    echo "[*] Starting experiment runner - Run $i of 10..."

    echo "Cleaning experiment results..."
    rm -rf benchmarks/experiments/benchmarks_energy_analysis/

    echo "Running experiment with RUN_NUMBER=$i..."
    export RUN_NUMBER=$i
    python experiment-runner/experiment-runner/_main_.py benchmarks/RunnerConfig.py
    
    # Wait 30 seconds between runs (skip wait after last run)
    if [ $i -lt 10 ]; then
        echo "[*] Waiting 30 seconds before next run..."
        sleep 30
    fi
done

echo "[*] Waking Macbook and stopping headless mode"
ssh macbook "nohup ./headless.sh stop > /dev/null 2>&1 &;"
ssh macbook "./headless.sh wake > /dev/null 2>&1 &"
echo "[*] Script completed - All 10 runs finished!"
