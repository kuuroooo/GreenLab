#!/bin/bash
EXPERIMENT_NUM=10
SRC_DIR="$(pwd)/benchmarks/Results/"
SRC_DIR_2="$(pwd)/benchmarks/experiments/"

echo "[*] Setting Macbook into headless mode..."
ssh macbook "nohup ./headless.sh start > /dev/null 2>&1 &"


echo "[*] Removing old experiment results..."
echo "Deleting original files from $SRC_DIR ..."
rm -rf "$SRC_DIR"/* 
rm -rf "$SRC_DIR_2"/* 

echo "[*] Running warmup configuration ..."
python WarmupConfig.py  #Run the experiment warmup
sleep 30

for ((i=1; i<=EXPERIMENT_NUM; i++)); do
    echo "[*] Starting experiment runner ..."
    echo "Running experiment..."
    python experiment-runner/experiment-runner/__main__.py benchmarks/RunnerConfig.py

    echo "Cleaning and migrating experiment results..."
    
    EXP_DIR="experiment_$i"
    EXP_DIR_2="results_$i"

    mkdir -p "$EXP_DIR"
    mkdir -p "$EXP_DIR_2"
    echo "Copying results to $EXP_DIR and $EXP_DIR_2 ..."
    cp -r "$SRC_DIR"/* "$EXP_DIR"/
    cp -r "$SRC_DIR_2"/* "$EXP_DIR_2"/

    echo "Deleting original files from $SRC_DIR ..."
    rm -rf "$SRC_DIR"/* 
    rm -rf "$SRC_DIR_2"/* 
    echo "Done! Archived to: $EXP_DIR"
    echo "Done! Archived to: $EXP_DIR_2"
    echo "Sleeping 30 seconds to cool down"
    sleep 30

done


echo "[*] Waking Macbook and stopping headless mode"
ssh macbook "nohup ./headless.sh stop > /dev/null 2>&1 &;"
ssh macbook "./headless.sh wake > /dev/null 2>&1 &"
echo "[*] Script completed "
