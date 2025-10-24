
from typing import Dict, List, Any, Optional
from pathlib import Path
from os.path import dirname, realpath
import os
import signal
import pandas as pd
import time
import subprocess
import shlex
import glob
from os.path import expanduser
import threading
from dotenv import load_dotenv
import os

load_dotenv()

LAPTOP_SSH = os.getenv("SSH_COMMAND")
REMOTE_REPO = os.getenv("REMOTE_REPO")

ROOT_DIR = Path(dirname(realpath(__file__)))
BENCHMARKS_DIR = ROOT_DIR / 'benchmarks'

def stream_output(proc):
    for line in proc.stdout:
        print(f"[REMOTE] {line}", end="")


def discover_python_files() -> List[Path]:
        """Discover all Python files in the benchmarks directory"""
        python_files = []
        
        # Find all .py files recursively in benchmarks directory
        for py_file in BENCHMARKS_DIR.rglob("*.py"):
            # Skip __pycache__ and other system files
            if "__pycache__" not in str(py_file) and "RunnerConfig.py" not in str(py_file):
                python_files.append(py_file)
        
        # Sort files for consistent ordering
        python_files.sort()
        return python_files



def start_warmup() -> None:
    ROOT_DIR = Path(dirname(realpath(__file__)))
    BENCHMARKS_DIR = ROOT_DIR
    
    # Discover all Python files in benchmarks directory
    python_files = discover_python_files()
    #print(f"{python_files}")

    for pyFile in python_files:

        rel_path = pyFile.relative_to(BENCHMARKS_DIR)
        #print(f'Warming up running: {rel_path}')
        
        remote_script = f"{REMOTE_REPO}/{rel_path}"

        #remote_cmd = f"source /Documents/compSci/p1/greenLab/GreenLab/venv/bin/activate"

        #print(REMOTE_REPO)
        remote_cmd = (
            f"cd / &&"
            f"cd {REMOTE_REPO} && "
            f"source venv/bin/activate && "
            f'export PATH="$HOME/.cargo/bin:$PATH" && '
            f"python3 {shlex.quote(remote_script)}"
        )


        print(f"[Warming_Script] Running script: {pyFile.name}")
        #print(f"[Executing command ] SSH command: {remote_cmd}")

        # start_measurement
        profiler = subprocess.Popen(
            [
                "ssh",
                "macbook",
                "bash", "-c", remote_cmd
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            bufsize=1,
            text=True,
        )
        #t = threading.Thread(target=stream_output, args=(profiler,))
        #t.start()

        profiler.wait() #Wait for the script to finish running
        #t.join()
        print(f"[Warming_Script] Completed script: {pyFile.name}")
    

if __name__ == "__main__":
    start_warmup()