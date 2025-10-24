from EventManager.Models.RunnerEvents import RunnerEvents
from EventManager.EventSubscriptionController import EventSubscriptionController
from ConfigValidator.Config.Models.RunTableModel import RunTableModel
from ConfigValidator.Config.Models.FactorModel import FactorModel
from ConfigValidator.Config.Models.RunnerContext import RunnerContext
from ConfigValidator.Config.Models.OperationType import OperationType
from ProgressManager.Output.OutputProcedure import OutputProcedure as output

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

LAPTOP_SSH = "ssh -i ~/.ssh/kwang kellywang@192.168.178.96"
REMOTE_REPO = "/Users/kellywang/Documents/compSci/p1/greenLab/GreenLab"

class RunnerConfig:
    ROOT_DIR = Path(dirname(realpath(__file__)))
    BENCHMARKS_DIR = ROOT_DIR

    # ================================ USER SPECIFIC CONFIG ================================
    """The name of the experiment."""
    name: str = "benchmarks_energy_analysis"

    """The path in which Experiment Runner will create a folder with the name `self.name`, in order to store the
    results from this experiment. (Path does not need to exist - it will be created if necessary.)
    Output path defaults to the config file's path, inside the folder 'experiments'"""
    results_output_path: Path = ROOT_DIR / 'experiments'

    """Experiment operation type. Unless you manually want to initiate each run, use `OperationType.AUTO`."""
    operation_type: OperationType = OperationType.AUTO

    """The time Experiment Runner will wait after a run completes.
    This can be essential to accommodate for cooldown periods on some systems."""
    time_between_runs_in_ms: int = 2000

    def __init__(self):
        """Executes immediately after program start, on config load"""
        
        # Discover all Python files in benchmarks directory
        self.python_files = self._discover_python_files()
        
        EventSubscriptionController.subscribe_to_multiple_events([
            (RunnerEvents.BEFORE_EXPERIMENT, self.before_experiment),
            (RunnerEvents.BEFORE_RUN, self.before_run),
            (RunnerEvents.START_RUN, self.start_run),
            (RunnerEvents.START_MEASUREMENT, self.start_measurement),
            (RunnerEvents.INTERACT, self.interact),
            (RunnerEvents.STOP_MEASUREMENT, self.stop_measurement),
            (RunnerEvents.STOP_RUN, self.stop_run),
            (RunnerEvents.POPULATE_RUN_DATA, self.populate_run_data),
            (RunnerEvents.AFTER_EXPERIMENT, self.after_experiment)
        ])
        self.run_table_model = None  # Initialized later
        output.console_log(f"Custom config loaded. Found {len(self.python_files)} Python files to benchmark.")

    def _discover_python_files(self) -> List[Path]:
        """Discover all Python files in the benchmarks directory"""
        python_files = []
        
        # Find all .py files recursively in benchmarks directory
        for py_file in self.BENCHMARKS_DIR.rglob("*.py"):
            # Skip __pycache__ and other system files
            if "__pycache__" not in str(py_file) and "RunnerConfig.py" not in str(py_file):
                python_files.append(py_file)
        
        # Sort files for consistent ordering
        python_files.sort()
        return python_files

    def create_run_table_model(self) -> RunTableModel:
        """Create and return the run_table model here. A run_table is a List (rows) of tuples (columns),
        representing each run performed"""
        
        # Create a single factor for the Python file to run
        file_paths = [str(py_file.relative_to(self.BENCHMARKS_DIR)) for py_file in self.python_files]
        file_factor = FactorModel("python_file", file_paths)
        
        self.run_table_model = RunTableModel(
            factors=[file_factor],
            data_columns=['total_energy_joules', 'execution_time_seconds', 'avg_power_watts', 
                          'max_power_watts', 'min_power_watts', 'cpu_usage_avg', 'memory_usage_avg']
        )
        return self.run_table_model

    def before_experiment(self) -> None:
        """Perform any activity required before starting the experiment here
        Invoked only once during the lifetime of the program."""
        output.console_log(f"Starting energy analysis for {len(self.python_files)} benchmark files")
        
        # Create results directories for each experiment type
        self._create_results_directories()

    def _create_results_directories(self):
        """Create local Results/<experiment_name> directory on the Pi"""
        results_dir = self.ROOT_DIR / "Results" / self.name
        results_dir.mkdir(parents=True, exist_ok=True)
        output.console_log(f"Results will be stored locally in: {results_dir}")

    def before_run(self) -> None:
        """Perform any activity required before starting a run.
        No context is available here as the run is not yet active (BEFORE RUN)"""
        pass

    def start_run(self, context: RunnerContext) -> None:
        """Perform any activity required for starting the run here.
        For example, starting the target system to measure.
        Activities after starting the run should also be performed here."""
        pass

    def start_measurement(self, context: RunnerContext) -> None:
        """Start the measurement on Kelly’s laptop via SSH"""
        py_file = None
        for factor_name, factor_value in context.execute_run.items():
            if factor_name == "python_file":
                py_file = Path(self.BENCHMARKS_DIR / factor_value)
                break
        if not py_file or not py_file.exists():
            raise FileNotFoundError(f"Python file not found: {py_file}")

        rel_path = py_file.relative_to(self.BENCHMARKS_DIR)  # e.g. ant_colony_optimization/.../foo.py

        # --- Correct remote paths (keep benchmarks/ prefix) ---
        remote_script = f"{REMOTE_REPO}/benchmarks/{rel_path}"
        remote_results_dir = f"{REMOTE_REPO}/Results/{self.name}"
        remote_output = f"{remote_results_dir}/{py_file.stem}.csv"
        venv_python = f"{REMOTE_REPO}/venv/bin/python3"
        activate = f"{REMOTE_REPO}/venv/bin/activate"

        # Build robust remote command
        remote_cmd = (
            f"cd {shlex.quote(REMOTE_REPO)} && "
            f'export PATH=\"$HOME/.cargo/bin:$PATH\" && '
            f"if [ -x {shlex.quote(venv_python)} ]; then :; "
            f"elif [ -f {shlex.quote(activate)} ]; then . {shlex.quote(activate)}; "
            f"else echo 'ERROR: venv not found at {activate}' >&2; exit 1; fi && "
            f"sudo -n true || (echo 'ERROR: sudo -n not permitted; add NOPASSWD for /usr/bin/powermetrics' >&2; exit 1); "
            f"mkdir -p {shlex.quote(remote_results_dir)} && "
            f"sudo -n energibridge --summary -o {shlex.quote(remote_output)} "
            f"{shlex.quote(venv_python)} {shlex.quote(remote_script)}"
        )

        output.console_log(f"[start_measurement] Running benchmark: {py_file.name}")
        output.console_log(f"[start_measurement] SSH command: {remote_cmd}")

        # Use your 'macbook' SSH alias so no -i path issues
        self.profiler = subprocess.Popen(
            ["ssh", "macbook", "bash", "-lc", remote_cmd],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            bufsize=1,
            text=True,
        )

        with open(f'{context.run_dir}/energibridge.log', 'w') as logf:
            for line in self.profiler.stdout:
                print(f"[macbook] {line.strip()}", flush=True)
                logf.write(line)



    def interact(self, context: RunnerContext) -> None:
        """Wait for the energibridge process to finish"""
        output.console_log("[interact] Waiting for energibridge measurement to complete...")
        self.profiler.wait()
        output.console_log("[interact] Measurement completed.")

    def stop_measurement(self, context): 
        pass

    def stop_run(self, context: RunnerContext) -> None:
        """After the run, pull results from Kelly's Mac into Pi's Results/<experiment>/"""
        py_file = None
        for factor_name, factor_value in context.execute_run.items():
            if factor_name == "python_file":
                py_file = Path(self.BENCHMARKS_DIR / factor_value)
                break

        if not py_file:
            return

        # Remote + local paths
        remote_output = f"{REMOTE_REPO}/Results/{self.name}/{py_file.stem}.csv"
        local_results_dir = self.ROOT_DIR / "Results" / self.name
        local_results_dir.mkdir(parents=True, exist_ok=True)

        # Rsync from laptop to Pi
        subprocess.run(
            [
                "rsync", "-av",
                f"macbook:{remote_output}",
                str(local_results_dir / f"{py_file.stem}.csv"),
            ],
            check=True
        )
        output.console_log(f"[stop_run] Results copied to {local_results_dir / (py_file.stem + '.csv')}")



    
    def populate_run_data(self, context: RunnerContext) -> Optional[Dict[str, Any]]:
        """Parse and process measurement data from Results/<experiment>/<script>.csv"""
        py_file = None
        for factor_name, factor_value in context.execute_run.items():
            if factor_name == "python_file":
                py_file = Path(self.BENCHMARKS_DIR / factor_value)
                break

        if not py_file:
            return None

        # Path to CSV on the Pi
        output_file = self.ROOT_DIR / "Results" / self.name / f"{py_file.stem}.csv"

        if not output_file.exists():
            output.console_log(f"Warning: Output file not found: {output_file}")
            return {
                'total_energy_joules': 0,
                'execution_time_seconds': 0,
                'avg_power_watts': 0,
                'max_power_watts': 0,
                'min_power_watts': 0,
                'cpu_usage_avg': 0,
                'memory_usage_avg': 0
            }

        try:
            df = pd.read_csv(output_file)
            output.console_log(f"[populate_run_data] CSV loaded with {len(df)} rows and {len(df.columns)} columns.")

            # Calculate power stats
            if 'SYSTEM_POWER (Watts)' in df.columns:
                power_data = df['SYSTEM_POWER (Watts)'].dropna()
                avg_power = power_data.mean() if len(power_data) > 0 else 0
                max_power = power_data.max() if len(power_data) > 0 else 0
                min_power = power_data.min() if len(power_data) > 0 else 0
            else:
                avg_power = max_power = min_power = 0

            # CPU usage avg
            cpu_columns = [col for col in df.columns if col.startswith('CPU_USAGE_')]
            if cpu_columns:
                cpu_data = df[cpu_columns].mean(axis=1).dropna()
                cpu_usage_avg = cpu_data.mean() if len(cpu_data) > 0 else 0
            else:
                cpu_usage_avg = 0

            # Memory usage avg
            if 'USED_MEMORY' in df.columns and 'TOTAL_MEMORY' in df.columns:
                memory_usage = (df['USED_MEMORY'] / df['TOTAL_MEMORY'] * 100).dropna()
                memory_usage_avg = memory_usage.mean() if len(memory_usage) > 0 else 0
            else:
                memory_usage_avg = 0

            # Parse energibridge log for total energy & time
            total_energy = 0
            execution_time = 0
            log_file = context.run_dir / "energibridge.log"
            if log_file.exists():
                with open(log_file, 'r') as f:
                    for line in f:
                        if "Energy consumption in joules:" in line:
                            try:
                                total_energy = float(line.split(':')[1].split()[0])
                            except:
                                pass
                        elif "sec of execution" in line:
                            try:
                                execution_time = float(line.split()[0])
                            except:
                                pass

            run_data = {
                'total_energy_joules': round(total_energy, 3),
                'execution_time_seconds': round(execution_time, 3),
                'avg_power_watts': round(avg_power, 3),
                'max_power_watts': round(max_power, 3),
                'min_power_watts': round(min_power, 3),
                'cpu_usage_avg': round(cpu_usage_avg, 3),
                'memory_usage_avg': round(memory_usage_avg, 3)
            }

            output.console_log(
                f"Processed data for {py_file.name}: {total_energy}J, {execution_time}s, "
                f"avgP={avg_power}, cpu={cpu_usage_avg}, mem={memory_usage_avg}"
            )
            return run_data

        except Exception as e:
            output.console_log(f"Error processing data for {py_file.name}: {e}")
            return {
                'total_energy_joules': 0,
                'execution_time_seconds': 0,
                'avg_power_watts': 0,
                'max_power_watts': 0,
                'min_power_watts': 0,
                'cpu_usage_avg': 0,
                'memory_usage_avg': 0
            }



    def after_experiment(self) -> None:
        """Perform any activity required after stopping the experiment here
        Invoked only once during the lifetime of the program."""
        output.console_log("Benchmark energy analysis completed!")
        
        # Create summary report
        self._create_summary_report()

    def _create_summary_report(self):
        """Create a summary report of all experiments (current experiment only)"""
        summary_root = self.ROOT_DIR / "Results" / self.name
        summary_file = summary_root / "energy_analysis_summary.md"
        summary_root.mkdir(parents=True, exist_ok=True)

        with open(summary_file, 'w') as f:
            f.write(f"# Energy Analysis Summary — {self.name}\n\n")
            f.write("This report summarizes the energy consumption analysis for this experiment.\n\n")

            csv_files = sorted(summary_root.glob("*.csv"))
            if csv_files:
                f.write(f"Found {len(csv_files)} measurement files:\n")
                for csv_file in csv_files:
                    f.write(f"- {csv_file.name}\n")
            else:
                f.write("No measurement files found.\n")
            f.write("\n")


    # ================================ DO NOT ALTER BELOW THIS LINE ================================
    experiment_path: Path = None