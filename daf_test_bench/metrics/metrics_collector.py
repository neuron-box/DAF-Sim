"""
Metrics Collection and Reporting Module

Handles collection, storage, and export of benchmark metrics (FR-5).
Provides utilities to aggregate results from multiple engine runs and
export to various formats for analysis.
"""

import json
import csv
from pathlib import Path
from typing import List, Dict, Any, Optional
from datetime import datetime
import pandas as pd

from ..data_models import BenchmarkResult, ScientificMetrics, ComputationalMetrics


class MetricsCollector:
    """
    Collects and manages benchmark results from multiple engine runs.

    This class aggregates BenchmarkResult objects from the Test Harness
    and provides export functionality (FR-5.3).
    """

    def __init__(self):
        """Initialize the metrics collector."""
        self.results: List[BenchmarkResult] = []
        self._results_by_test: Dict[str, List[BenchmarkResult]] = {}
        self._results_by_engine: Dict[str, List[BenchmarkResult]] = {}

    def add_result(self, result: BenchmarkResult) -> None:
        """
        Add a benchmark result to the collection.

        Args:
            result: BenchmarkResult from an engine run
        """
        self.results.append(result)

        # Index by test name
        if result.test_name not in self._results_by_test:
            self._results_by_test[result.test_name] = []
        self._results_by_test[result.test_name].append(result)

        # Index by engine name
        if result.engine_name not in self._results_by_engine:
            self._results_by_engine[result.engine_name] = []
        self._results_by_engine[result.engine_name].append(result)

    def get_results_for_test(self, test_name: str) -> List[BenchmarkResult]:
        """
        Get all results for a specific test case.

        Args:
            test_name: Name of test case

        Returns:
            List of BenchmarkResult objects for this test
        """
        return self._results_by_test.get(test_name, [])

    def get_results_for_engine(self, engine_name: str) -> List[BenchmarkResult]:
        """
        Get all results for a specific engine.

        Args:
            engine_name: Name of engine

        Returns:
            List of BenchmarkResult objects for this engine
        """
        return self._results_by_engine.get(engine_name, [])

    def export_to_json(self, filepath: str, pretty: bool = True) -> None:
        """
        Export all results to a JSON file.

        Args:
            filepath: Output file path
            pretty: If True, format with indentation
        """
        data = {
            "export_timestamp": datetime.utcnow().isoformat() + "Z",
            "num_results": len(self.results),
            "results": [result.to_dict() for result in self.results]
        }

        with open(filepath, 'w') as f:
            json.dump(data, f, indent=2 if pretty else None)

    def export_to_csv(self, filepath: str) -> None:
        """
        Export all results to a CSV file.

        Creates a flattened table with one row per result. Scientific and
        computational metrics are expanded into separate columns.

        Args:
            filepath: Output CSV file path
        """
        if not self.results:
            # Create empty CSV with headers
            with open(filepath, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow([
                    'engine_name', 'test_name', 'timestamp', 'run_status',
                    'particle_removal_eff', 'velocity_error_l2', 'tke_error_l2',
                    'outlet_psd_error_wasserstein', 'mean_floc_size_d50',
                    'mean_residence_time', 'mass_conservation_error_pct',
                    'wall_clock_sec', 'cpu_hours', 'peak_ram_gb', 'converged',
                    'num_iterations'
                ])
            return

        rows = []
        for result in self.results:
            row = {
                'engine_name': result.engine_name,
                'test_name': result.test_name,
                'timestamp': result.timestamp,
                'run_status': result.run_status,
                # Scientific metrics
                'particle_removal_eff': result.scientific_metrics.particle_removal_eff,
                'velocity_error_l2': result.scientific_metrics.velocity_error_l2,
                'tke_error_l2': result.scientific_metrics.tke_error_l2,
                'outlet_psd_error_wasserstein': result.scientific_metrics.outlet_psd_error_wasserstein,
                'mean_floc_size_d50': result.scientific_metrics.mean_floc_size_d50,
                'mean_residence_time': result.scientific_metrics.mean_residence_time,
                'mass_conservation_error_pct': result.scientific_metrics.mass_conservation_error_pct,
                # Computational metrics
                'wall_clock_sec': result.computational_metrics.wall_clock_sec,
                'cpu_hours': result.computational_metrics.cpu_hours,
                'peak_ram_gb': result.computational_metrics.peak_ram_gb,
                'converged': result.computational_metrics.converged,
                'num_iterations': result.computational_metrics.num_iterations,
            }
            rows.append(row)

        # Write using pandas for clean CSV formatting
        df = pd.DataFrame(rows)
        df.to_csv(filepath, index=False)

    def export_summary_table(self, filepath: str) -> None:
        """
        Export a summary comparison table.

        Creates a pivot table with engines as rows and metrics as columns,
        useful for quick comparison across engines for a single test.

        Args:
            filepath: Output CSV file path
        """
        if not self.results:
            return

        rows = []
        for result in self.results:
            row = {
                'Engine': result.engine_name,
                'Test': result.test_name,
                'Status': result.run_status,
                'Removal_Eff': result.scientific_metrics.particle_removal_eff,
                'Vel_Error_L2': result.scientific_metrics.velocity_error_l2,
                'TKE_Error_L2': result.scientific_metrics.tke_error_l2,
                'PSD_Error': result.scientific_metrics.outlet_psd_error_wasserstein,
                'Wall_Clock_s': result.computational_metrics.wall_clock_sec,
                'Peak_RAM_GB': result.computational_metrics.peak_ram_gb,
                'Converged': result.computational_metrics.converged,
            }
            rows.append(row)

        df = pd.DataFrame(rows)
        df.to_csv(filepath, index=False)

    def generate_report(self, output_dir: str) -> None:
        """
        Generate a complete benchmark report.

        Creates multiple output files in the specified directory:
        - full_results.json: Complete results in JSON format
        - results_table.csv: Flattened results table
        - summary.csv: Summary comparison table
        - report.txt: Human-readable text report

        Args:
            output_dir: Directory to save report files
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        # Export JSON
        self.export_to_json(str(output_path / "full_results.json"))

        # Export CSV tables
        self.export_to_csv(str(output_path / "results_table.csv"))
        self.export_summary_table(str(output_path / "summary.csv"))

        # Generate text report
        self._generate_text_report(str(output_path / "report.txt"))

    def _generate_text_report(self, filepath: str) -> None:
        """Generate human-readable text report."""
        with open(filepath, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("DAF Test Bench - Benchmark Report\n")
            f.write("=" * 80 + "\n\n")
            f.write(f"Generated: {datetime.utcnow().isoformat()}Z\n")
            f.write(f"Total Results: {len(self.results)}\n\n")

            # Summary by test
            f.write("-" * 80 + "\n")
            f.write("Results by Test Case\n")
            f.write("-" * 80 + "\n\n")

            for test_name, results in self._results_by_test.items():
                f.write(f"Test: {test_name}\n")
                f.write(f"  Engines Run: {len(results)}\n")
                successful = sum(1 for r in results if r.run_status == "Success")
                f.write(f"  Successful: {successful}/{len(results)}\n\n")

                for result in results:
                    f.write(f"  - {result.engine_name}: {result.run_status}\n")
                    if result.run_status == "Success":
                        sci = result.scientific_metrics
                        comp = result.computational_metrics
                        if sci.particle_removal_eff is not None:
                            f.write(f"      Removal Eff: {sci.particle_removal_eff:.4f}\n")
                        if comp.wall_clock_sec is not None:
                            f.write(f"      Wall Clock: {comp.wall_clock_sec:.1f} s\n")
                f.write("\n")

            # Summary by engine
            f.write("-" * 80 + "\n")
            f.write("Results by Engine\n")
            f.write("-" * 80 + "\n\n")

            for engine_name, results in self._results_by_engine.items():
                f.write(f"Engine: {engine_name}\n")
                f.write(f"  Tests Run: {len(results)}\n")
                successful = sum(1 for r in results if r.run_status == "Success")
                f.write(f"  Successful: {successful}/{len(results)}\n")
                if successful > 0:
                    avg_time = sum(
                        r.computational_metrics.wall_clock_sec or 0
                        for r in results if r.run_status == "Success"
                    ) / successful
                    f.write(f"  Avg Wall Clock: {avg_time:.1f} s\n")
                f.write("\n")

            f.write("=" * 80 + "\n")

    def clear(self) -> None:
        """Clear all collected results."""
        self.results.clear()
        self._results_by_test.clear()
        self._results_by_engine.clear()
