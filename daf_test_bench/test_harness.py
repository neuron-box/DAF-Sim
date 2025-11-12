"""
Test Harness - Main Orchestration Module

The Test Harness is the central driver that executes benchmarks across
multiple DAF Plant engines using the standardized IDAFPlant interface.

Implements FR-2 from DAF-TB-SRR-v1.0.
"""

import time
from pathlib import Path
from typing import List, Dict, Any, Optional, Type
import json

from .interfaces.idaf_plant import IDAFPlant
from .data_models import TestConfiguration, BenchmarkResult
from .metrics.metrics_collector import MetricsCollector


class TestHarness:
    """
    Main test orchestration class.

    The TestHarness loads test configurations, instantiates engine wrappers,
    executes simulations, and collects results for analysis.
    """

    def __init__(self):
        """Initialize the Test Harness."""
        self.metrics_collector = MetricsCollector()
        self._registered_engines: Dict[str, Type[IDAFPlant]] = {}
        self._current_test: Optional[str] = None

    def register_engine(self, engine_class: Type[IDAFPlant]) -> None:
        """
        Register a DAF Plant engine class.

        Args:
            engine_class: Class implementing IDAFPlant interface
        """
        # Get engine name from a temporary instance
        temp_instance = engine_class()
        engine_name = temp_instance.engine_name
        self._registered_engines[engine_name] = engine_class
        print(f"Registered engine: {engine_name}")

    def list_engines(self) -> List[str]:
        """
        Get list of registered engine names.

        Returns:
            List of engine names
        """
        return list(self._registered_engines.keys())

    def run_benchmark(
        self,
        config_file: str,
        engines: Optional[List[str]] = None,
        verbose: bool = True
    ) -> MetricsCollector:
        """
        Run a benchmark test on specified engines.

        This is the main entry point for executing benchmarks. It follows
        the workflow defined in Section 5.0 of DAF-TB-HLD-v1.0:
        1. Load TestConfiguration
        2. For each engine:
            a. Instantiate engine wrapper
            b. Call setup()
            c. Call initialize()
            d. Call run()
            e. Call get_metrics()
            f. Call finalize()
        3. Aggregate and return results

        Args:
            config_file: Path to config.json file
            engines: List of engine names to run (None = all registered)
            verbose: Print progress messages

        Returns:
            MetricsCollector with all benchmark results
        """
        # Load configuration
        if verbose:
            print("=" * 80)
            print("DAF Test Bench - Starting Benchmark")
            print("=" * 80)
            print(f"Loading configuration: {config_file}")

        config = TestConfiguration.from_json_file(config_file)
        self._current_test = config.test_name

        if verbose:
            print(f"Test: {config.test_name}")
            print(f"Description: {config.description}")
            print()

        # Determine which engines to run
        if engines is None:
            engines_to_run = self.list_engines()
        else:
            engines_to_run = engines

        if not engines_to_run:
            raise ValueError("No engines to run. Register engines first.")

        # Run each engine
        for engine_name in engines_to_run:
            if engine_name not in self._registered_engines:
                print(f"WARNING: Engine '{engine_name}' not registered. Skipping.")
                continue

            if verbose:
                print("-" * 80)
                print(f"Running: {engine_name}")
                print("-" * 80)

            result = self._run_single_engine(
                engine_name,
                config,
                verbose=verbose
            )

            # Add to metrics collector
            result_obj = BenchmarkResult.from_dict(result)
            self.metrics_collector.add_result(result_obj)

            if verbose:
                print(f"Status: {result['run_status']}")
                if result['run_status'] == 'Success':
                    comp = result['computational_metrics']
                    if comp.get('wall_clock_sec'):
                        print(f"Wall Clock: {comp['wall_clock_sec']:.3f} s")
                print()

        if verbose:
            print("=" * 80)
            print("Benchmark Complete")
            print("=" * 80)
            print(f"Total runs: {len(self.metrics_collector.results)}")
            successful = sum(
                1 for r in self.metrics_collector.results
                if r.run_status == "Success"
            )
            print(f"Successful: {successful}/{len(self.metrics_collector.results)}")
            print()

        return self.metrics_collector

    def _run_single_engine(
        self,
        engine_name: str,
        config: TestConfiguration,
        verbose: bool = True
    ) -> Dict[str, Any]:
        """
        Run a single engine on a test configuration.

        Args:
            engine_name: Name of engine to run
            config: Test configuration
            verbose: Print progress

        Returns:
            BenchmarkResult dictionary
        """
        start_time = time.time()

        try:
            # Instantiate engine
            engine_class = self._registered_engines[engine_name]
            engine = engine_class()

            if verbose:
                print(f"  1. Setting up...")

            # Setup
            setup_success = engine.setup(config.to_dict())
            if not setup_success:
                raise RuntimeError(f"Setup failed for {engine_name}")

            if verbose:
                print(f"  2. Initializing...")

            # Initialize
            init_success = engine.initialize()
            if not init_success:
                raise RuntimeError(f"Initialization failed for {engine_name}")

            if verbose:
                print(f"  3. Running simulation...")

            # Run
            run_success = engine.run()
            if not run_success:
                raise RuntimeError(f"Simulation failed for {engine_name}")

            if verbose:
                print(f"  4. Collecting metrics...")

            # Get metrics
            metrics = engine.get_metrics()

            if verbose:
                print(f"  5. Finalizing...")

            # Finalize
            engine.finalize()

            return metrics

        except Exception as e:
            # Create failure result
            error_msg = f"Engine execution failed: {str(e)}"
            if verbose:
                print(f"  ERROR: {error_msg}")

            elapsed = time.time() - start_time

            result = BenchmarkResult.create(
                engine_name=engine_name,
                test_name=config.test_name,
                run_status="Failed",
                run_log=error_msg
            )

            # Add elapsed time
            result.computational_metrics.wall_clock_sec = elapsed

            return result.to_dict()

    def run_test_suite(
        self,
        test_cases_dir: str,
        engines: Optional[List[str]] = None,
        test_filter: Optional[str] = None,
        verbose: bool = True
    ) -> MetricsCollector:
        """
        Run a complete test suite (T1-T5).

        Args:
            test_cases_dir: Directory containing test case subdirectories
            engines: List of engine names to run (None = all)
            test_filter: Filter pattern for test names (e.g., "T3_" for only T3 tests)
            verbose: Print progress

        Returns:
            MetricsCollector with all results
        """
        test_dir = Path(test_cases_dir)

        if not test_dir.exists():
            raise FileNotFoundError(f"Test cases directory not found: {test_cases_dir}")

        # Find all config.json files
        config_files = list(test_dir.glob("*/config.json"))

        if not config_files:
            raise FileNotFoundError(f"No config.json files found in {test_cases_dir}")

        # Apply filter
        if test_filter:
            config_files = [f for f in config_files if test_filter in f.parent.name]

        if verbose:
            print(f"Found {len(config_files)} test case(s)")
            print()

        # Run each test case
        for config_file in sorted(config_files):
            self.run_benchmark(
                str(config_file),
                engines=engines,
                verbose=verbose
            )

        return self.metrics_collector

    def save_results(self, output_dir: str) -> None:
        """
        Save all benchmark results to files.

        Creates:
        - full_results.json: Complete results in JSON
        - results_table.csv: Tabular results
        - summary.csv: Summary table
        - report.txt: Human-readable report

        Args:
            output_dir: Output directory path
        """
        self.metrics_collector.generate_report(output_dir)
        print(f"Results saved to: {output_dir}")

    def clear_results(self) -> None:
        """Clear all collected results."""
        self.metrics_collector.clear()
