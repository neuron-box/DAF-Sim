#!/usr/bin/env python3
"""
DAF Benchmark - Main Driver Script

Command-line interface for the DAF Test Bench framework.
Executes standardized benchmarks across multiple DAF Plant engines.

Usage:
    python daf_benchmark.py --config test_cases/T3_PBM/config.json
    python daf_benchmark.py --config test_cases/T2_Multiphase/config.json --engines Pillar3_Physics_Model
    python daf_benchmark.py --suite test_cases/ --output results/

Examples:
    # Run single test on all engines
    python daf_benchmark.py --config daf_test_bench/test_cases/T3_PBM/config.json

    # Run single test on specific engine
    python daf_benchmark.py --config daf_test_bench/test_cases/T3_PBM/config.json --engines Floc_Kinetics_PBM

    # Run full test suite
    python daf_benchmark.py --suite daf_test_bench/test_cases/ --output results/

    # Run only T3 tests
    python daf_benchmark.py --suite daf_test_bench/test_cases/ --filter T3_
"""

import argparse
import sys
import traceback
from pathlib import Path

from daf_test_bench import TestHarness
from daf_test_bench.engines import Pillar3PhysicsWrapper, FlocPBMWrapper


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="DAF Test Bench - Benchmark multiple DAF simulation engines",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    # Execution mode
    mode_group = parser.add_mutually_exclusive_group(required=True)
    mode_group.add_argument(
        '--config',
        type=str,
        help='Path to single test configuration file (config.json)'
    )
    mode_group.add_argument(
        '--suite',
        type=str,
        help='Path to test suite directory (runs all test cases)'
    )

    # Engine selection
    parser.add_argument(
        '--engines',
        type=str,
        nargs='+',
        help='List of engines to run (default: all registered engines)'
    )

    # Output options
    parser.add_argument(
        '--output',
        type=str,
        default='results',
        help='Output directory for results (default: results/)'
    )

    # Suite options
    parser.add_argument(
        '--filter',
        type=str,
        help='Filter pattern for test names (e.g., "T3_" for T3 tests only)'
    )

    # Verbosity
    parser.add_argument(
        '--quiet',
        action='store_true',
        help='Suppress progress output'
    )

    # List available engines
    parser.add_argument(
        '--list-engines',
        action='store_true',
        help='List available engines and exit'
    )

    args = parser.parse_args()

    # Create test harness
    harness = TestHarness()

    # Register available engines
    harness.register_engine(Pillar3PhysicsWrapper)
    harness.register_engine(FlocPBMWrapper)
    # Add more engines here as they are developed

    # List engines and exit if requested
    if args.list_engines:
        print("Available DAF Plant Engines:")
        for engine_name in harness.list_engines():
            print(f"  - {engine_name}")
        sys.exit(0)

    verbose = not args.quiet

    try:
        # Run benchmark(s)
        if args.config:
            # Single test case
            harness.run_benchmark(
                config_file=args.config,
                engines=args.engines,
                verbose=verbose
            )
        else:
            # Full test suite
            harness.run_test_suite(
                test_cases_dir=args.suite,
                engines=args.engines,
                test_filter=args.filter,
                verbose=verbose
            )

        # Save results
        if verbose:
            print()
        harness.save_results(args.output)

        if verbose:
            print()
            print("Benchmark complete!")
            print(f"Results saved to: {args.output}/")
            print()
            print("Output files:")
            print(f"  - {args.output}/full_results.json  (Complete JSON results)")
            print(f"  - {args.output}/results_table.csv  (Tabular format)")
            print(f"  - {args.output}/summary.csv        (Summary table)")
            print(f"  - {args.output}/report.txt         (Human-readable report)")

        sys.exit(0)

    except Exception as e:
        print(f"ERROR: {str(e)}", file=sys.stderr)
        if verbose:
            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
