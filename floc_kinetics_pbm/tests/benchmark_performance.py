"""
Performance benchmark for PBM solver optimizations.

This script benchmarks the performance improvement from the O(n³) to O(n²)
optimization in the aggregation_birth calculation.
"""

import numpy as np
import time
from src.pbe_solver import PopulationBalanceModel
from src.properties import FlocProperties

def benchmark_pbe_solve(n_bins: int, n_time_steps: int = 10) -> float:
    """
    Benchmark PBM solve time for given bin count.

    Args:
        n_bins: Number of size bins
        n_time_steps: Number of time steps to simulate

    Returns:
        Average time per time step [seconds]
    """
    # Initialize PBM
    props = FlocProperties(collision_efficiency=0.1)
    pbm = PopulationBalanceModel(
        n_bins=n_bins,
        d_min=1e-6,
        d_max=500e-6,
        properties=props
    )

    # Initial condition
    n_initial = np.zeros(pbm.n_bins)
    n_initial[n_bins//4:n_bins//3] = 1e9

    # Time span
    t_span = np.linspace(0, 10, n_time_steps)

    # Benchmark
    G = 100.0
    start = time.perf_counter()
    t, n_solution = pbm.solve(n_initial, t_span, G)
    elapsed = time.perf_counter() - start

    avg_time_per_step = elapsed / n_time_steps
    return avg_time_per_step


def benchmark_suite():
    """
    Run comprehensive performance benchmarks.
    """
    print("=" * 70)
    print(" PERFORMANCE BENCHMARK - PBM Solver Optimization")
    print("=" * 70)
    print()
    print("Testing the O(n³) → O(n²) optimization in aggregation_birth()")
    print()

    # Test different bin counts
    bin_counts = [10, 15, 20, 25, 30, 40, 50]
    results = []

    print(f"{'N_bins':<10} {'Time/step (ms)':<20} {'Total Time (s)':<20} {'Speedup':<10}")
    print("-" * 70)

    baseline = None
    for n_bins in bin_counts:
        time_per_step = benchmark_pbe_solve(n_bins, n_time_steps=11)
        total_time = time_per_step * 11

        # Calculate expected O(n²) scaling
        if baseline is None:
            baseline = (n_bins, time_per_step)
            expected_ratio = 1.0
        else:
            # For O(n²), time should scale as (n/n0)²
            expected_ratio = (n_bins / baseline[0])**2

        actual_ratio = time_per_step / baseline[1]

        print(f"{n_bins:<10} {time_per_step*1000:<20.3f} {total_time:<20.3f} {actual_ratio:<10.2f}x")

        results.append({
            'n_bins': n_bins,
            'time_per_step': time_per_step,
            'total_time': total_time,
            'actual_ratio': actual_ratio,
            'expected_ratio': expected_ratio
        })

    print()
    print("Analysis:")
    print("-" * 70)

    # Calculate average scaling exponent
    if len(results) > 2:
        # Fit power law: time = a * n_bins^b
        # log(time) = log(a) + b*log(n_bins)
        log_n = np.log([r['n_bins'] for r in results])
        log_t = np.log([r['time_per_step'] for r in results])

        # Linear regression
        A = np.vstack([log_n, np.ones(len(log_n))]).T
        b, log_a = np.linalg.lstsq(A, log_t, rcond=None)[0]

        print(f"Measured scaling exponent: {b:.2f}")
        print(f"Expected for O(n²): 2.00")
        print(f"Expected for O(n³): 3.00")
        print()

        if abs(b - 2.0) < 0.3:
            print("✓ Performance scales as O(n²) - Optimization SUCCESSFUL!")
        elif abs(b - 3.0) < 0.3:
            print("✗ Performance scales as O(n³) - Optimization NOT applied")
        else:
            print(f"⚠ Performance scales as O(n^{b:.2f}) - Between quadratic and cubic")

    print()
    print("Performance for realistic simulation:")
    print("-" * 70)

    # Realistic simulation: 25 bins, 60 seconds with 1-second time steps
    realistic_bins = 25
    realistic_steps = 61
    time_per_step = benchmark_pbe_solve(realistic_bins, n_time_steps=realistic_steps)
    total_time = time_per_step * realistic_steps

    print(f"Configuration: {realistic_bins} bins, 60 seconds, 1-second steps")
    print(f"Time per time step: {time_per_step*1000:.3f} ms")
    print(f"Total simulation time: {total_time:.3f} seconds")
    print()

    # Estimate for larger simulations
    large_bins = 100
    estimated_time_per_step = time_per_step * (large_bins / realistic_bins)**2
    estimated_total = estimated_time_per_step * realistic_steps

    print(f"Estimated for {large_bins} bins:")
    print(f"Time per time step: {estimated_time_per_step*1000:.3f} ms")
    print(f"Total simulation time: {estimated_total:.3f} seconds")
    print()

    print("=" * 70)
    print(" BENCHMARK COMPLETE")
    print("=" * 70)


if __name__ == "__main__":
    benchmark_suite()
