"""
Example simulation demonstrating the Floc Kinetics PBM module.

This script shows how to:
1. Calculate aggregation and breakage kernels
2. Set up and solve the Population Balance Equation
3. Analyze the results
4. (Optional) Visualize the particle size distribution evolution
"""

import sys
import os
import numpy as np

# Add parent directory to path to import the module
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.kernels import calculate_floc_kernels
from src.pbe_solver import PopulationBalanceModel
from src.properties import FlocProperties


def example_kernel_calculation():
    """
    Example 1: Calculate aggregation and breakage kernels for given conditions.
    """
    print("=" * 70)
    print("Example 1: Kernel Calculations")
    print("=" * 70)

    # Define floc properties
    props = FlocProperties(
        collision_efficiency=0.1,
        binding_strength=1e3,
        fractal_dimension=2.3
    )

    # Flocculation conditions
    G = 100.0  # Shear rate [1/s]
    d_i = 50e-6  # Particle i: 50 μm
    d_j = 30e-6  # Particle j: 30 μm

    # Calculate kernels
    beta, S_i = calculate_floc_kernels(G, d_i, d_j, properties=props)

    print(f"\nConditions:")
    print(f"  Shear rate (G): {G:.1f} 1/s")
    print(f"  Particle i diameter: {d_i*1e6:.1f} μm")
    print(f"  Particle j diameter: {d_j*1e6:.1f} μm")
    print(f"  Collision efficiency (α): {props.collision_efficiency}")

    print(f"\nResults:")
    print(f"  Aggregation kernel β(i,j): {beta:.3e} m³/s")
    print(f"  Breakage rate S(i): {S_i:.3e} 1/s")

    # Calculate individual components
    from src.kernels import FlocKernels
    kernels = FlocKernels(props)

    beta_ortho = kernels.beta_orthokinetic(d_i, d_j, G)
    beta_peri = kernels.beta_perikinetic(d_i, d_j)
    beta_ds = kernels.beta_differential_sedimentation(d_i, d_j)

    print(f"\nAggregation kernel breakdown:")
    print(f"  Orthokinetic (turbulent shear): {beta_ortho:.3e} m³/s ({beta_ortho/(beta_ortho+beta_peri+beta_ds)*100:.1f}%)")
    print(f"  Perikinetic (Brownian motion): {beta_peri:.3e} m³/s ({beta_peri/(beta_ortho+beta_peri+beta_ds)*100:.1f}%)")
    print(f"  Differential sedimentation: {beta_ds:.3e} m³/s ({beta_ds/(beta_ortho+beta_peri+beta_ds)*100:.1f}%)")

    print()


def example_pbe_simulation():
    """
    Example 2: Solve the PBE for a flocculation scenario.
    """
    print("=" * 70)
    print("Example 2: PBE Simulation - Flocculation in Rapid Mixing")
    print("=" * 70)

    # Set up floc properties
    props = FlocProperties(
        collision_efficiency=0.1,
        binding_strength=1e3,
        fractal_dimension=2.3
    )

    # Initialize PBM solver
    pbm = PopulationBalanceModel(
        n_bins=25,
        d_min=1e-6,   # Minimum diameter: 1 μm
        d_max=500e-6,  # Maximum diameter: 500 μm
        properties=props
    )

    # Initial condition: narrow distribution of small particles
    print(f"\nInitialization:")
    print(f"  Number of size bins: {pbm.n_bins}")
    print(f"  Size range: {pbm.d_min*1e6:.1f} - {pbm.d_max*1e6:.1f} μm")

    n_initial = np.zeros(pbm.n_bins)
    n_initial[2:5] = 1e10  # 10^10 particles/m³ in bins 2-4

    stats_initial = pbm.summary_statistics(n_initial)
    print(f"\nInitial particle size distribution:")
    print(f"  Total number concentration: {stats_initial['total_number']:.3e} #/m³")
    print(f"  Volume fraction: {stats_initial['volume_fraction']:.3e}")
    print(f"  Mean diameter (d_43): {stats_initial['mean_diameter']*1e6:.2f} μm")
    print(f"  d10: {stats_initial['d10']*1e6:.2f} μm")
    print(f"  d50 (median): {stats_initial['d50']*1e6:.2f} μm")
    print(f"  d90: {stats_initial['d90']*1e6:.2f} μm")

    # Simulation parameters
    G = 100.0  # Shear rate [1/s] (typical for rapid mixing)
    t_end = 60.0  # Simulation time [s]
    n_time_points = 61

    print(f"\nSimulation parameters:")
    print(f"  Shear rate (G): {G:.1f} 1/s")
    print(f"  Simulation time: {t_end:.1f} s")

    # Solve PBE
    print(f"\nSolving PBE...")
    t_span = np.linspace(0, t_end, n_time_points)
    t, n_solution = pbm.solve(n_initial, t_span, G)
    print(f"  Simulation complete!")

    # Analyze final state
    n_final = n_solution[-1, :]
    stats_final = pbm.summary_statistics(n_final)

    print(f"\nFinal particle size distribution (t = {t_end:.1f} s):")
    print(f"  Total number concentration: {stats_final['total_number']:.3e} #/m³")
    print(f"  Volume fraction: {stats_final['volume_fraction']:.3e}")
    print(f"  Mean diameter (d_43): {stats_final['mean_diameter']*1e6:.2f} μm")
    print(f"  d10: {stats_final['d10']*1e6:.2f} μm")
    print(f"  d50 (median): {stats_final['d50']*1e6:.2f} μm")
    print(f"  d90: {stats_final['d90']*1e6:.2f} μm")

    # Calculate changes
    print(f"\nChanges during flocculation:")
    print(f"  Mean diameter increase: {(stats_final['mean_diameter']/stats_initial['mean_diameter'] - 1)*100:.1f}%")
    print(f"  d50 increase: {(stats_final['d50']/stats_initial['d50'] - 1)*100:.1f}%")
    print(f"  Number concentration decrease: {(1 - stats_final['total_number']/stats_initial['total_number'])*100:.1f}%")

    # Return data for optional plotting
    return t, n_solution, pbm


def example_parametric_study():
    """
    Example 3: Study the effect of shear rate on flocculation.
    """
    print("\n" + "=" * 70)
    print("Example 3: Parametric Study - Effect of Shear Rate")
    print("=" * 70)

    props = FlocProperties(collision_efficiency=0.1)
    pbm = PopulationBalanceModel(n_bins=20, d_min=1e-6, d_max=200e-6, properties=props)

    # Initial condition
    n_initial = np.zeros(pbm.n_bins)
    n_initial[4] = 1e10

    # Test different shear rates
    shear_rates = [10.0, 50.0, 100.0, 200.0, 500.0]
    t_span = np.linspace(0, 30, 31)

    print(f"\nTesting shear rates: {shear_rates} 1/s")
    print(f"Simulation time: {t_span[-1]:.1f} s")
    print(f"\nResults:")
    print(f"{'G [1/s]':<10} {'d_mean [μm]':<15} {'d50 [μm]':<15} {'N_total [#/m³]':<20}")
    print("-" * 70)

    results = {}
    for G in shear_rates:
        t, n_sol = pbm.solve(n_initial, t_span, G, recompute_kernels=True)
        stats = pbm.summary_statistics(n_sol[-1, :])
        results[G] = stats

        print(f"{G:<10.1f} {stats['mean_diameter']*1e6:<15.2f} "
              f"{stats['d50']*1e6:<15.2f} {stats['total_number']:<20.3e}")

    print("\nObservations:")
    print("  - At low shear (< 50 1/s): Slow aggregation dominates")
    print("  - At moderate shear (50-200 1/s): Fast aggregation increases floc size")
    print("  - At high shear (> 200 1/s): Breakage becomes significant")

    return results


def visualize_results(t, n_solution, pbm):
    """
    Optional: Visualize the results using matplotlib.
    """
    try:
        import matplotlib.pyplot as plt

        print("\n" + "=" * 70)
        print("Generating visualizations...")
        print("=" * 70)

        fig, axes = plt.subplots(2, 2, figsize=(12, 10))

        # Plot 1: Size distribution at different times
        ax = axes[0, 0]
        time_indices = [0, len(t)//4, len(t)//2, 3*len(t)//4, -1]
        colors = ['blue', 'green', 'orange', 'red', 'purple']

        for idx, color in zip(time_indices, colors):
            ax.semilogx(pbm.diameters*1e6, n_solution[idx, :],
                       label=f't = {t[idx]:.1f} s', color=color, linewidth=2)

        ax.set_xlabel('Particle Diameter [μm]', fontsize=12)
        ax.set_ylabel('Number Density [#/m³]', fontsize=12)
        ax.set_title('Evolution of Particle Size Distribution', fontsize=13, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)

        # Plot 2: Mean diameter evolution
        ax = axes[0, 1]
        mean_diameters = [pbm.mean_diameter(n_solution[i, :]) for i in range(len(t))]
        ax.plot(t, np.array(mean_diameters)*1e6, linewidth=2, color='darkblue')
        ax.set_xlabel('Time [s]', fontsize=12)
        ax.set_ylabel('Mean Diameter [μm]', fontsize=12)
        ax.set_title('Mean Diameter vs Time', fontsize=13, fontweight='bold')
        ax.grid(True, alpha=0.3)

        # Plot 3: Total number concentration
        ax = axes[1, 0]
        total_numbers = [pbm.total_number_concentration(n_solution[i, :]) for i in range(len(t))]
        ax.semilogy(t, total_numbers, linewidth=2, color='darkred')
        ax.set_xlabel('Time [s]', fontsize=12)
        ax.set_ylabel('Total Number Concentration [#/m³]', fontsize=12)
        ax.set_title('Total Number Concentration vs Time', fontsize=13, fontweight='bold')
        ax.grid(True, alpha=0.3)

        # Plot 4: Volume fraction
        ax = axes[1, 1]
        volume_fractions = [pbm.total_volume_fraction(n_solution[i, :]) for i in range(len(t))]
        ax.plot(t, volume_fractions, linewidth=2, color='darkgreen')
        ax.set_xlabel('Time [s]', fontsize=12)
        ax.set_ylabel('Volume Fraction [-]', fontsize=12)
        ax.set_title('Volume Fraction vs Time (Conservation Check)', fontsize=13, fontweight='bold')
        ax.grid(True, alpha=0.3)

        plt.tight_layout()

        # Save figure
        output_file = os.path.join(os.path.dirname(__file__), 'flocculation_results.png')
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"\nFigure saved to: {output_file}")

        # Try to show the plot
        plt.show()

    except ImportError:
        print("\nMatplotlib not available. Skipping visualization.")
        print("Install matplotlib to enable plotting: pip install matplotlib")


def main():
    """
    Run all examples.
    """
    print("\n" + "=" * 70)
    print(" FLOC KINETICS PBM - EXAMPLE SIMULATIONS")
    print("=" * 70)
    print()

    # Example 1: Kernel calculations
    example_kernel_calculation()

    # Example 2: PBE simulation
    t, n_solution, pbm = example_pbe_simulation()

    # Example 3: Parametric study
    example_parametric_study()

    # Optional: Visualization
    try:
        visualize_results(t, n_solution, pbm)
    except Exception as e:
        print(f"\nVisualization skipped: {e}")

    print("\n" + "=" * 70)
    print(" ALL EXAMPLES COMPLETED SUCCESSFULLY")
    print("=" * 70)
    print()


if __name__ == "__main__":
    main()
