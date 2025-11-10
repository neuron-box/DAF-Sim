"""
Example usage of the Core Physics Model (Pillar 3)

This script demonstrates how to calculate collision efficiency factors
for various DAF scenarios.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from pillar3_physics_model import (
    Particle, Bubble, FluidProperties,
    calculate_attachment_efficiency
)


def example_1_basic():
    """Basic example with typical DAF parameters."""
    print("\n" + "="*60)
    print("Example 1: Basic DAF Scenario")
    print("="*60)

    # Define particle
    particle = Particle(
        diameter=10e-6,           # 10 μm
        zeta_potential=-0.025,    # -25 mV
        density=2650.0,           # Silica-like
        hamaker_constant=5e-21    # Typical for hydrophobic particles
    )

    # Define bubble
    bubble = Bubble(
        diameter=50e-6,           # 50 μm
        zeta_potential=-0.030     # -30 mV (air-water interface)
    )

    # Define fluid
    fluid = FluidProperties.water_at_20C(
        ionic_strength=0.001,     # 0.001 M
        ph=7.0
    )

    # Calculate collision efficiency
    result = calculate_attachment_efficiency(particle, bubble, fluid)

    # Display results
    print(f"\nParticle: {particle.diameter*1e6:.1f} μm, ζ = {particle.zeta_potential*1e3:.1f} mV")
    print(f"Bubble: {bubble.diameter*1e6:.1f} μm, ζ = {bubble.zeta_potential*1e3:.1f} mV")
    print(f"Ionic strength: {fluid.ionic_strength:.4f} M")
    print(f"\nResults:")
    print(f"  α_bp (total) = {result['alpha_bp']:.4f}")
    print(f"  η_collision  = {result['eta_collision']:.4f}")
    print(f"  α_attachment = {result['alpha_attachment']:.4f}")
    print(f"  Energy barrier = {result['energy_barrier_kT']:.1f} kT")
    print(f"  Critical offset = {result['critical_offset']*1e6:.2f} μm")
    print(f"  Debye length = {result['debye_length']*1e9:.2f} nm")
    print(f"  Attachment favorable: {result['attachment_favorable']}")


def example_2_after_coagulation():
    """Example after coagulation treatment (reduced zeta potential)."""
    print("\n" + "="*60)
    print("Example 2: After Coagulation Treatment")
    print("="*60)

    # After coagulation: lower zeta potential, hydrophobic surface
    particle = Particle(
        diameter=15e-6,           # Larger floc
        zeta_potential=-0.010,    # Reduced to -10 mV
        density=1200.0,           # Lower density (floc structure)
        hamaker_constant=8e-21    # More hydrophobic
    )

    bubble = Bubble(
        diameter=40e-6,
        zeta_potential=-0.020
    )

    fluid = FluidProperties.water_at_20C(
        ionic_strength=0.01,      # Higher after coagulant addition
        ph=7.0
    )

    result = calculate_attachment_efficiency(particle, bubble, fluid)

    print(f"\nParticle: {particle.diameter*1e6:.1f} μm, ζ = {particle.zeta_potential*1e3:.1f} mV")
    print(f"Bubble: {bubble.diameter*1e6:.1f} μm, ζ = {bubble.zeta_potential*1e3:.1f} mV")
    print(f"Ionic strength: {fluid.ionic_strength:.4f} M")
    print(f"\nResults:")
    print(f"  α_bp (total) = {result['alpha_bp']:.4f}")
    print(f"  η_collision  = {result['eta_collision']:.4f}")
    print(f"  α_attachment = {result['alpha_attachment']:.4f}")
    print(f"  Energy barrier = {result['energy_barrier_kT']:.1f} kT")
    print(f"  Critical offset = {result['critical_offset']*1e6:.2f} μm")
    print(f"  Debye length = {result['debye_length']*1e9:.2f} nm")


def example_3_opposite_charges():
    """Example with opposite charges (particle positive, bubble negative)."""
    print("\n" + "="*60)
    print("Example 3: Opposite Charges (Favorable Attachment)")
    print("="*60)

    # Positive particle (e.g., after cationic polymer addition)
    particle = Particle(
        diameter=12e-6,
        zeta_potential=0.020,     # Positive!
        density=1500.0,
        hamaker_constant=8e-21
    )

    bubble = Bubble(
        diameter=45e-6,
        zeta_potential=-0.025     # Negative
    )

    fluid = FluidProperties.water_at_20C(
        ionic_strength=0.005,
        ph=7.0
    )

    result = calculate_attachment_efficiency(particle, bubble, fluid)

    print(f"\nParticle: {particle.diameter*1e6:.1f} μm, ζ = {particle.zeta_potential*1e3:.1f} mV")
    print(f"Bubble: {bubble.diameter*1e6:.1f} μm, ζ = {bubble.zeta_potential*1e3:.1f} mV")
    print(f"Ionic strength: {fluid.ionic_strength:.4f} M")
    print(f"\nResults:")
    print(f"  α_bp (total) = {result['alpha_bp']:.4f}")
    print(f"  η_collision  = {result['eta_collision']:.4f}")
    print(f"  α_attachment = {result['alpha_attachment']:.4f}")
    print(f"  Energy barrier = {result['energy_barrier_kT']:.1f} kT")
    print(f"  Critical offset = {result['critical_offset']*1e6:.2f} μm")
    print(f"  Debye length = {result['debye_length']*1e9:.2f} nm")


def example_4_parametric_study():
    """Parametric study: effect of particle size."""
    print("\n" + "="*60)
    print("Example 4: Parametric Study - Effect of Particle Size")
    print("="*60)

    bubble = Bubble(diameter=50e-6, zeta_potential=-0.030)
    fluid = FluidProperties.water_at_20C(ionic_strength=0.001)

    particle_sizes = [5e-6, 10e-6, 20e-6, 30e-6]  # 5, 10, 20, 30 μm

    print(f"\nBubble: {bubble.diameter*1e6:.1f} μm")
    print(f"Ionic strength: {fluid.ionic_strength:.4f} M")
    print(f"\n{'Size (μm)':<12} {'α_bp':<10} {'η_coll':<10} {'α_attach':<10} {'Barrier (kT)':<12}")
    print("-" * 60)

    for d_p in particle_sizes:
        particle = Particle(
            diameter=d_p,
            zeta_potential=-0.025,
            density=2650.0,
            hamaker_constant=5e-21
        )

        result = calculate_attachment_efficiency(particle, bubble, fluid)

        print(f"{d_p*1e6:<12.1f} {result['alpha_bp']:<10.4f} "
              f"{result['eta_collision']:<10.4f} {result['alpha_attachment']:<10.4f} "
              f"{result['energy_barrier_kT']:<12.1f}")


def example_5_ionic_strength_effect():
    """Parametric study: effect of ionic strength."""
    print("\n" + "="*60)
    print("Example 5: Parametric Study - Effect of Ionic Strength")
    print("="*60)

    particle = Particle(
        diameter=10e-6,
        zeta_potential=-0.025,
        density=2650.0,
        hamaker_constant=5e-21
    )

    bubble = Bubble(diameter=50e-6, zeta_potential=-0.030)

    ionic_strengths = [0.0001, 0.001, 0.01, 0.1]  # M

    print(f"\nParticle: {particle.diameter*1e6:.1f} μm, ζ = {particle.zeta_potential*1e3:.1f} mV")
    print(f"Bubble: {bubble.diameter*1e6:.1f} μm, ζ = {bubble.zeta_potential*1e3:.1f} mV")
    print(f"\n{'I (M)':<12} {'α_bp':<10} {'Barrier (kT)':<12} {'λ_D (nm)':<12}")
    print("-" * 50)

    for I in ionic_strengths:
        fluid = FluidProperties.water_at_20C(ionic_strength=I)
        result = calculate_attachment_efficiency(particle, bubble, fluid)

        print(f"{I:<12.4f} {result['alpha_bp']:<10.4f} "
              f"{result['energy_barrier_kT']:<12.1f} {result['debye_length']*1e9:<12.2f}")


if __name__ == "__main__":
    print("\n" + "="*60)
    print("PILLAR 3: CORE PHYSICS MODEL - EXAMPLE USAGE")
    print("Collision Efficiency Factor (α_bp) Calculation")
    print("Based on Han et al. (2001)")
    print("="*60)

    example_1_basic()
    example_2_after_coagulation()
    example_3_opposite_charges()
    example_4_parametric_study()
    example_5_ionic_strength_effect()

    print("\n" + "="*60)
    print("Examples completed successfully!")
    print("="*60 + "\n")
