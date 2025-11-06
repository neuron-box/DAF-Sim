"""Quick debug test"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from pillar3_physics_model import (
    Particle, Bubble, FluidProperties,
    calculate_attachment_efficiency
)

# Simple test
particle = Particle(
    diameter=10e-6,
    zeta_potential=-0.025,
    density=2650.0,
    hamaker_constant=5e-21
)

bubble = Bubble(diameter=50e-6, zeta_potential=-0.030)
fluid = FluidProperties.water_at_20C(ionic_strength=0.001)

print("Testing basic calculation...")
result = calculate_attachment_efficiency(particle, bubble, fluid, tolerance=1e-6)

print(f"alpha_bp: {result['alpha_bp']}")
print(f"eta_collision: {result['eta_collision']}")
print(f"alpha_attachment: {result['alpha_attachment']}")
print(f"critical_offset: {result['critical_offset']*1e6:.3f} um")
print(f"energy_barrier_kT: {result['energy_barrier_kT']:.2f}")
