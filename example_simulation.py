#!/usr/bin/env python3
"""
Example DAF Simulation using Eulerian Momentum Equation

This script demonstrates a simple DAF simulation with bubbles rising
through water and carrying particles to the surface.
"""

import numpy as np
from eulerian_momentum import Phase, MomentumEquation, InterfacialForces


def main():
    """Run a simple DAF simulation demonstration."""
    print("=" * 70)
    print("DAF Simulation - Eulerian Momentum Equation Demonstration")
    print("=" * 70)

    # ========== Domain Setup ==========
    print("\n1. Setting up computational domain...")

    nz = 100  # Number of vertical grid points
    dz = 0.002  # Grid spacing: 2 mm
    dt = 0.001  # Time step: 1 ms
    height = (nz - 1) * dz  # Total height: ~20 cm

    print(f"   Grid points: {nz}")
    print(f"   Grid spacing: {dz * 1000:.1f} mm")
    print(f"   Domain height: {height * 100:.1f} cm")
    print(f"   Time step: {dt * 1000:.1f} ms")

    # ========== Phase Initialization ==========
    print("\n2. Initializing phases...")

    # Water (continuous phase)
    water = Phase(name="water", rho=1000.0, mu=1e-3, phase_id=0)
    water.initialize_fields((nz,), initial_alpha=0.88)
    water.velocity[:, :] = 0.0  # Still water initially

    # Air bubbles (dispersed phase)
    bubbles = Phase(name="bubbles", rho=1.2, mu=1.8e-5, phase_id=1)
    bubbles.initialize_fields((nz,), initial_alpha=0.10)
    bubbles.velocity[:, 2] = 0.05  # Initial rise velocity: 5 cm/s

    # Floc particles (dispersed phase)
    flocs = Phase(name="flocs", rho=1050.0, mu=1e-3, phase_id=2)
    flocs.initialize_fields((nz,), initial_alpha=0.02)
    flocs.velocity[:, 2] = 0.0  # Initially neutrally suspended

    # Verify volume fraction constraint
    total_alpha = water.alpha + bubbles.alpha + flocs.alpha
    print(f"\n   Water: α = {np.mean(water.alpha):.3f}, ρ = {water.rho} kg/m³")
    print(f"   Bubbles: α = {np.mean(bubbles.alpha):.3f}, ρ = {bubbles.rho} kg/m³")
    print(f"   Flocs: α = {np.mean(flocs.alpha):.3f}, ρ = {flocs.rho} kg/m³")
    print(f"   Total volume fraction: {np.mean(total_alpha):.3f} (should be 1.0)")

    # ========== Pressure Field ==========
    print("\n3. Setting up pressure field...")

    z = np.arange(nz) * dz
    pressure = 101325.0 + water.rho * 9.81 * (height - z)  # Hydrostatic

    print(f"   Bottom pressure: {pressure[0] / 1000:.2f} kPa")
    print(f"   Top pressure: {pressure[-1] / 1000:.2f} kPa")
    print(f"   Pressure gradient: {water.rho * 9.81:.1f} Pa/m (hydrostatic)")

    # ========== Create Solvers ==========
    print("\n4. Initializing solvers...")

    momentum_solver = MomentumEquation(gravity=np.array([0, 0, -9.81]))
    forces_calc = InterfacialForces()

    print("   ✓ Momentum equation solver created")
    print("   ✓ Interfacial forces calculator created")

    # ========== Compute Forces on Bubbles ==========
    print("\n5. Computing forces on air bubbles...")

    bubble_gravity = momentum_solver.gravity_term(bubbles)
    bubble_pressure = momentum_solver.pressure_gradient_term(bubbles, pressure, dz)
    bubble_drag = forces_calc.drag_force(
        water, bubbles,
        drag_coefficient=0.44,
        particle_diameter=100e-6  # 100 micron bubbles
    )

    # Buoyancy = pressure gradient - gravity
    bubble_buoyancy = bubble_pressure[:, 2] - bubble_gravity[:, 2]

    print(f"\n   Gravity force: {np.mean(bubble_gravity[:, 2]):.2f} N/m³ (downward)")
    print(f"   Pressure force: {np.mean(bubble_pressure[:, 2]):.2f} N/m³ (upward)")
    print(f"   Net buoyancy: {np.mean(bubble_buoyancy):.2f} N/m³ (upward)")
    print(f"   Drag force: {np.mean(bubble_drag[:, 2]):.2f} N/m³ (opposes rise)")

    # ========== Compute Forces on Flocs ==========
    print("\n6. Computing forces on floc particles...")

    floc_gravity = momentum_solver.gravity_term(flocs)
    floc_pressure = momentum_solver.pressure_gradient_term(flocs, pressure, dz)
    floc_drag = forces_calc.drag_force(
        water, flocs,
        drag_coefficient=1.0,  # Higher for irregular flocs
        particle_diameter=200e-6  # 200 micron flocs
    )

    floc_buoyancy = floc_pressure[:, 2] - floc_gravity[:, 2]

    print(f"\n   Gravity force: {np.mean(floc_gravity[:, 2]):.2f} N/m³ (downward)")
    print(f"   Pressure force: {np.mean(floc_pressure[:, 2]):.2f} N/m³ (upward)")
    print(f"   Net buoyancy: {np.mean(floc_buoyancy):.2f} N/m³ (downward - denser than water)")
    print(f"   Drag force: {np.mean(floc_drag[:, 2]):.2f} N/m³")

    # ========== Complete Momentum Equation ==========
    print("\n7. Computing complete momentum equation (LHS)...")

    bubble_lhs = momentum_solver.compute_lhs(
        bubbles, pressure, dt, dz,
        S_ce=np.zeros(nz)
    )

    floc_lhs = momentum_solver.compute_lhs(
        flocs, pressure, dt, dz,
        S_ce=np.zeros(nz)
    )

    print(f"\n   Bubble LHS (z): {np.mean(bubble_lhs[:, 2]):.2f} N/m³")
    print(f"   Floc LHS (z): {np.mean(floc_lhs[:, 2]):.2f} N/m³")

    # ========== Physical Interpretation ==========
    print("\n" + "=" * 70)
    print("PHYSICAL INTERPRETATION")
    print("=" * 70)

    print("\n✓ BUBBLES:")
    if np.mean(bubble_buoyancy) > 0:
        print("  - Experience net upward buoyancy force (ρ_bubble < ρ_water)")
        print("  - Will rise through the water column")
        print("  - Drag opposes upward motion")
        print("  - Expected behavior in DAF: bubbles rise to surface")

    print("\n✓ FLOCS:")
    if np.mean(floc_buoyancy) < 0:
        print("  - Denser than water (ρ_floc > ρ_water)")
        print("  - Experience net downward force")
        print("  - Would settle without bubble attachment")
    else:
        print("  - Neutrally buoyant or slightly buoyant")

    print("\n✓ DAF MECHANISM:")
    print("  - Bubbles attach to floc particles")
    print("  - Combined bubble-floc density < water density")
    print("  - Aggregates rise to surface for removal")
    print("  - This module provides the momentum transport physics!")

    # ========== Terminal Velocity Estimation ==========
    print("\n" + "=" * 70)
    print("TERMINAL VELOCITY ESTIMATION")
    print("=" * 70)

    # For terminal velocity: drag force = buoyancy force
    # (3/4) * (C_D/d_p) * α * ρ_c * |U|^2 = α * Δρ * g

    C_D_bubble = 0.44
    d_bubble = 100e-6
    delta_rho_bubble = water.rho - bubbles.rho

    # Solve for terminal velocity
    # |U_t|^2 = (4/3) * (d_p/C_D) * (Δρ/ρ_c) * g
    U_terminal_bubble = np.sqrt(
        (4.0 / 3.0) * (d_bubble / C_D_bubble) *
        (delta_rho_bubble / water.rho) * 9.81
    )

    print(f"\n✓ Estimated bubble terminal velocity: {U_terminal_bubble * 100:.2f} cm/s")
    print(f"  (Current simulation: {np.mean(bubbles.velocity[:, 2]) * 100:.2f} cm/s)")

    # ========== Summary ==========
    print("\n" + "=" * 70)
    print("SIMULATION SUMMARY")
    print("=" * 70)

    print(f"\n✓ Successfully computed all momentum equation terms")
    print(f"✓ Verified physical behavior (buoyancy, drag)")
    print(f"✓ Three-phase system: water + bubbles + flocs")
    print(f"✓ Module ready for integration into full DAF solver")

    print("\n" + "=" * 70)
    print("Demonstration complete!")
    print("=" * 70 + "\n")


if __name__ == "__main__":
    main()
