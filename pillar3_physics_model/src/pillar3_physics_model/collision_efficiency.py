"""
Collision Efficiency Factor (α_bp) Calculation

This module implements the main function to calculate the collision efficiency
factor using a semi-analytical approach that combines classical single-collector
theory with DLVO force corrections.

The methodology is based on Han et al. (2001) but uses analytical approximations
for improved computational efficiency and numerical stability.

References:
    - Han, M.Y., Kim, W., & Dockko, S. (2001). Water Science and Technology, 43(8), 139-144.
    - Yoon, R. H., & Luttrell, G. H. (1989). Mineral Processing and Extractive Metallurgy Review.
"""

import math
from typing import Optional, Dict, Any
from .particle import Particle
from .bubble import Bubble
from .fluid_properties import FluidProperties
from .trajectory_solver import TrajectorySolver
from .dlvo_forces import DLVOForces


def calculate_attachment_efficiency(
    particle: Particle,
    bubble: Bubble,
    fluid_properties: FluidProperties,
    method: str = "trajectory",
    tolerance: float = 1e-8
) -> Dict[str, Any]:
    """
    Calculate the collision efficiency factor (α_bp) for bubble-particle attachment.

    The collision efficiency factor represents the fraction of particles that
    successfully attach to the bubble upon collision, accounting for both
    hydrodynamic and physicochemical (DLVO) forces.

    Methodology - Semi-Analytical Approach:
    ---------------------------------------
    This implementation uses a semi-analytical method that combines classical
    single-collector collision theory with DLVO force corrections. This approach
    provides computational efficiency and numerical stability compared to full
    trajectory integration.

    1. **Base Collision Efficiency** (η_base):
       - Interception mechanism: η_int = 1.5 × (R_p/R_b)²
       - Gravity settling: η_grav = N_G × (R_p/R_b)
       - Combined: η_base = η_int + η_grav

    2. **DLVO Energy Barrier Analysis**:
       - Van der Waals attraction (Hamaker constant)
       - Electrostatic repulsion (zeta potentials, Debye length)
       - Total interaction energy profile calculated

    3. **DLVO Correction Factor** (f_DLVO):
       - If barrier > 0: f_DLVO = exp(-E_barrier/3kT) (reduces efficiency)
       - If barrier < 0: f_DLVO > 1 (attractive, enhances efficiency)

    4. **Collision Efficiency**:
       η_c = η_base × f_DLVO

    5. **Attachment Efficiency** (α):
       - Based on energy barrier height
       - If E_barrier < 0: α = 1.0 (favorable)
       - If 0 < E_barrier < 10 kT: α = 1 - exp(-10/E_barrier_kT)
       - If E_barrier > 10 kT: α = exp(-E_barrier_kT/2)

    6. **Total Collision Efficiency Factor**:
       α_bp = η_c × α

    Parameters:
    -----------
    particle : Particle
        Particle object with properties:
        - diameter [m]
        - zeta_potential [V]
        - density [kg/m³]
        - hamaker_constant [J]

    bubble : Bubble
        Bubble object with properties:
        - diameter [m]
        - zeta_potential [V]
        - rise_velocity [m/s] (optional)

    fluid_properties : FluidProperties
        Fluid properties object with:
        - temperature [K]
        - dynamic_viscosity [Pa·s]
        - density [kg/m³]
        - ionic_strength [M]
        - dielectric_constant [-]
        - ph [-]

    method : str, optional
        Calculation method. Default is "trajectory" which uses the semi-analytical
        approach described above. This parameter is reserved for potential future
        implementations (e.g., full numerical trajectory integration).

    tolerance : float, optional
        Numerical tolerance parameter (reserved for future use).

    Returns:
    --------
    dict
        Dictionary containing:
        - 'alpha_bp': Total collision efficiency factor (0 to 1)
        - 'eta_collision': Collision efficiency from trajectory analysis (0 to 1)
        - 'alpha_attachment': Attachment efficiency from DLVO analysis (0 to 1)
        - 'critical_offset': Critical radial offset for collision [m]
        - 'energy_barrier': Maximum DLVO energy barrier [J]
        - 'energy_barrier_kT': Energy barrier in units of kT
        - 'attachment_favorable': Boolean indicating if attachment is favorable
        - 'debye_length': Debye length [m]
        - 'particle_reynolds': Particle Reynolds number
        - 'stokes_regime': Boolean indicating if Stokes regime applies

    Physical Equations:
    ------------------

    1. **Van der Waals Force** (Derjaguin approximation):

       F_vdW = -A * R_p * R_b / (6 * h² * (R_p + R_b))

       Where:
       - A: Hamaker constant [J]
       - R_p: Particle radius [m]
       - R_b: Bubble radius [m]
       - h: Surface-to-surface separation [m]

    2. **Electrostatic Double Layer Force** (constant potential):

       F_EDL = 2π * ε₀ * εᵣ * κ * R_p * R_b / (R_p + R_b) *
               [(ζ_p² + ζ_b²) * exp(-κh) - 2ζ_pζ_b * exp(-2κh)]

       Where:
       - ε₀: Vacuum permittivity = 8.854×10⁻¹² F/m
       - εᵣ: Relative dielectric constant (≈78.5 for water)
       - κ: Debye parameter = 1/λ_D [1/m]
       - ζ_p, ζ_b: Zeta potentials [V]
       - h: Surface-to-surface separation [m]

    3. **Debye Length** (temperature-dependent):

       λ_D = √(ε₀ε_r k_B T / (2 N_A e² I))

       Where:
       - ε₀ = vacuum permittivity
       - ε_r = dielectric constant
       - k_B = Boltzmann constant
       - T = temperature
       - I = ionic strength [mol/L]

       (Simplified form at 25°C: λ_D ≈ 0.304/√I [nm])

    4. **Hydrodynamic Drag Force** (Stokes):

       F_drag = 6π * μ * R_p * v_rel

       Where:
       - μ: Dynamic viscosity [Pa·s]
       - R_p: Particle radius [m]
       - v_rel: Relative velocity between particle and fluid [m/s]

    5. **Collision Efficiency**:

       η_c = (y_c / R_b)²

       Where y_c is the critical radial offset

    6. **Total Efficiency**:

       α_bp = η_c * α_attachment

    Example:
    --------
    >>> from pillar3_physics_model import (
    ...     Particle, Bubble, FluidProperties,
    ...     calculate_attachment_efficiency
    ... )
    >>>
    >>> # Define particle (10 μm, negative zeta potential)
    >>> particle = Particle(
    ...     diameter=10e-6,
    ...     zeta_potential=-0.025,
    ...     density=2650.0,
    ...     hamaker_constant=5e-21
    ... )
    >>>
    >>> # Define bubble (50 μm, negative zeta potential)
    >>> bubble = Bubble(
    ...     diameter=50e-6,
    ...     zeta_potential=-0.030
    ... )
    >>>
    >>> # Define fluid properties (water at 20°C, 0.001 M ionic strength)
    >>> fluid = FluidProperties.water_at_20C(ionic_strength=0.001)
    >>>
    >>> # Calculate collision efficiency
    >>> result = calculate_attachment_efficiency(particle, bubble, fluid)
    >>> print(f"α_bp = {result['alpha_bp']:.4f}")
    >>> print(f"Energy barrier = {result['energy_barrier_kT']:.1f} kT")

    Notes:
    ------
    - This implementation uses a **semi-analytical approach** rather than full
      numerical trajectory integration for improved robustness and efficiency
    - The method assumes Stokes flow regime (low Reynolds number, Re < 0.1)
    - DLVO theory assumes additivity of van der Waals and electrostatic forces
    - Zeta potential is approximated as surface potential
    - Results are most accurate for particles and bubbles in the 1-100 μm range
    - The TrajectorySolver class is available for future trajectory-based methods
      but is not used in the current implementation
    """

    if method != "trajectory":
        raise ValueError(f"Unknown method: {method}. Only 'trajectory' is supported.")

    # Initialize DLVO force calculator
    dlvo = DLVOForces(particle, bubble, fluid_properties)

    # Note: TrajectorySolver is not used in the current semi-analytical implementation
    # It is available for future trajectory-based methods if needed
    # solver = TrajectorySolver(particle, bubble, fluid_properties)

    # Check if attachment is energetically favorable
    attachment_favorable = dlvo.is_favorable_for_attachment()

    # Calculate energy barrier
    h_range = [1e-9 + i * 99e-9 / 100 for i in range(101)]  # 1 to 100 nm
    energies = [dlvo.interaction_energy(h)[2] for h in h_range]
    max_energy = max(energies)
    kT = fluid_properties.boltzmann_temperature
    energy_barrier_kT = max_energy / kT

    # Determine attachment efficiency based on energy barrier
    if energy_barrier_kT > 10.0:
        # High energy barrier -> low attachment probability
        alpha_attachment = math.exp(-energy_barrier_kT / 2.0)  # Approximate
    else:
        # Low energy barrier -> high attachment probability
        alpha_attachment = 1.0 - math.exp(-10.0 / energy_barrier_kT) if energy_barrier_kT > 0 else 1.0

    # Calculate collision efficiency using semi-analytical approach
    # This combines the classic single-collector collision efficiency
    # with DLVO corrections
    try:
        # Use simplified but robust analytical approximation
        # Based on Yoon & Luttrell (1989) and Han et al. (2001)

        # Size ratio
        N_R = particle.radius / bubble.radius

        # Stokes number (ratio of inertial to viscous forces)
        v_bubble = bubble.rise_velocity or bubble.calculate_rise_velocity(
            fluid_properties.dynamic_viscosity, fluid_properties.density
        )
        St = (2.0 / 9.0) * (particle.density - fluid_properties.density) * \
             particle.diameter * v_bubble / (fluid_properties.dynamic_viscosity * bubble.diameter)

        # Base collision efficiency from interception and inertia
        # For small particles in DAF, interception dominates
        eta_interception = 1.5 * N_R**2  # Interception mechanism

        # Gravity settling contribution (usually small for micron particles)
        N_G = (particle.density - fluid_properties.density) * particle.diameter**2 * 9.81 / \
              (18.0 * fluid_properties.dynamic_viscosity * v_bubble)
        eta_gravity = N_G * N_R

        # Combined collision efficiency (addition for independent mechanisms)
        eta_base = eta_interception + eta_gravity

        # Apply DLVO correction factor
        # If energy barrier is high, reduce collision efficiency
        if energy_barrier_kT > 0:
            # Barrier reduction factor
            dlvo_factor = math.exp(-energy_barrier_kT / 3.0)
        else:
            # Attractive interaction enhances collision
            dlvo_factor = 1.0 + 0.5 * abs(energy_barrier_kT) / 10.0

        eta_collision = eta_base * dlvo_factor

        # Collision efficiency should not exceed 1.0
        eta_collision = min(max(eta_collision, 0.0), 1.0)

        # Calculate approximate critical offset
        critical_offset = bubble.radius * math.sqrt(eta_collision)

    except Exception as e:
        # If calculation fails, use simplest estimate
        print(f"Warning: Collision efficiency calculation failed: {e}")
        print("Using simplified collision efficiency estimate.")

        # Simplified estimate based on interception
        eta_collision = 1.5 * (particle.radius / bubble.radius) ** 2
        eta_collision = min(eta_collision, 1.0)
        # Calculate critical_offset consistent with η_c = (y_c / R_b)²
        critical_offset = bubble.radius * math.sqrt(eta_collision)

    # Total collision efficiency factor
    alpha_bp = eta_collision * alpha_attachment

    # Calculate Reynolds number for validation
    v_bubble = bubble.rise_velocity or bubble.calculate_rise_velocity(
        fluid_properties.dynamic_viscosity, fluid_properties.density
    )
    Re_p = (fluid_properties.density * particle.diameter * v_bubble /
            fluid_properties.dynamic_viscosity)

    # Return comprehensive results
    return {
        'alpha_bp': alpha_bp,
        'eta_collision': eta_collision,
        'alpha_attachment': alpha_attachment,
        'critical_offset': critical_offset,
        'energy_barrier': max_energy,
        'energy_barrier_kT': energy_barrier_kT,
        'attachment_favorable': attachment_favorable,
        'debye_length': fluid_properties.debye_length,
        'particle_reynolds': Re_p,
        'stokes_regime': Re_p < 0.1,
        'bubble_rise_velocity': v_bubble,
    }
