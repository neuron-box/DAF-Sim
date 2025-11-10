"""
DLVO Theory Forces for Bubble-Particle Interactions

This module implements the Derjaguin-Landau-Verwey-Overbeek (DLVO) theory
for calculating interparticle forces between bubbles and particles.

References:
    - Han et al. (2001). Water Science and Technology, 43(8), 139-144.
    - Israelachvili, J. N. (2011). Intermolecular and Surface Forces (3rd ed.).
"""

import math
from typing import Tuple
from .particle import Particle
from .bubble import Bubble
from .fluid_properties import FluidProperties


class DLVOForces:
    """
    Calculate DLVO forces between a bubble and a particle.

    The DLVO theory combines:
    1. van der Waals attractive forces
    2. Electrostatic double layer repulsive forces
    """

    # Physical constants
    VACUUM_PERMITTIVITY = 8.854187817e-12  # ε₀ [F/m or C²/(N·m²)]
    ELEMENTARY_CHARGE = 1.602176634e-19  # e [C]

    def __init__(self, particle: Particle, bubble: Bubble, fluid: FluidProperties):
        """
        Initialize DLVO force calculator.

        Args:
            particle: Particle object
            bubble: Bubble object
            fluid: FluidProperties object
        """
        self.particle = particle
        self.bubble = bubble
        self.fluid = fluid

    def van_der_waals_force(self, h: float) -> float:
        """
        Calculate van der Waals force between a bubble and particle.

        For two spheres of radii R₁ and R₂ separated by distance h:

        F_vdW = -A * R₁ * R₂ / (6 * h² * (R₁ + R₂))

        Where:
        - A is the Hamaker constant [J]
        - h is the surface-to-surface separation [m]
        - Negative sign indicates attractive force

        Args:
            h: Surface-to-surface separation distance [m]

        Returns:
            van der Waals force [N] (negative = attractive)

        Note:
            This is the Derjaguin approximation valid for h << R₁, R₂
        """
        if h <= 0:
            # Prevent division by zero; at contact, use a minimum separation
            h = 1e-10  # 0.1 nm minimum separation

        R1 = self.particle.radius
        R2 = self.bubble.radius
        A = self.particle.hamaker_constant

        # Derjaguin approximation for sphere-sphere interaction
        F_vdW = -A * R1 * R2 / (6.0 * h**2 * (R1 + R2))

        return F_vdW

    def electrostatic_force(self, h: float) -> float:
        """
        Calculate electrostatic double layer force between bubble and particle.

        For two spheres with constant surface potentials ψ₁ and ψ₂:

        F_EDL = 2π * ε₀ * εᵣ * κ * R₁ * R₂ / (R₁ + R₂) *
                [(ψ₁² + ψ₂²) * exp(-κh) - 2ψ₁ψ₂ * exp(-2κh)]

        For the constant potential boundary condition (simplified form):
        F_EDL ≈ 2π * ε₀ * εᵣ * κ * R₁ * R₂ / (R₁ + R₂) *
                (ψ₁² + ψ₂²) * exp(-κh)

        Where:
        - ε₀ is vacuum permittivity [F/m]
        - εᵣ is relative dielectric constant
        - κ is the Debye parameter (inverse Debye length) [1/m]
        - ψ₁, ψ₂ are the zeta potentials [V]
        - h is the surface-to-surface separation [m]
        - Positive sign indicates repulsive force (for like charges)

        Args:
            h: Surface-to-surface separation distance [m]

        Returns:
            Electrostatic force [N] (positive = repulsive)

        Note:
            This uses the linear superposition approximation (Derjaguin)
            and assumes zeta potential ≈ surface potential
        """
        if h <= 0:
            h = 1e-10  # Minimum separation

        R1 = self.particle.radius
        R2 = self.bubble.radius
        epsilon_0 = self.VACUUM_PERMITTIVITY
        epsilon_r = self.fluid.dielectric_constant
        kappa = self.fluid.debye_parameter
        psi_1 = self.particle.zeta_potential
        psi_2 = self.bubble.zeta_potential

        # Geometric factor
        geom_factor = (2.0 * math.pi * epsilon_0 * epsilon_r * kappa *
                      R1 * R2 / (R1 + R2))

        # For constant potential:
        # Full expression with cross-term
        term1 = (psi_1**2 + psi_2**2) * math.exp(-kappa * h)
        term2 = -2.0 * psi_1 * psi_2 * math.exp(-2.0 * kappa * h)

        F_EDL = geom_factor * (term1 + term2)

        return F_EDL

    def total_dlvo_force(self, h: float) -> float:
        """
        Calculate total DLVO force (van der Waals + electrostatic).

        F_DLVO = F_vdW + F_EDL

        Args:
            h: Surface-to-surface separation distance [m]

        Returns:
            Total DLVO force [N]
            (negative = net attractive, positive = net repulsive)
        """
        F_vdW = self.van_der_waals_force(h)
        F_EDL = self.electrostatic_force(h)

        return F_vdW + F_EDL

    def interaction_energy(self, h: float) -> Tuple[float, float, float]:
        """
        Calculate DLVO interaction energies (integral of forces).

        V_vdW = -A * R₁ * R₂ / (6 * h * (R₁ + R₂))

        V_EDL = 2π * ε₀ * εᵣ * R₁ * R₂ / (R₁ + R₂) *
                [(ψ₁² + ψ₂²) * (1/κ) * exp(-κh) + 2ψ₁ψ₂ * (1/2κ) * exp(-2κh)]

        Args:
            h: Surface-to-surface separation distance [m]

        Returns:
            Tuple of (V_vdW, V_EDL, V_total) [J]
        """
        if h <= 0:
            h = 1e-10

        R1 = self.particle.radius
        R2 = self.bubble.radius
        A = self.particle.hamaker_constant
        epsilon_0 = self.VACUUM_PERMITTIVITY
        epsilon_r = self.fluid.dielectric_constant
        kappa = self.fluid.debye_parameter
        psi_1 = self.particle.zeta_potential
        psi_2 = self.bubble.zeta_potential

        # van der Waals energy
        V_vdW = -A * R1 * R2 / (6.0 * h * (R1 + R2))

        # Electrostatic energy
        geom_factor = (2.0 * math.pi * epsilon_0 * epsilon_r *
                      R1 * R2 / (R1 + R2))
        term1 = (psi_1**2 + psi_2**2) * (1.0 / kappa) * math.exp(-kappa * h)
        term2 = 2.0 * psi_1 * psi_2 * (1.0 / (2.0 * kappa)) * math.exp(-2.0 * kappa * h)
        V_EDL = geom_factor * (term1 + term2)

        V_total = V_vdW + V_EDL

        return V_vdW, V_EDL, V_total

    def is_favorable_for_attachment(self, h_range: Tuple[float, float] = (1e-9, 100e-9)) -> bool:
        """
        Determine if attachment is energetically favorable.

        Checks if there's an energy barrier in the specified range.

        Args:
            h_range: Range of separations to check [m]

        Returns:
            True if attachment is favorable (no significant barrier)
        """
        h_min, h_max = h_range
        n_points = 100
        h_values = [h_min + i * (h_max - h_min) / n_points for i in range(n_points)]

        energies = [self.interaction_energy(h)[2] for h in h_values]

        # Check if there's a significant energy barrier
        max_energy = max(energies)
        kT = self.fluid.boltzmann_temperature

        # If maximum energy barrier > 10 kT, attachment is unfavorable
        return max_energy < 10.0 * kT
