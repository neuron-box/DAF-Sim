"""
Hydrodynamic Forces for Bubble-Particle Interactions

This module implements hydrodynamic forces acting on particles approaching bubbles,
including drag force in the bubble's flow field.

References:
    - Han et al. (2001). Water Science and Technology, 43(8), 139-144.
    - Yoon, R. H., & Luttrell, G. H. (1989). Mineral Processing and Extractive Metallurgy Review, 5(1-4), 101-122.
"""

import math
from typing import Tuple
from .particle import Particle
from .bubble import Bubble
from .fluid_properties import FluidProperties


class HydrodynamicForces:
    """
    Calculate hydrodynamic forces on a particle in the flow field around a rising bubble.
    """

    def __init__(self, particle: Particle, bubble: Bubble, fluid: FluidProperties):
        """
        Initialize hydrodynamic force calculator.

        Args:
            particle: Particle object
            bubble: Bubble object
            fluid: FluidProperties object
        """
        self.particle = particle
        self.bubble = bubble
        self.fluid = fluid

        # Calculate bubble rise velocity if not provided
        if self.bubble.rise_velocity is None:
            self.bubble.rise_velocity = self.bubble.calculate_rise_velocity(
                fluid.dynamic_viscosity, fluid.density
            )

    def stokes_drag_force(self, velocity: float) -> float:
        """
        Calculate Stokes drag force on a particle.

        For a spherical particle in Stokes flow:
        F_d = 6π * μ * R * v

        Where:
        - μ is the dynamic viscosity [Pa·s]
        - R is the particle radius [m]
        - v is the relative velocity [m/s]

        Args:
            velocity: Relative velocity between particle and fluid [m/s]

        Returns:
            Drag force magnitude [N]
        """
        mu = self.fluid.dynamic_viscosity
        R = self.particle.radius

        F_drag = 6.0 * math.pi * mu * R * abs(velocity)

        # Return signed force (opposite to velocity direction)
        return -math.copysign(F_drag, velocity)

    def bubble_flow_velocity(self, r: float, theta: float) -> Tuple[float, float]:
        """
        Calculate the flow velocity field around a rising bubble.

        For potential flow around a sphere (bubble), the velocity field is:

        v_r = U * (1 - (R_b/r)³) * cos(θ)
        v_θ = -U * (1 + 0.5 * (R_b/r)³) * sin(θ)

        Where:
        - U is the bubble rise velocity
        - R_b is the bubble radius
        - r is the radial distance from bubble center
        - θ is the angle from the vertical axis (0 = top, π = bottom)

        Note: This uses potential flow with the correction factor for a bubble
        with partially mobile interface (factor 0.5 vs 1.0 for rigid sphere).

        Args:
            r: Radial distance from bubble center [m]
            theta: Angle from vertical axis [rad]

        Returns:
            Tuple of (v_r, v_theta) in [m/s]
        """
        U = self.bubble.rise_velocity
        R_b = self.bubble.radius

        if r < R_b:
            # Inside bubble - should not happen, but return zero
            return 0.0, 0.0

        # Radial velocity component
        v_r = U * (1.0 - (R_b / r)**3) * math.cos(theta)

        # Tangential velocity component
        # Factor 0.5 for mobile interface (vs 1.0 for rigid sphere)
        v_theta = -U * (1.0 + 0.5 * (R_b / r)**3) * math.sin(theta)

        return v_r, v_theta

    def relative_velocity(self, r: float, theta: float, v_particle: Tuple[float, float]) -> Tuple[float, float]:
        """
        Calculate relative velocity between fluid and particle.

        Args:
            r: Radial distance from bubble center [m]
            theta: Angle from vertical axis [rad]
            v_particle: Particle velocity (v_r, v_theta) [m/s]

        Returns:
            Relative velocity (v_rel_r, v_rel_theta) [m/s]
        """
        v_fluid_r, v_fluid_theta = self.bubble_flow_velocity(r, theta)
        v_p_r, v_p_theta = v_particle

        v_rel_r = v_fluid_r - v_p_r
        v_rel_theta = v_fluid_theta - v_p_theta

        return v_rel_r, v_rel_theta

    def drag_force_spherical(self, r: float, theta: float,
                           v_particle: Tuple[float, float]) -> Tuple[float, float]:
        """
        Calculate drag force on particle in spherical coordinates.

        F_drag = 6π * μ * R_p * v_rel

        Args:
            r: Radial distance from bubble center [m]
            theta: Angle from vertical axis [rad]
            v_particle: Particle velocity (v_r, v_theta) [m/s]

        Returns:
            Drag force (F_r, F_theta) [N]
        """
        v_rel_r, v_rel_theta = self.relative_velocity(r, theta, v_particle)

        mu = self.fluid.dynamic_viscosity
        R_p = self.particle.radius

        # Stokes drag coefficient
        drag_coeff = 6.0 * math.pi * mu * R_p

        F_r = drag_coeff * v_rel_r
        F_theta = drag_coeff * v_rel_theta

        return F_r, F_theta

    def gravitational_force(self) -> float:
        """
        Calculate gravitational force on particle (buoyancy corrected).

        F_g = (ρ_p - ρ_f) * V_p * g

        Where:
        - ρ_p is particle density
        - ρ_f is fluid density
        - V_p is particle volume
        - g is gravitational acceleration

        Returns:
            Gravitational force [N] (positive = downward)
        """
        g = 9.81  # [m/s²]
        rho_p = self.particle.density
        rho_f = self.fluid.density
        V_p = self.particle.volume

        F_g = (rho_p - rho_f) * V_p * g

        return F_g

    def added_mass_coefficient(self) -> float:
        """
        Calculate the added mass coefficient for a sphere.

        For a sphere in fluid, the added mass coefficient is 0.5.

        Returns:
            Added mass coefficient (dimensionless)
        """
        return 0.5

    def reynolds_number(self, velocity: float) -> float:
        """
        Calculate particle Reynolds number.

        Re = ρ * d * v / μ

        Args:
            velocity: Characteristic velocity [m/s]

        Returns:
            Reynolds number (dimensionless)
        """
        rho = self.fluid.density
        d = self.particle.diameter
        mu = self.fluid.dynamic_viscosity

        Re = rho * d * abs(velocity) / mu

        return Re

    def is_stokes_regime(self, velocity: float) -> bool:
        """
        Check if particle is in Stokes flow regime (Re << 1).

        Args:
            velocity: Characteristic velocity [m/s]

        Returns:
            True if Re < 0.1 (Stokes regime)
        """
        Re = self.reynolds_number(velocity)
        return Re < 0.1
