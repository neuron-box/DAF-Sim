"""
Physical properties for floc kinetics modeling.

This module defines the physical properties of flocs, water, and related
parameters used in the Population Balance Model.
"""

import numpy as np
from dataclasses import dataclass
from typing import Optional


@dataclass
class FlocProperties:
    """
    Physical properties of flocs and the surrounding fluid.

    Attributes:
        rho_water: Water density [kg/m³]
        rho_floc: Floc density [kg/m³]
        mu_water: Dynamic viscosity of water [Pa·s]
        k_B: Boltzmann constant [J/K]
        temperature: Temperature [K]
        fractal_dimension: Fractal dimension of flocs [-]
        collision_efficiency: Collision efficiency factor (0-1) [-]
        binding_strength: Floc binding strength [N/m²]
        primary_particle_diameter: Diameter of primary particles [m]
    """

    # Fluid properties
    rho_water: float = 998.0  # kg/m³ at 20°C
    mu_water: float = 1.002e-3  # Pa·s at 20°C
    temperature: float = 293.15  # K (20°C)

    # Floc properties
    rho_floc: float = 1050.0  # kg/m³ (typical for wastewater flocs)
    fractal_dimension: float = 2.3  # Typical fractal dimension for flocs
    primary_particle_diameter: float = 1e-6  # m (1 µm)

    # Interaction parameters
    collision_efficiency: float = 0.1  # Sticking probability (0-1)
    binding_strength: float = 1e3  # N/m² (floc cohesive strength)

    # Physical constants
    k_B: float = 1.380649e-23  # Boltzmann constant [J/K]
    g: float = 9.81  # Gravitational acceleration [m/s²]

    def __post_init__(self):
        """Validate property values."""
        if self.rho_water <= 0:
            raise ValueError("Water density must be positive")
        if self.mu_water <= 0:
            raise ValueError("Dynamic viscosity must be positive")
        if self.temperature <= 0:
            raise ValueError("Temperature must be positive")
        if not 0 <= self.collision_efficiency <= 1:
            raise ValueError("Collision efficiency must be between 0 and 1")
        if not 1 <= self.fractal_dimension <= 3:
            raise ValueError("Fractal dimension must be between 1 and 3")

    @property
    def nu_water(self) -> float:
        """Kinematic viscosity [m²/s]."""
        return self.mu_water / self.rho_water

    def floc_density(self, diameter: float) -> float:
        """
        Calculate effective floc density based on fractal dimension.

        Larger flocs have lower effective density due to fractal structure.

        Args:
            diameter: Floc diameter [m]

        Returns:
            Effective floc density [kg/m³]
        """
        d_ratio = diameter / self.primary_particle_diameter

        # Fractal relationship: ρ_eff = ρ_water + (ρ_primary - ρ_water) * (d/d_0)^(Df-3)
        rho_eff = self.rho_water + (self.rho_floc - self.rho_water) * \
                  d_ratio**(self.fractal_dimension - 3.0)

        return np.clip(rho_eff, self.rho_water, self.rho_floc)

    def floc_mass(self, diameter: float) -> float:
        """
        Calculate floc mass based on fractal structure.

        Args:
            diameter: Floc diameter [m]

        Returns:
            Floc mass [kg]
        """
        volume = np.pi / 6.0 * diameter**3
        return self.floc_density(diameter) * volume

    def kolmogorov_length_scale(self, epsilon: float) -> float:
        """
        Calculate Kolmogorov length scale.

        Args:
            epsilon: Turbulent energy dissipation rate [m²/s³]

        Returns:
            Kolmogorov length scale [m]
        """
        return (self.nu_water**3 / epsilon)**0.25

    def kolmogorov_time_scale(self, epsilon: float) -> float:
        """
        Calculate Kolmogorov time scale.

        Args:
            epsilon: Turbulent energy dissipation rate [m²/s³]

        Returns:
            Kolmogorov time scale [s]
        """
        return (self.nu_water / epsilon)**0.5
