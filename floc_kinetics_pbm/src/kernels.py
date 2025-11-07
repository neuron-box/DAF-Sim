"""
Aggregation and breakage kernels for Population Balance Model.

This module implements the standard kernels for floc aggregation and breakage
based on the literature (Smoluchowski, Saffman-Turner, etc.), which form the
foundation of CFD-PBM coupling in flocculation modeling.

Mathematical Formulations (LaTeX):

1. Population Balance Equation (PBE):
   $$\\frac{\\partial n(v,t)}{\\partial t} + \\nabla \\cdot [\\mathbf{u}(\\mathbf{x},t) n(v,t)] =
   B_{agg} - D_{agg} + B_{break} - D_{break}$$

   Where:
   - Aggregation birth: $B_{agg} = \\frac{1}{2} \\int_0^v \\beta(v-\\lambda, \\lambda) n(v-\\lambda) n(\\lambda) d\\lambda$
   - Aggregation death: $D_{agg} = n(v) \\int_0^{\\infty} \\beta(v, \\lambda) n(\\lambda) d\\lambda$
   - Breakage birth: $B_{break} = \\int_v^{\\infty} S(\\lambda) \\Gamma(v|\\lambda) n(\\lambda) d\\lambda$
   - Breakage death: $D_{break} = S(v) n(v)$

2. Aggregation Kernel β(v_i, v_j):
   $$\\beta(v_i, v_j) = \\alpha [\\beta_{ortho}(v_i, v_j) + \\beta_{peri}(v_i, v_j) + \\beta_{ds}(v_i, v_j)]$$

   a) Orthokinetic (Turbulent Shear) - Saffman-Turner:
   $$\\beta_{ortho} = 1.3 \\sqrt{\\frac{\\varepsilon}{\\nu}} (d_i + d_j)^3$$

   b) Perikinetic (Brownian Motion) - Smoluchowski:
   $$\\beta_{peri} = \\frac{2 k_B T}{3 \\mu} \\frac{(d_i + d_j)^2}{d_i d_j}$$

   c) Differential Sedimentation:
   $$\\beta_{ds} = \\frac{\\pi g}{72 \\mu} |\\rho_i - \\rho_w| (d_i + d_j)^3 |d_i - d_j|$$

3. Breakage Kernel S(v):
   $$S(v_i) = k_b \\left(\\frac{\\varepsilon d_i^2}{\\sigma_f}\\right)^n \\exp\\left(-\\frac{C \\sigma_f}{\\rho_f \\varepsilon^{2/3} d_i^{2/3}}\\right)$$

   Where:
   - G: Velocity gradient (shear rate) [1/s]
   - ε: Turbulent energy dissipation rate [m²/s³]
   - ν: Kinematic viscosity [m²/s]
   - α: Collision efficiency factor [-]
   - σ_f: Floc binding strength [N/m²]
   - k_b, n, C: Empirical constants
"""

import numpy as np
from typing import Tuple, Union, Optional
from .properties import FlocProperties


class FlocKernels:
    """
    Calculate aggregation and breakage kernels for floc kinetics.

    This class implements the standard formulations used in CFD-PBM coupling
    for flocculation processes.
    """

    def __init__(self, properties: Optional[FlocProperties] = None):
        """
        Initialize kernel calculator with floc properties.

        Args:
            properties: FlocProperties object containing physical parameters.
                       If None, default properties are used.
        """
        self.props = properties if properties is not None else FlocProperties()

        # Breakage kernel empirical constants (typical values from literature)
        self.k_b = 1e-2  # Breakage rate constant [1/s]
        self.n_break = 1.0  # Breakage exponent [-]
        self.C_break = 0.1  # Breakage constant [-]

    def diameter_from_volume(self, volume: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Convert particle volume to equivalent spherical diameter.

        Args:
            volume: Particle volume [m³]

        Returns:
            Equivalent spherical diameter [m]
        """
        return (6.0 * volume / np.pi)**(1.0 / 3.0)

    def volume_from_diameter(self, diameter: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Convert particle diameter to volume.

        Args:
            diameter: Particle diameter [m]

        Returns:
            Particle volume [m³]
        """
        return np.pi / 6.0 * diameter**3

    def shear_rate_from_dissipation(self, epsilon: float) -> float:
        """
        Calculate shear rate from turbulent energy dissipation rate.

        Args:
            epsilon: Turbulent energy dissipation rate [m²/s³]

        Returns:
            Shear rate G [1/s]
        """
        return np.sqrt(epsilon / self.props.nu_water)

    def dissipation_from_shear_rate(self, G: float) -> float:
        """
        Calculate turbulent energy dissipation rate from shear rate.

        Args:
            G: Shear rate [1/s]

        Returns:
            Turbulent energy dissipation rate ε [m²/s³]
        """
        return G**2 * self.props.nu_water

    def beta_orthokinetic(self, d_i: float, d_j: float, G: float) -> float:
        """
        Calculate orthokinetic (turbulent shear) aggregation kernel.

        Based on Saffman-Turner formulation for turbulent collision.

        Args:
            d_i: Diameter of particle i [m]
            d_j: Diameter of particle j [m]
            G: Velocity gradient (shear rate) [1/s]

        Returns:
            Orthokinetic aggregation kernel [m³/s]
        """
        # Saffman-Turner kernel: β_ortho = 1.3 * sqrt(ε/ν) * (d_i + d_j)³
        # Or equivalently: β_ortho = 1.3 * G * (d_i + d_j)³
        return 1.3 * G * (d_i + d_j)**3

    def beta_perikinetic(self, d_i: float, d_j: float) -> float:
        """
        Calculate perikinetic (Brownian motion) aggregation kernel.

        Based on Smoluchowski's formulation for Brownian collisions.

        Args:
            d_i: Diameter of particle i [m]
            d_j: Diameter of particle j [m]

        Returns:
            Perikinetic aggregation kernel [m³/s]
        """
        # Smoluchowski kernel: β_peri = (2 k_B T)/(3 μ) * (d_i + d_j)²/(d_i * d_j)
        k_B = self.props.k_B
        T = self.props.temperature
        mu = self.props.mu_water

        return (2.0 * k_B * T) / (3.0 * mu) * (d_i + d_j)**2 / (d_i * d_j)

    def beta_differential_sedimentation(self, d_i: float, d_j: float) -> float:
        """
        Calculate differential sedimentation aggregation kernel.

        Based on gravity-driven collisions due to settling velocity differences.

        Args:
            d_i: Diameter of particle i [m]
            d_j: Diameter of particle j [m]

        Returns:
            Differential sedimentation aggregation kernel [m³/s]
        """
        # β_ds = (π g)/(72 μ) * |ρ_i - ρ_w| * (d_i + d_j)³ * |d_i - d_j|
        g = self.props.g
        mu = self.props.mu_water
        rho_w = self.props.rho_water

        # Get effective densities based on fractal structure
        rho_i = self.props.floc_density(d_i)
        rho_j = self.props.floc_density(d_j)

        # Average density difference
        delta_rho = 0.5 * (abs(rho_i - rho_w) + abs(rho_j - rho_w))

        return (np.pi * g) / (72.0 * mu) * delta_rho * \
               (d_i + d_j)**3 * abs(d_i - d_j)

    def beta_total(self, d_i: float, d_j: float, G: float,
                   include_perikinetic: bool = True,
                   include_sedimentation: bool = True) -> float:
        """
        Calculate total aggregation kernel.

        β_total = α * (β_ortho + β_peri + β_ds)

        Args:
            d_i: Diameter of particle i [m]
            d_j: Diameter of particle j [m]
            G: Velocity gradient (shear rate) [1/s]
            include_perikinetic: Include Brownian motion term
            include_sedimentation: Include differential sedimentation term

        Returns:
            Total aggregation kernel [m³/s]
        """
        alpha = self.props.collision_efficiency

        # Orthokinetic term (always included)
        beta = self.beta_orthokinetic(d_i, d_j, G)

        # Perikinetic term (important for small particles < 1 μm)
        if include_perikinetic:
            beta += self.beta_perikinetic(d_i, d_j)

        # Differential sedimentation (important for large particles > 40 μm)
        if include_sedimentation:
            beta += self.beta_differential_sedimentation(d_i, d_j)

        return alpha * beta

    def breakage_rate(self, d_i: float, G: float) -> float:
        """
        Calculate breakage rate (breakage kernel).

        Based on the balance between turbulent stress and floc strength.

        S(d_i) = k_b * (ε * d_i² / σ_f)^n * exp(-C * σ_f / (ρ_f * ε^(2/3) * d_i^(2/3)))

        Args:
            d_i: Diameter of particle i [m]
            G: Velocity gradient (shear rate) [1/s]

        Returns:
            Breakage rate [1/s]
        """
        # Convert shear rate to dissipation rate
        epsilon = self.dissipation_from_shear_rate(G)

        sigma_f = self.props.binding_strength
        rho_f = self.props.floc_density(d_i)

        # Turbulent stress factor
        stress_factor = (epsilon * d_i**2) / sigma_f

        # Exponential term (breakage resistance)
        exp_term = np.exp(-self.C_break * sigma_f /
                          (rho_f * epsilon**(2.0/3.0) * d_i**(2.0/3.0)))

        # Total breakage rate
        S = self.k_b * stress_factor**self.n_break * exp_term

        return S

    def daughter_distribution(self, v_parent: float, v_daughter: float) -> float:
        """
        Calculate daughter particle distribution from breakage.

        Γ(v_daughter | v_parent) - probability density that a particle of volume
        v_parent breaks into a daughter particle of volume v_daughter.

        Simple binary breakage model: two equal-sized daughters.

        Args:
            v_parent: Volume of parent particle [m³]
            v_daughter: Volume of daughter particle [m³]

        Returns:
            Daughter distribution function [1/m³]
        """
        # Binary breakage: two particles of volume v_parent/2
        if v_daughter <= v_parent:
            # Normalized distribution (two equal daughters)
            return 2.0 / v_parent
        else:
            return 0.0


def calculate_floc_kernels(
    G: float,
    particle_size_i: float,
    particle_size_j: Optional[float] = None,
    properties: Optional[FlocProperties] = None,
    size_units: str = 'diameter'
) -> Tuple[float, float]:
    """
    Calculate aggregation and breakage rates for floc kinetics.

    This is the main function interface for computing kernels as specified
    in the task requirements.

    Args:
        G: Local velocity gradient (shear rate) [1/s]
        particle_size_i: Size of particle i [m] (diameter or volume)
        particle_size_j: Size of particle j [m] (diameter or volume).
                        If None, only breakage rate is calculated.
        properties: FlocProperties object with physical parameters.
                   If None, default properties are used.
        size_units: 'diameter' or 'volume' to specify input size units

    Returns:
        Tuple of (aggregation_kernel, breakage_rate):
            - aggregation_kernel: β(i,j) [m³/s] (0 if particle_j is None)
            - breakage_rate: S(i) [1/s]

    Example:
        >>> props = FlocProperties(collision_efficiency=0.1)
        >>> G = 100.0  # 1/s (typical for flocculation)
        >>> d_i = 50e-6  # 50 μm
        >>> d_j = 30e-6  # 30 μm
        >>> beta, S = calculate_floc_kernels(G, d_i, d_j, props)
        >>> print(f"Aggregation kernel: {beta:.3e} m³/s")
        >>> print(f"Breakage rate: {S:.3e} 1/s")
    """
    if G < 0:
        raise ValueError("Shear rate G must be non-negative")
    if particle_size_i <= 0:
        raise ValueError("Particle size must be positive")
    if particle_size_j is not None and particle_size_j <= 0:
        raise ValueError("Particle size must be positive")

    # Initialize kernel calculator
    kernels = FlocKernels(properties)

    # Convert to diameter if needed
    if size_units == 'volume':
        d_i = kernels.diameter_from_volume(particle_size_i)
        d_j = kernels.diameter_from_volume(particle_size_j) if particle_size_j is not None else None
    elif size_units == 'diameter':
        d_i = particle_size_i
        d_j = particle_size_j
    else:
        raise ValueError("size_units must be 'diameter' or 'volume'")

    # Calculate breakage rate
    breakage_rate = kernels.breakage_rate(d_i, G)

    # Calculate aggregation kernel if second particle is provided
    if d_j is not None:
        aggregation_kernel = kernels.beta_total(d_i, d_j, G)
    else:
        aggregation_kernel = 0.0

    return aggregation_kernel, breakage_rate
