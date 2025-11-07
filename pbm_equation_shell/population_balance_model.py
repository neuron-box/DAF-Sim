"""
Population Balance Model (PBM) Equation Shell

This module implements the discretized Population Balance Equation (PBE)
for flocculation modeling based on Abreu et al. (2021).

Reference:
    Abreu, D. A. d. S., et al. (2021). "Integration of first principle models
    and machine learning in a modeling framework: An application to flocculation"
    Chemical Engineering Science.

The discretized PBE (Equation 2) describes the time evolution of particle
number concentrations through aggregation and breakage mechanisms.

This implementation supports both linear and geometric particle volume grids
using the fixed pivot technique (Kumar & Ramkrishna, 1996) for geometric grids.
"""

import numpy as np
from typing import Optional, Literal
import warnings


class PopulationBalanceModel:
    """
    Population Balance Model for flocculation processes.

    This class implements the discretized Population Balance Equation (PBE)
    that describes particle aggregation and breakage in a flocculation system.

    Supports both linear and geometric particle volume grids with proper
    mass conservation.

    Mathematical Formulation:
    -------------------------
    The discretized PBE for each size class i is:

    dN_i/dt = (Birth by aggregation) - (Death by aggregation)
              + (Birth by breakage) - (Death by breakage)

    For LINEAR grids (v_i = i·Δv):

    dN_i/dt = (1/2) * Σ_{j=1}^{i-1} β_{j,i-j} * α_{j,i-j} * N_j * N_{i-j}
              - N_i * Σ_{j=1}^{n} β_{i,j} * α_{i,j} * N_j
              + Σ_{j=i+1}^{n} S_j * γ_{j,i} * N_j
              - S_i * N_i

    For GEOMETRIC grids (v_i = v_0·r^i):
    Uses fixed pivot technique to distribute mass between adjacent bins.

    Attributes:
    -----------
    n_classes : int
        Number of discrete size classes
    particle_volumes : np.ndarray
        Volume of particles in each size class [m³]
    grid_type : str
        Type of volume grid: 'linear' or 'geometric'
    grid_ratio : float or None
        Ratio r for geometric grids (v_i = v_0 * r^i)
    validate_inputs : bool
        Whether to validate kernel inputs for consistency
    """

    def __init__(
        self,
        n_classes: int,
        particle_volumes: np.ndarray,
        validate_inputs: bool = True,
        grid_type: Optional[Literal['linear', 'geometric']] = None
    ):
        """
        Initialize the Population Balance Model.

        Parameters:
        -----------
        n_classes : int
            Number of discrete size classes (bins) for particle size distribution
        particle_volumes : np.ndarray
            Volume of particles in each size class [m³]
            Shape: (n_classes,)
            For geometric grids: v_i = v_0 * r^i
            For linear grids: v_i = v_0 + i * Δv
        validate_inputs : bool, optional
            Whether to validate kernel matrices/vectors for physical consistency
            Default is True
        grid_type : {'linear', 'geometric'}, optional
            Type of volume grid. If None, auto-detected from particle_volumes.
            Default is None (auto-detect)

        Raises:
        -------
        ValueError
            If n_classes < 2 or particle_volumes has wrong shape
        """
        if n_classes < 2:
            raise ValueError("Number of size classes must be at least 2")

        if particle_volumes.shape != (n_classes,):
            raise ValueError(
                f"particle_volumes must have shape ({n_classes},), "
                f"got {particle_volumes.shape}"
            )

        if np.any(particle_volumes <= 0):
            raise ValueError("particle_volumes must contain positive values")

        if not np.all(np.diff(particle_volumes) > 0):
            raise ValueError("particle_volumes must be strictly increasing")

        self.n_classes = n_classes
        self.particle_volumes = particle_volumes.copy()
        self.validate_inputs = validate_inputs

        # Detect or validate grid type
        if grid_type is None:
            self.grid_type = self._detect_grid_type()
        else:
            if grid_type not in ['linear', 'geometric']:
                raise ValueError("grid_type must be 'linear' or 'geometric'")
            self.grid_type = grid_type

        # For geometric grids, calculate and store the grid ratio
        if self.grid_type == 'geometric':
            self.grid_ratio = self._calculate_grid_ratio()
        else:
            self.grid_ratio = None

    def _detect_grid_type(self) -> str:
        """
        Auto-detect whether the grid is linear or geometric.

        Returns:
        --------
        grid_type : str
            'linear' or 'geometric'
        """
        # Check if ratios are approximately constant (geometric grid)
        ratios = self.particle_volumes[1:] / self.particle_volumes[:-1]

        # If ratios are approximately constant, it's geometric
        if np.std(ratios) / np.mean(ratios) < 0.01:  # 1% coefficient of variation
            return 'geometric'

        # Check if differences are approximately constant (linear grid)
        diffs = np.diff(self.particle_volumes)

        if np.std(diffs) / np.mean(diffs) < 0.01:
            return 'linear'

        # Default to geometric (more common in practice)
        warnings.warn(
            "Grid type could not be definitively determined. "
            "Defaulting to 'geometric'. Specify grid_type explicitly to avoid ambiguity."
        )
        return 'geometric'

    def _calculate_grid_ratio(self) -> float:
        """
        Calculate the grid ratio r for geometric grids.

        Returns:
        --------
        r : float
            Grid ratio (v_{i+1} / v_i)
        """
        ratios = self.particle_volumes[1:] / self.particle_volumes[:-1]
        return np.mean(ratios)

    def _find_target_bin(self, v_new: float) -> int:
        """
        Find the bin index where v_new falls.

        Uses binary search for efficiency.

        Parameters:
        -----------
        v_new : float
            Volume of newly formed particle [m³]

        Returns:
        --------
        i : int
            Bin index where v_i <= v_new < v_{i+1}
            Returns n_classes-1 if v_new >= v_max (overflow)
            Returns -1 if v_new < v_min (underflow, shouldn't happen)
        """
        # Handle edge cases
        if v_new >= self.particle_volumes[-1]:
            return self.n_classes - 1  # Overflow to largest bin

        if v_new < self.particle_volumes[0]:
            return -1  # Underflow (shouldn't happen in aggregation)

        # Binary search
        i = np.searchsorted(self.particle_volumes, v_new, side='right') - 1
        return i

    def _calculate_distribution_factors(self, v_new: float, i: int) -> tuple[float, float]:
        """
        Calculate distribution factors ξ and η for geometric grids.

        Distributes mass between bins i and i+1 to conserve number and volume.

        Parameters:
        -----------
        v_new : float
            Volume of newly formed particle [m³]
        i : int
            Lower bin index (v_i <= v_new < v_{i+1})

        Returns:
        --------
        xi : float
            Fraction going to bin i
        eta : float
            Fraction going to bin i+1

        Notes:
        ------
        ξ + η = 1 (number conservation)
        ξ*v_i + η*v_{i+1} = v_new (volume conservation)
        """
        # Handle overflow to largest bin
        if i >= self.n_classes - 1:
            return 1.0, 0.0

        # Handle underflow (shouldn't happen)
        if i < 0:
            return 0.0, 0.0

        v_i = self.particle_volumes[i]
        v_ip1 = self.particle_volumes[i + 1]

        # Distribution factors
        eta = (v_new - v_i) / (v_ip1 - v_i)
        xi = 1.0 - eta

        return xi, eta

    def calculate_dndt(
        self,
        N: np.ndarray,
        alpha_matrix: np.ndarray,
        beta_matrix: np.ndarray,
        S_vector: np.ndarray,
        gamma_matrix: np.ndarray
    ) -> np.ndarray:
        """
        Calculate the time derivative of particle number concentrations.

        This method implements the discretized Population Balance Equation (Equation 2
        from Abreu et al., 2021) to compute dN/dt for all size classes.

        For geometric grids, uses the fixed pivot technique (Kumar & Ramkrishna, 1996)
        to properly conserve mass during aggregation.

        Parameters:
        -----------
        N : np.ndarray
            Current particle number concentration vector [#/m³]
            Shape: (n_classes,)
            N[i] represents the number concentration in size class i

        alpha_matrix : np.ndarray
            Collision efficiency matrix [-]
            Shape: (n_classes, n_classes)
            alpha_matrix[i,j] is the collision efficiency between size classes i and j
            Dimensionless, typically ranges from 0 to 1

        beta_matrix : np.ndarray
            Collision frequency matrix [m³/s]
            Shape: (n_classes, n_classes)
            beta_matrix[i,j] is the collision frequency between size classes i and j
            Typically constructed from perikinetic, differential sedimentation,
            and orthokinetic mechanisms

        S_vector : np.ndarray
            Breakage rate vector [1/s]
            Shape: (n_classes,)
            S_vector[i] is the breakage rate for size class i

        gamma_matrix : np.ndarray
            Breakage distribution function matrix [-]
            Shape: (n_classes, n_classes)
            gamma_matrix[j,i] is the probability that breakage of particle j
            produces particle i
            Each row should sum to 1: Σ_i γ_{j,i} = 1

        Returns:
        --------
        dndt : np.ndarray
            Time derivative vector [#/(m³·s)]
            Shape: (n_classes,)
            dndt[i] = dN_i/dt for each size class i

        Raises:
        -------
        ValueError
            If input dimensions are inconsistent or values are invalid

        Notes:
        ------
        The equation consists of four terms:
        1. Birth by aggregation: Formation of size i from smaller particles
        2. Death by aggregation: Loss of size i through collision with other particles
        3. Birth by breakage: Formation of size i from breakage of larger particles
        4. Death by breakage: Loss of size i through particle breakage
        """
        # Input validation
        self._validate_inputs(N, alpha_matrix, beta_matrix, S_vector, gamma_matrix)

        # Initialize the derivative vector
        dndt = np.zeros(self.n_classes, dtype=np.float64)

        # Choose aggregation method based on grid type
        if self.grid_type == 'linear':
            birth_aggregation = self._calculate_birth_aggregation_linear(
                N, alpha_matrix, beta_matrix
            )
        else:  # geometric
            birth_aggregation = self._calculate_birth_aggregation_geometric(
                N, alpha_matrix, beta_matrix
            )

        # Add birth by aggregation
        dndt += birth_aggregation

        # Loop through each size class for remaining terms
        for i in range(self.n_classes):
            # Term 2: Death by aggregation
            # -N_i * Σ_{j=1}^{n} β_{i,j} * α_{i,j} * N_j
            # Particles of size i are lost through collision with any other particle j
            death_aggregation = 0.0
            for j in range(self.n_classes):
                death_aggregation += beta_matrix[i, j] * alpha_matrix[i, j] * N[j]
            dndt[i] -= death_aggregation * N[i]

            # Term 3: Birth by breakage
            # Σ_{j=i+1}^{n} S_j * γ_{j,i} * N_j
            # Particles of size i are formed by breakage of larger particles j
            birth_breakage = 0.0
            for j in range(i + 1, self.n_classes):
                birth_breakage += S_vector[j] * gamma_matrix[j, i] * N[j]
            dndt[i] += birth_breakage

            # Term 4: Death by breakage
            # -S_i * N_i
            # Particles of size i are lost through their own breakage
            dndt[i] -= S_vector[i] * N[i]

        return dndt

    def _calculate_birth_aggregation_linear(
        self,
        N: np.ndarray,
        alpha_matrix: np.ndarray,
        beta_matrix: np.ndarray
    ) -> np.ndarray:
        """
        Calculate birth by aggregation for LINEAR volume grids.

        For linear grids: v_i = v_0 + i*Δv, so v_j + v_k = v_{j+k}

        Parameters:
        -----------
        N : np.ndarray
            Particle number concentrations [#/m³]
        alpha_matrix : np.ndarray
            Collision efficiency matrix [-]
        beta_matrix : np.ndarray
            Collision frequency matrix [m³/s]

        Returns:
        --------
        birth : np.ndarray
            Birth rate by aggregation for each size class [#/(m³·s)]
        """
        birth = np.zeros(self.n_classes, dtype=np.float64)

        for i in range(self.n_classes):
            # (1/2) * Σ_{j=1}^{i-1} β_{j,i-j} * α_{j,i-j} * N_j * N_{i-j}
            # Particles of size i are formed by collision of smaller particles j and (i-j)
            for j in range(i):
                k = i - j - 1  # Corresponding size class (0-indexed)
                if k >= 0 and k < self.n_classes:
                    birth[i] += (
                        beta_matrix[j, k] * alpha_matrix[j, k] * N[j] * N[k]
                    )
            birth[i] *= 0.5  # Factor of 1/2 to avoid double counting

        return birth

    def _calculate_birth_aggregation_geometric(
        self,
        N: np.ndarray,
        alpha_matrix: np.ndarray,
        beta_matrix: np.ndarray
    ) -> np.ndarray:
        """
        Calculate birth by aggregation for GEOMETRIC volume grids.

        Uses fixed pivot technique (Kumar & Ramkrishna, 1996) to distribute
        mass between adjacent bins.

        For geometric grids: v_i = v_0 * r^i, so v_j + v_k doesn't match any grid point.
        We distribute the birth between bins i and i+1 where v_i <= v_j + v_k < v_{i+1}.

        Parameters:
        -----------
        N : np.ndarray
            Particle number concentrations [#/m³]
        alpha_matrix : np.ndarray
            Collision efficiency matrix [-]
        beta_matrix : np.ndarray
            Collision frequency matrix [m³/s]

        Returns:
        --------
        birth : np.ndarray
            Birth rate by aggregation for each size class [#/(m³·s)]
        """
        birth = np.zeros(self.n_classes, dtype=np.float64)

        # Loop over all possible collision pairs
        for j in range(self.n_classes):
            for k in range(j, self.n_classes):  # k >= j to avoid double counting
                # Volume of newly formed particle
                v_new = self.particle_volumes[j] + self.particle_volumes[k]

                # Find target bin
                i = self._find_target_bin(v_new)

                if i < 0:
                    # Underflow (shouldn't happen for aggregation)
                    continue

                # Calculate distribution factors
                xi, eta = self._calculate_distribution_factors(v_new, i)

                # Collision rate
                if j == k:
                    # Same size class: factor of 1/2 to avoid double counting
                    collision_rate = 0.5 * beta_matrix[j, k] * alpha_matrix[j, k] * N[j] * N[k]
                else:
                    # Different size classes: full rate (symmetry already accounts for double counting)
                    collision_rate = beta_matrix[j, k] * alpha_matrix[j, k] * N[j] * N[k]

                # Distribute to bins i and i+1
                birth[i] += xi * collision_rate

                if i + 1 < self.n_classes:
                    birth[i + 1] += eta * collision_rate

        return birth

    def _validate_inputs(
        self,
        N: np.ndarray,
        alpha_matrix: np.ndarray,
        beta_matrix: np.ndarray,
        S_vector: np.ndarray,
        gamma_matrix: np.ndarray
    ) -> None:
        """
        Validate input dimensions and physical consistency.

        Parameters:
        -----------
        N, alpha_matrix, beta_matrix, S_vector, gamma_matrix : np.ndarray
            Input arrays to validate

        Raises:
        -------
        ValueError
            If any validation check fails
        """
        # Check N vector
        if N.shape != (self.n_classes,):
            raise ValueError(
                f"N must have shape ({self.n_classes},), got {N.shape}"
            )
        if np.any(N < 0):
            raise ValueError("N must contain non-negative values")

        # Check alpha matrix
        if alpha_matrix.shape != (self.n_classes, self.n_classes):
            raise ValueError(
                f"alpha_matrix must have shape ({self.n_classes}, {self.n_classes}), "
                f"got {alpha_matrix.shape}"
            )
        if self.validate_inputs:
            if np.any(alpha_matrix < 0) or np.any(alpha_matrix > 1):
                raise ValueError("alpha_matrix must contain values in range [0, 1]")

        # Check beta matrix
        if beta_matrix.shape != (self.n_classes, self.n_classes):
            raise ValueError(
                f"beta_matrix must have shape ({self.n_classes}, {self.n_classes}), "
                f"got {beta_matrix.shape}"
            )
        if np.any(beta_matrix < 0):
            raise ValueError("beta_matrix must contain non-negative values")

        # Check S vector
        if S_vector.shape != (self.n_classes,):
            raise ValueError(
                f"S_vector must have shape ({self.n_classes},), got {S_vector.shape}"
            )
        if np.any(S_vector < 0):
            raise ValueError("S_vector must contain non-negative values")

        # Check gamma matrix
        if gamma_matrix.shape != (self.n_classes, self.n_classes):
            raise ValueError(
                f"gamma_matrix must have shape ({self.n_classes}, {self.n_classes}), "
                f"got {gamma_matrix.shape}"
            )
        if np.any(gamma_matrix < 0):
            raise ValueError("gamma_matrix must contain non-negative values")

    def validate_mass_conservation(
        self,
        N: np.ndarray,
        dndt: np.ndarray,
        tolerance: float = 1e-10
    ) -> tuple[bool, float]:
        """
        Check mass conservation in the PBM system.

        For a system with only aggregation (no breakage), the total mass should
        be conserved. This method checks if the rate of change of total mass
        is approximately zero.

        Parameters:
        -----------
        N : np.ndarray
            Current particle number concentration vector [#/m³]
        dndt : np.ndarray
            Time derivative vector [#/(m³·s)]
        tolerance : float, optional
            Tolerance for mass conservation check
            Default is 1e-10

        Returns:
        --------
        is_conserved : bool
            True if mass is conserved within tolerance
        mass_rate : float
            Rate of change of total mass [m³/(m³·s)]

        Notes:
        ------
        Total mass: M = Σ_i N_i * v_i
        Rate of change: dM/dt = Σ_i (dN_i/dt) * v_i
        For pure aggregation: dM/dt ≈ 0
        """
        # Calculate rate of change of total mass
        mass_rate = np.sum(dndt * self.particle_volumes)

        # Check if approximately zero
        is_conserved = abs(mass_rate) < tolerance

        return is_conserved, mass_rate


def construct_beta_matrix(
    n_classes: int,
    particle_diameters: np.ndarray,
    temperature: float = 298.15,
    viscosity: float = 8.9e-4,
    density_fluid: float = 998.2,
    density_particle: float = 1050.0,
    shear_rate: float = 0.0
) -> np.ndarray:
    """
    Construct the collision frequency matrix β from physical mechanisms.

    The collision frequency β_{i,j} is constructed from three mechanisms:
    1. Perikinetic (Brownian motion)
    2. Differential sedimentation
    3. Orthokinetic (shear-induced)

    β_{i,j} = β_perikinetic + β_differential_sedimentation + β_orthokinetic

    Parameters:
    -----------
    n_classes : int
        Number of discrete size classes
    particle_diameters : np.ndarray
        Diameter of particles in each size class [m]
        Shape: (n_classes,)
    temperature : float, optional
        Temperature [K], default 298.15 K (25°C)
    viscosity : float, optional
        Dynamic viscosity [Pa·s], default 8.9e-4 (water at 25°C)
    density_fluid : float, optional
        Fluid density [kg/m³], default 998.2 (water at 25°C)
    density_particle : float, optional
        Particle density [kg/m³], default 1050.0
    shear_rate : float, optional
        Shear rate [1/s], default 0.0 (no shear)

    Returns:
    --------
    beta_matrix : np.ndarray
        Collision frequency matrix [m³/s]
        Shape: (n_classes, n_classes)

    Notes:
    ------
    Perikinetic: β_P = (2 k_B T / 3 μ) * (d_i + d_j)² / (d_i * d_j)
    Differential sedimentation: β_DS = (π g |Δρ| / 72 μ) * (d_i + d_j)³ * |d_i - d_j|
    Orthokinetic: β_O = (4/3) * G * (d_i + d_j)³

    Where:
    - k_B: Boltzmann constant
    - T: temperature
    - μ: dynamic viscosity
    - g: gravitational acceleration
    - Δρ: density difference
    - G: shear rate

    Limitations:
    ------------
    - Assumes constant shear rate (no spatial/temporal variation)
    - Does not include hydrodynamic interactions
    - Uses point-particle approximation (no fractal dimensions)
    - Does not account for Hamaker constant in perikinetic collisions
    - For production use, consider implementing more sophisticated models
    """
    k_B = 1.380649e-23  # Boltzmann constant [J/K]
    g = 9.81  # Gravitational acceleration [m/s²]

    beta_matrix = np.zeros((n_classes, n_classes), dtype=np.float64)

    for i in range(n_classes):
        for j in range(n_classes):
            d_i = particle_diameters[i]
            d_j = particle_diameters[j]

            # Perikinetic (Brownian motion)
            beta_perikinetic = (
                (2.0 * k_B * temperature / (3.0 * viscosity))
                * (d_i + d_j)**2 / (d_i * d_j)
            )

            # Differential sedimentation
            delta_rho = abs(density_particle - density_fluid)
            beta_diff_sed = (
                (np.pi * g * delta_rho / (72.0 * viscosity))
                * (d_i + d_j)**3 * abs(d_i - d_j)
            )

            # Orthokinetic (shear-induced)
            beta_orthokinetic = (4.0 / 3.0) * shear_rate * (d_i + d_j)**3

            # Total collision frequency
            beta_matrix[i, j] = beta_perikinetic + beta_diff_sed + beta_orthokinetic

    return beta_matrix
