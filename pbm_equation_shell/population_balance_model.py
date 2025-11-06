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
"""

import numpy as np
from typing import Optional


class PopulationBalanceModel:
    """
    Population Balance Model for flocculation processes.

    This class implements the discretized Population Balance Equation (PBE)
    that describes particle aggregation and breakage in a flocculation system.

    Mathematical Formulation:
    -------------------------
    The discretized PBE for each size class i is:

    dN_i/dt = (Birth by aggregation) - (Death by aggregation)
              + (Birth by breakage) - (Death by breakage)

    Explicitly:

    dN_i/dt = (1/2) * Σ_{j=1}^{i-1} β_{j,i-j} * α_{j,i-j} * N_j * N_{i-j}
              - N_i * Σ_{j=1}^{n} β_{i,j} * α_{i,j} * N_j
              + Σ_{j=i+1}^{n} S_j * γ_{j,i} * N_j
              - S_i * N_i

    Where:
    ------
    N_i : float
        Number concentration of particles in size class i [#/m³]

    β_{i,j} : float
        Collision frequency between particles of size class i and j [m³/s]
        Constructed from three mechanisms:
        - Perikinetic (Brownian motion)
        - Differential sedimentation
        - Orthokinetic (shear-induced)

    α_{i,j} : float
        Collision efficiency between particles of size class i and j [-]
        Dimensionless, ranges from 0 to 1

    S_i : float
        Breakage rate of particles in size class i [1/s]

    γ_{j,i} : float
        Breakage distribution function [-]
        Probability that breakage of particle j produces particle i
        Satisfies: Σ_i γ_{j,i} = 1

    Attributes:
    -----------
    n_classes : int
        Number of discrete size classes
    validate_inputs : bool
        Whether to validate kernel inputs for consistency
    """

    def __init__(self, n_classes: int, validate_inputs: bool = True):
        """
        Initialize the Population Balance Model.

        Parameters:
        -----------
        n_classes : int
            Number of discrete size classes (bins) for particle size distribution
        validate_inputs : bool, optional
            Whether to validate kernel matrices/vectors for physical consistency
            Default is True
        """
        if n_classes < 2:
            raise ValueError("Number of size classes must be at least 2")

        self.n_classes = n_classes
        self.validate_inputs = validate_inputs

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

        # Loop through each size class
        for i in range(self.n_classes):
            # Term 1: Birth by aggregation
            # (1/2) * Σ_{j=1}^{i-1} β_{j,i-j} * α_{j,i-j} * N_j * N_{i-j}
            # Particles of size i are formed by collision of smaller particles j and (i-j)
            birth_aggregation = 0.0
            for j in range(i):  # j from 0 to i-1 (Python 0-indexed)
                k = i - j - 1  # Corresponding size class (i-j in 1-indexed becomes i-j-1 in 0-indexed)
                if k >= 0 and k < self.n_classes:
                    birth_aggregation += (
                        beta_matrix[j, k] * alpha_matrix[j, k] * N[j] * N[k]
                    )
            birth_aggregation *= 0.5  # Factor of 1/2 to avoid double counting

            # Term 2: Death by aggregation
            # -N_i * Σ_{j=1}^{n} β_{i,j} * α_{i,j} * N_j
            # Particles of size i are lost through collision with any other particle j
            death_aggregation = 0.0
            for j in range(self.n_classes):
                death_aggregation += beta_matrix[i, j] * alpha_matrix[i, j] * N[j]
            death_aggregation *= N[i]

            # Term 3: Birth by breakage
            # Σ_{j=i+1}^{n} S_j * γ_{j,i} * N_j
            # Particles of size i are formed by breakage of larger particles j
            birth_breakage = 0.0
            for j in range(i + 1, self.n_classes):
                birth_breakage += S_vector[j] * gamma_matrix[j, i] * N[j]

            # Term 4: Death by breakage
            # -S_i * N_i
            # Particles of size i are lost through their own breakage
            death_breakage = S_vector[i] * N[i]

            # Combine all terms
            dndt[i] = (
                birth_aggregation
                - death_aggregation
                + birth_breakage
                - death_breakage
            )

        return dndt

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

        # Optional: Check that gamma rows sum to approximately 1 (within tolerance)
        # This is commented out as it may be too strict for some applications
        # if self.validate_inputs:
        #     row_sums = np.sum(gamma_matrix, axis=1)
        #     if not np.allclose(row_sums, 1.0, atol=1e-6):
        #         raise ValueError(
        #             "Each row of gamma_matrix should sum to 1 (breakage distribution)"
        #         )

    def validate_mass_conservation(
        self,
        N: np.ndarray,
        dndt: np.ndarray,
        particle_volumes: np.ndarray,
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
        particle_volumes : np.ndarray
            Volume of particles in each size class [m³]
            Shape: (n_classes,)
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
        if particle_volumes.shape != (self.n_classes,):
            raise ValueError(
                f"particle_volumes must have shape ({self.n_classes},), "
                f"got {particle_volumes.shape}"
            )

        # Calculate rate of change of total mass
        mass_rate = np.sum(dndt * particle_volumes)

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
