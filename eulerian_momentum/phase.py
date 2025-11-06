"""
Phase Class - Represents a single phase in multiphase flow
"""

import numpy as np
from typing import Optional


class Phase:
    """
    Represents a single phase (liquid, gas, or solid) in an Eulerian multiphase model.

    Attributes:
        name (str): Phase name (e.g., 'water', 'air', 'particles')
        alpha (np.ndarray): Volume fraction field (dimensionless, 0 ≤ α ≤ 1)
        rho (float): Phase density [kg/m³]
        velocity (np.ndarray): Velocity field [m/s], shape (nx, ny, nz, 3) or (n, 3)
        mu (float): Dynamic viscosity [Pa·s]
        phase_id (int): Unique phase identifier
    """

    def __init__(
        self,
        name: str,
        rho: float,
        mu: float,
        alpha: Optional[np.ndarray] = None,
        velocity: Optional[np.ndarray] = None,
        phase_id: int = 0
    ):
        """
        Initialize a phase.

        Args:
            name: Phase name
            rho: Density [kg/m³]
            mu: Dynamic viscosity [Pa·s]
            alpha: Volume fraction field (if None, will be initialized as needed)
            velocity: Velocity field [m/s] (if None, will be initialized as needed)
            phase_id: Unique identifier for the phase
        """
        self.name = name
        self.rho = rho
        self.mu = mu
        self.alpha = alpha
        self.velocity = velocity
        self.phase_id = phase_id

        # Validate inputs
        if self.rho <= 0:
            raise ValueError(f"Density must be positive, got {self.rho}")
        if self.mu < 0:
            raise ValueError(f"Viscosity must be non-negative, got {self.mu}")

    def initialize_fields(self, shape: tuple, initial_alpha: float = 0.0):
        """
        Initialize volume fraction and velocity fields with given shape.

        Args:
            shape: Shape of the computational domain (nx, ny, nz) or (n,)
            initial_alpha: Initial volume fraction value
        """
        if len(shape) == 1:
            # 1D case
            self.alpha = np.full(shape[0], initial_alpha, dtype=np.float64)
            self.velocity = np.zeros((shape[0], 3), dtype=np.float64)
        elif len(shape) == 3:
            # 3D case
            self.alpha = np.full(shape, initial_alpha, dtype=np.float64)
            self.velocity = np.zeros((*shape, 3), dtype=np.float64)
        else:
            raise ValueError(f"Shape must be 1D or 3D, got {len(shape)}D")

    def get_momentum_density(self) -> np.ndarray:
        """
        Calculate momentum density: α_φ * ρ_φ * U_φ

        Returns:
            Momentum density field [kg/(m²·s)]
        """
        if self.alpha is None or self.velocity is None:
            raise ValueError("Fields not initialized. Call initialize_fields() first.")

        # Expand alpha to match velocity shape for broadcasting
        if len(self.alpha.shape) == 1:
            alpha_expanded = self.alpha[:, np.newaxis]
        else:
            alpha_expanded = self.alpha[..., np.newaxis]

        return alpha_expanded * self.rho * self.velocity

    def get_kinematic_viscosity(self) -> float:
        """
        Calculate kinematic viscosity: ν = μ / ρ

        Returns:
            Kinematic viscosity [m²/s]
        """
        return self.mu / self.rho

    def __repr__(self) -> str:
        return (f"Phase(name='{self.name}', rho={self.rho} kg/m³, "
                f"mu={self.mu} Pa·s, phase_id={self.phase_id})")
