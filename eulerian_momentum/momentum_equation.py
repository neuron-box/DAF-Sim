"""
Eulerian Momentum Equation Implementation

Implements the VTT (2023) Eulerian phase momentum equation:
∂(α_φ ρ_φ U_φ)/∂t + ∇·(α_φ ρ_φ U_φ U_φ) − S_ce,φ U_φ + R_φ + α_φ ∇p − ρ_φ g = M_φ
"""

import numpy as np
from typing import Optional, Tuple
from .phase import Phase


class MomentumEquation:
    """
    Computes all terms of the Eulerian phase momentum equation.

    This class provides methods to calculate each term in the momentum equation
    for a given phase in an Eulerian-Eulerian multiphase flow framework.
    """

    def __init__(self, gravity: np.ndarray = np.array([0, 0, -9.81])):
        """
        Initialize the momentum equation solver.

        Args:
            gravity: Gravitational acceleration vector [m/s²], default is [0, 0, -9.81]
        """
        self.gravity = np.array(gravity, dtype=np.float64)

    def transient_term(
        self,
        phase: Phase,
        dt: float,
        momentum_old: Optional[np.ndarray] = None
    ) -> np.ndarray:
        """
        Calculate the transient term: ∂(α_φ ρ_φ U_φ)/∂t

        Args:
            phase: Phase object containing current state
            dt: Time step [s]
            momentum_old: Previous time step momentum density (if None, assumes zero)

        Returns:
            Transient term [kg/(m²·s²)] or [N/m³]
        """
        momentum_new = phase.get_momentum_density()

        if momentum_old is None:
            # First time step - assume zero initial momentum
            momentum_old = np.zeros_like(momentum_new)

        # Forward Euler time derivative
        d_momentum_dt = (momentum_new - momentum_old) / dt

        return d_momentum_dt

    def convection_term(
        self,
        phase: Phase,
        dx: float,
        dy: float = None,
        dz: float = None
    ) -> np.ndarray:
        """
        Calculate the convection term: ∇·(α_φ ρ_φ U_φ U_φ)

        Uses finite difference approximation for the divergence of the momentum flux.

        Args:
            phase: Phase object
            dx: Grid spacing in x-direction [m]
            dy: Grid spacing in y-direction [m] (if None, assumes 1D)
            dz: Grid spacing in z-direction [m] (if None, assumes 1D or 2D)

        Returns:
            Convection term [kg/(m²·s²)] or [N/m³]
        """
        alpha = phase.alpha
        rho = phase.rho
        U = phase.velocity

        # Calculate momentum density
        if len(alpha.shape) == 1:
            # 1D case
            alpha_exp = alpha[:, np.newaxis]
        else:
            # 3D case
            alpha_exp = alpha[..., np.newaxis]

        momentum_density = alpha_exp * rho * U

        # Calculate flux tensor: (α ρ U) ⊗ U
        # This is a rank-2 tensor (vector × vector)
        if len(alpha.shape) == 1:
            # 1D: shape (n, 3) → flux is (n, 3, 3)
            n = alpha.shape[0]
            flux = np.zeros((n, 3, 3), dtype=np.float64)
            for i in range(3):
                for j in range(3):
                    flux[:, i, j] = momentum_density[:, i] * U[:, j]
        else:
            # 3D: shape (nx, ny, nz, 3)
            nx, ny, nz = alpha.shape
            flux = np.zeros((nx, ny, nz, 3, 3), dtype=np.float64)
            for i in range(3):
                for j in range(3):
                    flux[:, :, :, i, j] = momentum_density[:, :, :, i] * U[:, :, :, j]

        # Calculate divergence using central differences
        div_flux = self._calculate_divergence_tensor(flux, dx, dy, dz)

        return div_flux

    def compression_expansion_term(
        self,
        phase: Phase,
        S_ce: np.ndarray
    ) -> np.ndarray:
        """
        Calculate compression-expansion source term: S_ce,φ * U_φ

        Args:
            phase: Phase object
            S_ce: Compression-expansion source term [kg/(m³·s)]

        Returns:
            CE source term contribution [kg/(m²·s²)] or [N/m³]
        """
        if len(S_ce.shape) == 1:
            S_ce_exp = S_ce[:, np.newaxis]
        else:
            S_ce_exp = S_ce[..., np.newaxis]

        return S_ce_exp * phase.velocity

    def reynolds_stress_term(
        self,
        phase: Phase,
        dx: float,
        dy: float = None,
        dz: float = None,
        turbulent_viscosity: Optional[np.ndarray] = None
    ) -> np.ndarray:
        """
        Calculate Reynolds/turbulent stress term: R_φ = ∇·τ_φ

        For laminar flow: τ_φ = μ_φ (∇U_φ + (∇U_φ)ᵀ)
        For turbulent flow: τ_φ = (μ_φ + μ_t) (∇U_φ + (∇U_φ)ᵀ)

        Args:
            phase: Phase object
            dx: Grid spacing in x-direction [m]
            dy: Grid spacing in y-direction [m]
            dz: Grid spacing in z-direction [m]
            turbulent_viscosity: Turbulent viscosity field [Pa·s] (optional)

        Returns:
            Stress term [kg/(m²·s²)] or [N/m³]
        """
        # Effective viscosity
        mu_eff = phase.mu

        if turbulent_viscosity is not None:
            mu_eff = phase.mu + turbulent_viscosity

        # Calculate velocity gradient tensor
        grad_U = self._calculate_velocity_gradient(phase.velocity, dx, dy, dz)

        # Calculate stress tensor: μ_eff * (∇U + (∇U)ᵀ)
        if len(phase.alpha.shape) == 1:
            # 1D case
            n = phase.alpha.shape[0]
            stress = np.zeros((n, 3, 3), dtype=np.float64)
            mu_eff_scalar = mu_eff if np.isscalar(mu_eff) else mu_eff

            for i in range(3):
                for j in range(3):
                    if np.isscalar(mu_eff_scalar):
                        stress[:, i, j] = mu_eff_scalar * (grad_U[:, i, j] + grad_U[:, j, i])
                    else:
                        stress[:, i, j] = mu_eff_scalar * (grad_U[:, i, j] + grad_U[:, j, i])
        else:
            # 3D case
            nx, ny, nz = phase.alpha.shape
            stress = np.zeros((nx, ny, nz, 3, 3), dtype=np.float64)
            mu_eff_scalar = mu_eff if np.isscalar(mu_eff) else mu_eff

            for i in range(3):
                for j in range(3):
                    if np.isscalar(mu_eff_scalar):
                        stress[:, :, :, i, j] = mu_eff_scalar * (
                            grad_U[:, :, :, i, j] + grad_U[:, :, :, j, i]
                        )
                    else:
                        stress[:, :, :, i, j] = mu_eff_scalar[..., np.newaxis, np.newaxis] * (
                            grad_U[:, :, :, i, j] + grad_U[:, :, :, j, i]
                        )

        # Calculate divergence of stress tensor
        div_stress = self._calculate_divergence_tensor(stress, dx, dy, dz)

        return div_stress

    def pressure_gradient_term(
        self,
        phase: Phase,
        pressure: np.ndarray,
        dx: float,
        dy: float = None,
        dz: float = None
    ) -> np.ndarray:
        """
        Calculate pressure gradient term: α_φ ∇p

        Args:
            phase: Phase object
            pressure: Pressure field [Pa]
            dx: Grid spacing in x-direction [m]
            dy: Grid spacing in y-direction [m]
            dz: Grid spacing in z-direction [m]

        Returns:
            Pressure gradient term [kg/(m²·s²)] or [N/m³]
        """
        # Calculate pressure gradient
        grad_p = self._calculate_gradient(pressure, dx, dy, dz)

        # Multiply by volume fraction
        if len(phase.alpha.shape) == 1:
            alpha_exp = phase.alpha[:, np.newaxis]
        else:
            alpha_exp = phase.alpha[..., np.newaxis]

        return alpha_exp * grad_p

    def gravity_term(self, phase: Phase) -> np.ndarray:
        """
        Calculate gravitational body force term: ρ_φ g

        Note: The VTT equation shows ρ_φ g, but physically it should be α_φ ρ_φ g.
        This implementation uses α_φ ρ_φ g for physical correctness.

        Args:
            phase: Phase object

        Returns:
            Gravity term [kg/(m²·s²)] or [N/m³]
        """
        if len(phase.alpha.shape) == 1:
            n = phase.alpha.shape[0]
            gravity_expanded = np.tile(self.gravity, (n, 1))
            alpha_exp = phase.alpha[:, np.newaxis]
        else:
            nx, ny, nz = phase.alpha.shape
            gravity_expanded = np.tile(self.gravity, (*phase.alpha.shape, 1))
            alpha_exp = phase.alpha[..., np.newaxis]

        return alpha_exp * phase.rho * gravity_expanded

    def compute_lhs(
        self,
        phase: Phase,
        pressure: np.ndarray,
        dt: float,
        dx: float,
        dy: float = None,
        dz: float = None,
        momentum_old: Optional[np.ndarray] = None,
        S_ce: Optional[np.ndarray] = None,
        turbulent_viscosity: Optional[np.ndarray] = None
    ) -> np.ndarray:
        """
        Compute the left-hand side of the momentum equation:
        ∂(α_φ ρ_φ U_φ)/∂t + ∇·(α_φ ρ_φ U_φ U_φ) − S_ce,φ U_φ + R_φ + α_φ ∇p − ρ_φ g

        Args:
            phase: Phase object
            pressure: Pressure field [Pa]
            dt: Time step [s]
            dx: Grid spacing in x [m]
            dy: Grid spacing in y [m]
            dz: Grid spacing in z [m]
            momentum_old: Previous momentum density
            S_ce: Compression-expansion source [kg/(m³·s)]
            turbulent_viscosity: Turbulent viscosity [Pa·s]

        Returns:
            LHS of momentum equation [kg/(m²·s²)] or [N/m³]
        """
        # Initialize LHS
        lhs = np.zeros_like(phase.velocity, dtype=np.float64)

        # 1. Transient term
        lhs += self.transient_term(phase, dt, momentum_old)

        # 2. Convection term
        lhs += self.convection_term(phase, dx, dy, dz)

        # 3. Compression-expansion term (negative sign)
        if S_ce is not None:
            lhs -= self.compression_expansion_term(phase, S_ce)

        # 4. Reynolds stress term
        lhs += self.reynolds_stress_term(phase, dx, dy, dz, turbulent_viscosity)

        # 5. Pressure gradient term
        lhs += self.pressure_gradient_term(phase, pressure, dx, dy, dz)

        # 6. Gravity term (negative sign)
        lhs -= self.gravity_term(phase)

        return lhs

    # ==================== Helper Methods ====================

    def _calculate_gradient(
        self,
        field: np.ndarray,
        dx: float,
        dy: float = None,
        dz: float = None
    ) -> np.ndarray:
        """Calculate gradient of a scalar field using central differences."""
        if len(field.shape) == 1:
            # 1D case - assuming vertical (z-direction) domain for DAF context
            # DAF simulations are inherently vertical (bubbles rise, particles settle)
            n = field.shape[0]
            grad = np.zeros((n, 3), dtype=np.float64)
            # Use dz for 1D vertical gradient, fallback to dx if dz is None
            grid_spacing = dz if dz is not None else dx
            # Store gradient in z-component (index 2) for vertical domain
            grad[1:-1, 2] = (field[2:] - field[:-2]) / (2 * grid_spacing)
            grad[0, 2] = (field[1] - field[0]) / grid_spacing
            grad[-1, 2] = (field[-1] - field[-2]) / grid_spacing
        else:
            # 3D case
            nx, ny, nz = field.shape
            grad = np.zeros((nx, ny, nz, 3), dtype=np.float64)

            # x-direction
            grad[1:-1, :, :, 0] = (field[2:, :, :] - field[:-2, :, :]) / (2 * dx)
            grad[0, :, :, 0] = (field[1, :, :] - field[0, :, :]) / dx
            grad[-1, :, :, 0] = (field[-1, :, :] - field[-2, :, :]) / dx

            if dy is not None:
                # y-direction
                grad[:, 1:-1, :, 1] = (field[:, 2:, :] - field[:, :-2, :]) / (2 * dy)
                grad[:, 0, :, 1] = (field[:, 1, :] - field[:, 0, :]) / dy
                grad[:, -1, :, 1] = (field[:, -1, :] - field[:, -2, :]) / dy

            if dz is not None:
                # z-direction
                grad[:, :, 1:-1, 2] = (field[:, :, 2:] - field[:, :, :-2]) / (2 * dz)
                grad[:, :, 0, 2] = (field[:, :, 1] - field[:, :, 0]) / dz
                grad[:, :, -1, 2] = (field[:, :, -1] - field[:, :, -2]) / dz

        return grad

    def _calculate_velocity_gradient(
        self,
        velocity: np.ndarray,
        dx: float,
        dy: float = None,
        dz: float = None
    ) -> np.ndarray:
        """Calculate velocity gradient tensor ∇U."""
        if len(velocity.shape) == 2:
            # 1D case: velocity shape (n, 3)
            n = velocity.shape[0]
            grad_U = np.zeros((n, 3, 3), dtype=np.float64)

            for component in range(3):
                grad = self._calculate_gradient(velocity[:, component], dx, dy, dz)
                grad_U[:, component, :] = grad
        else:
            # 3D case: velocity shape (nx, ny, nz, 3)
            nx, ny, nz, _ = velocity.shape
            grad_U = np.zeros((nx, ny, nz, 3, 3), dtype=np.float64)

            for component in range(3):
                grad = self._calculate_gradient(velocity[:, :, :, component], dx, dy, dz)
                grad_U[:, :, :, component, :] = grad

        return grad_U

    def _calculate_divergence_tensor(
        self,
        tensor: np.ndarray,
        dx: float,
        dy: float = None,
        dz: float = None
    ) -> np.ndarray:
        """Calculate divergence of a rank-2 tensor field."""
        if len(tensor.shape) == 3:
            # 1D case: tensor shape (n, 3, 3)
            n = tensor.shape[0]
            div = np.zeros((n, 3), dtype=np.float64)

            for i in range(3):
                # Divergence in x-direction only for 1D
                div[1:-1, i] = (tensor[2:, i, 0] - tensor[:-2, i, 0]) / (2 * dx)
                div[0, i] = (tensor[1, i, 0] - tensor[0, i, 0]) / dx
                div[-1, i] = (tensor[-1, i, 0] - tensor[-2, i, 0]) / dx
        else:
            # 3D case: tensor shape (nx, ny, nz, 3, 3)
            nx, ny, nz = tensor.shape[:3]
            div = np.zeros((nx, ny, nz, 3), dtype=np.float64)

            for i in range(3):
                # x-direction contribution
                div[1:-1, :, :, i] += (tensor[2:, :, :, i, 0] - tensor[:-2, :, :, i, 0]) / (2 * dx)
                div[0, :, :, i] += (tensor[1, :, :, i, 0] - tensor[0, :, :, i, 0]) / dx
                div[-1, :, :, i] += (tensor[-1, :, :, i, 0] - tensor[-2, :, :, i, 0]) / dx

                if dy is not None:
                    # y-direction contribution
                    div[:, 1:-1, :, i] += (tensor[:, 2:, :, i, 1] - tensor[:, :-2, :, i, 1]) / (2 * dy)
                    div[:, 0, :, i] += (tensor[:, 1, :, i, 1] - tensor[:, 0, :, i, 1]) / dy
                    div[:, -1, :, i] += (tensor[:, -1, :, i, 1] - tensor[:, -2, :, i, 1]) / dy

                if dz is not None:
                    # z-direction contribution
                    div[:, :, 1:-1, i] += (tensor[:, :, 2:, i, 2] - tensor[:, :, :-2, i, 2]) / (2 * dz)
                    div[:, :, 0, i] += (tensor[:, :, 1, i, 2] - tensor[:, :, 0, i, 2]) / dz
                    div[:, :, -1, i] += (tensor[:, :, -1, i, 2] - tensor[:, :, -2, i, 2]) / dz

        return div
