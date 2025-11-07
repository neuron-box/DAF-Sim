"""
Interfacial Forces Module

Implements interfacial momentum transfer term M_φ, including:
- Drag force
- Lift force
- Virtual mass force
- Wall lubrication force
- Turbulent dispersion force
"""

import numpy as np
from typing import List, Dict, Optional
from .phase import Phase


class InterfacialForces:
    """
    Computes interfacial momentum transfer between phases.

    The total interfacial force M_φ includes contributions from:
    - Drag: F_D
    - Lift: F_L
    - Virtual mass: F_VM
    - Wall lubrication: F_WL
    - Turbulent dispersion: F_TD
    """

    def __init__(self):
        """Initialize interfacial forces calculator."""
        pass

    def drag_force(
        self,
        phase_continuous: Phase,
        phase_dispersed: Phase,
        drag_coefficient: float = 0.44,
        particle_diameter: float = 1e-4
    ) -> np.ndarray:
        """
        Calculate drag force between continuous and dispersed phases.

        F_D = (3/4) * (C_D / d_p) * α_d * ρ_c * |U_c - U_d| * (U_c - U_d)

        Args:
            phase_continuous: Continuous phase (e.g., water)
            phase_dispersed: Dispersed phase (e.g., bubbles, particles)
            drag_coefficient: Drag coefficient C_D (default: 0.44 for spheres)
            particle_diameter: Particle/bubble diameter [m]

        Returns:
            Drag force per unit volume [N/m³]
        """
        # Relative velocity
        U_rel = phase_continuous.velocity - phase_dispersed.velocity

        # Magnitude of relative velocity
        U_rel_mag = np.linalg.norm(U_rel, axis=-1, keepdims=True)

        # Avoid division by zero
        U_rel_mag = np.where(U_rel_mag < 1e-12, 1e-12, U_rel_mag)

        # Volume fraction of dispersed phase
        if len(phase_dispersed.alpha.shape) == 1:
            alpha_d = phase_dispersed.alpha[:, np.newaxis]
        else:
            alpha_d = phase_dispersed.alpha[..., np.newaxis]

        # Drag force
        drag_coeff = (3.0 / 4.0) * (drag_coefficient / particle_diameter)
        F_drag = drag_coeff * alpha_d * phase_continuous.rho * U_rel_mag * U_rel

        return F_drag

    def lift_force(
        self,
        phase_continuous: Phase,
        phase_dispersed: Phase,
        lift_coefficient: float = 0.5,
        dx: float = 0.01,
        dy: float = None,
        dz: float = None
    ) -> np.ndarray:
        """
        Calculate lift force (e.g., Saffman lift).

        F_L = C_L * α_d * ρ_c * (U_c - U_d) × (∇ × U_c)

        Args:
            phase_continuous: Continuous phase
            phase_dispersed: Dispersed phase
            lift_coefficient: Lift coefficient C_L
            dx, dy, dz: Grid spacings for computing curl

        Returns:
            Lift force per unit volume [N/m³]
        """
        # Relative velocity
        U_rel = phase_continuous.velocity - phase_dispersed.velocity

        # Compute vorticity (curl of continuous phase velocity)
        curl_Uc = self._calculate_curl(phase_continuous.velocity, dx, dy, dz)

        # Cross product: U_rel × curl_Uc
        cross_product = np.cross(U_rel, curl_Uc)

        # Volume fraction
        if len(phase_dispersed.alpha.shape) == 1:
            alpha_d = phase_dispersed.alpha[:, np.newaxis]
        else:
            alpha_d = phase_dispersed.alpha[..., np.newaxis]

        # Lift force
        F_lift = lift_coefficient * alpha_d * phase_continuous.rho * cross_product

        return F_lift

    def virtual_mass_force(
        self,
        phase_continuous: Phase,
        phase_dispersed: Phase,
        dt: float,
        virtual_mass_coefficient: float = 0.5,
        U_dispersed_old: Optional[np.ndarray] = None,
        U_continuous_old: Optional[np.ndarray] = None
    ) -> np.ndarray:
        """
        Calculate virtual mass force (added mass effect).

        F_VM = C_VM * α_d * ρ_c * (DU_c/Dt - DU_d/Dt)

        where D/Dt is the material derivative.

        Args:
            phase_continuous: Continuous phase
            phase_dispersed: Dispersed phase
            dt: Time step [s]
            virtual_mass_coefficient: Virtual mass coefficient C_VM (default: 0.5)
            U_dispersed_old: Previous dispersed phase velocity
            U_continuous_old: Previous continuous phase velocity

        Returns:
            Virtual mass force per unit volume [N/m³]
        """
        # Calculate acceleration of each phase
        if U_continuous_old is None:
            accel_c = np.zeros_like(phase_continuous.velocity)
        else:
            accel_c = (phase_continuous.velocity - U_continuous_old) / dt

        if U_dispersed_old is None:
            accel_d = np.zeros_like(phase_dispersed.velocity)
        else:
            accel_d = (phase_dispersed.velocity - U_dispersed_old) / dt

        # Relative acceleration
        accel_rel = accel_c - accel_d

        # Volume fraction
        if len(phase_dispersed.alpha.shape) == 1:
            alpha_d = phase_dispersed.alpha[:, np.newaxis]
        else:
            alpha_d = phase_dispersed.alpha[..., np.newaxis]

        # Virtual mass force
        F_vm = virtual_mass_coefficient * alpha_d * phase_continuous.rho * accel_rel

        return F_vm

    def wall_lubrication_force(
        self,
        phase_dispersed: Phase,
        wall_normal: np.ndarray,
        distance_to_wall: np.ndarray,
        particle_diameter: float = 1e-4,
        C_wl: float = 0.1
    ) -> np.ndarray:
        """
        Calculate wall lubrication force (prevents particles from touching walls).

        This is a simplified model. More sophisticated models exist.

        Args:
            phase_dispersed: Dispersed phase
            wall_normal: Normal vector to the wall (unit vector)
            distance_to_wall: Distance from each point to nearest wall [m]
            particle_diameter: Particle diameter [m]
            C_wl: Wall lubrication coefficient

        Returns:
            Wall lubrication force per unit volume [N/m³]
        """
        # Expand arrays for proper broadcasting
        if len(phase_dispersed.alpha.shape) == 1:
            d_wall = distance_to_wall[:, np.newaxis]
        else:
            d_wall = distance_to_wall[..., np.newaxis]

        # Force magnitude (inversely proportional to distance)
        # Avoid singularity at wall
        d_wall_safe = np.maximum(d_wall, particle_diameter / 10)

        force_magnitude = C_wl * (particle_diameter / d_wall_safe) ** 2

        # Force direction (away from wall)
        F_wl = force_magnitude * wall_normal

        return F_wl

    def turbulent_dispersion_force(
        self,
        phase_continuous: Phase,
        phase_dispersed: Phase,
        turbulent_diffusivity: float,
        dx: float,
        dy: float = None,
        dz: float = None
    ) -> np.ndarray:
        """
        Calculate turbulent dispersion force (drives particles down concentration gradient).

        F_TD = -C_TD * ρ_c * D_t * ∇α_d

        Args:
            phase_continuous: Continuous phase
            phase_dispersed: Dispersed phase
            turbulent_diffusivity: Turbulent diffusivity D_t [m²/s]
            dx, dy, dz: Grid spacings

        Returns:
            Turbulent dispersion force per unit volume [N/m³]
        """
        # Calculate gradient of dispersed phase volume fraction
        grad_alpha = self._calculate_gradient(phase_dispersed.alpha, dx, dy, dz)

        # Turbulent dispersion force
        F_td = -phase_continuous.rho * turbulent_diffusivity * grad_alpha

        return F_td

    def total_interfacial_force(
        self,
        phase_continuous: Phase,
        phase_dispersed: Phase,
        dt: float,
        dx: float,
        dy: float = None,
        dz: float = None,
        include_drag: bool = True,
        include_lift: bool = False,
        include_virtual_mass: bool = False,
        include_wall_lubrication: bool = False,
        include_turbulent_dispersion: bool = False,
        **kwargs
    ) -> np.ndarray:
        """
        Calculate total interfacial momentum transfer M_φ.

        Args:
            phase_continuous: Continuous phase
            phase_dispersed: Dispersed phase
            dt: Time step [s]
            dx, dy, dz: Grid spacings [m]
            include_*: Flags to include specific force contributions
            **kwargs: Additional parameters for force models

        Returns:
            Total interfacial force M_φ [N/m³]
        """
        M_phi = np.zeros_like(phase_dispersed.velocity, dtype=np.float64)

        if include_drag:
            drag_coeff = kwargs.get('drag_coefficient', 0.44)
            diameter = kwargs.get('particle_diameter', 1e-4)
            M_phi += self.drag_force(phase_continuous, phase_dispersed, drag_coeff, diameter)

        if include_lift:
            lift_coeff = kwargs.get('lift_coefficient', 0.5)
            M_phi += self.lift_force(phase_continuous, phase_dispersed, lift_coeff, dx, dy, dz)

        if include_virtual_mass:
            vm_coeff = kwargs.get('virtual_mass_coefficient', 0.5)
            U_d_old = kwargs.get('U_dispersed_old', None)
            U_c_old = kwargs.get('U_continuous_old', None)
            M_phi += self.virtual_mass_force(
                phase_continuous, phase_dispersed, dt, vm_coeff, U_d_old, U_c_old
            )

        if include_wall_lubrication:
            wall_normal = kwargs.get('wall_normal', np.array([0, 0, 1]))
            distance = kwargs.get('distance_to_wall', np.ones_like(phase_dispersed.alpha))
            diameter = kwargs.get('particle_diameter', 1e-4)
            C_wl = kwargs.get('C_wl', 0.1)
            M_phi += self.wall_lubrication_force(phase_dispersed, wall_normal, distance, diameter, C_wl)

        if include_turbulent_dispersion:
            D_t = kwargs.get('turbulent_diffusivity', 1e-5)
            M_phi += self.turbulent_dispersion_force(
                phase_continuous, phase_dispersed, D_t, dx, dy, dz
            )

        return M_phi

    # ==================== Helper Methods ====================

    def _calculate_gradient(
        self,
        field: np.ndarray,
        dx: float,
        dy: float = None,
        dz: float = None
    ) -> np.ndarray:
        """Calculate gradient of a scalar field."""
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

            grad[1:-1, :, :, 0] = (field[2:, :, :] - field[:-2, :, :]) / (2 * dx)
            grad[0, :, :, 0] = (field[1, :, :] - field[0, :, :]) / dx
            grad[-1, :, :, 0] = (field[-1, :, :] - field[-2, :, :]) / dx

            if dy is not None:
                grad[:, 1:-1, :, 1] = (field[:, 2:, :] - field[:, :-2, :]) / (2 * dy)
                grad[:, 0, :, 1] = (field[:, 1, :] - field[:, 0, :]) / dy
                grad[:, -1, :, 1] = (field[:, -1, :] - field[:, -2, :]) / dy

            if dz is not None:
                grad[:, :, 1:-1, 2] = (field[:, :, 2:] - field[:, :, :-2]) / (2 * dz)
                grad[:, :, 0, 2] = (field[:, :, 1] - field[:, :, 0]) / dz
                grad[:, :, -1, 2] = (field[:, :, -1] - field[:, :, -2]) / dz

        return grad

    def _calculate_curl(
        self,
        velocity: np.ndarray,
        dx: float,
        dy: float = None,
        dz: float = None
    ) -> np.ndarray:
        """Calculate curl (vorticity) of a vector field."""
        if len(velocity.shape) == 2:
            # 1D case - simplified (curl not fully defined in 1D)
            # Return zero curl
            return np.zeros_like(velocity)
        else:
            # 3D case
            nx, ny, nz, _ = velocity.shape
            curl = np.zeros((nx, ny, nz, 3), dtype=np.float64)

            # ω_x = ∂w/∂y - ∂v/∂z
            if dy is not None:
                curl[:, 1:-1, :, 0] += (velocity[:, 2:, :, 2] - velocity[:, :-2, :, 2]) / (2 * dy)
            if dz is not None:
                curl[:, :, 1:-1, 0] -= (velocity[:, :, 2:, 1] - velocity[:, :, :-2, 1]) / (2 * dz)

            # ω_y = ∂u/∂z - ∂w/∂x
            if dz is not None:
                curl[:, :, 1:-1, 1] += (velocity[:, :, 2:, 0] - velocity[:, :, :-2, 0]) / (2 * dz)
            curl[1:-1, :, :, 1] -= (velocity[2:, :, :, 2] - velocity[:-2, :, :, 2]) / (2 * dx)

            # ω_z = ∂v/∂x - ∂u/∂y
            curl[1:-1, :, :, 2] += (velocity[2:, :, :, 1] - velocity[:-2, :, :, 1]) / (2 * dx)
            if dy is not None:
                curl[:, 1:-1, :, 2] -= (velocity[:, 2:, :, 0] - velocity[:, :-2, :, 0]) / (2 * dy)

            return curl
