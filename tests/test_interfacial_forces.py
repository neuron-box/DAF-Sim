"""
Unit tests for InterfacialForces class
"""

import unittest
import numpy as np
from eulerian_momentum.phase import Phase
from eulerian_momentum.interfacial_forces import InterfacialForces


class TestInterfacialForces(unittest.TestCase):
    """Test cases for InterfacialForces class."""

    def setUp(self):
        """Set up test fixtures."""
        self.n = 100
        self.dx = 0.01

        # Continuous phase (water)
        self.water = Phase(name="water", rho=1000.0, mu=1e-3, phase_id=0)
        self.water.initialize_fields((self.n,), initial_alpha=0.9)
        self.water.velocity[:, 2] = 0.0  # Stationary water

        # Dispersed phase (air bubbles)
        self.bubbles = Phase(name="bubbles", rho=1.2, mu=1.8e-5, phase_id=1)
        self.bubbles.initialize_fields((self.n,), initial_alpha=0.1)
        self.bubbles.velocity[:, 2] = 0.1  # Rising at 0.1 m/s

        # Interfacial forces calculator
        self.forces = InterfacialForces()

    def test_initialization(self):
        """Test interfacial forces initialization."""
        forces = InterfacialForces()
        self.assertIsNotNone(forces)

    def test_drag_force_zero_relative_velocity(self):
        """Test drag force with zero relative velocity."""
        # Set same velocity for both phases
        self.water.velocity[:, :] = 0.0
        self.bubbles.velocity[:, :] = 0.0

        drag = self.forces.drag_force(
            self.water,
            self.bubbles,
            drag_coefficient=0.44,
            particle_diameter=1e-3
        )

        # Should be nearly zero
        self.assertTrue(np.allclose(drag, 0.0, atol=1e-6))

    def test_drag_force_direction(self):
        """Test that drag force opposes relative motion."""
        # Bubbles rising, water stationary
        self.water.velocity[:, 2] = 0.0
        self.bubbles.velocity[:, 2] = 0.2  # Rising at 0.2 m/s

        drag = self.forces.drag_force(
            self.water,
            self.bubbles,
            drag_coefficient=0.44,
            particle_diameter=1e-3
        )

        # Drag on dispersed phase should be negative (downward)
        # since bubbles are moving upward relative to water
        # F_drag points in direction of (U_c - U_d) = (0 - 0.2) = -0.2
        self.assertTrue(np.all(drag[:, 2] < 0))

    def test_drag_force_magnitude(self):
        """Test drag force magnitude calculation."""
        self.water.velocity[:, 2] = 0.0
        self.bubbles.velocity[:, 2] = 0.1

        C_D = 0.44
        d_p = 1e-3
        U_rel = 0.1

        drag = self.forces.drag_force(self.water, self.bubbles, C_D, d_p)

        # Expected: (3/4) * (C_D / d_p) * α_d * ρ_c * |U_rel| * U_rel
        expected = (3.0 / 4.0) * (C_D / d_p) * 0.1 * 1000.0 * U_rel * (-U_rel)

        self.assertTrue(np.allclose(drag[:, 2], expected, rtol=0.01))

    def test_lift_force_zero_vorticity(self):
        """Test lift force with zero vorticity (uniform flow)."""
        # Uniform velocity field - no rotation
        self.water.velocity[:, 0] = 1.0
        self.bubbles.velocity[:, 0] = 1.1

        lift = self.forces.lift_force(
            self.water,
            self.bubbles,
            lift_coefficient=0.5,
            dx=self.dx
        )

        # Should be nearly zero for uniform flow (no curl)
        self.assertTrue(np.allclose(lift[2:-2, :], 0.0, atol=1e-6))

    def test_virtual_mass_force_zero_acceleration(self):
        """Test virtual mass force with constant velocity."""
        # Both phases at constant velocity
        self.water.velocity[:, 0] = 1.0
        self.bubbles.velocity[:, 0] = 1.5

        # Previous velocities same as current
        U_water_old = self.water.velocity.copy()
        U_bubbles_old = self.bubbles.velocity.copy()

        dt = 0.001

        vm_force = self.forces.virtual_mass_force(
            self.water,
            self.bubbles,
            dt,
            virtual_mass_coefficient=0.5,
            U_dispersed_old=U_bubbles_old,
            U_continuous_old=U_water_old
        )

        # Should be zero for no acceleration
        self.assertTrue(np.allclose(vm_force, 0.0, atol=1e-6))

    def test_virtual_mass_force_acceleration(self):
        """Test virtual mass force with acceleration."""
        # Water accelerating
        self.water.velocity[:, 0] = 2.0
        U_water_old = np.zeros((self.n, 3))
        U_water_old[:, 0] = 1.0

        # Bubbles at constant velocity
        self.bubbles.velocity[:, 0] = 0.5
        U_bubbles_old = self.bubbles.velocity.copy()

        dt = 0.01

        vm_force = self.forces.virtual_mass_force(
            self.water,
            self.bubbles,
            dt,
            virtual_mass_coefficient=0.5,
            U_dispersed_old=U_bubbles_old,
            U_continuous_old=U_water_old
        )

        # Water acceleration = (2-1)/0.01 = 100 m/s²
        # Expected: C_VM * α_d * ρ_c * (accel_c - accel_d)
        # = 0.5 * 0.1 * 1000 * 100 = 5000 N/m³
        expected = 0.5 * 0.1 * 1000.0 * 100.0

        self.assertTrue(np.allclose(vm_force[:, 0], expected, rtol=0.01))

    def test_wall_lubrication_force_direction(self):
        """Test that wall lubrication force points away from wall."""
        wall_normal = np.array([0, 0, 1])  # Wall at z=0, normal pointing up
        distance = np.linspace(0.01, 0.1, self.n)  # Distance from wall

        wl_force = self.forces.wall_lubrication_force(
            self.bubbles,
            wall_normal,
            distance,
            particle_diameter=1e-3,
            C_wl=0.1
        )

        # Force should be in direction of wall normal (away from wall)
        self.assertTrue(np.all(wl_force[:, 2] > 0))

    def test_turbulent_dispersion_force_uniform_alpha(self):
        """Test turbulent dispersion with uniform volume fraction."""
        # Uniform volume fraction - no gradient
        self.bubbles.alpha[:] = 0.1

        td_force = self.forces.turbulent_dispersion_force(
            self.water,
            self.bubbles,
            turbulent_diffusivity=1e-5,
            dx=self.dx
        )

        # Should be zero for uniform concentration
        self.assertTrue(np.allclose(td_force[2:-2, :], 0.0, atol=1e-6))

    def test_turbulent_dispersion_force_gradient(self):
        """Test turbulent dispersion with concentration gradient."""
        # Create gradient in volume fraction
        self.bubbles.alpha = np.linspace(0.05, 0.15, self.n)

        td_force = self.forces.turbulent_dispersion_force(
            self.water,
            self.bubbles,
            turbulent_diffusivity=1e-5,
            dx=self.dx
        )

        # Force should oppose gradient (Fick's law)
        # Gradient is positive (increasing), so force should be negative
        # Actually, the formula is F_TD = -ρ_c * D_t * ∇α
        # So if ∇α > 0, then F_TD < 0 in that direction
        interior_force = td_force[5:-5, 0]
        self.assertTrue(np.mean(interior_force) < 0)

    def test_total_interfacial_force_drag_only(self):
        """Test total force with only drag enabled."""
        dt = 0.001

        M_phi = self.forces.total_interfacial_force(
            self.water,
            self.bubbles,
            dt,
            self.dx,
            include_drag=True,
            include_lift=False,
            include_virtual_mass=False,
            particle_diameter=1e-3,
            drag_coefficient=0.44
        )

        # Should be same as drag force alone
        drag_only = self.forces.drag_force(
            self.water,
            self.bubbles,
            drag_coefficient=0.44,
            particle_diameter=1e-3
        )

        self.assertTrue(np.allclose(M_phi, drag_only))

    def test_total_interfacial_force_multiple_contributions(self):
        """Test total force with multiple contributions."""
        dt = 0.001

        # Set previous velocities
        U_water_old = self.water.velocity.copy()
        U_bubbles_old = self.bubbles.velocity.copy()
        U_bubbles_old[:, 2] = 0.05  # Was rising slower

        M_phi = self.forces.total_interfacial_force(
            self.water,
            self.bubbles,
            dt,
            self.dx,
            include_drag=True,
            include_virtual_mass=True,
            particle_diameter=1e-3,
            drag_coefficient=0.44,
            virtual_mass_coefficient=0.5,
            U_dispersed_old=U_bubbles_old,
            U_continuous_old=U_water_old
        )

        # Should have contributions from both drag and virtual mass
        self.assertFalse(np.allclose(M_phi, 0.0))


class TestInterfacialForces3D(unittest.TestCase):
    """Test 3D interfacial forces."""

    def setUp(self):
        """Set up 3D test case."""
        self.nx, self.ny, self.nz = 10, 10, 10
        self.dx = 0.01

        self.water = Phase(name="water", rho=1000.0, mu=1e-3)
        self.water.initialize_fields((self.nx, self.ny, self.nz), initial_alpha=0.8)

        self.bubbles = Phase(name="bubbles", rho=1.2, mu=1.8e-5)
        self.bubbles.initialize_fields((self.nx, self.ny, self.nz), initial_alpha=0.2)

        # Create rotational flow for testing lift
        y, z = np.meshgrid(
            np.arange(self.ny) * self.dx,
            np.arange(self.nz) * self.dx,
            indexing='ij'
        )
        for i in range(self.nx):
            # Circular flow in y-z plane
            self.water.velocity[i, :, :, 1] = -z.T
            self.water.velocity[i, :, :, 2] = y.T

        self.forces = InterfacialForces()

    def test_3d_drag_force_shape(self):
        """Test 3D drag force returns correct shape."""
        drag = self.forces.drag_force(
            self.water,
            self.bubbles,
            drag_coefficient=0.44,
            particle_diameter=1e-3
        )

        self.assertEqual(drag.shape, (self.nx, self.ny, self.nz, 3))

    def test_3d_lift_force_with_rotation(self):
        """Test lift force in 3D with rotational flow."""
        # Create a stronger rotational flow with vorticity
        # Set up a simple shear flow with non-zero curl
        for i in range(self.nx):
            for j in range(self.ny):
                for k in range(self.nz):
                    # Shear flow: u_y varies with z
                    self.water.velocity[i, j, k, 1] = 0.1 * (k * self.dx)

        # Bubbles moving perpendicular to the shear
        self.bubbles.velocity[:, :, :, 0] = 0.05

        lift = self.forces.lift_force(
            self.water,
            self.bubbles,
            lift_coefficient=0.5,
            dx=self.dx,
            dy=self.dx,
            dz=self.dx
        )

        # Lift force should be computed without errors
        # Check that the shape is correct
        self.assertEqual(lift.shape, (self.nx, self.ny, self.nz, 3))

        # For this simple shear flow, lift force may be small but should be computable
        # Just verify that we get reasonable finite values
        self.assertTrue(np.all(np.isfinite(lift)))


if __name__ == '__main__':
    unittest.main()
