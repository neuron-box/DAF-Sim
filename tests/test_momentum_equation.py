"""
Unit tests for MomentumEquation class
"""

import unittest
import numpy as np
from eulerian_momentum.phase import Phase
from eulerian_momentum.momentum_equation import MomentumEquation


class TestMomentumEquation(unittest.TestCase):
    """Test cases for MomentumEquation class."""

    def setUp(self):
        """Set up test fixtures."""
        # Create a simple 1D test case
        self.n = 50
        self.dx = 0.01  # 1 cm spacing
        self.dt = 0.001  # 1 ms time step

        # Create water phase
        self.phase = Phase(name="water", rho=1000.0, mu=1e-3)
        self.phase.initialize_fields((self.n,), initial_alpha=0.8)

        # Set up a simple velocity profile (linear)
        x = np.linspace(0, (self.n - 1) * self.dx, self.n)
        self.phase.velocity[:, 0] = 0.1 * x  # u = 0.1 * x

        # Create momentum equation solver
        self.solver = MomentumEquation(gravity=np.array([0, 0, -9.81]))

        # Pressure field (linear gradient)
        self.pressure = 101325.0 - 1000.0 * x  # 1000 Pa/m gradient

    def test_initialization(self):
        """Test momentum equation initialization."""
        solver = MomentumEquation()
        self.assertTrue(np.allclose(solver.gravity, [0, 0, -9.81]))

        custom_gravity = np.array([0, 0, -10.0])
        solver_custom = MomentumEquation(gravity=custom_gravity)
        self.assertTrue(np.allclose(solver_custom.gravity, custom_gravity))

    def test_transient_term_zero_velocity(self):
        """Test transient term with zero initial velocity."""
        self.phase.velocity[:, :] = 0.0
        momentum_old = np.zeros_like(self.phase.velocity)

        transient = self.solver.transient_term(self.phase, self.dt, momentum_old)

        # Should be zero since velocity hasn't changed
        self.assertTrue(np.allclose(transient, 0.0))

    def test_transient_term_acceleration(self):
        """Test transient term with uniform acceleration."""
        # Initial velocity = 0, final velocity = 1 m/s
        momentum_old = np.zeros((self.n, 3))
        self.phase.velocity[:, 0] = 1.0

        transient = self.solver.transient_term(self.phase, self.dt, momentum_old)

        # Expected: (α * ρ * 1.0 - 0) / dt = 0.8 * 1000 * 1.0 / 0.001 = 800000
        expected = 0.8 * 1000.0 * 1.0 / self.dt
        self.assertTrue(np.allclose(transient[:, 0], expected))

    def test_pressure_gradient_term_uniform_pressure(self):
        """Test pressure gradient with uniform pressure."""
        pressure_uniform = np.ones(self.n) * 101325.0

        grad_p = self.solver.pressure_gradient_term(
            self.phase, pressure_uniform, self.dx
        )

        # Gradient should be zero for uniform pressure
        # (except at boundaries due to one-sided differences)
        self.assertTrue(np.allclose(grad_p[1:-1, :], 0.0, atol=1e-6))

    def test_pressure_gradient_term_linear_gradient(self):
        """Test pressure gradient with linear pressure field."""
        # Pressure with constant gradient: dp/dx = -1000 Pa/m
        x = np.linspace(0, (self.n - 1) * self.dx, self.n)
        pressure = 101325.0 - 1000.0 * x

        grad_p = self.solver.pressure_gradient_term(self.phase, pressure, self.dx)

        # Expected: α * (-1000) in x-direction
        expected_x = 0.8 * (-1000.0)

        # Check interior points (boundaries may have different values due to one-sided differences)
        self.assertTrue(np.allclose(grad_p[2:-2, 0], expected_x, rtol=0.1))

    def test_gravity_term(self):
        """Test gravitational body force term."""
        gravity_force = self.solver.gravity_term(self.phase)

        # Expected: α * ρ * g = 0.8 * 1000 * [0, 0, -9.81]
        expected_z = 0.8 * 1000.0 * (-9.81)

        self.assertTrue(np.allclose(gravity_force[:, 0], 0.0))
        self.assertTrue(np.allclose(gravity_force[:, 1], 0.0))
        self.assertTrue(np.allclose(gravity_force[:, 2], expected_z))

    def test_compression_expansion_term(self):
        """Test compression-expansion source term."""
        # Create a source term
        S_ce = np.ones(self.n) * 10.0  # 10 kg/(m³·s)
        self.phase.velocity[:, 0] = 2.0  # u = 2 m/s

        ce_term = self.solver.compression_expansion_term(self.phase, S_ce)

        # Expected: S_ce * U = 10 * 2 = 20 in x-direction
        expected = 10.0 * 2.0
        self.assertTrue(np.allclose(ce_term[:, 0], expected))

    def test_reynolds_stress_term_zero_gradient(self):
        """Test Reynolds stress term with zero velocity gradient."""
        # Uniform velocity
        self.phase.velocity[:, 0] = 5.0

        stress = self.solver.reynolds_stress_term(self.phase, self.dx)

        # Should be near zero for uniform velocity (except boundaries)
        self.assertTrue(np.allclose(stress[2:-2, :], 0.0, atol=1e-6))

    def test_convection_term_uniform_velocity(self):
        """Test convection term with uniform velocity."""
        # Set uniform velocity
        self.phase.velocity[:, 0] = 3.0
        self.phase.alpha[:] = 0.8

        convection = self.solver.convection_term(self.phase, self.dx)

        # For uniform velocity and volume fraction, convection should be zero
        # (except at boundaries)
        self.assertTrue(np.allclose(convection[2:-2, :], 0.0, atol=1e-3))

    def test_compute_lhs_dimensions(self):
        """Test that LHS computation returns correct dimensions."""
        lhs = self.solver.compute_lhs(
            self.phase,
            self.pressure,
            self.dt,
            self.dx
        )

        self.assertEqual(lhs.shape, self.phase.velocity.shape)

    def test_momentum_conservation_1d(self):
        """Test basic momentum conservation in a simple scenario."""
        # Create two time steps with small velocity change
        self.phase.velocity[:, 0] = 1.0
        momentum_old = self.phase.get_momentum_density()

        # Advance velocity slightly
        self.phase.velocity[:, 0] = 1.01

        transient = self.solver.transient_term(self.phase, self.dt, momentum_old)

        # Check that transient term is proportional to velocity change
        delta_u = 0.01
        expected_magnitude = self.phase.alpha[0] * self.phase.rho * delta_u / self.dt

        self.assertAlmostEqual(
            np.mean(np.abs(transient[:, 0])),
            expected_magnitude,
            delta=expected_magnitude * 0.01
        )


class TestMomentumEquation3D(unittest.TestCase):
    """Test cases for 3D momentum equation."""

    def setUp(self):
        """Set up 3D test case."""
        self.nx, self.ny, self.nz = 10, 10, 10
        self.dx = self.dy = self.dz = 0.01

        self.phase = Phase(name="water", rho=1000.0, mu=1e-3)
        self.phase.initialize_fields((self.nx, self.ny, self.nz), initial_alpha=0.9)

        self.solver = MomentumEquation()
        self.pressure = np.ones((self.nx, self.ny, self.nz)) * 101325.0

    def test_3d_gravity_term(self):
        """Test gravity in 3D."""
        gravity_force = self.solver.gravity_term(self.phase)

        self.assertEqual(gravity_force.shape, (self.nx, self.ny, self.nz, 3))

        expected_z = 0.9 * 1000.0 * (-9.81)
        self.assertTrue(np.allclose(gravity_force[:, :, :, 2], expected_z))

    def test_3d_pressure_gradient(self):
        """Test 3D pressure gradient."""
        # Create pressure gradient in z-direction
        z = np.arange(self.nz) * self.dz
        for i in range(self.nx):
            for j in range(self.ny):
                self.pressure[i, j, :] = 101325.0 + 1000.0 * 9.81 * z

        grad_p = self.solver.pressure_gradient_term(
            self.phase, self.pressure, self.dx, self.dy, self.dz
        )

        # Should have gradient in z-direction
        expected_z = 0.9 * 1000.0 * 9.81

        # Check interior points
        interior_grad_z = grad_p[2:-2, 2:-2, 2:-2, 2]
        self.assertTrue(np.allclose(interior_grad_z, expected_z, rtol=0.1))


if __name__ == '__main__':
    unittest.main()
