"""
Unit tests for Phase class
"""

import unittest
import numpy as np
from eulerian_momentum.phase import Phase


class TestPhase(unittest.TestCase):
    """Test cases for Phase class."""

    def setUp(self):
        """Set up test fixtures."""
        self.phase_water = Phase(
            name="water",
            rho=1000.0,
            mu=1e-3,
            phase_id=0
        )

        self.phase_air = Phase(
            name="air",
            rho=1.2,
            mu=1.8e-5,
            phase_id=1
        )

    def test_initialization(self):
        """Test phase initialization."""
        self.assertEqual(self.phase_water.name, "water")
        self.assertEqual(self.phase_water.rho, 1000.0)
        self.assertEqual(self.phase_water.mu, 1e-3)
        self.assertEqual(self.phase_water.phase_id, 0)

    def test_invalid_density(self):
        """Test that negative density raises error."""
        with self.assertRaises(ValueError):
            Phase(name="test", rho=-100.0, mu=1e-3)

    def test_invalid_viscosity(self):
        """Test that negative viscosity raises error."""
        with self.assertRaises(ValueError):
            Phase(name="test", rho=1000.0, mu=-1e-3)

    def test_initialize_fields_1d(self):
        """Test 1D field initialization."""
        n = 100
        self.phase_water.initialize_fields((n,), initial_alpha=0.8)

        self.assertEqual(self.phase_water.alpha.shape, (n,))
        self.assertEqual(self.phase_water.velocity.shape, (n, 3))
        self.assertTrue(np.allclose(self.phase_water.alpha, 0.8))
        self.assertTrue(np.allclose(self.phase_water.velocity, 0.0))

    def test_initialize_fields_3d(self):
        """Test 3D field initialization."""
        nx, ny, nz = 10, 15, 20
        self.phase_water.initialize_fields((nx, ny, nz), initial_alpha=0.5)

        self.assertEqual(self.phase_water.alpha.shape, (nx, ny, nz))
        self.assertEqual(self.phase_water.velocity.shape, (nx, ny, nz, 3))
        self.assertTrue(np.allclose(self.phase_water.alpha, 0.5))
        self.assertTrue(np.allclose(self.phase_water.velocity, 0.0))

    def test_get_momentum_density_1d(self):
        """Test momentum density calculation in 1D."""
        n = 50
        self.phase_water.initialize_fields((n,), initial_alpha=0.9)
        self.phase_water.velocity[:, 0] = 1.0  # x-velocity = 1 m/s

        momentum = self.phase_water.get_momentum_density()

        expected = 0.9 * 1000.0 * 1.0  # α * ρ * U
        self.assertTrue(np.allclose(momentum[:, 0], expected))
        self.assertTrue(np.allclose(momentum[:, 1], 0.0))
        self.assertTrue(np.allclose(momentum[:, 2], 0.0))

    def test_get_momentum_density_3d(self):
        """Test momentum density calculation in 3D."""
        nx, ny, nz = 5, 5, 5
        self.phase_water.initialize_fields((nx, ny, nz), initial_alpha=0.7)
        self.phase_water.velocity[:, :, :, 2] = -2.0  # z-velocity = -2 m/s

        momentum = self.phase_water.get_momentum_density()

        expected = 0.7 * 1000.0 * (-2.0)
        self.assertTrue(np.allclose(momentum[:, :, :, 2], expected))

    def test_kinematic_viscosity(self):
        """Test kinematic viscosity calculation."""
        nu_water = self.phase_water.get_kinematic_viscosity()
        self.assertAlmostEqual(nu_water, 1e-6, places=9)

        nu_air = self.phase_air.get_kinematic_viscosity()
        self.assertAlmostEqual(nu_air, 1.5e-5, places=8)

    def test_repr(self):
        """Test string representation."""
        repr_str = repr(self.phase_water)
        self.assertIn("water", repr_str)
        self.assertIn("1000", repr_str)


if __name__ == '__main__':
    unittest.main()
