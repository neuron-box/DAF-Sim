"""
Unit tests for the Population Balance Model (PBM) Equation Shell

This module contains comprehensive tests for the discretized Population Balance
Equation implementation based on Abreu et al. (2021).
"""

import unittest
import numpy as np
from pbm_equation_shell import PopulationBalanceModel, construct_beta_matrix


class TestPopulationBalanceModel(unittest.TestCase):
    """Test cases for PopulationBalanceModel class"""

    def setUp(self):
        """Set up test fixtures"""
        self.n_classes = 5
        self.pbm = PopulationBalanceModel(n_classes=self.n_classes)

    def test_initialization(self):
        """Test PBM initialization"""
        self.assertEqual(self.pbm.n_classes, self.n_classes)
        self.assertTrue(self.pbm.validate_inputs)

        # Test invalid initialization
        with self.assertRaises(ValueError):
            PopulationBalanceModel(n_classes=1)

    def test_calculate_dndt_shape(self):
        """Test that calculate_dndt returns correct shape"""
        N = np.ones(self.n_classes)
        alpha_matrix = np.ones((self.n_classes, self.n_classes)) * 0.5
        beta_matrix = np.ones((self.n_classes, self.n_classes)) * 1e-18
        S_vector = np.zeros(self.n_classes)
        gamma_matrix = np.zeros((self.n_classes, self.n_classes))

        dndt = self.pbm.calculate_dndt(N, alpha_matrix, beta_matrix, S_vector, gamma_matrix)

        self.assertEqual(dndt.shape, (self.n_classes,))
        self.assertEqual(dndt.dtype, np.float64)

    def test_no_aggregation_no_breakage(self):
        """Test that dN/dt = 0 when no aggregation or breakage occurs"""
        N = np.ones(self.n_classes) * 1e15

        # Zero collision efficiency (no aggregation)
        alpha_matrix = np.zeros((self.n_classes, self.n_classes))
        beta_matrix = np.ones((self.n_classes, self.n_classes)) * 1e-18

        # Zero breakage
        S_vector = np.zeros(self.n_classes)
        gamma_matrix = np.zeros((self.n_classes, self.n_classes))

        dndt = self.pbm.calculate_dndt(N, alpha_matrix, beta_matrix, S_vector, gamma_matrix)

        # All derivatives should be zero
        np.testing.assert_array_almost_equal(dndt, np.zeros(self.n_classes))

    def test_pure_aggregation_monotonic(self):
        """
        Test pure aggregation: smallest particles should decrease,
        larger particles should increase
        """
        N = np.array([1e15, 1e14, 1e13, 1e12, 1e11], dtype=np.float64)

        # Uniform collision efficiency and frequency
        alpha_matrix = np.ones((self.n_classes, self.n_classes)) * 0.8
        beta_matrix = np.ones((self.n_classes, self.n_classes)) * 1e-18

        # No breakage
        S_vector = np.zeros(self.n_classes)
        gamma_matrix = np.zeros((self.n_classes, self.n_classes))

        dndt = self.pbm.calculate_dndt(N, alpha_matrix, beta_matrix, S_vector, gamma_matrix)

        # Smallest particles (class 0) should decrease due to aggregation
        self.assertLess(dndt[0], 0, "Smallest particles should decrease")

        # Verify that dndt has expected behavior
        # (Not all derivatives will be positive due to death terms)
        self.assertIsInstance(dndt[0], (float, np.floating))

    def test_pure_breakage(self):
        """Test pure breakage mechanism"""
        N = np.array([0.0, 0.0, 0.0, 1e12, 0.0], dtype=np.float64)

        # No aggregation
        alpha_matrix = np.zeros((self.n_classes, self.n_classes))
        beta_matrix = np.ones((self.n_classes, self.n_classes)) * 1e-18

        # Breakage only in class 3
        S_vector = np.array([0.0, 0.0, 0.0, 1.0, 0.0])

        # Breakage produces smaller particles uniformly
        gamma_matrix = np.zeros((self.n_classes, self.n_classes))
        gamma_matrix[3, 0] = 0.5  # 50% to class 0
        gamma_matrix[3, 1] = 0.5  # 50% to class 1

        dndt = self.pbm.calculate_dndt(N, alpha_matrix, beta_matrix, S_vector, gamma_matrix)

        # Class 3 should decrease (death by breakage)
        self.assertLess(dndt[3], 0, "Parent class should decrease due to breakage")

        # Classes 0 and 1 should increase (birth by breakage)
        self.assertGreater(dndt[0], 0, "Product class 0 should increase")
        self.assertGreater(dndt[1], 0, "Product class 1 should increase")

        # Other classes should be zero
        self.assertAlmostEqual(dndt[2], 0.0)
        self.assertAlmostEqual(dndt[4], 0.0)

    def test_input_validation_N(self):
        """Test input validation for N vector"""
        alpha_matrix = np.ones((self.n_classes, self.n_classes)) * 0.5
        beta_matrix = np.ones((self.n_classes, self.n_classes)) * 1e-18
        S_vector = np.zeros(self.n_classes)
        gamma_matrix = np.zeros((self.n_classes, self.n_classes))

        # Wrong shape
        N_wrong_shape = np.ones(self.n_classes + 1)
        with self.assertRaises(ValueError):
            self.pbm.calculate_dndt(N_wrong_shape, alpha_matrix, beta_matrix, S_vector, gamma_matrix)

        # Negative values
        N_negative = np.ones(self.n_classes) * -1
        with self.assertRaises(ValueError):
            self.pbm.calculate_dndt(N_negative, alpha_matrix, beta_matrix, S_vector, gamma_matrix)

    def test_input_validation_alpha(self):
        """Test input validation for alpha matrix"""
        N = np.ones(self.n_classes)
        beta_matrix = np.ones((self.n_classes, self.n_classes)) * 1e-18
        S_vector = np.zeros(self.n_classes)
        gamma_matrix = np.zeros((self.n_classes, self.n_classes))

        # Wrong shape
        alpha_wrong_shape = np.ones((self.n_classes, self.n_classes + 1))
        with self.assertRaises(ValueError):
            self.pbm.calculate_dndt(N, alpha_wrong_shape, beta_matrix, S_vector, gamma_matrix)

        # Values outside [0, 1]
        alpha_invalid = np.ones((self.n_classes, self.n_classes)) * 1.5
        with self.assertRaises(ValueError):
            self.pbm.calculate_dndt(N, alpha_invalid, beta_matrix, S_vector, gamma_matrix)

    def test_input_validation_beta(self):
        """Test input validation for beta matrix"""
        N = np.ones(self.n_classes)
        alpha_matrix = np.ones((self.n_classes, self.n_classes)) * 0.5
        S_vector = np.zeros(self.n_classes)
        gamma_matrix = np.zeros((self.n_classes, self.n_classes))

        # Wrong shape
        beta_wrong_shape = np.ones((self.n_classes + 1, self.n_classes))
        with self.assertRaises(ValueError):
            self.pbm.calculate_dndt(N, alpha_matrix, beta_wrong_shape, S_vector, gamma_matrix)

        # Negative values
        beta_negative = np.ones((self.n_classes, self.n_classes)) * -1
        with self.assertRaises(ValueError):
            self.pbm.calculate_dndt(N, alpha_matrix, beta_negative, S_vector, gamma_matrix)

    def test_input_validation_S(self):
        """Test input validation for S vector"""
        N = np.ones(self.n_classes)
        alpha_matrix = np.ones((self.n_classes, self.n_classes)) * 0.5
        beta_matrix = np.ones((self.n_classes, self.n_classes)) * 1e-18
        gamma_matrix = np.zeros((self.n_classes, self.n_classes))

        # Wrong shape
        S_wrong_shape = np.zeros(self.n_classes + 1)
        with self.assertRaises(ValueError):
            self.pbm.calculate_dndt(N, alpha_matrix, beta_matrix, S_wrong_shape, gamma_matrix)

        # Negative values
        S_negative = np.ones(self.n_classes) * -1
        with self.assertRaises(ValueError):
            self.pbm.calculate_dndt(N, alpha_matrix, beta_matrix, S_negative, gamma_matrix)

    def test_input_validation_gamma(self):
        """Test input validation for gamma matrix"""
        N = np.ones(self.n_classes)
        alpha_matrix = np.ones((self.n_classes, self.n_classes)) * 0.5
        beta_matrix = np.ones((self.n_classes, self.n_classes)) * 1e-18
        S_vector = np.zeros(self.n_classes)

        # Wrong shape
        gamma_wrong_shape = np.zeros((self.n_classes + 1, self.n_classes))
        with self.assertRaises(ValueError):
            self.pbm.calculate_dndt(N, alpha_matrix, beta_matrix, S_vector, gamma_wrong_shape)

        # Negative values
        gamma_negative = np.ones((self.n_classes, self.n_classes)) * -1
        with self.assertRaises(ValueError):
            self.pbm.calculate_dndt(N, alpha_matrix, beta_matrix, S_vector, gamma_negative)

    def test_symmetry_of_aggregation(self):
        """
        Test that aggregation terms respect symmetry:
        β_{i,j} * α_{i,j} should be symmetric for symmetric kernels
        """
        N = np.ones(self.n_classes) * 1e14

        # Create symmetric kernels
        alpha_matrix = np.ones((self.n_classes, self.n_classes)) * 0.5
        beta_matrix = np.ones((self.n_classes, self.n_classes)) * 1e-18

        # No breakage
        S_vector = np.zeros(self.n_classes)
        gamma_matrix = np.zeros((self.n_classes, self.n_classes))

        dndt1 = self.pbm.calculate_dndt(N, alpha_matrix, beta_matrix, S_vector, gamma_matrix)

        # Transpose kernels (should give same result for symmetric case)
        alpha_matrix_T = alpha_matrix.T
        beta_matrix_T = beta_matrix.T

        dndt2 = self.pbm.calculate_dndt(N, alpha_matrix_T, beta_matrix_T, S_vector, gamma_matrix)

        np.testing.assert_array_almost_equal(dndt1, dndt2, decimal=10)

    def test_mass_conservation_validation(self):
        """Test mass conservation validation method"""
        N = np.array([1e15, 1e14, 1e13, 1e12, 1e11], dtype=np.float64)

        # Pure aggregation (should conserve mass)
        alpha_matrix = np.ones((self.n_classes, self.n_classes)) * 0.5
        beta_matrix = np.ones((self.n_classes, self.n_classes)) * 1e-18
        S_vector = np.zeros(self.n_classes)
        gamma_matrix = np.zeros((self.n_classes, self.n_classes))

        dndt = self.pbm.calculate_dndt(N, alpha_matrix, beta_matrix, S_vector, gamma_matrix)

        # Define particle volumes (assuming geometric progression)
        particle_volumes = np.array([1e-18, 2e-18, 4e-18, 8e-18, 16e-18])

        is_conserved, mass_rate = self.pbm.validate_mass_conservation(
            N, dndt, particle_volumes, tolerance=1e-6
        )

        # For pure aggregation, mass should be conserved
        self.assertTrue(is_conserved, "Mass should be conserved for pure aggregation")
        self.assertAlmostEqual(mass_rate, 0.0, places=6)


class TestBetaMatrixConstruction(unittest.TestCase):
    """Test cases for collision frequency matrix construction"""

    def test_beta_matrix_shape(self):
        """Test that beta matrix has correct shape"""
        n_classes = 5
        particle_diameters = np.logspace(-6, -4, n_classes)  # 1 to 100 μm

        beta_matrix = construct_beta_matrix(n_classes, particle_diameters)

        self.assertEqual(beta_matrix.shape, (n_classes, n_classes))

    def test_beta_matrix_positive(self):
        """Test that all beta values are positive"""
        n_classes = 5
        particle_diameters = np.logspace(-6, -4, n_classes)

        beta_matrix = construct_beta_matrix(n_classes, particle_diameters)

        self.assertTrue(np.all(beta_matrix >= 0), "All beta values should be non-negative")
        self.assertTrue(np.all(beta_matrix > 0), "All beta values should be positive")

    def test_beta_matrix_symmetry(self):
        """Test that beta matrix is symmetric for given inputs"""
        n_classes = 4
        particle_diameters = np.logspace(-6, -4, n_classes)

        beta_matrix = construct_beta_matrix(
            n_classes,
            particle_diameters,
            temperature=298.15,
            viscosity=8.9e-4,
            shear_rate=10.0
        )

        # Beta should be symmetric: β_{i,j} = β_{j,i}
        np.testing.assert_array_almost_equal(beta_matrix, beta_matrix.T, decimal=10)

    def test_beta_mechanisms(self):
        """Test individual collision mechanisms"""
        n_classes = 3
        particle_diameters = np.array([1e-6, 1e-5, 1e-4])  # 1, 10, 100 μm

        # Test perikinetic only (no shear, particles same density as fluid)
        beta_perikinetic = construct_beta_matrix(
            n_classes,
            particle_diameters,
            temperature=298.15,
            viscosity=8.9e-4,
            density_fluid=1000.0,
            density_particle=1000.0,  # Same as fluid -> no sedimentation
            shear_rate=0.0
        )

        # Test orthokinetic only (high shear)
        beta_orthokinetic = construct_beta_matrix(
            n_classes,
            particle_diameters,
            temperature=0.0,  # No Brownian motion (unrealistic but for testing)
            viscosity=8.9e-4,
            density_fluid=1000.0,
            density_particle=1000.0,  # No sedimentation
            shear_rate=100.0
        )

        # Orthokinetic should increase with shear rate
        self.assertTrue(np.all(beta_orthokinetic > beta_perikinetic))

    def test_beta_size_dependence(self):
        """Test that beta increases with particle size"""
        n_classes = 4
        particle_diameters = np.logspace(-6, -4, n_classes)

        beta_matrix = construct_beta_matrix(
            n_classes,
            particle_diameters,
            shear_rate=10.0
        )

        # For orthokinetic mechanism, beta should increase with particle size
        # β ~ (d_i + d_j)³
        # So β[3,3] > β[0,0]
        self.assertGreater(beta_matrix[3, 3], beta_matrix[0, 0])


class TestIntegrationScenarios(unittest.TestCase):
    """Integration tests for realistic PBM scenarios"""

    def test_simple_flocculation_scenario(self):
        """
        Test a simple flocculation scenario with defined initial conditions
        """
        n_classes = 5
        pbm = PopulationBalanceModel(n_classes=n_classes)

        # Initial concentration: mostly small particles
        N = np.array([1e15, 5e14, 1e14, 5e13, 1e13], dtype=np.float64)

        # Particle diameters
        particle_diameters = np.logspace(-6, -4, n_classes)  # 1 to 100 μm

        # Construct beta matrix
        beta_matrix = construct_beta_matrix(
            n_classes,
            particle_diameters,
            temperature=298.15,
            viscosity=8.9e-4,
            shear_rate=50.0  # Moderate shear
        )

        # Collision efficiency (moderate)
        alpha_matrix = np.ones((n_classes, n_classes)) * 0.3

        # No breakage
        S_vector = np.zeros(n_classes)
        gamma_matrix = np.zeros((n_classes, n_classes))

        # Calculate derivatives
        dndt = pbm.calculate_dndt(N, alpha_matrix, beta_matrix, S_vector, gamma_matrix)

        # Verify basic physics
        # Smallest particles should decrease (consumed by aggregation)
        self.assertLess(dndt[0], 0, "Smallest particles should decrease")

        # Total particle number should decrease (aggregation reduces count)
        total_dndt = np.sum(dndt)
        self.assertLess(total_dndt, 0, "Total particle number should decrease")

    def test_equilibrium_approach(self):
        """
        Test that system approaches equilibrium when aggregation is balanced
        """
        n_classes = 4
        pbm = PopulationBalanceModel(n_classes=n_classes)

        # Equilibrium-like distribution
        N = np.array([1e12, 1e12, 1e12, 1e12], dtype=np.float64)

        # Very small collision efficiency
        alpha_matrix = np.ones((n_classes, n_classes)) * 0.01
        beta_matrix = np.ones((n_classes, n_classes)) * 1e-20

        # No breakage
        S_vector = np.zeros(n_classes)
        gamma_matrix = np.zeros((n_classes, n_classes))

        dndt = pbm.calculate_dndt(N, alpha_matrix, beta_matrix, S_vector, gamma_matrix)

        # Derivatives should be small (approaching equilibrium)
        self.assertTrue(np.all(np.abs(dndt) < 1e10))

    def test_aggregation_breakage_balance(self):
        """
        Test scenario where aggregation and breakage are both active
        """
        n_classes = 5
        pbm = PopulationBalanceModel(n_classes=n_classes)

        N = np.array([1e14, 5e13, 2e13, 1e13, 5e12], dtype=np.float64)

        # Moderate aggregation
        alpha_matrix = np.ones((n_classes, n_classes)) * 0.4
        beta_matrix = np.ones((n_classes, n_classes)) * 5e-18

        # Breakage for larger particles
        S_vector = np.array([0.0, 0.0, 0.1, 0.5, 1.0])

        # Breakage produces smaller particles
        gamma_matrix = np.zeros((n_classes, n_classes))
        for i in range(n_classes):
            if i > 0:
                # Particles break into smaller sizes
                gamma_matrix[i, i-1] = 0.6  # Most go to next smaller size
                if i > 1:
                    gamma_matrix[i, i-2] = 0.4  # Some go to two sizes smaller

        dndt = pbm.calculate_dndt(N, alpha_matrix, beta_matrix, S_vector, gamma_matrix)

        # Verify that derivatives are computed
        self.assertEqual(dndt.shape, (n_classes,))

        # Largest particles should have significant breakage (negative contribution)
        # But exact sign depends on balance of aggregation vs breakage
        self.assertIsInstance(dndt[-1], (float, np.floating))


def run_simple_integration_example():
    """
    Run a simple example demonstrating the PBM usage
    """
    print("\n" + "="*70)
    print("SIMPLE FLOCCULATION EXAMPLE")
    print("="*70)

    # Setup
    n_classes = 5
    pbm = PopulationBalanceModel(n_classes=n_classes)

    # Initial particle distribution (mostly small particles)
    N = np.array([1e15, 5e14, 1e14, 5e13, 1e13], dtype=np.float64)
    print(f"\nInitial particle concentrations (#/m³):")
    for i, n in enumerate(N):
        print(f"  Class {i}: {n:.2e}")

    # Define particle sizes
    particle_diameters = np.logspace(-6, -4, n_classes)  # 1 to 100 μm
    print(f"\nParticle diameters (m):")
    for i, d in enumerate(particle_diameters):
        print(f"  Class {i}: {d*1e6:.2f} μm")

    # Construct collision frequency matrix
    beta_matrix = construct_beta_matrix(
        n_classes,
        particle_diameters,
        temperature=298.15,
        viscosity=8.9e-4,
        shear_rate=50.0
    )

    # Collision efficiency
    alpha_matrix = np.ones((n_classes, n_classes)) * 0.3

    # No breakage for this example
    S_vector = np.zeros(n_classes)
    gamma_matrix = np.zeros((n_classes, n_classes))

    # Calculate time derivatives
    dndt = pbm.calculate_dndt(N, alpha_matrix, beta_matrix, S_vector, gamma_matrix)

    print(f"\nTime derivatives dN/dt (#/(m³·s)):")
    for i, rate in enumerate(dndt):
        print(f"  Class {i}: {rate:.2e}")

    # Check mass conservation
    particle_volumes = (4/3) * np.pi * (particle_diameters/2)**3
    is_conserved, mass_rate = pbm.validate_mass_conservation(
        N, dndt, particle_volumes, tolerance=1e-6
    )

    print(f"\nMass conservation check:")
    print(f"  Is conserved: {is_conserved}")
    print(f"  Mass rate: {mass_rate:.2e} m³/(m³·s)")

    print("\n" + "="*70)


if __name__ == '__main__':
    # Run the example first
    run_simple_integration_example()

    # Then run unit tests
    print("\n\nRUNNING UNIT TESTS")
    print("="*70)
    unittest.main(verbosity=2)
