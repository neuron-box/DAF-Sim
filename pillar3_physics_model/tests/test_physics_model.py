"""
Comprehensive unit tests for the Core Physics Model (Pillar 3)

Tests cover:
1. Data classes (Particle, Bubble, FluidProperties)
2. DLVO forces (van der Waals, electrostatic)
3. Hydrodynamic forces
4. Trajectory solver
5. Collision efficiency calculation
"""

import unittest
import math
import numpy as np
import sys
import os

# Add the src directory to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))

from pillar3_physics_model import (
    Particle, Bubble, FluidProperties,
    calculate_attachment_efficiency
)
from pillar3_physics_model.dlvo_forces import DLVOForces
from pillar3_physics_model.hydrodynamic_forces import HydrodynamicForces
from pillar3_physics_model.trajectory_solver import TrajectorySolver


class TestParticle(unittest.TestCase):
    """Test Particle class."""

    def test_particle_initialization(self):
        """Test particle creation with valid parameters."""
        particle = Particle(
            diameter=10e-6,
            zeta_potential=-0.025,
            density=2650.0,
            hamaker_constant=5e-21
        )
        self.assertEqual(particle.diameter, 10e-6)
        self.assertEqual(particle.zeta_potential, -0.025)
        self.assertEqual(particle.density, 2650.0)
        self.assertEqual(particle.hamaker_constant, 5e-21)

    def test_particle_properties(self):
        """Test calculated particle properties."""
        particle = Particle(
            diameter=10e-6,
            zeta_potential=-0.025,
            density=2650.0
        )
        self.assertAlmostEqual(particle.radius, 5e-6)
        self.assertAlmostEqual(particle.volume, (4/3) * math.pi * (5e-6)**3, places=20)
        self.assertAlmostEqual(particle.mass, 2650.0 * particle.volume, places=20)

    def test_particle_validation(self):
        """Test particle validation."""
        with self.assertRaises(ValueError):
            Particle(diameter=-10e-6, zeta_potential=-0.025, density=2650.0)

        with self.assertRaises(ValueError):
            Particle(diameter=10e-6, zeta_potential=-0.5, density=2650.0)

        with self.assertRaises(ValueError):
            Particle(diameter=10e-6, zeta_potential=-0.025, density=-2650.0)


class TestBubble(unittest.TestCase):
    """Test Bubble class."""

    def test_bubble_initialization(self):
        """Test bubble creation."""
        bubble = Bubble(diameter=50e-6, zeta_potential=-0.030)
        self.assertEqual(bubble.diameter, 50e-6)
        self.assertEqual(bubble.zeta_potential, -0.030)

    def test_bubble_rise_velocity(self):
        """Test bubble rise velocity calculation."""
        bubble = Bubble(diameter=50e-6, zeta_potential=-0.030)
        fluid = FluidProperties.water_at_20C()

        v_rise = bubble.calculate_rise_velocity(
            fluid.dynamic_viscosity,
            fluid.density
        )

        # Check that rise velocity is positive and reasonable
        self.assertGreater(v_rise, 0)
        self.assertLess(v_rise, 0.1)  # Should be < 10 cm/s for small bubbles


class TestFluidProperties(unittest.TestCase):
    """Test FluidProperties class."""

    def test_fluid_initialization(self):
        """Test fluid properties initialization."""
        fluid = FluidProperties.water_at_20C(ionic_strength=0.001)
        self.assertAlmostEqual(fluid.temperature, 293.15)
        self.assertAlmostEqual(fluid.ionic_strength, 0.001)

    def test_debye_length(self):
        """Test Debye length calculation."""
        fluid = FluidProperties.water_at_20C(ionic_strength=0.001)
        debye_length = fluid.debye_length

        # For I = 0.001 M, Debye length should be ~9.6 nm
        expected = 0.304 / math.sqrt(0.001) * 1e-9
        self.assertAlmostEqual(debye_length, expected, places=12)

    def test_thermal_voltage(self):
        """Test thermal voltage calculation."""
        fluid = FluidProperties.water_at_20C()
        V_T = fluid.thermal_voltage

        # At 20°C, thermal voltage should be ~25.3 mV
        self.assertAlmostEqual(V_T, 0.0253, places=3)


class TestDLVOForces(unittest.TestCase):
    """Test DLVO force calculations."""

    def setUp(self):
        """Set up test fixtures."""
        self.particle = Particle(
            diameter=10e-6,
            zeta_potential=-0.025,
            density=2650.0,
            hamaker_constant=5e-21
        )
        self.bubble = Bubble(diameter=50e-6, zeta_potential=-0.030)
        self.fluid = FluidProperties.water_at_20C(ionic_strength=0.001)
        self.dlvo = DLVOForces(self.particle, self.bubble, self.fluid)

    def test_van_der_waals_force(self):
        """Test van der Waals force calculation."""
        h = 10e-9  # 10 nm separation

        F_vdw = self.dlvo.van_der_waals_force(h)

        # van der Waals force should be attractive (negative)
        self.assertLess(F_vdw, 0)

        # Force should increase (become more negative) as separation decreases
        F_vdw_close = self.dlvo.van_der_waals_force(5e-9)
        self.assertLess(F_vdw_close, F_vdw)

    def test_electrostatic_force(self):
        """Test electrostatic force calculation."""
        h = 10e-9  # 10 nm separation

        F_edl = self.dlvo.electrostatic_force(h)

        # For like-charged particles (both negative), force should be repulsive (positive)
        self.assertGreater(F_edl, 0)

        # Force should decrease as separation increases
        F_edl_far = self.dlvo.electrostatic_force(20e-9)
        self.assertLess(F_edl_far, F_edl)

    def test_total_dlvo_force(self):
        """Test total DLVO force."""
        h = 10e-9

        F_total = self.dlvo.total_dlvo_force(h)
        F_vdw = self.dlvo.van_der_waals_force(h)
        F_edl = self.dlvo.electrostatic_force(h)

        # Total should be sum of components
        self.assertAlmostEqual(F_total, F_vdw + F_edl, places=15)

    def test_interaction_energy(self):
        """Test interaction energy calculation."""
        h = 10e-9

        V_vdw, V_edl, V_total = self.dlvo.interaction_energy(h)

        # van der Waals energy should be negative (attractive)
        self.assertLess(V_vdw, 0)

        # Electrostatic energy should be positive (repulsive for like charges)
        self.assertGreater(V_edl, 0)

        # Total should be sum
        self.assertAlmostEqual(V_total, V_vdw + V_edl, places=15)


class TestHydrodynamicForces(unittest.TestCase):
    """Test hydrodynamic force calculations."""

    def setUp(self):
        """Set up test fixtures."""
        self.particle = Particle(
            diameter=10e-6,
            zeta_potential=-0.025,
            density=2650.0
        )
        self.bubble = Bubble(diameter=50e-6, zeta_potential=-0.030)
        self.fluid = FluidProperties.water_at_20C()
        self.hydro = HydrodynamicForces(self.particle, self.bubble, self.fluid)

    def test_stokes_drag_force(self):
        """Test Stokes drag force."""
        velocity = 0.001  # 1 mm/s

        F_drag = self.hydro.stokes_drag_force(velocity)

        # Drag should oppose motion (negative for positive velocity)
        self.assertLess(F_drag, 0)

        # Drag magnitude should be proportional to velocity
        F_drag_2v = self.hydro.stokes_drag_force(2 * velocity)
        self.assertAlmostEqual(abs(F_drag_2v), 2 * abs(F_drag), places=15)

    def test_bubble_flow_velocity(self):
        """Test bubble flow field calculation."""
        R_b = self.bubble.radius
        r = 2 * R_b  # 2 bubble radii from center
        theta = math.pi / 4  # 45 degrees

        v_r, v_theta = self.hydro.bubble_flow_velocity(r, theta)

        # Velocities should be on the order of bubble rise velocity
        U = self.bubble.rise_velocity or self.bubble.calculate_rise_velocity(
            self.fluid.dynamic_viscosity, self.fluid.density
        )
        self.assertLess(abs(v_r), 2 * U)
        self.assertLess(abs(v_theta), 2 * U)

    def test_reynolds_number(self):
        """Test Reynolds number calculation."""
        velocity = 0.001  # 1 mm/s

        Re = self.hydro.reynolds_number(velocity)

        # For small particles at low velocities, Re should be << 1
        self.assertLess(Re, 1.0)

    def test_gravitational_force(self):
        """Test gravitational force."""
        F_g = self.hydro.gravitational_force()

        # For particles denser than water, force should be positive (downward)
        self.assertGreater(F_g, 0)


class TestTrajectorySolver(unittest.TestCase):
    """Test trajectory solver."""

    def setUp(self):
        """Set up test fixtures."""
        self.particle = Particle(
            diameter=5e-6,
            zeta_potential=-0.020,
            density=2650.0,
            hamaker_constant=5e-21
        )
        self.bubble = Bubble(diameter=50e-6, zeta_potential=-0.030)
        self.fluid = FluidProperties.water_at_20C(ionic_strength=0.001)
        self.solver = TrajectorySolver(self.particle, self.bubble, self.fluid)

    def test_trajectory_solver_initialization(self):
        """Test trajectory solver initialization."""
        self.assertIsNotNone(self.solver)
        self.assertIsNotNone(self.solver.dlvo)
        self.assertIsNotNone(self.solver.hydro)

    def test_collision_detection(self):
        """Test collision detection."""
        contact_distance = self.particle.radius + self.bubble.radius

        # Just at contact
        self.assertTrue(self.solver.collision_detected(contact_distance))

        # Well beyond contact
        self.assertFalse(self.solver.collision_detected(contact_distance * 2))

    def test_solve_trajectory_centerline(self):
        """Test trajectory solution for centerline approach."""
        # Particle starting on centerline should collide
        result = self.solver.solve_trajectory(initial_offset=0.0, max_time=0.1)

        self.assertTrue(result.success)
        # For centerline, collision should occur (in most cases)
        # This may depend on the specific force balance


class TestCollisionEfficiency(unittest.TestCase):
    """Test collision efficiency calculation."""

    def test_calculate_attachment_efficiency_basic(self):
        """Test basic collision efficiency calculation."""
        particle = Particle(
            diameter=10e-6,
            zeta_potential=-0.025,
            density=2650.0,
            hamaker_constant=5e-21
        )
        bubble = Bubble(diameter=50e-6, zeta_potential=-0.030)
        fluid = FluidProperties.water_at_20C(ionic_strength=0.001)

        result = calculate_attachment_efficiency(particle, bubble, fluid)

        # Check that all expected keys are present
        self.assertIn('alpha_bp', result)
        self.assertIn('eta_collision', result)
        self.assertIn('alpha_attachment', result)
        self.assertIn('critical_offset', result)
        self.assertIn('energy_barrier', result)
        self.assertIn('energy_barrier_kT', result)
        self.assertIn('debye_length', result)

        # Check that values are in reasonable ranges
        self.assertGreaterEqual(result['alpha_bp'], 0.0)
        self.assertLessEqual(result['alpha_bp'], 1.0)
        self.assertGreaterEqual(result['eta_collision'], 0.0)
        self.assertLessEqual(result['eta_collision'], 1.0)
        self.assertGreaterEqual(result['alpha_attachment'], 0.0)
        self.assertLessEqual(result['alpha_attachment'], 1.0)

    def test_opposite_charges(self):
        """Test collision efficiency with opposite charges."""
        # Particle positive, bubble negative -> should have high efficiency
        particle = Particle(
            diameter=10e-6,
            zeta_potential=0.025,  # Positive
            density=2650.0,
            hamaker_constant=8e-21  # Higher Hamaker for hydrophobic
        )
        bubble = Bubble(diameter=50e-6, zeta_potential=-0.030)
        fluid = FluidProperties.water_at_20C(ionic_strength=0.001)

        result = calculate_attachment_efficiency(particle, bubble, fluid)

        # Opposite charges should result in favorable attachment
        self.assertTrue(result['attachment_favorable'])
        self.assertGreater(result['alpha_attachment'], 0.5)

    def test_like_charges_high_ionic_strength(self):
        """Test collision efficiency with like charges at high ionic strength."""
        particle = Particle(
            diameter=10e-6,
            zeta_potential=-0.025,
            density=2650.0
        )
        bubble = Bubble(diameter=50e-6, zeta_potential=-0.030)
        # High ionic strength reduces electrostatic repulsion
        fluid = FluidProperties.water_at_20C(ionic_strength=0.1)

        result = calculate_attachment_efficiency(particle, bubble, fluid)

        # Higher ionic strength should reduce energy barrier
        self.assertLess(result['energy_barrier_kT'], 20.0)

    def test_size_dependence(self):
        """Test that collision efficiency depends on particle size."""
        bubble = Bubble(diameter=50e-6, zeta_potential=-0.030)
        fluid = FluidProperties.water_at_20C(ionic_strength=0.001)

        # Small particle
        particle_small = Particle(
            diameter=5e-6,
            zeta_potential=-0.025,
            density=2650.0
        )
        result_small = calculate_attachment_efficiency(particle_small, bubble, fluid)

        # Large particle
        particle_large = Particle(
            diameter=20e-6,
            zeta_potential=-0.025,
            density=2650.0
        )
        result_large = calculate_attachment_efficiency(particle_large, bubble, fluid)

        # Collision efficiencies should be different
        self.assertNotEqual(result_small['eta_collision'], result_large['eta_collision'])


class TestIntegration(unittest.TestCase):
    """Integration tests with realistic scenarios."""

    def test_typical_daf_conditions(self):
        """Test under typical DAF conditions."""
        # Typical DAF scenario: coagulated particles with reduced zeta potential
        particle = Particle(
            diameter=15e-6,  # Flocculated particle
            zeta_potential=-0.010,  # Reduced by coagulant
            density=1200.0,  # Lower density due to floc structure
            hamaker_constant=8e-21  # Hydrophobic after coagulation
        )

        bubble = Bubble(
            diameter=40e-6,  # Typical DAF bubble
            zeta_potential=-0.020
        )

        fluid = FluidProperties.water_at_20C(
            ionic_strength=0.01,  # After coagulant addition
            ph=7.0
        )

        result = calculate_attachment_efficiency(particle, bubble, fluid)

        # Under good DAF conditions, efficiency should be reasonable
        self.assertGreater(result['alpha_bp'], 0.01)
        self.assertLess(result['energy_barrier_kT'], 30.0)

        # Print results for inspection
        print("\n=== Typical DAF Conditions ===")
        print(f"α_bp = {result['alpha_bp']:.4f}")
        print(f"η_collision = {result['eta_collision']:.4f}")
        print(f"α_attachment = {result['alpha_attachment']:.4f}")
        print(f"Energy barrier = {result['energy_barrier_kT']:.1f} kT")
        print(f"Critical offset = {result['critical_offset']*1e6:.2f} μm")
        print(f"Debye length = {result['debye_length']*1e9:.2f} nm")

    def test_poor_coagulation(self):
        """Test with poor coagulation (high zeta potential)."""
        # Poor coagulation: high zeta potential with lower Hamaker constant
        # (more hydrophilic particles)
        particle = Particle(
            diameter=10e-6,
            zeta_potential=-0.040,  # High magnitude
            density=2650.0,
            hamaker_constant=1e-21  # Lower (more hydrophilic)
        )

        bubble = Bubble(diameter=50e-6, zeta_potential=-0.035)
        fluid = FluidProperties.water_at_20C(ionic_strength=0.0001)  # Lower ionic strength

        result = calculate_attachment_efficiency(particle, bubble, fluid)

        # Check that results are in valid range
        self.assertGreaterEqual(result['alpha_bp'], 0.0)
        self.assertLessEqual(result['alpha_bp'], 1.0)

        # Attachment may be unfavorable
        print("\n=== Poor Coagulation Conditions ===")
        print(f"α_bp = {result['alpha_bp']:.4f}")
        print(f"Energy barrier = {result['energy_barrier_kT']:.1f} kT")
        print(f"Attachment favorable: {result['attachment_favorable']}")


def run_tests():
    """Run all unit tests."""
    # Create test suite
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()

    # Add all test classes
    suite.addTests(loader.loadTestsFromTestCase(TestParticle))
    suite.addTests(loader.loadTestsFromTestCase(TestBubble))
    suite.addTests(loader.loadTestsFromTestCase(TestFluidProperties))
    suite.addTests(loader.loadTestsFromTestCase(TestDLVOForces))
    suite.addTests(loader.loadTestsFromTestCase(TestHydrodynamicForces))
    suite.addTests(loader.loadTestsFromTestCase(TestTrajectorySolver))
    suite.addTests(loader.loadTestsFromTestCase(TestCollisionEfficiency))
    suite.addTests(loader.loadTestsFromTestCase(TestIntegration))

    # Run tests
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)

    return result


if __name__ == '__main__':
    result = run_tests()
    sys.exit(0 if result.wasSuccessful() else 1)
