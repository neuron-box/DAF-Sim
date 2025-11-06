"""
Unit tests for floc properties module.
"""

import pytest
import numpy as np
from src.properties import FlocProperties


class TestFlocProperties:
    """Test suite for FlocProperties class."""

    def test_default_initialization(self):
        """Test that default properties are initialized correctly."""
        props = FlocProperties()

        assert props.rho_water == 998.0
        assert props.mu_water == 1.002e-3
        assert props.temperature == 293.15
        assert props.fractal_dimension == 2.3
        assert 0 <= props.collision_efficiency <= 1

    def test_custom_properties(self):
        """Test initialization with custom properties."""
        props = FlocProperties(
            rho_water=1000.0,
            mu_water=1.0e-3,
            temperature=298.15,
            collision_efficiency=0.5
        )

        assert props.rho_water == 1000.0
        assert props.mu_water == 1.0e-3
        assert props.temperature == 298.15
        assert props.collision_efficiency == 0.5

    def test_kinematic_viscosity(self):
        """Test kinematic viscosity calculation."""
        props = FlocProperties(rho_water=1000.0, mu_water=1.0e-3)
        nu = props.nu_water

        assert nu == pytest.approx(1.0e-6, rel=1e-6)

    def test_floc_density_decreases_with_size(self):
        """Test that floc density decreases with increasing size (fractal structure)."""
        props = FlocProperties(fractal_dimension=2.3)

        d_small = 1e-6  # 1 μm
        d_large = 100e-6  # 100 μm

        rho_small = props.floc_density(d_small)
        rho_large = props.floc_density(d_large)

        # Larger flocs should have lower density
        assert rho_small > rho_large
        # Density should be between water and floc material density
        assert props.rho_water <= rho_large <= rho_small <= props.rho_floc

    def test_floc_mass(self):
        """Test floc mass calculation."""
        props = FlocProperties()
        d = 50e-6  # 50 μm

        mass = props.floc_mass(d)

        # Mass should be positive
        assert mass > 0
        # Order of magnitude check
        volume = np.pi / 6.0 * d**3
        assert mass < props.rho_floc * volume  # Less than solid sphere

    def test_kolmogorov_scales(self):
        """Test Kolmogorov length and time scale calculations."""
        props = FlocProperties()
        epsilon = 0.1  # m²/s³ (typical for flocculation)

        eta = props.kolmogorov_length_scale(epsilon)
        tau = props.kolmogorov_time_scale(epsilon)

        # Reasonable values for flocculation
        assert 1e-6 < eta < 1e-3  # Between 1 μm and 1 mm
        assert 0.001 < tau < 1.0  # Between 1 ms and 1 s

    def test_validation_positive_density(self):
        """Test that negative density raises error."""
        with pytest.raises(ValueError, match="Water density must be positive"):
            FlocProperties(rho_water=-1.0)

    def test_validation_positive_viscosity(self):
        """Test that negative viscosity raises error."""
        with pytest.raises(ValueError, match="Dynamic viscosity must be positive"):
            FlocProperties(mu_water=-1.0)

    def test_validation_collision_efficiency_range(self):
        """Test that collision efficiency must be between 0 and 1."""
        with pytest.raises(ValueError, match="Collision efficiency must be between 0 and 1"):
            FlocProperties(collision_efficiency=1.5)

        with pytest.raises(ValueError, match="Collision efficiency must be between 0 and 1"):
            FlocProperties(collision_efficiency=-0.1)

    def test_validation_fractal_dimension_range(self):
        """Test that fractal dimension must be between 1 and 3."""
        with pytest.raises(ValueError, match="Fractal dimension must be between 1 and 3"):
            FlocProperties(fractal_dimension=3.5)

        with pytest.raises(ValueError, match="Fractal dimension must be between 1 and 3"):
            FlocProperties(fractal_dimension=0.5)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
