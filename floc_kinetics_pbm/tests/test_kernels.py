"""
Unit tests for aggregation and breakage kernels.
"""

import pytest
import numpy as np
from src.kernels import FlocKernels, calculate_floc_kernels
from src.properties import FlocProperties


class TestFlocKernels:
    """Test suite for FlocKernels class."""

    def setup_method(self):
        """Set up test fixtures."""
        self.props = FlocProperties()
        self.kernels = FlocKernels(self.props)
        self.G = 100.0  # 1/s (typical for flocculation)

    def test_volume_diameter_conversion(self):
        """Test volume to diameter and diameter to volume conversions."""
        d = 50e-6  # 50 μm
        v = self.kernels.volume_from_diameter(d)
        d_back = self.kernels.diameter_from_volume(v)

        assert d_back == pytest.approx(d, rel=1e-10)

    def test_shear_dissipation_conversion(self):
        """Test shear rate and dissipation rate conversions."""
        G = 100.0  # 1/s
        epsilon = self.kernels.dissipation_from_shear_rate(G)
        G_back = self.kernels.shear_rate_from_dissipation(epsilon)

        assert G_back == pytest.approx(G, rel=1e-10)
        # G² = ε/ν
        assert epsilon == pytest.approx(G**2 * self.props.nu_water, rel=1e-10)

    def test_orthokinetic_kernel_positive(self):
        """Test that orthokinetic kernel is positive."""
        d_i = 50e-6  # 50 μm
        d_j = 30e-6  # 30 μm

        beta = self.kernels.beta_orthokinetic(d_i, d_j, self.G)

        assert beta > 0

    def test_orthokinetic_kernel_increases_with_shear(self):
        """Test that orthokinetic kernel increases with shear rate."""
        d_i = 50e-6
        d_j = 30e-6

        beta_low = self.kernels.beta_orthokinetic(d_i, d_j, G=10.0)
        beta_high = self.kernels.beta_orthokinetic(d_i, d_j, G=100.0)

        assert beta_high > beta_low

    def test_orthokinetic_kernel_symmetric(self):
        """Test that orthokinetic kernel is symmetric: β(i,j) = β(j,i)."""
        d_i = 50e-6
        d_j = 30e-6

        beta_ij = self.kernels.beta_orthokinetic(d_i, d_j, self.G)
        beta_ji = self.kernels.beta_orthokinetic(d_j, d_i, self.G)

        assert beta_ij == pytest.approx(beta_ji, rel=1e-10)

    def test_orthokinetic_kernel_increases_with_size(self):
        """Test that orthokinetic kernel increases with particle size (cubic dependence)."""
        d_small = 10e-6
        d_large = 100e-6

        beta_small = self.kernels.beta_orthokinetic(d_small, d_small, self.G)
        beta_large = self.kernels.beta_orthokinetic(d_large, d_large, self.G)

        # Should scale roughly as (d_i + d_j)³
        ratio_expected = (2 * d_large / (2 * d_small))**3
        ratio_actual = beta_large / beta_small

        assert ratio_actual == pytest.approx(ratio_expected, rel=0.01)

    def test_perikinetic_kernel_positive(self):
        """Test that perikinetic kernel is positive."""
        d_i = 1e-6  # 1 μm (small particles where Brownian motion matters)
        d_j = 0.5e-6

        beta = self.kernels.beta_perikinetic(d_i, d_j)

        assert beta > 0

    def test_perikinetic_kernel_symmetric(self):
        """Test that perikinetic kernel is symmetric."""
        d_i = 1e-6
        d_j = 0.5e-6

        beta_ij = self.kernels.beta_perikinetic(d_i, d_j)
        beta_ji = self.kernels.beta_perikinetic(d_j, d_i)

        assert beta_ij == pytest.approx(beta_ji, rel=1e-10)

    def test_perikinetic_more_important_for_small_particles(self):
        """Test that perikinetic is relatively more important for small particles."""
        d_small = 1e-6  # 1 μm
        d_large = 50e-6  # 50 μm

        beta_peri_small = self.kernels.beta_perikinetic(d_small, d_small)
        beta_ortho_small = self.kernels.beta_orthokinetic(d_small, d_small, self.G)

        beta_peri_large = self.kernels.beta_perikinetic(d_large, d_large)
        beta_ortho_large = self.kernels.beta_orthokinetic(d_large, d_large, self.G)

        # Perikinetic ratio should be larger for small particles
        ratio_small = beta_peri_small / beta_ortho_small
        ratio_large = beta_peri_large / beta_ortho_large

        assert ratio_small > ratio_large

    def test_differential_sedimentation_kernel_positive(self):
        """Test that differential sedimentation kernel is positive."""
        d_i = 100e-6  # 100 μm (large particles where settling matters)
        d_j = 50e-6

        beta = self.kernels.beta_differential_sedimentation(d_i, d_j)

        assert beta >= 0

    def test_differential_sedimentation_kernel_symmetric(self):
        """Test that differential sedimentation kernel is symmetric."""
        d_i = 100e-6
        d_j = 50e-6

        beta_ij = self.kernels.beta_differential_sedimentation(d_i, d_j)
        beta_ji = self.kernels.beta_differential_sedimentation(d_j, d_i)

        assert beta_ij == pytest.approx(beta_ji, rel=1e-10)

    def test_differential_sedimentation_zero_for_equal_sizes(self):
        """Test that differential sedimentation is zero for equal-sized particles."""
        d = 50e-6

        beta = self.kernels.beta_differential_sedimentation(d, d)

        assert beta == pytest.approx(0.0, abs=1e-20)

    def test_total_aggregation_kernel(self):
        """Test total aggregation kernel combines all mechanisms."""
        d_i = 50e-6
        d_j = 30e-6

        beta_ortho = self.kernels.beta_orthokinetic(d_i, d_j, self.G)
        beta_peri = self.kernels.beta_perikinetic(d_i, d_j)
        beta_ds = self.kernels.beta_differential_sedimentation(d_i, d_j)

        beta_total = self.kernels.beta_total(d_i, d_j, self.G)

        expected = self.props.collision_efficiency * (beta_ortho + beta_peri + beta_ds)
        assert beta_total == pytest.approx(expected, rel=1e-10)

    def test_collision_efficiency_effect(self):
        """Test that collision efficiency scales the aggregation kernel."""
        props_low_alpha = FlocProperties(collision_efficiency=0.1)
        props_high_alpha = FlocProperties(collision_efficiency=0.5)

        kernels_low = FlocKernels(props_low_alpha)
        kernels_high = FlocKernels(props_high_alpha)

        d_i = 50e-6
        d_j = 30e-6

        beta_low = kernels_low.beta_total(d_i, d_j, self.G)
        beta_high = kernels_high.beta_total(d_i, d_j, self.G)

        ratio = beta_high / beta_low
        expected_ratio = 0.5 / 0.1

        assert ratio == pytest.approx(expected_ratio, rel=1e-6)

    def test_breakage_rate_positive(self):
        """Test that breakage rate is positive."""
        d_i = 50e-6

        S = self.kernels.breakage_rate(d_i, self.G)

        assert S >= 0

    def test_breakage_rate_increases_with_shear(self):
        """Test that breakage rate increases with shear rate."""
        # Use larger particle where breakage is more significant
        d_i = 200e-6  # 200 μm

        S_low = self.kernels.breakage_rate(d_i, G=50.0)
        S_high = self.kernels.breakage_rate(d_i, G=500.0)

        # For large particles and high shear, breakage should be observable
        assert S_high > S_low
        assert S_high > 0  # Should have non-zero breakage

    def test_breakage_rate_increases_with_size(self):
        """Test that larger flocs break more easily."""
        # Use high shear and larger size range
        d_small = 100e-6  # 100 μm
        d_large = 400e-6  # 400 μm
        G_high = 500.0  # High shear rate

        S_small = self.kernels.breakage_rate(d_small, G_high)
        S_large = self.kernels.breakage_rate(d_large, G_high)

        # At high shear, larger particles should break more
        assert S_large > S_small
        assert S_large > 0  # Should have non-zero breakage

    def test_daughter_distribution(self):
        """Test daughter particle distribution for binary breakage."""
        v_parent = 1e-15  # m³

        # Daughter volume (half of parent for binary breakage)
        v_daughter = 0.5 * v_parent

        gamma = self.kernels.daughter_distribution(v_parent, v_daughter)

        # Should be positive
        assert gamma > 0

        # For binary breakage, should be 2/v_parent
        assert gamma == pytest.approx(2.0 / v_parent, rel=1e-6)

    def test_daughter_distribution_zero_for_larger_daughter(self):
        """Test that daughter distribution is zero if daughter > parent."""
        v_parent = 1e-15
        v_daughter = 2 * v_parent

        gamma = self.kernels.daughter_distribution(v_parent, v_daughter)

        assert gamma == 0.0


class TestCalculateFlocKernels:
    """Test suite for calculate_floc_kernels function."""

    def test_basic_calculation(self):
        """Test basic kernel calculation."""
        G = 100.0
        d_i = 50e-6
        d_j = 30e-6

        beta, S = calculate_floc_kernels(G, d_i, d_j)

        assert beta > 0
        assert S >= 0

    def test_without_second_particle(self):
        """Test that aggregation kernel is zero when only one particle given."""
        G = 100.0
        d_i = 50e-6

        beta, S = calculate_floc_kernels(G, d_i)

        assert beta == 0.0
        assert S >= 0

    def test_volume_units(self):
        """Test calculation with volume units."""
        G = 100.0
        v_i = 1e-15  # m³
        v_j = 0.5e-15

        beta, S = calculate_floc_kernels(G, v_i, v_j, size_units='volume')

        assert beta > 0
        assert S >= 0

    def test_custom_properties(self):
        """Test calculation with custom properties."""
        props = FlocProperties(collision_efficiency=0.5)
        G = 100.0
        d_i = 50e-6
        d_j = 30e-6

        beta, S = calculate_floc_kernels(G, d_i, d_j, properties=props)

        assert beta > 0
        assert S >= 0

    def test_validation_negative_shear(self):
        """Test that negative shear rate raises error."""
        with pytest.raises(ValueError, match="Shear rate G must be non-negative"):
            calculate_floc_kernels(-1.0, 50e-6, 30e-6)

    def test_validation_negative_size(self):
        """Test that negative particle size raises error."""
        with pytest.raises(ValueError, match="Particle size must be positive"):
            calculate_floc_kernels(100.0, -50e-6, 30e-6)

        with pytest.raises(ValueError, match="Particle size must be positive"):
            calculate_floc_kernels(100.0, 50e-6, -30e-6)

    def test_validation_invalid_units(self):
        """Test that invalid size units raise error."""
        with pytest.raises(ValueError, match="size_units must be 'diameter' or 'volume'"):
            calculate_floc_kernels(100.0, 50e-6, 30e-6, size_units='invalid')

    def test_zero_shear_rate(self):
        """Test behavior at zero shear rate (only Brownian motion)."""
        G = 0.0
        d_i = 1e-6
        d_j = 0.5e-6

        beta, S = calculate_floc_kernels(G, d_i, d_j)

        # Should still have Brownian aggregation
        assert beta > 0
        # No breakage at zero shear
        assert S == pytest.approx(0.0, abs=1e-20)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
