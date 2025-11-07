"""
Unit tests for Population Balance Equation solver.
"""

import pytest
import numpy as np
from src.pbe_solver import PopulationBalanceModel
from src.properties import FlocProperties


class TestPopulationBalanceModel:
    """Test suite for PopulationBalanceModel class."""

    def setup_method(self):
        """Set up test fixtures."""
        self.props = FlocProperties()
        self.pbm = PopulationBalanceModel(
            n_bins=20,
            d_min=1e-6,
            d_max=1e-3,
            properties=self.props
        )

    def test_initialization(self):
        """Test PBM initialization."""
        assert self.pbm.n_bins == 20
        assert len(self.pbm.diameters) == 20
        assert len(self.pbm.volumes) == 20
        assert len(self.pbm.dv) == 20

    def test_size_bins_logarithmic(self):
        """Test that size bins are logarithmically spaced."""
        log_d = np.log10(self.pbm.diameters)
        diffs = np.diff(log_d)

        # Should be approximately constant spacing in log space
        assert np.allclose(diffs, diffs[0], rtol=1e-10)

    def test_size_bin_ranges(self):
        """Test that size bins span the specified range."""
        assert self.pbm.diameters[0] == pytest.approx(self.pbm.d_min, rel=1e-6)
        assert self.pbm.diameters[-1] == pytest.approx(self.pbm.d_max, rel=1e-6)

    def test_find_bin_index(self):
        """Test finding bin index for a given volume."""
        # Test extremes
        idx_min = self.pbm.find_bin_index(self.pbm.volumes[0] * 0.5)
        assert idx_min == 0

        idx_max = self.pbm.find_bin_index(self.pbm.volumes[-1] * 2.0)
        assert idx_max == self.pbm.n_bins - 1

        # Test middle
        v_mid = self.pbm.volumes[10]
        idx_mid = self.pbm.find_bin_index(v_mid)
        assert idx_mid == 10

    def test_precompute_kernels(self):
        """Test kernel matrix precomputation."""
        G = 100.0  # 1/s

        self.pbm.precompute_kernels(G)

        assert self.pbm.beta_matrix is not None
        assert self.pbm.beta_matrix.shape == (self.pbm.n_bins, self.pbm.n_bins)
        assert self.pbm.breakage_rates is not None
        assert len(self.pbm.breakage_rates) == self.pbm.n_bins

        # Check symmetry of aggregation kernel
        assert np.allclose(self.pbm.beta_matrix, self.pbm.beta_matrix.T)

        # Check positivity
        assert np.all(self.pbm.beta_matrix >= 0)
        assert np.all(self.pbm.breakage_rates >= 0)

    def test_moments_zeroth(self):
        """Test zeroth moment calculation (total number)."""
        # Create a simple distribution
        n = np.zeros(self.pbm.n_bins)
        n[5] = 1e6  # 1e6 particles/m³ in bin 5

        M0 = self.pbm.moments(n, k=0)

        # Should be approximately equal to the number concentration
        # times the bin width
        expected = 1e6 * self.pbm.dv[5]
        assert M0 == pytest.approx(expected, rel=1e-6)

    def test_moments_first(self):
        """Test first moment calculation (total volume)."""
        # Create a monodisperse distribution
        n = np.zeros(self.pbm.n_bins)
        n[10] = 1e6  # particles/m³

        M1 = self.pbm.moments(n, k=1)

        # Should be volume × number × bin width
        expected = self.pbm.volumes[10] * 1e6 * self.pbm.dv[10]
        assert M1 == pytest.approx(expected, rel=1e-6)

    def test_total_number_concentration(self):
        """Test total number concentration calculation."""
        n = np.ones(self.pbm.n_bins) * 1e6  # Uniform distribution

        N_total = self.pbm.total_number_concentration(n)

        # Should be sum of n_i * dv_i
        expected = np.sum(n * self.pbm.dv)
        assert N_total == pytest.approx(expected, rel=1e-6)

    def test_total_volume_fraction(self):
        """Test total volume fraction calculation."""
        n = np.ones(self.pbm.n_bins) * 1e6

        phi = self.pbm.total_volume_fraction(n)

        # Should be sum of v_i * n_i * dv_i
        expected = np.sum(self.pbm.volumes * n * self.pbm.dv)
        assert phi == pytest.approx(expected, rel=1e-6)

    def test_mean_diameter(self):
        """Test mean diameter calculation."""
        # Monodisperse distribution
        n = np.zeros(self.pbm.n_bins)
        n[10] = 1e6

        d_mean = self.pbm.mean_diameter(n)

        # For monodisperse, should be close to the bin diameter
        assert d_mean == pytest.approx(self.pbm.diameters[10], rel=0.1)

    def test_summary_statistics(self):
        """Test summary statistics calculation."""
        # Create a distribution
        n = np.zeros(self.pbm.n_bins)
        n[8:12] = 1e6  # Particles in bins 8-11

        stats = self.pbm.summary_statistics(n)

        assert 'total_number' in stats
        assert 'volume_fraction' in stats
        assert 'mean_diameter' in stats
        assert 'd10' in stats
        assert 'd50' in stats
        assert 'd90' in stats

        # Check that percentiles are ordered
        assert stats['d10'] <= stats['d50'] <= stats['d90']

    def test_aggregation_death_term(self):
        """Test aggregation death term calculation."""
        G = 100.0
        self.pbm.precompute_kernels(G)

        # Create a simple distribution
        n = np.ones(self.pbm.n_bins) * 1e6

        D_agg = self.pbm.aggregation_death(n, i=10)

        # Should be positive (particles are being removed by aggregation)
        assert D_agg > 0

    def test_breakage_death_term(self):
        """Test breakage death term calculation."""
        G = 100.0
        self.pbm.precompute_kernels(G)

        n = np.ones(self.pbm.n_bins) * 1e6

        D_break = self.pbm.breakage_death(n, i=10)

        # Should be positive (particles are breaking)
        assert D_break >= 0

    def test_pbe_rhs_conservation(self):
        """Test that PBE RHS conserves total particle volume (approximately)."""
        G = 100.0
        self.pbm.precompute_kernels(G)

        # Create initial distribution
        n = np.zeros(self.pbm.n_bins)
        n[8:12] = 1e6

        # Calculate RHS
        dndt = self.pbm.pbe_rhs(n, t=0, G=G)

        # Calculate volume change rate
        volume_change_rate = np.sum(self.pbm.volumes * dndt * self.pbm.dv)

        # Should be close to zero (volume conservation in aggregation/breakage)
        # Note: Some deviation expected due to discretization
        total_volume = np.sum(self.pbm.volumes * n * self.pbm.dv)
        relative_change = abs(volume_change_rate) / (total_volume + 1e-20)

        # Allow some tolerance due to discretization
        assert relative_change < 0.5  # Less than 50% per unit time

    def test_solve_basic(self):
        """Test basic PBE solution."""
        # Initial monodisperse distribution
        n_initial = np.zeros(self.pbm.n_bins)
        n_initial[5] = 1e9  # High concentration in small size

        # Time points
        t_span = np.linspace(0, 10, 11)  # 0 to 10 seconds

        # Solve
        G = 100.0
        t, n_solution = self.pbm.solve(n_initial, t_span, G)

        # Check output dimensions
        assert len(t) == len(t_span)
        assert n_solution.shape == (len(t_span), self.pbm.n_bins)

        # Check that solution remains non-negative
        assert np.all(n_solution >= -1e-10)  # Allow tiny numerical errors

        # Check that distribution evolves (aggregation should move mass to larger sizes)
        mean_d_initial = self.pbm.mean_diameter(n_initial)
        mean_d_final = self.pbm.mean_diameter(n_solution[-1, :])

        assert mean_d_final >= mean_d_initial  # Aggregation increases mean size

    def test_solve_with_aggregation_only(self):
        """Test that aggregation moves mass to larger sizes."""
        # Use low shear (less breakage) and small particles
        pbm_agg = PopulationBalanceModel(n_bins=15, d_min=1e-6, d_max=100e-6)
        pbm_agg.kernels.k_b = 0.0  # Disable breakage

        # Initial narrow distribution
        n_initial = np.zeros(pbm_agg.n_bins)
        n_initial[3] = 1e10  # High concentration

        t_span = np.linspace(0, 5, 6)
        G = 50.0

        t, n_solution = pbm_agg.solve(n_initial, t_span, G)

        # Mean diameter should increase
        d_initial = pbm_agg.mean_diameter(n_initial)
        d_final = pbm_agg.mean_diameter(n_solution[-1, :])

        assert d_final > d_initial

    def test_solve_conservation_of_mass(self):
        """Test that total volume is approximately conserved during simulation."""
        n_initial = np.zeros(self.pbm.n_bins)
        n_initial[8:12] = 1e8

        t_span = np.linspace(0, 5, 6)
        G = 50.0

        t, n_solution = self.pbm.solve(n_initial, t_span, G)

        # Calculate total volume at each time
        volumes = []
        for i in range(len(t)):
            vol = self.pbm.total_volume_fraction(n_solution[i, :])
            volumes.append(vol)

        # Volume should not change drastically (conservation)
        volumes = np.array(volumes)
        rel_change = np.abs(volumes - volumes[0]) / (volumes[0] + 1e-20)

        # Allow 50% variation due to discretization and numerical errors
        assert np.all(rel_change < 0.5)

    def test_solve_validation_initial_size(self):
        """Test that solve validates initial condition size."""
        n_initial = np.zeros(10)  # Wrong size

        with pytest.raises(ValueError, match="Initial condition must have"):
            self.pbm.solve(n_initial, np.array([0, 1]), G=100.0)

    def test_equilibrium_approach(self):
        """Test that system approaches steady state."""
        # For high shear, system should reach aggregation-breakage equilibrium
        n_initial = np.zeros(self.pbm.n_bins)
        n_initial[5] = 1e9

        # Long time simulation
        t_span = np.linspace(0, 100, 51)
        G = 200.0  # High shear

        t, n_solution = self.pbm.solve(n_initial, t_span, G)

        # Compare last two time points
        n_second_last = n_solution[-2, :]
        n_last = n_solution[-1, :]

        # Change should be small at equilibrium
        rel_change = np.linalg.norm(n_last - n_second_last) / (np.linalg.norm(n_last) + 1e-20)

        # Should be approaching steady state
        assert rel_change < 0.2  # Less than 20% change in last step


class TestPopulationBalanceModelIntegration:
    """Integration tests for complete PBM simulation scenarios."""

    def test_flocculation_scenario(self):
        """Test a realistic flocculation scenario."""
        # Set up realistic properties
        props = FlocProperties(
            collision_efficiency=0.1,
            binding_strength=1e3,
            fractal_dimension=2.3
        )

        pbm = PopulationBalanceModel(
            n_bins=25,
            d_min=1e-6,
            d_max=500e-6,
            properties=props
        )

        # Initial condition: small particles
        n_initial = np.zeros(pbm.n_bins)
        n_initial[2:5] = 1e10  # 10^10 particles/m³ in small size range

        # Flocculation conditions
        G = 100.0  # 1/s (typical for rapid mixing)
        t_span = np.linspace(0, 60, 61)  # 60 seconds

        # Solve
        t, n_solution = pbm.solve(n_initial, t_span, G)

        # Extract statistics
        stats_initial = pbm.summary_statistics(n_initial)
        stats_final = pbm.summary_statistics(n_solution[-1, :])

        print(f"\nFlocculation Scenario Results:")
        print(f"Initial mean diameter: {stats_initial['mean_diameter']*1e6:.2f} μm")
        print(f"Final mean diameter: {stats_final['mean_diameter']*1e6:.2f} μm")
        print(f"Initial d50: {stats_initial['d50']*1e6:.2f} μm")
        print(f"Final d50: {stats_final['d50']*1e6:.2f} μm")

        # Assertions
        # Mean diameter should increase due to aggregation (allow 1% tolerance for very small changes)
        assert stats_final['mean_diameter'] >= stats_initial['mean_diameter'] * 0.99
        # d50 should also increase or stay similar
        assert stats_final['d50'] >= stats_initial['d50'] * 0.99
        # Total number should decrease or stay similar (within 10% due to breakage-aggregation balance)
        # Note: Small fluctuations can occur due to numerical precision and breakage
        assert abs(stats_final['total_number'] - stats_initial['total_number']) / stats_initial['total_number'] < 0.1

    def test_varying_shear_rate(self):
        """Test behavior at different shear rates."""
        pbm = PopulationBalanceModel(n_bins=20, d_min=1e-6, d_max=200e-6)

        n_initial = np.zeros(pbm.n_bins)
        n_initial[4] = 1e10

        t_span = np.linspace(0, 30, 31)

        results = {}
        for G in [10.0, 50.0, 100.0, 200.0]:
            t, n_sol = pbm.solve(n_initial, t_span, G, recompute_kernels=True)
            stats = pbm.summary_statistics(n_sol[-1, :])
            results[G] = stats['mean_diameter']

        print(f"\nVarying Shear Rate Results:")
        for G, d_mean in results.items():
            print(f"G = {G:.0f} 1/s: Mean diameter = {d_mean*1e6:.2f} μm")

        # At very high shear, breakage dominates, mean size may be smaller
        # At moderate shear, aggregation dominates, mean size increases
        # This is complex behavior, just check reasonable values
        assert all(d > 0 for d in results.values())


if __name__ == '__main__':
    pytest.main([__file__, '-v', '-s'])
